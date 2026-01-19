import sys
import argparse
import pysam
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from typing import Optional


def parse_args(argv):
    """
    Commandline parser

    :param argv: Command line arguments
    :type argv: List
    """
    usage = "Command line interface to telemap"
    parser = argparse.ArgumentParser(
        description=usage, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        "--hs", dest="HS", required=True, help="Hotspot"
    )
    parser.add_argument(
        "--bam", dest="BAM", required=True, help="Input .bam file"
    )
    parser.add_argument(
        "--bai", dest="BAI", required=True, help="Input .bai file"
    )
    parser.add_argument(
        "--vcf", dest="VCF", required=True, help="VCF file with called variants"
    )
    parser.add_argument(
        "--tbi", dest="TBI", required=True, help="Index to .vcf file"
    )
    parser.add_argument(
        "--positions", dest="POS", required=True, help="Positions file"
    )
    parser.add_argument(
        "-o", "--output", dest="OUTPUT", required=True, help="Output folder"
    )

    args = parser.parse_args(argv)

    return args


def find_base_in_alignment(alignment: pysam.AlignedSegment,
                           pos: int, 
                           bam_stores_revcomp: bool = False) -> Optional[str]:
    idx_q = 0
    idx_r = pos - alignment.reference_start
    if bam_stores_revcomp:
        seq = alignment.query_sequence
    else:
        seq = alignment.get_forward_sequence()
    
    if seq is None:
        return [None, None]
    
    # Use CIGAR string from bam file
    # cigar alignment is represented by list of tuples (operation,length)
    for op, l in alignment.cigartuples:
        ref_consumed = op in {0, 2, 3, 7, 8} # matches, deletions, ref skip, equal or diff
        query_consumed = op in {0, 1, 4, 7, 8} # matches, insertions, soft clip, equal or diff
        
        if ref_consumed:
            idx_r -= l
        if query_consumed:
            idx_q += l
        
        if idx_r <= 0:
            if query_consumed:
                # base is in query between idx_q-l , idx_q
                base = seq[idx_q + idx_r - 1]

                quality = alignment.query_qualities[idx_q + idx_r - 1]
                if quality >= 20:
                    return [base, quality]
                else:
                    return [base, quality]
            else:
                # position has been deleted
                return ["D", None]

def get_bases(bam, vcf_file, tbi_file, pos_file,out_dir):
    # Load VCF file
    vcf_snps = pysam.VariantFile(vcf_file, index_filename=tbi_file)

    # Open the positions file
    bases = []

    with open(pos_file, "r") as positions_file:
        for line in positions_file:
            chrom, pos = line.strip().split()
            pos = int(pos)

            vcf_records = list(vcf_snps.fetch('chr' + chrom, pos-1, pos))
            if len(vcf_records) != 1:
                print(f'VCF records for pos {chrom}:{pos} = {len(vcf_records)}')

            for rec in vcf_records:
                if len(rec.alleles) != 2:
                    print(f'More Alt Alleles for {chrom}:{pos} -- {rec.alleles}')
                ref, alt = rec.alleles[0], rec.alleles[1]

            for read in bam.fetch('chr' + chrom, pos-1, pos+1):
                if not read.is_supplementary and not read.is_unmapped and not read.is_secondary:
                    # minimap2 always store forward strand (even when original sequence is aligned to reverse strand)
                    base, qual = find_base_in_alignment(read, pos, bam_stores_revcomp = True) #if find_base_in_alignment(read, pos, bam_stores_revcomp = True) is not None else ("0")

                    if base == ref:
                        base = "B"
                    elif base == alt:
                        base = "P"
                    elif base in ['A', 'C', 'T', 'G']:
                        base = 'X'
                    if base not in ['B', 'P', 'X', 'N', 'D']:
                        print(f'Unknown assigments in the sequence: {base}')
                    bases.append({
                        "read": read.query_name,
                        "position": f'{chrom}:{pos}',
                        "base": base,
                        "qual": qual
                        })
                else:
                    continue

    # Prepare output
    df = pd.DataFrame(bases)

    df_long = pd.pivot_table(df, 
                         index='read', 
                         columns='position', 
                         values=['base','qual'], 
                         aggfunc=lambda x: '$'.join(x.astype(str)), 
                         fill_value='.')
    df_long.columns = [f"{level0}_{level1}" for level0, level1 in df_long.columns]
    df_long.to_csv(out_dir + '_counts.csv', 
                   sep=',', index=True)

    B_seqs = filtered_df = df_long.loc[(df_long.filter(regex='^base') == 'B').all(axis=1)].index.tolist()
    
    return df_long, B_seqs


def get_penalties(bam, B_seqs):
    total_len = 0
    D_len = 0
    N_len = 0

    for read in bam.fetch():
        if read.query_name not in B_seqs:
            continue
        if read.is_secondary or read.is_supplementary:
            continue
        else:
            cigar = read.cigartuples
            length = read.query_alignment_length
            quals = read.query_qualities

            D = 0
            for op, l in cigar:
                if op == 2:
                    D += l
            
            frontS = cigar[0][1] if cigar[0][0] == '4' else 0
            rearS = cigar[-1][1] if cigar[-1][0] == '4' else None
            N = 0
            for q in quals[frontS:rearS]:
                if q < 20:
                    N += 1

            total_len += length
            D_len += D
            N_len += N
    
    D_pen = D_len / total_len
    N_pen = N_len / total_len

    print(f'\tNumber of seqs contr to penalties: {len(B_seqs)}')
    print(f'\tT/D/N - {total_len}/{D_len}/{N_len}')
    print(f'\tD penalty: {D_pen}')
    print(f'\tN penalty: {N_pen}')

    return D_pen, N_pen


def get_ratios(row, D_pen, N_pen):
    base_cols = [c for c in row.index if str(c).startswith('base')]
    qual_cols = [c for c in row.index if str(c).startswith('qual')]

    base_sr = row[base_cols]
    qual_sr = pd.to_numeric(row[qual_cols], errors='coerce')

    occurances_dict = base_sr.value_counts().to_dict()
   
    B = occurances_dict.get('B', 0)
    P = occurances_dict.get('P', 0)
    N = (qual_sr < 20).sum()
    X = occurances_dict.get('X', 0)
    M = occurances_dict.get('M', 0)
    D = occurances_dict.get('D', 0)
    dot = occurances_dict.get('.', 0)
    SNPs = len(base_sr) - dot
    # 1 - B6 strand, 0 - PWD strand
    # should I count others as well?
    if B == 0:
        ratioBP = 0
    else:
        ratioBP = B / (B + P)
    ratioN = N / SNPs
    ratioD = D / SNPs
   
    minority_allele = 'P' if P < B else 'B'
    majority_allele = 'P' if P >= B else 'B'

#    # Calculate separate freq-based penalties
#    s_n = 0.75 * ( ratioN ** 3 )
#    s_d = 0.8 * ( ratioD ** 3 )

#    # Initial prob incorporating N and D penalties
#    prob_all_correct = ((1.0 - s_n) ** N) * ((1.0 - s_d) ** D)

#    for idx, allele in enumerate(base_sr):
#       if allele in ['N', 'D']:
#         continue
#       if allele == minority_allele or allele in ['X', 'M']:
#         q_score = int(qual_sr.iloc[idx])
#         q_factor = 10 ** (-q_score / 10.0)

#         if allele in ['P', 'B']:
#             s_i = (1 / 3.0) * q_factor
#         if allele == 'X':
#             s_i = 0.75 * q_factor
#         if allele == 'M':
#             s_i = 0.5 * q_factor

#         prob_all_correct *= (1.0 - s_i)

    # Calc separate D and N penalties
    s_n = N_pen
    s_d = D_pen

    # Initial prob incorporating N and D penalties
    prob_all_correct = ((1.0 - s_n) ** N) * ((1.0 - s_d) ** D)

    for idx, allele in enumerate(base_sr):
        if allele in ['N', 'D']:
            continue
        if allele in ['P', 'B', 'X']:
            q_score = int(qual_sr.iloc[idx])
            q_factor = 10 ** (-q_score / 10.0)

            if allele == minority_allele:
                s_i = (1 / 3.0) * q_factor
            elif allele == majority_allele:
                s_i = q_factor
            elif allele == 'X':
                s_i = (2 / 3.0) * q_factor
            else:
                print(f'Unknown allele: {allele}, minor {minority_allele}, q_score {q_score}')

            prob_all_correct *= (1.0 - s_i)

    error = 1.0 - prob_all_correct

    return pd.Series([round(ratioBP,2), round(ratioN,2), round(error, 5), int(B), int(P), int(N), int(X), int(M), int(D)]) 


def create_barplot(percentages, counts, out_dir, chimera_perc, chimeras, hs):
    plt.figure(figsize=(7, 5))
    plt.bar(percentages, counts, width=0.05, align='center')
    plt.xlim(0, 1)
    plt.ylim(0, max(counts) + 100)
    plt.xlabel('Percentage')
    plt.ylabel('No. of clusters')
    plt.suptitle(f'{hs} - B6 ratio',  fontsize = 18)
    plt.title('Number of reads between [0.2, 0.8]:' + str(chimeras) + ' out of ' + str(sum(counts)) + ' (' + str(chimera_perc) + '%)', fontsize = 12)
    plt.savefig(out_dir + '_barplot.png', bbox_inches = 'tight')
    plt.close()
    
def create_lineplot(perc_counts, col, out_dir, hs):
    plt.figure(figsize=(7, 5))
    sns.lineplot(x = 'Percentage', y = 'Count', data = perc_counts)

    # Set plot labels and title
    plt.xlabel('Percentage')
    plt.ylabel('No. of clusters')
    if col == 'perc_B':
        middle = sum(perc_counts[perc_counts.Percentage.between(0.2, 0.8, "both")]['Count'])
        chimera_perc = round(middle/sum(perc_counts['Count'])*100, 2)
        
        # Create a bar plot using Matplot
        print('Plotting barplot...')
        create_barplot(perc_counts['Percentage'], perc_counts['Count'], out_dir, chimera_perc, middle, hs)
        
        plt.suptitle(f'{hs} - B6 ratio',  fontsize = 18)
        plt.title('Number of reads between [0.2, 0.8]:' + str(middle) + ' out of ' + str(sum(perc_counts['Count'])) + ' (' + str(chimera_perc) + '%)', fontsize = 12)
    elif col == 'perc_N':
        plt.suptitle('Percentage of low quality SNPs',  fontsize = 18)
    plt.xticks(rotation=45)  # Rotate x-axis labels for better readability
    plt.xlim(0, 1)
    plt.ylim(0, max(perc_counts['Count']) + 100)

    # Save plot
    plt.savefig(out_dir + '_plot_'+ col +'.png', bbox_inches='tight')
    plt.close()

def get_BP_graph(df, col, out_dir, hs):
    # Count occurrences of each rounded percentage
    percentage_counts = df[col].value_counts().reset_index()
    percentage_counts.columns = ['Percentage', 'Count']

    # Create lineplot
    print('Plotting lineplot...')
    create_lineplot(percentage_counts, col, out_dir, hs)

def get_tables(df, out_dir):
    f_stats = open(out_dir + '_stats.txt', 'w')
    print(f'Total seqs: {len(df)}', file=f_stats)
    print(f'Error seqs: {len(df[df.Err > 0.05])}', file=f_stats)
    print(f'High confidence seqs: {len(df[df.Err <= 0.05])}', file=f_stats)

    df[df.Err > 0.05].to_csv(out_dir + '_errors.csv', sep='\t', index=True)

    # Filter for Error threshold
    df = df[df.Err <= 0.05]
    
    df[df.perc_B == 0].to_csv(out_dir + '_allPWD.csv', sep='\t', index=True)
    print(f'PWD seqs: {len(df[df.perc_B == 0])}', file=f_stats)

    df[df.perc_B == 1].to_csv(out_dir + '_allB6.csv', sep='\t', index=True)
    print(f'B6 seqs: {len(df[df.perc_B == 1])}', file=f_stats)

    df[df.perc_B.between(0, 1, "neither")].to_csv(out_dir + '_recombo.csv', sep='\t', index=True)
    print(f'Recombo seqs: {len(df[df.perc_B.between(0, 1, "neither")])}', file=f_stats)

    f_stats.close()

    with open(out_dir + '_seqs_PWD.txt', 'w') as pwd_file:
        for name in list(df[df.perc_B == 0].index):
            print(name, file=pwd_file)

    with open(out_dir + '_seqs_B6.txt', 'w') as b6_file:
        for name in list(df[df.perc_B == 1].index):
            print(name, file=b6_file)

def get_cluster_type(df, out_dir):
    df['seq_type'] = 'X'
    df.loc[df.Err > 0.05, 'seq_type'] = 'Low_Confidence'
    df.loc[(df.Err <= 0.05) & (df.perc_B == 1), 'seq_type'] = 'Complete_B6'
    df.loc[(df.Err <= 0.05) & (df.perc_B == 0), 'seq_type'] = 'Complete_PWD'
    df.loc[(df.Err <= 0.05) & (df.perc_B.between(0, 1, 'neither')), 'seq_type'] = 'Recombo'

    df.index = df.index.map(lambda x: x.split('=')[1].split('_')[0])  # Simplify cluster names

    df['seq_type'].to_csv(out_dir + '_cluster_types.tsv', sep='\t', index=True)                                    


def run_pipeline(hs, input_bam, input_bai, vcf, tbi, positions, output_dir):
    
    print('Reading BAM...')
    bam = pysam.AlignmentFile(input_bam, "rb", index_filename=input_bai)

    print('Getting bases...')
    df, B_seqs = get_bases(bam, vcf, tbi, positions, output_dir)

    print('Getting penalties...')
    D_pen, N_pen = get_penalties(bam, B_seqs)

    df[['perc_B', 'perc_N', 'Err', 'cB', 'cP', 'cN', 'cX', 'cM', 'cD']] = df.apply(get_ratios, axis=1, args=(D_pen, N_pen))
    int_cols = ['cB', 'cP', 'cN', 'cX', 'cM', 'cD']
    df[int_cols] = df[int_cols].astype(int)

    positions = [col.split('_')[1] for col in df.columns if col.startswith('base_')]
    for pos in positions:
        base_col = f'base_{pos}'
        qual_col = f'qual_{pos}'
        df[qual_col] = pd.to_numeric(df[qual_col], errors='coerce')
        df.loc[df[qual_col] < 20, base_col] = 'N'

    base_cols = [c for c in df.columns if str(c).startswith('base')]
    qual_cols = [c for c in df.columns if str(c).startswith('qual')]
    df[qual_cols] = df[qual_cols].fillna(0).astype(int)
    df['assign'] = df[base_cols].astype(str).agg(''.join, axis=1)
    df['quals'] = df[qual_cols].astype(str).agg(';'.join, axis=1)
    df.drop(columns=base_cols, inplace=True)
    df.drop(columns=qual_cols, inplace=True)
    #df.columns = [str(c).split('_')[1] if str(c).startswith('base') else str(c) for c in df.columns]

    get_BP_graph(df, 'perc_B', output_dir, hs)
    get_BP_graph(df, 'perc_N', output_dir, hs)

    get_tables(df, output_dir)
    get_cluster_type(df, output_dir)

def main(argv=sys.argv[1:]):
    """
    Basic command line interface.
    """
    args = parse_args(argv=argv)

    out_file = f'{args.OUTPUT}/analysis_{args.HS}'

    run_pipeline(args.HS, args.BAM, args.BAI, args.VCF, args.TBI, args.POS, out_file)

if __name__ == "__main__":
    main()
