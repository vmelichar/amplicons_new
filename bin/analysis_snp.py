import sys
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
        return None
    
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
                if quality >= 15:
                    return base
                else:
                    return "N"
            else:
                # position has been deleted
                return "D"

def get_bases(bam, vcf_file, tbi_file, pos_file,out_dir):
    # Load VCF file
    vcf_snps = pysam.VariantFile(vcf_file, index_filename=tbi_file)

    # Open the positions file
    bases = []

    with open(pos_file, "r") as positions_file:
        for line in positions_file:
            chrom, pos = line.strip().split()
            pos = int(pos)

            for rec in vcf_snps.fetch('chr' + chrom, pos-1, pos):
                ref, alt = rec.alleles[0], rec.alleles[1]

            for read in bam.fetch('chr' + chrom, pos-1, pos+1):
                if not read.is_supplementary and not read.is_unmapped and not read.is_secondary:
                    # minimap2 always store forward strand (even when original sequence is aligned to reverse strand)
                    base = find_base_in_alignment(read, pos, bam_stores_revcomp = True) #if find_base_in_alignment(read, pos, bam_stores_revcomp = True) is not None else ("0")

                    if base == ref:
                        base = "B"
                    elif base == alt:
                        base = "P"
                    elif base.islower():
                        base == 'M'
                    elif base in ['A', 'C', 'T', 'G']:
                        base = 'X'
                    bases.append({
                        "read": read.query_name,
                        "position": f'{chrom}:{pos}',
                        "base": base,
                        })
                else: continue

    # Prepare output
    df = pd.DataFrame(bases)

    df_long = pd.pivot_table(df, 
                         index='read', 
                         columns='position', 
                         values='base', 
                         aggfunc=lambda x: ' '.join(x), 
                         fill_value='.')
    df_long.to_csv(out_dir + '_counts.csv', 
                   sep=',', index=True)
    return df_long


def get_ratios(row):
   occurances_dict = row.value_counts().to_dict()
   B = occurances_dict.get('B', 0)
   P = occurances_dict.get('P', 0)
   N = occurances_dict.get('N', 0)
   dot = occurances_dict.get('.', 0)
   SNPs = len(row) - dot
   # 1 - B6 strand, 0 - PWD strand
   # should I count others as well?
   if B == 0:
      ratioBP = 0
   else:
      ratioBP = B / (B + P)
   ratioN = N / SNPs
   #other = SNPs - B - P - N
   return pd.Series([round(ratioBP,2), round(ratioN,2)]) 

def create_barplot(percentages, counts, out_dir, chimera_perc, chimeras):
    plt.figure(figsize=(10, 6))
    plt.bar(percentages, counts, width=0.05, align='center')
    plt.xlim(0, 1)
    plt.ylim(0, max(counts) + 100)
    plt.xlabel('Percentage')
    plt.ylabel('No. of reads')
    plt.suptitle('B6 ratio',  fontsize = 18)
    plt.title('Number of reads between [0.2, 0.8]:' + str(chimeras) + ' out of ' + str(sum(counts)) + ' (' + str(chimera_perc) + '%)', fontsize = 12)
    plt.savefig(out_dir + '_barplot.png', bbox_inches = 'tight')
    plt.close()
    
def create_lineplot(perc_counts, col, out_dir):
    sns.lineplot(x = 'Percentage', y = 'Count', data = perc_counts)

    # Set plot labels and title
    plt.xlabel('Percentage')
    plt.ylabel('Count')
    if col == 'perc_B':
        middle = sum(perc_counts[perc_counts.Percentage.between(0.2, 0.8, "both")]['Count'])
        chimera_perc = round(middle/sum(perc_counts['Count'])*100, 2)
        
        # Create a bar plot using Matplot
        print('Plotting barplot...')
        create_barplot(perc_counts['Percentage'], perc_counts['Count'], out_dir, chimera_perc, middle)
        
        plt.suptitle('B6 ratio',  fontsize = 18)
        plt.title('Number of reads between [0.2, 0.8]:' + str(middle) + ' out of ' + str(sum(perc_counts['Count'])) + ' (' + str(chimera_perc) + '%)', fontsize = 12)
    elif col == 'perc_N':
        plt.suptitle('Percentage of low quality SNPs',  fontsize = 18)
    plt.xticks(rotation=45)  # Rotate x-axis labels for better readability
    plt.xlim(0, 1)
    plt.ylim(0, max(perc_counts['Count']) + 100)

    # Save plot
    plt.savefig(out_dir + '_plot_'+ col +'.png', bbox_inches='tight')
    plt.close()

def get_BP_graph(df, col, out_dir):
    # Count occurrences of each rounded percentage
    percentage_counts = df[col].value_counts().reset_index()
    percentage_counts.columns = ['Percentage', 'Count']

    df[df.perc_B.between(0.2, 0.8, "both")].to_csv(out_dir + '_perc_counts.csv', 
                   sep=',', index=True)

    # Create lineplot
    print('Plotting lineplot...')
    create_lineplot(percentage_counts, col, out_dir)

def get_rest(df, out_dir):
    df[df.perc_B.between(0, 0.2, "neither")].to_csv(out_dir + '_perc_counts_upto2.csv',
                   sep=',', index=True)

    df[df.perc_B.between(0.8, 1, "neither")].to_csv(out_dir + '_perc_counts_from8.csv',
                   sep=',', index=True)
    
    df[df.perc_B == 0].to_csv(out_dir + '_perc_counts_allPWD.csv',
                   sep=',', index=True)

    df[df.perc_B == 1].to_csv(out_dir + '_perc_counts_allB6.csv',
                   sep=',', index=True)

    with open(out_dir + '_seqs_PWD.txt', 'w') as pwd_file:
        for name in list(df[df.perc_B == 0].index):
            print(name, file=pwd_file)

    with open(out_dir + '_seqs_B6.txt', 'w') as b6_file:
        for name in list(df[df.perc_B == 1].index):
            print(name, file=b6_file)


def run_pipeline(input_bam, input_bai, vcf, tbi, positions, output_dir):
    
    bam = pysam.AlignmentFile(input_bam, "rb", index_filename=input_bai)

    df = get_bases(bam, vcf, tbi, positions, output_dir)
    df[['perc_B', 'perc_N']] = df.apply(get_ratios, axis=1)

    get_BP_graph(df, 'perc_B', output_dir)
    get_BP_graph(df, 'perc_N', output_dir)
    get_rest(df, output_dir)


def main(argv=sys.argv[1:]):
    """
    Basic command line interface.
    """
    args = parse_args(argv=argv)

    run_pipeline(args.BAM, args.BAI, args.VCF, args.TBI, args.POS, args.OUTPUT)

if __name__ == "__main__":
    main()