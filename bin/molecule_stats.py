import sys
import argparse
import pandas as pd

def parse_args(argv):
    """
    Commandline parser

    :param argv: Command line arguments
    :type argv: List
    """
    usage = "Command line interface to molecule statistics"
    parser = argparse.ArgumentParser(
        description=usage, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        "--hotspot", dest="HOTSPOT", required=True, help="Hotspot identifier"
    )
    parser.add_argument(
        "--extraction_raw_syn_file", dest="EXTRACTION_RAW_SYN_FILE", required=True, help="Input extraction raw synthetic tsv file"
    )
    parser.add_argument(
        "--extraction_raw_umi_file", dest="EXTRACTION_RAW_UMI_FILE", required=True, help="Input extraction raw umi tsv file"
    )
    parser.add_argument(
        "--extraction_cons_syn_file", dest="EXTRACTION_CONS_SYN_FILE", required=True, help="Input extraction consensus synthetic tsv file"
    )
    parser.add_argument(
        "--extraction_cons_umi_file", dest="EXTRACTION_CONS_UMI_FILE", required=True, help="Input extraction consensus umi tsv file"
    )
    parser.add_argument(
        "--cluster_hash_file", dest="CLUSTER_HASH_FILE", required=True, help="Input cluster hash tsv file"
    )
    parser.add_argument(
        "--seq_type_hash_file", dest="SEQ_TYPE_HASH_FILE", required=True, help="Input sequence type hash tsv file"
    )
    parser.add_argument(
        "-o", "--output", dest="OUTPUT", required=True, help="Output folder"
    )
    args = parser.parse_args(argv)

    return args


def main(argv=sys.argv[1:]):
    args = parse_args(argv)

    # Load input files
    extraction_raw_syn_df = pd.read_csv(args.EXTRACTION_RAW_SYN_FILE, sep='\t', header=0)
    extraction_raw_syn_df.drop(columns=['umi_pattern', 'umi_seq'], inplace=True)
    extraction_raw_syn_df = extraction_raw_syn_df.pivot(index='hash', columns='orientation')
    extraction_raw_syn_df.columns = [f"{col[0]}_{col[1]}" for col in extraction_raw_syn_df.columns]
    extraction_raw_syn_df.reset_index(inplace=True)
    extraction_raw_syn_df.drop(columns=['strand_rev'], inplace=True)
    extraction_raw_syn_df.rename(columns={'strand_fwd': 'strand'}, inplace=True)

    extraction_raw_umi_df = pd.read_csv(args.EXTRACTION_RAW_UMI_FILE, sep='\t', header=0)
    extraction_raw_umi_df.drop(columns=['syn_pattern', 'syn_seq'], inplace=True)
    extraction_raw_umi_df = extraction_raw_umi_df.pivot(index='hash', columns='orientation')
    extraction_raw_umi_df.columns = [f"{col[0]}_{col[1]}" for col in extraction_raw_umi_df.columns]
    extraction_raw_umi_df.reset_index(inplace=True)
    extraction_raw_umi_df.drop(columns=['umi_pattern_rev', 'strand_fwd', 'strand_rev'], inplace=True)
    extraction_raw_umi_df.rename(columns={'umi_pattern_fwd': 'umi_pattern'}, inplace=True)

    extraction_cons_syn_df = pd.read_csv(args.EXTRACTION_CONS_SYN_FILE, sep='\t', header=0)
    extraction_cons_syn_df.drop(columns=['umi_pattern', 'umi_seq'], inplace=True)
    extraction_cons_syn_df = extraction_cons_syn_df.pivot(index='hash', columns='orientation')
    extraction_cons_syn_df.columns = [f"{col[0]}_{col[1]}" for col in extraction_cons_syn_df.columns]
    extraction_cons_syn_df.reset_index(inplace=True)
    extraction_cons_syn_df.drop(columns=['strand_rev'], inplace=True)
    extraction_cons_syn_df.rename(columns={'strand_fwd': 'strand'}, inplace=True)

    extraction_cons_umi_df = pd.read_csv(args.EXTRACTION_CONS_UMI_FILE, sep='\t', header=0)
    extraction_cons_umi_df.drop(columns=['syn_pattern', 'syn_seq'], inplace=True)
    extraction_cons_umi_df = extraction_cons_umi_df.pivot(index='hash', columns='orientation')
    extraction_cons_umi_df.columns = [f"{col[0]}_{col[1]}" for col in extraction_cons_umi_df.columns]
    extraction_cons_umi_df.reset_index(inplace=True)
    extraction_cons_umi_df.drop(columns=['umi_pattern_fwd', 'umi_pattern_rev', 'strand_fwd', 'strand_rev'], inplace=True)

    cluster_hash_df = pd.read_csv(args.CLUSTER_HASH_FILE, sep='\t', header=0)
    seq_type_hash_df = pd.read_csv(args.SEQ_TYPE_HASH_FILE, sep='\t', header=0)
    seq_type_hash_df.columns = ['cluster_id', 'seq_type']

    # Merge extraction raw
    print('Merging extraction raw data...')
    merged_extr_raw = pd.merge(
        extraction_raw_syn_df, 
        extraction_raw_umi_df,
        on='hash',
        how='left',
        suffixes=('_raw_syn', '_raw_umi')
    )

    extraction_cons_syn_df['cluster_id'] = extraction_cons_syn_df['hash'].apply(lambda x: x.split('_')[0])
    extraction_cons_syn_df.drop(columns=['hash'], inplace=True)
    extraction_cons_umi_df['cluster_id'] = extraction_cons_umi_df['hash'].apply(lambda x: x.split('_')[0])
    extraction_cons_umi_df.drop(columns=['hash'], inplace=True)
    
    # Merge extraction consensus
    print('Merging extraction consensus data...')
    merged_extr_cons = pd.merge(
        extraction_cons_syn_df, 
        extraction_cons_umi_df,
        on='cluster_id',
        how='left',
        suffixes=('_cons_syn', '_cons_umi')
    )

    # Merge with cluster hash
    print('Merging with cluster hash data...')
    merged_raw_ext_cluster = pd.merge(
        merged_extr_raw,
        cluster_hash_df,
        on='hash',
        how='left'
    )

    # Merge raw extraction and cluster with consensus extraction
    print('Merging all data together...')
    merged_cluster = pd.merge(
        merged_raw_ext_cluster,
        merged_extr_cons,
        on='cluster_id',
        how='left',
        suffixes=('_raw', '_cons')
    )

    print(merged_cluster.columns)
    print(seq_type_hash_df.columns)

    # Merge with seq type hash
    print('Merging with sequence type hash data...')
    merged_final = pd.merge(
        merged_cluster,
        seq_type_hash_df,
        on='cluster_id',
        how='left'
    )

    merged_final = merged_final.fillna('')

    # Put hash, cluster_id, seq_type at the front
    cols = merged_final.columns.tolist()
    cols.insert(0, cols.pop(cols.index('seq_type')))
    cols.insert(0, cols.pop(cols.index('cluster_id')))
    cols.insert(0, cols.pop(cols.index('hash')))
    merged_data = merged_final[cols]

    # Save the processed data to output file
    output_file = f"{args.OUTPUT}/{args.HOTSPOT}_molecule_statistics.tsv"
    merged_data.to_csv(output_file, sep='\t', index=False)
    print(f"Molecule statistics saved to {output_file}")


if __name__ == "__main__":
    main()