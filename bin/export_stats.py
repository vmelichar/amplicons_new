"""
Merge filtering tsv file from different chunks
"""

import pandas as pd
import argparse
import sys
import os


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
        "-i", "--input", dest="INPUT", required=True, help="Input folder"
    )
    parser.add_argument(
        "-l", "--low_clusters", dest="LOW_CL", required=True, nargs="+", help="Counts for low clusters"
    )
    parser.add_argument(
        "-x", "--hs_index", dest="HS_IDX", required=True, nargs="+", help="HS index for counts of low clusters"
    )

    args = parser.parse_args(argv)

    return args

# BASIC STATS

def get_basic_stats_to_dict(barcode, input):

    file = f'{input}/barcode01/HS{barcode}/stats/raw/filtering_stats.tsv'
    
    if os.path.exists(file):
        with open(file) as f:
            lines = f.readlines()

        count_line = lines[1]

        counts = count_line.split('\t')

        dict = {    'reads_found': counts[2],
                    'reads_unmapped': counts[3],
                    'reads_secondary': counts[4],
                    'reads_supplementary': counts[5],
                    'reads_on_target': counts[6],
                    'reads_concatamer': counts[7],
                    'reads_short': counts[8],
                    'reads_long': counts[9],
                    'reads_filtered': counts[10]}
        return dict

    else:
        dict = {    'reads_found': 0,
                    'reads_unmapped': 0,
                    'reads_secondary': 0,
                    'reads_supplementary': 0,
                    'reads_on_target': 0,
                    'reads_concatamer': 0,
                    'reads_short': 0,
                    'reads_long': 0,
                    'reads_filtered': 0}
        return dict

def get_basic_stats_to_dict_perc(dict, part=False):
    output_dict = {}

    if not part:
        for key, value in dict.items():
            output_dict[key] = count_perc_string(value, dict['reads_found'])
        return output_dict

    high_level_keys = ['reads_found', 'reads_unmapped', 'reads_secondary', 'reads_supplementary', 'reads_on_target']
    low_level_keys = ['reads_concatamer', 'reads_short', 'reads_long', 'reads_filtered']
    for key,value in dict.items():
        if key in high_level_keys:
            output_dict[key] = count_perc_string(value, dict['reads_found'])
        if key in low_level_keys:
            output_dict[key] = count_perc_string(value, dict['reads_on_target'])
    return output_dict       


def get_basic_strand_stats_to_dict(barcode, strand, input):
    
    file = f'{input}/barcode01/HS{barcode}/stats/raw/strand_filter.txt'

    if os.path.exists(file):
        with open(file) as f:
            lines = f.readlines()

        lines = [int(i) for i in lines]
        
        if strand == '+':
            counts = lines[0::2]
            target = sum(counts)
            dict = {    'reads_found': '-',
                        'reads_unmapped': '-',
                        'reads_secondary': '-',
                        'reads_supplementary': '-',
                        'reads_on_target': target,
                        'reads_concatamer': counts[0],
                        'reads_short': counts[1],
                        'reads_long': counts[2],
                        'reads_filtered': counts[3]}
            return dict
        if strand == '-':
            counts = lines[1::2]
            target = sum(counts)
            dict = {    'reads_found': '-',
                        'reads_unmapped': '-',
                        'reads_secondary': '-',
                        'reads_supplementary': '-',
                        'reads_on_target': target,
                        'reads_concatamer': counts[0],
                        'reads_short': counts[1],
                        'reads_long': counts[2],
                        'reads_filtered': counts[3]}
            return dict
    
    else:
        dict = {    'reads_found': '-',
                    'reads_unmapped': '-',
                    'reads_secondary': '-',
                    'reads_supplementary': '-',
                    'reads_on_target': 0,
                    'reads_concatamer': 0,
                    'reads_short': 0,
                    'reads_long': 0,
                    'reads_filtered': 0}
        return dict


def get_detected_umis_stats(barcode, df, input):
    hs = f'hs{barcode}'
    file = f'{input}/barcode01/HS{barcode}/stats/raw/detected_umi_stats.tsv'
    file_strands = f'{input}/barcode01/HS{barcode}/stats/raw/umi_strand_filter.txt'

    if os.path.exists(file):

        with open(file, 'r') as f:
            line = f.readlines()[1]

        reads_used = int(line.split('\t')[6].strip())
        reads_used_perc_part = f'{round(float(line.split('\t')[7].strip()), 1)}%'
        reads_used_perc_tot = f"{round(reads_used / int(df.at['reads_found', f'{hs}_count']) * 100, 1)}%"
        
        if os.path.exists(file_strands):

            with open(file_strands, 'r') as f:
                lines = f.readlines()

            plus = int(lines[0])
            minus = int(lines[1])

            return {f'hs{barcode}_count': reads_used, 
                    f'hs{barcode}_percPart': reads_used_perc_part, 
                    f'hs{barcode}_percTot': reads_used_perc_tot,
                    f'hs{barcode}_plus': plus, 
                    f'hs{barcode}_minus': minus}
        else:

            return {f'hs{barcode}_count': reads_used, 
                    f'hs{barcode}_percPart': reads_used_perc_part, 
                    f'hs{barcode}_percTot': reads_used_perc_tot,
                    f'hs{barcode}_plus': '--',
                    f'hs{barcode}_minus': '--'}
    else:
        return {f'hs{barcode}_count': 0, 
                f'hs{barcode}_percPart': 0, 
                f'hs{barcode}_percTot': 0,
                f'hs{barcode}_plus': 0, 
                f'hs{barcode}_minus': 0}


def get_csv_stats(input):
    df = pd.DataFrame()
    for barcode in range(1,9):
        df = pd.concat([df, pd.DataFrame(get_basic_stats_to_dict(barcode, input), index=[f'hs{barcode}_count'])])
        df = pd.concat([df, pd.DataFrame(get_basic_stats_to_dict_perc(get_basic_stats_to_dict(barcode, input), part=True), index=[f'hs{barcode}_percPart'])])
        df = pd.concat([df, pd.DataFrame(get_basic_stats_to_dict_perc(get_basic_stats_to_dict(barcode, input), part=False), index=[f'hs{barcode}_percTot'])])
        df = pd.concat([df, pd.DataFrame(get_basic_strand_stats_to_dict(barcode, '+', input), index=[f'hs{barcode}_plus'])])
        df = pd.concat([df, pd.DataFrame(get_basic_strand_stats_to_dict(barcode, '-', input), index=[f'hs{barcode}_minus'])])
        
    df = df.transpose()

    df_umi = pd.DataFrame()
    for barcode in range(1,9):
        df_umi = pd.merge(df_umi, pd.DataFrame(get_detected_umis_stats(barcode, df, input), index=['detected_both_umis']), left_index=True, right_index=True, how='outer')
    
    df = pd.concat([df, df_umi])

    new_cols = []
    for col in df.columns:
        values = col.split('_')
        new_cols.append((values[0], values[1]))

    df.columns = pd.MultiIndex.from_tuples(new_cols, names=['Hotspot','Type_of_Value'])

    df.to_csv('basic_stats.csv')

    return df


# EXTRACTION STATS

def get_extraction_stats_to_list(hs, mode, extraction, input, func=['med','mean']):
    file = f'{input}/barcode01/HS{hs[2]}/stats/{mode}/extraction_{extraction}_stats.tsv'

    if os.path.exists(file):
        df = pd.read_csv(file, sep='\t', skiprows=1, header=None)
        df.columns = ['patt_syn', 'seq_syn', 'patt_umi', 'seq_umi', 'orient', 'strand', 'ed', 'len', 'strart', 'end']

        fwd_plus = df[(df['strand'] == '+') & (df['orient'] == 'fwd') & (df['ed'] != -1)][['ed', 'len', 'strart', 'end']]
        rev_plus = df[(df['strand'] == '+') & (df['orient'] == 'rev') & (df['ed'] != -1)][['ed', 'len', 'strart', 'end']]
        fwd_minus = df[(df['strand'] == '-') & (df['orient'] == 'fwd') & (df['ed'] != -1)][['ed', 'len', 'strart', 'end']]
        rev_minus = df[(df['strand'] == '-') & (df['orient'] == 'rev') & (df['ed'] != -1)][['ed', 'len', 'strart', 'end']]

        if func == 'med':
            fwd_plus = [int(i) if pd.notna(i) else 0 for i in fwd_plus.median().tolist()]
            rev_plus = [int(i) if pd.notna(i) else 0 for i in rev_plus.median().tolist()]
            fwd_minus = [int(i) if pd.notna(i) else 0 for i in fwd_minus.median().tolist()]
            rev_minus = [int(i) if pd.notna(i) else 0 for i in rev_minus.median().tolist()]
        elif func == 'mean':
            fwd_plus = [round(i,2) for i in fwd_plus.mean().tolist()]
            rev_plus = [round(i,2) for i in rev_plus.mean().tolist()]
            fwd_minus = [round(i,2) for i in fwd_minus.mean().tolist()]
            rev_minus = [round(i,2) for i in rev_minus.mean().tolist()]

        return fwd_plus + rev_plus + fwd_minus + rev_minus
    else:
        return [0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0]


def get_csv_extraction(input):
    row_tupples = []
    columns_tupples = []

    for hs in ['hs1', 'hs2', 'hs3', 'hs4', 'hs5', 'hs6', 'hs7', 'hs8']:
        for v in ['synthetic', 'umi']:
            row_tupples.append((hs,v))

    for strand in ['fwd', 'rev']:
        for umi in ['head', 'tail']:
            for v in ['ED', 'Len', 'S', 'E']:
                columns_tupples.append((strand,umi,v))

    # Define the multiindex for rows and columns
    rows = pd.MultiIndex.from_tuples(row_tupples, names=['Hotspot', 'Extraction'])
    columns = pd.MultiIndex.from_tuples(columns_tupples, names=['Strand', 'Orientation', 'Value'])

    # Create an empty dataframe with the multiindex rows and columns
    df_med = pd.DataFrame(index=rows, columns=columns)
    df_med_cons = pd.DataFrame(index=rows, columns=columns)

    for hs, extr in df_med.index:
        df_to_add = pd.DataFrame(get_extraction_stats_to_list(hs,'raw',extr,input,'med')).T
        df_to_add.columns = columns
        df_med.loc[(hs,extr)] = df_to_add.loc[0]

    for hs, extr in df_med_cons.index:
        df_to_add = pd.DataFrame(get_extraction_stats_to_list(hs,'consensus',extr,input,'med')).T
        df_to_add.columns = columns
        df_med_cons.loc[(hs,extr)] = df_to_add.loc[0]

    # Create an empty dataframe with the multiindex rows and columns
    df_mean = pd.DataFrame(index=rows, columns=columns)
    df_mean_cons = pd.DataFrame(index=rows, columns=columns)

    for hs, extr in df_mean.index:
        df_to_add = pd.DataFrame(get_extraction_stats_to_list(hs,'raw',extr,input,'mean')).T
        df_to_add.columns = columns
        df_mean.loc[(hs,extr)] = df_to_add.loc[0]

    for hs, extr in df_mean_cons.index:
        df_to_add = pd.DataFrame(get_extraction_stats_to_list(hs,'consensus',extr,input,'mean')).T
        df_to_add.columns = columns
        df_mean_cons.loc[(hs,extr)] = df_to_add.loc[0]

    df_med.to_csv('extraction_stats_med.csv')
    df_mean.to_csv('extraction_stats_mean.csv')
    df_med_cons.to_csv('extraction_stats_med_cons.csv')
    df_mean_cons.to_csv('extraction_stats_mean_cons.csv')

    return df_med, df_mean, df_med_cons, df_mean_cons


# EXTRACTION BASICS STATS

def get_extraction_basics_to_list(hs, mode, extraction, val, input):
    file = f'{input}/barcode01/HS{hs[2]}/stats/{mode}/extraction_{extraction}_stats.tsv'

    if os.path.exists(file):
        df = pd.read_csv(file, sep='\t', skiprows=1, header=None)
        df.columns = ['patt_syn', 'seq_syn', 'patt_umi', 'seq_umi', 'orient', 'strand', 'ed', 'len', 'strart', 'end']

        total = len(df)
        total_perc = count_perc_string(total, total)
        
        if extraction == 'synthetic':
            no_match = len(df[df.ed == -1])
            no_match_perc = count_perc_string(no_match, total)
            extracted = len(df[df.ed > -1])
            extracted_perc = count_perc_string(extracted, total)

        if extraction == 'umi':
            no_match = len(df[df.seq_umi.isnull()])
            no_match_perc = count_perc_string(no_match, total)
            extracted = len(df[(df.ed > -1) & (df.seq_umi.notnull())])
            extracted_perc = count_perc_string(extracted, total)
        
        if val == 'count':
            return [total, no_match, extracted]
        if val == 'perc':
            return [total_perc, no_match_perc, extracted_perc]
    else:
        return [0,0,0]


def get_csv_extraction_basics(input):
    columns_tupples = []
    rows = ['total', 'no_match', 'extracted']

    for hs in ['hs1', 'hs2', 'hs3', 'hs4', 'hs5', 'hs6', 'hs7', 'hs8']:
        for v in ['count', 'perc']:
            columns_tupples.append((hs,v))
        
    # Define the multiindex for columns
    columns = pd.MultiIndex.from_tuples(columns_tupples, names=['Hotspot', 'Value'])

    # Create an empty dataframe with the multiindex rows and columns
    df_synthetic = pd.DataFrame(index=rows, columns=columns)
    df_umi = pd.DataFrame(index=rows, columns=columns)
    df_synthetic_cons = pd.DataFrame(index=rows, columns=columns)
    df_umi_cons = pd.DataFrame(index=rows, columns=columns)

    for hs, val in df_synthetic.columns:
        df_to_add = pd.DataFrame(get_extraction_basics_to_list(hs, 'raw', 'synthetic', val, input))
        df_to_add.index = rows
        df_synthetic.loc[:,(hs,val)] = df_to_add[0]

    for hs, val in df_synthetic_cons.columns:
        df_to_add = pd.DataFrame(get_extraction_basics_to_list(hs, 'consensus', 'synthetic', val, input))
        df_to_add.index = rows
        df_synthetic_cons.loc[:,(hs,val)] = df_to_add[0]

    for hs, val in df_umi.columns:
        df_to_add = pd.DataFrame(get_extraction_basics_to_list(hs, 'raw', 'umi', val, input))
        df_to_add.index = rows
        df_umi.loc[:,(hs,val)] = df_to_add[0]

    for hs, val in df_umi_cons.columns:
        df_to_add = pd.DataFrame(get_extraction_basics_to_list(hs, 'consensus', 'umi', val, input))
        df_to_add.index = rows
        df_umi_cons.loc[:,(hs,val)] = df_to_add[0]

    df_synthetic.to_csv('extraction_basics_synthetic.csv')
    df_umi.to_csv('extraction_basics_umi.csv')

    df_synthetic_cons.to_csv('extraction_basics_synthetic_cons.csv')
    df_umi_cons.to_csv('extraction_basics_umi_cons.csv')

    return df_synthetic, df_umi, df_synthetic_cons, df_umi_cons


# CLUSTER STATS

def get_cluster_stats_to_list(hs, input, func=['med','mean']):
    file = f'{input}/barcode01/HS{hs[2]}/stats/raw/split_cluster_stats.tsv'

    if os.path.exists(file):
        df = pd.read_csv(file, sep='\t', header=0)

        df_0 = df[df['cluster_written'] == 0][['reads_found', 'reads_found_fwd', 'reads_found_rev']]
        df_1 = df[df['cluster_written'] == 1][['reads_found', 'reads_found_fwd', 'reads_found_rev']]
        df_tot = df[['reads_found', 'reads_found_fwd', 'reads_found_rev']]

        if func == 'med':
            df_0_list = df_0.median().tolist()
            df_1_list = df_1.median().tolist()
            df_0_list = [int(i) if pd.notna(i) else 0 for i in df_0_list]
            df_1_list = [int(i) if pd.notna(i) else 0 for i in df_1_list]
        elif func == 'mean':
            df_0_list = df_0.mean().tolist()
            df_1_list = df_1.mean().tolist()
            df_0_list = [round(i,2) for i in df_0_list]
            df_1_list = [round(i,2) for i in df_1_list]

        file_det = f'{input}/barcode01/HS{hs[2]}/stats/raw/detected_umi_stats.tsv'

        with open(file_det, 'r') as f_det:
            line_det = f_det.readlines()[1]

        sum_0_cl = len(df_0)
        sum_1_cl = len(df_1)
        sum_tot_cl = len(df_tot)
        perc_1_cl = count_perc_string(sum_1_cl, sum_tot_cl)
        sum_tot_rd = df_tot['reads_found'].sum()
        sum_1_rd = df_1['reads_found'].sum()
        rd_per_cl = round(sum_tot_rd / sum_tot_cl, 1)
        rd_det = int(line_det.split('\t')[6].strip())
        perc_det = count_perc_string(sum_tot_rd, rd_det)
        perc_used = count_perc_string(sum_1_rd, rd_det)

        df_tot = [sum_tot_cl, sum_0_cl, sum_1_cl, perc_1_cl, sum_tot_rd, rd_per_cl, perc_det, perc_used]

        return df_tot + df_0_list + df_1_list
    else:
        return [0,0,0,0,0,0,0,0, 0,0,0, 0,0,0]


def get_csv_cluster(input):
    columns_tupples = []

    rows = ['hs1', 'hs2', 'hs3', 'hs4', 'hs5', 'hs6', 'hs7', 'hs8']

    for tot in ['Tot_clusters', '0Tot_clusters', '1Tot_clusters', 'perc1Tot_clusters', 'Tot_reads', 'reads_per_cluster', 'perc_of_detected', 'perc_used_polishing']:
        columns_tupples.append(('Total',tot))
    for written in ['0', '1']:
        for v in ['reads_found', 'reads_found_plus', 'reads_found_minus']:
            columns_tupples.append((written,v))

    # Define the multiindex for rows and columns
    columns = pd.MultiIndex.from_tuples(columns_tupples, names=['Written', 'Value'])

    # Create an empty dataframe with the multiindex rows and columns
    df_med = pd.DataFrame(index=rows, columns=columns)

    for hs in df_med.index:
        df_to_add = pd.DataFrame(get_cluster_stats_to_list(hs,input,'med')).T
        df_to_add.columns = columns
        df_med.loc[(hs)] = df_to_add.loc[0]

    # Create an empty dataframe with the multiindex rows and columns
    df_mean = pd.DataFrame(index=rows, columns=columns)

    for hs in df_mean.index:
        df_to_add = pd.DataFrame(get_cluster_stats_to_list(hs,input,'mean')).T
        df_to_add.columns = columns
        df_mean.loc[(hs)] = df_to_add.loc[0]

    df_med.to_csv('cluster_stats_med.csv')
    df_mean.to_csv('cluster_stats_mean.csv')

    return df_med, df_mean

# SANKEY

def get_sankey_values(hs,input,low_clusters,hs_index):

    filter_values = get_basic_stats_to_dict(hs, input)

    file_detection = f'{input}/barcode01/HS{hs}/stats/raw/detected_umi_stats.tsv'

    if os.path.exists(file_detection):
        with open(file_detection, 'r') as f:
            line = f.readlines()[1]

        det_umi = int(line.split('\t')[6].strip())

        undetected_umi = int(filter_values['reads_filtered']) - det_umi
    else:
        det_umi = 0

        undetected_umi = int(filter_values['reads_filtered']) - 0

    file_clustering = f'{input}/barcode01/HS{hs}/stats/raw/split_cluster_stats.tsv'
    
    if os.path.exists(file_clustering):
        df = pd.read_csv(file_clustering, sep='\t', header=0)

        cl_0 = int(df[df['cluster_written'] == 0]['reads_found'].sum())
        cl_1 = int(df[df['cluster_written'] == 1]['reads_found'].sum())

        hs_idx = hs_index.index(f'HS{hs}')
        low_cluster_count = int(low_clusters[hs_idx])

        singletons = det_umi - low_cluster_count - cl_0 - cl_1

        recombo = df[df['cluster_id'].isin(get_recombo_clusters_names(hs,input))]['reads_found'].sum()
    else:
        cl_0 = 0
        cl_1 = 0
        if f'HS{hs}' not in hs_index:
            low_cluster_count = 0
        else:
            hs_idx = hs_index.index(f'HS{hs}')
            low_cluster_count = int(low_clusters[hs_idx])
        singletons = det_umi - low_cluster_count - cl_0 - cl_1
        recombo = 0

    return [
            filter_values['reads_unmapped'], 
            filter_values['reads_secondary'],
            filter_values['reads_supplementary'],
            filter_values['reads_on_target'],
            filter_values['reads_concatamer'],
            filter_values['reads_short'],
            filter_values['reads_long'],
            filter_values['reads_filtered'],
            undetected_umi,
            det_umi,
            singletons,
            low_cluster_count,
            cl_0,
            cl_1,
            recombo
            ]


def get_sankey_all(input,low_clusters,hs_index):

    sources = [0, 0, 0, 0, 4, 4, 4, 4, 8, 8, 10, 10, 10, 10, 14]
    targets = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]

    for hs in [1,2,3,4,5,6,7,8]:
        values = get_sankey_values(hs, input, low_clusters, hs_index)

        df = pd.DataFrame({'source': sources, 'target': targets, 'value': values})
        df.to_csv(f'sankey_HS{hs}.csv', index=False)


# MISCELLANEOUS

def get_recombo_clusters_names(hs, input):
    file = f'{input}/barcode01/HS{hs}/analysis/analysis_HS{hs}_perc_counts.csv'
    if os.path.exists(file):
        df = pd.read_csv(file)
        names = [i.split('=')[1].replace('_', '_sub') for i in list(df['read'])]

        return names
    else:
        return []


def check_strands(df):
    n_err = 0
    for index, row in df.iterrows():
        if index in ['reads_on_target','reads_concatamer', 'reads_short', 'reads_long', 'reads_filtered', 'detected_both_umis']:
            for i in ['1','2','3', '4', '5', '6', '7', '8']:
                if int(row.loc[(f'hs{i}', 'count')]) != int(row.loc[(f'hs{i}', 'plus')]) + int(row.loc[(f'hs{i}', 'minus')]):
                    n_err += 1
                    print('Strands do not match!')
    if n_err == 0:
        print('No errors in strand counts.')


def count_perc_string(x,y):
    perc = round(int(x) / int(y) * 100, 1)
    return f'{perc}%'


def main(argv=sys.argv[1:]):
    args = parse_args(argv=argv)

    print('Parsing stats...')


    print('Parsing filtering stats...')
    df_basic = get_csv_stats(args.INPUT)
    check_strands(df_basic)

    print('Parsing extraction stats...')
    df_extraction_med, df_extraction_mean, _, _ = get_csv_extraction(args.INPUT)
    df_extraction_basics_synthetic, df_extraction_basics_umi, _, _ = get_csv_extraction_basics(args.INPUT)

    print('Parsing cluster stats...')
    df_cluster_med, df_cluster_mean = get_csv_cluster(args.INPUT)

    print('Sankey values calculation...')
    get_sankey_all(args.INPUT, args.LOW_CL, args.HS_IDX)


if __name__ == '__main__':
    main()