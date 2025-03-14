"""
Merge stats tsv file from different chunks
"""

import argparse
import os
import sys


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
        "--detected_tsv", dest="DETECTED_TSV", nargs="+", required=True, type=str, help="TSV reports of detected UMIs"
    )
    parser.add_argument(
        "--synthetic_tsv", dest="SYNTHETIC_TSV", nargs="+", required=True, type=str, help="TSV reports of synthetic extraction"
    )
    parser.add_argument(
        "--umi_tsv", dest="UMI_TSV", nargs="+", required=True, type=str, help="TSV reports of UMI extraction"
    )
    parser.add_argument(
        "-o", "--output", dest="OUTPUT", required=True, help="Output folder"
    )
    parser.add_argument(
        "--cons", dest="CONS", action="store_true", help="Mark consensus round of extraction"
    )
    parser.add_argument(
        "--tsv", dest="TSV", action="store_true", help="Write stats tsv"
    )

    args = parser.parse_args(argv)

    return args


def merge_detected(files, out_dir):
    output_file = os.path.join(out_dir, "detected_umi_stats.tsv")

    with open(output_file, 'w') as out:
        print(
            "output_file",
            "max_pattern_distance",
            "detected_forward_strands",
            "detected_reverse_strands",
            "ratio_fwd_rvs",
            "total_reads",
            "included_reads",
            "percent_of_total_reads",
            sep="\t",
            file=out
        )

        name = "detected_umi_stats.tsv"
        distance = 0
        fwd = 0
        rev = 0
        ratio = 0
        total = 0
        included = 0
        perc_total = 0
        error = 0

        for i, file in enumerate(files):
            with open(file, 'r') as tsv_file:
                header = tsv_file.readline()
                values = tsv_file.readline().split('\t')
                n_error = int(tsv_file.readline().split("'E': ")[1].split('}')[0])

            if i == 0:
                distance = values[1]
            
            fwd += int(values[2])
            rev += int(values[3])
            
            total += int(values[5])
            included += int(values[6])

            error += n_error
        
        ratio = fwd / rev
        perc_total = included / total

        print(
            name,
            distance,
            fwd,
            rev,
            ratio,
            total,
            included,
            perc_total,
            sep="\t",
            file=out
        )
        print(f'Strand detection error: {error}', file=out)


def merge_extraction(files, out_dir, type):
    output_file = os.path.join(out_dir, f'extraction_{type}_stats.tsv')

    with open(output_file, 'w') as out_f:
        header_written = False
        
        for file in files:
            with open(file, 'r') as in_f:
                lines = in_f.readlines()
                if not lines:
                    continue
                
                if not header_written:
                    out_f.write(lines[0])  # Write the header from the first file
                    header_written = True
                
                out_f.writelines(lines[1:])


def main(argv=sys.argv[1:]):
    """
    Basic command line interface.
    """
    args = parse_args(argv=argv)

    if args.TSV:
        merge_detected(args.DETECTED_TSV, args.OUTPUT)
        merge_extraction(args.SYNTHETIC_TSV, args.OUTPUT, "synthetic")
        merge_extraction(args.UMI_TSV, args.OUTPUT, "umi")


if __name__ == "__main__":
    main()