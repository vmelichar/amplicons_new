"""
Merge filtering tsv file from different chunks
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
        "--filter_tsv", dest="FILTER_TSV", nargs="+", required=True, type=str, help="TSV reports of filtering"
    )
    parser.add_argument(
        "-o", "--output", dest="OUTPUT", required=True, help="Output folder"
    )
    parser.add_argument(
        "--tsv", dest="TSV", action="store_true", help="Write stats tsv"
    )

    args = parser.parse_args(argv)

    return args


def merge_filter(files, out_dir):
    output_file = os.path.join(out_dir, 'filtering_stats.tsv')

    with open(output_file, 'w') as out_f:
        region = ''
        found = 0
        unmapped = 0
        second = 0
        suppl = 0
        target = 0
        conca = 0
        short = 0
        long_r = 0
        filtered = 0
        incl = ''
        
        for i,file in enumerate(files):
            with open(file, 'r') as in_f:
                lines = in_f.readlines()
                if not lines:
                    continue

                values = lines[1].split('\t')
                
                if i == 0:
                    out_f.write(lines[0])  # Write the header from the first file
                    region += values[1]
                    incl += values[11].strip()
                
                found += int(values[2])
                unmapped += int(values[3])
                second += int(values[4])
                suppl += int(values[5])
                target += int(values[6])
                conca += int(values[7])
                short += int(values[8])
                long_r += int(values[9])
                filtered += int(values[10])
        
        print(
            'count',
            region,
            found,
            unmapped,
            second,
            suppl,
            target,
            conca,
            short,
            long_r,
            filtered,
            incl,
            sep='\t',
            file=out_f
        )
        print(
            '%',
            region,
            round(found / found * 100),
            round(unmapped / found * 100),
            round(second / found * 100),
            round(suppl / found * 100),
            round(target / found * 100),
            round(conca / target * 100),
            round(short / target * 100),
            round(long_r / target * 100),
            round(filtered / found * 100),
            incl,
            sep='\t',
            file=out_f
        )


def main(argv=sys.argv[1:]):
    """
    Basic command line interface.
    """
    args = parse_args(argv=argv)

    if args.TSV:
        merge_filter(args.FILTER_TSV, args.OUTPUT)


if __name__ == "__main__":
    main()