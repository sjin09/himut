#!/usr/bin/env python3

import sys
import pyfastx
import argparse
import statistics


def parse_args(args):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="FASTQ file to read",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="file to write",
    )
    args = args[1:]
    return parser.parse_args(args)


def fastq2stat(
    seq_file: str,
    out_file: str
): 

    qlen_lst = [len(seq.seq) for seq in pyfastx.Fastq(seq_file)]
    base_sum = sum(qlen_lst) 
    read_count = len(qlen_lst)
    qsmall = min(qlen_lst)
    qbig = max(qlen_lst)
    qlen_std = statistics.stdev(qlen_lst)
    qlen_mean = base_sum/float(read_count)

    o = open(out_file, "w")
    o.write("read_count: {}\n".format(read_count))
    o.write("base_sum: {}\n".format(base_sum))
    o.write("smallest_length: {}\n".format(qsmall))
    o.write("biggest_length: {}\n".format(qbig))
    o.write("mean_length: {}\n".format(qlen_mean))
    o.write("read_length_std: {}\n".format(qlen_std))
    o.close()


def main():
    options = parse_args(sys.argv)
    fastq2stat(
        options.input, 
        options.output
    )
    sys.exit(0)


if __name__ == "__main__":
    main()
