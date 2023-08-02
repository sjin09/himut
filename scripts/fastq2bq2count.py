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
        help="FASTA file to read",
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


def fastq2bq2count(
    seq_file: str,
    out_file: str
): 

    bq2count = {i: 0 for i in range(1,94)}
    for seq in pyfastx.Fastq(seq_file):
        for bq_str in seq.qual:
            bq_int = ord(bq_str) - 33
            bq2count[bq_int] += 1

    o = open(out_file, "w")
    o.write("{}\t{}\n".format("bq", "count"))
    for bq in range(1, 94):
        o.write("{}\t{}\n".format(bq, bq2count[bq]))
    o.close()
  


def main():
    options = parse_args(sys.argv)
    fastq2bq2count(
        options.input, 
        options.output
    )
    sys.exit(0)


if __name__ == "__main__":
    main()
