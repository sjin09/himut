#!/usr/bin/env python3

import sys
import pyfastx
import argparse


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


def fastq2bq2count(
    seq_file: str,
    out_file: str
): 

    o = open(out_file, "w")
    o.write("{}\t{}\n".format("qname", "qlen"))
    for seq in pyfastx.Fastq(seq_file):
        o.write("{}\t{}\n".format(seq.name, len(seq.seq)))
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
