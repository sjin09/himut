#!/ur/bin/env python

import re
import sys
import gzip
import natsort
import pyfastx
import argparse

def parse_args(args):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-i",
        "--fofn",
        type=str,
        required=True,
        help="FASTQ file of file names",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="FASTQ file to write",
    )
    args = args[1:]
    return parser.parse_args(args)


def cat_fastq(
    fofn_file: str,
    outfile: str,
):

    o = open(outfile, "w")
    fq_lst = natsort.natsorted([line.strip() for line in open(fofn_file)])
    for fq in fq_lst:
        for seq in pyfastx.Fastq(fq):
            o.wrseqte("@{}\n{}\n+\n{}\n".format(seq.name, seq.seq, seq.qual))
    o.close()


def main():
    options = parse_args(sys.argv)
    cat_fastq(
        options.fofn,
        options.output
    )
    sys.exit(0)


if __name__ == "__main__":
    main()

