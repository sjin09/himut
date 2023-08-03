#!/usr/bin/env python3

import sys
import pysam
import argparse
import himut.bamlib

def parse_args(args):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="BAM file to read",
    )
    parser.add_argument(
        "-n",
        "--count",
        type=int,
        default=1000,
        required=False,
        help="number of reads",
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


def bam2bq2count(
    bam_file: str,
    count_threshold: int,
    out_file: str
): 

    counter = 0
    bq2count = {i: 0 for i in range(1,94)}
    alignments = pysam.AlignmentFile(bam_file, "rb")
    for line in alignments:
        ccs = himut.bamlib.BAM(line)
        if not ccs.is_primary:
            continue
        for bq_int in ccs.bq_int_lst:
            bq2count[bq_int] += 1
        counter += 1 
        if counter > count_threshold:
            break

    o = open(out_file, "w")
    o.write("{}\t{}\n".format("bq", "count"))
    for bq_int in range(1, 94):
        bq_str = chr(bq_int + 33)
        o.write("{}:{}\t{}\n".format(bq_int, bq_str, bq2count[bq_int]))
    o.close()
  

def main():
    options = parse_args(sys.argv)
    bam2bq2count(
        options.input, 
        options.count, 
        options.output
    )
    sys.exit(0)


if __name__ == "__main__":
    main()
