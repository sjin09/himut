#!/usr/bin/env python

import argparse
import sys
from pathlib import Path
from typing import List

import pysam


def parse_args(args):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="BAM file to read"
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="file to write"
    )
    args = args[1:]
    return parser.parse_args(args)


def get_chromosome_length_lookup(bam_path: Path) -> List[int]:
    sequence_length_lookup = {}
    alignments = pysam.AlignmentFile(bam_path, "rb")
    header = str(alignments.header).strip().split("\n")
    for line in header:
        if not line.startswith("@SQ"):
            continue
        fields = line.rstrip().split() 
        chrom = fields[1].split(":")[1]
        chrom_len = fields[2].split(":")[1]
        sequence_length_lookup[chrom] = chrom_len
    return sequence_length_lookup


def write_chromosome_lengths(bam_path: Path, out_path: Path):
    sequence_length_lookup = get_chromosome_length_lookup(bam_path)
    with open(out_path, "w") as outfile:
        for (chrom, chrom_length) in sequence_length_lookup.items():
            outfile.write("{}\t{}\n".format(
                chrom,
                chrom_length
            ))


def main():
    options = parse_args(sys.argv)
    write_chromosome_lengths(options.input, options.output)
    return 0


if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
