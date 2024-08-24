#!/usr/bin/env python

import argparse
import os
import sys
import argparse
from collections import defaultdict
from pathlib import Path
from typing import List

import pysam
import natsort
from Bio import SeqIO

NTS = ["A", "C", "G", "T"]
NTS_SET = set(NTS)


def parse_args(args):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="FASTA file to read"
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


def get_assembly_n50(seq_lens: List[int]):
    cumsum = 0
    seq_total = sum(len_lst)
    for i in len_lst:
        cumsum += i
        if cumsum >= seq_total / 2:
            return i


# def gap_length():
# def get_number_of_gaps_bases()
# def get_total_bases()
# def get_number_of_contigs()
# def get_number_of_scaffolds()
# def get_number_of_chromosomes()
 

def write_assembly_statistics(seq_path: Path, out_path: Path):
    seqfile = pysam.FastaFile(seq_path)
    # len_lst = [len(seq.seq) for seq in seqfile]
    # n50 = get_n50(sorted(len_lst, reverse=True))
    # seq_cnt = len(len_lst)
    # seq_total = sum(len_lst)
    # seq_min = min(len_lst)
    # seq_max = max(len_lst)
    # seq_mean = seq_total / seq_cnt

    # write
    with open(out_path, "w") as outfile:
        outfile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
            "SAMPLE",
            "NUMBER_OF_CHROMOSOMES",
            "NUMBER_OF_SCAFFOLDS",
            "NUMBER_OF_CONTIGS",
            "NUMBER_OF_GAPS",
            "NUMBER_OF_GAP_BASES",
            "TOTAL_BASES"
        ))
        # o.write("number_of_sequences: {}\n".format(seq_cnt))
        # o.write("N50: {}\n".format(n50))
        # o.write("min: {}\n".format(seq_min))
        # o.write("mean: {}\n".format(seq_mean))
        # o.write("max: {}\n".format(seq_max))
        # o.write("total (bp): {}\n\n".format(seq_total))


def main():
    options = parse_args(sys.argv)
    write_assembly_statistics(options.input, options.output)
    return 0


if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
