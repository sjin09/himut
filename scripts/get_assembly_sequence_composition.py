#!/usr/bin/env python

import argparse
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List

import pysam


NTS = ["A", "C", "G", "T", "N"]


def parse_args(args):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        required=True,
        help="FASTA file to read"
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="file to write"
    )
    args = args[1:]
    return parser.parse_args(args)


def get_sequence_composition(seq_path: Path) -> Dict[str, int]:
    sequence_lookup = pysam.FastaFile(seq_path)
    count_per_nt = defaultdict(lambda: 0)
    for seqid in sequence_lookup.references:
        seq = sequence_lookup[seqid]
        for nt in NTS:
            count_per_nt[nt] += seq.count(nt)
    return count_per_nt


def write_assembly_sequence_composition(seq_path: Path, out_path: Path):
    count_per_nt = get_sequence_composition(seq_path)
    nts_sum = sum(count_per_nt.values())
    with open(out_path, "w") as outfile:
        outfile.write("SAMPLE\tA\tC\tG\tT\tN\n")
        cols = [seq_path.stem]
        for nt in NTS:
            nt_proportion = (count_per_nt[nt]/float(nts_sum)) * 100
            nt_proportion = "{:.2f}".format(nt_proportion)
            cols.append(str(nt_proportion))
        outfile.write("\t".join(cols))


def main():
    options = parse_args(sys.argv)
    write_assembly_sequence_composition(options.input, options.output)
    return 0


if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)

