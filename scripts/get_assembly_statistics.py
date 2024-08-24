#!/usr/bin/env python

import argparse
import sys
from pathlib import Path
from typing import List

import pysam

NTS = ["A", "C", "G", "T"]
NTS_SET = set(NTS)
SEX_CHROM_SET = set(["U", "W", "X", "Y", "Z"])  # X,Y in mammals, # "U" in plants # W, Z in insects


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
        "--sample",
        type=str,
        required=True,
        help="sample"
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


def get_chromosomes(sequence_lookup: pysam.FastaFile) -> List[str]:
    chroms = []
    for seqid in sequence_lookup.references:
        if seqid.isdigit():  # autosomes
            chroms.append(seqid)
        else:  # sex chromosomes
            if seqid in SEX_CHROM_SET:
                chroms.append(seqid)
    return chroms


def get_total_bases(sequence_lookup: pysam.FastaFile, chroms: List[str]) -> int:
    total_bases = sum([len(sequence_lookup.fetch(chrom)) for chrom in chroms])
    return total_bases


def get_gap_count(sequence_lookup: pysam.FastaFile, chroms: List[str]) -> int:
    def count_gaps(sequence: str) -> int:
        gap_count = 0
        in_gap = False
        for base in sequence:
            if base == "N":
                if not in_gap:  # new gap
                    gap_count += 1
                    in_gap = True
            else:
                in_gap = False
        return gap_count

    total_gap_count = 0
    for chrom in chroms:
        seq = sequence_lookup.fetch(chrom)
        total_gap_count += count_gaps(seq)
    return total_gap_count


def get_number_of_gap_bases(sequence_lookup: pysam.FastaFile, chroms: List[str]):
    gap_bases = 0
    for chrom in chroms:
        seq = sequence_lookup.fetch(chrom)
        gap_bases += seq.count("N")
    return gap_bases


def write_assembly_statistics(seq_path: Path, sample: str, out_path: Path):
    sequence_lookup = pysam.FastaFile(seq_path)
    chroms = get_chromosomes(sequence_lookup)
    chrom_count = len(chroms)
    total_bases = get_total_bases(sequence_lookup, chroms)
    total_gap_bases = get_number_of_gap_bases(sequence_lookup, chroms)
    total_gap_count = get_gap_count(sequence_lookup, chroms)
    if chrom_count == 0:
        raise ValueError("No chromosomes found in {} FASTA file".format(seq_path))
    with open(out_path, "w") as outfile:
        outfile.write("{}\t{}\t{}\t{}\t{}\n".format(
            "SAMPLE",
            "NUMBER_OF_CHROMOSOMES",
            "CHROMOSOME_BASES",
            "NUMBER_OF_CHROMOSOME_GAPS",
            "NUMBER_OF_GAP_CHROMOSOME_BASES",
        ))
        outfile.write("{}\t{}\t{}\t{}\t{}\n".format(
            sample,
            chrom_count,
            total_bases,
            total_gap_count,
            total_gap_bases
        ))


def main():
    options = parse_args(sys.argv)
    write_assembly_statistics(options.input, options.sample, options.output)
    return 0


if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
