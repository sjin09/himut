#!/usr/bin/env python

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List

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
        help="file to read samples.fofn"
    )
    parser.add_argument(
        "--sequence-length-fofn",
        type=str,
        required=True,
        help="file of file names to read sequence length lookup from BAM file"
    )
    parser.add_argument(
        "--dtol-assemblies",
        type=str,
        required=True,
        help="file to read dtol_assemblies.json"
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


def load_sequence_length_lookup_per_sample(sequence_length_lookup_fofn_path: Path) -> Dict[str, Dict[str, str]]:
    sequence_length_lookup_per_sample = defaultdict(dict)
    for line in open(sequence_length_lookup_fofn_path).readlines():
        sample_sequence_length_lookup_path = line.rstrip()
        sample_sequence_length_lookup = line.rstrip().split("/")[-1]
        fields = sample_sequence_length_lookup.split(".")
        if len(fields) == 3:
            sample = fields[0]
        elif len(fields) == 4:
            sample = "{}.{}".format(fields[0], fields[1])
        sample_length_lookup = dict(
            [line.rstrip().split() for line in open(sample_sequence_length_lookup_path).readlines()]
        )
        sequence_length_lookup_per_sample[sample] = sample_length_lookup
    return sequence_length_lookup_per_sample


def load_fasta_sequence_lengths_per_sample(dtol_assemblies_json_path: Path) -> Dict[str, Dict[str, List[int]]]:
    released_assembly_per_sample = json.load(open(dtol_assemblies_json_path))
    assembly_lengths_per_assembly_per_sample = defaultdict(dict)
    for sample in released_assembly_per_sample:
        for released_assembly in released_assembly_per_sample[sample]:
            reference_sequence_lookup = pysam.FastaFile(released_assembly)
            reference_sequence_lengths = [
                reference_sequence_lookup.get_reference_length(seqid)for seqid in reference_sequence_lookup.references
            ]
            reference_sequence_lengths.sort()
            assembly_lengths_per_assembly_per_sample[sample][released_assembly] = reference_sequence_lengths
    return assembly_lengths_per_assembly_per_sample


def get_correct_released_assembly_per_sample(
    samples_path: Path,
    dtol_assemblies_json_path: Path,
    sequence_length_lookup_fofn_path: Path,
    out_path: Path
):
    assembly_lengths_per_assembly_per_sample = load_fasta_sequence_lengths_per_sample(dtol_assemblies_json_path)
    sequence_length_lookup_per_sample = load_sequence_length_lookup_per_sample(sequence_length_lookup_fofn_path)
    with open(out_path, "w") as outfile:
        for line in open(samples_path).readlines():
            fields = line.rstrip().split(".")
            if len(fields) == 1:
                ref_sample = fields[0]
                sample = ref_sample
            else:
                ref_sample = fields[0]
                sample = ".".join(fields)
            contig_lengths = []
            for contig, contig_length in sequence_length_lookup_per_sample[sample].items():
                contig_lengths.append(int(contig_length))
            contig_lengths.sort()
            for (assembly, assembly_lengths) in assembly_lengths_per_assembly_per_sample[ref_sample].items():
                if contig_lengths == assembly_lengths:
                    outfile.write("{}\t{}\t{}\tTrue\n".format(ref_sample, sample, assembly))
                else:
                    outfile.write("{}\t{}\t{}\tFalse\n".format(ref_sample, sample, assembly))

                


def main():
    options = parse_args(sys.argv)
    get_correct_released_assembly_per_sample(
        options.input, 
        options.dtol_assemblies, 
        options.sequence_length_fofn,
        options.output
    )
    return 0


if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
