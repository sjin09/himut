#!/usr/bin/env python

import argparse
import sys
from pathlib import Path


def parse_args(args):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--hdp-matrix",
        type=Path,
        required=True,
        help="file to read matrix for HDP mutational signature extraction"
    )
    parser.add_argument(
        "--rtol-filtered-samples",
        type=Path,
        required=True,
        help="file to read HDP mutational signatures"
    )
    parser.add_argument(
        "--rtol-filtered-hdp-matrix",
        type=Path,
        required=False,
        help="file to write rTOL filtered HDP matrix for mutational signature extraction"
    )
    args = args[1:]
    return parser.parse_args(args)


def load_samples(rtol_filtered_sample_path):
    samples = [line.rstrip().split()[0] for line in open(rtol_filtered_sample_path).readlines()]
    return set(samples)


def write_rTOL_filtered_hdp_matrix(
    hdp_matrix_path: Path,
    rtol_filtered_sample_path: Path,
    rtol_filtered_hdp_matrix_path: Path
):
    samples = load_samples(rtol_filtered_sample_path)
    with open(rtol_filtered_hdp_matrix_path, "w") as outfile:
        for line in open(hdp_matrix_path):
            if line.startswith("C>A"):
                outfile.write(line)
            else:
                fields = line.rstrip().split()
                sample = fields[0]
                if sample in samples:
                    outfile.write(line)


def main():
    options = parse_args(sys.argv)
    write_rTOL_filtered_hdp_matrix(options.hdp_matrix, options.rtol_filtered_samples, options.rtol_filtered_hdp_matrix)
    sys.exit(0)


if __name__ == "__main__":
    main()
