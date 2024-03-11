#!/usr/bin/env python

import argparse
import sys
from pathlib import Path

import natsort


def parse_args(args):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--hdp-mutsig-exposure",
        type=str,
        required=True,
        help="file to read mutational signature exposure per sample"
    )
    parser.add_argument(
        "--rTOL-mutsig",
        type=Path,
        required=False,
        help="list of rTOL mutational signatures"
    )
    parser.add_argument(
        "--rTOL-mutsig-attribution-threshold",
        type=float,
        required=False,
        help="rTOL mutational signature attribution threshold"
    )
    args = args[1:]
    return parser.parse_args(args)


def get_rTOL_filtered_samples(
    mutsig_exposure_path: Path,
    rTOL_mutsig_path: Path,
    rTOL_mutsig_attribution_threshold: int
):
    samples = []
    rTOL_mutsigs = [line.strip() for line in open(rTOL_mutsig_path).readlines()]
    header = open(mutsig_exposure_path).readline().rstrip().split(",")
    mutsigs = ["HDP{}".format(i) for i in header]
    for line in open(mutsig_exposure_path).readlines()[1:]:
        rTOL_sum = 0
        fields = line.rstrip().split(",")
        for rTOL in rTOL_mutsigs:
            rTOL_sum += float(fields[mutsigs.index(rTOL)])
        if rTOL_sum > rTOL_mutsig_attribution_threshold:
            continue
        sample = fields[0]
        samples.append((sample, rTOL_sum))
    samples = natsort.natsorted(samples)
    return samples


def write_rTOL_filtered_matrix(
    mutsig_exposure_path: Path,
    rTOL_mutsig_path: Path,
    rTOL_mutsig_attribution_threshold: float,
):
    rTOL_filtered_samples = get_rTOL_filtered_samples(
        mutsig_exposure_path,
        rTOL_mutsig_path,
        rTOL_mutsig_attribution_threshold
    )
    outfile_path = "dtol.rTOL_filtered.maximum_rTOL_attribution_{}.txt".format(rTOL_mutsig_attribution_threshold)
    with open(outfile_path, "w") as outfile:
        for (sample, rTOL_attribution) in rTOL_filtered_samples:
            print(sample, rTOL_attribution, sep="\t", file=outfile)


def main():
    options = parse_args(sys.argv)
    write_rTOL_filtered_matrix(
        options.hdp_mutsig_exposure,
        options.rTOL_mutsig,
        options.rTOL_mutsig_attribution_threshold
    )
    sys.exit(0)


if __name__ == "__main__":
    main()
