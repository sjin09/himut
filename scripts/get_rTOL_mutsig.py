#!/usr/bin/env python

import argparse
import sys
from pathlib import Path
from typing import Dict, List

from numpy import dot
from numpy.linalg import norm


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
        "--hdp-mutsig",
        type=Path,
        required=True,
        help="file to read HDP mutational signatures"
    )
    parser.add_argument(
        "--cord-blood",
        type=str,
        required=True,
        help="file to read cord blood mutational spectrum"
    )
    parser.add_argument(
        "--library-error",
        type=Path,
        required=False,
        help="file to read mutational spectrum from DToL sample with library errors"
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=False,
        help="file to write rTOL mutsigs"
    )
    args = args[1:]
    return parser.parse_args(args)


def load_sbs96(hdp_matrix_file_path: Path) -> List[str]:
    sbs96_lst = []
    fields = open(hdp_matrix_file_path).readline().rstrip().split()
    for i in fields:
        sub, tri = i.split(",")
        ubase, _, dbase = list(tri)
        sbs96 = "{}[{}]{}".format(ubase, sub, dbase)
        sbs96_lst.append(sbs96)
    return sbs96_lst


def load_mutsigs(hdp_mutsig_file_path: Path) -> Dict[str, List[float]]:
    header = open(hdp_mutsig_file_path).readline().rstrip().split(",")
    mutsigs = ["HDP{}".format(i.replace("X", "")) for i in header[1:]]
    mutprobs_per_mutsig = {mutsig: [] for mutsig in mutsigs}
    for line in open(hdp_mutsig_file_path).readlines()[1:]:
        fields = line.rstrip().split(",")
        for mutsig_idx, mutprob in enumerate(fields[1:]):
            mutsig = mutsigs[mutsig_idx]
            mutprobs_per_mutsig[mutsig].append(float(mutprob))
    return mutprobs_per_mutsig


def load_mutcounts(mutcount_file_path: Path, sbs96_lst: List[str]) -> List[float]:
    mutcount_per_sbs96 = {}
    for line in open(mutcount_file_path).readlines():
        if line.startswith("#"):
            continue
        elif line.startswith("sub"):
            continue
        fields = line.rstrip().split()
        sbs96 = fields[2]
        mutcount = fields[4]
        mutcount_per_sbs96[sbs96] = float(mutcount)
    mutcounts = [mutcount_per_sbs96[sbs96] for sbs96 in sbs96_lst]
    return mutcounts


def get_cosine_similarity(va, vb) -> float:
    cosine_similarity = dot(va, vb)/(norm(va)*norm(vb))
    return cosine_similarity


def get_rTOL_mutsigs(
    mutprobs_per_hdp_mutsig: Dict[str, List[float]],
    cord_blood_mutcounts: List[float],
    dtol_mutcounts: List[float]
) -> List[str]:
    rTOL_mutsigs = []
    for hdp_mutsig, mutprobs in mutprobs_per_hdp_mutsig.items():
        cord_blood_cos_sim = get_cosine_similarity(mutprobs, cord_blood_mutcounts)
        dtol_cos_sim = get_cosine_similarity(mutprobs, dtol_mutcounts)
        if cord_blood_cos_sim > 0.85:
            rTOL_mutsigs.append(hdp_mutsig)
        if dtol_cos_sim > 0.85:
            rTOL_mutsigs.append(hdp_mutsig)
    return rTOL_mutsigs


def write_rTOL_mutsig(
    hdp_matrix_file_path: Path,
    hdp_mutsig_file_path: Path,
    cord_blood_mutcount_file_path: Path,
    dtol_mutcount_file_path: Path,
    rtol_mutsig_file_path: Path
):
    sbs96_lst = load_sbs96(hdp_matrix_file_path)
    mutprobs_per_hdp_mutsig = load_mutsigs(hdp_mutsig_file_path)
    cord_blood_mutcounts = load_mutcounts(cord_blood_mutcount_file_path, sbs96_lst)
    dtol_mutcounts = load_mutcounts(dtol_mutcount_file_path, sbs96_lst)
    rTOL_mutsigs = get_rTOL_mutsigs(mutprobs_per_hdp_mutsig, cord_blood_mutcounts, dtol_mutcounts)
    with open(rtol_mutsig_file_path, "w") as outfile:
        print("\n".join(rTOL_mutsigs), file=outfile)


def main():
    options = parse_args(sys.argv)
    write_rTOL_mutsig(options.hdp_matrix, options.hdp_mutsig, options.cord_blood, options.library_error, options.output)
    sys.exit(0)


if __name__ == "__main__":
    main()
