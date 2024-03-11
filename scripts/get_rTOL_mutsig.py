#!/usr/bin/env python

import argparse
import itertools
import sys
import natsort
from pathlib import Path
from collections import defaultdict
from typing import Set, Dict, List, Tuple

from numpy import dot
from numpy.linalg import norm


PUR_SET = set(["A", "G"])
NTS = ["A", "C", "G", "T"]
SUB_LST = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
COMPLEMENTARY_BASE_LOOKUP = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}

SBS96_TRI_SET = set()
SORTED_SBS96_LST = []
SUB_TO_SBS96_LST = {sub: [] for sub in SUB_LST}
for sub in SUB_LST:
    ref, alt = sub.split(">")
    for nti in NTS:
        for ntj in NTS:
            tri = f"{nti}{ref}{ntj}"
            sbs96 = f"{nti}[{sub}]{ntj}"
            SBS96_TRI_SET.add(tri)
            SUB_TO_SBS96_LST[sub].append(sbs96)
    SUB_TO_SBS96_LST[sub] = natsort.natsorted(SUB_TO_SBS96_LST[sub])
    SORTED_SBS96_LST.extend(SUB_TO_SBS96_LST[sub])
SORTED_SBS96_LST.append("N[N>N]N")
SBS96_TRI_LST = natsort.natsorted(list(SBS96_TRI_SET))
SBS96_TRI_COUNT = len(SBS96_TRI_LST)
SBS96_TRI_WEIGHT = 1/float(SBS96_TRI_COUNT)
SBS96_MUTSIG_FILL_COLOURS = ("#98D7EC", "#212121", "#FF003A", "#A6A6A6", "#83A603", "#F5ABCC")


def parse_args(args):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--hdp_matrix",
        type=Path,
        required=True,
        help="file to read matrix for HDP mutational signatures"
    )
    parser.add_argument(
        "--hdp_mutsig",
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


def get_rTOL_filtered_samples(mutsig_exposure_path: Path, artefactual_mutsig_path: Path):
    samples = set()
    filtered_samples = set()
    rTOL_mutsigs = [line.strip() for line in open(artefactual_mutsig_path).readlines()]
    header = open(mutsig_exposure_path).readline().rstrip().split(",")
    mutsigs = ["HDP{}".format(i) for i in header]
    for line in open(mutsig_exposure_path).readlines()[1:]:
        rTOL_sum = 0
        fields = line.rstrip().split(",")
        for rTOL in rTOL_mutsigs:
            rTOL_sum += float(fields[mutsigs.index(rTOL)])
        sample = fields[0]
        print(sample, rTOL_sum)
        if rTOL_sum > 0.5:
            filtered_samples.add(sample)
            continue
        samples.add(sample)
    return samples


def write_rTOL_filtered_matrix(
    matrix_infile_path: Path,
    mutsig_exposure_path: Path,
    artefactual_mutsig_path: Path,
    matrix_outfile_path: Path,
):
    samples = get_rTOL_filtered_samples(mutsig_exposure_path, artefactual_mutsig_path)
    with open(matrix_outfile_path, "w") as outfile:
        for line in open(matrix_infile_path).readlines():
            if line.startswith("C>A"):
                outfile.write(line)
            else:
                fields = line.rstrip().split()
                sample = fields[0]
                if sample in samples:
                    outfile.write(line)


def get_cosine_similarity(va, vb) -> float:
    cosine_similarity = dot(va, vb)/(norm(va)*norm(vb))
    return cosine_similarity


# def load_sbs96_cnt(infile):

#     sbs2cnt = {}
#     for line in open(infile).readlines():
#         if line.startswith("sub"):
#             continue
#         arr = line.strip().split()
#         sbs2cnt[arr[2]] = arr[3]
#     cnt_lst = [float(sbs2cnt[sbs]) for sbs in sbs_lst]
#     return cnt_lst


# def dump_cos_sim(target, query):
#     """
#     measures and returns cosine similarity
#     """
    
#     query_sbs96_cnt = load_sbs96_cnt(query)
#     target_sbs96_cnt = load_sbs96_cnt(target)
#     cos_sim = get_cos_sim(target_sbs96_cnt, query_sbs96_cnt)
#     print("target:{}\tquery:{}\tsim:{}".format(target, query, cos_sim))


def write_rTOL_mutsig(
    hdp_matrix_file_path: Path,
    hdp_mutsig_file_path: Path,
    cord_blood_mutspectrum_file_path: Path,
    dtol_mutspectrum_file_path: Path,
    rtol_mutsig_file_path: Path 
):
    


def main():
    options = parse_args(sys.argv)
    write_rTOL_mutsig(options.mutsig, options.cord_blood, options.library_error, options.output)
    sys.exit(0)


if __name__ == "__main__":
    main()

