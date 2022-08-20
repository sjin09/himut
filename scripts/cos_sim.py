#!/usr/bin/env python

import os
import sys
import natsort
import argparse
import itertools
from numpy import dot
from numpy.linalg import norm
from collections import defaultdict
from typing import Set, Dict, List, Tuple


def parse_args(args):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--target",
        type=str,
        required=True,
        help="path to target SBS96 count (reference)",
    )
    parser.add_argument(
        "--query",
        type=str,
        required=True,
        help="path to query SBS96 count",
    )
    args = args[1:]
    return parser.parse_args(args)


sbs_lst = []
tri_set = set()
sbs2sub = {}
sbs2tri = {}
purine = set(["A", "G"])
pyrimidine = set(["T", "C"])
sub_lst = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
purine2pyrimidine = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
for sub in sub_lst:
    for upstream, downstream in itertools.product(list("ATGC"), repeat=2):
        ref, alt = sub.split(">")
        tri = upstream + ref + downstream
        sbs = "{}[{}]{}".format(upstream, sub, downstream)
        tri_set.add(tri)
        sbs_lst.append(sbs)
        sbs2sub[sbs] = sub
        sbs2tri[sbs] = tri
sbs_lst = natsort.natsorted(sbs_lst)
tri_lst = natsort.natsorted(list(tri_set))


def get_cos_sim(va, vb) -> float:
## def get_cos_sim(va: List[int], vb:List[int]) -> float:

    cos_sim = dot(va, vb)/(norm(va)*norm(vb))
    return cos_sim


def load_sbs96_cnt(infile):

    sbs2cnt = {}
    for line in open(infile).readlines():
        if line.startswith("sub"):
            continue
        arr = line.strip().split()
        sbs2cnt[arr[2]] = arr[3]
    cnt_lst = [float(sbs2cnt[sbs]) for sbs in sbs_lst]
    return cnt_lst


def dump_cos_sim(target, query):
    """
    measures and returns cosine similarity
    """
    
    query_sbs96_cnt = load_sbs96_cnt(query)
    target_sbs96_cnt = load_sbs96_cnt(target)
    cos_sim = get_cos_sim(target_sbs96_cnt, query_sbs96_cnt)
    print("target:{}\tquery:{}\tsim:{}".format(target, query, cos_sim))


def main():
    options = parse_args(sys.argv)
    dump_cos_sim(options.target, options.query)
    sys.exit(0)


if __name__ == "__main__":
    main()
