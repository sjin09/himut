#!/usr/bin/env python

import os
import sys
import math
import gzip
import natsort
import argparse
import itertools
import statistics
from collections import defaultdict 

sub2sbs = defaultdict(list)
sub_lst = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
for sub in sub_lst:
    for upstream, downstream in itertools.product(list("ATGC"), repeat=2):
        sbs = "{},{}-{}".format(sub, upstream, downstream)
        sub2sbs[sub].append(sbs) 

sbs_lst = []
for sub in sub_lst:
    sub2sbs[sub] = natsort.natsorted(list(set(sub2sbs[sub])))
    for sbs in sub2sbs[sub]:
        sbs_lst.append(sbs)


def parse_args(args):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="normcounts.fofn"
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="SBS96 matrix"
    )

    args = args[1:]
    return parser.parse_args(args)


def get_sbs96_matrix(infile, outfile):

    sample_sbs96_count_hsh = defaultdict(str)
    sample_lst = [i.strip().split("/")[-1].split(".")[0] for i in open(infile).readlines()]
    for i,j in enumerate(open(infile).readlines()):
        sample = sample_lst[i]
        for line in open(j.strip()).readlines():
            if line.startswith("sub"):
                continue
            arr = line.strip().split()
            sub, tri, _, normcounts = arr[:4] 
            upstream, _, downstream = list(tri)
            sbs = "{},{}-{}".format(sub, upstream, downstream)
            sample_sbs96_count_hsh["{}:{}".format(sample, sbs)] = normcounts

    o = open(outfile, "w")
    o.write("{}\n".format("\t".join(sbs_lst))) 
    for sample in sample_lst:
        count_lst = [str(math.ceil(float(sample_sbs96_count_hsh["{}:{}".format(sample,sbs)]))) for sbs in sbs_lst]
        o.write("{}\n".format("\t".join([sample] + count_lst)))


def main():
    options = parse_args(sys.argv)
    get_sbs96_matrix(options.input, options.output)
    sys.exit(0)


if __name__ == "__main__":
    main()
