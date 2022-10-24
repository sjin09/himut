#!/usr/bin/env python

from re import U
import sys
import math
import argparse
from collections import defaultdict 


def parse_args(args):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="hdp matrix"
    )
    parser.add_argument(
        "--sample",
        type=str,
        required=True,
        help="sample id"
    )
    parser.add_argument(
        "--counts",
        type=str,
        required=True,
        help="table with SBS96 normcounts"
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="updated hdp SBS96 matrix"
    )
    args = args[1:]
    return parser.parse_args(args)


def update_sbs96_matrix(infile, cntfile, nsample, outfile):

    sample_lst = []
    sample_sbs96_count_hsh = defaultdict(str)
    for line in open(infile):
        if line.startswith("C>A"):
            sbs_lst = line.strip().split()
        else:
            arr = line.strip().split()
            sample = arr[0]
            sample_lst.append(sample)
            sample_count_lst = arr[1:]
            for i, counts in enumerate(sample_count_lst):
                sbs = sbs_lst[i] 
                sample_sbs96_count_hsh["{}:{}".format(sample, sbs)] = counts
    sample_lst.append(nsample)
    
    for line in open(cntfile):
        if line.startswith("sub"):
            continue
        arr = line.strip().split()
        sub, tri, _, normcounts = arr[:4] 
        upstream, _, downstream = list(tri)
        sbs = "{},{}-{}".format(sub, upstream, downstream)
        sample_sbs96_count_hsh["{}:{}".format(nsample, sbs)] = normcounts

    o = open(outfile, "w")
    o.write("{}\n".format("\t".join(sbs_lst))) 
    for sample in sample_lst:
        row_lst = [sample] + [str(math.ceil(float(sample_sbs96_count_hsh["{}:{}".format(sample,sbs)]))) for sbs in sbs_lst]
        o.write("{}\n".format("\t".join(row_lst)))
    o.close()


def main():
    options = parse_args(sys.argv)
    update_sbs96_matrix(
        options.input, 
        options.counts, 
        options.sample,
        options.output)
    sys.exit(0)


if __name__ == "__main__":
    main()

