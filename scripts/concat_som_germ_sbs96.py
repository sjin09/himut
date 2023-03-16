#!/usr/bin/env python

import sys
import argparse
from collections import defaultdict 
from himut.mutlib import sbs96_lst, sbs96_to_tri, sbs96_to_sub



def parse_args(args):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--somatic",
        type=str,
        required=True,
        help="file to read somatic SBS96 counts"
    )
    parser.add_argument(
        "--somatic_error_normalised",
        type=str,
        required=True,
        help="file to read somatic SBS96 counts"
    )
    parser.add_argument(
        "--somatic_error_normalised_tri_equal_weight",
        type=str,
        required=True,
        help="file to read somatic SBS96 counts"
    )
    parser.add_argument(
        "--germline",
        type=str,
        required=True,
        help="file to read germline SBS96 counts"
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="file to return SBS48 counts"
    )
    args = args[1:]
    return parser.parse_args(args)


def load_sbs96_counts(tsv_file: str):

    sbs96_counts = defaultdict(lambda: 0)
    for line in open(tsv_file).readlines():
        if line.startswith("#"):
            continue
        elif line.startswith("sub"):
            continue
        _, _, sbs96, _, normcounts, _, _, _, _, _ = line.strip().split()
        sbs96_counts[sbs96] += float(normcounts)
    return sbs96_counts


def load_germline_sbs96_counts(tsv_file: str):

    sbs96_counts = defaultdict(lambda: 0)
    for line in open(tsv_file).readlines():
        if line.startswith("#"):
            continue
        elif line.startswith("sub"):
            continue
        _, _, sbs96, counts = line.strip().split()
        sbs96_counts[sbs96] += int(counts)
    return sbs96_counts


def dump_sbs48_counts(som_file, som_error_file, som_error_tri_file, germ_file, outfile):

    o = open(outfile, "w")
    o.write("sub\ttri\tsbs48\tcounts\tsource\n")
    som_sbs96_counts = load_sbs96_counts(som_file)
    som_error_sbs96_counts = load_sbs96_counts(som_error_file)
    som_error_tri_sbs96_counts = load_sbs96_counts(som_error_tri_file)
    germ_sbs96_counts = load_germline_sbs96_counts(germ_file)
    for sbs96 in sbs96_lst:
        tri = sbs96_to_tri[sbs96]
        sub = sbs96_to_sub[sbs96]
        o.write("{}\t{}\t{}\t{:.0f}\tGermline\n".format(sub, tri, sbs96, germ_sbs96_counts[sbs96], "Germline")) 
    for sbs96 in sbs96_lst:
        tri = sbs96_to_tri[sbs96]
        sub = sbs96_to_sub[sbs96]
        o.write("{}\t{}\t{}\t{}\tSomatic:0\n".format(sub, tri, sbs96, som_sbs96_counts[sbs96], "Somatic:0")) 
    for sbs96 in sbs96_lst:
        tri = sbs96_to_tri[sbs96]
        sub = sbs96_to_sub[sbs96]
        o.write("{}\t{}\t{}\t{}\tSomatic:1\n".format(sub, tri, sbs96, som_error_sbs96_counts[sbs96], "Somatic:1")) 
    for sbs96 in sbs96_lst:
        tri = sbs96_to_tri[sbs96]
        sub = sbs96_to_sub[sbs96]
        o.write("{}\t{}\t{}\t{}\tSomatic:2\n".format(sub, tri, sbs96, som_error_tri_sbs96_counts[sbs96], "Somatic:2")) 

    o.close()

        
def main():
    options = parse_args(sys.argv)
    dump_sbs48_counts(options.somatic, options.somatic_error_normalised, options.somatic_error_normalised_tri_equal_weight,  options.germline, options.output)
    sys.exit(0)

if __name__ == "__main__":
    main()
