#!/usr/bin/env python

import sys
import argparse
from collections import defaultdict 
from himut.mutlib import sbs48_lst, sbs96_to_sbs48, sbs96_to_tri, sbs96_to_sub



def parse_args(args):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--somatic",
        type=str,
        required=True,
        help="file to read somatic SBS48/SBS96 counts"
    )
    parser.add_argument(
        "--germline",
        type=str,
        required=True,
        help="file to read germline SBS48/SBS96 counts"
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


def load_sbs48_counts(tsv_file: str):

    sbs48_counts = defaultdict(lambda: 0)
    for line in open(tsv_file).readlines():
        if line.startswith("#"):
            continue
        elif line.startswith("sub"):
            continue
        _, _, sbs96, _, normcounts, _, _, _, _, _ = line.strip().split()
        sbs48 = sbs96_to_sbs48[sbs96]
        sbs48_counts[sbs48] += float(normcounts)
    return sbs48_counts


def dump_sbs48_counts(som_file, germ_file, outfile):

    o = open(outfile, "w")
    o.write("sub\ttri\tsbs48\tcounts\tsource\n")
    som_sbs48_counts = load_sbs48_counts(som_file)
    germ_sbs48_counts = load_sbs48_counts(germ_file)
    for sbs48 in sbs48_lst:
        tri = sbs96_to_tri[sbs48]
        sub = sbs96_to_sub[sbs48]
        o.write("{}\t{}\t{}\t{}\tSomatic\n".format(sub, tri, sbs48, som_sbs48_counts[sbs48], "Somatic")) 
    for sbs48 in sbs48_lst:
        tri = sbs96_to_tri[sbs48]
        sub = sbs96_to_sub[sbs48]
        o.write("{}\t{}\t{}\t{:.0f}\tGermline\n".format(sub, tri, sbs48, germ_sbs48_counts[sbs48], "Germline")) 
    o.close()

        
def main():
    options = parse_args(sys.argv)
    dump_sbs48_counts(options.somatic, options.germline, options.output)
    sys.exit(0)


if __name__ == "__main__":
    main()
