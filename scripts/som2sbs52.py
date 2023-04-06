#!/usr/bin/env python

import sys
import natsort
import argparse
from collections import defaultdict 
from himut.mutlib import sbs52_lst, sbs96_to_sbs52, sbs96_to_tri, sbs96_to_sub


def parse_args(args):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="file to read SBS96 counts"
    )
    parser.add_argument(
        "--tri",
        type=str,
        required=True,
        help="file to return SBS52 reference trinucleotide sequence contexts"
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="file to return SBS52 counts"
    )
    args = args[1:]
    return parser.parse_args(args)


def load_sbs52_counts(tsv_file: str):

    sbs52_counts = defaultdict(lambda: 0)
    sbs52_tri_counts = defaultdict(lambda: 0)
    for line in open(tsv_file).readlines():
        if line.startswith("#"):
            continue
        elif line.startswith("sub"):
            continue
        _, _, sbs96, _, normcounts, _, _, ref_tri_count, _, _ = line.strip().split()
        sbs52 = sbs96_to_sbs52[sbs96]
        tri = sbs96_to_tri[sbs52]
        sbs52_counts[sbs52] += float(normcounts)
        sbs52_tri_counts[tri] += int(ref_tri_count)
    return sbs52_counts, sbs52_tri_counts


def dump_sbs52_counts(infile, trifile, outfile):

    o = open(outfile, "w")
    sbs52_counts, sbs52_tri_counts = load_sbs52_counts(infile)
    o.write("sub\ttri\tsbs48\tcounts\tnormcounts\tref_tri_ratio\tref_ccs_tri_ratio\tref_tri_count\tref_callable_tri_count\tccs_callable_tri_count\n")
    for sbs52 in sbs52_lst:
        tri = sbs96_to_tri[sbs52]
        sub = sbs96_to_sub[sbs52]
        o.write("{}\t{}\t{}\t.\t{}\t.\t.\t.\t.\t.\n".format(sub, tri, sbs52, sbs52_counts[sbs52])) 
    o.close()

    p = open(trifile, "w")
    sbs52_tri_lst = natsort.natsorted(list(sbs52_tri_counts.keys()))
    for tri in sbs52_tri_lst:
        p.write("{}\t{}\n".format(tri, sbs52_tri_counts[tri]))
    p.close()
        
        
def main():
    options = parse_args(sys.argv)
    dump_sbs52_counts(options.input, options.tri, options.output)
    sys.exit(0)


if __name__ == "__main__":
    main()
