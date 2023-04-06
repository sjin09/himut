
import sys
import argparse
from collections import defaultdict
from himut.mutlib import sbs52_lst, sbs96_to_tri

def parse_args(args):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="file to read SBS952 counts"
    )
    parser.add_argument(
        "--tri",
        type=str,
        required=True,
        help="file to target reference genome trinucleotide sequence context counts"
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="file to return SBS52 with equal trinucleotide sequence context weights"
    )
    args = args[1:]
    return parser.parse_args(args)


def load_tri_counts(trifile):

    tri2counts = defaultdict(lambda: 0)
    for line in open(trifile).readlines():
        tri, count = line.split()
        tri2counts[tri] = int(count)
    tri_sum = sum(tri2counts.values())
    return tri_sum, tri2counts


def dump_sbs52_tri_equal_weights(
    infile: str,
    trifile: str,
    outfile: str,
) -> None:

    tri_weight = 1/26
    tri_sum, tri2counts = load_tri_counts(trifile)
    tri2freq = {tri: tri_count/float(tri_sum) for tri, tri_count in tri2counts.items()}
    
    o = open(outfile, "w")
    o.write("sub\ttri\tsbs96\tcounts\tnormcounts\tref_tri_ratio\tref_ccs_tri_ratio\tref_tri_count\tref_callable_tri_count\tccs_callable_tri_count\n")
    for line in open(infile).readlines():
        if line.startswith("sub"):
            continue
        sub, tri, sbs52, _, normcounts, _, _, _, _, _ = line.strip().split()
        tri_count = tri2counts[tri]
        tri_equal_weight = tri2freq[tri]/tri_weight
        tri_equal_weight_normcounts = float(normcounts)/tri_equal_weight
        o.write("{}\t{}\t{}\t.\t{}\t.\t.\t{}\t.\t.\n".format(sub, tri, sbs52, tri_equal_weight_normcounts, tri_count))
    o.close()


def main():
    options = parse_args(sys.argv)
    dump_sbs52_tri_equal_weights(options.input, options.tri, options.output)
    sys.exit(0)


if __name__ == "__main__":
    main()

