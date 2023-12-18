
import sys
import argparse
from collections import defaultdict


def parse_args(args):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-i",
        "--som",
        type=str,
        required=True,
        help="file to read normalised somatic mutations counts"
    )
    parser.add_argument(
        "--tri",
        type=str,
        required=True,
        help="file to reference FASTA file trinucleotide sequence contexts"
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="file to return SBS48 counts with equal trinucleotide weights"
    )
    args = args[1:]
    return parser.parse_args(args)


def load_tri_counts(tri_file):

    tri2counts = defaultdict(lambda: 0)
    for line in open(tri_file).readlines():
        tri, count = line.strip().split()
        tri2counts[tri] = int(count)
    return tri2counts


def dump_sbs96_tri_equal_weights(
    som_file: str,
    tri_file: str,
    out_file: str,
) -> None:


    tri2freq = {}
    tri_weight = 1/32
    tri2counts = load_tri_counts(tri_file)
    tri_sum = sum(tri2counts.values())
    for tri, tri_count in tri2counts.items():
        tri2freq[tri] = tri_count/float(tri_sum)
    
    o = open(out_file, "w")
    o.write("sub\ttri\tsbs96\tcounts\tnormcounts\tref_tri_ratio\tref_ccs_tri_ratio\tref_tri_count\tref_callable_tri_count\tccs_callable_tri_count\n")
    for line in open(som_file).readlines():
        if line.startswith("#"):
            continue
        elif line.startswith("sub"):
            continue
        sub, tri, sbs48, _, normcounts, _, _, _, _, _ = line.strip().split()
        tri_weight_normaliser = tri2freq[tri]/tri_weight
        tri_weight_normcounts = float(normcounts)/tri_weight_normaliser
        o.write("{}\t{}\t{}\t.\t{}\t.\t.\t.\t.\t.\n".format(sub, tri, sbs48, tri_weight_normcounts))
    o.close()


def main():
    options = parse_args(sys.argv)
    dump_sbs96_tri_equal_weights(options.som, options.tri, options.output)
    sys.exit(0)


if __name__ == "__main__":
    main()

