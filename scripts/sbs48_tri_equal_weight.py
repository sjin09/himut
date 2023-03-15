
import sys
import argparse
from collections import defaultdict


def parse_args(args):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--som",
        type=str,
        required=True,
        help="file to read somatic mutations normcounts"
    )
    parser.add_argument(
        "--germ",
        type=str,
        required=True,
        help="file to read germline mutations SBS48 counts"
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


def load_tri_counts(som_file):

    tri2counts = defaultdict(lambda: 0)
    for line in open(som_file).readlines():
        if line.startswith("#"):
            continue
        elif line.startswith("sub"):
            continue
        _, tri, _, _, _, _, _, ref_tri_count, _,  _ = line.strip().split()
        upstream, ref, downstream = list(tri)
        if ref == "T":
            tri = "{}C{}".format(upstream,downstream)
        tri2counts[tri] = int(ref_tri_count)
    return tri2counts


def dump_sbs48_tri_equal_weights(
    som_file: str,
    germ_file: str,
    out_file: str,
) -> None:


    tri2freq = {}
    tri_weight = 1/16
    tri2counts = load_tri_counts(som_file)
    tri_sum = sum(tri2counts.values())
    for tri, tri_count in tri2counts.items():
        tri2freq[tri] = tri_count/float(tri_sum)
    
    o = open(out_file, "w")
    o.write("sub\ttri\tsbs96\tcounts\tnormcounts\tref_tri_ratio\tref_ccs_tri_ratio\tref_tri_count\tref_callable_tri_count\tccs_callable_tri_count\n")
    for line in open(germ_file).readlines():
        if line.startswith("#"):
            continue
        elif line.startswith("sub"):
            continue
        sub, tri, sbs48, _, normcounts, _, _, _, _, _ = line.strip().split()
        ref_tri_count = tri2counts[tri]
        tri_weight_normaliser = tri2freq[tri]/tri_weight
        tri_weight_normcounts = float(normcounts)/tri_weight_normaliser
        o.write("{}\t{}\t{}\t.\t{}\t.\t.\t{}\t.\t.\n".format(sub, tri, sbs48, tri_weight_normcounts, ref_tri_count))
    o.close()


def main():
    options = parse_args(sys.argv)
    dump_sbs48_tri_equal_weights(options.som, options.germ, options.output)
    sys.exit(0)


if __name__ == "__main__":
    main()

