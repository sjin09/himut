import sys
import natsort
import argparse
import pandas as pd
from plotnine import *
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
        help="file to read SBS96 normalised counts"
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="file to return SBS96 normalised counts with equal trinucleotide weight"
    )
    args = args[1:]
    return parser.parse_args(args)


NTS = list("ATGC")
PYR_LST = list("CT")
TRI_LST = [f"{nt1}{pyr}{nt2}" for nt1 in NTS for pyr in PYR_LST for nt2 in NTS]
TRI_LST = natsort.natsorted(TRI_LST)
TRI_COUNT = len(TRI_LST)
TRI_WEIGHT = 1/TRI_COUNT


def load_tri_freq(sbs96_file):

    tri_count = defaultdict(lambda: 0)
    for line in open(sbs96_file):
        if line.startswith("#"):
            continue
        elif line.startswith("sub"):
            continue
        (
            _sub,
            tri,
            _sbs96,
            _count,
            _normcount,
            _ref_tri_ratio,
            _ref_ccs_tri_ratio,
            ref_tri_count,
            _ref_callable_tri_count,
            _ccs_callable_tri_count,
        ) = line.strip().split()
        tri_count[tri] = int(ref_tri_count)
    tri_sum = sum(tri_count.values())
    tri_freq = {tri: tri_count/float(tri_sum) for tri, tri_count in tri_count.items()}
    return tri_freq



def dump_sbs96_plt(
    infile: str, 
    outfile: str
):

    df = pd.read_csv(infile, sep="\t")
    plot = (
        ggplot(df, aes(x="TRI", y="WEIGHTED_COUNT", fill="SUB"))
        + geom_bar(stat="identity")
        + theme_bw()
        + facet_grid(". ~ SUB", scales="free")
        + scale_fill_manual(
            values=("#98D7EC", "#212121", "#FF003A", "#A6A6A6", "#83A603", "#F5ABCC")
        )
        + labs(x="\nTrinucleotide Context\n", y="\nCounts\n")
        + theme(
            text=element_text(size=10),
            legend_title=element_blank(),
            axis_text_x=element_text(
                family="monospace", angle=90, ha="center"
            ),
        )
    )
    plot.save(outfile, width=22, height=12)


def trinucleotide_equal_weights_normcounts(
    sbs96_file_path: str,
    weighted_sbs96_file_path: str,
) -> None:


    tri_freq = load_tri_freq(sbs96_file_path)
    with open(weighted_sbs96_file_path, "w") as outfile:
        print(
            "SUB",
            "TRI",
            "SBS96",
            "WEIGHTED_COUNT",
            "COUNT",
            "WEIGHT",
            sep="\t",
            file=outfile
        )
        for line in open(sbs96_file_path).readlines():
            if line.startswith("#"):
                continue
            elif line.startswith("sub"):
                continue
            (
                sub,
                tri,
                sbs96,
                _count,
                normcount,
                _ref_tri_ratio,
                _ref_ccs_tri_ratio,
                _ref_tri_count,
                _ref_callable_tri_count,
                _ccs_callable_tri_count,
            ) = line.strip().split()
            tri_equal_weight = tri_freq[tri]/float(TRI_WEIGHT)
            tri_equal_weight_normcount = float(normcount)/tri_equal_weight
            print(
                sub,
                tri,
                sbs96,
                tri_equal_weight_normcount,
                normcount,
                tri_equal_weight,
                sep="\t",
                file=outfile
            )

    if weighted_sbs96_file_path.endswith(".tsv"):
        weighted_sbs96_pdf_file_path = weighted_sbs96_file_path.replace(".tsv", ".pdf")
    else:
        weighted_sbs96_pdf_file_path = weighted_sbs96_file_path + ".pdf"
    dump_sbs96_plt(weighted_sbs96_file_path, weighted_sbs96_pdf_file_path)

def main():
    options = parse_args(sys.argv)
    trinucleotide_equal_weights_normcounts(options.input, options.output)
    sys.exit(0)


if __name__ == "__main__":
    main()

