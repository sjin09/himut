
#!/usr/bin/env python

import argparse
import itertools
import pandas as pd
import sys

from plotnine import *


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
        help="target chromosome trinucleotide sequence context counts"
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="file to return trinucleotide sequence context weighted SBS96 counts"
    )
    args = args[1:]
    return parser.parse_args(args)


NTS = ["A", "T", "C", "G"]
PURINE = set(["A", "G"])
PYRIMIDINE = set(["C", "T"])
PURINE2PYRIMIDINE = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
SBS96_TRI_LST = ["".join(nt) for nt in itertools.product(NTS, PYRIMIDINE, NTS)]
SBS96_TRI_WEIGHT = 1/len(SBS96_TRI_LST)
SBS96_LST = [
    "A[C>A]A",
    "A[C>A]C",
    "A[C>A]G",
    "A[C>A]T",
    "A[C>G]A",
    "A[C>G]C",
    "A[C>G]G",
    "A[C>G]T",
    "A[C>T]A",
    "A[C>T]C",
    "A[C>T]G",
    "A[C>T]T",
    "A[T>A]A",
    "A[T>A]C",
    "A[T>A]G",
    "A[T>A]T",
    "A[T>C]A",
    "A[T>C]C",
    "A[T>C]G",
    "A[T>C]T",
    "A[T>G]A",
    "A[T>G]C",
    "A[T>G]G",
    "A[T>G]T",
    "C[C>A]A",
    "C[C>A]C",
    "C[C>A]G",
    "C[C>A]T",
    "C[C>G]A",
    "C[C>G]C",
    "C[C>G]G",
    "C[C>G]T",
    "C[C>T]A",
    "C[C>T]C",
    "C[C>T]G",
    "C[C>T]T",
    "C[T>A]A",
    "C[T>A]C",
    "C[T>A]G",
    "C[T>A]T",
    "C[T>C]A",
    "C[T>C]C",
    "C[T>C]G",
    "C[T>C]T",
    "C[T>G]A",
    "C[T>G]C",
    "C[T>G]G",
    "C[T>G]T",
    "G[C>A]A",
    "G[C>A]C",
    "G[C>A]G",
    "G[C>A]T",
    "G[C>G]A",
    "G[C>G]C",
    "G[C>G]G",
    "G[C>G]T",
    "G[C>T]A",
    "G[C>T]C",
    "G[C>T]G",
    "G[C>T]T",
    "G[T>A]A",
    "G[T>A]C",
    "G[T>A]G",
    "G[T>A]T",
    "G[T>C]A",
    "G[T>C]C",
    "G[T>C]G",
    "G[T>C]T",
    "G[T>G]A",
    "G[T>G]C",
    "G[T>G]G",
    "G[T>G]T",
    "T[C>A]A",
    "T[C>A]C",
    "T[C>A]G",
    "T[C>A]T",
    "T[C>G]A",
    "T[C>G]C",
    "T[C>G]G",
    "T[C>G]T",
    "T[C>T]A",
    "T[C>T]C",
    "T[C>T]G",
    "T[C>T]T",
    "T[T>A]A",
    "T[T>A]C",
    "T[T>A]G",
    "T[T>A]T",
    "T[T>C]A",
    "T[T>C]C",
    "T[T>C]G",
    "T[T>C]T",
    "T[T>G]A",
    "T[T>G]C",
    "T[T>G]G",
    "T[T>G]T",
]


def load_sbs96_counts(
    sbs96_count_file: str
):
    sbs96_counts = {sbs96: 0 for sbs96 in SBS96_LST}
    for line in open(sbs96_count_file).readlines():
        if line.startswith("SUB"):
            continue
        _sub, _tri, sbs96, count = line.strip().split()
        sbs96_counts[sbs96] = int(count)
    return sbs96_counts 


def load_trinucleotide_frequencies(
    tri_file: str
):
    tri_counts = {tri: 0 for tri in SBS96_TRI_LST}
    for line in open(tri_file).readlines():
        tri, count = line.strip().split()
        tri_counts[tri] = int(count) 
    tri_sum = sum(tri_counts.values())
    tri_freq = {tri: count/tri_sum for tri, count in tri_counts.items()}
    return tri_freq


def draw_sbs96_barplot(
    sbs96_file: str, 
    sbs96_pdf_file: str
):

    df = pd.read_csv(sbs96_file, sep="\t")
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
    plot.save(sbs96_pdf_file, width=22, height=12)


def dump_weighted_sbs96_counts(
    sbs96_count_file: str, 
    tri_file: str,
    weighted_sbs96_count_file: str,
):

    sbs96_counts = load_sbs96_counts(sbs96_count_file)
    trinucleotide_frequencies = load_trinucleotide_frequencies(tri_file)
    with open(weighted_sbs96_count_file, "w") as outfile:
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
        for sbs96 in SBS96_LST:
            ubase, _, ref, _, alt, _, dbase = list(sbs96)
            sub = f"{ref}>{alt}"
            tri = f"{ubase}{ref}{dbase}"
            tri_equal_weight = 1/(trinucleotide_frequencies[tri]/float(SBS96_TRI_WEIGHT))
            sbs96_count = sbs96_counts[sbs96] 
            weighted_sbs96_count =  sbs96_count * tri_equal_weight
            print(
                sub,
                tri,
                sbs96,
                weighted_sbs96_count,
                sbs96_count,
                tri_equal_weight,
                sep="\t",
                file=outfile
            )
    if weighted_sbs96_count_file.endswith(".tsv"):
        weighted_sbs96_count_barplot_file = weighted_sbs96_count_file.replace(".tsv", ".pdf")
    else:
        weighted_sbs96_count_barplot_file = weighted_sbs96_count_file + ".pdf" 
    draw_sbs96_barplot(weighted_sbs96_count_file, weighted_sbs96_count_barplot_file)
       
        
def main():
    options = parse_args(sys.argv)
    dump_weighted_sbs96_counts(options.input, options.tri, options.output)
    sys.exit(0)


if __name__ == "__main__":
    main()
