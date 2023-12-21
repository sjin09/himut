#!/usr/bin/env python

import sys
import argparse
import natsort
import pandas as pd

import natsort
import pandas as pd

from collections import defaultdict 
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
        help="file to return trinucleotide sequence context weighted SBS52 counts"
    )
    args = args[1:]
    return parser.parse_args(args)


NTS = ["A", "T", "C", "G"]
PURINE = set(["A", "G"])
PYRIMIDINE = set(["C", "T"])
PURINE2PYRIMIDINE = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
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
SBS96_TO_SBS52 = {
    "A[C>A]A": "T[T>G]T",
    "T[T>G]T": "T[T>G]T",
    "A[C>A]C": "A[C>A]C",
    "G[T>G]T": "A[C>A]C",
    "A[C>A]G": "C[T>G]T",
    "C[T>G]T": "C[T>G]T",
    "A[C>A]T": "A[C>A]T",
    "A[T>G]T": "A[C>A]T",
    "C[C>A]A": "C[C>A]A",
    "T[T>G]G": "C[C>A]A",
    "C[C>A]C": "C[C>A]C",
    "G[T>G]G": "C[C>A]C",
    "C[C>A]G": "C[C>A]G",
    "C[T>G]G": "C[C>A]G",
    "C[C>A]T": "C[C>A]T",
    "A[T>G]G": "C[C>A]T",
    "G[C>A]A": "T[T>G]C",
    "T[T>G]C": "T[T>G]C",
    "G[C>A]C": "G[C>A]C",
    "G[T>G]C": "G[C>A]C",
    "G[C>A]G": "C[T>G]C",
    "C[T>G]C": "C[T>G]C",
    "G[C>A]T": "A[T>G]C",
    "A[T>G]C": "A[T>G]C",
    "T[C>A]A": "T[C>A]A",
    "T[T>G]A": "T[C>A]A",
    "T[C>A]C": "T[C>A]C",
    "G[T>G]A": "T[C>A]C",
    "T[C>A]G": "C[T>G]A",
    "C[T>G]A": "C[T>G]A",
    "T[C>A]T": "T[C>A]T",
    "A[T>G]A": "T[C>A]T",
    "A[C>T]A": "A[C>T]A",
    "A[T>C]A": "A[C>T]A",
    "A[C>T]C": "A[C>T]C",
    "A[T>C]C": "A[C>T]C",
    "A[C>T]G": "A[C>T]G",
    "A[T>C]G": "A[C>T]G",
    "A[C>T]T": "A[C>T]T",
    "A[T>C]T": "A[C>T]T",
    "C[C>T]A": "C[C>T]A",
    "C[T>C]A": "C[C>T]A",
    "C[C>T]C": "C[C>T]C",
    "C[T>C]C": "C[C>T]C",
    "C[C>T]G": "C[C>T]G",
    "C[T>C]G": "C[C>T]G",
    "C[C>T]T": "C[C>T]T",
    "C[T>C]T": "C[C>T]T",
    "G[C>T]A": "G[C>T]A",
    "G[T>C]A": "G[C>T]A",
    "G[C>T]C": "G[C>T]C",
    "G[T>C]C": "G[C>T]C",
    "G[C>T]G": "G[C>T]G",
    "G[T>C]G": "G[C>T]G",
    "G[C>T]T": "G[C>T]T",
    "G[T>C]T": "G[C>T]T",
    "T[C>T]A": "T[C>T]A",
    "T[T>C]A": "T[C>T]A",
    "T[C>T]C": "T[C>T]C",
    "T[T>C]C": "T[C>T]C",
    "T[C>T]G": "T[C>T]G",
    "T[T>C]G": "T[C>T]G",
    "T[C>T]T": "T[C>T]T",
    "T[T>C]T": "T[C>T]T",
    "A[C>G]A": "T[C>G]T",
    "T[C>G]T": "T[C>G]T",
    "A[C>G]C": "A[C>G]C",
    "G[C>G]T": "A[C>G]C",
    "A[C>G]G": "C[C>G]T",
    "C[C>G]T": "C[C>G]T",
    "A[C>G]T": "A[C>G]T",
    "C[C>G]A": "C[C>G]A",
    "T[C>G]G": "C[C>G]A",
    "C[C>G]C": "C[C>G]C",
    "G[C>G]G": "C[C>G]C",
    "C[C>G]G": "C[C>G]G",
    "G[C>G]A": "T[C>G]C",
    "T[C>G]C": "T[C>G]C",
    "G[C>G]C": "G[C>G]C",
    "T[C>G]A": "T[C>G]A",
    "A[T>A]A": "T[T>A]T",
    "T[T>A]T": "T[T>A]T",
    "A[T>A]C": "A[T>A]C",
    "G[T>A]T": "A[T>A]C",
    "A[T>A]G": "C[T>A]T",
    "C[T>A]T": "C[T>A]T",
    "A[T>A]T": "A[T>A]T",
    "C[T>A]A": "C[T>A]A",
    "T[T>A]G": "C[T>A]A",
    "C[T>A]C": "C[T>A]C",
    "G[T>A]G": "C[T>A]C",
    "C[T>A]G": "C[T>A]G",
    "G[T>A]A": "T[T>A]C",
    "T[T>A]C": "T[T>A]C",
    "G[T>A]C": "G[T>A]C",
    "T[T>A]A": "T[T>A]A",
}
SBS96_TO_TRI = {
    "A[C>A]A": "ACA",
    "A[C>A]T": "ACT",
    "A[C>A]G": "ACG",
    "A[C>A]C": "ACC",
    "T[C>A]A": "TCA",
    "T[C>A]T": "TCT",
    "T[C>A]G": "TCG",
    "T[C>A]C": "TCC",
    "G[C>A]A": "GCA",
    "G[C>A]T": "GCT",
    "G[C>A]G": "GCG",
    "G[C>A]C": "GCC",
    "C[C>A]A": "CCA",
    "C[C>A]T": "CCT",
    "C[C>A]G": "CCG",
    "C[C>A]C": "CCC",
    "A[C>G]A": "ACA",
    "A[C>G]T": "ACT",
    "A[C>G]G": "ACG",
    "A[C>G]C": "ACC",
    "T[C>G]A": "TCA",
    "T[C>G]T": "TCT",
    "T[C>G]G": "TCG",
    "T[C>G]C": "TCC",
    "G[C>G]A": "GCA",
    "G[C>G]T": "GCT",
    "G[C>G]G": "GCG",
    "G[C>G]C": "GCC",
    "C[C>G]A": "CCA",
    "C[C>G]T": "CCT",
    "C[C>G]G": "CCG",
    "C[C>G]C": "CCC",
    "A[C>T]A": "ACA",
    "A[C>T]T": "ACT",
    "A[C>T]G": "ACG",
    "A[C>T]C": "ACC",
    "T[C>T]A": "TCA",
    "T[C>T]T": "TCT",
    "T[C>T]G": "TCG",
    "T[C>T]C": "TCC",
    "G[C>T]A": "GCA",
    "G[C>T]T": "GCT",
    "G[C>T]G": "GCG",
    "G[C>T]C": "GCC",
    "C[C>T]A": "CCA",
    "C[C>T]T": "CCT",
    "C[C>T]G": "CCG",
    "C[C>T]C": "CCC",
    "A[T>A]A": "ATA",
    "A[T>A]T": "ATT",
    "A[T>A]G": "ATG",
    "A[T>A]C": "ATC",
    "T[T>A]A": "TTA",
    "T[T>A]T": "TTT",
    "T[T>A]G": "TTG",
    "T[T>A]C": "TTC",
    "G[T>A]A": "GTA",
    "G[T>A]T": "GTT",
    "G[T>A]G": "GTG",
    "G[T>A]C": "GTC",
    "C[T>A]A": "CTA",
    "C[T>A]T": "CTT",
    "C[T>A]G": "CTG",
    "C[T>A]C": "CTC",
    "A[T>C]A": "ATA",
    "A[T>C]T": "ATT",
    "A[T>C]G": "ATG",
    "A[T>C]C": "ATC",
    "T[T>C]A": "TTA",
    "T[T>C]T": "TTT",
    "T[T>C]G": "TTG",
    "T[T>C]C": "TTC",
    "G[T>C]A": "GTA",
    "G[T>C]T": "GTT",
    "G[T>C]G": "GTG",
    "G[T>C]C": "GTC",
    "C[T>C]A": "CTA",
    "C[T>C]T": "CTT",
    "C[T>C]G": "CTG",
    "C[T>C]C": "CTC",
    "A[T>G]A": "ATA",
    "A[T>G]T": "ATT",
    "A[T>G]G": "ATG",
    "A[T>G]C": "ATC",
    "T[T>G]A": "TTA",
    "T[T>G]T": "TTT",
    "T[T>G]G": "TTG",
    "T[T>G]C": "TTC",
    "G[T>G]A": "GTA",
    "G[T>G]T": "GTT",
    "G[T>G]G": "GTG",
    "G[T>G]C": "GTC",
    "C[T>G]A": "CTA",
    "C[T>G]T": "CTT",
    "C[T>G]G": "CTG",
    "C[T>G]C": "CTC",
}
SBS96_TRI_LST = list(set(SBS96_TO_TRI.values()))
SBS96_TRI_COUNT = len(SBS96_TRI_LST)
TRI_WEIGHT = 1/SBS96_TRI_COUNT
SBS52_LST = natsort.natsorted(list(set([sbs52 for (_sbs96, sbs52) in SBS96_TO_SBS52.items()])))
SBS52_TO_TRI_LST = {
    "T[T>G]T": ["ACA", "TTT"],
    "A[C>A]C": ["ACC", "GTT"],
    "C[T>G]T": ["ACG", "CTT"],
    "A[C>A]T": ["ACT", "ATT"],
    "C[C>A]A": ["CCA", "TTG"],
    "C[C>A]C": ["CCC", "GTG"],
    "C[C>A]G": ["CCG", "CTG"],
    "C[C>A]T": ["ATG", "CCT"],
    "T[T>G]C": ["GCA", "TTC"],
    "G[C>A]C": ["GCC", "GTC"],
    "C[T>G]C": ["CTC", "GCG"],
    "A[T>G]C": ["ATC", "GCT"],
    "T[C>A]A": ["TCA", "TTA"],
    "T[C>A]C": ["GTA", "TCC"],
    "C[T>G]A": ["CTA", "TCG"],
    "T[C>A]T": ["ATA", "TCT"],
    "A[C>T]A": ["ACA", "ATA"],
    "A[C>T]C": ["ACC", "ATC"],
    "A[C>T]G": ["ACG", "ATG"],
    "A[C>T]T": ["ACT", "ATT"],
    "C[C>T]A": ["CCA", "CTA"],
    "C[C>T]C": ["CCC", "CTC"],
    "C[C>T]G": ["CCG", "CTG"],
    "C[C>T]T": ["CCT", "CTT"],
    "G[C>T]A": ["GCA", "GTA"],
    "G[C>T]C": ["GCC", "GTC"],
    "G[C>T]G": ["GCG", "GTG"],
    "G[C>T]T": ["GCT", "GTT"],
    "T[C>T]A": ["TCA", "TTA"],
    "T[C>T]C": ["TCC", "TTC"],
    "T[C>T]G": ["TCG", "TTG"],
    "T[C>T]T": ["TCT", "TTT"],
    "T[C>G]T": ["ACA", "TCT"],
    "A[C>G]C": ["ACC", "GCT"],
    "C[C>G]T": ["ACG", "CCT"],
    "A[C>G]T": ["ACT"],
    "C[C>G]A": ["CCA", "TCG"],
    "C[C>G]C": ["CCC", "GCG"],
    "C[C>G]G": ["CCG"],
    "T[C>G]C": ["GCA", "TCC"],
    "G[C>G]C": ["GCC"],
    "T[C>G]A": ["TCA"],
    "T[T>A]T": ["ATA", "TTT"],
    "A[T>A]C": ["ATC", "GTT"],
    "C[T>A]T": ["ATG", "CTT"],
    "A[T>A]T": ["ATT"],
    "C[T>A]A": ["CTA", "TTG"],
    "C[T>A]C": ["CTC", "GTG"],
    "C[T>A]G": ["CTG"],
    "T[T>A]C": ["GTA", "TTC"],
    "G[T>A]C": ["GTC"],
    "T[T>A]A": ["TTA"]
}


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


def load_sbs52_counts(
    sbs96_count_file: str
):
    sbs52_counts = {sbs52: 0 for sbs52 in SBS52_LST}
    sbs96_counts = load_sbs96_counts(sbs96_count_file) 
    for (sbs96, count) in sbs96_counts.items():
        sbs52_counts[SBS96_TO_SBS52[sbs96]] += count
    return sbs52_counts


def load_trinucleotide_counts(
    tri_file: str
):
    trinucleotide_counts = {tri: 0 for tri in SBS96_TRI_LST}
    for line in open(tri_file).readlines():
        tri, count = line.strip().split()
        trinucleotide_counts[tri] = int(count) 
    trinucleotide_sum = sum(trinucleotide_counts.values())
    return trinucleotide_sum, trinucleotide_counts


def draw_sbs52_barplot(
    sbs52_file: str, 
    sbs52_pdf_file: str
):

    df = pd.read_csv(sbs52_file, sep="\t")
    plot = (
        ggplot(df, aes(x="TRI", y="WEIGHTED_COUNT", fill="SUB"))
        + geom_bar(stat="identity")
        + theme_bw()
        + facet_grid(". ~ SUB", scales="free")
        + scale_fill_manual(
            values=("#98D7EC", "#212121", "#FF003A", "#A6A6A6", "#F5ABCC")
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
    plot.save(sbs52_pdf_file, width=22, height=12)


def dump_weighted_sbs52_counts(
    sbs96_count_file: str, 
    tri_file: str,
    weighted_sbs52_count_file: str,
):

    sbs52_counts = load_sbs52_counts(sbs96_count_file)
    trinucleotide_sum, trinucleotide_counts = load_trinucleotide_counts(tri_file)
    with open(weighted_sbs52_count_file, "w") as outfile:
        print(
            "SUB",
            "TRI",
            "SBS52",
            "WEIGHTED_COUNT",
            "COUNT",
            "WEIGHT",
            sep="\t",
            file=outfile
        )
        for sbs52 in SBS52_LST:
            ubase, _, ref, _, alt, _, dbase = list(sbs52)
            tri = f"{ubase}{ref}{dbase}"
            sub = f"{ref}>{alt}"
            sbs52_tri_lst = SBS52_TO_TRI_LST[sbs52]
            sbs52_tri_sum = sum([trinucleotide_counts[tri] for tri in sbs52_tri_lst])
            sbs52_tri_freq = sbs52_tri_sum/float(trinucleotide_sum)
            sbs52_tri_weight = 1/(sbs52_tri_freq/TRI_WEIGHT)
            sbs52_count = sbs52_counts[sbs52]
            weighted_sbs52_count = sbs52_count * sbs52_tri_weight
            print(
                sub,
                tri,
                sbs52,
                weighted_sbs52_count,
                sbs52_count,
                sbs52_tri_weight,
                sep="\t",
                file=outfile
            )
    if weighted_sbs52_count_file.endswith(".tsv"):
        weighted_sbs52_count_barplot_file = weighted_sbs52_count_file.replace(".tsv", ".pdf")
    else:
        weighted_sbs52_count_barplot_file = weighted_sbs52_count_file + ".pdf" 
    draw_sbs52_barplot(weighted_sbs52_count_file, weighted_sbs52_count_barplot_file)
    
def main():
    options = parse_args(sys.argv)
    dump_weighted_sbs52_counts(options.input, options.tri, options.output)
    sys.exit(0)


if __name__ == "__main__":
    main()
