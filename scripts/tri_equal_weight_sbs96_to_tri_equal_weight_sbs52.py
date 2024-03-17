#!/usr/bin/env python

import argparse
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict

import natsort
import pandas as pd
import plotnine as p9


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

SBS52_TO_SBS96_LST = defaultdict(list)
SBS52_SUB_LST = ["C>A", "C>G", "C>T", "T>A", "T>G"]
for SBS96, SBS52 in SBS96_TO_SBS52.items():
    SBS52_TO_SBS96_LST[SBS52].append(SBS96)

SUB_TO_SBS52_LST = defaultdict(list)
for SBS52 in set(list(SBS96_TO_SBS52.values())):
    ubase, _, ref, _, alt, _, dbase = list(SBS52)
    sub = "{}>{}".format(ref, alt)
    tri = "{}{}{}".format(ubase, ref, dbase)
    SUB_TO_SBS52_LST[sub].append(SBS52)

SBS52_LST = []
for sbs52_sub in SBS52_SUB_LST:
    SUB_TO_SBS52_LST[sbs52_sub] = natsort.natsorted(SUB_TO_SBS52_LST[sbs52_sub])
    SBS52_LST.extend(SUB_TO_SBS52_LST[sbs52_sub])

for SBS52 in SBS52_TO_SBS96_LST:
    SBS52_TO_SBS96_LST[SBS52] = natsort.natsorted(list(set(SBS52_TO_SBS96_LST[SBS52])))

TRI_WITH_ANNOTATION = []
for sbs52_sub in SBS52_SUB_LST:
    for sbs52 in SUB_TO_SBS52_LST[sbs52_sub]:
        ubase, _, ref, _, alt, _, dbase = list(sbs52)
        tri = f"{ubase}{ref}{dbase}"
        TRI_WITH_ANNOTATION.append("({}) {}".format(";".join(SBS52_TO_SBS96_LST[sbs52]), tri))

SBS52_MUTSIG_FILL_COLOURS = ("#98D7EC", "#212121", "#FF003A", "#A6A6A6", "#F5ABCC")


def parse_args(args):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        required=True,
        help="file to read trinucleotide normalised SBS96 counts"
    )
    parser.add_argument(
        "--sample",
        type=str,
        required=True,
        help="sample id"
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="file to read trinucleotide normalised SBS52 counts"
    )
    args = args[1:]
    return parser.parse_args(args)


def load_sbs96_counts(sbs92_file_path: Path) -> Dict[str, int]:
    sbs96_counts = {}
    for line in open(sbs92_file_path).readlines():
        if line.startswith("SUB"):
            continue
        _sub, _tri, sbs96, normcount, _count, _weight = line.rstrip().split("\t")
        sbs96_counts[sbs96] = float(normcount)
    return sbs96_counts


def load_sbs52_counts(sbs92_file_path: Path) -> Dict[str, int]:
    sbs52_counts = {sbs52: 0 for sbs52 in SBS52_LST}
    sbs96_counts = load_sbs96_counts(sbs92_file_path)
    for sbs96, count in sbs96_counts.items():
        if sbs96.count("N") != 0:
            continue
        sbs52_counts[SBS96_TO_SBS52[sbs96]] += count
    return sbs52_counts


def draw_sbs52_barplot(
    sbs52_file: Path,
    sample: str,
    sbs52_pdf_file: str
):

    df = pd.read_csv(sbs52_file, sep="\t")
    df["TRI"] = pd.Categorical(
        df["TRI"],
        categories=TRI_WITH_ANNOTATION,
        ordered=True
    )
    plot = (
        p9.ggplot(df, p9.aes(x="TRI", y="NORMCOUNT", fill="SUB"))
        + p9.geom_bar(stat="identity")
        + p9.theme_bw()
        + p9.facet_grid(". ~ SUB", scales="free")
        + p9.scale_fill_manual(
            values=SBS52_MUTSIG_FILL_COLOURS
        )
        + p9.labs(x="\nTrinucleotide Context\n", y="\nCounts\n")
        + p9.ggtitle("\n{}\n".format(sample))
        + p9.theme(
            text=p9.element_text(size=10),
            legend_title=p9.element_blank(),
            axis_text_x=p9.element_text(
                family="monospace", angle=90, ha="center"
            ),
        )
    )
    plot.save(sbs52_pdf_file, width=22, height=12)


def write_tri_equal_weight_sbs52_counts(
    input_path: Path,
    sample: str,
    output_path: Path
):
    sbs52_counts = load_sbs52_counts(input_path)
    sbs52_pdf_file_path = Path("{}.pdf".format(output_path.stem))
    with open(output_path, "w") as outfile:
        print(
            "SUB",
            "TRI",
            "SBS52",
            "NORMCOUNT",
            sep="\t",
            file=outfile
        )
        for sbs52 in SBS52_LST:
            ubase, _, ref, _, alt, _, dbase = list(sbs52)
            sub = f"{ref}>{alt}"
            tri = f"{ubase}{ref}{dbase}"
            tri_with_annotation = "({}) {}".format(";".join(SBS52_TO_SBS96_LST[sbs52]), tri)
            sbs52_count = sbs52_counts[sbs52]
            print(
                sub,
                tri_with_annotation,
                sbs52,
                sbs52_count,
                sep="\t",
                file=outfile
            )
    draw_sbs52_barplot(output_path, sample, sbs52_pdf_file_path)


def main():
    options = parse_args(sys.argv)
    write_tri_equal_weight_sbs52_counts(options.input, options.sample, options.output)
    sys.exit(0)


if __name__ == "__main__":
    main()
