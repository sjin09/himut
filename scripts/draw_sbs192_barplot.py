#!/usr/bin/env python

import argparse
from collections import defaultdict
from typing import Dict
import sys

import natsort
import pandas as pd
from pathlib import Path
import plotnine as p9


PUR_SET = set(["A", "G"])
COMPLEMENTARY_BASE_LOOKUP = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}

NTS = ["A", "C", "G", "T"]
SUB_LST = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
SUB_TO_SBS96_LST = defaultdict(list)
SORTED_SBS96_LST = []
for sub in SUB_LST:
    for nti in NTS:
        for ntj in NTS:
            sbs96 = f"{nti}[{sub}]{ntj}"
            SUB_TO_SBS96_LST[sub].append(sbs96)
    SORTED_SBS96_LST.extend(natsort.natsorted(SUB_TO_SBS96_LST[sub]))

SBS192_FILL_COLOURS = (
    '#72CCE8', '#004F94',  # Lighter shade, darker shade, original color: "#98D7EC"
    '#595959', '#262626',
    '#FFAAA5', '#DB5665',  # Lighter shade, darker shade, original color: "#212121"
    '#BFBFBF', '#737373',  # Lighter shade, darker shade, original color: "#A6A6A6"
    '#AAD11E', '#769803',  # Lighter shade, darker shade, original color: "#83A603"
    '#F2ACC6', '#F279BC',  # Lighter shade, darker shade, original color: "#F5ABCC"
)

TRANSCRIPTIONAL_STRAND_BIAS_ORDER = [
    "T:C>A",
    "U:C>A",
    "T:C>G",
    "U:C>G",
    "T:C>T",
    "U:C>T",
    "T:T>A",
    "U:T>A",
    "T:T>C",
    "U:T>C",
    "T:T>G",
    "U:T>G",
]


def parse_args(args):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        required=True,
        help="SigProfiler SBS192 file to read"
    )
    parser.add_argument(
        "--sample",
        type=Path,
        required=True,
        help="sample id"
    )
    args = args[1:]
    return parser.parse_args(args)


def load_sbs192_counts(
    sbs192_file_path: Path,
) -> Dict[str, int]:

    sbs192_counts = {sbs192_state: defaultdict(lambda: 0) for sbs192_state in ["T", "U", "N"]}
    for line in open(sbs192_file_path, "r").readlines()[1:]:
        sbs192, count = line.split()
        sbs192_state, sbs96 = sbs192.split(":")
        if sbs192_state == "N":
            continue
        sbs192_counts[sbs192_state][sbs96] = int(count)
    return sbs192_counts


def draw_sbs192_barplot(
    sbs96_file: Path,
    sample: str,
    sbs96_pdf_file: str
):

    df = pd.read_csv(sbs96_file, sep="\t")
    df["TRANSCRIPTIONAL_STRAND_BIAS"] = pd.Categorical(
        df["TRANSCRIPTIONAL_STRAND_BIAS"],
        categories=TRANSCRIPTIONAL_STRAND_BIAS_ORDER,
        ordered=True
    )
    plot = (
        p9.ggplot(df, p9.aes(x="TRI", y="COUNT", fill="TRANSCRIPTIONAL_STRAND_BIAS"))
        + p9.geom_bar(stat="identity")
        + p9.theme_bw()
        + p9.facet_grid(". ~ SUB", scales="free")
        + p9.scale_fill_manual(
            values=SBS192_FILL_COLOURS
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
    plot.save(sbs96_pdf_file, width=22, height=12)


def write_sbs192_counts(
    sigprofiler_sbs192_file_path: Path,
    sample: str,
):
    sbs192_tsv_file_path = Path(str(sigprofiler_sbs192_file_path).replace(".SBS288.all", ".himut_unphased.sbs192.tsv"))
    sbs192_pdf_file_path = Path(str(sigprofiler_sbs192_file_path).replace(".SBS288.all", ".himut_unphased.sbs192.pdf"))
    sbs192_counts = load_sbs192_counts(sigprofiler_sbs192_file_path)
    with open(sbs192_tsv_file_path, "w") as outfile:
        print(
            "SUB",
            "TRI",
            "SBS96",
            "COUNT",
            "TRANSCRIPTIONAL_STRAND_BIAS",
            sep="\t",
            file=outfile
        )
        for sbs192_state in ["T", "U"]:
            sbs96_counts = sbs192_counts[sbs192_state]
            for sbs96 in SORTED_SBS96_LST:
                ubase, _, ref, _, alt, _, dbase = list(sbs96)
                sub = f"{ref}>{alt}"
                transcriptional_strand_bias = f"{sbs192_state}:{ref}>{alt}"
                tri = f"{ubase}{ref}{dbase}"
                print(
                    sub,
                    tri,
                    sbs96,
                    sbs96_counts[sbs96],
                    transcriptional_strand_bias,
                    sep="\t",
                    file=outfile
                )
    draw_sbs192_barplot(sbs192_tsv_file_path, sample, sbs192_pdf_file_path)


def main():
    options = parse_args(sys.argv)
    write_sbs192_counts(options.input, options.sample)
    sys.exit(0)


if __name__ == "__main__":
    main()
