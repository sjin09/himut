import sys
import argparse
from pathlib import Path
from collections import defaultdict

import natsort
import numpy as np
import pandas as pd
import plotnine as p9

NTS = ["A", "C", "G", "T"]
SUB_LST = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
SUB_TO_SBS96_LST = defaultdict(list)
for sub in SUB_LST:
    for nti in NTS:
        for ntj in NTS:
            sbs96 = f"{nti}[{sub}]{ntj}"
            SUB_TO_SBS96_LST[sub].append(sbs96)

SBS96_LST = []
for sub in SUB_LST:
    SUB_TO_SBS96_LST[sub] = natsort.natsorted(SUB_TO_SBS96_LST[sub])
    SBS96_LST.extend(SUB_TO_SBS96_LST[sub])

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
for SBS96, SBS52 in SBS96_TO_SBS52.items():
    SBS52_TO_SBS96_LST[SBS52].append(SBS96)

for SBS52 in SBS52_TO_SBS96_LST:
    SBS52_TO_SBS96_LST[SBS52] = natsort.natsorted(list(set(SBS52_TO_SBS96_LST[SBS52])))


SUB_LST = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
SUB_TO_SBS52_LST = {sub: [] for sub in SUB_LST}
for sbs52 in list(set(SBS96_TO_SBS52.values())):
    ubase, _, ref, _, alt, _, dbase = list(sbs52)
    sub = "{}>{}".format(ref, alt)
    SUB_TO_SBS52_LST[sub].append(sbs52)

SBS52_LST = []
for sub in SUB_LST:
    SUB_TO_SBS52_LST[sub] = natsort.natsorted(SUB_TO_SBS52_LST[sub])
    SBS52_LST.extend(SUB_TO_SBS52_LST[sub])

TRI_WITH_ANNOTATION = []
for sub in SUB_LST:
    for sbs52 in SUB_TO_SBS52_LST[sub]:
        ubase, _, ref, _, alt, _, dbase = list(sbs52)
        tri = f"{ubase}{ref}{dbase}"
        TRI_WITH_ANNOTATION.append("({}) {}".format(";".join(SBS52_TO_SBS96_LST[sbs52]), tri))


def parse_args(args):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        required=True,
        help="file to read HDP signatures",
    )
    args = args[1:]
    return parser.parse_args(args)


def get_hdp_signature_p9_plot(
    title: str,
    hdp_signature_tsv_path: Path
):

    df = pd.read_csv(hdp_signature_tsv_path, sep="\t")
    df["TRI"] = pd.Categorical(
        df["TRI"],
        categories=TRI_WITH_ANNOTATION,
        ordered=True
    )
    plot = (
        p9.ggplot(df, p9.aes(x="TRI", y="PROBABILITY", fill="SUB"))
        + p9.geom_bar(stat="identity")
        + p9.theme_bw()
        + p9.facet_grid(". ~ SUB", scales="free")
        + p9.scale_fill_manual(
            values=(
                "#98D7EC",
                "#212121",
                "#FF003A",
                "#A6A6A6",
                "#F5ABCC"
                )
        )
        + p9.labs(
            title=title,
            x="\nTrinucleotide Context\n",
            y="\nMutational Probability (%)\n"
        )
        + p9.theme(
            figure_size=(22, 12),
            text=p9.element_text(size=10),
            legend_title=p9.element_blank(),
            axis_text_x=p9.element_text(
                family="monospace", angle=90, ha="center"
            ),
        )
    )
    return plot


def load_hdp_sbs96_signatures(hdp_signature_csv_path: Path):
    hdp_idxes = open(hdp_signature_csv_path).readline().strip().split(",")[1:]
    hdp_signatures = [hdp_idx.replace("X", "HDP") for hdp_idx in hdp_idxes]
    hdp_sig_mutational_probabilities = {
        hdp_sig: [] for hdp_sig in hdp_signatures
    }
    for line in open(hdp_signature_csv_path).readlines()[1:]:
        fields = line.rstrip().split(",")
        for op_idx, mutational_probability in enumerate(fields[1:]):
            hdp_sig = hdp_signatures[op_idx]
            hdp_sig_mutational_probabilities[hdp_sig].append(
                "{:.3}".format(float(mutational_probability)*100)
            )
    return hdp_sig_mutational_probabilities


def load_hdp_sbs52_signatures(hdp_signature_csv_path: Path):
    hdp_sbs52_sig_mutational_probabilities = {}
    hdp_sbs96_sig_mutational_probabilities = load_hdp_sbs96_signatures(
        hdp_signature_csv_path
    )
    for hdp_sig, sbs96_mprobs in hdp_sbs96_sig_mutational_probabilities.items():
        sbs52_to_mprob = np.zeros((len(SBS52_LST), 1))
        for op_idx, mprob in enumerate(sbs96_mprobs):
            sbs96 = SBS96_LST[op_idx]
            sbs52 = SBS96_TO_SBS52[sbs96]
            sbs52_idx = SBS52_LST.index(sbs52)
            sbs52_to_mprob[sbs52_idx] += float(mprob)
        hdp_sbs52_sig_mutational_probabilities[hdp_sig] = sbs52_to_mprob
    return hdp_sbs52_sig_mutational_probabilities


def write_and_draw_transformed_hdp_signatures(
    hdp_signature_csv_path: Path,
):
    hdp_signature_plots = []
    hdp_sbs52_sig_mutational_probabilities = load_hdp_sbs52_signatures(
        hdp_signature_csv_path
    )
    hdp_signature_plots_pdf_path = f"{hdp_signature_csv_path.stem}.sbs52.pdf"
    for hdp_sig, mprobs in hdp_sbs52_sig_mutational_probabilities.items():
        hdp_signature_tsv_path = f"{hdp_sig}.sbs52.mutational_probabilities.tsv"
        with open(hdp_signature_tsv_path, "w") as outfile:
            print(
                "SUB",
                "TRI",
                "SBS52",
                "PROBABILITY",
                sep="\t",
                file=outfile
            )
            for op_idx, mprob in enumerate(mprobs):
                sbs52 = SBS52_LST[op_idx]
                ubase, _, ref, _, alt, _, dbase = list(sbs52)
                sub = f"{ref}>{alt}"
                tri = f"{ubase}{ref}{dbase}"
                tri_with_annotation = "({}) {}".format(";".join(SBS52_TO_SBS96_LST[sbs52]), tri)
                print(
                    sub,
                    tri_with_annotation,
                    sbs52,
                    mprob[0],
                    sep="\t",
                    file=outfile
                )
        hdp_signature_plot = get_hdp_signature_p9_plot(
            hdp_sig,
            hdp_signature_tsv_path
        )
        hdp_signature_plots.append(hdp_signature_plot)
    p9.save_as_pdf_pages(hdp_signature_plots, hdp_signature_plots_pdf_path)


def main():
    options = parse_args(sys.argv)
    write_and_draw_transformed_hdp_signatures(options.input)
    sys.exit(0)


if __name__ == "__main__":
    main()
