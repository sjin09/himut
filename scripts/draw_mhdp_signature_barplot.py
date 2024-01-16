import sys
import argparse
from pathlib import Path
from collections import defaultdict

import natsort
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
        help="file to read mHDP signatures",
    )
    args = args[1:]
    return parser.parse_args(args)


def load_mhdp_signatures(hdp_signature_csv_path: Path):
    mhdp_signatures = open(
        hdp_signature_csv_path
    ).readline().rstrip().split(",")
    mhdp_signatures = [
        "mHDP{}".format(mhdp_sig.split(".")[1]) for mhdp_sig in mhdp_signatures
    ]
    mhdp_signature_probabilities = {
        mhdp_sig: [] for mhdp_sig in mhdp_signatures
    }
    for line in open(hdp_signature_csv_path).readlines()[1:]:
        fields = line.rstrip().split(",")
        for op_idx, mutational_probability in enumerate(fields):
            mhdp_signature_probabilities[mhdp_signatures[op_idx]].append(
                "{:.3}".format(float(mutational_probability)*100)
            )
    return mhdp_signature_probabilities


def get_mhdp_signature_p9_barplot(
    title: str,
    mhdp_signature_tsv_path: Path
):
    df = pd.read_csv(mhdp_signature_tsv_path, sep="\t")
    plot = (
        p9.ggplot(df, p9.aes(x="TRI", y="MUTATIONAL_PROBABILITY", fill="SUB"))
        + p9.geom_bar(stat="identity")
        + p9.theme_bw()
        + p9.facet_grid(". ~ SUB", scales="free")
        + p9.scale_fill_manual(
            values=(
                "#98D7EC",
                "#212121",
                "#FF003A",
                "#A6A6A6",
                "#83A603",
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


def draw_mhdp_signature_barplot(
    mhdp_signature_csv_path: Path
):
    mhdp_signature_plots = []
    mhdp_signature_probabilities = load_mhdp_signatures(
        mhdp_signature_csv_path
    )
    for mhdp_signature in mhdp_signature_probabilities:
        mutational_probabilities = mhdp_signature_probabilities[mhdp_signature]
        mhdp_signature_tsv_path = Path(f"{mhdp_signature}.tsv")
        with open(mhdp_signature_tsv_path, "w") as outfile:
            print(
                "SUB",
                "TRI",
                "SBS96",
                "MUTATIONAL_PROBABILITY",
                sep="\t",
                file=outfile
            )
            for op_idx, mprob in enumerate(mutational_probabilities):
                sbs96 = SBS96_LST[op_idx]
                ubase, _, ref, _, alt, _, dbase = list(sbs96)
                sub = f"{ref}>{alt}"
                tri = f"{ubase}{ref}{dbase}"
                print(
                    sub,
                    tri,
                    sbs96,
                    mprob,
                    sep="\t",
                    file=outfile
                )
        mhdp_signature_plot = get_mhdp_signature_p9_barplot(
            mhdp_signature,
            mhdp_signature_tsv_path
        )
        mhdp_signature_plots.append(mhdp_signature_plot)
    mhdp_signature_plots_pdf_path = f"{mhdp_signature_csv_path.stem}.pdf"
    p9.save_as_pdf_pages(mhdp_signature_plots, mhdp_signature_plots_pdf_path)


def main():
    options = parse_args(sys.argv)
    draw_mhdp_signature_barplot(options.input)
    sys.exit(0)


if __name__ == "__main__":
    main()
