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
        help="file to read HDP signatures",
    )
    args = args[1:]
    return parser.parse_args(args)


def get_hdp_signature_p9_plot(
    title: str,
    hdp_signature_tsv_path: Path
):

    df = pd.read_csv(hdp_signature_tsv_path, sep="\t")
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


def load_hdp_signatures(hdp_signature_csv_path: Path):
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


def write_and_draw_transformed_hdp_signatures(
    hdp_signature_csv_path: Path,
):
    hdp_signature_plots = []
    hdp_sig_mutational_probabilities = load_hdp_signatures(
        hdp_signature_csv_path
    )
    hdp_signature_plots_pdf_path = f"{hdp_signature_csv_path.stem}.sbs96.pdf"
    for hdp_sig, mprobs in hdp_sig_mutational_probabilities.items():
        hdp_signature_tsv_path = f"{hdp_sig}.sbs96.mutational_probabilities.tsv"
        with open(hdp_signature_tsv_path, "w") as outfile:
            print(
                "SUB",
                "TRI",
                "SBS96",
                "PROBABILITY",
                sep="\t",
                file=outfile
            )
            for op_idx, mprob in enumerate(mprobs):
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
