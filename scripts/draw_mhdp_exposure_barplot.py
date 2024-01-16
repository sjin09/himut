import argparse
from collections import defaultdict
from pathlib import Path
import sys

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

BARPLOT_FILL_VALUE = [
    "#FF0000",
    "#00FF00",
    "#0000FF",
    "#FFFF00",
    "#FF00FF",
    "#00FFFF",
    "#800000",
    "#008000",
    "#000080",
    "#808000",
    "#800080",
    "#008080",
    "#FF4500",
    "#DA70D6",
    "#FF8C00",
    "#FFD700",
    "#ADFF2F",
    "#32CD32",
    "#00FF7F",
    "#40E0D0",
    "#4682B4",
    "#87CEEB",
    "#7B68EE",
    "#8A2BE2",
    "#FF69B4",
    "#FF1493",
    "#00BFFF",
    "#1E90FF",
    "#20B2AA",
    "#98FB98",
    "#8FBC8F",
    "#2E8B57",
    "#6B8E23",
    "#556B2F",
    "#008B8B",
    "#800000",
    "#8B4513",
    "#A0522D",
    "#D2691E",
    "#CD853F",
    "#B8860B",
    "#DAA520",
    "#FFD700",
    "#FFA07A",
    "#FF7F50",
    "#FF6347",
    "#FF4500",
    "#FF8C00",
    "#FFD700",
    "#FFDAB9",
    "#FFE4B5",
    "#FFFACD",
    "#FAFAD2",
    "#FFEFD5",
    "#FFEBCD",
    "#FFE4C4",
    "#FFDEAD",
    "#F5DEB3",
    "#DEB887",
    "#D2B48C",
    "#BC8F8F",
    "#F08080",
    "#CD5C5C",
    "#8B0000",
    "#A52A2A",
    "#B22222",
    "#DC143C",
    "#FF0000",
    "#FF4500",
    "#FF6347",
    "#FF7F50",
    "#CD5C5C",
    "#DC143C",
    "#B22222",
    "#8B0000",
    "#A52A2A",
    "#D2691E",
    "#FF8C00",
    "#FFA07A",
    "#FFD700",
    "#FFDAB9",
    "#F5DEB3",
    "#DEB887",
    "#D2B48C",
    "#BC8F8F",
    "#F08080",
    "#CD853F",
    "#DAA520",
    "#FF6347",
    "#FF4500",
    "#FF8C00",
    "#FFA07A",
    "#FA8072",
    "#E9967A",
    "#F08080",
    "#CD5C5C",
    "#DC143C",
    "#B22222",
    "#8B0000",
    "#800000",
]

RANGE_STEP = 15


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
        help="file to read mhdp signature exposure per sample",
    )
    args = args[1:]
    return parser.parse_args(args)


def get_samples(
    mhdp_exposure_csv_path: Path,
):

    samples = set()
    samples = open(mhdp_exposure_csv_path).readline().rstrip().split(",")[1:]
    samples = [sample.replace('"', '') for sample in samples]
    return samples


def get_mhdp_signatures(
    mhdp_exposure_csv_path: Path,
):
    mhdp_signatures = []
    for line in open(mhdp_exposure_csv_path).readlines()[1:]:
        fields = line.rstrip().split(",")
        mhdp_sig = fields[0]
        mhdp_sig = mhdp_sig.replace('"', '')
        mhdp_signatures.append(mhdp_sig)
    return mhdp_signatures


def write_mhdp_exposures(
    mhdp_exposure_csv_path: Path,
    mhdp_exposure_tsv_path: Path
):
    # init
    samples = get_samples(mhdp_exposure_csv_path)
    mhdp_signatures = get_mhdp_signatures(mhdp_exposure_csv_path)
    mhdp_signature_burden_per_sample = {}
    for sample in samples:
        mhdp_signature_burden_per_sample[sample] = {}
        for mhdp_sig in mhdp_signatures:
            mhdp_signature_burden_per_sample[sample][mhdp_sig] = 0
    # collect
    for line in open(mhdp_exposure_csv_path).readlines()[1:]:
        fields = line.rstrip().split(",")
        mhdp_sig = fields[0].replace('"', '')
        mutation_burdens = fields[1:]
        for op_idx, mburden in enumerate(mutation_burdens):
            sample = samples[op_idx]
            mhdp_signature_burden_per_sample[sample][mhdp_sig] = mburden
    # write
    with open(mhdp_exposure_tsv_path, "w") as outfile:
        print(
            "SAMPLE",
            "mHDP_SIGNATURE",
            "MUTATION_BURDEN",
            sep="\t",
            file=outfile
        )
        for sample in samples:
            for mhdp_sig in mhdp_signatures:
                mburden = mhdp_signature_burden_per_sample[sample][mhdp_sig]
                mhdp_sig = "mHDP{}".format(mhdp_sig.split(".")[1])
                print(
                    sample,
                    mhdp_sig,
                    mburden,
                    sep="\t",
                    file=outfile
                )


def load_mhdp_exposures(
    mhdp_exposure_tsv_path: Path,
):
    df = pd.read_csv(mhdp_exposure_tsv_path, sep="\t")
    mhdp_signatures = natsort.natsorted(list(set(df["mHDP_SIGNATURE"])))
    df["mHDP_SIGNATURE"] = pd.Categorical(
        df["mHDP_SIGNATURE"],
        categories=mhdp_signatures,
        ordered=True
    )
    return df


def get_mhdp_exposure_p9_barplot(
    df: pd.DataFrame,
):
    plot = (
        p9.ggplot(
            df,
            p9.aes(
                x="SAMPLE",
                y="MUTATION_BURDEN",
                fill="mHDP_SIGNATURE"
            )
        )
        + p9.scale_fill_manual(values=BARPLOT_FILL_VALUE)
        + p9.geom_bar(stat="identity")
        + p9.theme_bw()
        + p9.labs(x="\nSamples\n", y="\nmHDP Mutation Burden\n")
        + p9.geom_label(
            p9.aes(label="mHDP_SIGNATURE"),
            position="stack",
            size=8,
            format_string='{}'
        )
        + p9.theme(
            figure_size=(24, 14),
            text=p9.element_text(size=10),
            legend_title=p9.element_blank(),
        )
    )
    return plot


def draw_mhdp_exposure_p9_barplot(mhdp_exposure_csv_path: Path):
    # init
    samples = get_samples(mhdp_exposure_csv_path)
    sample_count = len(samples)
    mhdp_exposure_tsv_path = f"{mhdp_exposure_csv_path.stem}.tsv"
    mhdp_exposure_pdf_path = f"{mhdp_exposure_csv_path.stem}.pdf"
    write_mhdp_exposures(mhdp_exposure_csv_path, mhdp_exposure_tsv_path)
    df = load_mhdp_exposures(mhdp_exposure_tsv_path)
    # collect
    mhdp_exposure_plots = []
    if sample_count < RANGE_STEP:
        hdp_exposure_plot = get_mhdp_exposure_p9_barplot(
            df,
        )
        mhdp_exposure_plots.append(hdp_exposure_plot)
    else:
        for op_idx in range(0, sample_count, 15)[:-1]:
            sample_subset = samples[op_idx:(op_idx + 15)]
            df_subset = df[df["SAMPLE"].isin(sample_subset)]
            mhdp_exposure_plot = get_mhdp_exposure_p9_barplot(
                df_subset,
            )
            mhdp_exposure_plots.append(mhdp_exposure_plot)
        sample_subset = samples[(op_idx + 15):sample_count]
        df_subset = df[df["SAMPLE"].isin(sample_subset)]
        mhdp_exposure_plot = get_mhdp_exposure_p9_barplot(df_subset)
        mhdp_exposure_plots.append(mhdp_exposure_plot)
    # draw
    p9.save_as_pdf_pages(mhdp_exposure_plots, mhdp_exposure_pdf_path)


def main():
    options = parse_args(sys.argv)
    draw_mhdp_exposure_p9_barplot(options.input)
    sys.exit(0)


if __name__ == "__main__":
    main()
