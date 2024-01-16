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
        help="file to read hdp signature exposure per sample",
    )
    args = args[1:]
    return parser.parse_args(args)


def get_samples(
    hdp_exposure_csv_path: Path,
):

    samples = set()
    for line in open(hdp_exposure_csv_path).readlines()[1:]:
        fields = line.rstrip().split(",")
        sample = fields[0]
        samples.add(sample)
    samples = natsort.natsorted(list(samples))
    return samples


def get_hdp_signatures(
    hdp_exposure_csv_path: Path,
):
    hdp_signatures = open(
        hdp_exposure_csv_path
    ).readline().strip().split(",")[1:]
    hdp_signatures = [
        hdp_sig.replace("X", "HDP") for hdp_sig in hdp_signatures
    ]
    return hdp_signatures


def write_hdp_exposures(
    hdp_exposure_csv_path: Path,
    hdp_exposure_tsv_path: Path
):

    samples = []
    hdp_signatures = open(
        hdp_exposure_csv_path
    ).readline().strip().split(",")[1:]
    hdp_signatures = [
        f"HDP{hdp_sig}" for hdp_sig in hdp_signatures
    ]
    with open(hdp_exposure_tsv_path, "w") as outfile:
        print(
            "SAMPLE",
            "HDP_SIGNATURE",
            "SIGNATURE_EXPOSURE",
            sep="\t",
            file=outfile
        )
        for line in open(hdp_exposure_csv_path).readlines()[1:]:
            fields = line.rstrip().split(",")
            sample = fields[0]
            samples.append(sample)
            signature_exposures = fields[1:]
            for op_idx, sig_exposure in enumerate(signature_exposures):
                hdp_sig = hdp_signatures[op_idx]
                print(sample, hdp_sig, sig_exposure, sep="\t", file=outfile)


def load_hdp_exposures(
    hdp_exposure_tsv_path: Path,
):
    df = pd.read_csv(hdp_exposure_tsv_path, sep="\t")
    hdp_signatures = natsort.natsorted(list(set(df["HDP_SIGNATURE"])))
    df["HDP_SIGNATURE"] = pd.Categorical(
        df["HDP_SIGNATURE"],
        categories=hdp_signatures,
        ordered=True
    )
    df["SIGNATURE_EXPOSURE"] *= 100   # TODO
    return df


def get_hdp_exposure_p9_barplot(
    df: pd.DataFrame,
):
    plot = (
        p9.ggplot(
            df,
            p9.aes(
                x="SAMPLE",
                y="SIGNATURE_EXPOSURE",
                fill="HDP_SIGNATURE"
            )
        )
        + p9.scale_fill_manual(values=BARPLOT_FILL_VALUE)
        + p9.geom_bar(stat="identity")
        + p9.theme_bw()
        + p9.labs(x="\nSamples\n", y="\nHDP Signature Exposure (%)\n")
        + p9.geom_label(
            p9.aes(label="HDP_SIGNATURE"),
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


def draw_hdp_exposure_p9_barplot(hdp_exposure_csv_path: Path):
    # init
    samples = get_samples(hdp_exposure_csv_path)
    sample_count = len(samples)
    hdp_exposure_tsv_path = f"{hdp_exposure_csv_path.stem}.tsv"
    hdp_exposure_pdf_path = f"{hdp_exposure_csv_path.stem}.pdf"
    write_hdp_exposures(hdp_exposure_csv_path, hdp_exposure_tsv_path)
    df = load_hdp_exposures(hdp_exposure_tsv_path)
    # collect
    hdp_exposure_plots = []
    for op_jdx in range(0, sample_count, 15)[:-1]:
        sample_subset = samples[op_jdx:(op_jdx + 15)]
        df_subset = df[df["SAMPLE"].isin(sample_subset)]
        hdp_exposure_plot = get_hdp_exposure_p9_barplot(
            df_subset,
        )
        hdp_exposure_plots.append(hdp_exposure_plot)
    sample_subset = samples[(op_jdx + 15):sample_count]
    df_subset = df[df["SAMPLE"].isin(sample_subset)]
    hdp_exposure_plot = get_hdp_exposure_p9_barplot(
        df_subset,
    )
    hdp_exposure_plots.append(hdp_exposure_plot)
    # draw
    p9.save_as_pdf_pages(hdp_exposure_plots, hdp_exposure_pdf_path)


def main():
    options = parse_args(sys.argv)
    draw_hdp_exposure_p9_barplot(options.input)
    sys.exit(0)


if __name__ == "__main__":
    main()
