import sys
import natsort
import argparse
import itertools
import pandas as pd
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
        help="file to read somatic and germline SBS48 counts"
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
        type=str,
        required=True,
        help="file to generate mutational spectrum plot file"
    )
    args = args[1:]
    return parser.parse_args(args)

    
def dump_sbs48_plt(infile: str, sample: str, outfile: str) -> None:

    df = pd.read_csv(infile, sep="\t")
    plot = (ggplot(df, aes(x="tri", y="counts", fill="sub")) +
      geom_bar(stat="identity")  +
      theme_bw() +
      facet_grid("source ~ sub", scales = "free") +
      scale_fill_manual(values = ("#98D7EC","#212121","#FF003A","#A6A6A6","#83A603","#F5ABCC")) +
      labs(x = "\nTrinucleotide sequence context\n", y = "\nCounts\n") +
      ggtitle("\n{}\n".format(sample)) +
      theme(
          text = element_text(size=10),
          legend_title = element_blank(),
          axis_text_x = element_text(family = "monospace", angle = 90, ha="center")
        )
    )
    plot.save(outfile, width = 22, height = 12)


def main():
    options = parse_args(sys.argv)
    dump_sbs48_plt(options.input, options.sample, options.output)
    sys.exit(0)


if __name__ == "__main__":
    main()

