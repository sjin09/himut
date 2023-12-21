import sys
import natsort
import argparse
import pandas as pd
from plotnine import *
from collections import defaultdict


def parse_args(args):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-i",
        "--sbs96",
        type=str,
        required=True,
        help="file to read SBS96 normalised counts",
    )
    parser.add_argument(
        "-o",
        "--sbs52",
        type=str,
        required=True,
        help="file to return normalised SBS52 counts",
    )
    args = args[1:]
    return parser.parse_args(args)


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
SBS52_LST = natsort.natsorted(list(set([sbs52 for (_sbs96, sbs52) in SBS96_TO_SBS52.items()])))



def load_sbs52_counts(sbs96_file):

    sbs52_counts = {sbs52: 0 for sbs52 in SBS52_LST}
    for line in open(sbs96_file):
        if line.startswith("#"):
            continue
        elif line.startswith("sub"):
            continue
        (
            _sub,
            _tri,
            sbs96,
            _count,
            normcount,
            _ref_tri_ratio,
            _ref_ccs_tri_ratio,
            _ref_tri_count,
            _ref_callable_tri_count,
            _ccs_callable_tri_count,
        ) = line.strip().split()
        sbs52 = SBS96_TO_SBS52[sbs96]
        sbs52_counts[sbs52] += float(normcount)
    return sbs52_counts


def dump_sbs52_plt(
    infile: str, 
    outfile: str
):

    df = pd.read_csv(infile, sep="\t")
    plot = (
        ggplot(df, aes(x="TRI", y="COUNT", fill="SUB"))
        + geom_bar(stat="identity")
        + theme_bw()
        + facet_grid(". ~ SUB", scales="free")
        + scale_fill_manual(
            values=("#98D7EC", "#212121", "#FF003A", "#A6A6A6", "#83A603", "#F5ABCC")
        )
        + labs(x="\nTrinucleotide Context\n", y="\nWeighted Count\n")
        + theme(
            text=element_text(size=10),
            legend_title=element_blank(),
            axis_text_x=element_text(
                family="monospace", angle=90, ha="center"
            ),
        )
    )
    plot.save(outfile, width=22, height=12)


def dump_normalised_sbs52_counts(
    sbs96_file: str,
    sbs52_file: str,
):

    sbs52_counts= load_sbs52_counts(sbs96_file)
    with open(sbs52_file, "w") as outfile:
        print(
            "SUB",
            "TRI",
            "SBS52",
            "COUNT",
            sep="\t",
            file=outfile
        )
       
        for sbs52 in SBS52_LST:
            _ubase, _, ref, _, alt, _, _dbase = list(sbs52)
            sub = f"{ref}>{alt}"
            sbs52_normcount = sbs52_counts[sbs52]
            print(
                sub,
                SBS96_TO_TRI[sbs52],
                sbs52,
                sbs52_normcount,
                sep="\t",
                file=outfile
            )
            
    if sbs52_file.endswith(".tsv"):
        sbs52_pdf_file = sbs52_file.replace(".tsv", ".pdf")
    else:
        sbs52_pdf_file = sbs52_file + ".pdf"
    dump_sbs52_plt(sbs52_file, sbs52_pdf_file)


def main():
    options = parse_args(sys.argv)
    dump_normalised_sbs52_counts(options.sbs96, options.sbs52)
    sys.exit(0)


if __name__ == "__main__":
    main()
