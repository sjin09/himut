#!/usr/bin/env python

import argparse
from collections import defaultdict
from typing import Dict, List
import sys

import natsort
import pandas as pd
from pathlib import Path
import plotnine as p9
import pysam

PUR_SET = set(["A", "G"])
PUR2PYR = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}

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
SBS96_LST.append("N[N>N]N")

SBS52_TO_SBS96_LST = defaultdict(list)
for SBS96, SBS52 in SBS96_TO_SBS52.items():
    SBS52_TO_SBS96_LST[SBS52].append(SBS96)

for SBS52 in SBS52_TO_SBS96_LST:
    SBS52_TO_SBS96_LST[SBS52] = natsort.natsorted(list(set(SBS52_TO_SBS96_LST[SBS52])))

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

SBS52_MUTSIG_FILL_COLOURS = ("#98D7EC", "#212121", "#FF003A", "#A6A6A6", "#F5ABCC")


class ExpandedVariantRecord:

    def __init__(self, variant: pysam.VariantRecord):
        self.variant = variant

    @property
    def chrom(self) -> str:
        return self.variant.chrom

    @property
    def pos(self) -> str:
        return self.variant.pos

    @property
    def ref(self) -> str:
        return self.variant.ref

    @property
    def alt(self) -> str:
        return self.alts[0]

    @property
    def alts(self) -> str:
        return self.variant.alts

    @property
    def filter(self) -> str:
        return str(list(self.variant.filter)[0])

    def is_snp(self, alt) -> bool:
        return True if len(self.ref) == 1 and len(alt) == 1 else False

    @property
    def is_pass(self) -> bool:
        return True if self.filter == "PASS" else False

    @property
    def is_biallelic(self) -> bool:
        return True if len(self.alts) == 1 else False

    @property
    def is_biallelic_snp(self) -> bool:
        if self.is_biallelic:
            return self.is_snp(self.alts[0])
        return False


def parse_args(args):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        required=True,
        help="VCF file to read"
    )
    parser.add_argument(
        "--ref",
        type=str,
        required=True,
        help="reference FASTA file to read"
    )
    parser.add_argument(
        "--region",
        type=Path,
        required=False,
        help="target chromosome"
    )
    parser.add_argument(
        "--region_list",
        type=Path,
        required=False,
        help="list of target chromosomes separated by new line"
    )
    args = args[1:]
    return parser.parse_args(args)


def get_sample(vcf_file_path: Path):

    vcf_file = pysam.VariantFile(vcf_file_path)
    sample = vcf_file.header.samples[0]
    vcf_file.close()
    return sample


def load_loci(
    region: str,
    region_list: str
):

    chrom_lst = []
    if region is None and region_list is not None:
        for line in open(region_list).readlines():
            chrom_lst.append(line.strip())
    elif region is not None and region_list is None:
        chrom_lst.append(region)
    elif region is not None and region_list is not None:
        for line in open(region_list).readlines():
            chrom_lst.append(line.strip())
    else:
        sys.exit()
    return natsort.natsorted(chrom_lst)


def get_sbs96(
    xvariant: ExpandedVariantRecord,
    reference_sequence_lookup: pysam.FastaFile,
):
    tri = reference_sequence_lookup.fetch(xvariant.chrom, xvariant.pos - 2, xvariant.pos + 1)
    if xvariant.ref in PUR_SET:
        ubase, _, dbase = tri[::-1]
        sbs96 = "{}[{}>{}]{}".format(
            PUR2PYR.get(ubase, "N"),
            PUR2PYR.get(xvariant.ref, "N"),
            PUR2PYR.get(xvariant.alt, "N"),
            PUR2PYR.get(dbase, "N"),
        )
    else:
        ubase, _, dbase = tri
        sbs96 = "{}[{}>{}]{}".format(ubase, xvariant.ref, xvariant.alt, dbase)
    if sbs96.count("N") == 0:
        return sbs96
    return "N[N>N]N"


def load_sbs96_counts(
    chrom_lst: List[str],
    vcf_file_path: Path,
    ref_file_path: Path,
) -> Dict[str, int]:
    vcf_file = pysam.VariantFile(vcf_file_path)
    reference_sequence_lookup = pysam.FastaFile(ref_file_path)
    if vcf_file_path.suffix == ".vcf":
        sbs96_per_chrom = defaultdict(lambda: {sbs96: 0 for sbs96 in SBS96_LST})
        for variant in vcf_file:
            xvariant = ExpandedVariantRecord(variant)
            if xvariant.is_pass and xvariant.is_biallelic_snp:
                sbs96 = get_sbs96(xvariant, reference_sequence_lookup)
                sbs96_per_chrom[variant.chrom][sbs96] += 1
    elif vcf_file_path.suffix == ".bgz":
        sbs96_per_chrom = {chrom: {sbs96: 0 for sbs96 in SBS96_LST} for chrom in chrom_lst}
        for chrom in chrom_lst:
            for variant in vcf_file.fetch(chrom):
                xvariant = ExpandedVariantRecord(variant)
                if xvariant.is_pass and xvariant.is_biallelic_snp:
                    sbs96 = get_sbs96(xvariant, reference_sequence_lookup)
                    sbs96_per_chrom[chrom][sbs96] += 1
    sbs96_counts = {sbs96: 0 for sbs96 in SBS96_LST}
    for chrom in chrom_lst:
        for sbs96, count in sbs96_per_chrom[chrom].items():
            sbs96_counts[sbs96] += count
    return sbs96_counts


def load_sbs52_counts(
    chrom_lst: List[str],
    vcf_file_path: Path,
    ref_file_path: Path,
):
    sbs52_counts = {sbs52: 0 for sbs52 in SBS52_LST}
    sbs96_counts = load_sbs96_counts(chrom_lst, vcf_file_path, ref_file_path)
    for sbs96, count in sbs96_counts.items():
        if sbs96.count("N") == 0:
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
        p9.ggplot(df, p9.aes(x="TRI", y="COUNT", fill="SUB"))
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


def dump_sbs52_counts(
    vcf_file_path: Path,
    ref_file_path: Path,
    region: Path,
    region_list: Path,
):

    sample = get_sample(vcf_file_path)
    chrom_lst = load_loci(region, region_list)
    sbs52_file_prefix = str(vcf_file_path).replace(".vcf", "").replace(".bgz", "")
    sbs52_file_path = Path(f"{sbs52_file_prefix}.sbs52.tsv")
    sbs52_pdf_file_path = Path(f"{sbs52_file_prefix}.sbs52.pdf")
    sbs52_counts = load_sbs52_counts(chrom_lst, vcf_file_path, ref_file_path)
    with open(sbs52_file_path, "w") as outfile:
        print(
            "SUB",
            "TRI",
            "SBS52",
            "COUNT",
            sep="\t",
            file=outfile
        )
        for sbs52 in SBS52_LST:
            ubase, _, ref, _, alt, _, dbase = list(sbs52)
            sub = f"{ref}>{alt}"
            tri = f"{ubase}{ref}{dbase}"
            tri_with_annotation = "({}) {}".format(";".join(SBS52_TO_SBS96_LST[sbs52]), tri)
            print(
                sub,
                tri_with_annotation,
                sbs52,
                sbs52_counts[sbs52],
                sep="\t",
                file=outfile
            )
    draw_sbs52_barplot(sbs52_file_path, sample, sbs52_pdf_file_path)


def main():
    options = parse_args(sys.argv)
    dump_sbs52_counts(options.input, options.ref, options.region, options.region_list)
    sys.exit(0)


if __name__ == "__main__":
    main()
