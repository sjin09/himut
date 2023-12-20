#!/usr/bin/env python

import argparse
import cyvcf2
import natsort
import pandas as pd
import pyfastx
import sys

from collections import defaultdict 
from plotnine import *
from typing import Dict, List, Tuple


PURINE = set(["A", "G"])
PURINE2PYRIMIDINE = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
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
SBS96_LST = list(SBS96_TO_SBS52.keys())
SBS52_LST = list(set(SBS96_TO_SBS52.values()))


class VCF:
    def __init__(self, line):
        arr = line.strip().split()
        self.chrom = arr[0]
        self.pos = int(arr[1])
        self.id = arr[2]
        self.ref = arr[3]
        self.alt_lst = arr[4].split(",")
        self.qual = arr[5]
        self.qual = float(self.qual) if self.qual != "." else self.qual
        self.is_pass = True if arr[6] == "PASS" else False
        self.info = arr[7]
        self.format_lst = arr[8].split(":")
        self.sample_format_lst = arr[9].split(":")
        hsh = {i: j for i, j in zip(self.format_lst, self.sample_format_lst)}
        if "GT" in hsh:
            self.sample_gt = hsh["GT"]
        if "PS" in hsh:
            self.sample_phase_set = hsh["PS"]
        if "AD" in hsh:
            arr = hsh["AD"].split(",")
            self.ref_count = arr[0]
            self.alt_count_arr = arr[1:]

        self.is_snp = False
        self.is_dbs = False
        self.is_indel = False
        if len(self.alt_lst) == 1:  # bi-allelic
            self.is_biallelic = True
            self.alt = self.alt_lst[0]
            if len(self.ref) == 1 and len(self.alt) == 1:  # snp
                self.is_snp = True
            elif len(self.ref) == len(self.alt) == 2:
                self.is_dbs = True
            elif len(self.ref) > len(self.alt):  # del
                self.is_indel = True
            elif len(self.ref) < len(self.alt):  # ins
                self.is_indel = True
        else:
            self.is_biallelic = False


def parse_args(args):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-i",
        "--input",
        type=str,
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
        type=str,
        required=False,
        help="target chromosome"
    )
    parser.add_argument(
        "--region_list",
        type=str,
        required=False,
        help="list of target chromosomes separated by new line"
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="file to return SBS96 counts"
    )
    args = args[1:]
    return parser.parse_args(args)


def get_sample(vcf_file: str):

    if vcf_file.endswith(".vcf"):
        for line in open(vcf_file).readlines():
            if line.startswith("##"):
                continue
            elif line.startswith("#CHROM"):
                sample = line.strip().split()[-1]
                break
    elif vcf_file.endswith(".bgz"):
        v = cyvcf2.VCF(vcf_file)
        sample = v.samples[0]
    return sample


def get_tname2tsize(vcf_file: str):

    tname2tsize = {}
    if vcf_file.endswith(".vcf"):
        for line in open(vcf_file).readlines():
            if line.startswith("##"):
                if line.startswith("##contig"):
                    arr = line.strip().replace("##contig=<ID=", "").split(",")
                    tname = arr[0]
                    tsize = int(arr[1].replace("length=", "").replace(">", ""))
                    tname2tsize[tname] = tsize
            elif line.startswith("#CHROM"):
                break
    elif vcf_file.endswith(".bgz"):
        for line in cyvcf2.VCF(vcf_file).raw_header.split("\n"):
            if line.startswith("##"):
                if line.startswith("##contig"):
                    arr = line.strip().replace("##contig=<ID=", "").split(",")
                    tname = arr[0]
                    tsize = int(arr[1].replace("length=", "").replace(">", ""))
                    tname2tsize[tname] = tsize
            elif line.startswith("#CHROM"):
                break
    return tname2tsize


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
    chrom: str,
    pos: str,
    ref: str,
    alt: str,
    refseq,
):
    if ref in PURINE:
        ubase = PURINE2PYRIMIDINE.get(refseq[chrom][pos + 1], "N")
        dbase = PURINE2PYRIMIDINE.get(refseq[chrom][pos - 1], "N")
        sbs96 = "{}[{}>{}]{}".format(
            ubase,
            PURINE2PYRIMIDINE.get(ref, "N"),
            PURINE2PYRIMIDINE.get(alt, "N"),
            dbase,
        )
    else:
        ubase = refseq[chrom][pos - 1]
        dbase = refseq[chrom][pos + 1]
        sbs96 = "{}[{}>{}]{}".format(ubase, ref, alt, dbase)
    return sbs96


def load_sbs96_counts(
    vcf_file: str, ref_file: str, chrom_lst: List[str]
) -> Dict[str, int]:

    tname_lst = []
    tname2sbs96_counts = defaultdict()
    refseq = pyfastx.Fasta(ref_file)
    if vcf_file.endswith(".vcf"):
        for line in open(vcf_file).readlines():
            if line.startswith("##"):
                if line.startswith("##contig"):
                    tname = line.strip().replace("##contig=<ID=", "").split(",")[0]
                    tname_lst.append(tname)
                continue
            elif line.startswith("#CHROM"):
                for tname in tname_lst:
                    tname2sbs96_counts[tname] = defaultdict(lambda: 0)
                continue
            v = VCF(line)
            if v.is_snp and v.is_pass and v.is_biallelic:
                sbs96 = get_sbs96(v.chrom, int(v.pos) - 1, v.ref, v.alt, refseq)
                tname2sbs96_counts[v.chrom][sbs96] += 1
    elif vcf_file.endswith(".vcf.bgz"):
        for line in cyvcf2.VCF(vcf_file).raw_header.split("\n"):
            if line.startswith("##"):
                if line.startswith("##contig"):
                    tname = line.replace("##contig=<ID=", "").split(",")[0]
                    tname_lst.append(tname)
                continue
            elif line.startswith("#CHROM"):
                for tname in tname_lst:
                    tname2sbs96_counts[tname] = defaultdict(lambda: 0)
        for i in cyvcf2.VCF(vcf_file):
            v = VCF(str(i))
            if v.is_snp and v.is_pass and v.is_biallelic:
                sbs96 = get_sbs96(v.chrom, int(v.pos) - 1, v.ref, v.alt, refseq)
                tname2sbs96_counts[v.chrom][sbs96] += 1

    sbs96_counts = {sbs96: 0 for sbs96 in SBS96_LST}
    for chrom in chrom_lst:
        for sbs96, count in tname2sbs96_counts[chrom].items():
            if sbs96.count("N") > 0:
                continue
            sbs96_counts[sbs96] += count
    return sbs96_counts


def load_sbs52_counts(
    sbs96_counts: Dict[str, int]
):
    sbs52_counts = {sbs52: 0 for sbs52 in SBS52_LST}
    for sbs96 in SBS96_LST:
        sbs52 = SBS96_TO_SBS52[sbs96]
        sbs96_count = sbs96_counts[sbs96]
        sbs52_counts[sbs52] += sbs96_count
    return sbs52_counts


def draw_sbs52_barplot(
    sbs52_file: str, 
    sample: str,
    sbs52_pdf_file: str
):

    df = pd.read_csv(sbs52_file, sep="\t")
    plot = (
        ggplot(df, aes(x="TRI", y="COUNT", fill="SUB"))
        + geom_bar(stat="identity")
        + theme_bw()
        + facet_grid(". ~ SUB", scales="free")
        + scale_fill_manual(
            values=("#98D7EC", "#212121", "#FF003A", "#A6A6A6", "#F5ABCC")
        )
        + labs(x="\nTrinucleotide Context\n", y="\nCounts\n")
        + ggtitle("\n{}\n".format(sample))
        + theme(
            text=element_text(size=10),
            legend_title=element_blank(),
            axis_text_x=element_text(
                family="monospace", angle=90, ha="center"
            ),
        )
    )
    plot.save(sbs52_pdf_file, width=22, height=12)


def dump_sbs52_counts(
    vcf_file, 
    ref_file,
    region,
    region_list,
    sbs52_file,
):

    sample = get_sample(vcf_file)
    chrom_lst = load_loci(region, region_list)
    sbs96_counts = load_sbs96_counts(vcf_file, ref_file, chrom_lst)
    sbs52_counts = load_sbs52_counts(sbs96_counts)
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
            ubase, _, ref, _, alt, _, dbase = list(sbs52)
            sub = f"{ref}>{alt}"
            tri = f"{ubase}{ref}{dbase}"
            print(
                sub,
                tri,
                sbs52,
                sbs52_counts[sbs52], 
                sep="\t",
                file=outfile
        )
    if sbs52_file.endswith(".tsv"):
        sbs52_pdf_file = sbs52_file.replace(".tsv", ".pdf")
    else:
        sbs52_pdf_file = sbs52_file + ".pdf" 
    draw_sbs52_barplot(sbs52_file, sample, sbs52_pdf_file)
       
        
def main():
    options = parse_args(sys.argv)
    dump_sbs52_counts(options.input, options.ref, options.region, options.region_list, options.output)
    sys.exit(0)


if __name__ == "__main__":
    main()
