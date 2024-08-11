#!/usr/bin/env python

import argparse
import sys
from pathlib import Path
from typing import List

import pysam


TRANSITIONS = set(["A>G", "G>A", "C>T", "T>C"])
TRANSVERSIONS = set(["A>C", "C>A", "C>G", "G>C", "A>T", "T>A", "G>T", "T>G"])


def parse_args(args):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--vcf",
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
        "--tgt",
        type=str,
        required=True,
        help="list of chromosomes separated by new line"
    )
    parser.add_argument(
        "--hymenoptera-sample",
        required=False,
        action="store_true",
        help="hymenoptera sample"
    )
    parser.add_argument(
        "--sample-is-reference-sample",
        required=False,
        action="store_true",
        help="sample is reference sample"
    )
    parser.add_argument(
        "-o",
        "--out",
        type=str,
        required=True,
        help="file to write"
    )
    args = args[1:]
    return parser.parse_args(args)


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

    @property
    def format(self) -> str:
        return self.variant.format

    @property
    def is_snp(self) -> bool:
        return True if len(self.ref) == 1 and len(self.alt) == 1 else False

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


def load_chromosomes(tgt_path: Path):
    chroms = [line.strip() for line in open(tgt_path)]
    return chroms


def load_mutation_counts(chroms: List[str], vcf_path: Path, sample: str):
    ti_count = 0
    tv_count = 0
    snp_count = 0
    del_count = 0
    ins_count = 0
    variant_records = pysam.VariantFile(vcf_path)
    for chrom in chroms:
        for variant in variant_records.fetch(chrom):
            xvariant = ExpandedVariantRecord(variant)
            if not xvariant.is_pass:
                continue
            if not xvariant.is_biallelic:
                continue
            if xvariant.is_snp:
                snp_count += 1
                sub = "{}>{}".format(xvariant.ref, xvariant.alt)
                if sub in TRANSITIONS:
                    ti_count += 1
                elif sub in TRANSVERSIONS:
                    tv_count += 1
            else:
                if len(xvariant.ref) > len(xvariant.alt):  # deletion
                    del_count += 1
                elif len(xvariant.ref) < len(xvariant.alt):  # insertion
                    ins_count += 1
    titv_ratio = ti_count/float(tv_count) if tv_count != 0 else 0.0
    indel_count = del_count + ins_count
    indel_ratio = ins_count/float(del_count) if del_count != 0 else 0.0
    return snp_count, ti_count, tv_count, titv_ratio, del_count, ins_count, indel_count, indel_ratio


def load_het_mutation_counts(chroms: List[str], vcf_path: Path, sample: str):

    ti_count = 0
    tv_count = 0
    snp_count = 0
    del_count = 0
    ins_count = 0
    variant_records = pysam.VariantFile(vcf_path)
    for chrom in chroms:
        for variant in variant_records.fetch(chrom):
            xvariant = ExpandedVariantRecord(variant)
            if not xvariant.is_pass:
                continue
            sample_gt = xvariant.format[sample]["GT"]
            if sample_gt != "0/1":
                continue
            if xvariant.is_snp:
                snp_count += 1
                sub = "{}>{}".format(xvariant.ref, xvariant.alt)
                if sub in TRANSITIONS:
                    ti_count += 1
                elif sub in TRANSVERSIONS:
                    tv_count += 1
            else:
                if len(xvariant.ref) > len(xvariant.alt):  # deletion
                    del_count += 1
                elif len(xvariant.ref) < len(xvariant.alt):  # insertion
                    ins_count += 1
    titv_ratio = ti_count/float(tv_count) if tv_count != 0 else 0.0
    indel_count = del_count + ins_count
    indel_ratio = ins_count/float(del_count) if del_count != 0 else 0.0
    return snp_count, ti_count, tv_count, titv_ratio, del_count, ins_count, indel_count, indel_ratio


def load_hom_mutation_counts(chroms: List[str], vcf_path: Path, sample: str):
    ti_count = 0
    tv_count = 0
    snp_count = 0
    del_count = 0
    ins_count = 0
    variant_records = pysam.VariantFile(vcf_path)
    for chrom in chroms:
        for variant in variant_records.fetch(chrom):
            xvariant = ExpandedVariantRecord(variant)
            if not xvariant.is_pass:
                continue
            sample_gt = xvariant.format[sample]["GT"]
            if sample_gt != "0/1":
                continue
            if xvariant.is_snp:
                snp_count += 1
                sub = "{}>{}".format(xvariant.ref, xvariant.alt)
                if sub in TRANSITIONS:
                    ti_count += 1
                elif sub in TRANSVERSIONS:
                    tv_count += 1
            else:
                if len(xvariant.ref) > len(xvariant.alt):  # deletion
                    del_count += 1
                elif len(xvariant.ref) < len(xvariant.alt):  # insertion
                    ins_count += 1
    titv_ratio = ti_count/float(tv_count) if tv_count != 0 else 0.0
    indel_count = del_count + ins_count
    indel_ratio = ins_count/float(del_count) if del_count != 0 else 0.0
    return snp_count, ti_count, tv_count, titv_ratio, del_count, ins_count, indel_count, indel_ratio


def write_germline_mutation_statistics(
    vcf_path: Path,
    ref_path: Path,
    tgt_path: Path,
    hymenoptera_sample: bool,
    sample_is_reference_sample: bool,
    out_path: Path
):
    reference_sequence_lookup = pysam.FastaFile(ref_path)
    chroms = load_chromosomes(tgt_path)
    genome_sum = sum([len(reference_sequence_lookup[chrom]) for chrom in chroms])
    variant_records = pysam.VariantFile(vcf_path)
    sample = variant_records.header.samples[0]
    with open(out_path, "w") as outfile:
        # outfile.write("SAMPLE\tSNP_COUNT\tSNP_DENSIY\tTI_COUNT\tTV_COUNT\tTITV\tDEL_COUNT\tDEL_DENSITY\tINS_COUNT\tINS_DENSITY\tINDEL_COUNT\tINDEL_DENSITY\n")
        if hymenoptera_sample and sample_is_reference_sample:
            (
                het_snp_count,
                het_ti_count,
                het_tv_count,
                het_titv_ratio,
                het_del_count,
                het_ins_count,
                het_indel_count,
                het_indel_ratio
            ) = load_het_mutation_counts(chroms, vcf_path, sample)
            (
                hom_snp_count,
                hom_ti_count,
                hom_tv_count,
                hom_titv_ratio,
                hom_del_count,
                hom_ins_count,
                hom_indel_count,
                hom_indel_ratio
            ) = load_hom_mutation_counts(chroms, vcf_path, sample)
            het_snp_density = het_snp_count/genome_sum
            het_del_density = het_del_count/genome_sum
            het_ins_density = het_ins_count/genome_sum
            het_indel_density = het_indel_count/genome_sum
            hom_snp_density = hom_snp_count/genome_sum
            hom_del_density = hom_del_count/genome_sum
            hom_ins_density = hom_ins_count/genome_sum
            hom_indel_density = hom_indel_count/genome_sum
            outfile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                sample,
                het_snp_count,
                het_snp_density,
                het_ti_count,
                het_tv_count,
                het_titv_ratio,
                het_del_count,
                het_del_density,
                het_ins_count,
                het_ins_density,
                het_indel_count,
                het_indel_density,
                het_indel_ratio
            ))
            outfile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                sample,
                hom_snp_count,
                hom_snp_density,
                hom_ti_count,
                hom_tv_count,
                hom_titv_ratio,
                hom_del_count,
                hom_del_density,
                hom_ins_count,
                hom_ins_density,
                hom_indel_count,
                hom_indel_density,
                hom_indel_ratio
            ))
        elif not hymenoptera_sample and sample_is_reference_sample:
            (
                het_snp_count,
                het_ti_count,
                het_tv_count,
                het_titv_ratio,
                het_del_count,
                het_ins_count,
                het_indel_count,
                het_indel_ratio
            ) = load_het_mutation_counts(chroms, vcf_path, sample)
            het_snp_density = het_snp_count/genome_sum
            het_del_density = het_del_count/genome_sum
            het_ins_density = het_ins_count/genome_sum
            het_indel_density = het_indel_count/genome_sum
            outfile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                sample,
                het_snp_count,
                het_snp_density,
                het_ti_count,
                het_tv_count,
                het_titv_ratio,
                het_del_count,
                het_del_density,
                het_ins_count,
                het_ins_density,
                het_indel_count,
                het_indel_density,
                het_indel_ratio
            ))
        elif hymenoptera_sample and not sample_is_reference_sample:
            (
                het_snp_count,
                het_ti_count,
                het_tv_count,
                het_titv_ratio,
                het_del_count,
                het_ins_count,
                het_indel_count,
                het_indel_ratio
            ) = load_het_mutation_counts(chroms, vcf_path, sample)
            (
                hom_snp_count,
                hom_ti_count,
                hom_tv_count,
                hom_titv_ratio,
                hom_del_count,
                hom_ins_count,
                hom_indel_count,
                hom_indel_ratio
            ) = load_hom_mutation_counts(chroms, vcf_path, sample)
            het_snp_density = het_snp_count/genome_sum
            het_del_density = het_del_count/genome_sum
            het_ins_density = het_ins_count/genome_sum
            het_indel_density = het_indel_count/genome_sum
            hom_snp_density = hom_snp_count/genome_sum
            hom_del_density = hom_del_count/genome_sum
            hom_ins_density = hom_ins_count/genome_sum
            hom_indel_density = hom_indel_count/genome_sum
            outfile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                sample,
                het_snp_count,
                het_snp_density,
                het_ti_count,
                het_tv_count,
                het_titv_ratio,
                het_del_count,
                het_del_density,
                het_ins_count,
                het_ins_density,
                het_indel_count,
                het_indel_density,
                het_indel_ratio
            ))
            outfile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                sample,
                hom_snp_count,
                hom_snp_density,
                hom_ti_count,
                hom_tv_count,
                hom_titv_ratio,
                hom_del_count,
                hom_del_density,
                hom_ins_count,
                hom_ins_density,
                hom_indel_count,
                hom_indel_density,
                hom_indel_ratio
            ))
        elif not hymenoptera_sample and not sample_is_reference_sample:
            (
                snp_count,
                ti_count,
                tv_count,
                titv_ratio,
                del_count,
                ins_count,
                indel_count,
                indel_ratio
            ) = load_mutation_counts(chroms, vcf_path, sample)
            snp_density = snp_count/genome_sum
            del_density = del_count/genome_sum
            ins_density = ins_count/genome_sum
            indel_density = indel_count/genome_sum
            outfile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                sample,
                snp_count,
                snp_density,
                ti_count,
                tv_count,
                titv_ratio,
                del_count,
                del_density,
                ins_count,
                ins_density,
                indel_count,
                indel_density,
                indel_ratio
            ))


def main():
    options = parse_args(sys.argv)
    write_germline_mutation_statistics(
        options.vcf,
        options.ref,
        options.tgt,
        options.hymenoptera_sample,
        options.sample_is_reference_sample,
        options.out
    )
    sys.exit(0)


if __name__ == "__main__":
    main()
