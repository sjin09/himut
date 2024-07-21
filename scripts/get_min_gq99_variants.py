#!/usr/bin/env python

import argparse
import sys

from pathlib import Path
import pysam


MIN_GQ = 99


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
        # return str(list(self.variant.filter)[0])

    @property
    def samples(self) -> str:
        return self.variant.samples

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
        "-o",
        "--output",
        type=Path,
        required=True,
        help="VCF file to read"
    )
    args = args[1:]
    return parser.parse_args(args)


def get_samples(vcf_header: str):
    samples = vcf_header.rstrip().split("\n")[-1].split()[9:]
    return samples


def get_min_gq99_variants(input_path: Path, output_path: Path):
    with open(output_path, "w") as outfile:
        vcf_file = pysam.VariantFile(input_path)
        vcf_header = str(vcf_file.header)
        outfile.write("{}".format(vcf_header))
        samples = get_samples(vcf_header)
        sample = samples[0]
        for variant in vcf_file:
            xvariant = ExpandedVariantRecord(variant)
            if not xvariant.is_pass:
                continue
            if not xvariant.is_biallelic_snp:
                continue
            gq = int(variant.samples[sample]["GQ"])
            if gq >= MIN_GQ:
                outfile.write("{}".format(str(variant)))


def main():
    options = parse_args(sys.argv)
    get_min_gq99_variants(options.input, options.output)
    sys.exit(0)


if __name__ == "__main__":
    main()
