#!/usr/bin/env python

import sys
import pysam
import argparse
from pathlib import Path


def parse_args(args):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-i",
        "--vcf",
        type=str,
        required=True,
        help="VCF file to read"
    )
    args = args[1:]
    return parser.parse_args(args)


def get_samples(vcf_file: pysam.VariantFile):
    samples = vcf_file.header.samples
    return samples


def get_variant_per_sample(vcf_file_path: Path):
    vcf_file = pysam.VariantFile(vcf_file_path)
    samples = get_samples(vcf_file)
    header_lines = str(vcf_file.header).split("\n")[:-2]
    sample_files = {
        sample: open(f"{sample}.deepvariant.phased.vcf", 'w')
        for sample in samples
    }
    for sample in samples:
        last_line = f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}"
        for line in header_lines:
            sample_files[sample].write(f"{line}\n")
        sample_files[sample].write(f"{last_line}\n")
    for record in vcf_file:
        fields = str(record).rstrip().split("\t")
        constant_fields = fields[0:9]
        sample_records = fields[9:]
        for sample, sample_record in zip(samples, sample_records):
            sample_line = "\t".join(constant_fields + [sample_record])
            sample_files[sample].write(f"{sample_line}\n")
    for sample_file in sample_files.values():
        sample_file.close()
    vcf_file.close()


def main():
    options = parse_args(sys.argv)
    get_variant_per_sample(options.vcf)
    sys.exit(0)


if __name__ == "__main__":
    main()
