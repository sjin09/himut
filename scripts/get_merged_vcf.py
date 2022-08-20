#!/usr/bin/env python3

import re
import os
import sys
import time
import pysam
import natsort
import argparse
from collections import defaultdict
from typing import Set, Dict, List, Tuple


def parse_args(args):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="path to VCF files separated by new line",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="file to return",
    )
    args = args[1:]
    return parser.parse_args(args)


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
        self.format = arr[8]
        self.sample_format = arr[9]
        self.sample_format_lst = self.sample_format.split(":")
        self.sample_gt = self.sample_format_lst[0]
        self.is_snp = False
        self.is_indel = False
        if len(self.alt_lst) == 1:  # bi-allelic
            self.is_biallelic = True
            self.alt = self.alt_lst[0]
            if len(self.ref) == 1 and len(self.alt) == 1:  # snp
                self.is_snp = True
            else:  # indel
                self.is_indel = True
        else:
            self.is_biallelic = False


def get_merged_vcf(infile: str, outfile: str) -> None:
    """
    merge VCF file from each chromosome from the same sample to form one VCF file 
    """

    o = open(outfile, "w") 
    vcf_lst = natsort.natsorted([line.strip() for line in open(infile).readlines()])
    for line in open(vcf_lst[0]).readlines():
        if line.startswith("#"):
            o.write("{}\n".format(line.strip()))
        else:
            break

    for vcf in vcf_lst:
        for line in open(vcf).readlines():
            if line.startswith("#"):
                continue
            else:
                o.write("{}\n".format(line.strip()))

def main():
    options = parse_args(sys.argv)
    get_merged_vcf(options.input, options.output)
    sys.exit(0)


if __name__ == "__main__":
    main()
