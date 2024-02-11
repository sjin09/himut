#!/usr/bin/env python

import argparse
import gzip
from pathlib import Path
import sys


def parse_args(args):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        required=True,
        help="gff3 gene annotation"
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="file to return transcripts for identification of mutational specturms with transcriptional-strand bias"
    )
    args = args[1:]
    return parser.parse_args(args)


class GFF3:
    def __init__(self, line):
        arr = line.strip().split()
        self.chrom = arr[0]
        self.source = arr[1]
        self.type = arr[2]
        self.start = arr[3]
        self.end = arr[4]
        self.score = arr[5]
        self.strand = 1 if arr[6] == "+" else "-1"
        self.phase = arr[7]
        self.attribute_arr = arr[8].split(";")


def gff2transcripts(gff_file_path, transcripts_file_path):
    with open(transcripts_file_path, "w") as outfile:
        for line in gzip.open(gff_file_path, "rt").readlines():
            if line.startswith("#"):
                continue
            gff = GFF3(line)
            if gff.type == "gene":
                continue
            elif gff.type == "mRNA":
                for attribute in gff.attribute_arr:
                    if attribute.startswith("ID"):
                        transcript = attribute.split("=")[1].split(":")[1]
                    elif attribute.startswith("Parent"):
                        gene = attribute.split("=")[1].split(":")[1]
                    else:
                        continue
                if gff.chrom.isdigit():
                    print(
                        gene,
                        transcript,
                        gff.chrom,
                        gff.strand,
                        gff.start,
                        gff.end,
                        gene,
                        "protein_coding",
                        sep="\t",
                        file=outfile
                    )
            elif gff.type == "three_prime_UTR":
                continue
            elif gff.type == "exon":
                continue
            elif gff.type == "CDS":
                continue
            elif gff.type == "five_prime_UTR":
                continue


def main():
    options = parse_args(sys.argv)
    gff2transcripts(options.input, options.output)
    sys.exit(0)


if __name__ == "__main__":
    main()
