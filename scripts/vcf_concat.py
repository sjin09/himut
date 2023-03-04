#!/usr/bin/env python3

import sys
import natsort
import argparse
import himut.vcflib
from collections import defaultdict

def parse_args(args):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="file containing path to VCF files separated by new line",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="file to write",
    )
    args = args[1:]
    return parser.parse_args(args)


def vcf_concat(fofn: str, out_file: str) -> None:

    vcf_header_lst = []
    vcf_lst = [i.strip() for i in open(fofn).readlines()]
    for line in open(vcf_lst[0]):
        if line.startswith("##"):
            vcf_header_lst.append(line.strip())

    mut_set = set()
    sample_lst = []
    sample2mut_set = defaultdict(set)
    sample2mut2format = defaultdict(dict)
    for vcf in vcf_lst:
        for line in open(vcf).readlines():
            if line.startswith("##"):
                continue
            elif line.startswith("#CHROM"):
                sample = line.split()[-1]
                sample_lst.append(sample)
            else:
                v = himut.vcflib.VCF(line)
                if v.is_pass:
                    mut = (v.chrom, v.pos, v.ref, v.alt)
                    mut_set.add(mut)
                    sample2mut_set[sample].add(mut)
                    sample2mut2format[sample][mut] = ":".join(v.sample_format_lst)
    sample_lst = natsort.natsorted(sample_lst)
    mut_lst = natsort.natsorted(list(mut_set))

    o = open(out_file, "w") 
    vcf_header_lst.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}".format("\t".join(sample_lst)))
    o.write("{}\n".format("\n".join(vcf_header_lst)))
    for mut in mut_lst:
        chrom, pos, ref, alt = mut
        sample_mut_lst = [sample2mut2format[sample].get(mut, ".:.:.:.:.:.") for sample in sample_lst]
        o.write("{}\t{}\t.\t{}\t{}\t.\tPASS\t.\tGT:BQ:DP:AD:VAF:PS\t{}\n".format(chrom, pos, ref, alt, "\t".join(sample_mut_lst)))
    o.close()


def main():
    options = parse_args(sys.argv)
    vcf_concat(options.input, options.output)
    sys.exit(0)


if __name__ == "__main__":
    main()
