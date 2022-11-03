#!/usr/bin/env python3

import re
import os
import sys
import time
import pysam
import natsort
import argparse
import itertools
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
        help="VCF fofn (file of file names) to read",
    )
    parser.add_argument(
        "--region",
        type=str,
        required=False,
        help="target chromosome",
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
        help="file to write",
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


purine = set(["A", "G"])
pyrimidine = set(["T", "C"])
transitions = set(["A>G", "G>A", "C>T", "T>C"])
transversions = set(["A>C", "C>A", "C>G", "G>C", "A>T", "T>A", "G>T", "T>G"])


def load_loci(
    region: str, 
    region_list: str, 
) -> Tuple[List[str], List[Tuple[str, int, int]]]:
    
    chrom_lst = []
    if region is None and region_list is not None:
        for line in open(region_list).readlines():
            arr = line.strip().split()
            chrom_lst.append(arr[0])
    elif region is not None and region_list is None:
        chrom_lst.append(region)
    elif region is not None and region_list is not None:
        for line in open(region_list).readlines():
            arr = line.strip().split()
            chrom_lst.append(arr[0])
    else:
        pass
    return chrom_lst


def get_vcf_statistics(vcf_file: str, chrom_lst: List[str]):

    chrom2ti_count = defaultdict(lambda: 0)
    chrom2tv_count = defaultdict(lambda: 0)
    chrom2ins_count = defaultdict(lambda: 0)
    chrom2del_count = defaultdict(lambda: 0)
    chrom2hetsnp_count = defaultdict(lambda: 0)
    chrom2homsnp_count = defaultdict(lambda: 0)
    chrom2hetindel_count = defaultdict(lambda: 0)
    chrom2homindel_count = defaultdict(lambda: 0)
    if vcf_file.endswith(".vcf"):
        for line in open(vcf_file).readlines():
            if line.startswith("#"):
                continue
            v = VCF(line)
            if v.is_pass and v.is_biallelic:
                if v.is_snp:
                    if v.sample_gt == "0/1":
                        chrom2hetsnp_count[v.chrom] += 1
                    elif v.sample_gt == "1/1":
                        chrom2homsnp_count[v.chrom] += 1
                    sub = "{}>{}".format(v.ref, v.alt)
                    if sub in transitions:
                        chrom2ti_count[v.chrom] += 1
                    elif sub in transversions:
                        chrom2tv_count[v.chrom] += 1
                elif v.is_indel:
                    if len(v.ref) > len(v.alt):
                        chrom2del_count[v.chrom] += 1
                    elif len(v.ref) < len(v.alt):
                        chrom2ins_count[v.chrom] += 1 
                    if v.sample_gt == "0/1":
                        chrom2hetindel_count[v.chrom] += 1
                    elif v.sample_gt == "1/1":
                        chrom2homindel_count[v.chrom] += 1
    ti_count = 0
    tv_count = 0
    ins_count = 0
    del_count = 0
    hetsnp_count = 0
    homsnp_count = 0
    hetindel_count = 0
    homindel_count = 0
    for chrom in chrom_lst:
        ti_count += chrom2ti_count[chrom] 
        tv_count += chrom2tv_count[chrom] 
        ins_count += chrom2ins_count[chrom] 
        del_count += chrom2del_count[chrom] 
        hetsnp_count += chrom2hetsnp_count[chrom] 
        homsnp_count += chrom2homsnp_count[chrom] 
        hetindel_count += chrom2hetindel_count[chrom] 
        homindel_count += chrom2homindel_count[chrom] 
        
    titv_ratio = ti_count/float(tv_count)
    snp_hethom_ratio = hetsnp_count/float(homsnp_count)
    indel_hethom_ratio = hetindel_count/float(homindel_count)
    indel_ratio = ins_count/float(del_count)
    return titv_ratio, hetsnp_count, homsnp_count, snp_hethom_ratio, hetindel_count, homindel_count, indel_hethom_ratio, indel_ratio


def dump_vcf_statistics(infile: str, region: str, region_list: str, outfile: str):

    o = open(outfile, "w")
    o.write(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
            "vcf_file",
            "titv_ratio",
            "hetsnp",
            "homsnp",
            "snp het/hom ratio",
            "hetindel",
            "homindel",
            "indel het/hom ratio",
            "indel ratio",
        )
    )
    chrom_lst = load_loci(region, region_list)
    vcf_file_lst = natsort.natsorted(
        [line.strip() for line in open(infile).readlines()]
    )
 
    for vcf_file in vcf_file_lst:
        (
            titv_ratio,
            hetsnp_count,
            homsnp_count,
            snp_hethom_ratio,
            hetindel_count,
            homindel_count,
            indel_hethom_ratio,
            indel_ratio,
        ) = get_vcf_statistics(vcf_file, chrom_lst)
        o.write(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                vcf_file,
                titv_ratio,
                hetsnp_count,
                homsnp_count,
                snp_hethom_ratio,
                hetindel_count,
                homindel_count,
                indel_hethom_ratio,
                indel_ratio,
            )
        )


def main():
    options = parse_args(sys.argv)
    dump_vcf_statistics(options.input, options.region, options.region_list, options.output)
    sys.exit(0)


if __name__ == "__main__":
    main()
