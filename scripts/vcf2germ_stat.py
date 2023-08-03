#!/usr/bin/env python3

import sys
import tabix
import cyvcf2
import natsort
import argparse
import numpy as np 
import himut.vcflib
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
        help="VCF file to read",
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
        print("--region or --region_list parameter is required")
        sys.exit()
    return chrom_lst


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
        v = cyvcf2.VCF(vcf_file)
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


def get_vcf_statistics(
    vcf_file: str, 
    chrom_lst: List[str],
    tname2tsize: Dict[str, int],
):

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
            v = himut.vcflib.VCF(line)
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
    else:
        tb = tabix.open(vcf_file)
        for chrom in chrom_lst:
            chrom_len = tname2tsize[chrom]
            records = tb.query(chrom, 0, chrom_len)
            for record in records:
                v = himut.vcflib.VCF("\t".join(record))
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


def dump_vcf_statistics(vcf_file: str, region: str, region_list: str, outfile: str):


    chrom_lst = load_loci(region, region_list)
    tname2tsize = get_tname2tsize(vcf_file)
    o = open(outfile, "w")
    o.write(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
            "titv_ratio",
            "het_snp",
            "hom_snp",
            "snp_het/hom_ratio",
            "het_indel",
            "hom_indel",
            "indel_het/hom_ratio",
            "indel_ratio",
        )
    )
    
    (
        titv_ratio,
        hetsnp_count,
        homsnp_count,
        snp_hethom_ratio,
        hetindel_count,
        homindel_count,
        indel_hethom_ratio,
        indel_ratio,
    ) = get_vcf_statistics(vcf_file, chrom_lst, tname2tsize)
    o.write(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
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
