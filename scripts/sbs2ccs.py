#!/usr/bin/env python3

import re
import sys
import time
import pysam
import cyvcf2
import argparse
import multiprocessing as mp
from collections import defaultdict
from typing import Dict, List, Tuple


def parse_args(args):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-i",
        "--bam",
        type=str,
        required=True,
        help="BAM filed to read",
    )
    parser.add_argument(
        "--vcf",
        type=str,
        required=True,
        help="single base substitution VCF file to read",
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
        "--min_mapq",
        type=int,
        default=60,
        required=False,
        help="minimum mapping quality score",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        required=False,
        help="number of threads"
    )
    parser.add_argument(
        "-o",
        "--out",
        type=str,
        required=True,
        help="file to return",
    )
    args = args[1:]
    return parser.parse_args(args)


class BAM:
    def __init__(self, line):
        # target
        self.tname = line.reference_name
        self.tstart = line.reference_start
        self.tend = line.reference_end
        self.target_alignment_length = self.tend - self.tstart
        # query
        self.qname = line.query_name
        self.qstart = line.query_alignment_start
        self.qend = line.query_alignment_end
        self.qseq = line.query_sequence
        self.qlen = len(self.qseq)
        self.mapq = line.mapping_quality
        self.bq_int_lst = line.query_qualities
        self.query_alignment_length = self.qend - self.qstart
        self.query_alignment_proportion = self.query_alignment_length/float(self.qlen)
        self.cs_tag = line.get_tag("cs") if line.has_tag("cs") else "."
        if line.has_tag("tp"):
            if line.get_tag("tp") == "P":
                self.is_primary = True
            else:
                self.is_primary = False
        else:
            self.is_primary = False

    def cs2lst(self):
        cslst = [cs for cs in re.split("(:[0-9]+|\*[a-z][a-z]|[=\+\-][A-Za-z]+)", self.cs_tag)]
        cslst = [cs.upper() for cs in cslst if cs != ""]
        return cslst
            
    def cs2tuple(self) -> List[Tuple[int, str, str, int, int]]:
        qpos = self.qstart
        self.cstuple_lst = []
        cs_lst = self.cs2lst()
        for cs in cs_lst:
            m = cs[1:]
            mlen = len(m)
            qstart = qpos
            if cs.startswith("="):  # match # --cs=long
                cs = ":{}".format(mlen)
                t = (1, m, m, mlen, mlen)
            elif cs.startswith(":"):  # match # --cs=short
                mlen = int(m)
                qend = qpos + mlen
                m = self.qseq[qstart:qend]
                t = (1, m, m, mlen, mlen)
            elif cs.startswith("*"):  # snp # target and query
                mlen = 1
                ref, alt = list(m)
                t = (2, ref, alt, 1, 1)
            elif cs.startswith("+"):  # insertion # query
                ref = self.qseq[qpos - 1]
                alt = ref + m
                t = (3, ref, alt, 0, mlen)
            elif cs.startswith("-"):  # deletion # target
                alt = self.qseq[qpos - 1]
                ref = alt + m
                t = (4, ref, alt, mlen, 0)
                mlen = 0
            self.cstuple_lst.append(t)
            qpos += mlen

    def cs2subindel(self):
        self.cs2tuple()
        tpos = self.tstart
        qpos = self.qstart
        self.tsbs_lst = []
        self.qsbs_lst = []
        self.qsbs_bq_lst = []
        self.tindel_lst = []
        for cstuple in self.cstuple_lst:
            state, ref, alt, ref_len, alt_len, = cstuple
            if state == 2:  # snp 
                if ref == "N": 
                    continue
                self.qsbs_bq_lst.append(self.bq_int_lst[qpos])
                self.tsbs_lst.append((tpos + 1, ref, alt))
                self.qsbs_lst.append((qpos + 1, alt, ref))
            elif state == 3 or state == 4:  # insertion
                self.tindel_lst.append((tpos, ref, alt))
            tpos += ref_len 
            qpos += alt_len

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


def load_chrom_lst(
    region: str, 
    region_file: str, 
) -> Tuple[List[str], List[Tuple[str, int, int]]]:

    chrom_lst = []
    if region is None and region_file is not None:
        for line in open(region_file).readlines():
            chrom_lst.append(line.strip().split()[0])
    elif region is not None and region_file is None:
        chrom_lst = [region]
    elif region is not None and region_file is not None:
        for line in open(region_file).readlines():
            chrom_lst.append(line.strip().split()[0])
    else:
        print("--region or --region_list parameter is required")
        sys.exit()
    return chrom_lst


def load_vcf_file(vcf_file: str):

    chrom2sbs_lst = defaultdict(list)
    if vcf_file.endswith(".vcf"):
        for line in open(vcf_file).readlines():
            if line.startswith("#"):
                continue
            v = VCF(line)
            if v.is_pass:
                chrom2sbs_lst[v.chrom].append((v.pos, v.ref, v.alt))
    elif vcf_file.endswith(".bgz"):
        for i in cyvcf2.VCF(vcf_file):
            v = VCF(str(i))
            if v.is_pass:
                chrom2sbs_lst[v.chrom].append((v.pos, v.ref, v.alt))
    return chrom2sbs_lst


def get_ccs(
    chrom: str,
    sbs_lst: List[Tuple[int, int, str]],
    bam_file: str,
    min_mapq: int,
    chrom2sbs2ccs: Dict[str, Dict[str, Tuple[str, int, int, int, int]]]
):
    sbs2ccs = {}
    alignments = pysam.AlignmentFile(bam_file, "rb")
    for sbs in sbs_lst:
        pos = sbs[0]
        for line in alignments.fetch(chrom, pos, pos + 1):
            ccs = BAM(line)
            if ccs.mapq >= min_mapq and ccs.is_primary:
                ccs.cs2subindel()
                tsbs_set = set(ccs.tsbs_lst)
                if sbs in tsbs_set:
                    sbs2ccs[sbs] = [ccs.qname, ccs.qlen, ccs.mapq, len(ccs.tsbs_lst), len(ccs.tindel_lst)]
    chrom2sbs2ccs[chrom] = sbs2ccs

def sbs2ccs(
    bam_file: str, 
    vcf_file: str,
    region: str,
    region_file: str,
    min_mapq: int,
    threads: int,
    out_file: str
): 

    chrom_lst = load_chrom_lst(region, region_file)
    chrom2sbs_lst = load_vcf_file(vcf_file)
    p = mp.Pool(threads)
    manager = mp.Manager()
    chrom2sbs2ccs = manager.dict()
    get_ccs_arg_lst = [
        (
            chrom,
            chrom2sbs_lst[chrom],
            bam_file,
            min_mapq,
            chrom2sbs2ccs,
        )
        for chrom in chrom_lst
    ]
    print("sbs2ccs started retrieving CCS read for each single base substitution")
    p.starmap(
        get_ccs, get_ccs_arg_lst,
    )
    p.close()
    p.join()
    print("sbs2ccs finished retrieving CCS read for each single base substitution")
   
    # def dump_sbs2ccs()
    o = open(out_file, "w") 
    o.write("{}\n".format("\t".join(["sbs", "ccs", "ccs_len", "mapq", "sbs_cnt", "indel_cnt", "cell", "zmw"])))
    for chrom in chrom2sbs2ccs:
        for (pos, ref, alt), (qname, qlen, mapq, tsbs_cnt, tindel_cnt) in chrom2sbs2ccs[chrom].items():
            cell, zmw, _ = qname.split("/")
            o.write("{}:{}_{}/{}\t{}\t{}\t60\t{}\t{}\t{}\t{}\n".format(chrom, pos, ref, alt, qname, qlen, tsbs_cnt, tindel_cnt, cell, zmw))
    o.close() 


def main():
    options = parse_args(sys.argv)
    sbs2ccs(
        options.bam, 
        options.vcf, 
        options.region, 
        options.region_list, 
        options.min_mapq, 
        options.threads,
        options.out
    )
    sys.exit(0)


if __name__ == "__main__":
    main()
