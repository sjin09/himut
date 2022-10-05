#!/usr/bin/env python3

import re
import sys
import pysam
import tabix
import cyvcf2
import natsort
import argparse
import multiprocessing as mp
from collections import defaultdict
from typing import Dict, List, Tuple, Set


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
        "--vcf",
        type=str,
        required=True,
        help="deepvariant VCF file to read",
    )
    parser.add_argument(
        "--min_bq",
        type=int,
        default=93,
        required=False,
        help="minimum base quality score",
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=1,
        required=False,
        help="maximum number of threads",
    )
    parser.add_argument(
        "-o",
        "--out",
        type=str,
        required=True,
        help="file to write",
    )
    args = args[1:]
    return parser.parse_args(args)


class BAM:
    def __init__(self, line):
        # target
        self.tname = line.reference_name
        self.tstart = line.reference_start
        self.tend = line.reference_end
        self.tcoord = "{}:{}-{}".format(self.tname, self.tstart, self.tend)
        # query
        self.qname = line.query_name.replace("/ccs", "")
        self.qstart = line.query_alignment_start
        self.qend = line.query_alignment_end
        self.qseq = line.query_sequence
        self.qlen = len(self.qseq)
        self.mapq = line.mapping_quality
        self.bq_int_lst = line.query_qualities
        self.qv = sum(self.bq_int_lst)/self.qlen
        self.hbq_proportion = self.bq_int_lst.count(93)/float(self.qlen)
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
        self.format_lst = arr[8].split(":")
        self.sample_format_lst = arr[9].split(":")
        hsh = {i:j for i,j in zip(self.format_lst, self.sample_format_lst)}
        if "GT" in hsh: self.sample_gt = hsh["GT"] 
        if "PS" in hsh: self.sample_phase_set = hsh["PS"] 
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
            elif len(self.ref) > len(self.alt): # del
                self.is_indel = True
            elif len(self.ref) < len(self.alt): # ins
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


def get_tname2tsize(bam_file: str) -> Tuple[List[str], Dict[str, int]]:

    tname2tsize = {}
    alignments = pysam.AlignmentFile(bam_file, "rb")
    bam_header_lst = str(alignments.header).strip().split("\n")
    for h in bam_header_lst:
        if h.startswith("@SQ"):
            _tag, tname, tsize = h.split("\t")
            tname = tname.replace("SN:", "")
            tsize = tsize.replace("LN:", "")
            tname2tsize[tname] = int(tsize)
    alignments.close()
    return tname2tsize


def load_vcf_file(
    chrom: str,
    chrom_len: int,
    vcf_file: str
) -> Set[Tuple[int, str, str]]:

    snp_set = set()
    if vcf_file.endswith(".vcf"):
        for line in open(vcf_file):
            if line.startswith("#"):
                continue
            v = VCF(line)
            if v.chrom != chrom:
                continue
            if v.is_pass:
                if v.is_biallelic:
                    if v.is_snp:
                        snp_set.add((v.pos, v.ref, v.alt))
                else:
                    for alt in v.alt_lst:
                        if len(v.ref) == 1 and len(alt) == 1:
                            snp_set.add((v.pos, v.ref, alt))
    elif vcf_file.endswith(".bgz"):
        tb = tabix.open(vcf_file)
        records = tb.query(chrom, 0, chrom_len)
        for record in records:
            v = VCF("\t".join(record))            
            if v.is_snp and v.is_pass:
                if v.is_biallelic:
                    if v.is_snp:
                        snp_set.add((v.pos, v.ref, v.alt))
                else:
                    for alt in v.alt_lst:
                        if len(v.ref) == 1 and len(alt) == 1:
                            snp_set.add((v.pos, v.ref, alt))
    return snp_set

def get_ccs_statistics(
    chrom: str,
    chrom_len: int,
    bam_file: str,
    vcf_file: str,
    min_bq: int,
    chrom2ccs_lst: Dict[str, str] 
):
  
    ccs_lst = [] 
    alignments = pysam.AlignmentFile(bam_file, "rb")
    deepvariant_set = load_vcf_file(chrom, chrom_len, vcf_file)
    for line in alignments.fetch(chrom):
        ccs = BAM(line)
        if not ccs.is_primary:
            continue
        
        if ccs.mapq > 0:
            ccs.cs2subindel()
            deepvariant_snp_lst = []                    
            himut_candidate_lst = []
            himut_q93_candidate_lst = []
            for tsbs, qbq in zip(ccs.tsbs_lst, ccs.qsbs_bq_lst): 
                if tsbs in deepvariant_set:
                    deepvariant_snp_lst.append(tsbs)
                else:
                    himut_candidate_lst.append(tsbs)
                    if qbq >= min_bq:
                        himut_q93_candidate_lst.append(tsbs)
            ccs_lst.append(
                "\t".join(
                    [
                        ccs.qname, 
                        ccs.tcoord, 
                        str(ccs.qv), 
                        str(ccs.qlen), 
                        str(ccs.mapq), 
                        str(ccs.hbq_proportion), 
                        str(ccs.query_alignment_proportion), 
                        str(len(himut_q93_candidate_lst)), 
                        str(len(himut_candidate_lst)),
                        str(len(deepvariant_snp_lst)),
                        str(len(ccs.tsbs_lst))
                    ]
                )
            )
    chrom2ccs_lst[chrom] = ccs_lst

def dump_ccs_statistics(
    bam_file: str,
    vcf_file: str, 
    region: str,
    region_file: str,
    min_bq: int,
    threads: int,
    out_file: str
): 
    
    tname2tsize = get_tname2tsize(bam_file)
    chrom_lst = load_chrom_lst(region, region_file)

    p = mp.Pool(threads)
    manager = mp.Manager()
    chrom2ccs_lst = manager.dict()
    get_ccs_statistics_arg_lst = [
        (
            chrom,
            tname2tsize[chrom],
            bam_file,
            vcf_file,
            min_bq,
            chrom2ccs_lst
        )
        for chrom in chrom_lst
    ]
    p.starmap(
        get_ccs_statistics, get_ccs_statistics_arg_lst,
    )
    p.close()
    p.join()
  
    o = open(out_file, "w")
    o.write("{}\n".format("\t".join(["qname", "tcoord", "qv", "qlen", "mapq", "hbq_proportion", "alignment_proportion", "himut_q93_candidate_count", "himut_candidate_count", "snp_count", "sbs_count"])))
    for chrom in chrom_lst:
        for p in chrom2ccs_lst[chrom]: 
            o.write("{}\n".format(p))
    o.close() 


def main():
    options = parse_args(sys.argv)
    dump_ccs_statistics(
        options.bam, 
        options.vcf,
        options.region, 
        options.region_list,
        options.min_bq,
        options.threads, 
        options.out
    )
    sys.exit(0)


if __name__ == "__main__":
    main()
