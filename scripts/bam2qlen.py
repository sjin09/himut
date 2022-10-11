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
        "-t",
        "--threads",
        type=int,
        default=1,
        required=False,
        help="number of threads",
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


def load_chrom_lst(
    bam_file: str,
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
        alignments = pysam.AlignmentFile(bam_file, "rb")
        bam_header_lst = str(alignments.header).strip().split("\n")
        for h in bam_header_lst:
            if h.startswith("@SQ"):
                tname = h.split("\t")[1].replace("SN:", "")
                chrom_lst.append(tname)
        alignments.close()
    return chrom_lst 


def get_ccs_statistics(
    chrom: str,
    bam_file: str,
    chrom2qlen_lst: Dict[str, List[int]] 
):
  
    qlen_lst = [] 
    alignments = pysam.AlignmentFile(bam_file, "rb")
    for line in alignments.fetch(chrom):
        ccs = BAM(line)
        if not ccs.is_primary:
            continue
        qlen_lst.append((ccs.qname, ccs.qlen))
    chrom2qlen_lst[chrom] = qlen_lst


def dump_ccs_statistics(
    bam_file: str,
    region: str,
    region_file: str,
    threads: int,
    out_file: str
): 
    
    p = mp.Pool(threads)
    manager = mp.Manager()
    chrom2qlen_lst = manager.dict()
    chrom_lst = load_chrom_lst(bam_file, region, region_file)
    get_ccs_statistics_arg_lst = [
        (
            chrom,
            bam_file,
            chrom2qlen_lst
        )
        for chrom in chrom_lst
    ]
    p.starmap(
        get_ccs_statistics, get_ccs_statistics_arg_lst,
    )
    p.close()
    p.join()
  
    o = open(out_file, "w")
    o.write("{}\t{}\n".format("qname", "qlen"))
    for chrom in chrom_lst:
        for (qname, qlen) in chrom2qlen_lst[chrom]: 
            o.write("{}\t{}\n".format(qname, qlen))
    o.close() 


def main():
    options = parse_args(sys.argv)
    dump_ccs_statistics(
        options.bam, 
        options.region, 
        options.region_list,
        options.threads, 
        options.out
    )
    sys.exit(0)


if __name__ == "__main__":
    main()
