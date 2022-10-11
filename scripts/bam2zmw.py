#!/usr/bin/env python

import sys
import pysam
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
        "--input",
        type=str,
        required=True,
        help="BAM file to read"
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
        required=False,
        help="number of threads"
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

def get_qname(
    chrom: str,
    bam_file: str,
    chrom2qname_lst: Dict[str, List[str]]
):

    qname_lst = [] 
    alignments = pysam.AlignmentFile(bam_file, "rb")
    for line in alignments.fetch(chrom):
        ccs = BAM(line)
        if not ccs.is_primary:
            continue
        qname_lst.append("/".join(ccs.qname.split("/")[0:2]))
    chrom2qname_lst[chrom] = qname_lst


def dump_zmw(bam_file, region, region_file, threads):

    p = mp.Pool(threads)
    manager = mp.Manager()
    chrom2qname_lst = manager.dict()
    chrom_lst = load_chrom_lst(bam_file, region, region_file)
    get_qname_arg_lst = [
        (
            chrom,
            bam_file,
            chrom2qname_lst
        )
        for chrom in chrom_lst
    ]
    p.starmap(
        get_qname, get_qname_arg_lst,
    )
    p.close()
    p.join()
 
    movie2zmw_lst = defaultdict(list)
    for chrom in chrom_lst:
        for qname in chrom2qname_lst[chrom]:
            movie, zmw = qname.split("/")
            movie2zmw_lst[movie].append(zmw) 

    for movie, zmw_lst in movie2zmw_lst.items():
        o = open("{}.zmw".format(movie), "w")
        for zmw in natsort.natsorted(zmw_lst):
            o.write("{}\n".format(zmw))
        o.close()
  

def main():
    options = parse_args(sys.argv)
    dump_zmw(
        options.input, 
        options.region, 
        options.region_list, 
        options.threads, 
    )
    sys.exit(0)


if __name__ == "__main__":
    main()
