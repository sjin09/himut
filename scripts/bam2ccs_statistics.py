#!/usr/bin/env python3

import sys
import pysam
import natsort
import argparse
import numpy as np
import himut.bamlib
import multiprocessing as mp
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
        print("--region or --region_list parameter has not been provided")
        print("defaulting to all the contigs and chromosomes")
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


def get_ccs_statistics(
    chrom: str,
    bam_file: str,
    chrom2qlen_lst: Dict[str, List[int]] 
):

    qlen_lst = []
    alignments = pysam.AlignmentFile(bam_file, "rb")
    for line in alignments.fetch(chrom):
        ccs = himut.bamlib.BAM(line)
        if not ccs.is_primary:
            continue
        qlen_lst.append(ccs.qlen)
    chrom2qlen_lst[chrom] = qlen_lst


def dump_ccs_statistics(
    bam_file: str,
    region: str,
    region_file: str,
    threads: int,
    out_file: str
): 
    
    tname2tsize = get_tname2tsize(bam_file)
    chrom_lst = load_chrom_lst(region, region_file)
    if len(chrom_lst) == 0:
        tgt_lst = natsort.natsorted(tname2tsize.keys())
    else:
        tgt_lst = chrom_lst

    p = mp.Pool(threads)
    manager = mp.Manager()
    chrom2qlen_lst = manager.dict()
    get_ccs_statistics_arg_lst = [
        (
            tgt,
            bam_file,
            chrom2qlen_lst
        )
        for tgt in tgt_lst
    ]
    p.starmap(
        get_ccs_statistics, get_ccs_statistics_arg_lst,
    )
    p.close()
    p.join()
 
    qsum = 0
    genome_qlen_lst = []
    o = open(out_file, "w")
    tsum = sum([tname2tsize[tgt] for tgt in tgt_lst])
    o.write("{}\n".format("\t".join(["ccs_base", "assembly_base", "coverage", "qlen_mean", "qlen_std"])))
    for tgt in tgt_lst:
        for p in chrom2qlen_lst[tgt]: 
            chrom_qlen_lst = chrom2qlen_lst[tgt]
            genome_qlen_lst.extend(chrom_qlen_lst)
            qsum += sum(chrom_qlen_lst)
    coverage = qsum/tsum
    qlen_mean = np.mean(genome_qlen_lst)
    qlen_std = np.std(genome_qlen_lst)
    o.write("{}\t{}\t{}\t{}\t{}\n".format(qsum, tsum, coverage, qlen_mean, qlen_std))
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
