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
    chrom2ccs_len_lst: Dict[str, List[int]],
    chrom2q93_base_count: Dict[str, int],
):

    # counter = 0
    ccs_len_lst = []
    q93_base_count = 0
    alignments = pysam.AlignmentFile(bam_file, "rb")
    for line in alignments.fetch(chrom):
        ccs = himut.bamlib.BAM(line)
        if not ccs.is_primary:
            continue
        # counter += 1
        # if counter > 100:
        #     break
        ccs_len_lst.append(ccs.qlen)
        q93_base_count += ccs.bq_int_lst.count(93)
    chrom2ccs_len_lst[chrom] = ccs_len_lst
    chrom2q93_base_count[chrom] = q93_base_count


def dump_seq_cov(
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
    chrom2ccs_len_lst = manager.dict()
    chrom2q93_base_count = manager.dict()
    get_ccs_statistics_arg_lst = [
        (
            tgt,
            bam_file,
            chrom2ccs_len_lst,
            chrom2q93_base_count
        )
        for tgt in tgt_lst
    ]
    p.starmap(
        get_ccs_statistics, get_ccs_statistics_arg_lst,
    )
    p.close()
    p.join()
 
    ccs_sum = 0
    ccs_count = 0
    q93_count = 0
    genome_ccs_len_lst = []
    o = open(out_file, "w")
    tgt_sum = sum([tname2tsize[tgt] for tgt in tgt_lst])
    for tgt in tgt_lst:
        ccs_len_lst = chrom2ccs_len_lst[tgt]
        genome_ccs_len_lst.extend(ccs_len_lst)
        q93_count += chrom2q93_base_count[tgt]
        ccs_count += len(ccs_len_lst)
        ccs_sum += sum(ccs_len_lst)
    tgt_cov = ccs_sum/tgt_sum
    q93_proportion = (q93_count/ccs_sum)
    ccs_len_std = np.std(genome_ccs_len_lst)
    ccs_mean_len = np.mean(genome_ccs_len_lst)
    o.write("{:,}\t{:,}\t{:.2f} ± {:.2f}\t{:.2%}\t{:,}\t{:.2f}\n".format(ccs_count, ccs_sum, ccs_mean_len, ccs_len_std, q93_proportion, tgt_sum, tgt_cov))
    o.close() 


def main():
    options = parse_args(sys.argv)
    dump_seq_cov(
        options.bam, 
        options.region, 
        options.region_list,
        options.threads, 
        options.out
    )
    sys.exit(0)


if __name__ == "__main__":
    main()
