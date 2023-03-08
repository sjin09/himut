
#!/usr/bin/env python3

import sys
import math
import pysam
import pyfastx
import argparse
import numpy as np
import himut.util
import himut.gtlib
import himut.bamlib
import himut.caller
import himut.normcounts
import multiprocessing as mp
from collections import defaultdict
from typing import Dict, List, Tuple, Set


def init_allelecounts():
    rpos2allelecounts = defaultdict(lambda: np.zeros(6))
    rpos2allele2bq_lst = defaultdict(lambda: {0: [], 1: [], 2: [], 3: []})
    return rpos2allelecounts, rpos2allele2bq_lst


def parse_args(args):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-i",
        "--bam",
        type=str,
        required=True,
        help="BAM file to read",
    )
    parser.add_argument(
        "--ref",
        type=str,
        required=True,
        help="ref FASTA file to read",
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
        help="list of target chromosomes separated by new line",
    )
    parser.add_argument(
        "--min_gq",
        type=int,
        default=20,
        required=False,
        help="minimum germline genotype quality (GQ) score ",
    )
    parser.add_argument(
        "--germline_snv_prior",
        type=float,
        default=1/(10**3),
        required=False,
        help="germline snv prior",
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=1,
        required=False,
        help="number of threads to use",
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


def get_match_mismatch_count_per_bq(
    chrom: str,
    chunkloci_list: List[Tuple[str, int, int]],
    bam_file: str,
    chrom_seq: str,
    min_gq: int,
    germline_snv_prior: float,
    chrom2bq2match_count: Dict[str, Dict[int, int]],
    chrom2bq2mismatch_count: Dict[str, Dict[int, int]], 
):

    himut.gtlib.init(germline_snv_prior)
    bq2match_count = {i: 0 for i in range(1,94)}
    bq2mismatch_count = {i: 0 for i in range(1,94)}
    alignments = pysam.AlignmentFile(bam_file, "rb")
    for (chrom, chunk_start, chunk_end) in chunkloci_list: 
        rpos2allelecounts, rpos2allele2bq_lst= init_allelecounts() 
        for i in alignments.fetch(chrom, chunk_start, chunk_end): # iterate through reads
            ccs = himut.bamlib.BAM(i)
            if not ccs.is_primary:
                continue
            himut.normcounts.update_allelecounts(ccs, rpos2allelecounts, rpos2allele2bq_lst)

        for rpos in range(chunk_start, chunk_end): # iterate through reference positions
            ref = chrom_seq[rpos]
            if not ref in himut.util.base_set:
                continue

            allelecounts = rpos2allelecounts[rpos]
            allele2bq_lst = rpos2allele2bq_lst[rpos]
            del_count, ins_count, read_depth = himut.bamlib.get_read_depth(allelecounts)
            if (del_count != 0 or ins_count != 0):
                continue
            
            germ_gt, germ_gq, germ_gt_state, _ = himut.gtlib.get_germ_gt(ref, allele2bq_lst)
            if himut.caller.is_low_gq(germ_gq, min_gq):
                continue
            if germ_gt_state == "het" or germ_gt_state == "hetalt":
                base_sum = sum([allelecounts[himut.util.base2idx[base]] for base in list(germ_gt)])
                if base_sum == read_depth:
                    for base in list(germ_gt):
                        for bq in allele2bq_lst[himut.util.base2idx[base]]:
                            bq2match_count[bq] += 1
                    continue 
            elif germ_gt_state == "homalt" or germ_gt_state == "homref":
                base_sum = allelecounts[himut.util.base2idx[germ_gt[0]]]
                if base_sum == read_depth:
                    for bq in allele2bq_lst[himut.util.base2idx[germ_gt[0]]]:
                        bq2match_count[bq] += 1
                    continue
            germ_gt_base_set = set(list(germ_gt))
            for base in list(himut.util.base_set.difference(germ_gt_base_set)):
                for bq in allele2bq_lst[himut.util.base2idx[base]]:
                    bq2mismatch_count[bq] += 1
    chrom2bq2match_count[chrom] = bq2match_count
    chrom2bq2mismatch_count[chrom] = bq2mismatch_count


def dump_match_mismatch_count_per_bq(
    bam_file: str,
    ref_file: str,
    region: str,
    region_list: str,
    min_gq: int,
    germline_snv_prior: float,
    threads: int,
    out_file: str
): 

    p = mp.Pool(threads)
    manager = mp.Manager()
    ref_seq = pyfastx.Fasta(ref_file)
    chrom2bq2match_count = manager.dict()
    chrom2bq2mismatch_count = manager.dict()
    _, tname2tsize = himut.bamlib.get_tname2tsize(bam_file)
    chrom_lst, chrom2chunkloci_lst = himut.util.load_loci(region, region_list, tname2tsize)
    match_mismatch_count_arg_lst = [
        (
            chrom,
            chrom2chunkloci_lst[chrom],
            bam_file,
            str(ref_seq[chrom]),
            min_gq,
            germline_snv_prior,
            chrom2bq2match_count,
            chrom2bq2mismatch_count,
        )
        for chrom in chrom_lst
    ]
    p.starmap(
        get_match_mismatch_count_per_bq, match_mismatch_count_arg_lst,
    )
    p.close()
    p.join()

    o = open(out_file, "w")
    bq2match_count = defaultdict(lambda: 0)
    bq2mismatch_count = defaultdict(lambda: 0)
    o.write("{}\t{}\t{}\t{}\n".format("bq", "mismatch", "match", "pq"))
    for chrom in chrom_lst:
        for bq in range(1, 94): 
            bq2match_count[bq] += chrom2bq2match_count[chrom][bq]
            bq2mismatch_count[bq] += chrom2bq2mismatch_count[chrom][bq]
            
    for bq in range(1, 94): 
        match_count = bq2match_count[bq]
        mismatch_count = bq2mismatch_count[bq]
        if match_count != 0 and mismatch_count != 0:
            pq = -10*math.log10(mismatch_count/float(match_count))
            o.write("{}\t{}\t{}\t{}\n".format(bq, mismatch_count, match_count, pq))
        else: 
            o.write("{}\t{}\t{}\t{}\n".format(bq, mismatch_count, match_count, "NA"))
    o.close()


def main():
    options = parse_args(sys.argv)
    dump_match_mismatch_count_per_bq(
        options.bam,
        options.ref,
        options.region,
        options.region_list, 
        options.min_gq,
        options.germline_snv_prior,
        options.threads, 
        options.out
    )
    sys.exit(0)


if __name__ == "__main__":
    main()
