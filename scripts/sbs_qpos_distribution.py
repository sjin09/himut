
#!/usr/bin/env python3

import sys
import pysam
import natsort
import argparse
import numpy as np
import himut.util
import himut.gtlib
import himut.bamlib
import himut.caller
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
        help="BAM file to read",
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
        "--min_bq",
        type=int,
        default=93,
        required=False,
        help="minimum base quality (BQ) score ",
    )
    parser.add_argument(
        "--min_gq",
        type=int,
        default=20,
        required=False,
        help="minimum germline genotype quality (GQ) score ",
    )
    parser.add_argument(
        "--min_mapq",
        type=int,
        default=60,
        required=False,
        help="minimum mapping quality (MAPQ) score",
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


def init_allelecounts():
    rpos2allelecounts = defaultdict(lambda: np.zeros(6))
    rpos2allele2bq_lst = defaultdict(lambda: {0: [], 1: [], 2: [], 3: []})
    return rpos2allelecounts, rpos2allele2bq_lst


def update_allelecounts(
    ccs,
    rpos2allelecounts: Dict[int, np.ndarray],
    rpos2allele2bq_lst: Dict[int, Dict[int, List[int]]],
):

    tpos = ccs.tstart
    qpos = ccs.qstart
    for (state, ref, alt, ref_len, alt_len) in ccs.cstuple_lst:
        if state == 1:  # match
            for i, alt_base in enumerate(alt):
                epos = tpos + i
                bidx = himut.util.base2idx[alt_base]
                rpos2allelecounts[epos][bidx] += 1
                rpos2allele2bq_lst[epos][bidx].append(ccs.bq_int_lst[qpos + i])
        elif state == 2:  # sub
            bidx = himut.util.base2idx[alt]
            rpos2allelecounts[tpos][bidx] += 1
            rpos2allele2bq_lst[tpos][bidx].append(ccs.bq_int_lst[qpos])
        elif state == 3:  # insertion
            rpos2allelecounts[tpos][4] += 1
        elif state == 4:  # deletion
            for j in range(len(ref[1:])):
                rpos2allelecounts[tpos + j][5] += 1
        tpos += ref_len
        qpos += alt_len


def get_qpos_distribution(
    chrom: str,
    chunkloci_list: List[Tuple[str, int, int]],
    bam_file: str,
    min_bq: int,
    min_gq: int,
    min_mapq: int,
    germline_snv_prior: float,
    chrom2qpos2count: Dict[str, Dict[float, int]],
):

    som_seen = set() 
    himut.gtlib.init(germline_snv_prior)
    qpos2count = {(i/100): 0 for i in range(101)}
    alignments = pysam.AlignmentFile(bam_file, "rb")
    for (chrom, chunk_start, chunk_end) in chunkloci_list: 
        chunk_tsbs_lst = []
        tsbs2qpos = defaultdict(list)
        rpos2allelecounts, rpos2allele2bq_lst= init_allelecounts() 
        for i in alignments.fetch(chrom, chunk_start, chunk_end): # iterate through reads
            ccs = himut.bamlib.BAM(i)
            if not ccs.is_primary:
                continue
            update_allelecounts(ccs, rpos2allelecounts, rpos2allele2bq_lst)
            if himut.caller.is_low_mapq(ccs.mapq, min_mapq):
                continue
            ccs.cs2subindel()
            chunk_tsbs_lst.extend(ccs.tsbs_lst)
            for (tsbs), (qpos, _, _ ) in zip(ccs.tsbs_lst, ccs.qsbs_lst):
                if tsbs in som_seen:
                    continue
                normalised_qpos = round(qpos/ccs.qlen, 2)
                tsbs2qpos[tsbs].append(normalised_qpos)
                
        for (tpos, ref, alt) in natsort.natsorted(list(set(chunk_tsbs_lst))): # iterate through substitutions
            if (tpos, ref, alt) in som_seen: 
                continue
            
            if not himut.caller.is_chunk(tpos, chunk_start, chunk_end):
                continue

            rpos = tpos - 1
            tsbs = (tpos, ref, alt)
            som_gt = "{}{}".format(ref, alt)
            allelecounts = rpos2allelecounts[rpos]
            allele2bq_lst = rpos2allele2bq_lst[rpos]
            _, _, germ_gt_state, gt2gt_state = himut.gtlib.get_germ_gt(ref, allele2bq_lst)
            if germ_gt_state != "homref":
                continue
           
            germ_gq = himut.gtlib.get_germ_gq(som_gt, gt2gt_state, allele2bq_lst)
            if himut.caller.is_low_gq(germ_gq, min_gq): 
                continue
            
            if himut.caller.is_low_bq(alt, min_bq, allele2bq_lst):
                continue
          
            alt_count = allelecounts[himut.util.base2idx[alt]] 
            if alt_count != 1:
                continue
            som_seen.add((tpos, ref, alt))
            for qpos in tsbs2qpos[tsbs]:
                qpos2count[qpos] += 1

    chrom2qpos2count[chrom] = qpos2count


def dump_qpos_distribution(
    bam_file: str,
    region: str,
    region_list: str,
    min_bq: int,
    min_gq: int,
    min_mapq: int,
    germline_snv_prior: float,
    threads: int,
    out_file: str
): 

    p = mp.Pool(threads)
    manager = mp.Manager()
    chrom2qpos2count = manager.dict()
    _, tname2tsize = himut.bamlib.get_tname2tsize(bam_file)
    chrom_lst, chrom2chunkloci_lst = himut.util.load_loci(region, region_list, tname2tsize)
    get_qpos_distribution_arg_lst = [
        (
            chrom,
            chrom2chunkloci_lst[chrom],
            bam_file,
            min_bq,
            min_gq,
            min_mapq,
            germline_snv_prior,
            chrom2qpos2count,
        )
        for chrom in chrom_lst
    ]
    p.starmap(
        get_qpos_distribution, get_qpos_distribution_arg_lst,
    )
    p.close()
    p.join()

    qpos2count = defaultdict(lambda: 0)
    for chrom in chrom_lst:
        for qpos, count in chrom2qpos2count[chrom].items():
            qpos2count[qpos] += count  

    o = open(out_file, "w")
    o.write("{}\t{}\n".format("qpos", "count"))
    for qpos, count in qpos2count.items():
        o.write("{}\t{}\n".format(qpos, count))    
    o.close()


def main():
    options = parse_args(sys.argv)
    dump_qpos_distribution(
        options.bam,
        options.region,
        options.region_list, 
        options.min_bq,
        options.min_gq,
        options.min_mapq,
        options.germline_snv_prior,
        options.threads, 
        options.out
    )
    sys.exit(0)


if __name__ == "__main__":
    main()
