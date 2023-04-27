
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


def get_count_per_qname(
    chrom: str,
    chunkloci_list: List[Tuple[str, int, int]],
    bam_file: str,
    min_gq: int,
    min_mapq: int,
    germline_snv_prior: float,
    chrom2qname2qv: Dict[str, Dict[str, float]],
    chrom2qname2count: Dict[str, Dict[str, int]],
):

    qname2qv = defaultdict()
    qname2count = defaultdict(lambda: 0)
    himut.gtlib.init(germline_snv_prior)
    alignments = pysam.AlignmentFile(bam_file, "rb")
    for (chrom, chunk_start, chunk_end) in chunkloci_list: 
        chunk_tsbs_lst = []
        qname2tsbs_lst = {}
        rpos2allelecounts, rpos2allele2bq_lst= init_allelecounts() 
        for i in alignments.fetch(chrom, chunk_start, chunk_end): # iterate through reads
            ccs = himut.bamlib.BAM(i)
            if not ccs.is_primary:
                continue
            update_allelecounts(ccs, rpos2allelecounts, rpos2allele2bq_lst)
            if himut.caller.is_low_mapq(ccs.mapq, min_mapq):
                continue
            ccs.cs2subindel()
            qname2count[ccs.qname] = 0
            qname2qv[ccs.qname] = ccs.get_qv()
            chunk_tsbs_lst.extend(ccs.tsbs_lst)
            qname2tsbs_lst[ccs.qname] = ccs.tsbs_lst
            
        germ_tsbs_set = set() 
        for (tpos, ref, alt) in natsort.natsorted(list(set(chunk_tsbs_lst))): # iterate through substitutions
            rpos = tpos - 1
            tsbs = tpos, ref, alt
            som_gt = "{}{}".format(ref, alt)
            allele2bq_lst = rpos2allele2bq_lst[rpos]
            _, _, germ_gt_state, gt2gt_state = himut.gtlib.get_germ_gt(ref, allele2bq_lst)
            if germ_gt_state == "homref":
                continue
            germ_gq = himut.gtlib.get_germ_gq(som_gt, gt2gt_state, allele2bq_lst)
            if himut.caller.is_low_gq(germ_gq, min_gq): 
                continue
            germ_tsbs_set.add(tsbs)
        
        for qname, tsbs_lst in qname2tsbs_lst.items(): 
            for tsbs in tsbs_lst:
                if tsbs in germ_tsbs_set:
                    continue
                qname2count[qname] += 1

    chrom2qname2qv[chrom] = dict(qname2qv)
    chrom2qname2count[chrom] = dict(qname2count)


def dump_count_per_qname(
    bam_file: str,
    region: str,
    region_list: str,
    min_gq: int,
    min_mapq: int,
    germline_snv_prior: float,
    threads: int,
    out_file: str
): 

    p = mp.Pool(threads)
    manager = mp.Manager()
    chrom2qname2qv = manager.dict()
    chrom2qname2count = manager.dict()
    _, tname2tsize = himut.bamlib.get_tname2tsize(bam_file)
    chrom_lst, chrom2chunkloci_lst = himut.util.load_loci(region, region_list, tname2tsize)
    get_count_per_qname_arg_lst = [
        (
            chrom,
            chrom2chunkloci_lst[chrom],
            bam_file,
            min_gq,
            min_mapq,
            germline_snv_prior,
            chrom2qname2qv,
            chrom2qname2count,
        )
        for chrom in chrom_lst
    ]
    p.starmap(
        get_count_per_qname, get_count_per_qname_arg_lst,
    )
    p.close()
    p.join()

    o = open(out_file, "w") # return
    o.write("{}\t{}\t{}\n".format("qname", "qv", "count"))
    for chrom in chrom_lst:
        for qname, qv in chrom2qname2qv[chrom].items():
            count = chrom2qname2count[chrom][qname] 
            o.write("{}\t{}\t{}\n".format(qname, qv, count))    
    o.close()


def main():
    options = parse_args(sys.argv)
    dump_count_per_qname(
        options.bam,
        options.region,
        options.region_list, 
        options.min_gq,
        options.min_mapq,
        options.germline_snv_prior,
        options.threads, 
        options.out
    )
    sys.exit(0)


if __name__ == "__main__":
    main()
