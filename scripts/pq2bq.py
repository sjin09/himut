#!/usr/bin/env python3

import re  # packages
import sys
import math
import time
import pysam
import bisect
import random
import natsort
import pyfastx
import argparse
import numpy as np
import multiprocessing as mp
from collections import defaultdict
from typing import Dict, List, Tuple, Set


purine = set(["A", "G"])
base_lst = list("ATGC")  # global
base_set = set(base_lst)
allele_lst = list("ATGC+-")
sub_lst = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
purine2pyrimidine = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
base2idx = {base: idx for idx, base in enumerate(base_lst)}
idx2base = {idx: base for idx, base in enumerate(base_lst)}
allele2idx = {allele: idx for idx, allele in enumerate(allele_lst)}
idx2allele = {idx: allele for idx, allele in enumerate(allele_lst)}
gt_lst = ["AA", "TA", "CA", "GA", "TT", "CT", "GT", "CC", "GC", "GG"]
sub_lst = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
sub2idx = {sub: idx for idx, sub in enumerate(sub_lst)}
tri_lst = [
    "ACA",
    "ACC",
    "ACG",
    "ACT",
    "ATA",
    "ATC",
    "ATG",
    "ATT",
    "CCA",
    "CCC",
    "CCG",
    "CCT",
    "CTA",
    "CTC",
    "CTG",
    "CTT",
    "GCA",
    "GCC",
    "GCG",
    "GCT",
    "GTA",
    "GTC",
    "GTG",
    "GTT",
    "TCA",
    "TCC",
    "TCG",
    "TCT",
    "TTA",
    "TTC",
    "TTG",
    "TTT",
    "NNN",
]
tri_count = len(tri_lst)
tri2idx = {tri: idx for idx, tri in enumerate(tri_lst)}
sbs96_lst = [
    "A[C>A]A",
    "A[C>A]C",
    "A[C>A]G",
    "A[C>A]T",
    "A[C>G]A",
    "A[C>G]C",
    "A[C>G]G",
    "A[C>G]T",
    "A[C>T]A",
    "A[C>T]C",
    "A[C>T]G",
    "A[C>T]T",
    "A[T>A]A",
    "A[T>A]C",
    "A[T>A]G",
    "A[T>A]T",
    "A[T>C]A",
    "A[T>C]C",
    "A[T>C]G",
    "A[T>C]T",
    "A[T>G]A",
    "A[T>G]C",
    "A[T>G]G",
    "A[T>G]T",
    "C[C>A]A",
    "C[C>A]C",
    "C[C>A]G",
    "C[C>A]T",
    "C[C>G]A",
    "C[C>G]C",
    "C[C>G]G",
    "C[C>G]T",
    "C[C>T]A",
    "C[C>T]C",
    "C[C>T]G",
    "C[C>T]T",
    "C[T>A]A",
    "C[T>A]C",
    "C[T>A]G",
    "C[T>A]T",
    "C[T>C]A",
    "C[T>C]C",
    "C[T>C]G",
    "C[T>C]T",
    "C[T>G]A",
    "C[T>G]C",
    "C[T>G]G",
    "C[T>G]T",
    "G[C>A]A",
    "G[C>A]C",
    "G[C>A]G",
    "G[C>A]T",
    "G[C>G]A",
    "G[C>G]C",
    "G[C>G]G",
    "G[C>G]T",
    "G[C>T]A",
    "G[C>T]C",
    "G[C>T]G",
    "G[C>T]T",
    "G[T>A]A",
    "G[T>A]C",
    "G[T>A]G",
    "G[T>A]T",
    "G[T>C]A",
    "G[T>C]C",
    "G[T>C]G",
    "G[T>C]T",
    "G[T>G]A",
    "G[T>G]C",
    "G[T>G]G",
    "G[T>G]T",
    "T[C>A]A",
    "T[C>A]C",
    "T[C>A]G",
    "T[C>A]T",
    "T[C>G]A",
    "T[C>G]C",
    "T[C>G]G",
    "T[C>G]T",
    "T[C>T]A",
    "T[C>T]C",
    "T[C>T]G",
    "T[C>T]T",
    "T[T>A]A",
    "T[T>A]C",
    "T[T>A]G",
    "T[T>A]T",
    "T[T>C]A",
    "T[T>C]C",
    "T[T>C]G",
    "T[T>C]T",
    "T[T>G]A",
    "T[T>G]C",
    "T[T>G]G",
    "T[T>G]T",
    "N[N>N]N",
]
sbs96_to_idx = {sbs96: idx for idx, sbs96 in enumerate(sbs96_lst)}
sbs96_to_sub = {
    "A[C>A]A": "C>A",
    "A[C>A]T": "C>A",
    "A[C>A]G": "C>A",
    "A[C>A]C": "C>A",
    "T[C>A]A": "C>A",
    "T[C>A]T": "C>A",
    "T[C>A]G": "C>A",
    "T[C>A]C": "C>A",
    "G[C>A]A": "C>A",
    "G[C>A]T": "C>A",
    "G[C>A]G": "C>A",
    "G[C>A]C": "C>A",
    "C[C>A]A": "C>A",
    "C[C>A]T": "C>A",
    "C[C>A]G": "C>A",
    "C[C>A]C": "C>A",
    "A[C>G]A": "C>G",
    "A[C>G]T": "C>G",
    "A[C>G]G": "C>G",
    "A[C>G]C": "C>G",
    "T[C>G]A": "C>G",
    "T[C>G]T": "C>G",
    "T[C>G]G": "C>G",
    "T[C>G]C": "C>G",
    "G[C>G]A": "C>G",
    "G[C>G]T": "C>G",
    "G[C>G]G": "C>G",
    "G[C>G]C": "C>G",
    "C[C>G]A": "C>G",
    "C[C>G]T": "C>G",
    "C[C>G]G": "C>G",
    "C[C>G]C": "C>G",
    "A[C>T]A": "C>T",
    "A[C>T]T": "C>T",
    "A[C>T]G": "C>T",
    "A[C>T]C": "C>T",
    "T[C>T]A": "C>T",
    "T[C>T]T": "C>T",
    "T[C>T]G": "C>T",
    "T[C>T]C": "C>T",
    "G[C>T]A": "C>T",
    "G[C>T]T": "C>T",
    "G[C>T]G": "C>T",
    "G[C>T]C": "C>T",
    "C[C>T]A": "C>T",
    "C[C>T]T": "C>T",
    "C[C>T]G": "C>T",
    "C[C>T]C": "C>T",
    "A[T>A]A": "T>A",
    "A[T>A]T": "T>A",
    "A[T>A]G": "T>A",
    "A[T>A]C": "T>A",
    "T[T>A]A": "T>A",
    "T[T>A]T": "T>A",
    "T[T>A]G": "T>A",
    "T[T>A]C": "T>A",
    "G[T>A]A": "T>A",
    "G[T>A]T": "T>A",
    "G[T>A]G": "T>A",
    "G[T>A]C": "T>A",
    "C[T>A]A": "T>A",
    "C[T>A]T": "T>A",
    "C[T>A]G": "T>A",
    "C[T>A]C": "T>A",
    "A[T>C]A": "T>C",
    "A[T>C]T": "T>C",
    "A[T>C]G": "T>C",
    "A[T>C]C": "T>C",
    "T[T>C]A": "T>C",
    "T[T>C]T": "T>C",
    "T[T>C]G": "T>C",
    "T[T>C]C": "T>C",
    "G[T>C]A": "T>C",
    "G[T>C]T": "T>C",
    "G[T>C]G": "T>C",
    "G[T>C]C": "T>C",
    "C[T>C]A": "T>C",
    "C[T>C]T": "T>C",
    "C[T>C]G": "T>C",
    "C[T>C]C": "T>C",
    "A[T>G]A": "T>G",
    "A[T>G]T": "T>G",
    "A[T>G]G": "T>G",
    "A[T>G]C": "T>G",
    "T[T>G]A": "T>G",
    "T[T>G]T": "T>G",
    "T[T>G]G": "T>G",
    "T[T>G]C": "T>G",
    "G[T>G]A": "T>G",
    "G[T>G]T": "T>G",
    "G[T>G]G": "T>G",
    "G[T>G]C": "T>G",
    "C[T>G]A": "T>G",
    "C[T>G]T": "T>G",
    "C[T>G]G": "T>G",
    "C[T>G]C": "T>G",
}
sbs96_to_tri = {
    "A[C>A]A": "ACA",
    "A[C>A]T": "ACT",
    "A[C>A]G": "ACG",
    "A[C>A]C": "ACC",
    "T[C>A]A": "TCA",
    "T[C>A]T": "TCT",
    "T[C>A]G": "TCG",
    "T[C>A]C": "TCC",
    "G[C>A]A": "GCA",
    "G[C>A]T": "GCT",
    "G[C>A]G": "GCG",
    "G[C>A]C": "GCC",
    "C[C>A]A": "CCA",
    "C[C>A]T": "CCT",
    "C[C>A]G": "CCG",
    "C[C>A]C": "CCC",
    "A[C>G]A": "ACA",
    "A[C>G]T": "ACT",
    "A[C>G]G": "ACG",
    "A[C>G]C": "ACC",
    "T[C>G]A": "TCA",
    "T[C>G]T": "TCT",
    "T[C>G]G": "TCG",
    "T[C>G]C": "TCC",
    "G[C>G]A": "GCA",
    "G[C>G]T": "GCT",
    "G[C>G]G": "GCG",
    "G[C>G]C": "GCC",
    "C[C>G]A": "CCA",
    "C[C>G]T": "CCT",
    "C[C>G]G": "CCG",
    "C[C>G]C": "CCC",
    "A[C>T]A": "ACA",
    "A[C>T]T": "ACT",
    "A[C>T]G": "ACG",
    "A[C>T]C": "ACC",
    "T[C>T]A": "TCA",
    "T[C>T]T": "TCT",
    "T[C>T]G": "TCG",
    "T[C>T]C": "TCC",
    "G[C>T]A": "GCA",
    "G[C>T]T": "GCT",
    "G[C>T]G": "GCG",
    "G[C>T]C": "GCC",
    "C[C>T]A": "CCA",
    "C[C>T]T": "CCT",
    "C[C>T]G": "CCG",
    "C[C>T]C": "CCC",
    "A[T>A]A": "ATA",
    "A[T>A]T": "ATT",
    "A[T>A]G": "ATG",
    "A[T>A]C": "ATC",
    "T[T>A]A": "TTA",
    "T[T>A]T": "TTT",
    "T[T>A]G": "TTG",
    "T[T>A]C": "TTC",
    "G[T>A]A": "GTA",
    "G[T>A]T": "GTT",
    "G[T>A]G": "GTG",
    "G[T>A]C": "GTC",
    "C[T>A]A": "CTA",
    "C[T>A]T": "CTT",
    "C[T>A]G": "CTG",
    "C[T>A]C": "CTC",
    "A[T>C]A": "ATA",
    "A[T>C]T": "ATT",
    "A[T>C]G": "ATG",
    "A[T>C]C": "ATC",
    "T[T>C]A": "TTA",
    "T[T>C]T": "TTT",
    "T[T>C]G": "TTG",
    "T[T>C]C": "TTC",
    "G[T>C]A": "GTA",
    "G[T>C]T": "GTT",
    "G[T>C]G": "GTG",
    "G[T>C]C": "GTC",
    "C[T>C]A": "CTA",
    "C[T>C]T": "CTT",
    "C[T>C]G": "CTG",
    "C[T>C]C": "CTC",
    "A[T>G]A": "ATA",
    "A[T>G]T": "ATT",
    "A[T>G]G": "ATG",
    "A[T>G]C": "ATC",
    "T[T>G]A": "TTA",
    "T[T>G]T": "TTT",
    "T[T>G]G": "TTG",
    "T[T>G]C": "TTC",
    "G[T>G]A": "GTA",
    "G[T>G]T": "GTT",
    "G[T>G]G": "GTG",
    "G[T>G]C": "GTC",
    "C[T>G]A": "CTA",
    "C[T>G]T": "CTT",
    "C[T>G]G": "CTG",
    "C[T>G]C": "CTC",
}

class BAM:
    def __init__(self, line):
        if line.is_secondary:
            self.is_primary = False
        else:
            self.is_primary = True
            self.tname = line.reference_name
            self.tstart = line.reference_start
            self.tend = line.reference_end
            # query
            self.qname = line.query_name
            self.qstart = line.query_alignment_start
            self.qend = line.query_alignment_end
            self.qseq = line.query_sequence
            self.qlen = len(self.qseq)
            self.mapq = line.mapping_quality
            self.pq_int_lst = line.query_qualities
            self.cs2tuple(line.get_tag("cs"))

    def cs2lst(self, cs_tag):
        cs_lst = [
            cs for cs in re.split("(:[0-9]+|\*[a-z][a-z]|[=\+\-][A-Za-z]+)", cs_tag)
        ]
        cs_lst = [cs.upper() for cs in cs_lst if cs != ""]
        return cs_lst

    def cs2tuple(self, cs_tag) -> List[Tuple[int, str, str, int, int]]:
        qpos = self.qstart
        self.cstuple_lst = []
        cs_lst = self.cs2lst(cs_tag)
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
            qpos += mlen
            self.cstuple_lst.append(t)

    def cs2subindel(self):

        tpos = self.tstart  # 0-coordinate
        qpos = self.qstart
        self.mismatch_tpos_lst = []
        for (state, ref, alt, ref_len, alt_len) in self.cstuple_lst:
            if ref != "N" and state == 2:  # snp
                self.mismatch_tpos_lst.append(tpos)
            elif state == 3 or state == 4:  # insertion or deletion
                self.mismatch_tpos_lst.append(tpos)
            tpos += ref_len
            qpos += alt_len

    def get_trimmed_range(self, min_trim):
        trimmed_qstart = int(self.qlen * min_trim)
        trimmed_qend = int(self.qlen * (1 - min_trim))
        return trimmed_qstart, trimmed_qend

    def is_trimmed_base(self, qpos: int, trimmed_qstart: int, trimmed_qend: int):
        if qpos < trimmed_qstart:
            return True
        elif qpos > trimmed_qend:
            return True
        else:
            return False

    def get_mismatch_range(self, tpos: int, qpos: int, qlen: int, window: int):
        qstart, qend = [qpos - window, qpos + window]
        if qstart < 0:
            urange = window + qstart
            drange = window + abs(qstart)
        elif qend > qlen:
            urange = window + abs(qend - qlen)
            drange = qlen - qpos
        else:
            urange = window
            drange = window
        tstart = tpos - urange
        tend = tpos + drange
        return tstart, tend

    def is_mismatch_base(self, mismatch_count, max_mismatch_count):
        if mismatch_count > max_mismatch_count:
            return True
        else:
            return False


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
        "--min_mapq",
        type=int,
        default=60,
        required=False,
        help="minimum mapping quality (MAPQ) score ",
    )
    parser.add_argument(
        "--min_trim",
        type=float,
        default=0.01,
        required=False,
        help="minimum proportion of bases to be trimmed from the start and end of the read",
    )
    parser.add_argument(
        "--mismatch_window",
        type=int,
        default=20,
        required=False,
        help="mismatch window size",
    )
    parser.add_argument(
        "--max_mismatch_count",
        type=int,
        default=0,
        required=False,
        help="maximum number of mismatches within the mismatch window",
    )
    parser.add_argument(
        "--germline_snv_prior",
        type=float,
        default=1 / (10**3),
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
        "-p",
        "--prefix",
        type=str,
        required=True,
        help="prefix for output file to write",
    )
    args = args[1:]
    return parser.parse_args(args)


def get_tname2tsize(bam_file: str) -> Tuple[List[str], Dict[str, int]]:

    tname2tsize = {}
    alignments = pysam.AlignmentFile(bam_file, "rb")
    bam_header_lst = str(alignments.header).strip().split("\n")
    for h in bam_header_lst:
        if h.startswith("@SQ"):
            arr = h.split("\t")
            tname = arr[1].replace("SN:", "")
            tsize = int(arr[2].replace("LN:", ""))
            tname2tsize[tname] = tsize
    alignments.close()
    return tname2tsize


def get_md_threshold(read_depth: int) -> int:
    md_threshold = math.ceil(read_depth + (4 * math.sqrt(read_depth)))
    return md_threshold


def get_qlen_threshold(qlen_lst: List[int]):

    qlen_std = np.std(qlen_lst)
    qlen_mean = math.ceil(np.mean(qlen_lst))
    qlen_lower_limit = (
        0
        if math.ceil(qlen_mean - 2 * qlen_std) < 0
        else math.ceil(qlen_mean - 2 * qlen_std)
    )
    qlen_upper_limit = math.ceil(qlen_mean + 2 * qlen_std)
    return qlen_lower_limit, qlen_upper_limit


def get_bam_thresholds(
    bam_file: str,
    chrom_lst: List[str],
    chrom2len: Dict[str, int],
) -> Tuple[int, int, int, int]:

    n = 100
    ccs_sum = 0
    qlen_lst = []
    random.seed(10)
    chunk_len = 100000
    chunk_sum = n * chunk_len * len(chrom_lst)
    alignments = pysam.AlignmentFile(bam_file, "rb")
    for chrom in chrom_lst:
        chrom_len = chrom2len[chrom]
        random_start_lst = random.sample(range(chrom_len), n)
        for start in random_start_lst:
            end = start + chunk_len
            for line in alignments.fetch(chrom, start, end):
                ccs = BAM(line)
                if not ccs.is_primary:
                    continue
                if is_low_mapq(ccs.mapq, 0):
                    continue
                if start < ccs.tstart and ccs.tend < end:  # within
                    ccs_sum += ccs.qlen
                elif start > ccs.tstart and ccs.tend < end:  # upstream
                    ccs_sum += ccs.qlen - (start - ccs.tstart)
                elif start < ccs.tstart and ccs.tend > end:  # downstream
                    ccs_sum += ccs.qlen - (ccs.tend - end)
                else:
                    continue
                qlen_lst.append(ccs.qlen)
    alignments.close()

    read_depth = ccs_sum / float(chunk_sum)
    qlen_lower_limit, qlen_upper_limit = get_qlen_threshold(qlen_lst)
    return read_depth, qlen_lower_limit, qlen_upper_limit


def chunkloci(loci: Tuple[str, int, int]) -> List[Tuple[str, int, int]]:

    chrom, start, end = loci
    chrom_len = end - start
    if chrom_len > 200000:
        chunk_loci_lst = [(chrom, 1, 200000)]
        chunk_start_lst = list(range(200000, end, 200000))
        for i, chunk_start in enumerate(chunk_start_lst[:-1]):
            chunk_loci_lst.append((chrom, chunk_start, chunk_start_lst[i + 1]))
        if (chrom, chunk_start_lst[-1], end) not in chunk_loci_lst:
            chunk_loci_lst.append((chrom, chunk_start_lst[-1], end - 2))
    else:
        chunk_loci_lst = [(chrom, start, end)]
    return chunk_loci_lst


def load_chunkloci(loci: Tuple[str, int, int]) -> List[Tuple[str, int, int]]:

    chrom, start, end = loci
    chrom_len = end - start
    if chrom_len > 200000:
        chunk_loci_lst = [(chrom, 1, 200000)]
        chunk_start_lst = list(range(200000, end, 200000))
        for i, chunk_start in enumerate(chunk_start_lst[:-1]):
            chunk_loci_lst.append((chrom, chunk_start, chunk_start_lst[i + 1]))
        if (chrom, chunk_start_lst[-1], end) not in chunk_loci_lst:
            chunk_loci_lst.append((chrom, chunk_start_lst[-1], end - 2))
    else:
        chunk_loci_lst = [(chrom, start, end)]
    return chunk_loci_lst


def load_loci(
    region: str, region_list: str, tname2tsize: Dict[str, int]
) -> Tuple[List[str], List[Tuple[str, int, int]]]:

    chrom2loci_lst = defaultdict(list)
    if region is None and region_list is not None:
        for line in open(region_list).readlines():
            arr = line.strip().split()
            if len(arr) == 1:
                chrom2loci_lst[arr[0]].append((arr[0], 0, tname2tsize[arr[0]]))
            else:
                chrom2loci_lst[arr[0]].append((arr[0], int(arr[1]), int(arr[2])))
    elif region is not None and region_list is None:
        if region in tname2tsize:
            chrom2loci_lst[region].append((region, 0, tname2tsize[region]))
        else:
            print("{} does not exist in the BAM file".format(region))
            exit()
    elif region is not None and region_list is not None:
        for line in open(region_list).readlines():
            arr = line.strip().split()
            if len(arr) == 1:
                chrom2loci_lst[arr[0]].append((arr[0], 0, tname2tsize[arr[0]]))
            else:
                chrom2loci_lst[arr[0]].append((arr[0], int(arr[1]), int(arr[2])))
    else:
        for tname, tsize in tname2tsize.items():
            chrom2loci_lst[tname].append((tname, 0, tsize))

    chrom2chunkloci_lst = defaultdict(list)
    chrom_lst = natsort.natsorted(list(chrom2loci_lst.keys()))
    for chrom, loci_lst in chrom2loci_lst.items():
        chunkloci_lst = [
            chunkloci for loci in loci_lst for chunkloci in load_chunkloci(loci)
        ]
        chrom2chunkloci_lst[chrom] = chunkloci_lst
    return chrom_lst, chrom2chunkloci_lst


def gt_init(germline_snv_prior: float):
    global gt_state2gt_prior
    gt_state2gt_prior = {}
    gt_state2gt_prior["het"] = germline_snv_prior
    gt_state2gt_prior["hetalt"] = germline_snv_prior * germline_snv_prior * 2
    gt_state2gt_prior["homref"] = 1 - (
        (1.5 * germline_snv_prior) + (germline_snv_prior * germline_snv_prior)
    )
    gt_state2gt_prior["homalt"] = germline_snv_prior / 2


def get_epsilon(bq: int) -> float:
    epsilon = 10 ** (-bq / 10)
    return epsilon


def get_log10_epsilon(bq: int) -> float:
    return math.log10(get_epsilon(bq))


def get_one_minus_epsilon(bq: int) -> float:
    return 1 - get_epsilon(bq)


def get_one_half_minus_epsilon(bq: int) -> float:
    return 0.5 - get_epsilon(bq) / 2.0


def get_log10_one_minus_epsilon(bq: int) -> float:
    return math.log10(get_one_minus_epsilon(bq))


def get_log10_one_half_minus_epsilon(bq: int) -> float:
    return math.log10(get_one_half_minus_epsilon(bq))


def get_germ_gt_state(
    b1: str,
    b2: str,
    ref: str,
) -> Tuple[str, float]:

    gt_state = 0
    if b1 == b2 == ref:  # homref
        gt_state = "homref"
    elif (b1 == ref and b2 != ref) or (b1 != ref and b2 == ref):  # het
        gt_state = "het"
    elif b1 != ref and b2 != ref and b1 != b2:  # hetalt
        gt_state = "hetalt"
    elif b1 != ref and b2 != ref and b1 == b2:  # homalt
        gt_state = "homalt"
    return gt_state


def get_log10_germ_gt_prior(gt_state: str) -> Tuple[str, float]:
    gt_prior = math.log10(gt_state2gt_prior[gt_state])
    return gt_prior


def get_log10_gt_pD(
    gt: str,
    ref: str,
    allele2bq_lst: Dict[int, List[int]],
) -> float:

    gt_pD = 0
    b1, b2 = list(gt)
    gt_state = get_germ_gt_state(b1, b2, ref)
    gt_prior = get_log10_germ_gt_prior(gt_state)
    for base in base_lst:
        base_bq_lst = allele2bq_lst[base2idx[base]]
        if (b1 == b2) and (base == b1 or base == b2):  ## hom
            gt_pD += sum(
                [get_log10_one_minus_epsilon(base_bq) for base_bq in base_bq_lst]
            )
        elif (b1 != b2) and (base == b1 or base == b2):  # het
            gt_pD += sum(
                [get_log10_one_half_minus_epsilon(base_bq) for base_bq in base_bq_lst]
            )
        else:  ## error
            gt_pD += sum([get_log10_epsilon(base_bq / 3) for base_bq in base_bq_lst])
    gt_pD += gt_prior
    gt_pl = -10 * gt_pD
    return gt_pl, gt_state


def get_germ_gt_pD(
    ref: str,
    allele2bq_lst: Dict[int, List[int]],
) -> Tuple[List[str], List[float], str]:

    gt_pl_lst = []
    gt2gt_state = {}
    for gt in gt_lst:
        gt_pl, gt_state = get_log10_gt_pD(gt, ref, allele2bq_lst)
        gt_pl_lst.append(gt_pl)
        gt2gt_state[gt] = gt_state
    return gt_pl_lst, gt2gt_state


def get_argmin_gt(gt_pl_lst):
    ilst = np.argsort(gt_pl_lst)
    gt = [gt_lst[i] for i in ilst][0]
    gt_pl_lst = [gt_pl_lst[i] for i in ilst]
    gq = (np.array(gt_pl_lst) - min(gt_pl_lst))[1]
    gq = int(gq) if gq < 99 else 99
    return gt, gq


def get_germ_gt(
    ref: str,
    allele2bq_lst: Dict[int, List[int]],
) -> Tuple[List[str], List[float], str]:

    gt_pl_lst, gt2gt_state = get_germ_gt_pD(
        ref,
        allele2bq_lst,
    )
    gt, gq = get_argmin_gt(gt_pl_lst)
    gt_state = gt2gt_state[gt]
    if gt[0] != ref and gt.count(ref) == 1:
        gt = gt[::-1]
    return gt, gq, gt_state


def init_allelecounts():
    rpos2allelecounts = defaultdict(lambda: np.zeros(6))
    rpos2allele2pq_lst = defaultdict(lambda: {0: [], 1: [], 2: [], 3: []})
    return rpos2allelecounts, rpos2allele2pq_lst


def update_allelecounts(
    ccs: BAM,
    rpos2allelecounts: Dict[int, np.ndarray],
    rpos2allele2pq_lst: Dict[int, Dict[int, List[int]]],
):

    tpos = ccs.tstart
    qpos = ccs.qstart
    for (state, ref, alt, ref_len, alt_len) in ccs.cstuple_lst:
        if state == 1:  # match
            for i, alt_base in enumerate(alt):
                epos = tpos + i
                bidx = base2idx[alt_base]
                rpos2allelecounts[epos][bidx] += 1
                rpos2allele2pq_lst[epos][bidx].append(ccs.pq_int_lst[qpos + i])
        elif state == 2:  # sub
            bidx = base2idx[alt]
            rpos2allelecounts[tpos][bidx] += 1
            rpos2allele2pq_lst[tpos][bidx].append(ccs.pq_int_lst[qpos])
        elif state == 3:  # insertion
            rpos2allelecounts[tpos][4] += 1
        elif state == 4:  # deletion
            for j in range(len(ref[1:])):
                rpos2allelecounts[tpos + j][5] += 1
        tpos += ref_len
        qpos += alt_len


def update_filtered_allelecounts(
    ccs: BAM,
    min_trim: float,
    mismatch_window: int,
    max_mismatch_count: int,
    rpos2allele2pq_lst: Dict[int, Dict[int, List[int]]],
):

    ccs.cs2subindel()
    rpos = ccs.tstart
    qpos = ccs.qstart
    trimmed_qstart, trimmed_qend = ccs.get_trimmed_range(min_trim)  ## done
    for cstuple in ccs.cstuple_lst:
        state, ref, alt, ref_len, alt_len = cstuple
        if state == 1:  # match
            mismatch_tstart, mismatch_tend = ccs.get_mismatch_range(
                rpos, qpos, ccs.qlen, mismatch_window
            )
            for j, ref_base in enumerate(ref):
                if ccs.is_trimmed_base(qpos + j, trimmed_qstart, trimmed_qend):
                    continue
                mismatch_count = bisect.bisect_right(
                    ccs.mismatch_tpos_lst, mismatch_tend + j
                ) - bisect.bisect_left(ccs.mismatch_tpos_lst, mismatch_tstart + j)
                if ccs.is_mismatch_base(mismatch_count, max_mismatch_count):
                    continue
                rpos2allele2pq_lst[rpos + j][base2idx[ref_base]].append(
                    ccs.pq_int_lst[qpos + j]
                )
        elif state == 2:
            mismatch_tstart, mismatch_tend = ccs.get_mismatch_range(
                rpos, qpos, ccs.qlen, mismatch_window
            )
            mismatch_count = (
                bisect.bisect_right(ccs.mismatch_tpos_lst, mismatch_tend)
                - bisect.bisect_left(ccs.mismatch_tpos_lst, mismatch_tstart)
                - 1
            )
            if ccs.is_trimmed_base(qpos, trimmed_qstart, trimmed_qend):
                pass
            if ccs.is_mismatch_base(mismatch_count, max_mismatch_count):
                pass
            rpos2allele2pq_lst[rpos][base2idx[alt]].append(ccs.pq_int_lst[qpos])
        rpos += ref_len
        qpos += alt_len


def get_base_counts(
    allelecounts: Dict[int, Dict[int, int]],
):
    ins_count = allelecounts[4]
    del_count = allelecounts[5]
    # base_sum = sum(allelecounts[0:4])
    read_depth = sum(allelecounts) - ins_count
    return del_count, ins_count, read_depth
    # return base_sum, del_count, ins_count, read_depth


def is_low_mapq(ccs_mapq: int, min_mapq: int):
    if ccs_mapq < min_mapq:
        return True
    else:
        return False


def is_outlier(qlen, qlen_lower_limit, qlen_upper_limit):
    if qlen_lower_limit < qlen and qlen < qlen_upper_limit:
        return False
    else:
        return True


def is_low_gq(
    germ_gq: float,
    min_gq: int,
) -> bool:
    if germ_gq < min_gq:
        return True
    else:
        return False


def get_tri(
    seq: str,
    pos: int,
):

    tri = seq[pos - 1 : pos + 2]
    if len(tri) == 3 and tri.count("N") == 0:
        if tri[1] in purine:
            tri_pyr = "".join([purine2pyrimidine.get(base, "N") for base in tri[::-1]])
            return tri_pyr
        return tri
    return "NNN"


def get_sbs6(
    ref: str,
    alt: str,
):

    if ref in purine:
        sbs6 = "{}>{}".format(
            purine2pyrimidine.get(ref, "N"),
            purine2pyrimidine.get(alt, "N"),
        )
    else:
        sbs6 = "{}>{}".format(ref, alt)
    return sbs6


def get_sbs96(
    ubase: str,
    ref: str,
    alt: str,
    dbase: str,
):

    if ref in purine:
        sbs96 = "{}[{}>{}]{}".format(
            ubase,
            purine2pyrimidine.get(ref, "N"),
            purine2pyrimidine.get(alt, "N"),
            dbase,
        )
    else:
        sbs96 = "{}[{}>{}]{}".format(ubase, ref, alt, dbase)

    if sbs96.count("N") == 0:
        return sbs96
    return "N[N>N]N"


def add_match_count(
    tri: str,
    hom_base: str,
    allele2pq_lst: Dict[int, List[int]],
    pq2tri_count: Dict[int, Dict[int, str]],
    pq2match_mismatch_count: Dict[int, List[int]],
):
    for pq in allele2pq_lst[base2idx[hom_base]]:
        pq2tri_count[pq][tri2idx[tri]] += 1
        pq2match_mismatch_count[pq][0] += 1


def add_mismatch_count(
    sub_sbs96_lst: List[str],
    allele2pq_lst: Dict[int, List[int]],
    pq2sbs96_count: Dict[int, Dict[int, int]],
    pq2match_mismatch_count: Dict[int, List[int]],
):
    for (sub, sbs96) in sub_sbs96_lst:
        for pq in allele2pq_lst[base2idx[sub]]:
            pq2match_mismatch_count[pq][1] += 1
            pq2sbs96_count[pq][sbs96_to_idx[sbs96]] += 1


def get_match_mismatch_count(
    seq: str,
    chrom: str,
    bam_file: str,
    chunkloci_list: List[Tuple[str, int, int]],
    min_gq: int,
    min_mapq: int,
    min_trim: float,
    md_threshold: float,
    qlen_lower_limit: int,
    qlen_upper_limit: int,
    mismatch_window: int,
    max_mismatch_count: int,
    germline_snv_prior: float,
    chrom2pq2tri_count: Dict[str, Dict[int, int]],
    chrom2pq2sbs96_count: Dict[str, Dict[int, int]],
    chrom2pq2match_mismatch_count: Dict[str, Dict[int, List[int]]],
    chrom2pq2filtered_tri_count: Dict[str, Dict[int, int]],
    chrom2pq2filtered_sbs96_count: Dict[str, Dict[int, int]],
    chrom2pq2filtered_match_mismatch_count: Dict[str, Dict[int, List[int]]],
):

    # counter = 0
    gt_init(germline_snv_prior)
    alignments = pysam.AlignmentFile(bam_file, "rb")
    pq2tri_count = {pq: np.zeros(len(tri_lst)) for pq in range(1, 94)}
    pq2sbs96_count = {pq: np.zeros(len(sbs96_lst)) for pq in range(1, 94)}
    pq2match_mismatch_count = {pq: np.zeros(2) for pq in range(1, 94)}
    pq2filtered_tri_count = {pq: np.zeros(len(tri_lst)) for pq in range(1, 94)}
    pq2filtered_sbs96_count = {pq: np.zeros(len(sbs96_lst)) for pq in range(1, 94)}
    pq2filtered_match_mismatch_count = {pq: np.zeros(2) for pq in range(1, 94)}
    for (chrom, chunk_start, chunk_end) in chunkloci_list:
        _, rpos2filtered_allele2pq_lst = init_allelecounts()
        rpos2allelecounts, rpos2allele2pq_lst = init_allelecounts()
        for i in alignments.fetch(chrom, chunk_start, chunk_end):  # iterate through reads
            ccs = BAM(i)
            if not ccs.is_primary:
                continue
            update_allelecounts(ccs, rpos2allelecounts, rpos2allele2pq_lst)
            if is_low_mapq(ccs.mapq, min_mapq):
                continue
            if is_outlier(ccs.qlen, qlen_lower_limit, qlen_upper_limit):
                continue
            update_filtered_allelecounts(
                ccs,
                min_trim,
                mismatch_window,
                max_mismatch_count,
                rpos2filtered_allele2pq_lst,
            )

        for rpos in range(chunk_start, chunk_end): 
            ref = seq[rpos]
            if ref not in base_set:
                continue

            tri = get_tri(seq, rpos)
            ubase = tri[0]
            dbase = tri[-1]
            allelecounts = rpos2allelecounts[rpos]
            allele2pq_lst = rpos2allele2pq_lst[rpos]
            filtered_allele2pq_lst = rpos2filtered_allele2pq_lst[rpos]
            germ_gt, germ_gq, germ_gt_state = get_germ_gt(ref, allele2pq_lst)
            if is_low_gq(germ_gq, min_gq):  ## proportional to sequence coverage
                continue
            if germ_gt_state != "homref":  # not homozygous reference
                continue

            sub_sbs96_lst = [
                (sub, get_sbs96(ubase, ref, sub, dbase))
                for sub in base_set.difference(ref)
                if allelecounts[base2idx[sub]] == 1
            ]
            add_match_count(
                tri, ref, allele2pq_lst, pq2tri_count, pq2match_mismatch_count
            )
            add_mismatch_count(
                sub_sbs96_lst,
                allele2pq_lst,
                pq2sbs96_count,
                pq2match_mismatch_count,
            )
                
            del_count, ins_count, read_depth = get_base_counts(allelecounts)
            if read_depth > md_threshold:  # high-depth
                continue
            if del_count > 0 or ins_count > 0:  # indel present
                continue

            add_match_count(
                tri,
                ref,
                filtered_allele2pq_lst,
                pq2filtered_tri_count,
                pq2filtered_match_mismatch_count,
            )
            add_mismatch_count(
                sub_sbs96_lst,
                filtered_allele2pq_lst,
                pq2filtered_sbs96_count,
                pq2filtered_match_mismatch_count,
            )
        # counter += 1
        # if counter > 0:
        #     break

    chrom2pq2tri_count[chrom] = pq2tri_count
    chrom2pq2sbs96_count[chrom] = pq2sbs96_count
    chrom2pq2match_mismatch_count[chrom] = pq2match_mismatch_count
    chrom2pq2filtered_tri_count[chrom] = pq2filtered_tri_count
    chrom2pq2filtered_sbs96_count[chrom] = pq2filtered_sbs96_count
    chrom2pq2filtered_match_mismatch_count[chrom] = pq2filtered_match_mismatch_count
    alignments.close()


def get_bq(match_count: int, mismatch_count: int):
    bq = -10 * math.log10(mismatch_count / float(match_count))
    return bq


def cumsum_pq2tri_count(
    chrom_lst: List[str],
    chrom2pq2tri_count: Dict[str, Dict[int, int]], 
):

    pq2tri_count = {pq: np.zeros(len(tri_lst)) for pq in range(1, 94)}
    for chrom in chrom_lst: 
        for pq in range(1,94):
            pq2tri_count[pq] += chrom2pq2tri_count[chrom][pq]
    return pq2tri_count


def cumsum_pq2sbs96_count(
    chrom_lst: List[str],
    chrom2pq2sbs96_count: Dict[str, Dict[int, int]]
):
    
    pq2sbs96_count = {pq: np.zeros(len(sbs96_lst)) for pq in range(1, 94)}
    for chrom in chrom_lst: 
        for pq in range(1,94):
            pq2sbs96_count[pq] += chrom2pq2sbs96_count[chrom][pq]
    return pq2sbs96_count
    

def cumsum_pq2match_mismatch_count(
    chrom_lst: List[str],
    chrom2pq2match_mismatch_count: Dict[str, Dict[int, int]],
):

    pq2match_count = defaultdict(lambda: 0)
    pq2mismatch_count = defaultdict(lambda: 0)
    for chrom in chrom_lst:
        for pq in range(1, 94):
            pq2match_count[pq] += chrom2pq2match_mismatch_count[chrom][pq][0]
            pq2mismatch_count[pq] += chrom2pq2match_mismatch_count[chrom][pq][1]
    return pq2match_count, pq2mismatch_count


def dump_pq2bq(
    chrom_lst: List[str],
    chrom2pq2match_mismatch_count: Dict[str, Dict[int, List[int]]],
    out_file: str,
):

    o = open(out_file, "w")
    o.write("{}\t{}\t{}\t{}\n".format("pq", "bq", "mismatch", "match"))
    pq2match_count, pq2mismatch_count = cumsum_pq2match_mismatch_count(
        chrom_lst, 
        chrom2pq2match_mismatch_count
    )
    for pq in range(1, 94):
        match_count = pq2match_count[pq]
        mismatch_count = pq2mismatch_count[pq]
        if match_count != 0 and mismatch_count != 0:
            bq = get_bq(match_count, mismatch_count)
        else:
            bq = "."
        o.write("{}\t{}\t{}\t{}\n".format(pq, bq, mismatch_count, match_count))
    o.close()


def dump_pq2sbs96(
    chrom_lst: List[str],
    chrom2pq2tri_count: Dict[str, Dict[int, List[int]]],
    chrom2pq2sbs96_count: Dict[str, Dict[int, List[int]]],
    out_file: str,
): 

    o = open(out_file, "w")
    o.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
        "pq", 
        "bq", 
        "sbs96", 
        "sbs96_count",
        "tri_count",
        "sub", 
        "tri", 
        "ubase", 
        "dbase"
        )
    )
    pq2tri_count = cumsum_pq2tri_count(
        chrom_lst, 
        chrom2pq2tri_count
    ) 
    pq2sbs96_count = cumsum_pq2sbs96_count(
        chrom_lst, 
        chrom2pq2sbs96_count
    )
    for pq in range(1, 94): 
        for sbs96 in sbs96_lst:
            if sbs96.count("N") != 0:
                continue
            sub = sbs96_to_sub[sbs96] 
            tri = sbs96_to_tri[sbs96] 
            tri_count = pq2tri_count[pq][tri2idx[tri]]
            sbs96_count = pq2sbs96_count[pq][sbs96_to_idx[sbs96]]
            ubase, _, dbase = list(tri)
            if tri_count != 0 and sbs96_count != 0:
                sbs96_bq = get_bq(tri_count, sbs96_count)
            else:
                sbs96_bq = "."

            o.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                pq,
                sbs96_bq,
                sbs96,
                int(sbs96_count),
                int(tri_count),
                sub,
                tri,
                ubase,
                dbase, 
                )
            ) 


def collect_counts(
    bam_file: str,
    ref_file: str,
    region: str,
    region_list: str,
    min_gq: int,
    min_mapq: int,
    min_trim: float,
    mismatch_window: int,
    max_mismatch_count: int,
    germline_snv_prior: float,
    threads: int,
    prefix: str,
):
    cpu_start = time.time() / 60
    print(
        "counting the number of matches and mismatches with {} threads".format(threads)
    )
    p = mp.Pool(threads)
    manager = mp.Manager()
    ref_seq = pyfastx.Fasta(ref_file)
    tname2tsize = get_tname2tsize(bam_file)
    chrom2pq2tri_count = manager.dict()
    chrom2pq2sbs96_count = manager.dict()
    chrom2pq2match_mismatch_count = manager.dict()
    chrom2pq2filtered_tri_count = manager.dict()
    chrom2pq2filtered_sbs96_count = manager.dict()
    chrom2pq2filtered_match_mismatch_count = manager.dict()
    chrom_lst, chrom2chunkloci_lst = load_loci(region, region_list, tname2tsize)
    read_depth, qlen_lower_limit, qlen_upper_limit = get_bam_thresholds(
        bam_file, chrom_lst, tname2tsize
    )
    md_threshold = get_md_threshold(read_depth)
    match_mismatch_count_arg_lst = [
        (
            str(ref_seq[chrom]),
            chrom,
            bam_file,
            chrom2chunkloci_lst[chrom],
            min_gq,
            min_mapq,
            min_trim,
            md_threshold,
            qlen_lower_limit,
            qlen_upper_limit,
            mismatch_window,
            max_mismatch_count,
            germline_snv_prior,
            chrom2pq2tri_count,
            chrom2pq2sbs96_count,
            chrom2pq2match_mismatch_count,
            chrom2pq2filtered_tri_count,
            chrom2pq2filtered_sbs96_count,
            chrom2pq2filtered_match_mismatch_count,
        )
        for chrom in chrom_lst
    ]
    p.starmap(
        get_match_mismatch_count,
        match_mismatch_count_arg_lst,
    )
    p.close()
    p.join()

    dump_pq2bq(chrom_lst, chrom2pq2match_mismatch_count, "{}.pq2bq.txt".format(prefix))
    dump_pq2bq(chrom_lst, chrom2pq2filtered_match_mismatch_count, "{}.filtered_pq2bq.txt".format(prefix))
    dump_pq2sbs96(chrom_lst, chrom2pq2tri_count, chrom2pq2sbs96_count, "{}.pq2sbs96.txt".format(prefix))
    dump_pq2sbs96(chrom_lst, chrom2pq2filtered_tri_count, chrom2pq2filtered_sbs96_count, "{}.filtered_pq2sbs96.txt".format(prefix))

    cpu_end = time.time() / 60
    duration = cpu_end - cpu_start
    print("pq2bq took {} minutes".format(duration))


def main():
    options = parse_args(sys.argv)
    collect_counts(
        options.bam,
        options.ref,
        options.region,
        options.region_list,
        options.min_gq,
        options.min_mapq,
        options.min_trim,
        options.mismatch_window,
        options.max_mismatch_count,
        options.germline_snv_prior,
        options.threads,
        options.prefix,
    )
    sys.exit(0)


if __name__ == "__main__":
    main()
