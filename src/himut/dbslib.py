import math
import bisect
import himut.util
import numpy as np
from typing import Set, Dict, List, Tuple
from collections import defaultdict, Counter


def get_dbs_candidates(
   read,
   pon_dbs_set: Set[Tuple[int, str, str]],
   trimmed_qstart: int,
   trimmed_qend: int,
   mpos_lst: List[int],
   mismatch_set: Tuple[str, int, str, str],
   mismatch_window: int,
   max_mismatch_count: int,
   tree_of_life_sample: bool,
   create_panel_of_normals: bool
):

    tdbs_candidate_lst = [] 
    for j, (tdbs, qdbs), in enumerate(zip(read.tdbs_lst, read.qdbs_lst)):
        bq_sum = sum(read.qdbs_bq_lst[j])
        if bq_sum != 186:
            continue
        
        if not tree_of_life_sample and not create_panel_of_normals:
            if tdbs in pon_dbs_set:
                continue

        tpos = tdbs[0] 
        qpos = qdbs[0]                   
        if qpos < trimmed_qstart:
            continue
        elif qpos > trimmed_qend:
            continue

        mismatch_start, mismatch_end = himut.util.get_mismatch_range(tpos, qpos, read.qlen, mismatch_window)
        idx = bisect.bisect_left(mpos_lst, mismatch_start)
        jdx = bisect.bisect_right(mpos_lst, mismatch_end)
        if tdbs in mismatch_set:
            mismatch_count = (jdx - idx) - 1 
        else:
            mismatch_count = jdx - idx 

        if mismatch_count > max_mismatch_count:
            continue
        tdbs_candidate_lst.append(tdbs)
    return tdbs_candidate_lst


def get_dbs_allelecounts(
    tpos: int,
    ref: str,
    alt: str,
    read2tpos2qbase: Dict[str, Dict[int,  str]],
    tpos2allelecounts: Dict[int, np.ndarray],
    tpos2qbase2read_lst: Dict[int, Dict[int, List[str]]],
) -> Tuple[int, int, str, float, int, int, int]:
 
    ref_count = 0
    alt_count = 0
    ins_count = tpos2allelecounts[tpos][4] + tpos2allelecounts[tpos+1][4]
    del_count = tpos2allelecounts[tpos][5] + tpos2allelecounts[tpos+1][5] 
    ref_read_lst = list(set(tpos2qbase2read_lst[tpos][himut.util.base2idx[ref[0]]]).intersection(set(tpos2qbase2read_lst[tpos + 1][himut.util.base2idx[ref[1]]])))
    alt_read_lst = list(set(tpos2qbase2read_lst[tpos][himut.util.base2idx[alt[0]]]).intersection(set(tpos2qbase2read_lst[tpos + 1][himut.util.base2idx[alt[1]]])))

    for ref_read in ref_read_lst:
        ref_read_allele = "{}{}".format(read2tpos2qbase[ref_read][tpos][0], read2tpos2qbase[ref_read][tpos + 1][0])
        if ref == ref_read_allele:
            ref_count += 1 
               
    for alt_read in alt_read_lst:
        alt_read_allele = "{}{}".format(read2tpos2qbase[alt_read][tpos][0], read2tpos2qbase[alt_read][tpos + 1][0])
        if alt == alt_read_allele:
            alt_count += 1 

    total_count = ref_count + alt_count + math.ceil(del_count/float(2)) 
    vaf = alt_count / float(total_count)
    return vaf, ref_count, alt_count, ins_count, del_count, total_count, ref_read_lst, alt_read_lst


def get_dbs(
    chrom: str, 
    tdbs_candidate_lst: List[Tuple[str, int, str, str]],
    read2tpos2qbase: Dict[str, Dict[int,  str]],
    tpos2allelecounts: Dict[int, np.ndarray],
    tpos2qbase2read_lst: Dict[int, Dict[int, List[str]]],
    md_threshold: float,
    min_ref_count: int,
    min_alt_count: int,
):

    tdbs_lst = []
    tdbs2count = Counter(tdbs_candidate_lst) 
    for (tpos, ref, alt), _ in tdbs2count.items():
        vaf, ref_count, alt_count, ins_count, del_count, total_count, _ , _ = get_dbs_allelecounts(tpos, ref, alt, read2tpos2qbase, tpos2allelecounts, tpos2qbase2read_lst)
        if del_count != 0 or ins_count != 0:
            continue

        if alt_count != 1:
            continue 
        
        if total_count > md_threshold:
            continue
        
        if ref_count >= min_ref_count and alt_count >= min_alt_count:
            tdbs_lst.append(
                (chrom, tpos, ref, alt, "93.0", total_count, ref_count, alt_count, vaf, ".")
            )
    return tdbs_lst


def get_phased_dbs(
    chrom: str,
    chrom_tdbs_candidate_lst: List[Tuple[int, str, str]],
    read2tpos2qbase: Dict[str, Dict[int,  str]],
    tpos2allelecounts: Dict[int, np.ndarray],
    tpos2qbase2read_lst: Dict[int, Dict[int, List[str]]],
    md_threshold: float,
    min_ref_count: int,
    min_alt_count: int,
    min_hap_count: int,
    hpos_lst: List[int],
    hetsnp_lst: List[Tuple[str, int, str, str]],
    hblock_lst: List[List[Tuple[int, int]]],
    hidx2hstate: Dict[int, str],
    hidx2hetsnp: Dict[int, Tuple[str, int, str, str]],
    hetsnp2bidx: Dict[Tuple[str, int, str, str], int],
    hetsnp2hidx: Dict[Tuple[str, int, str, str], int],
):
    
    tdbs_lst = []
    tdbs2count = Counter(chrom_tdbs_candidate_lst) 
    for (tpos, ref, alt), _ in tdbs2count.items():
        vaf, ref_count, alt_count, ins_count, del_count, total_count, ref_read_lst, alt_read_lst = get_dbs_allelecounts(tpos, ref, alt, read2tpos2qbase, tpos2allelecounts, tpos2qbase2read_lst)
        if del_count != 0 or ins_count != 0:
            continue

        if total_count > md_threshold:
            continue
        
        if ref_count >= min_ref_count and alt_count >= min_alt_count:
            bidx2count = defaultdict(lambda: 0)
            alt_read_hap_count = defaultdict(lambda: 0)
            for alt_read in alt_read_lst:
                bidx, alt_hap = himut.haplib.get_read_haplotype(
                    read2tpos2qbase[alt_read],
                    hpos_lst,
                    hetsnp_lst,
                    hidx2hstate,
                    hetsnp2bidx,
                    hetsnp2hidx,
                )
                bidx2count[bidx] += 1
                alt_read_hap_count[alt_hap] += 1

            if len(bidx2count.keys()) == 1 and len(alt_read_hap_count.keys()) == 1:
                bidx = list(bidx2count.keys())[0]
                if bidx == ".":
                    continue
                
                alt_hap = list(alt_read_hap_count.keys())[0]
                if alt_hap == ".":
                    continue
            else:
                continue

            ref_hap = "1" if alt_hap == "0" else "0"
            hap2count = himut.haplib.get_region_hap2count( 
                ref_read_lst, read2tpos2qbase, hblock_lst[bidx], hidx2hetsnp,
            )
            ref_hap_count = hap2count[ref_hap]
            alt_hap_count = hap2count[alt_hap]  
            if ref_hap_count >= min_hap_count and alt_hap_count >= min_hap_count:
                phase_set = hidx2hetsnp[hblock_lst[bidx][0][0]][0]
                tdbs_lst.append(
                    (chrom, tpos, ref, alt, "93.0", total_count, ref_count, alt_count, vaf, phase_set)
                )               
    return tdbs_lst