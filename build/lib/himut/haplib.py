import json
import bisect
from collections import defaultdict
from typing import Dict, List, Tuple


def get_hamming_distance(h0_lst: List[str], hbit_lst: List[str]) -> Tuple[int, int]:

    h0_distance = 0
    h1_distance = 0
    for (i, j) in zip(h0_lst, hbit_lst):
        if j == "-":
            continue
        if i != j:
            h0_distance += 1
        else:
            h1_distance += 1
    return h0_distance, h1_distance


def get_haplotype(h0_lst: List[str], hbit_lst: List[str]) -> str: 
    h0_distance, h1_distance = get_hamming_distance(h0_lst, hbit_lst)
    if h0_distance > h1_distance:
        return "1"
    elif h0_distance < h1_distance:
        return "0"
    elif h0_distance == h1_distance:
        return "."


def get_read_hbits(
    tpos2qbase: Dict[int, Tuple[str, int]], 
    hetsnp_subset_lst: List[Tuple[str, int, str, str]]
) -> List[str]:
    
    hbit_lst = []
    for hetsnp in hetsnp_subset_lst: # get read haplotype bits
        qbase, _ = tpos2qbase[hetsnp[1]]
        if qbase == hetsnp[2]:
            hbit_lst.append("0")
        elif qbase == hetsnp[3]:
            hbit_lst.append("1")
        else:
            hbit_lst.append("-")
    return hbit_lst


def get_read_haplotype(
    tpos2qbase: Dict[int, Tuple[str, int]],
    hpos_lst: List[int],
    hetsnp_lst: List[Tuple[str, int, str, str]],
    hidx2hstate: Dict[int, str],
    hetsnp2bidx: Dict[Tuple[str, int, str, str], int],
    hetsnp2hidx: Dict[Tuple[str, int, str, str], int]
) -> str:

    tpos_lst = list(tpos2qbase.keys())
    idx = bisect.bisect_right(hpos_lst, tpos_lst[0])
    jdx = bisect.bisect_right(hpos_lst, tpos_lst[-1])
    if idx == jdx:
        return ".", "."
    elif jdx - idx == 1:
        return ".", "."
    else:
        bidx2hetsnp_subset_lst = defaultdict(list)
        hetsnp_subset_lst = [hetsnp_lst[kdx] for kdx in range(idx, jdx)]
        for hetsnp in hetsnp_subset_lst:
            bidx2hetsnp_subset_lst[hetsnp2bidx[hetsnp]].append(hetsnp)

        bidx = 0
        hetsnp_subset_lst = 0 
        hetsnp_subset_count = 0
        for bidx_candidate, hetsnp_subset_candidate_lst in bidx2hetsnp_subset_lst.items():
            hetsnp_candidate_count = len(hetsnp_subset_candidate_lst)
            if hetsnp_subset_count < hetsnp_candidate_count:
                bidx = bidx_candidate
                hetsnp_subset_count = hetsnp_candidate_count
                hetsnp_subset_lst = hetsnp_subset_candidate_lst
        hbit_lst = get_read_hbits(tpos2qbase, hetsnp_subset_lst)
        h0_lst = [hidx2hstate[hetsnp2hidx[hetsnp]] for hetsnp in hetsnp_subset_lst]
        hap = get_haplotype(h0_lst, hbit_lst) 
        return bidx, hap


def get_region_hap2count(
    read_lst: List[str],
    read2tpos2qbase: Dict[str, Dict[int, Tuple[str, int]]],
    hblock: List[Tuple[int, str]],
    hidx2hetsnp: Dict[int, Tuple[str, int, str, str]],
) -> Tuple[int, int]:

    h0_lst = []
    hetsnp_lst = []  
    hap2count = defaultdict(lambda: 0)
    for (hidx, hstate) in hblock:
        h0_lst.append(hstate)
        hetsnp_lst.append(hidx2hetsnp[hidx])
        
    hpos_lst = [hetsnp[1] for hetsnp in hetsnp_lst]
    for read in read_lst:
        tpos_lst = list(read2tpos2qbase[read].keys())
        read_tstart = tpos_lst[0]
        read_tend = tpos_lst[-1]
        idx = bisect.bisect_right(hpos_lst, read_tstart)
        jdx = bisect.bisect_right(hpos_lst, read_tend)
        if idx == jdx:
            continue
        elif jdx - idx == 1:
            continue

        h0_subset_lst = []
        hetsnp_subset_lst = []
        kdx_lst = list(range(idx, jdx))
        for kdx in kdx_lst:
            h0_subset_lst.append(h0_lst[kdx])
            hetsnp_subset_lst.append(hetsnp_lst[kdx])

        hbit_lst = get_read_hbits(read2tpos2qbase[read], hetsnp_subset_lst)
        hap = get_haplotype(h0_subset_lst, hbit_lst)
        hap2count[hap] += 1
    return hap2count
   