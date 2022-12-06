import json
import bisect
import himut.util
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


def get_closest_haplotype(h0_lst: List[str], hbit_lst: List[str]) -> str: 
    h0_distance, h1_distance = get_hamming_distance(h0_lst, hbit_lst)
    if h0_distance > h1_distance:
        return "1"
    elif h0_distance < h1_distance:
        return "0"
    elif h0_distance == h1_distance:
        return "."


def get_ccs_hbits(
    tpos2qbase: Dict[int, Tuple[str, int]], 
    hetsnp_subset_lst: List[Tuple[str, int, str, str]]
) -> List[str]:
    
    hbit_lst = []
    for hetsnp in hetsnp_subset_lst: # get read haplotype bits
        qbase, _ = tpos2qbase[hetsnp[0]]
        if qbase == hetsnp[1]:
            hbit_lst.append("0")
        elif qbase == hetsnp[2]:
            hbit_lst.append("1")
        else:
            hbit_lst.append("-")
    return hbit_lst


def get_ccs_haplotype(
    h0_lst: List[str],
    hpos_lst: List[int], 
    hetsnp_lst: List[Tuple[int, str, str]],
    tpos2qbase: Dict[str, Tuple[str, int]],  
):
    tpos_lst = list(tpos2qbase.keys())
    ref_ccs_start = tpos_lst[0]
    ref_ccs_end = tpos_lst[-1]
    idx = bisect.bisect_right(hpos_lst, ref_ccs_start)
    jdx = bisect.bisect_right(hpos_lst, ref_ccs_end)
    if idx == jdx:
        return "."
    elif jdx - idx == 1:
        return "."
    else:
        h0_subset_lst = []
        hetsnp_subset_lst = []
        kdx_lst = list(range(idx, jdx))
        for kdx in kdx_lst:
            h0_subset_lst.append(h0_lst[kdx])
            hetsnp_subset_lst.append(hetsnp_lst[kdx])
        ccs_hbit_lst = get_ccs_hbits(tpos2qbase, hetsnp_subset_lst)
        ccs_hap = get_closest_haplotype(h0_subset_lst, ccs_hbit_lst)
        return ccs_hap


# def get_read_haplotype(
#     tpos2qbase: Dict[int, Tuple[str, int]],
#     hpos_lst: List[int],
#     hetsnp_lst: List[Tuple[str, int, str, str]],
#     hidx2hstate: Dict[int, str],
#     hetsnp2bidx: Dict[Tuple[str, int, str, str], int],
#     hetsnp2hidx: Dict[Tuple[str, int, str, str], int]
# ) -> str:

#     tpos_lst = list(tpos2qbase.keys())
#     idx = bisect.bisect_right(hpos_lst, tpos_lst[0])
#     jdx = bisect.bisect_right(hpos_lst, tpos_lst[-1])
#     if idx == jdx:
#         return ".", "."
#     elif jdx - idx == 1:
#         return ".", "."
#     else:
#         bidx2hetsnp_subset_lst = defaultdict(list)
#         hetsnp_subset_lst = [hetsnp_lst[kdx] for kdx in range(idx, jdx)]
#         for hetsnp in hetsnp_subset_lst:
#             bidx2hetsnp_subset_lst[hetsnp2bidx[hetsnp]].append(hetsnp)

#         bidx = 0
#         hetsnp_subset_lst = 0 
#         hetsnp_subset_count = 0
#         for bidx_candidate, hetsnp_subset_candidate_lst in bidx2hetsnp_subset_lst.items():
#             hetsnp_candidate_count = len(hetsnp_subset_candidate_lst)
#             if hetsnp_subset_count < hetsnp_candidate_count:
#                 bidx = bidx_candidate
#                 hetsnp_subset_count = hetsnp_candidate_count
#                 hetsnp_subset_lst = hetsnp_subset_candidate_lst
#         hbit_lst = get_ccs_hbits(tpos2qbase, hetsnp_subset_lst)
#         h0_lst = [hidx2hstate[hetsnp2hidx[hetsnp]] for hetsnp in hetsnp_subset_lst]
#         hap = get_haplotype(h0_lst, hbit_lst) 
#         return bidx, hap


# def get_region_hap2count(
#     read_lst: List[str],
#     read2tpos2qbase: Dict[str, Dict[int, Tuple[str, int]]],
#     hblock: List[Tuple[int, str]],
#     hidx2hetsnp: Dict[int, Tuple[str, int, str, str]],
# ) -> Tuple[int, int]:

#     h0_lst = []
#     hetsnp_lst = []  
#     hap2count = defaultdict(lambda: 0)
#     for (hidx, hstate) in hblock:
#         h0_lst.append(hstate)
#         hetsnp_lst.append(hidx2hetsnp[hidx])
        
#     hpos_lst = [hetsnp[0] for hetsnp in hetsnp_lst]
#     for read in read_lst:
#         tpos_lst = list(read2tpos2qbase[read].keys())
#         read_tstart = tpos_lst[0]
#         read_tend = tpos_lst[-1]
#         idx = bisect.bisect_right(hpos_lst, read_tstart)
#         jdx = bisect.bisect_right(hpos_lst, read_tend)
#         if idx == jdx:
#             continue
#         elif jdx - idx == 1:
#             continue

#         h0_subset_lst = []
#         hetsnp_subset_lst = []
#         kdx_lst = list(range(idx, jdx))
#         for kdx in kdx_lst:
#             h0_subset_lst.append(h0_lst[kdx])
#             hetsnp_subset_lst.append(hetsnp_lst[kdx])

#         hbit_lst = get_read_hbits(read2tpos2qbase[read], hetsnp_subset_lst)
#         hap = get_haplotype(h0_subset_lst, hbit_lst)
#         hap2count[hap] += 1
#     return hap2count


def get_bidx(
    start: int,
    end: float,
    hpos_lst: List[int],
    hetsnp_lst: List[int],
    hetsnp2bidx: Dict[Tuple[int, str, str], int]
):
    idx = bisect.bisect_right(hpos_lst, start)
    jdx = bisect.bisect_right(hpos_lst, end)
    if idx == jdx:
        return "."
    elif jdx - idx == 1:
        return "."
    else:
        hetsnp_subset_lst = [hetsnp_lst[kdx] for kdx in range(idx, jdx)]
        bidx_lst = list(set([hetsnp2bidx[hetsnp] for hetsnp in hetsnp_subset_lst]).difference("."))
        if len(bidx_lst) == 1:
            bidx = bidx_lst[0]
        else:
            bidx = "." 
    return bidx


def get_hap_counts(
    ref: str,
    alt: str,
    hblock: List[Tuple[str, int]],
    hidx2hetsnp: Dict[int, List[Tuple[int, str, str]]],
    ccs2tpos2qbase: Dict[str, Dict[int, str]],
    allele2ccs_lst: Dict[int, Dict[str, List[str]]],
):
  
    h0_lst = []
    hpos_lst = []
    hetsnp_lst = []  
    wt_ccs_lst = allele2ccs_lst[himut.util.base2idx[ref]]
    som_ccs_lst = allele2ccs_lst[himut.util.base2idx[alt]]
    for (hidx, hstate) in hblock:
        hetsnp = hidx2hetsnp[hidx]
        hpos_lst.append(hetsnp[0])
        hetsnp_lst.append(hetsnp)
        h0_lst.append(hstate)

    hap2count = defaultdict(lambda: 0)
    for wt_ccs in wt_ccs_lst:
        wt_ccs_hap = get_ccs_haplotype(
           h0_lst,
           hpos_lst, 
           hetsnp_lst,
           ccs2tpos2qbase[wt_ccs],
        )
        hap2count[wt_ccs_hap] += 1

    som_ccs_hap_lst = []
    for som_ccs in som_ccs_lst: 
        som_ccs_hap = get_ccs_haplotype(
           h0_lst,
           hpos_lst, 
           hetsnp_lst,
           ccs2tpos2qbase[som_ccs],
        )
        if som_ccs_hap != ".":
            som_ccs_hap_lst.append(som_ccs_hap)

    h0_count = hap2count["0"]  
    h1_count = hap2count["1"] 
    if len(set(som_ccs_hap_lst)) == 1:
        som_hap = som_ccs_hap_lst[0]
        som_hap_count = som_ccs_hap_lst.count(som_hap) 
        return h0_count, h1_count, som_hap, som_hap_count 
    else:
        return h0_count, h1_count, ".", ".", 
