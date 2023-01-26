import bisect
import himut.bamlib
from collections import defaultdict
from typing import Set, Dict, List, Tuple
bit_complement_hsh = {"0": "1", "1": "0", "-": "-"}


def phase_set2chunkloci_lst(chrom: int, phase_set2hpos_lst: Dict[str, List[int]]):
    chunkloci_lst = [
        (chrom, hpos_lst[0], hpos_lst[-1]) for hpos_lst in phase_set2hpos_lst.values()
    ]
    return chunkloci_lst


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


def get_phase_set(
    chunk_start: int,
    chunk_end: int,
    hpos_lst: List[int],
    hpos2phase_set: Dict[int, str],
):
    idx = bisect.bisect_right(hpos_lst, chunk_start)
    jdx = bisect.bisect_right(hpos_lst, chunk_end)
    hpos_subset_lst = [hpos_lst[kdx] for kdx in range(idx, jdx)]
    phase_set_lst = list(set([hpos2phase_set[hpos] for hpos in hpos_subset_lst]))
    if len(phase_set_lst) == 1:
        phase_set = phase_set_lst[0]
    else:
        phase_set = "."
    return phase_set


def get_ccs_hbit(ccs, hetsnp_subset_lst: List[Tuple[str, int, str, str]]) -> List[str]:

    ccs_hbit = ""
    ccs.cs2tpos2qbase()
    for hetsnp in hetsnp_subset_lst:  # get read haplotype bits
        qbase = ccs.tpos2qbase[hetsnp[0]][0]
        if qbase == hetsnp[1]:  # ref
            ccs_hbit += "0"
        elif qbase == hetsnp[2]:  # alt
            ccs_hbit += "1"
        else:
            ccs_hbit += "-"
    return ccs_hbit


def get_ccs_hap(
    ccs,
    hbit_lst: Dict[str, List[str]],
    hpos_lst: Dict[str, List[int]],
    hetsnp_lst: List[Tuple[int, str, str]],
):

    idx = bisect.bisect_right(hpos_lst, ccs.tstart)
    jdx = bisect.bisect_right(hpos_lst, ccs.tend)
    if (jdx - idx) < 2:
        return "."
    else:
        hetsnp_subset_lst = [hetsnp_lst[kdx] for kdx in range(idx, jdx)]
        h0_hbit = "".join([hbit_lst[kdx] for kdx in range(idx, jdx)])
        h1_hbit = "".join([bit_complement_hsh[hbit] for hbit in h0_hbit])
        ccs_hbit = get_ccs_hbit(ccs, hetsnp_subset_lst)
        if h0_hbit == ccs_hbit:
            ccs_hap = "0"
        elif h1_hbit == ccs_hbit:
            ccs_hap = "1"
        else:
            ccs_hap = "."
        return ccs_hap

def get_tpos_hap(
    alignments,
    loci: Tuple[str, int, int],
    hbit_lst: List[str],
    hpos_lst: List[int],
    hetsnp_lst: List[Tuple[int, str, str]],
) -> bool:

    ccs_hap_lst = []
    for i in alignments.fetch(*loci):
        ccs = himut.bamlib.BAM(i)
        ccs_hap = himut.haplib.get_ccs_hap(
            ccs,
            hbit_lst,
            hpos_lst,
            hetsnp_lst,
        )
        ccs_hap_lst.append((ccs.qname, ccs_hap))
    return ccs_hap_lst
