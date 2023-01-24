import pyfastx
import himut.mutlib
import multiprocessing as mp
from collections import defaultdict
from typing import Dict, List, Tuple
from himut.mutlib import purine, purine2pyrimidine


def load_ref_tricounts(
    chrom: str,
    chrom_seq: str,
    chunkloci_lst: List[Tuple[str, int, int]],
    chrom2tri2count: Dict[str, Dict[str, int]],
) -> Dict[str, Dict[str, int]]:

    tri2count = defaultdict(lambda: 0)
    for (_chrom, chunk_start, chunk_end) in chunkloci_lst:
        for i in range(chunk_start, chunk_end):
            tri = chrom_seq[i-1:i+2]
            if tri[1] in purine:
                tri_pyr = "".join(
                    [purine2pyrimidine.get(base, "N") for base in tri[::-1]]
                )
                tri2count[tri_pyr] += 1
            else:
                tri2count[tri] += 1
    chrom2tri2count[chrom] = dict(tri2count)


def get_ref_tricounts(
    refseq: str, chrom2chunkloci_lst: List[str], threads: int
) -> Dict[str, Dict[str, int]]:

    p = mp.Pool(threads)
    manager = mp.Manager()
    chrom2tri2count = manager.dict()
    load_ref_tricount_arg_lst = [
        (chrom, str(refseq[chrom]), chrom2chunkloci_lst[chrom], chrom2tri2count)
        for chrom in chrom2chunkloci_lst
    ]
    p.starmap(load_ref_tricounts, load_ref_tricount_arg_lst)
    p.close()
    p.join()
    return chrom2tri2count

