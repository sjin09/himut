import pyfastx
import himut.mutlib
import multiprocessing as mp
from collections import defaultdict
from typing import Dict, List, Tuple
from himut.mutlib import tri_lst, purine, purine2pyrimidine


def get_chrom_tricount(
    chrom: str,
    seq: str,
    chrom2tri2count: Dict[str, Dict[str, int]],
) -> Dict[str, Dict[str, int]]:

    tri2count = defaultdict(lambda: 0)
    for i in range(len(seq)-2):
        base = seq[i]
        if base == "N":
            continue
        tri = seq[i:i+3]
        if tri[1] in purine:
            tri_pyr = "".join(
                [purine2pyrimidine.get(base, "N") for base in tri[::-1]]
            )
            tri2count[tri_pyr] += 1
        else:
            tri2count[tri] += 1
    chrom2tri2count[chrom] = dict(tri2count)


def get_genome_tricounts(
    refseq: str, 
    chrom_lst: List[str], 
    threads: int
) -> Dict[str, Dict[str, int]]:

    p = mp.Pool(threads)
    manager = mp.Manager()
    chrom2tri2count = manager.dict()
    get_chrom_tricount_arg_lst = [
        (
            chrom, 
            str(refseq[chrom]), 
            chrom2tri2count
        )
        for chrom in chrom_lst
    ]
    p.starmap(get_chrom_tricount, get_chrom_tricount_arg_lst)
    p.close()
    p.join()
    
    genome_sum = 0
    tri2count = defaultdict(lambda: 0)
    for chrom in chrom_lst:
        genome_sum += len(refseq[chrom])
        for tri in tri_lst:
            tri2count[tri] += chrom2tri2count[chrom][tri]        
    return genome_sum, dict(tri2count)


