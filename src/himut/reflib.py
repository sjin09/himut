import pyfastx
import natsort
import himut.util
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
    
    tri2count = defaultdict(lambda: 0)
    for chrom in chrom_lst:
        for tri in tri_lst:
            tri2count[tri] += chrom2tri2count[chrom][tri]        
    return dict(tri2count)


def get_ref_tricount(
    ref_file: str, 
    region: str,
    region_list: str,
    threads: int,
    out_file: str
) -> Dict[str, Dict[str, int]]:

    chrom_lst = []
    refseq = pyfastx.Fasta(ref_file)
    if region is None and region_list is not None:
        for line in open(region_list).readlines():
            chrom = line.strip()
            chrom_lst.append(chrom)
    elif region is not None and region_list is None:
        chrom_lst.append(region)
    elif region is not None and region_list is not None:
        for line in open(region_list).readlines():
            chrom = line.strip()
            chrom_lst.append(chrom)
    else:
        print("Please provide --region or --region_list")
        himut.util.exit()
   
    o = open(out_file, "w")  
    tri2count = get_genome_tricounts(refseq, chrom_lst, threads) 
    tri_lst = natsort.natsorted(list(tri2count.keys()))
    for tri in tri_lst:
        o.write("{}\t{}\n".format(tri, tri2count[tri]))
    o.close()
