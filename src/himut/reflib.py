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
            tri = chrom_seq[i : i + 3]
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


def load_seq_loci(
    ref_file: str, chrom_lst: List[str]
) -> Dict[str, List[Tuple[str, int, int]]]:

    refseq = pyfastx.Fasta(ref_file)
    chrom2loci_lst = defaultdict(list)
    for chrom in chrom_lst:
        status = 0
        chrom_seq = str(refseq[chrom])
        chrom_len = len(chrom_seq) - 1
        for i, j in enumerate(chrom_seq):
            if status == 0:
                if j == "N":
                    continue
                else:
                    start = i
                    status = 1
            elif status == 1:
                if j == "N":
                    end = i
                    status = 0
                    chrom2loci_lst[chrom].append((chrom, start, end))
                else:
                    if i == chrom_len:
                        end = i - 1
                        status = 0
                        chrom2loci_lst[chrom].append((chrom, start, end))
                    else:
                        continue
    chrom2chunkloci_lst = {}
    for chrom in chrom_lst:
        chunkloci_lst = [
            chunkloci
            for loci in chrom2loci_lst[chrom]
            for chunkloci in himut.util.chunkloci(loci)
        ]
        chrom2chunkloci_lst[chrom] = chunkloci_lst
    return chrom2chunkloci_lst
