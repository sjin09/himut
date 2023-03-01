#!/usr/bin/env python

import sys
import natsort
import pyfastx
import argparse
from typing import Dict
import multiprocessing as mp
from collections import defaultdict


def parse_args(args):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="FASTA file to read"
    )
    parser.add_argument(
        "--target",
        type=str,
        required=True,
        help="target chromosomes separated by new line"
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        required=True,
        help="number of threads"
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="file to return trinucleotide sequence context counts"
    )
    args = args[1:]
    return parser.parse_args(args)


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
]
purine = set(["A", "G"])
purine2pyrimidine = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}

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


def ref2tri(seqfile, tgtfile, threads, outfile):

    p = mp.Pool(threads)
    manager = mp.Manager()
    manager = mp.Manager()
    refseq = pyfastx.Fasta(seqfile)
    chrom2tri2count = manager.dict()
    chrom_lst = natsort.natsorted([line.strip() for line in open(tgtfile)])
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

    o = open(outfile, "w")
    for tri in tri_lst:
        o.write("{}\t{}\n".format(tri, tri2count[tri]))
    o.close()

    
def main():
    options = parse_args(sys.argv)
    ref2tri(options.input, options.target, options.threads, options.output) 
    sys.exit(0)


if __name__ == "__main__":
    main()
