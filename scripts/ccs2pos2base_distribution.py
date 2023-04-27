
#!/usr/bin/env python3

import sys
import pysam
import pyfastx
import natsort
import argparse
import numpy as np
import himut.util
import himut.gtlib
import himut.bamlib
import himut.caller
import multiprocessing as mp
from collections import defaultdict
from typing import Dict, List, Tuple


def parse_args(args):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-i",
        "--seq",
        type=str,
        required=True,
        help="FASTA file to read",
    )
    parser.add_argument(
        "-o",
        "--out",
        type=str,
        required=True,
        help="file to write",
    )
    args = args[1:]
    return parser.parse_args(args)


def dump_pos_base_distribution(
    seq_file: str,
    out_file: str
): 

    base_lst = list("ATGC")
    fafile = pyfastx.Fasta(seq_file)
    rpos2base2count = {(i/100): defaultdict(lambda: 0) for i in range(101)}
    for ccs in fafile:
        qlen = len(ccs.seq)
        for i, base in enumerate(ccs.seq):
            pos = i
            rpos = round(pos/qlen, 2)
            rpos2base2count[rpos][base] += 1

    o = open(out_file, "w")
    o.write("{}\t{}\t{}\t{}\n".format("qpos", "base", "count", "fraction"))
    for rpos in rpos2base2count:
        base_sum = 0
        for base in base_lst:
            count = rpos2base2count[rpos][base]    
            base_sum += count
        for base in base_lst:
            count = rpos2base2count[rpos][base]    
            count_fraction = count/float(base_sum)
            o.write("{}\t{}\t{}\t{}\n".format(rpos, base, count, count_fraction))   
    o.close()            

def main():
    options = parse_args(sys.argv)
    dump_pos_base_distribution(
        options.seq,
        options.out
    )
    sys.exit(0)


if __name__ == "__main__":
    main()
