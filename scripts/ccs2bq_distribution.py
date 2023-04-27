
#!/usr/bin/env python3

import sys
import pyfastx
import argparse
import numpy as np
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
        help="FASTQ file to read",
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


def dump_bq_distribution(
    seq_file: str,
    out_file: str
): 

    fqfile = pyfastx.Fastq(seq_file)
    bq2count = defaultdict(lambda: 0)
    for bq in range(94):
        bq2count[bq] = 0
    for ccs in fqfile:
        for ascii_bq in ccs.qual:
            bq = ord(ascii_bq)-33
            bq2count[bq] += 1
    o = open(out_file, "w")
    o.write("{}\t{}\n".format("bq", "count"))
    for bq in range(94):
        o.write("{}\t{}\n".format(bq, bq2count[bq]))   
    o.close()            


def main():
    options = parse_args(sys.argv)
    dump_bq_distribution(
        options.seq,
        options.out
    )
    sys.exit(0)


if __name__ == "__main__":
    main()
