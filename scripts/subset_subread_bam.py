#!/usr/bin/env python3

import re
import sys
import pysam
import natsort
import argparse
import statistics
from collections import defaultdict
from typing import Dict, List, Tuple, Set


def parse_args(args):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="BAM filed to read",
    )
    parser.add_argument(
        "--count",
        type=int,
        required=False,
        help="maximum number of full-length subreads"
    ) 
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="file to write",
    )
    args = args[1:]
    return parser.parse_args(args)


def dump_subreads(o, subread_lst, max_subread_count):

    qlen_lst = [len(subread.query_sequence) for subread in subread_lst]
    full_length_qlen_lst = qlen_lst[1:-1]
    median_len = statistics.median(full_length_qlen_lst)
    upper_len_limit = 1.2 * median_len
    lower_len_limit = 0.8 * median_len

    counter = 0
    for i, qlen in enumerate(full_length_qlen_lst):
        if lower_len_limit < qlen and qlen < upper_len_limit:
            counter += 1 

    full_length_subread_count = len(full_length_qlen_lst)
    if counter == len(full_length_qlen_lst):
        if max_subread_count <= full_length_subread_count:
            for subread in subread_lst[1:1+max_subread_count]:
                o.write(subread)


def subset_subreads(
    bam_file: str,
    max_subread_count: int, 
    out_file: str
): 
    
    state = 0
    alignments = pysam.AlignmentFile(bam_file, "r", check_sq=False)
    o = pysam.AlignmentFile(out_file, "wb", template=alignments)
    for line in alignments:
        qname = "/".join(line.query_name.split("/")[0:2])
        if state == 0:
            state = 1
            subread_lst = [line]
            current_qname = qname
        else:
            if qname != current_qname:
                dump_subreads(o, subread_lst, max_subread_count)
                current_qname = qname
                subread_lst = [line]
                state = 1
            else:
                subread_lst.append(line)
    dump_subreads(o, subread_lst, max_subread_count)


def main():
    options = parse_args(sys.argv)
    subset_subreads(
        options.input, 
        options.count,
        options.output
    )
    sys.exit(0)


if __name__ == "__main__":
    main()
