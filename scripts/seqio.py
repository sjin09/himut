#!/usr/bin/env python3

import sys
import gzip
import argparse
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
        help="FASTQ file to read",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="FASTQ file to write",
    )
    args = args[1:]
    return parser.parse_args(args)


def seqio(
    in_file: str,
    out_file: str
):

    o = open(out_file, "w")
    if in_file.endswith(".gz"):
        for i,j in enumerate(gzip.open(in_file).readlines()):
            k = i%4 
            if k == 0:
                qname = j.strip().decode("utf-8")
            elif k == 1:
                qseq = j.strip().decode("utf-8")
            elif k == 3:
                qual = j.strip().decode("utf-8")
                if qual.count("!") == 0:
                    o.write("{}\n{}\n+\n{}\n".format(qname, qseq, qual))
                else:
                    pq_cat = ""
                    for i, pq_str in enumerate(qual):
                        pq_int = ord(pq_str) - 33
                        if pq_int < 1:
                            pq_cat += '"'
                        else:
                            pq_cat += pq_str
                    o.write("{}\n{}\n+\n{}\n".format(qname, qseq, pq_cat))
    else:
        for i,j in enumerate(open(in_file).readlines()):
            k = i%4 
            if k == 0:
                qname = j.strip()
            elif k == 1:
                qseq = j.strip()
            elif k == 3:
                qual = j.strip()
                if qual.count("!") == 0:
                    o.write("{}\n{}\n+\n{}\n".format(qname, qseq, qual))
                else:
                    pq_cat = ""
                    for i, pq_str in enumerate(qual):
                        pq_int = ord(pq_str) - 33
                        if pq_int < 1:
                            pq_cat += '"'
                        else:
                            pq_cat += pq_str
                    o.write("{}\n{}\n+\n{}\n".format(qname, qseq, pq_cat))
    o.close()


def main():
    options = parse_args(sys.argv)
    seqio(
        options.input, 
        options.output
    )
    sys.exit(0)


if __name__ == "__main__":
    main()
