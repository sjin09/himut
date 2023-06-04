#!/usr/bin/env python3

import sys
import pyfastx
import argparse
import statistics


def parse_args(args):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="BAM file to read",
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


def dump_fastx_statistics(
    infile: str,
    outfile: str
): 
    
    hbq_sum = 0
    qlen_lst = []
    if infile.endswith((".fasta", ".fasta.gz")):
        for (name, seq, comment) in pyfastx.Fastx(infile):
            qlen = len(seq)
            qlen_lst.append(qlen) 
    elif infile.endswith((".fastq", ".fastq.gz")):
        for i in pyfastx.Fastx(infile):
            qlen = len(i.seq)
            qlen_lst.append(qlen)
            hbq_sum += i.qual.count("~")
    seq_sum = sum(qlen_lst) 
    seq_cnt = len(qlen_lst)
    qlen_min = min(qlen_lst)
    qlen_max = max(qlen_lst)
    qlen_std = statistics.stdev(qlen_lst)
    qlen_mean = seq_sum/float(seq_cnt)
    hbq_proportion = (hbq_sum/float(seq_sum))*100

    o = open(outfile, "w")
    o.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2f}%\n".format(infile, seq_cnt, qlen_min, qlen_mean, qlen_std, qlen_max, seq_sum, hbq_proportion))
    o.close()
  


def main():
    options = parse_args(sys.argv)
    dump_fastx_statistics(
        options.input, 
        options.output
    )
    sys.exit(0)


if __name__ == "__main__":
    main()
