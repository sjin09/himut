#!/usr/bin/python

import sys
import gzip
import pyfastx
import argparse
import pyabpoa as poa
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
        help="subreads_to_ccs.fasta file to read"
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="file to return CCS fastq"
    )
    args = args[1:]
    return parser.parse_args(args)


def seq2msa(seq_file):

    aln = poa.msa_aligner()
    seq_lst = [i.seq for i in seq_file]
    msa = aln.msa(seq_lst, out_cons=False, out_msa=True)
    return msa.msa_seq


def msa2ccs_bq(msa_lst):


    ccs_msa = msa_lst[0]
    idx2base2count = {i: defaultdict(lambda: 0) for i in range(len(ccs_msa))}
    for subread in msa_lst[1:]: # subreads
        for i, subread_base in enumerate(subread):
            idx2base2count[i][subread_base] += 1
   
    ccs_bq_lst = []
    subread_count = len(msa_lst[1:])
    for i, j in enumerate(ccs_msa):
        if j == "-":
            ccs_bq_lst.append("-")
            continue
        subread_base_count = idx2base2count[i][j] 
        if subread_count == subread_base_count:
            ccs_bq_lst.append("~")
        else:
            ccs_bq_lst.append("!")
    ccs_bq = "".join(ccs_bq_lst).replace("-", "")
    return ccs_bq


def dump_fastq(ccs_header, ccs_seq, ccs_bq, outfile):

    o = open(outfile, "w") 
    o.write("@{}\n{}\n+\n{}\n".format(ccs_header, ccs_seq, ccs_bq))
    o.close()


def seq2ccs(infile, outfile):

    seq_file = pyfastx.Fasta(infile)
    msa_lst = seq2msa(seq_file)
    ccs_header = seq_file[0].name
    ccs_seq = seq_file[0].seq
    ccs_bq = msa2ccs_bq(msa_lst)
    dump_fastq(ccs_header, ccs_seq, ccs_bq, outfile)


def main():

    options = parse_args(sys.argv)
    seq2ccs(options.input, options.output)
    sys.exit(0)


if __name__ == "__main__":
    main()


