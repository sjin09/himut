#!/usr/bin/python

import sys
import gzip
import math
import scipy
import pyfastx
import argparse
import pyabpoa as poa
from scipy.special import binom
from collections import defaultdict
prior_p = 0.1

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
    for i, ccs_base in enumerate(ccs_msa):
        if ccs_base == "-":
            ccs_bq_lst.append("-")
            continue
        subread_base_count = idx2base2count[i][ccs_base] 
        prior_likelihood  = prior_p ** subread_base_count
        if subread_count == subread_base_count:
            ccs_bq_lst.append("~")
        else:
            bq = -10 * math.log10(prior_likelihood * binom(subread_count, subread_base_count))
            bq = math.floor(bq) 
            if bq <= 0:
                ccs_bq_lst.append('"')
            else:
                ccs_bq_lst.append(chr(bq+33))
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


