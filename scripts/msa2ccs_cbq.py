import sys
import gzip
import math
import argparse
import scipy.special
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
        help="fasta.fofn (file of file names)"
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="fastq to return"
    )
    args = args[1:]
    return parser.parse_args(args)

def msa2ccs(infile, outfile):

    ccs_seq = None
    ccs_header = None
    subread_lst = []
    for i, line in enumerate(open(infile).readlines()):
        j = i % 2
        if j == 0:
            seq_header = line.strip()
        elif j == 1:
            seq = line.strip()
            if seq_header.endswith("ccs"):
                ccs_seq = seq
                ccs_header = seq_header.replace(">", "")
            else:
                subread_lst.append(seq)

    bq_lst = []
    prior_p = 0.1
    seq_count = len(subread_lst) 
    for x, y in enumerate(ccs_seq):
        if y == "-":
            continue

        base2count = defaultdict(lambda: 0) 
        for subread in subread_lst:
            base2count[subread[x]] += 1
        ref_count = base2count[y]
        prior_likelihood  = prior_p ** ref_count
        if prior_likelihood == 0:
            bq = 93
        else:
            bq = math.floor(-10 * math.log10(prior_likelihood * scipy.special.binom(seq_count, ref_count)))

        if bq >= 93:
            ascii_bq = "~"
        elif bq < 0:
            ascii_bq = "!"
        else:
            ascii_bq = chr(bq+33)
        bq_lst.append(ascii_bq)
    ccs_bq = "".join(bq_lst)

    if not outfile.endswith(".gz"):
        outfile = "{}.gz".format(outfile)
        
    if outfile.endswith(".gz"):
        o = gzip.open(outfile, "wb")
        o.write("@{}\n{}\n+\n{}\n".format(ccs_header, ccs_seq.replace("-", ""),  ccs_bq).encode('utf-8'))
        o.close()


def main():
    options = parse_args(sys.argv)
    msa2ccs(options.input, options.output)
    sys.exit(0)


if __name__ == "__main__":
    main()

