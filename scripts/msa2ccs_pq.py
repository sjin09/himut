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
        help="abpoa multiple sequence alignment"
    )
    parser.add_argument(
        "--min_subread_count",
        type=int,
        default=8,
        required=False,
        help="abpoa multiple sequence alignment"
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


def msa2ccs(infile, min_subread_cnt, outfile):

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
    seq_count = len(subread_lst)
    if seq_count >= min_subread_cnt:
        for x, y in enumerate(ccs_seq):
            if y == "-":
                continue

            base2count = defaultdict(lambda: 0) 
            for subread in subread_lst:
                base2count[subread[x]] += 1
            ccs_base_count = base2count[y]
            if seq_count == ccs_base_count:
                bq_lst.append("~") # 93
            else:
                bq_lst.append("!") # 1
        ccs_bq = "".join(bq_lst)
    else:
        ccs_bq = "!" * len(ccs_seq)

    if not outfile.endswith(".gz"):
        outfile = "{}.gz".format(outfile)
        
    if outfile.endswith(".gz"):
        o = gzip.open(outfile, "wb")
        o.write("@{}\n{}\n+\n{}\n".format(ccs_header, ccs_seq.replace("-", ""),  ccs_bq).encode('utf-8'))
        o.close()


def main():
    options = parse_args(sys.argv)
    msa2ccs(options.input, options.min_subread_count, options.output)
    sys.exit(0)


if __name__ == "__main__":
    main()

