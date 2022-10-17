import sys
import pysam
import argparse
from collections import defaultdict


class BAM:
    def __init__(self, line):
        self.qname = "/".join(line.query_name.split("/")[0:2])


def parse_args(args):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help=".subreads.bam file to read"
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="file to write"
    )
    args = args[1:]
    return parser.parse_args(args)


def dump_qname_counts(
    infile: str,
    outfile: str
):

    cnt_hsh = defaultdict(lambda: 0)
    alignments = pysam.AlignmentFile(infile, "rb", check_sq=False)
    for line in alignments:
        ccs = BAM(line)
        cnt_hsh[ccs.qname] += 1

    o = open(outfile, "w")
    o.write("qname\tsubread_count\n")
    for qname, subread_cnt in cnt_hsh.items():
        o.write("{}\t{}\n".format(qname, subread_cnt))
    o.close()


def main():
    options = parse_args(sys.argv)
    dump_qname_counts(options.input, options.output)
    sys.exit(0)


if __name__ == "__main__":
    main()
