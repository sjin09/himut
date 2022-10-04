import sys
import pysam
import argparse
from collections import defaultdict


class BAM:
    def __init__(self, line):
        self.qname = line.query_name
        self.zmw = self.qname.split("/")[1]


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


def dump_zmw_counts(
    infile: str,
    outfile: str
):

    cnt_hsh = defaultdict(lambda: 0)
    alignments = pysam.AlignmentFile(infile, "rb", check_sq=False)
    for line in alignments:
        ccs = BAM(line)
        cnt_hsh[ccs.zmw] += 1
        counter += 1

    o = open(outfile, "w")
    o.write("zmw\tcount\n")
    for zmw, cnt in cnt_hsh.items():
        o.write("{}\t{}\n".format(zmw, cnt))
    o.close()


def main():
    options = parse_args(sys.argv)
    dump_zmw_counts(options.input, options.output)
    sys.exit(0)


if __name__ == "__main__":
    main()
