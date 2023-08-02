#!/usr/bin/env python

import os
import sys
import pyfastx
import argparse


def parse_args(args):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="FASTA file to read"
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="directory to return FASTA files"
    )
    args = args[1:]
    return parser.parse_args(args)


def chunkstring(string):
    chunks = [string[i : i + 60] for i in range(0, len(string), 60)]
    return chunks


def faSplit(fa_file, out_dir):

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    for seq in pyfastx.Fasta(fa_file):
        o = open("{}/{}.fasta".format(out_dir, seq.name), "w")
        o.write(">{}\n".format(seq.id))
        for chunk in chunkstring(seq.seq):
            o.write("{}\n".format(chunk))
        o.close()    


def main():
    options = parse_args(sys.argv)
    faSplit(options.input, options.output)
    sys.exit(0)


if __name__ == "__main__":
    main()


