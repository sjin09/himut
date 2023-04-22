#!/usr/bin/env python3

import sys
import pysam
import argparse


def parse_args(args):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-i",
        "--bam",
        type=str,
        required=True,
        help="subreads_to_ccs.bam file to read",
    )
    parser.add_argument(
        "--threshold",
        type=int,
        default=10000,
        required=True,
        help="number of ZMWs per BAM file",
    )
    args = args[1:]
    return parser.parse_args(args)



class BAM:
    def __init__(self, line):
        # target
        if line.is_secondary:
            self.is_primary = False
        else:
            self.is_primary = True
            self.tname = line.reference_name
            self.zmw = self.tname.split("/")[1]
            self.tstart = line.reference_start
            self.tend = line.reference_end
            # query
            self.qname = line.query_name
            self.qstart = line.query_alignment_start
            self.qend = line.query_alignment_end
            self.qseq = line.query_sequence
            self.qlen = len(self.qseq)
            self.mapq = line.mapping_quality
            self.bq_int_lst = line.query_qualities


def bam_split(
    infile: str,
    threshold: int,
): 

    aln_lst = []
    current_zmw = 0
    zmw_counter = 0
    shard_counter = 0
    prefix = infile.split(".")[0] # outfile prefix
    alignments = pysam.AlignmentFile(infile, "r", check_sq=False) # load BAM file
    for line in alignments:
        i_aln = BAM(line)
        if i_aln.zmw == current_zmw:
            aln_lst.append(line)
        else:
            zmw_counter += 1
            if zmw_counter == threshold:
                o = pysam.AlignmentFile("{}.{}.subreads_to_ccs.bam".format(prefix, shard_counter), "wb", template=alignments) ## dump alignment
                for j_aln in aln_lst:
                    o.write(j_aln)
                o.close()

                aln_lst = [line]
                zmw_counter = 0
                shard_counter += 1
                current_zmw = i_aln.zmw
            else:
                aln_lst.append(line)
                current_zmw = i_aln.zmw

    o = pysam.AlignmentFile("{}.{}.subreads_to_ccs.bam".format(prefix, shard_counter), "wb", template=alignments) ## dump alignment
    for j_aln in aln_lst:
        o.write(j_aln)
    o.close()


def main():
    options = parse_args(sys.argv)
    bam_split(options.bam, options.seq, options.threhsold)
    sys.exit(0)


if __name__ == "__main__":
    main()

