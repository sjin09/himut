#!/usr/bin/env python

import os
import sys
import gzip
import pysam
import pyfastx
import argparse
import numpy as np
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}


class BAM:
    def __init__(self, line):
        if line.is_secondary:
            self.is_primary = False
        else:
            # alignment
            self.is_primary = True
            # target
            self.tname = line.reference_name
            # query
            self.flag = line.flag
            self.qname = line.query_name
            self.qseq = line.get_forward_sequence()
            self.qlen = len(self.qseq)
            ## zmw
            self.zmw = self.tname.split("/")[1]

def parse_args(args):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-i",
        "--bam",
        type=str,
        required=True,
        help="BAM file to read subreads to ccs alignments"
    )
    parser.add_argument(
        "--seq",
        type=str,
        required=True,
        help="FASTA file to read ccs sequences"
    )
    parser.add_argument(
        "--blk",
        type=str,
        required=True,
        help="ZMW blacklist"
    )
    parser.add_argument(
        "--min_subread_count",
        type=int,
        default=10,
        required=False,
        help="number of subreads per CCS read"
    )
    parser.add_argument(
        "--subread_lower_len_threshold",
        type=float,
        default=0.8,
        required=False,
        help="subreads below subread_lower_len_threshold times median subread length are discarded"
    )
    parser.add_argument(
        "--subread_upper_len_threshold",
        type=float,
        default=1.2,
        required=False,
        help="subreads above subread_upper_len_threshold times median subread length are discarded"
    )
    parser.add_argument(
        "--fofn",
        type=str,
        required=True,
        help="fofn"
    )
    args = args[1:]
    return parser.parse_args(args)


def load_blk_file(blk_file):
    blk_set = set([line.strip() for line in open(blk_file).readlines()])
    return blk_set


def get_complement(seq):
    seq_complement = "".join([complement[base] for base in seq])
    return seq_complement

def get_reverse_complement(seq):
    seq_rc = "".join([complement[base] for base in seq[::-1]])
    return seq_rc

    
def bam2seq(
    bam_file, 
    seq_file,
    blk_file, 
    fofn_file,
    min_subread_count, 
    lower_subread_len_threshold, 
    upper_subread_len_threshold, 
):

    if not os.path.exists("zmw"):
        os.mkdir("zmw")

    i_state = 0
    j_state = 0
    i_counter = 0
    j_counter = 0
    current_zmw = 0 
    subread_len_lst = []
    fofn = open(fofn_file, "w")
    lower_ccs_len_threshold = 10000
    upper_ccs_len_threshold = 30000
    blk_set = load_blk_file(blk_file) ## blklist
    ccs_seq = pyfastx.Fasta(seq_file) ## ccs FASTA file
    tname2ccs_idx = {j.name: i for i, j in enumerate(ccs_seq)}
    alignments = pysam.AlignmentFile(bam_file, "rb") ## subreads_to_ccs BAM file
    if not os.path.exists("zmw/{}".format(j_counter)):
        os.mkdir("zmw/{}".format(j_counter))

    for line in alignments: # iterate
        i_aln = BAM(line)
        if i_state == 0:
            if not i_aln.is_primary:
                continue
            i_state = 1
            zmw = i_aln.zmw
            aln_lst = [line]
            current_zmw = zmw
            current_tname = i_aln.tname
            subread_len_lst.append(i_aln.qlen)
        elif i_state == 1:
            if current_zmw == i_aln.zmw: ## second to last alignments
                aln_lst.append(line)
                subread_len_lst.append(i_aln.qlen)
            else: ## return ZMW ## init ZMW
                aln_cnt = len(aln_lst)
                if current_zmw in blk_set:
                    j_state = 2
                else:
                    j_state = 1
                if aln_cnt < min_subread_count:
                    j_state = 2
                else:
                    j_state = 1
                
                if j_state == 1:
                    median_subread_length = np.median(subread_len_lst) 
                    if not (lower_ccs_len_threshold < median_subread_length and median_subread_length < upper_ccs_len_threshold):
                        pass 
                    full_length_subread_count = 0
                    zmw_subread_lower_threshold = median_subread_length * lower_subread_len_threshold
                    zmw_subread_upper_threshold = median_subread_length * upper_subread_len_threshold
                    for line in aln_lst[1:-1]:
                        j_aln = BAM(line)
                        if zmw_subread_lower_threshold < j_aln.qlen and j_aln.qlen < zmw_subread_upper_threshold:
                            full_length_subread_count += 1
                            if full_length_subread_count == min_subread_count:
                                break
                    if full_length_subread_count == min_subread_count: ## at least 
                        j_state = 3    
                    else:
                        j_state = 2
                        
                if j_state == 2: ## init
                    ## init
                    i_counter += 1
                    aln_lst = [line]
                    current_zmw = i_aln.zmw
                    current_tname = i_aln.tname
                    if i_counter > 1000:
                        i_counter = 0  
                        j_counter += 1
                        if not os.path.exists("zmw/{}".format(j_counter)):
                            os.mkdir("zmw/{}".format(j_counter))
                elif j_state == 3: ## return and init
                    ## return
                    i_counter += 1
                    full_length_subread_count = 0
                    fofn.write("{}\n".format(zmw))
                    ccs = ccs_seq[tname2ccs_idx[current_tname]]
                    o = gzip.open("zmw/{}/{}.fasta.gz".format(j_counter, current_zmw), "wb")
                    o.write(">{}\n{}\n".format(ccs.name, ccs.seq).encode('utf-8')) ## return CCS read
                    for line in aln_lst[1:-1]: ## second to penultimate subreads ## return subreads
                        j_aln = BAM(line)
                        if zmw_subread_lower_threshold < j_aln.qlen and j_aln.qlen < zmw_subread_upper_threshold:
                            if j_aln.flag == 0:
                                o.write(">{}/{}\n{}\n".format(j_aln.qname, j_aln.flag, j_aln.qseq).encode('utf-8'))
                            elif j_aln.flag == 16:
                                qseq_rc = get_reverse_complement(j_aln.qseq)
                                o.write(">{}/{}\n{}\n".format(j_aln.qname, j_aln.flag, qseq_rc).encode('utf-8'))
                            full_length_subread_count += 1
                            if full_length_subread_count == min_subread_count:
                                break
                    o.close()
                    
                    ## init
                    aln_lst = [line]
                    current_zmw = i_aln.zmw
                    current_tname = i_aln.tname
                    if i_counter > 1000:
                        i_counter = 0  
                        j_counter += 1
                        if not os.path.exists("zmw/{}".format(j_counter)):
                            os.mkdir("zmw/{}".format(j_counter))

  
def main():
    options = parse_args(sys.argv)
    bam2seq(
        options.bam, 
        options.seq, 
        options.blk, 
        options.fofn,
        options.min_subread_count, 
        options.subread_lower_len_threshold, 
        options.subread_upper_len_threshold, 
    )
    sys.exit(0)


if __name__ == "__main__":
    main()


