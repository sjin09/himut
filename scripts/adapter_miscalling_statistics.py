#!/usr/bin/env python3

import sys
import pysam
import argparse
import numpy as np

class BAM:
    def __init__(self, line):
        # query
        self.qname = line.query_name
        self.zmw = self.qname.split("/")[1]
        self.qseq = line.query_sequence
        self.qlen = len(self.qseq)


def parse_args(args):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="raw subreads BAM file to read",
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


def is_misadapter_detection(subread_len_lst, subread_median_len):

    counter = 0
    upper_len_limit = 2 * subread_median_len
    lower_len_limit = 0.5 * subread_median_len
    for subread_len in subread_len_lst[1:-1]:
        if subread_len < lower_len_limit:
            counter += 1 
        elif subread_len > upper_len_limit:
            counter += 1 
    if counter > 0:
        return True
    else:
        return False


def dump_adapter_miscalling_statistics(
    bam_file: str,
    out_file: str
): 
  
    # init
    state = 0 
    counter = 0
    o = open(out_file, "w")    
    o.write("zmw\tmedian_length\tsubread_lengths\tmiscall\n")
    alignments = pysam.AlignmentFile(bam_file, "r", check_sq=False)
    for i in alignments:
        j = BAM(i)
        if state == 0: # init
            state = 1
            current_zmw = j.zmw
            subread_len_lst = [j.qlen]
        else:
            if j.zmw != current_zmw: # return
                subread_count = len(subread_len_lst)
                if subread_count > 3:
                    full_length_len_lst = subread_len_lst[1:-1]
                    subread_median_len = np.median(full_length_len_lst)
                    full_length_len_lst = [str(l) for l in full_length_len_lst]
                    if is_misadapter_detection(subread_len_lst, subread_median_len):
                        o.write("{}\t{}\t{}\t{}\n".format(current_zmw, str(subread_median_len), ",".join(full_length_len_lst), "TRUE"))
                    else: 
                        o.write("{}\t{}\t{}\t{}\n".format(current_zmw, str(subread_median_len), ",".join(full_length_len_lst), "FALSE"))
                else: 
                        o.write("{}\t{}\t{}\t{}\n".format(current_zmw, ".", ",".join(subread_len_lst), "."))
                # init
                state = 1
                current_zmw = j.zmw
                subread_len_lst = [j.qlen]
                counter += 1
                if counter > 99:
                    break
            else: # init
                subread_len_lst.append(j.qlen)
    o.close()
    
def main():
    options = parse_args(sys.argv)
    dump_adapter_miscalling_statistics(
        options.input,
        options.output
    )
    sys.exit(0)


if __name__ == "__main__":
    main()
