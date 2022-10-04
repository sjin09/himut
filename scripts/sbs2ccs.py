#!/usr/bin/env python3

import re
import sys
import pysam
import cyvcf2
import natsort
import argparse
from typing import List, Tuple
from collections import defaultdict


def parse_args(args):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-i",
        "--bam",
        type=str,
        required=True,
        help="BAM file to read",
    )
    parser.add_argument(
        "--vcf",
        type=str,
        required=True,
        help="single base substitution VCF file to read",
    )
    parser.add_argument(
        "--min_bq",
        type=int,
        default=93,
        required=False,
        help="minimum base quality score",
    )
    parser.add_argument(
        "--min_mapq",
        type=int,
        default=60,
        required=False,
        help="minimum mapping quality score",
    )
    parser.add_argument(
        "-o",
        "--out",
        type=str,
        required=True,
        help="file to write",
    )
    args = args[1:]
    return parser.parse_args(args)


class BAM:
    def __init__(self, line):
        # target
        self.tname = line.reference_name
        self.tstart = line.reference_start
        self.tend = line.reference_end
        self.target_alignment_length = self.tend - self.tstart
        # query
        self.qname = line.query_name.replace("/ccs", "")
        self.qstart = line.query_alignment_start
        self.qend = line.query_alignment_end
        self.qseq = line.query_sequence
        self.qlen = len(self.qseq)
        self.mapq = line.mapping_quality
        self.bq_int_lst = line.query_qualities
        self.qv = sum(self.bq_int_lst)/self.qlen
        self.query_alignment_length = self.qend - self.qstart
        self.query_alignment_proportion = self.query_alignment_length/float(self.qlen)
        self.cs_tag = line.get_tag("cs") if line.has_tag("cs") else "."
        if line.has_tag("tp"):
            if line.get_tag("tp") == "P":
                self.is_primary = True
            else:
                self.is_primary = False
        else:
            self.is_primary = False


    def cs2lst(self):
        cslst = [cs for cs in re.split("(:[0-9]+|\*[a-z][a-z]|[=\+\-][A-Za-z]+)", self.cs_tag)]
        cslst = [cs.upper() for cs in cslst if cs != ""]
        return cslst

            
    def cs2tuple(self) -> List[Tuple[int, str, str, int, int]]:
        qpos = self.qstart
        self.cstuple_lst = []
        cs_lst = self.cs2lst()
        for cs in cs_lst:
            m = cs[1:]
            mlen = len(m)
            qstart = qpos
            if cs.startswith("="):  # match # --cs=long
                cs = ":{}".format(mlen)
                t = (1, m, m, mlen, mlen)
            elif cs.startswith(":"):  # match # --cs=short
                mlen = int(m)
                qend = qpos + mlen
                m = self.qseq[qstart:qend]
                t = (1, m, m, mlen, mlen)
            elif cs.startswith("*"):  # snp # target and query
                mlen = 1
                ref, alt = list(m)
                t = (2, ref, alt, 1, 1)
            elif cs.startswith("+"):  # insertion # query
                ref = self.qseq[qpos - 1]
                alt = ref + m
                t = (3, ref, alt, 0, mlen)
            elif cs.startswith("-"):  # deletion # target
                alt = self.qseq[qpos - 1]
                ref = alt + m
                t = (4, ref, alt, mlen, 0)
                mlen = 0
            self.cstuple_lst.append(t)
            qpos += mlen


    def cs2subindel(self):
        self.cs2tuple()
        tpos = self.tstart
        qpos = self.qstart
        self.tsbs_lst = []
        self.qsbs_lst = []
        self.qsbs_bq_lst = []
        self.tindel_lst = []
        for cstuple in self.cstuple_lst:
            state, ref, alt, ref_len, alt_len, = cstuple
            if state == 2:  # snp 
                if ref == "N": 
                    continue
                self.qsbs_bq_lst.append(self.bq_int_lst[qpos])
                self.tsbs_lst.append((tpos + 1, ref, alt))
                self.qsbs_lst.append((qpos + 1, alt, ref))
            elif state == 3 or state == 4:  # insertion
                self.tindel_lst.append((tpos, ref, alt))
            tpos += ref_len 
            qpos += alt_len


class VCF:
    def __init__(self, line):
        arr = line.strip().split()
        self.chrom = arr[0]
        self.pos = int(arr[1])
        self.id = arr[2]
        self.ref = arr[3]
        self.alt_lst = arr[4].split(",")
        self.qual = arr[5]
        self.qual = float(self.qual) if self.qual != "." else self.qual
        self.is_pass = True if arr[6] == "PASS" else False
        if len(self.alt_lst) == 1:  # bi-allelic
            self.is_biallelic = True
            self.alt = self.alt_lst[0]
        else:
            self.is_biallelic = False


def load_vcf_file(vcf_file: str):

    chrom2sbs_lst = defaultdict(list)
    if vcf_file.endswith(".vcf"):
        for line in open(vcf_file).readlines():
            if line.startswith("#"):
                continue
            v = VCF(line)
            if v.is_pass:
                chrom2sbs_lst[v.chrom].append((v.pos, v.ref, v.alt))
    elif vcf_file.endswith(".bgz"):
        for i in cyvcf2.VCF(vcf_file):
            v = VCF(str(i))
            if v.is_pass:
                chrom2sbs_lst[v.chrom].append((v.pos, v.ref, v.alt))
    return chrom2sbs_lst


def sbs2ccs(
    bam_file: str, 
    vcf_file: str,
    min_bq: int,
    min_mapq: int,
    out_file: str
): 

    o = open(out_file, "w") 
    chrom2sbs_lst = load_vcf_file(vcf_file)
    alignments = pysam.AlignmentFile(bam_file, "rb")
    chrom_lst = natsort.natsorted(list(chrom2sbs_lst.keys()))
    o.write("{}\n".format("\t".join(["sbs", "qname", "qv", "qlen", "mapq", "sbs_count", "bq93_sbs_count", "indel_count"])))
    for chrom in chrom_lst:
        for sbs in chrom2sbs_lst[chrom]:
            pos, ref, alt = sbs
            for line in alignments.fetch(chrom, pos, pos + 1):
                ccs = BAM(line)
                if not ccs.is_primary:
                    continue
                
                if ccs.mapq >= min_mapq:
                    ccs.cs2subindel()
                    tsbs_set = set(ccs.tsbs_lst)
                    if sbs in tsbs_set:
                        passed_tsbs_lst = [tsbs for tsbs, qbq in zip(ccs.tsbs_lst, ccs.qsbs_bq_lst) if qbq >= min_bq]
                        o.write("{}:{}_{}/{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                            chrom, 
                            pos, 
                            ref, 
                            alt, 
                            ccs.qname, 
                            ccs.qv, 
                            ccs.qlen, 
                            ccs.mapq, 
                            len(ccs.tsbs_lst),
                            len(passed_tsbs_lst),
                            len(ccs.tindel_lst)
                            )
                        )
                        break
    o.close() 


def main():
    options = parse_args(sys.argv)
    sbs2ccs(
        options.bam, 
        options.vcf, 
        options.min_bq,
        options.min_mapq, 
        options.out
    )
    sys.exit(0)


if __name__ == "__main__":
    main()
