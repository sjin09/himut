import json
import math
import pysam
import random
import natsort
import numpy as np
import himut.util
import himut.cslib
import himut.vcflib
from typing import List, Dict, Tuple


class BAM:
    def __init__(self, line):
        # target
        self.tname = line.reference_name
        self.tstart = line.reference_start
        self.tend = line.reference_end
        self.target_alignment_length = self.tend - self.tstart
        # query
        self.qname = line.query_name
        self.qstart = line.query_alignment_start
        self.qend = line.query_alignment_end
        self.qseq = line.query_sequence
        self.qlen = len(self.qseq)
        self.mapq = line.mapping_quality
        self.bq_int_lst = line.query_qualities
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
        himut.cslib.cs2tuple(self)
        himut.cslib.cs2mut(self)

def get_sample(bam_file: str) -> str:
   
    state = 0
    alignments = pysam.AlignmentFile(bam_file, "rb")
    bam_header_lines = str(alignments.header).strip().split("\n")
    for line in bam_header_lines:
        if line.startswith("@RG"):
            for i in line.split():
                if i.startswith("SM"):
                    sample = i.split(":")[1]
                    return sample
    if state == 0:
        print("SM field is missing")
        print("Please provide a BAM file with @RG group")
        print(
            "samtools reheader in.header.sam in.bam > out.bam command can be used to insert a new header"
        )
        himut.util.exit() 


def get_tname2tsize(bam_file: str) -> Tuple[List[str], Dict[str, int]]:

    tname2tsize = {}
    alignments = pysam.AlignmentFile(bam_file, "rb")
    bam_header_lst = str(alignments.header).strip().split("\n")
    for h in bam_header_lst:
        if h.startswith("@SQ"):
            arr = h.split("\t")
            tname = arr[1].replace("SN:", "")
            tsize = int(arr[2].replace("LN:", ""))
            tname2tsize[tname] = tsize
    alignments.close()

    tname_lst = natsort.natsorted(list(tname2tsize.keys()))
    if len(tname_lst) == 0:
        print("@SQ header is missing from BAM file")
        print("Please use samtools reheader to insert approprirate header to your BAM file")
        himut.util.exit()
    return tname_lst, tname2tsize


def get_md_threshold(coverage: int) -> int:
    md_threshold = math.ceil(coverage + 4 * math.sqrt(coverage))
    return md_threshold


def get_thresholds(
    bam_file: str, 
    chrom_lst: List[str], 
    chrom2len: Dict[str, int],
) -> Tuple[int, int, int, int]:


    if len(chrom_lst) == 0:
        print("target is missing")
        print("Please check .vcf file or .target file")
        himut.util.exit() 

    qlen_lst = []
    random.seed(10)
    sample_count = 100
    sample_range = 100000
    genome_read_sum = 0
    genome_sample_sum = sample_count * sample_range * len(chrom_lst)
    alignments = pysam.AlignmentFile(bam_file, "rb")
    for chrom in chrom_lst:
        chrom_len = chrom2len[chrom]
        random_start_lst = random.sample(range(chrom_len), sample_count)
        for start in random_start_lst:
            end = start + 100000
            for read in alignments.fetch(chrom, start, end):
                mapq = int(read.mapping_quality)
                alignment_type = read.get_tag("tp") if read.has_tag("tp") else "."
                if mapq > 0 and alignment_type == "P":
                    qlen = len(read.query_sequence)
                    genome_read_sum += qlen
                    qlen_lst.append(qlen)
    alignments.close()
    qlen_std = np.std(qlen_lst)
    qlen_mean = math.ceil(np.mean(qlen_lst))
    qlen_lower_limit = math.ceil(qlen_mean - 2 * qlen_std)
    qlen_upper_limit = math.ceil(qlen_mean + 2 * qlen_std)
    if qlen_lower_limit < 0:
        qlen_lower_limit = 0
    coverage = genome_read_sum / float(genome_sample_sum)
    md_threshold = get_md_threshold(coverage)
    return qlen_mean, qlen_lower_limit, qlen_upper_limit, md_threshold


def get_hq_base_proportion(read):
    hq_base_proportion = read.bq_int_lst.count(93)/float(read.qlen)
    return hq_base_proportion


