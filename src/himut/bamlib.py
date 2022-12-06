import json
import math
import pysam
import random
import natsort
import numpy as np
import himut.util
import himut.cslib
import himut.vcflib
from collections import defaultdict
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
        self.strand = "+" if line.is_forward else "-"
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
    qlen_lower_limit = 0 if math.ceil(qlen_mean - 2 * qlen_std) < 0 else math.ceil(qlen_mean - 2 * qlen_std)
    qlen_upper_limit = math.ceil(qlen_mean + 2 * qlen_std)
    coverage = genome_read_sum / float(genome_sample_sum)
    md_threshold = get_md_threshold(coverage)
    return qlen_std, qlen_mean, qlen_lower_limit, qlen_upper_limit, md_threshold


def get_hq_base_proportion(read):
    hq_base_proportion = read.bq_int_lst.count(93)/float(read.qlen)
    return hq_base_proportion


def get_basecounts(
    alignments,
    chrom: str,
    start: int,
):
    
    basecounts = [0] * 4
    base2bq_lst = defaultdict(lambda: {0: [], 1:[], 2:[], 3:[]})
    for pileupcolumn in alignments.pileup(chrom, start, start+1, stepper="samtools", flag_filter=256, min_base_quality=0, min_mapping_quality=0):
        if pileupcolumn.pos == start - 1:
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip:
                    base = pileupread.alignment.query_sequence[pileupread.query_position]
                    base_bq = pileupread.alignment.query_qualities[pileupread.query_position]
                    basecounts[himut.util.base2idx[base]] += 1
                    base2bq_lst[himut.util.base2idx[base]].append(base_bq)
            break
    return basecounts, base2bq_lst


def get_sbs_allelecounts(
    ref: str,
    alt: str,
    allelecounts: Dict[int, Dict[int, int]],
    allele2bq_lst: Dict[int, Dict[str, List[int]]],
):
   
    read_depth = sum(allelecounts)
    indel_count = allelecounts[4] + allelecounts[5]
    ref_count = allelecounts[himut.util.base2idx[ref]] 
    alt_count = allelecounts[himut.util.base2idx[alt]]
    vaf = alt_count/float(read_depth) 
    bq = sum(allele2bq_lst[himut.util.base2idx[alt]])/float(alt_count)
    return bq, vaf, ref_count, alt_count, indel_count, read_depth


def get_dbs_allelecounts(
    ref: str,
    alt: str,
    allelecounts: Dict[int, Dict[int, int]],
    allele2bq_lst: Dict[int, Dict[str, List[int]]],
):
   
    read_depth = sum(allelecounts)
    indel_count = allelecounts[4] + allelecounts[5]
    ref_count = allelecounts[himut.util.base2idx[ref]] 
    alt_count = allelecounts[himut.util.base2idx[alt]]
    vaf = alt_count/float(read_depth) 
    bq = sum(allele2bq_lst[himut.util.base2idx[alt]])/float(alt_count)
    return bq, vaf, ref_count, alt_count, indel_count, read_depth


def get_dbs_allelecounts(
    pos: int,
    ref: str,
    alt: str,
    read2tpos2qbase: Dict[str, Dict[int,  str]],
    tpos2allelecounts: Dict[int, np.ndarray],
    tpos2qbase2read_lst: Dict[int, Dict[int, List[str]]],
) -> Tuple[int, int, str, float, int, int, int]:
 
    ref_count = 0
    alt_count = 0
    indel_count = sum(tpos2allelecounts[pos][4], tpos2allelecounts[pos+1][4], tpos2allelecounts[pos][5], tpos2allelecounts[pos+1][5])
    ref_read_lst = list(set(tpos2qbase2read_lst[pos][himut.util.base2idx[ref[0]]]).intersection(set(tpos2qbase2read_lst[pos + 1][himut.util.base2idx[ref[1]]])))
    alt_read_lst = list(set(tpos2qbase2read_lst[pos][himut.util.base2idx[alt[0]]]).intersection(set(tpos2qbase2read_lst[pos + 1][himut.util.base2idx[alt[1]]])))

    for ref_read in ref_read_lst:
        ref_read_allele = "{}{}".format(read2tpos2qbase[ref_read][pos][0], read2tpos2qbase[ref_read][pos + 1][0])
        if ref == ref_read_allele:
            ref_count += 1 
               
    for alt_read in alt_read_lst:
        alt_read_allele = "{}{}".format(read2tpos2qbase[alt_read][pos][0], read2tpos2qbase[alt_read][pos + 1][0])
        if alt == alt_read_allele:
            alt_count += 1 

    read_depth = len(ref_read_lst) + len(alt_read_lst)
    vaf = alt_count / float(read_depth)
    return vaf, ref_count, alt_count, indel_count, read_depth, ref_read_lst, alt_read_lst

