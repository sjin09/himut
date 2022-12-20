import math
import pysam
import random
import bisect
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
        if line.is_secondary:
            self.is_primary = False
        else: 
            self.is_primary = True
            self.tname = line.reference_name
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
            himut.cslib.cs2tuple(self, line.get_tag("cs"))

    def cs2mut(self):
        himut.cslib.cs2mut(self)

    def cs2subindel(self):
        himut.cslib.cs2subindel(self)
    
    def cs2tpos2qbase(self):
        himut.cslib.cs2tpos2qbase(self) 
        
    def get_hq_base_proportion(self):
        hq_base_proportion = self.bq_int_lst.count(93)/float(self.qlen)
        return hq_base_proportion

    def get_blast_sequence_identity(self):
        
        mismatch_count = 0
        target_alignment_len = self.tstart - self.tend 
        for cstuple in self.cstuple_lst:
            mstate, _, _, ref_len, alt_len = cstuple
            if mstate == 1:  # match
                continue
            elif mstate == 2:  # mismatch: snp
                mismatch_count += alt_len
            elif mstate == 3:  # mismatch: insertion
                mismatch_count += alt_len
            elif mstate == 4:  # mismatch: deletion
                mismatch_count += ref_len
        blast_sequence_identity = (target_alignment_len  - mismatch_count)/float(target_alignment_len)
        return blast_sequence_identity

    def get_query_alignment_proportion(self):
        query_alignment_proportion = (self.qend - self.qstart)/float(self.qlen) 
        return query_alignment_proportion


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
    md_threshold = math.ceil(coverage + (4*math.sqrt(coverage)))
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
    return qlen_lower_limit, qlen_upper_limit, md_threshold


def get_basecounts(
    chrom: str,
    start: int,
    alignments
):
    
    basecounts = [0] * 4
    base2bq_lst = {0: [], 1:[], 2:[], 3:[]}
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


def get_allelecounts(
    chrom: str,
    pos: int,
    alignments
):

    tpos2allelecounts = defaultdict(lambda: np.zeros(6)) 
    tpos2allele2bq_lst = defaultdict(lambda: {0: [], 1:[], 2:[], 3:[], 4:[], 5:[]})
    for i in alignments.fetch(chrom, pos-1, pos+1):
        ccs = himut.bamlib.BAM(i)
        tpos = ccs.tstart
        qpos = ccs.qstart
        for (state, ref, alt, ref_len, alt_len) in ccs.cstuple_lst:
            if state == 1:  # match
                for i, alt_base in enumerate(alt):
                    tpos2allelecounts[tpos + i + 1][himut.util.base2idx[alt_base]] += 1
                    tpos2allele2bq_lst[tpos + i + 1][himut.util.base2idx[alt_base]].append(ccs.bq_int_lst[qpos + i])
            elif state == 2:  # sub
                tpos2allelecounts[tpos + 1][himut.util.base2idx[alt]] += 1
                tpos2allele2bq_lst[tpos + 1][himut.util.base2idx[alt]].append(ccs.bq_int_lst[qpos])
            elif state == 3:  # insertion
                tpos2allelecounts[tpos + 1][4] += 1
                pass
            elif state == 4:  # deletion
                for j in range(len(ref[1:])):
                    tpos2allelecounts[tpos + j + 1][5] += 1
            tpos += ref_len
            qpos += alt_len
    allelecounts = tpos2allelecounts[pos] 
    allele2bq_lst = tpos2allele2bq_lst[pos]
    return allelecounts, allele2bq_lst


def get_sbs_allelecounts(
    ref: str,
    alt: str,
    allelecounts: Dict[int, Dict[int, int]],
    allele2bq_lst: Dict[int, Dict[str, List[int]]],
):

    ins_count = allelecounts[4] 
    del_count = allelecounts[5]
    indel_count = ins_count + del_count 
    read_depth = sum(allelecounts) - ins_count 
    ref_count = allelecounts[himut.util.base2idx[ref]]
    alt_count = allelecounts[himut.util.base2idx[alt]]
    alt_bq = sum(allele2bq_lst[himut.util.base2idx[alt]])/float(alt_count)
    alt_vaf = alt_count/float(read_depth) 
    return ref_count, alt_bq, alt_vaf, alt_count, indel_count, read_depth


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


def is_trimmed(
    ccs,
    qpos: int,
    min_trim: float,
) -> bool:

    trimmed_qstart = math.floor(min_trim*ccs.qlen)
    trimmed_qend = math.ceil((1 - min_trim)*ccs.qlen)
    if qpos < trimmed_qstart:
        return True
    elif qpos > trimmed_qend:
        return True
    else:
        return False
  

def get_mismatch_range(
    tpos: int, 
    qpos: int,
    qlen: int, 
    window: int
):
    qstart, qend = [qpos - window, qpos + window]
    if qstart < 0:
        urange = window + qstart
        drange = window + abs(qstart)
    elif qend > qlen:
        urange = window + abs(qend - qlen)
        drange = qlen - qpos
    else:
        urange = window
        drange = window
    tstart = tpos - urange
    tend = tpos + drange
    return tstart, tend 
   
    
def is_mismatch_conflict(
    ccs,
    tpos: int,
    qpos: int,
    mismatch_window: int,
    max_mismatch_count: int
) -> bool:

    mpos_lst = [mismatch[0] for mismatch in ccs.mismatch_lst]
    mismatch_start, mismatch_end = get_mismatch_range(tpos, qpos, ccs.qlen, mismatch_window)
    mismatch_count = bisect.bisect_right(mpos_lst, mismatch_end) - bisect.bisect_left(mpos_lst, mismatch_start) - 1
    if mismatch_count > max_mismatch_count:
        return True
    else:
        return False