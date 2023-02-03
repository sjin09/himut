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

    def get_qv(self):
        self.qv = np.mean(self.bq_int_lst)
        return self.qv

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
 
        match_count = 0
        mismatch_count = 0
        for cstuple in self.cstuple_lst:
            mstate, _, _, ref_len, alt_len = cstuple
            if mstate == 1:  # match
                match_count += ref_len
            elif mstate == 2:  # mismatch: snp
                mismatch_count += alt_len
            elif mstate == 3:  # mismatch: insertion
                mismatch_count += alt_len
            elif mstate == 4:  # mismatch: deletion
                mismatch_count += ref_len
        alignment_len = match_count + mismatch_count
        blast_sequence_identity = (match_count)/float(alignment_len)
        return blast_sequence_identity

    def get_query_alignment_proportion(self):
        qaln_proportion = (self.qend - self.qstart)/float(self.qlen)
        return qaln_proportion

    def get_tsbs_candidates(self, som_seen, min_trim, mismatch_window, max_mismatch_count):

        self.cs2subindel()
        ccs_somatic_tsbs_candidate_lst = []
        trimmed_qstart, trimmed_qend = himut.bamlib.get_trimmed_range(self.qlen, min_trim)
        for tsbs, qsbs in zip(self.tsbs_lst, self.qsbs_lst):
            tpos = tsbs[0]
            qpos = qsbs[0]
            if tpos in som_seen:
                continue
            if himut.bamlib.is_trimmed(qpos, trimmed_qstart, trimmed_qend):
                continue
            if himut.bamlib.is_mismatch_conflict(
                self, tpos, qpos, mismatch_window, max_mismatch_count
            ):
                continue 
            ccs_somatic_tsbs_candidate_lst.append(tsbs)
        return ccs_somatic_tsbs_candidate_lst
        

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
        print(
            "Please use samtools reheader to insert approprirate header to your BAM file"
        )
        himut.util.exit()
    return tname_lst, tname2tsize


def get_md_threshold(coverage: int) -> int:
    md_threshold = math.ceil(coverage + (4 * math.sqrt(coverage)))
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
    qlen_lower_limit = (
        0
        if math.ceil(qlen_mean - 2 * qlen_std) < 0
        else math.ceil(qlen_mean - 2 * qlen_std)
    )
    qlen_upper_limit = math.ceil(qlen_mean + 2 * qlen_std)
    coverage = genome_read_sum / float(genome_sample_sum)
    md_threshold = get_md_threshold(coverage)
    return qlen_lower_limit, qlen_upper_limit, md_threshold


def get_ref_counts(
    ref: str,
    read_depth: int,
    allelecounts: Dict[int, Dict[int, int]],
    allele2bq_lst: Dict[int, Dict[str, List[int]]],
):

    ref_count = allelecounts[himut.util.base2idx[ref]]
    ref_vaf = ref_count/float(read_depth)
    if ref_count != 0:
        ref_bq = sum(allele2bq_lst[himut.util.base2idx[ref]])/float(ref_count)
    else:
        ref_bq = 0.0
    return ref_bq, ref_vaf, ref_count


def get_alt_counts(
    alt: str,
    read_depth: int,
    allelecounts: Dict[int, Dict[int, int]],
    allele2bq_lst: Dict[int, Dict[str, List[int]]],
):

    alt_count = allelecounts[himut.util.base2idx[alt]]
    alt_vaf = alt_count/float(read_depth)
    if alt_count != 0:
        alt_bq = sum(allele2bq_lst[himut.util.base2idx[alt]])/float(alt_count)
    else:
        alt_bq = 0.0
    return alt_bq, alt_vaf, alt_count


def get_read_depth(
    allelecounts: Dict[int, Dict[int, int]],
):
    ins_count = allelecounts[4]
    del_count = allelecounts[5]
    read_depth = sum(allelecounts) - ins_count
    return del_count, ins_count, read_depth


def get_trimmed_range(
    qlen: int,
    min_trim: float,
):
    trimmed_qstart = math.floor(min_trim * qlen)
    trimmed_qend = math.ceil((1 - min_trim) * qlen)
    return trimmed_qstart, trimmed_qend
    

def is_trimmed(
    qpos: int,
    trimmed_qstart: float,
    trimmed_qend: float
) -> bool:

    if qpos < trimmed_qstart:
        return True
    elif qpos > trimmed_qend:
        return True
    else:
        return False


def get_mismatch_range(tpos: int, qpos: int, qlen: int, window: int):
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


def get_mismatch_positions(ccs):
    mismatch_tpos_lst = [mismatch[0] for mismatch in ccs.mismatch_lst] # 1-coordinate
    return mismatch_tpos_lst


def is_mismatch_conflict(
    ccs, tpos: int, qpos: int, mismatch_window: int, max_mismatch_count: int
) -> bool:

    mismatch_tpos_lst = himut.bamlib.get_mismatch_positions(ccs) 
    mismatch_start, mismatch_end = get_mismatch_range(
        tpos, qpos, ccs.qlen, mismatch_window
    )
    mismatch_count = (
        bisect.bisect_right(mismatch_tpos_lst, mismatch_end)
        - bisect.bisect_left(mismatch_tpos_lst, mismatch_start)
        - 1
    )
    if mismatch_count > max_mismatch_count:
        return True
    else:
        return False

