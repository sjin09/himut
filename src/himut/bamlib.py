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
        hq_base_proportion = self.bq_int_lst.count(93) / float(self.qlen)
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
        blast_sequence_identity = (target_alignment_len - mismatch_count) / float(
            target_alignment_len
        )
        return blast_sequence_identity

    def get_query_alignment_proportion(self):
        query_alignment_proportion = (self.qend - self.qstart) / float(self.qlen)
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
    alt_bq = sum(allele2bq_lst[himut.util.base2idx[alt]]) / float(alt_count)
    alt_vaf = alt_count / float(read_depth)
    return ref_count, alt_bq, alt_vaf, alt_count, indel_count, read_depth


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
    mismatch_tpos_lst = [mismatch[0] for mismatch in ccs.mismatch_lst]
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

