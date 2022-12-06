import json
import time
import math
import pysam
import bisect
import pyfastx
import natsort
import himut.util
import himut.cslib
import himut.haplib
import himut.bamlib
import himut.reflib
import numpy as np
import multiprocessing as mp
from dataclasses import dataclass
from collections import defaultdict
from typing import Dict, List, Tuple
from himut.mutlib import (
    purine,
    sbs2tri,
    tri_lst,
    purine2pyrimidine,
)


@dataclass
class METRICS:
    num_reads: int = 0
    num_lq_reads: int = 0
    num_hq_reads: int = 0
    num_hq_phased_reads: int = 0
    num_hq_unphased_reads: int = 0
    num_total_bases: int = 0
    num_lq_bases: int = 0
    num_hq_bases: int = 0
    num_phased_hq_bases: int = 0
    num_unphased_hq_bases: int = 0
    num_aligned_hq_bases: int = 0
    num_unaligned_hq_bases: int = 0
    num_bq_filtered_hq_bases: int = 0
    num_trimmed_hq_bases: int = 0
    num_snp_filtered_hq_bases: int = 0
    num_pon_filtered_hq_bases: int = 0  
    num_mismatch_filtered_hq_bases: int = 0
    num_uncallable_hq_bases: int = 0
    num_md_filtered_hq_bases: int = 0
    num_ab_filtered_hq_bases: int = 0
    num_callable_bases: int = 0


def update_tricounts(
    pos: int, 
    seq: str, 
    tri2count: Dict[str, int],
    basecounts: np.ndarray 
) -> None:
    for idx, count in enumerate(basecounts[pos]):
        base = himut.util.idx2base[idx]
        if count == 0:
            continue
        else:
            if base in purine:
                upstream = purine2pyrimidine.get(seq[pos + 1], "N")
                downstream = purine2pyrimidine.get(seq[pos - 1], "N")
                tri = (upstream + purine2pyrimidine.get(base, "N") + downstream).upper()
            else:
                upstream = seq[pos - 1]
                downstream = seq[pos + 1]
                tri = (upstream + base + downstream).upper()
            tri2count[tri] += count


def update_allelecounts(
    read,
    allelecounts: Dict[int, np.ndarray]
) -> None:

    tpos2qbase = {}
    tpos = read.tstart
    qpos = read.qstart
    himut.cslib.cs2tuple(read)
    if read.is_primary:
        for cstuple in read.cstuple_lst:
            state, ref, alt, ref_len, alt_len, = cstuple
            if state == 1:  # match
                for i, alt_base in enumerate(alt):
                    allelecounts[tpos + i + 1][himut.util.base2idx[alt_base]] += 1
                    tpos2qbase[tpos + i + 1] = (alt_base, read.bq_int_lst[qpos + i])
            elif state == 2:  # sub
                allelecounts[tpos + 1][himut.util.base2idx[alt]] += 1
                tpos2qbase[tpos + 1] = (alt, read.bq_int_lst[qpos])
            elif state == 3:  # insertion
                allelecounts[tpos + 1][4] += 1
            elif state == 4:  # deletion
                for j in range(len(ref[1:])):
                    allelecounts[tpos + j + 1][5] += 1
                    tpos2qbase[tpos + j + 1] = ("-", 0)
            tpos += ref_len
            qpos += alt_len
    else:
        tpos2qbase = himut.cslib.cs2tpos2qbase(read)
    return tpos2qbase


def get_cumsum_tricounts(chrom2tri2count: Dict[str, Dict[str, int]]):

    tri2count = defaultdict(lambda: 0)
    for chrom in chrom2tri2count:
        for tri in tri_lst:
            tri2count[tri] += chrom2tri2count[chrom][tri]            
    return tri2count 


def get_trinucleotide_frequency(tri2count: Dict[str, int]) -> Dict[str, float]:

    tri_freq = {}
    tri_sum = sum(tri2count.values())
    for tri, count in tri2count.items():
        tri_freq[tri] = count / float(tri_sum)
    return tri_freq


def get_trinucleotide_frequency_ratio(
    ref_tricounts: Dict[str, Dict[str, int]],
    ccs_tricounts: Dict[str, Dict[str, int]],
) -> Dict[str, float]:

    trifreq_ratio = {}
    ref_trifreq = get_trinucleotide_frequency(ref_tricounts)
    ccs_trifreq = get_trinucleotide_frequency(ccs_tricounts)
    for tri, ref_freq in ref_trifreq.items():
        trifreq_ratio[tri] = ref_freq / ccs_trifreq[tri]
    return trifreq_ratio


def get_normalised_sbs96_counts(
    sbs2counts: Dict[Tuple[str, str, str], int],
    tri_ratio: Dict[Tuple[str, str, str], float],
) -> Dict[Tuple[str, str, str], float]:

    norm_sbs96_counts = {}
    for sbs, count in sbs2counts.items():
        norm_sbs96_counts[sbs] = (count * tri_ratio[sbs2tri[sbs]])
    return norm_sbs96_counts


def get_callable_tricounts(
    chrom: str,
    chrom_len: int,
    chrom_seq: str,
    bam_file: str,
    vcf_file: str,
    phased_vcf_file: str,
    common_snps: str,
    panel_of_normals: str,
    loci_lst: List[Tuple[str, int, int]],
    ploidy: str,
    min_mapq: int,
    min_trim: float,
    qlen_mean: float, 
    qlen_lower_limit: int,
    qlen_upper_limit: int,
    min_sequence_identity: float,
    min_hq_base_proportion: float,
    min_alignment_proportion: float,
    min_bq: int,
    mismatch_window: int,
    max_mismatch_count: int,
    min_base_count: int,
    md_threshold: int,
    phase: bool,
    non_human_sample: bool,
    chrom2ccs_tri2count: Dict[str, Dict[str, int]],
    chrom2base_statistics: Dict[str, List[int]],
) -> Dict[str, int]:

    seen = set()
    m = METRICS()
    tri2count = defaultdict(lambda: 0)
    for tri in tri_lst:
        tri2count[tri] = 0

    if vcf_file.endswith(".vcf"):
        sample_snp_set = himut.vcflib.load_snp(chrom, vcf_file)

    common_snp_set = set()
    if not non_human_sample:
        if common_snps.endswith(".vcf"):
            common_snp_set = himut.vcflib.load_common_snp(chrom, vcf_file)

        if panel_of_normals.endswith(".vcf"):
            pon_sbs_set, _ = himut.vcflib.load_pon(chrom, panel_of_normals)

    if phase and ploidy == "diploid": 
        (
            hpos_lst,
            hblock_lst,
            hetsnp_lst,
            hidx2hetsnp,
            hidx2hstate,
            hetsnp2bidx,
            hetsnp2hidx,
        ) = himut.vcflib.get_phased_hetsnps(phased_vcf_file, chrom, chrom_len)
        del hidx2hetsnp
        del hblock_lst
   
    alignments = pysam.AlignmentFile(bam_file, "rb")
    for loci in loci_lst:
        chunkloci_lst = himut.util.chunkloci(loci)
        for chunkloci in chunkloci_lst: 
            chunk_start, chunk_end = chunkloci[1:]
            basecounts = defaultdict(lambda: np.zeros(4)) 
            allelecounts = defaultdict(lambda: np.zeros(6))
            if vcf_file.endswith(".bgz"):
                sample_snp_set = himut.vcflib.load_bgz_snp(chunkloci, vcf_file)
            
            if not non_human_sample:
                if common_snps.endswith(".bgz"):
                    common_snp_set = himut.vcflib.load_bgz_common_snp((chrom, chunk_start - qlen_mean, chunk_end + qlen_mean), vcf_file)
                
                if panel_of_normals.endswith(".bgz"):
                    pon_sbs_set, _ = himut.vcflib.load_bgz_pon((chrom, chunk_start - qlen_mean, chunk_end + qlen_mean), panel_of_normals) 
            
            for line in alignments.fetch(*chunkloci):
                read = himut.bamlib.BAM(line)
                if read.qname in seen:
                    tpos2qbase = update_allelecounts(read, allelecounts)
                    continue  
               
                m.num_reads += 1
                seen.add(read.qname)
                m.num_total_bases += read.qlen
                ccs_somatic_tsbs_candidate_lst = [] 
                for idx, tsbs in enumerate(read.tsbs_lst):
                    if tsbs in sample_snp_set:
                        continue
                    ccs_somatic_tsbs_candidate_lst.append(tsbs)
                
                tpos2qbase = update_allelecounts(read, allelecounts)
                if not read.is_primary:
                    m.num_lq_reads += 1
                    m.num_lq_bases += read.qlen
                    continue

                if read.mapq < min_mapq:
                    m.num_lq_reads += 1
                    m.num_lq_bases += read.qlen
                    continue

                if read.qlen < qlen_lower_limit or read.qlen > qlen_upper_limit:
                    m.num_lq_reads += 1
                    m.num_lq_bases += read.qlen
                    continue

                if read.query_alignment_proportion < min_alignment_proportion:
                    m.num_lq_reads += 1
                    m.num_lq_bases += read.qlen
                    continue 

                if himut.bamlib.get_hq_base_proportion(read) < min_hq_base_proportion:
                    m.num_lq_reads += 1
                    m.num_lq_bases += read.qlen
                    continue

                if himut.util.get_blast_sequence_identity(read) < min_sequence_identity:
                    m.num_lq_reads += 1
                    m.num_lq_bases += read.qlen
                    continue
                
                common_cnt = 0
                if not non_human_sample:
                    for tsbs in ccs_somatic_tsbs_candidate_lst:
                        if tsbs in common_snp_set:
                            common_cnt += 1
                            
                    if common_cnt > 0:
                        m.num_lq_reads += 1
                        m.num_lq_bases += read.qlen
                        continue
                m.num_hq_reads += 1
                m.num_hq_bases += read.qlen

                if phase:
                    _, hap = himut.haplib.get_read_haplotype(
                        tpos2qbase,
                        hpos_lst,
                        hetsnp_lst,
                        hidx2hstate,
                        hetsnp2bidx,
                        hetsnp2hidx,
                    )
                    if hap == ".":
                        m.num_hq_unphased_reads += 1
                        continue
                    else:
                        m.num_hq_phased_reads += 1

                read_ins_count = 0
                tpos = read.tstart
                qpos = read.qstart
                trimmed_qstart = math.floor(min_trim * read.qlen)
                trimmed_qend = math.ceil((1 - min_trim) * read.qlen)
                mismatch_lst = natsort.natsorted(ccs_somatic_tsbs_candidate_lst + read.mismatch_lst)
                mismatch_set = set(mismatch_lst)
                mpos_lst = [mismatch[0] for mismatch in mismatch_lst]
                for cstuple in read.cstuple_lst:
                    state, ref, alt, ref_len, alt_len, = cstuple
                    if state == 1:  # match
                        mismatch_start, mismatch_end = himut.util.get_mismatch_range(tpos, qpos, read.qlen, mismatch_window)
                        for i, alt_base in enumerate(alt):
                            bq = read.bq_int_lst[qpos + i] 
                            if bq < min_bq:
                                m.num_bq_filtered_hq_bases += 1
                                continue

                            if qpos < trimmed_qstart: 
                                m.num_trimmed_hq_bases += 1
                                continue 
                            elif qpos > trimmed_qend:
                                m.num_trimmed_hq_bases += 1
                                continue

                            idx = bisect.bisect_left(mpos_lst, mismatch_start + i)
                            jdx = bisect.bisect_right(mpos_lst, mismatch_end + i)
                            mismatch_count = jdx - idx 
                            if mismatch_count > max_mismatch_count:                
                                m.num_mismatch_filtered_hq_bases += 1
                                continue 
                            basecounts[tpos + i][himut.util.base2idx[alt_base]] += 1
                    elif state == 2:  # sub
                        bq = read.bq_int_lst[qpos]
                        tsbs = (chrom, tpos, ref, alt)
                        if bq < min_bq:
                            tpos += ref_len
                            qpos += alt_len 
                            m.num_bq_filtered_hq_bases += 1
                            continue

                        if qpos < trimmed_qstart or qpos > trimmed_qend:
                            tpos += ref_len
                            qpos += alt_len
                            m.num_trimmed_hq_bases += 1
                            continue

                        if not non_human_sample:
                            if tsbs in pon_sbs_set:
                                tpos += ref_len
                                qpos += alt_len
                                m.num_pon_filtered_hq_bases += 1
                                continue
                            
                        if tsbs in sample_snp_set:
                            tpos += ref_len
                            qpos += alt_len
                            m.num_snp_filtered_hq_bases += 1 
                            continue

                        mismatch_start, mismatch_end = himut.util.get_mismatch_range(tpos, qpos, read.qlen, mismatch_window)
                        idx = bisect.bisect_left(mpos_lst, mismatch_start)
                        jdx = bisect.bisect_right(mpos_lst, mismatch_end)
                        if tsbs in mismatch_set:
                            mismatch_count = (jdx - idx) - 1 
                        else:
                            mismatch_count = jdx - idx 

                        if mismatch_count > max_mismatch_count:                
                            tpos += ref_len
                            qpos += alt_len
                            m.num_mismatch_filtered_hq_bases += 1
                            continue 
                        basecounts[tpos][himut.util.base2idx[alt]] += 1
                    elif state == 3:  # insertion
                        read_ins_count += alt_len
                    tpos += ref_len
                    qpos += alt_len
                m.num_unaligned_hq_bases += read_ins_count
                m.num_unaligned_hq_bases += (read.qlen - read.query_alignment_length)
                m.num_aligned_hq_bases += (read.query_alignment_length - read_ins_count) 
                
            for pos in range(*chunkloci[1:]):
                ins_count = allelecounts[pos][4] 
                del_count = allelecounts[pos][5] 
                base_count = sum(basecounts[pos])
                total_count = sum(allelecounts[pos])
                if del_count != 0 or ins_count != 0:
                    m.num_uncallable_hq_bases += base_count 

                if total_count < min_base_count:
                    m.num_ab_filtered_hq_bases += base_count
                    continue

                if total_count > md_threshold:               
                    m.num_md_filtered_hq_bases += base_count 
                    continue
                m.num_callable_bases += base_count
                update_tricounts(pos, chrom_seq, tri2count, basecounts)
        
    chrom2ccs_tri2count[chrom] = dict(tri2count)
    chrom2base_statistics[chrom] = [
        m.num_reads,
        m.num_lq_reads,
        m.num_hq_reads,
        m.num_hq_phased_reads,
        m.num_hq_unphased_reads,
        m.num_total_bases,
        m.num_lq_bases,
        m.num_hq_bases,
        m.num_phased_hq_bases,
        m.num_unphased_hq_bases,
        m.num_aligned_hq_bases,
        m.num_unaligned_hq_bases,
        m.num_bq_filtered_hq_bases,
        m.num_trimmed_hq_bases,
        m.num_snp_filtered_hq_bases,
        m.num_pon_filtered_hq_bases,
        m.num_mismatch_filtered_hq_bases,
        m.num_uncallable_hq_bases,
        m.num_ab_filtered_hq_bases,
        m.num_md_filtered_hq_bases,
        m.num_callable_bases 
    ]
    alignments.close()


def get_normcounts(
    bam_file: str,
    ref_file: str,
    sbs_file: str,
    vcf_file: str,
    phased_vcf_file: str,
    ploidy: str,
    min_mapq: int,
    min_trim: float,
    min_sequence_identity: float,
    min_hq_base_proportion: float,
    min_alignment_proportion: float,
    common_snps: str,
    panel_of_normals: str,
    min_bq: int,
    mismatch_window: int,
    max_mismatch_count: int,
    min_ref_count: int,
    min_alt_count: int,
    threads: int,
    phase: bool,
    non_human_sample: bool,
    out_file: str,
) -> None:

    cpu_start = time.time() / 60
    himut.util.check_normcounts_input_exists(
        bam_file,
        ref_file,
        sbs_file,
        vcf_file,
        phased_vcf_file,
        common_snps,
        panel_of_normals, 
        phase,
        non_human_sample,
        out_file,
    )

    print(
        "starting himut SBS96 count normalisation with {} threads".format(threads)
    )
    p = mp.Pool(threads)
    manager = mp.Manager()
    refseq = pyfastx.Fasta(ref_file)
    min_base_count = min_ref_count + min_alt_count
    _, tname2tsize = himut.bamlib.get_tname2tsize(bam_file)
    chrom_lst, sbs2count, chrom2sbs2count = himut.mutlib.load_chrom2sbs2count(sbs_file, ref_file)
    chrom2loci_lst = himut.reflib.load_seq_loci(ref_file, chrom_lst)
    qlen_mean, qlen_lower_limit, qlen_upper_limit, md_threshold = himut.bamlib.get_thresholds(
        bam_file, chrom_lst, tname2tsize
    )
    chrom2ccs_tri2count = manager.dict()
    chrom2base_statistics = manager.dict()
    tricount_arg_lst = [
        (
            chrom,
            tname2tsize[chrom],
            str(refseq[chrom]),
            bam_file,
            vcf_file,
            phased_vcf_file,
            common_snps,
            panel_of_normals,
            chrom2loci_lst[chrom],
            ploidy,
            min_mapq,
            min_trim,
            qlen_mean,
            qlen_lower_limit,
            qlen_upper_limit,
            min_sequence_identity,
            min_hq_base_proportion,
            min_alignment_proportion,
            min_bq,
            mismatch_window,
            max_mismatch_count,
            min_base_count,
            md_threshold,
            phase,
            non_human_sample,
            chrom2ccs_tri2count,
            chrom2base_statistics,
        )
        for chrom in chrom_lst
    ]
    p.starmap(get_callable_tricounts, tricount_arg_lst)
    p.close()
    p.join()
    print(
        "finished himut SBS96 count normalisation with {} threads".format(threads)
    )
    chrom2ref_tri2count = himut.reflib.get_ref_tricounts(refseq, chrom_lst, threads) ## TODO: fix
    ref_tricounts = get_cumsum_tricounts(chrom2ref_tri2count)
    ccs_tricounts = get_cumsum_tricounts(chrom2ccs_tri2count)
    trifreq_ratio = get_trinucleotide_frequency_ratio(ref_tricounts, ccs_tricounts)
    sample = himut.bamlib.get_sample(bam_file)
    himut.mutlib.dump_normcounts(sbs2count, ref_tricounts, ccs_tricounts, trifreq_ratio, out_file)
    himut.mutlib.dump_sbs96_plt(out_file, sample, out_file.replace(".tsv", ".pdf")) 
    himut.mutlib.dump_normalisation_statistics(
        chrom_lst, 
        chrom2base_statistics, 
        chrom2sbs2count, 
        chrom2ref_tri2count, 
        chrom2ccs_tri2count, 
        out_file.replace(".tsv", ".log")
    )
    print("finished returning normcounts")
    cpu_end = time.time() / 60
    duration = cpu_end - cpu_start
    print("SBS96 count normalisation took {} minutes".format(duration))
    himut.util.exit()