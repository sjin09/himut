import time
import pysam
import bisect
import pyfastx
import himut.util
import himut.cslib
import himut.gtlib
import himut.bamlib
import himut.caller
import himut.haplib
import himut.somlib
import himut.reflib
import numpy as np
import multiprocessing as mp
from collections import defaultdict
from typing import Set, Dict, List, Tuple
from himut.mutlib import (
    purine,
    tri_lst,
    purine2pyrimidine,
)


def init_allelecounts():
    rpos2basecounts = defaultdict(lambda: np.zeros(4))
    rpos2allelecounts = defaultdict(lambda: np.zeros(6))
    rpos2allele2bq_lst = defaultdict(lambda: {0: [], 1: [], 2: [], 3: [], 4: [], 5: []})
    rpos2allele2ccs_lst = defaultdict(
        lambda: {0: [], 1: [], 2: [], 3: [], 4: [], 5: []}
    )
    return rpos2basecounts, rpos2allele2bq_lst, rpos2allele2ccs_lst, rpos2allelecounts


def is_low_bq(
    bq: int, 
    min_bq: int, 
) -> bool:

    if bq < min_bq:
        return True
    else:
        return False


def update_basecounts(
    ccs,
    min_bq: int,
    min_trim: float,
    mismatch_window: int,
    max_mismatch_count: int,
    non_human_sample: bool,
    pon_sbs_set: Set[Tuple[int, str, str]],
    common_snp_set: Set[Tuple[int, str, str]],
    tpos2basecounts: Dict[int, int]
):

    i = 0
    ccs.cs2mut()
    tpos = ccs.tstart
    qpos = ccs.qstart
    mismatch_tpos_lst = himut.bamlib.get_mismatch_positions(ccs) 
    trimmed_qstart, trimmed_qend = himut.bamlib.get_trimmed_range(ccs.qlen, min_trim)
    mismatch_tstart, mismatch_tend = himut.bamlib.get_mismatch_range(tpos, qpos, ccs.qlen, mismatch_window)
    for cstuple in ccs.cstuple_lst:
        state, ref, alt, ref_len, alt_len = cstuple
        if state == 1:  # match
            for j, alt_base in enumerate(alt):
                bq = ccs.bq_int_lst[qpos + j] 
                mismatch_count = (
                    bisect.bisect_right(mismatch_tpos_lst, mismatch_tend + i)
                    - bisect.bisect_left(mismatch_tpos_lst, mismatch_tstart + i)
                )
                if is_low_bq(bq, min_bq):
                    continue
                if himut.bamlib.is_trimmed(qpos, trimmed_qstart, trimmed_qend): 
                    continue
                if mismatch_count > max_mismatch_count:
                    continue
                tpos2basecounts[tpos + j][himut.util.base2idx[alt_base]] += 1
                i += 1
        elif state == 2:  
            sbs_state = 1
            tsbs = (tpos, ref, alt)
            bq = ccs.bq_int_lst[qpos]
            mismatch_count = (
                bisect.bisect_right(mismatch_tpos_lst, mismatch_tend + i)
                - bisect.bisect_left(mismatch_tpos_lst, mismatch_tstart + i)
                - 1
            )
            if is_low_bq(bq, min_bq):
                sbs_state = 0
            if mismatch_count > max_mismatch_count:                
                sbs_state = 0
            if (tsbs in pon_sbs_set or tsbs in common_snp_set) and not non_human_sample:
                sbs_state = 0
            if himut.bamlib.is_trimmed(qpos, trimmed_qstart, trimmed_qend):
                sbs_state = 0
            if sbs_state:
                tpos2basecounts[tpos][himut.util.base2idx[alt]] += 1
            i += 1
        elif state == 4: # deletion
            i += ref_len
        tpos += ref_len
        qpos += alt_len 
    

def is_phased(
   ccs, 
   hbit_lst: List[str], 
   hpos_lst: List[int], 
   hetsnp_lst: List[Tuple[str, int, int]]
) -> bool:

    hap = himut.haplib.get_ccs_hap(ccs, hbit_lst, hpos_lst, hetsnp_lst) 
    if hap == "0":
        return True
    elif hap == "1":
        return True
    else:
        return False


def get_read_depth(allelecounts: np.ndarray):
    read_depth = sum(allelecounts) - allelecounts[4]
    return read_depth


def get_indel_count(allelecounts: np.ndarray):
    indel_count = allelecounts[4] + allelecounts[5]
    return indel_count


def update_tricounts(
    pos: int,
    base: str,
    hq_base_count: int,
    chrom_seq: str,
    tri2count: Dict[str, int],
):
    if base in purine:
        upstream = purine2pyrimidine.get(chrom_seq[pos + 1], "N")
        downstream = purine2pyrimidine.get(chrom_seq[pos - 1], "N")
        tri = "{}{}{}".format(upstream, purine2pyrimidine.get(base, "N"), downstream)
    else:
        upstream = chrom_seq[pos - 1]
        downstream = chrom_seq[pos + 1]
        tri = "{}{}{}".format(upstream, base, downstream)
    tri2count[tri] += hq_base_count


def get_callable_tricounts(
    chrom: str,
    chrom_seq: str,
    bam_file: str,
    common_snps: str,
    panel_of_normals: str,
    chunkloci_lst: List[Tuple[str, int, int]],
    phase_set2hbit_lst: Dict[str, List[str]],
    phase_set2hpos_lst: Dict[int, List[int]],
    phase_set2hetsnp_lst: Dict[int, List[Tuple[int, str, str]]], 
    min_mapq: int,
    min_trim: float,
    qlen_lower_limit: int,
    qlen_upper_limit: int,
    min_sequence_identity: float,
    min_hq_base_proportion: float,
    min_alignment_proportion: float,
    min_gq: int,
    min_bq: int,
    mismatch_window: int,
    max_mismatch_count: int,
    min_ref_count,
    min_alt_count,
    md_threshold: int,
    somatic_snv_prior: float,
    germline_snv_prior: float,
    germline_indel_prior: float,
    phase: bool,
    non_human_sample: bool,
    chrom2ccs_tri2count: Dict[str, Dict[str, int]],
) -> Dict[str, int]:

    pon_sbs_set = set()
    common_snp_set = set()
    tri2count = defaultdict(lambda: 0)
    himut.gtlib.init(germline_snv_prior)
    himut.somlib.init(germline_snv_prior)
    for tri in tri_lst:
        tri2count[tri] = 0
    min_read_depth = min_ref_count + min_alt_count
    alignments = pysam.AlignmentFile(bam_file, "rb")
    if not non_human_sample:
        if common_snps.endswith(".vcf"):
            common_snp_set = himut.vcflib.load_common_snp(chrom, common_snps)
        if panel_of_normals.endswith(".vcf"):
            pon_sbs_set = himut.vcflib.load_pon(chrom, panel_of_normals)

    for (_chrom, chunk_start, chunk_end) in chunkloci_lst:  ## TODO
        if not non_human_sample:
            if common_snps.endswith(".bgz"):
                common_snp_set = himut.vcflib.load_bgz_common_snp(
                    (
                        chrom,
                        chunk_start - qlen_upper_limit,
                        chunk_end + qlen_upper_limit,
                    ),
                    common_snps,
                )
            if panel_of_normals.endswith(".bgz"):
                pon_sbs_set = himut.vcflib.load_bgz_pon(
                    (
                        chrom,
                        chunk_start - qlen_upper_limit,
                        chunk_end + qlen_upper_limit,
                    ),
                    panel_of_normals,
                )

        rpos2basecounts, rpos2allele2bq_lst, rpos2allele2ccs_lst, rpos2allelecounts = init_allelecounts()
        for i in alignments.fetch(chrom, chunk_start, chunk_end):
            ccs = himut.bamlib.BAM(i)
            phase_set = str(chunk_start)
            hbit_lst = phase_set2hbit_lst[phase_set] 
            hpos_lst = phase_set2hpos_lst[phase_set] 
            hetsnp_lst = phase_set2hetsnp_lst[phase_set] 
            if not ccs.is_primary:
                continue
            himut.cslib.update_allelecounts(
                ccs, rpos2allelecounts, rpos2allele2bq_lst, rpos2allele2ccs_lst
            )
            if himut.caller.is_low_mapq(ccs.mapq, min_mapq):
                continue
            if ccs.qlen < qlen_lower_limit or ccs.qlen > qlen_upper_limit:
                continue
            if ccs.get_hq_base_proportion() < min_hq_base_proportion:
                continue
            if ccs.get_blast_sequence_identity() < min_sequence_identity:
                continue
            if ccs.get_query_alignment_proportion() < min_alignment_proportion:
                continue

            if phase and is_phased(ccs, hbit_lst, hpos_lst, hetsnp_lst):
                update_basecounts(ccs, min_bq, min_trim, mismatch_window, max_mismatch_count, non_human_sample, pon_sbs_set, common_snp_set, rpos2basecounts)
            else:
                update_basecounts(ccs, min_bq, min_trim, mismatch_window, max_mismatch_count, non_human_sample, pon_sbs_set, common_snp_set, rpos2basecounts)

        for tpos in range(chunk_start, chunk_end):  # 1-based coordinate
            rpos = tpos - 1  # 0-based coordinate
            ref = chrom_seq[rpos]
            if ref == "N":
                continue
            basecounts = rpos2basecounts[rpos]
            allelecounts = rpos2allelecounts[rpos]
            allele2bq_lst = rpos2allele2bq_lst[rpos]
            read_depth = get_read_depth(allelecounts)
            indel_count = get_indel_count(allelecounts)
            _germ_gt, germ_gq, germ_gt_state, gt2gt_state = himut.gtlib.get_germ_gt(
                ref, allele2bq_lst
            )
            if germ_gt_state == "het":
                continue
            elif germ_gt_state == "hetalt":
                continue
            elif germ_gt_state == "homalt":
                continue
            if indel_count != 0:
                continue
            if read_depth < min_read_depth:
                continue
            elif read_depth > md_threshold:
                continue
            
            ref_sum = allelecounts[himut.util.base2idx[ref]]
            ref_count = basecounts[himut.util.base2idx[ref]]
            if ref_sum == read_depth: 
                if himut.caller.is_low_gq(germ_gq, min_gq):
                    continue
                update_tricounts(rpos, ref, ref_count, chrom_seq, tri2count)
                continue
            
            alt_lst = himut.util.base_set.difference(ref)
            alt_count_lst = [basecounts[himut.util.base2idx[alt_base]] for alt_base in alt_lst]
            for alt, alt_count in zip(alt_lst, alt_count_lst): 
                if alt_count == 0: 
                    continue
                germ_gq = himut.gtlib.get_germ_gq(alt, gt2gt_state, allele2bq_lst)
                if himut.caller.is_low_gq(germ_gq, min_gq):
                    continue
                update_tricounts(rpos, alt, alt_count, chrom_seq, tri2count)
            update_tricounts(rpos, ref, ref_count, chrom_seq, tri2count)
    chrom2ccs_tri2count[chrom] = dict(tri2count)
    alignments.close()


def get_normcounts(
    bam_file: str,
    ref_file: str,
    sbs_file: str,
    vcf_file: str,
    phased_vcf_file: str,
    common_snps: str,
    panel_of_normals: str,
    region: str,
    region_list: str,
    min_mapq: int,
    min_sequence_identity: float,
    min_hq_base_proportion: float,
    min_alignment_proportion: float,
    min_gq: int,
    min_bq: int,
    min_trim: float,
    mismatch_window: int,
    max_mismatch_count: int,
    min_ref_count: int,
    min_alt_count: int,
    somatic_snv_prior: float,
    germline_snv_prior: float,
    germline_indel_prior: float,
    threads: int,
    phase: bool,
    non_human_sample: bool,
    reference_sample: bool,
    out_file: str,
) -> None:

    cpu_start = time.time() / 60
    _, tname2tsize = himut.bamlib.get_tname2tsize(bam_file)
    chrom_lst, _ = himut.util.load_loci(region, region_list, tname2tsize)
    himut.util.check_normcounts_input_exists(
        bam_file,
        ref_file,
        sbs_file,
        vcf_file,
        phased_vcf_file,
        common_snps,
        panel_of_normals,
        chrom_lst,
        tname2tsize,
        phase,
        reference_sample,
        non_human_sample,
        out_file,
    )
    
    print("starting himut SBS96 count normalisation with {} threads".format(threads))
    if phase:
        (
            chrom2ps2hbit_lst,
            chrom2ps2hpos_lst,
            chrom2ps2hetsnp_lst,
            chrom2chunkloci_lst,
        ) = himut.vcflib.load_phased_hetsnps(phased_vcf_file, chrom_lst, tname2tsize)
    else:
        chrom2ps2hbit_lst = defaultdict(dict)
        chrom2ps2hpos_lst = defaultdict(dict)
        chrom2ps2hetsnp_lst = defaultdict(dict)
        chrom2chunkloci_lst = himut.reflib.load_seq_loci(ref_file, chrom_lst)
    qlen_lower_limit, qlen_upper_limit, md_threshold = himut.bamlib.get_thresholds(
        bam_file, chrom_lst, tname2tsize
    )
    if non_human_sample:
        germline_snv_prior, germline_indel_prior = himut.vcflib.get_germline_priors(
            chrom_lst, ref_file, vcf_file, reference_sample
        )

    p = mp.Pool(threads)
    manager = mp.Manager()
    refseq = pyfastx.Fasta(ref_file)
    chrom2ccs_tri2count = manager.dict()
    tricount_arg_lst = [
        (
            chrom,
            str(refseq[chrom]),
            bam_file,
            common_snps,
            panel_of_normals,
            chrom2chunkloci_lst[chrom],
            chrom2ps2hbit_lst[chrom],
            chrom2ps2hpos_lst[chrom],
            chrom2ps2hetsnp_lst[chrom],
            min_mapq,
            min_trim,
            qlen_lower_limit,
            qlen_upper_limit,
            min_sequence_identity,
            min_hq_base_proportion,
            min_alignment_proportion,
            min_gq,
            min_bq,
            mismatch_window,
            max_mismatch_count,
            min_ref_count,
            min_alt_count,
            md_threshold,
            somatic_snv_prior,
            germline_snv_prior,
            germline_indel_prior,
            phase,
            non_human_sample,
            chrom2ccs_tri2count,
        )
        for chrom in chrom_lst
    ]
    p.starmap(get_callable_tricounts, tricount_arg_lst) 
    p.close()
    p.join()
    print("finished himut SBS96 count normalisation with {} threads".format(threads))
    cmdline = himut.mutlib.get_normcounts_cmdline(
        bam_file,
        ref_file,
        sbs_file,
        vcf_file,
        phased_vcf_file,
        min_mapq,
        min_sequence_identity,
        min_hq_base_proportion,
        min_alignment_proportion,
        min_gq,
        min_bq,
        min_trim,
        mismatch_window,
        max_mismatch_count,
        min_ref_count,
        min_alt_count,
        common_snps,
        panel_of_normals,
        somatic_snv_prior,
        germline_snv_prior,
        germline_indel_prior,
        threads,
        phase,
        non_human_sample,
        reference_sample,
        out_file, 
    )
    sbs2count = himut.mutlib.load_sbs96_counts(sbs_file, ref_file, chrom_lst)
    chrom2ref_tri2count = himut.reflib.get_ref_tricounts(refseq, chrom2chunkloci_lst, threads)
    himut.mutlib.dump_normcounts(
        cmdline, 
        sbs2count,
        chrom2ref_tri2count, 
        chrom2ccs_tri2count, 
        phase,
        tname2tsize,
        chrom2chunkloci_lst,
        out_file
    )
    if out_file.endswith(".tsv"):
        pdf_file = out_file.replace(".tsv", ".pdf")
    else:
        pdf_file = "{}.pdf".format(out_file)
    himut.mutlib.dump_norm_sbs96_plt(out_file, himut.bamlib.get_sample(bam_file), pdf_file)
    print("finished returning normcounts")
    cpu_end = time.time() / 60
    duration = cpu_end - cpu_start
    print("SBS96 count normalisation took {} minutes".format(duration))
    himut.util.exit()
