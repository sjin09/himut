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
from dataclasses import dataclass
from collections import defaultdict
from typing import Set, Dict, List, Tuple
from himut.mutlib import (
    purine,
    tri_lst,
    purine2pyrimidine,
)


@dataclass
class METRICS:
    num_ccs: int = 0
    num_bases: int = 0
    num_unphased_bases: int = 0
    num_het_bases: int = 0 
    num_hetalt_bases: int = 0 
    num_homalt_bases: int = 0 
    num_homref_bases: int = 0
    num_uncallable_bases: int = 0
    num_md_filtered_bases: int = 0
    num_ab_filtered_bases: int = 0
    num_low_gq_bases: int = 0
    num_pon_filtered_bases: int = 0  
    num_pop_filtered_bases: int = 0  
    num_callable_bases: int = 0


def init_allelecounts():
    rpos2allelecounts = defaultdict(lambda: np.zeros(6))
    rpos2allele2bq_lst = defaultdict(lambda: {0: [], 1: [], 2: [], 3: []})
    return rpos2allelecounts, rpos2allele2bq_lst

        
def get_tri_context(
    seq: str,
    pos: int,
):

    tri = seq[pos-1:pos+2]    
    if len(tri) == 3: 
        if tri[1] in purine:
            tri_pyr = "".join(
                [purine2pyrimidine.get(base, "N") for base in tri[::-1]]
            )
            return tri_pyr
        return tri
    return "NNN"

         
def update_tri2count(
    ccs,
    min_bq: int,
    min_trim: float, 
    mismatch_window: int, 
    max_mismatch_count: int, 
    rpos2count: Dict[int, Dict[str, int]]
):
                
    ccs.cs2subindel()
    rpos = ccs.tstart
    qpos = ccs.qstart
    mismatch_tpos_lst = himut.bamlib.get_mismatch_positions(ccs) 
    trimmed_qstart, trimmed_qend = himut.bamlib.get_trimmed_range(ccs.qlen, min_trim)
    for cstuple in ccs.cstuple_lst:
        state, ref, _alt, ref_len, alt_len = cstuple
        if state == 1:  # match
            mismatch_tstart, mismatch_tend = himut.bamlib.get_mismatch_range(rpos, qpos, ccs.qlen, mismatch_window)
            for j, _ref_base in enumerate(ref):
                mismatch_count = (
                    bisect.bisect_right(mismatch_tpos_lst, mismatch_tend + j)
                    - bisect.bisect_left(mismatch_tpos_lst, mismatch_tstart + j)
                )
                if is_low_bq(ccs.bq_int_lst[qpos+j], min_bq):
                    continue
                if is_mismatch_conflict(mismatch_count, max_mismatch_count):
                    continue
                if himut.bamlib.is_trimmed(qpos + j, trimmed_qstart, trimmed_qend): # 
                    continue
                rpos2count[rpos+j] += 1
        elif state == 2:  
            mismatch_tstart, mismatch_tend = himut.bamlib.get_mismatch_range(rpos, qpos, ccs.qlen, mismatch_window)
            mismatch_count = (
                bisect.bisect_right(mismatch_tpos_lst, mismatch_tend)
                - bisect.bisect_left(mismatch_tpos_lst, mismatch_tstart)
                - 1
            )
            if is_low_bq(ccs.bq_int_lst[qpos], min_bq):
                pass
            if is_mismatch_conflict(mismatch_count, max_mismatch_count):
                pass
            if himut.bamlib.is_trimmed(qpos, trimmed_qstart, trimmed_qend):
                pass
            rpos2count[rpos] += 1
        rpos += ref_len
        qpos += alt_len 
            


def update_allelecounts(
    ccs,
    rpos2allelecounts: Dict[int, np.ndarray],
    rpos2allele2bq_lst: Dict[int, Dict[int, List[int]]],
):

    tpos = ccs.tstart
    qpos = ccs.qstart
    for (state, ref, alt, ref_len, alt_len) in ccs.cstuple_lst:
        if state == 1:  # match
            for i, alt_base in enumerate(alt):
                epos = tpos + i
                bidx = himut.util.base2idx[alt_base]
                rpos2allelecounts[epos][bidx] += 1
                rpos2allele2bq_lst[epos][bidx].append(ccs.bq_int_lst[qpos + i])
        elif state == 2:  # sub
            bidx = himut.util.base2idx[alt]
            rpos2allelecounts[tpos][bidx] += 1
            rpos2allele2bq_lst[tpos][bidx].append(ccs.bq_int_lst[qpos])
        elif state == 3:  # insertion
            rpos2allelecounts[tpos][4] += 1
        elif state == 4:  # deletion
            for j in range(len(ref[1:])):
                rpos2allelecounts[tpos + j][5] += 1
        tpos += ref_len
        qpos += alt_len



def update_phased_allelecounts(
    ccs,
    rpos2hapcount: Dict[str, int],
    rpos2allelecounts: Dict[int, np.ndarray],
    rpos2allele2bq_lst: Dict[int, Dict[int, List[int]]],
):

    tpos = ccs.tstart
    qpos = ccs.qstart
    for (state, ref, alt, ref_len, alt_len) in ccs.cstuple_lst:
        if state == 1:  # match
            for i, alt_base in enumerate(alt):
                epos = tpos + i
                bidx = himut.util.base2idx[alt_base]
                rpos2hapcount[epos][ccs.hap] += 1
                rpos2allelecounts[epos][bidx] += 1
                rpos2allele2bq_lst[epos][bidx].append(ccs.bq_int_lst[qpos + i])
        elif state == 2:  # sub
            bidx = himut.util.base2idx[alt]
            rpos2hapcount[tpos][ccs.hap] += 1
            rpos2allelecounts[tpos][bidx] += 1
            rpos2allele2bq_lst[tpos][bidx].append(ccs.bq_int_lst[qpos])
        elif state == 3:  # insertion
            rpos2allelecounts[tpos][4] += 1
        elif state == 4:  # deletion
            for j in range(len(ref[1:])):
                rpos2allelecounts[tpos + j][5] += 1
        tpos += ref_len
        qpos += alt_len


def is_low_bq(
    bq: int, 
    min_bq: int, 
) -> bool:

    if bq < min_bq:
        return True
    else:
        return False


def is_mismatch_conflict(mismatch_count: int, max_mismatch_count: int) -> bool:

    if mismatch_count > max_mismatch_count:
        return True
    else:
        return False


def is_rpos_phased(
    hap2count: Dict[str, int],
    min_hap_count: int,
) -> bool:

    h0_count = hap2count["0"]
    h1_count = hap2count["1"]
    if (h0_count >= min_hap_count and h1_count >= min_hap_count):
        return True
    else:
        return False


def get_callable_tricounts(
    chrom: str,
    seq: str,
    bam_file: str,
    common_snps: str,
    panel_of_normals: str,
    chunkloci_lst: List[Tuple[str, int, int]],
    phase_set2hbit_lst: Dict[str, List[str]],
    phase_set2hpos_lst: Dict[int, List[int]],
    phase_set2hetsnp_lst: Dict[int, List[Tuple[int, str, str]]], 
    min_qv: int,
    min_mapq: int,
    min_trim: float,
    qlen_lower_limit: int,
    qlen_upper_limit: int,
    min_sequence_identity: float,
    min_gq: int,
    min_bq: int,
    mismatch_window: int,
    max_mismatch_count: int,
    min_ref_count: int,
    min_alt_count: int,
    min_hap_count: int,
    md_threshold: int,
    somatic_snv_prior: float,
    germline_snv_prior: float,
    germline_indel_prior: float,
    phase: bool,
    non_human_sample: bool,
    chrom2ccs_callable_tri2count: Dict[str, Dict[str, int]],
    chrom2ref_callable_tri2count: Dict[str, Dict[str, int]],
    chrom2norm_log: Dict[str, List[int]],
) -> Dict[str, int]:

    seen = set() # init
    pon_sbs_set = set() 
    common_snp_set = set()
    himut.gtlib.init(germline_snv_prior)
    ccs_tri2count = defaultdict(lambda: 0)
    ref_tri2count = defaultdict(lambda: 0)
    for tri in tri_lst:
        ccs_tri2count[tri] = 0
        ref_tri2count[tri] = 0
    alignments = pysam.AlignmentFile(bam_file, "rb")
    if not non_human_sample:
        if common_snps.endswith(".vcf"):
            common_snp_set = himut.vcflib.load_common_snp(chrom, common_snps)
        if panel_of_normals.endswith(".vcf"):
            pon_sbs_set = himut.vcflib.load_pon(chrom, panel_of_normals)

    m = METRICS()
    for (_chrom, chunk_start, chunk_end) in chunkloci_lst: # traverse reads 
        if not non_human_sample: # load
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

        if phase:
            phase_set = str(chunk_start)
            hbit_lst = phase_set2hbit_lst[phase_set] 
            hpos_lst = phase_set2hpos_lst[phase_set] 
            hetsnp_lst = phase_set2hetsnp_lst[phase_set] 

        rpos2count = defaultdict(lambda: 0) 
        rpos2allelecounts, rpos2allele2bq_lst= init_allelecounts() 
        rpos2hap2count = defaultdict(lambda: defaultdict(lambda: 0))
        for i in alignments.fetch(chrom, chunk_start, chunk_end): # iterate through reads 
            ccs = himut.bamlib.BAM(i)
            if not ccs.is_primary:
                continue
            if phase: 
                ccs.hap = himut.haplib.get_ccs_hap(ccs, hbit_lst, hpos_lst, hetsnp_lst)  
                update_phased_allelecounts(ccs, rpos2hap2count, rpos2allelecounts, rpos2allele2bq_lst)
                if not himut.caller.is_ccs_phased(ccs.hap):
                    continue
            else:
                update_allelecounts(ccs, rpos2allelecounts, rpos2allele2bq_lst)
                
            if himut.caller.is_low_qv(ccs, min_qv):
                continue
            if himut.caller.is_low_mapq(ccs.mapq, min_mapq):
                continue
            if ccs.get_blast_sequence_identity() < min_sequence_identity:
                continue
            if not (qlen_lower_limit < ccs.qlen and ccs.qlen < qlen_upper_limit):
                continue
            if ccs.qname not in seen:
                m.num_ccs += 1
                seen.add(ccs.qname)
            update_tri2count(ccs, min_bq, min_trim, mismatch_window, max_mismatch_count, rpos2count)
                
        for rpos in range(chunk_start, chunk_end): # iterate through reference positions
            ref = seq[rpos]
            tri_sum = rpos2count[rpos]
            if not ref in himut.util.base_set:
                continue
            if tri_sum == 0:
                continue

            m.num_bases += tri_sum
            if phase:
                hap2count = rpos2hap2count[rpos]
                if not is_rpos_phased(hap2count, min_hap_count):
                    m.num_unphased_bases += tri_sum
                    continue
                
            ridx = himut.util.base2idx[ref]
            allelecounts = rpos2allelecounts[rpos]
            allele2bq_lst = rpos2allele2bq_lst[rpos]
            del_count, ins_count, read_depth = himut.bamlib.get_read_depth(allelecounts)
            _, germ_gq, germ_gt_state, gt2gt_state = himut.gtlib.get_germ_gt(ref, allele2bq_lst)
            if germ_gt_state == "het":
                m.num_het_bases += tri_sum
                continue
            elif germ_gt_state == "hetalt":
                m.num_hetalt_bases += tri_sum
                continue
            elif germ_gt_state == "homalt":
                m.num_homalt_bases += tri_sum
                continue

            m.num_homref_bases += tri_sum
            if (del_count != 0 or ins_count != 0):
                m.num_uncallable_bases += tri_sum 
                continue
            if read_depth > md_threshold:
                m.num_md_filtered_bases += tri_sum 
                continue
 
            ref_count = allelecounts[ridx]
            tri = get_tri_context(seq, rpos)     
            if read_depth == ref_count: # implict 
                if himut.caller.is_low_gq(germ_gq, min_gq):
                    m.num_low_gq_bases += tri_sum
                    continue
                if ref_count < min_ref_count:
                    m.num_ab_filtered_bases += tri_sum
                    continue
                ref_tri2count[tri] += 1
                ccs_tri2count[tri] += tri_sum
                m.num_callable_bases += tri_sum
            else:
                alt_state = 0
                alt_count_lst = []
                alt_lst = list(himut.util.base_set.difference(ref))
                aidx_lst = [himut.util.base2idx[alt] for alt in alt_lst]
                for alt, aidx in zip(alt_lst, aidx_lst):
                    tsbs = (rpos + 1, ref, alt)
                    alt_count = allelecounts[aidx]
                    alt_count_lst.append(alt_count)
                    if alt_count == 0:
                        continue
                    if tsbs in pon_sbs_set and not non_human_sample:
                        alt_state = 1
                        m.num_pon_filtered_bases += tri_sum
                        break
                    if tsbs in common_snp_set and not non_human_sample:
                        alt_state = 1
                        m.num_pop_filtered_bases += tri_sum
                        break
                if alt_state == 1:
                    continue
                alt = alt_lst[alt_count_lst.index(max(alt_count_lst))]
                germ_gq = himut.gtlib.get_germ_gq(alt, gt2gt_state, allele2bq_lst)
                if himut.caller.is_low_gq(germ_gq, min_gq):
                    m.num_low_gq_bases += tri_sum
                    continue

                alt_count = allelecounts[himut.util.base2idx[alt]]
                if not (ref_count >= min_ref_count and alt_count >= min_alt_count):
                    m.num_ab_filtered_bases += tri_sum
                    continue
                ref_tri2count[tri] += 1
                ccs_tri2count[tri] += tri_sum
                m.num_callable_bases += tri_sum
                
    chrom2ccs_callable_tri2count[chrom] = dict(ccs_tri2count) # return
    chrom2ref_callable_tri2count[chrom] = dict(ref_tri2count) # return
    chrom2norm_log[chrom] = [
        m.num_ccs,
        m.num_bases,
        m.num_unphased_bases,
        m.num_het_bases,
        m.num_hetalt_bases,
        m.num_homalt_bases,
        m.num_homref_bases,
        m.num_uncallable_bases,
        m.num_md_filtered_bases,
        m.num_ab_filtered_bases,
        m.num_low_gq_bases,
        m.num_pon_filtered_bases,  
        m.num_pop_filtered_bases,  
        m.num_callable_bases,
    ]
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
    min_qv: int,
    min_mapq: int,
    min_sequence_identity: float,
    min_gq: int,
    min_bq: int,
    min_trim: float,
    mismatch_window: int,
    max_mismatch_count: int,
    min_ref_count: int,
    min_alt_count: int,
    min_hap_count: int,
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
    chrom_lst, chrom2chunkloci_lst = himut.util.load_loci(region, region_list, tname2tsize)
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

    if non_human_sample:
        germline_snv_prior, germline_indel_prior = himut.vcflib.get_germline_priors(
            chrom_lst, ref_file, vcf_file, reference_sample
        )


    qlen_lower_limit, qlen_upper_limit, md_threshold = himut.vcflib.get_thresholds(sbs_file)       

    p = mp.Pool(threads)
    manager = mp.Manager()
    refseq = pyfastx.Fasta(ref_file)
    chrom2norm_log = manager.dict()
    chrom2ccs_callable_tri2count = manager.dict()
    chrom2ref_callable_tri2count = manager.dict()
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
            min_qv,
            min_mapq,
            min_trim,
            qlen_lower_limit,
            qlen_upper_limit,
            min_sequence_identity,
            min_gq,
            min_bq,
            mismatch_window,
            max_mismatch_count,
            min_ref_count,
            min_alt_count,
            min_hap_count,
            md_threshold,
            somatic_snv_prior,
            germline_snv_prior,
            germline_indel_prior,
            phase,
            non_human_sample,
            chrom2ccs_callable_tri2count,
            chrom2ref_callable_tri2count,
            chrom2norm_log
        )
        for chrom in chrom_lst
    ]
    p.starmap(get_callable_tricounts, tricount_arg_lst) 
    p.close()
    p.join()
    print("finished himut SBS96 count normalisation with {} threads".format(threads))

    cmdline = himut.mutlib.get_normcounts_cmdline( # return
        bam_file,
        ref_file,
        sbs_file,
        vcf_file,
        phased_vcf_file,
        min_qv,
        min_mapq,
        min_sequence_identity,
        min_gq,
        min_bq,
        min_trim,
        mismatch_window,
        max_mismatch_count,
        min_ref_count,
        min_alt_count,
        min_hap_count,
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
    ref_tri2count = himut.reflib.get_genome_tricounts(refseq, chrom_lst, threads)
    himut.mutlib.dump_normcounts(
        sbs2count,
        ref_tri2count,
        chrom2ref_callable_tri2count, 
        chrom2ccs_callable_tri2count, 
        cmdline, 
        out_file
    )
    himut.mutlib.dump_norm_log(
        chrom_lst, 
        chrom2norm_log,
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
