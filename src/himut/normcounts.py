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
    num_het_bases: int = 0 
    num_hetalt_bases: int = 0 
    num_homalt_bases: int = 0 
    num_homref_bases: int = 0
    num_uncallable_bases: int = 0
    num_low_gq_bases: int = 0
    num_pon_filtered_bases: int = 0  
    num_pop_filtered_bases: int = 0  
    num_ab_filtered_bases: int = 0
    num_md_filtered_bases: int = 0
    num_callable_bases: int = 0


def init_allelecounts():
    rpos2basecounts = defaultdict(lambda: np.zeros(4))
    rpos2allelecounts = defaultdict(lambda: np.zeros(6))
    rpos2allele2bq_lst = defaultdict(lambda: {0: [], 1: [], 2: [], 3: [], 4: [], 5: []})
    return rpos2basecounts, rpos2allelecounts, rpos2allele2bq_lst, 


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


def get_read_depth(allelecounts: np.ndarray):
    read_depth = sum(allelecounts) - allelecounts[4]
    return read_depth


def get_indel_count(allelecounts: np.ndarray):
    indel_count = allelecounts[4] + allelecounts[5]
    return indel_count


def update_tricounts(
    tri2count: Dict[str, int],
    ccs_tri2count: Dict[str, int],
):
    for tri, count in tri2count.items():
        ccs_tri2count[tri] += count


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
    min_ref_count: int,
    min_alt_count: int,
    min_hap_count: int,
    md_threshold: int,
    somatic_snv_prior: float,
    germline_snv_prior: float,
    germline_indel_prior: float,
    phase: bool,
    non_human_sample: bool,
    chrom2ccs_tri2count: Dict[str, Dict[str, int]],
    chrom2norm_log: Dict[str, List[int]],
) -> Dict[str, int]:

    seen = set()
    pon_sbs_set = set() # init
    common_snp_set = set()
    himut.gtlib.init(germline_snv_prior)
    ccs_tri2count = defaultdict(lambda: 0)
    for tri in tri_lst:
        ccs_tri2count[tri] = 0
    alignments = pysam.AlignmentFile(bam_file, "rb")
    if not non_human_sample:
        if common_snps.endswith(".vcf"):
            common_snp_set = himut.vcflib.load_common_snp(chrom, common_snps)
        if panel_of_normals.endswith(".vcf"):
            pon_sbs_set = himut.vcflib.load_pon(chrom, panel_of_normals)

    m = METRICS()
    for (_chrom, chunk_start, chunk_end) in chunkloci_lst: # traverse reads 
        print(chrom, chunk_start, chunk_end)
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

        hap2count = defaultdict(lambda: 0) ## TODO
        rpos2tri2count = defaultdict(lambda: defaultdict(lambda: 0)) 
        rpos2basecounts, rpos2allelecounts, rpos2allele2bq_lst = init_allelecounts() 
        for i in alignments.fetch(chrom, chunk_start, chunk_end): # traverse reads 
            ccs = himut.bamlib.BAM(i)
            if not ccs.is_primary:
                continue
            # if phase: ## TODO
            #     ccs_hap = himut.haplib.get_ccs_hap(ccs, hbit_lst, hpos_lst, hetsnp_lst)  
            #     hap2count[ccs_hap] += 1 ## TODO

            update_allelecounts(ccs, rpos2allelecounts, rpos2allele2bq_lst)
            if himut.caller.is_low_mapq(ccs.mapq, min_mapq):
                continue
            if not (qlen_lower_limit < ccs.qlen and ccs.qlen < qlen_upper_limit):
                continue
            if ccs.get_hq_base_proportion() < min_hq_base_proportion:
                continue
            if ccs.get_blast_sequence_identity() < min_sequence_identity:
                continue
            if ccs.get_query_alignment_proportion() < min_alignment_proportion:
                continue
            if phase:
                ccs_hap = himut.haplib.get_ccs_hap(ccs, hbit_lst, hpos_lst, hetsnp_lst)  
                if not himut.caller.is_ccs_phased(ccs_hap):
                    continue
            if ccs.qname not in seen:
                m.num_ccs += 1
                seen.add(ccs.qname)
               
          
            counter = 0 
            ccs.cs2subindel()
            rpos = ccs.tstart
            qpos = ccs.qstart
            mismatch_tpos_lst = himut.bamlib.get_mismatch_positions(ccs) 
            trimmed_qstart, trimmed_qend = himut.bamlib.get_trimmed_range(ccs.qlen, min_trim)
            for cstuple in ccs.cstuple_lst:
                state, ref, alt, ref_len, alt_len = cstuple
                if state == 1:  # match
                    mismatch_tstart, mismatch_tend = himut.bamlib.get_mismatch_range(rpos, qpos, ccs.qlen, mismatch_window)
                    for j, alt_base in enumerate(alt):
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
                        counter += 1
                        rpos2basecounts[rpos+j][himut.util.base2idx[alt_base]] += 1
                        tri = get_tri_context(seq, rpos+j)               
                        rpos2tri2count[rpos+j][tri] += 1
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
                    counter += 1
                    rpos2basecounts[rpos][himut.util.base2idx[alt]] += 1
                    tri = get_tri_context(seq, rpos)               
                    rpos2tri2count[rpos][tri] += 1
                rpos += ref_len
                qpos += alt_len 
            
        # if phase: ## TODO 
        #     if not himut.caller.is_chunk_phased(hap2count, min_hap_count): 
        #        continue 
           
        for rpos in range(chunk_start, chunk_end):  
            ref = seq[rpos]
            if ref == "N":
                continue
            if rpos == 1:
                continue
            base_sum = sum(rpos2basecounts[rpos])
            allele2bq_lst = rpos2allele2bq_lst[rpos]
            _, germ_gq, germ_gt_state, gt2gt_state = himut.gtlib.get_germ_gt(
                ref, allele2bq_lst
            )
            
            m.num_bases += base_sum
            if germ_gt_state == "het":
                m.num_het_bases += base_sum
                # print(chrom, rpos, "het", allele2bq_lst)
                continue
            elif germ_gt_state == "hetalt":
                m.num_hetalt_bases += base_sum
                # print(chrom, rpos, "hetalt", allele2bq_lst)
                continue
            elif germ_gt_state == "homalt":
                m.num_homalt_bases += base_sum
                # print(chrom, rpos, "homalt", allele2bq_lst)
                continue

            m.num_homref_bases += base_sum
            allelecounts = rpos2allelecounts[rpos]
            read_depth = get_read_depth(allelecounts)
            indel_count = get_indel_count(allelecounts)
            if indel_count != 0:
                m.num_uncallable_bases += base_sum 
                continue
            if himut.caller.is_low_gq(germ_gq, min_gq):
                m.num_low_gq_bases += base_sum
                continue
            if read_depth > md_threshold:
                m.num_md_filtered_bases += base_sum 
                continue
  
            ridx = himut.util.base2idx[ref]
            tri2count = rpos2tri2count[rpos]
            ref_allelecount = allelecounts[ridx]
            if read_depth == ref_allelecount: # implict 
                if ref_allelecount < min_ref_count:
                    # print(chrom, rpos, "allele imbalance", ref_allelecount)
                    m.num_ab_filtered_bases += base_sum
                    continue
                m.num_callable_bases += base_sum
                update_tricounts(tri2count, ccs_tri2count)
            else:
                alt_state = 0
                alt_basecount_lst = []
                basecounts = rpos2basecounts[rpos]
                alt_lst = list(himut.util.base_set.difference(ref))
                aidx_lst = [himut.util.base2idx[alt] for alt in alt_lst]
                for alt, aidx in zip(alt_lst, aidx_lst):
                    tsbs = (rpos + 1, ref, alt)
                    alt_basecount_lst.append(basecounts[aidx])
                    if tsbs in pon_sbs_set and not non_human_sample:
                        alt_state = 1
                        m.num_pon_filtered_bases += base_sum
                        break
                    if tsbs in common_snp_set and not non_human_sample:
                        alt_state = 1
                        m.num_pop_filtered_bases += base_sum
                        break
                if alt_state:
                    continue
               
                alt = alt_lst[alt_basecount_lst.index(max(alt_basecount_lst))]
                germ_gq = himut.gtlib.get_germ_gq(alt, gt2gt_state, allele2bq_lst)
                if himut.caller.is_low_gq(germ_gq, min_gq):
                    m.num_low_gq_bases += base_sum
                    continue

                alt_allelecount = allelecounts[himut.util.base2idx[alt]]
                if not (ref_allelecount >= min_ref_count and alt_allelecount >= min_alt_count):
                    m.num_ab_filtered_bases += 1
                    # print(chrom, rpos, "allele imbalance", ref_allelecount, alt_allelecount) ## TODO
                    continue
                # print(chrom, rpos, ref, alt, "PASS") ## TODO
                update_tricounts(tri2count, ccs_tri2count)
                
    chrom2ccs_tri2count[chrom] = dict(ccs_tri2count) # return
    chrom2norm_log[chrom] = [
        m.num_ccs,
        m.num_bases,
        m.num_het_bases,
        m.num_hetalt_bases,
        m.num_homalt_bases,
        m.num_homref_bases,
        m.num_uncallable_bases,
        m.num_low_gq_bases,
        m.num_md_filtered_bases,
        m.num_ab_filtered_bases,
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
    qlen_lower_limit, qlen_upper_limit, md_threshold = himut.vcflib.get_thresholds(sbs_file)       
    if non_human_sample:
        germline_snv_prior, germline_indel_prior = himut.vcflib.get_germline_priors(
            chrom_lst, ref_file, vcf_file, reference_sample
        )

    p = mp.Pool(threads)
    manager = mp.Manager()
    refseq = pyfastx.Fasta(ref_file)
    chrom2norm_log = manager.dict()
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
            min_hap_count,
            md_threshold,
            somatic_snv_prior,
            germline_snv_prior,
            germline_indel_prior,
            phase,
            non_human_sample,
            chrom2ccs_tri2count,
            chrom2norm_log
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
