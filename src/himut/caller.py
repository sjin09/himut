import time
import pysam
import natsort
import himut.util
import himut.cslib
import himut.gtlib
import himut.bamlib
import himut.haplib
import himut.vcflib
import numpy as np
import multiprocessing as mp
from collections import defaultdict
from typing import Dict, List, Tuple


class METRICS:
    num_ccs: int = 0
    num_sbs: int = 0
    num_het_sbs: int = 0
    num_hetalt_sbs: int = 0
    num_homalt_sbs: int = 0
    num_somrev_sbs: int = 0
    num_homref_sbs: int = 0
    num_uncallable_sbs: int = 0
    num_low_gq_sbs: int = 0 
    num_low_bq_sbs: int = 0 
    num_pon_filtered_sbs: int = 0  
    num_pop_filtered_sbs: int = 0  
    num_md_filtered_sbs: int = 0
    num_ab_filtered_sbs: int = 0
    num_som: int = 0
    


def init_allelecounts():
    rpos2allelecounts = defaultdict(lambda: np.zeros(6))
    rpos2allele2bq_lst = defaultdict(lambda: {0: [], 1: [], 2: [], 3: [], 4: [], 5: []})
    rpos2allele2ccs_lst = defaultdict(
        lambda: {0: [], 1: [], 2: [], 3: [], 4: [], 5: []}
    )
    return rpos2allelecounts, rpos2allele2bq_lst, rpos2allele2ccs_lst, 


def update_allelecounts(
    ccs,
    rpos2allelecounts: Dict[int, np.ndarray],
    rpos2allele2bq_lst: Dict[int, Dict[int, List[int]]],
    rpos2allele2ccs_lst: Dict[int, Dict[int, List[str]]],
):

    tpos = ccs.tstart
    qpos = ccs.qstart
    for (state, ref, alt, ref_len, alt_len) in ccs.cstuple_lst:
        if state == 1:  # match
            for i, alt_base in enumerate(alt):
                epos = tpos + i
                bidx = himut.util.base2idx[alt_base]
                rpos2allelecounts[epos][bidx] += 1
                rpos2allele2ccs_lst[epos][bidx].append(ccs.qname)
                rpos2allele2bq_lst[epos][bidx].append(ccs.bq_int_lst[qpos + i])
        elif state == 2:  # sub
            bidx = himut.util.base2idx[alt]
            rpos2allelecounts[tpos][bidx] += 1
            rpos2allele2ccs_lst[tpos][bidx].append(ccs.qname)
            rpos2allele2bq_lst[tpos][bidx].append(ccs.bq_int_lst[qpos])
        elif state == 3:  # insertion
            rpos2allelecounts[tpos][4] += 1
        elif state == 4:  # deletion
            for j in range(len(ref[1:])):
                rpos2allelecounts[tpos + j][5] += 1
        tpos += ref_len
        qpos += alt_len


def is_ccs_phased(hap) -> bool:

    if hap == "0":
        return True
    elif hap == "1":
        return True
    else:
        return False


def is_low_qv(
    ccs,
    min_qv: int,
) -> bool:

    ccs.get_qv()
    if ccs.qv < min_qv:
        return True
    else:
        return False


def is_low_mapq(ccs_mapq: int, min_mapq: int):
    if ccs_mapq < min_mapq:
        return True
    else:
        return False


def is_chunk(tpos: int, start: int, end: int) -> bool:
    if start <= tpos and tpos <= end:
        return True
    else:
        return False


def is_germ_gt(
    som_gt: str,
    germ_gt: str,
    germ_gt_state: str,
    allelecounts: Dict[int, int],
) -> bool:

    if germ_gt_state == "het":
        if som_gt == germ_gt:
            return True
        else:
            return False
    elif germ_gt_state == "hetalt":
        base_sum = sum(
            [allelecounts[himut.util.base2idx[base]] for base in himut.util.base_lst]
        )
        base_counts = (
            allelecounts[himut.util.base2idx[germ_gt[0]]]
            + allelecounts[himut.util.base2idx[germ_gt[1]]]
        )
        if base_sum == base_counts and (
            som_gt[1] == germ_gt[0] or som_gt[1] == germ_gt[1]
        ):
            return True
        else:
            return False
    elif germ_gt_state == "homalt":
        ref_count = allelecounts[himut.util.base2idx[som_gt[0]]]
        if ref_count == 0 and germ_gt.count(som_gt[1]) == 2:
            return True
        else:
            return False
    elif germ_gt_state == "homref":
        if som_gt[1] == germ_gt[0]:
            return True
        else:
            return False


def is_low_gq(
    germ_gq: float,
    min_gq: int,
) -> bool:
    if germ_gq < min_gq:
        return True
    else:
        return False


def is_low_bq(
    alt: str, 
    min_bq: int, 
    allele2bq_lst: Dict[str, List[int]]
) -> bool:
    alt_bq_lst = allele2bq_lst[himut.util.base2idx[alt]]
    if alt_bq_lst.count(min_bq) == 0:
        return True
    else:
        return False


def get_hetalt_counts(
    p: str,
    q: str,
    read_depth: int,
    allelecounts: Dict[int, int],
    allele2bq_lst: Dict[str, List[int]],
):
    pidx = himut.util.base2idx[p]
    qidx = himut.util.base2idx[q]
    p_count = allelecounts[pidx]
    q_count = allelecounts[qidx]
    pbq = sum(allele2bq_lst[pidx])/float(p_count)
    qbq = sum(allele2bq_lst[qidx])/float(q_count)
    alt_bq = "{:0.1f},{:0.1f}".format(pbq, qbq)
    alt_count = "{:0.0f},{:0.0f}".format(p_count, q_count)
    alt_vaf = "{:.2f},{:.2f}".format(
        p_count/float(read_depth), q_count/float(read_depth)
    )
    return alt_bq, alt_vaf, alt_count


def is_chunk_phased(
    hap2count: Dict[str, int],
    min_hap_count: int,
) -> bool:

    h0_count = hap2count["0"]
    h1_count = hap2count["1"]
    if (h0_count >= min_hap_count and h1_count >= min_hap_count):
        return True
    else:
        return False


def get_somatic_substitutions(
    chrom: str,
    bam_file: str,
    common_snps: str,
    panel_of_normals: str,
    chunkloci_lst: List[Tuple[str, int, int]],
    phase_set2hbit_lst: Dict[str, List[str]],
    phase_set2hpos_lst: Dict[int, List[int]],
    phase_set2hetsnp_lst: Dict[int, List[Tuple[int, str, str]]], 
    min_qv: int,
    min_mapq: int,
    qlen_lower_limit: int,
    qlen_upper_limit: int,
    min_sequence_identity: float,
    min_gq: int,
    min_bq: int,
    min_trim: float,
    mismatch_window: int,
    max_mismatch_count: int,
    md_threshold: int,
    min_ref_count: int,
    min_alt_count: int,
    min_hap_count: int,
    somatic_snv_prior: float,
    germline_snv_prior: float,
    germline_indel_prior: float,
    phase: bool,
    non_human_sample: bool,
    create_panel_of_normals: bool,
    chrom2tsbs_lst: Dict[
        str, List[Tuple[str, int, str, str, int, int, int, int, str]]
    ],
    chrom2tsbs_log: Dict[str, List[int]]
) -> List[Tuple[str, int, str, str, int, int, int, float, float]]:

    ccs_seen = set()
    som_seen = set()
    pon_sbs_set = set()
    common_snp_set = set()
    himut.gtlib.init(germline_snv_prior)
    if not non_human_sample and common_snps.endswith(".vcf"):
        common_snp_set = himut.vcflib.load_common_snp(chrom, common_snps)

    if (
        not create_panel_of_normals
        and not non_human_sample
        and panel_of_normals.endswith(".vcf")
    ):
        pon_sbs_set = himut.vcflib.load_pon(chrom, panel_of_normals)

    m = METRICS()
    somatic_tsbs_lst = []
    filtered_somatic_tsbs_lst = []
    alignments = pysam.AlignmentFile(bam_file, "rb")
    for (chrom, chunk_start, chunk_end) in chunkloci_lst: ## TODO
        if not non_human_sample and common_snps.endswith(".bgz"):
            common_snp_set = himut.vcflib.load_bgz_common_snp(
                (
                    chrom,
                    (chunk_start - qlen_upper_limit),
                    (chunk_end + qlen_upper_limit),
                ),
                common_snps,
            )
        if (
            not non_human_sample
            and not create_panel_of_normals
            and panel_of_normals.endswith(".bgz")
        ):
            pon_sbs_set = himut.vcflib.load_bgz_pon(
                (
                    chrom,
                    (chunk_start - qlen_upper_limit),
                    (chunk_end + qlen_upper_limit),
                ),
                panel_of_normals,
            )
        if phase:
            phase_set = str(chunk_start)
            hbit_lst = phase_set2hbit_lst[phase_set] 
            hpos_lst = phase_set2hpos_lst[phase_set] 
            hetsnp_lst = phase_set2hetsnp_lst[phase_set] 
           
        somatic_tsbs_candidate_lst = []
        rpos2allelecounts, rpos2allele2bq_lst, rpos2allele2ccs_lst = init_allelecounts()
        for i in alignments.fetch(chrom, chunk_start, chunk_end): # iterate through reads
            ccs = himut.bamlib.BAM(i)
            if not ccs.is_primary:
                continue
            update_allelecounts(
                ccs, rpos2allelecounts, rpos2allele2bq_lst, rpos2allele2ccs_lst
            )
            if phase:
                ccs_hap = himut.haplib.get_ccs_hap(ccs, hbit_lst, hpos_lst, hetsnp_lst)  
                if not is_ccs_phased(ccs_hap):
                    continue 
            if is_low_qv(ccs, min_qv):
                continue
            if is_low_mapq(ccs.mapq, min_mapq):
                continue
            if ccs.get_blast_sequence_identity() < min_sequence_identity:
                continue
            if not (qlen_lower_limit < ccs.qlen and ccs.qlen < qlen_upper_limit):
                continue
            if ccs.qname not in ccs_seen:
                m.num_ccs += 1 
                ccs_seen.add(ccs.qname)
            ccs_somatic_tsbs_candidate_lst = ccs.get_tsbs_candidates(som_seen, min_trim, mismatch_window, max_mismatch_count)
            somatic_tsbs_candidate_lst.extend(ccs_somatic_tsbs_candidate_lst)

        for (tpos, ref, alt) in set(somatic_tsbs_candidate_lst): # iterate through substitutions
            if not is_chunk(tpos, chunk_start, chunk_end):
                continue

            m.num_sbs += 1
            rpos = tpos - 1 # reference FASTA file position 
            tsbs = (tpos, ref, alt)
            som_gt = "{}{}".format(ref, alt)
            allelecounts = rpos2allelecounts[rpos]
            allele2bq_lst = rpos2allele2bq_lst[rpos]
            del_count, ins_count, read_depth = himut.bamlib.get_read_depth(allelecounts)
            _, _, ref_count = himut.bamlib.get_ref_counts(ref, read_depth, allelecounts, allele2bq_lst)
            alt_bq, alt_vaf, alt_count = himut.bamlib.get_alt_counts(alt, read_depth, allelecounts, allele2bq_lst)
            germ_gt, _, germ_gt_state, gt2gt_state = himut.gtlib.get_germ_gt(ref, allele2bq_lst)
            if is_germ_gt(som_gt, germ_gt, germ_gt_state, allelecounts):
                if germ_gt_state == "het":
                    m.num_het_sbs += 1
                elif germ_gt_state == "hetalt":
                    m.num_hetalt_sbs += 1
                elif germ_gt_state == "homalt":
                    m.num_homalt_sbs += 1
                continue

            som_seen.add(tpos)
            germ_gq = himut.gtlib.get_germ_gq(som_gt, gt2gt_state, allele2bq_lst)
            if germ_gt_state == "het":
                m.num_somrev_sbs += 1
                filtered_somatic_tsbs_lst.append(
                    (
                        chrom,
                        tpos,
                        ref,
                        alt,
                        "HetSite",
                        germ_gq,
                        alt_bq,
                        read_depth,
                        ref_count,
                        alt_count,
                        alt_vaf,
                        ".",
                    )
                )
                continue
            elif germ_gt_state == "hetalt":
                m.num_somrev_sbs += 1
                a1, a2 = list(germ_gt)
                alt = "{},{}".format(a1, a2)
                alt_bq, alt_vaf, alt_count = get_hetalt_counts(
                    a1, a2, read_depth, allelecounts, allele2bq_lst
                )
                filtered_somatic_tsbs_lst.append(
                    (
                        chrom,
                        tpos,
                        ref,
                        alt,
                        "HetAltSite",
                        germ_gq,
                        alt_bq,
                        read_depth,
                        ref_count,
                        alt_count,
                        alt_vaf,
                        ".",
                    )
                )
                continue
            elif germ_gt_state == "homalt":
                m.num_somrev_sbs += 1
                filtered_somatic_tsbs_lst.append(
                    (
                        chrom,
                        tpos,
                        ref,
                        alt,
                        "HomAltSite",
                        germ_gq,
                        alt_bq,
                        read_depth,
                        ref_count,
                        alt_count,
                        alt_vaf,
                        ".",
                    )
                )
                continue

            if (del_count != 0 or ins_count != 0):
                m.num_uncallable_sbs += 1
                filtered_somatic_tsbs_lst.append(
                    (
                        chrom,
                        tpos,
                        ref,
                        alt,
                        "IndelSite",
                        germ_gq,
                        alt_bq,
                        read_depth,
                        ref_count,
                        alt_count,
                        alt_vaf,
                        ".",
                    )
                )
                continue

            m.num_homref_sbs += 1
            if is_low_gq(germ_gq, min_gq): # homref
                m.num_low_gq_sbs += 1
                filtered_somatic_tsbs_lst.append(
                    (
                        chrom,
                        tpos,
                        ref,
                        alt,
                        "LowGQ",
                        germ_gq,
                        alt_bq,
                        read_depth,
                        ref_count,
                        alt_count,
                        alt_vaf,
                        ".",
                    )
                )
                continue
            if is_low_bq(alt, min_bq, allele2bq_lst):
                m.num_low_bq_sbs += 1
                filtered_somatic_tsbs_lst.append(
                    (
                        chrom,
                        tpos,
                        ref,
                        alt,
                        "LowBQ",
                        germ_gq,
                        alt_bq,
                        read_depth,
                        ref_count,
                        alt_count,
                        alt_vaf,
                        ".",
                    )
                )
                continue
            if (
                tsbs in pon_sbs_set
                and not non_human_sample
                and not create_panel_of_normals
            ):
                m.num_pon_filtered_sbs += 1
                filtered_somatic_tsbs_lst.append(
                    (
                        chrom,
                        tpos,
                        ref,
                        alt,
                        "PanelOfNormal",
                        germ_gq,
                        alt_bq,
                        read_depth,
                        ref_count,
                        alt_count,
                        alt_vaf,
                        ".",
                    )
                )
                continue
            if tsbs in common_snp_set and not non_human_sample:
                m.num_pop_filtered_sbs += 1
                filtered_somatic_tsbs_lst.append(
                    (
                        chrom,
                        tpos,
                        ref,
                        alt,
                        "ComSnp",
                        germ_gq,
                        alt_bq,
                        read_depth,
                        ref_count,
                        alt_count,
                        alt_vaf,
                        ".",
                    )
                )
                continue
            if not (ref_count >= min_ref_count and alt_count >= min_alt_count):
                m.num_ab_filtered_sbs += 1
                filtered_somatic_tsbs_lst.append(
                    (
                        chrom,
                        tpos,
                        ref,
                        alt,
                        "LowDepth",
                        germ_gq,
                        alt_bq,
                        read_depth,
                        ref_count,
                        alt_count,
                        alt_vaf,
                        ".",
                    )
                )
                continue
            if read_depth > md_threshold:
                m.num_md_filtered_sbs += 1
                filtered_somatic_tsbs_lst.append(
                    (
                        chrom,
                        tpos,
                        ref,
                        alt,
                        "HighDepth",
                        germ_gq,
                        alt_bq,
                        read_depth,
                        ref_count,
                        alt_count,
                        alt_vaf,
                        ".",
                    )
                )
                continue

            if phase:
                m.num_som += 1
                som_hap_set = set()
                hap2count = defaultdict(lambda: 0)
                wt_ccs_set = set(rpos2allele2ccs_lst[rpos][himut.util.base2idx[ref]])
                alt_ccs_set = set(rpos2allele2ccs_lst[rpos][himut.util.base2idx[alt]])
                for j in alignments.fetch(chrom, tpos, tpos + 1):
                    ccs = himut.bamlib.BAM(j)
                    if not ccs.is_primary:
                        continue
                    if ccs.qname in wt_ccs_set:
                        ccs_hap = himut.haplib.get_ccs_hap(ccs, hbit_lst, hpos_lst, hetsnp_lst)  
                        hap2count[ccs_hap] += 1
                    elif ccs.qname in alt_ccs_set:
                        ccs_hap = himut.haplib.get_ccs_hap(ccs, hbit_lst, hpos_lst, hetsnp_lst)  
                        som_hap_set.add(ccs_hap)
                        
                som_hap_set = som_hap_set.difference(set(["."]))
                if is_chunk_phased(hap2count, min_hap_count) and len(som_hap_set) == 1:
                    somatic_tsbs_lst.append(
                        (
                            chrom,
                            tpos,
                            ref,
                            alt,
                            "PASS",
                            germ_gq,
                            alt_bq,
                            read_depth,
                            ref_count,
                            alt_count,
                            alt_vaf,
                            phase_set,
                        )
                    )
                else:
                    filtered_somatic_tsbs_lst.append(
                        (
                            chrom,
                            tpos,
                            ref,
                            alt,
                            "Unphased",
                            germ_gq,
                            alt_bq,
                            read_depth,
                            ref_count,
                            alt_count,
                            alt_vaf,
                            ".",
                        )
                    )
            else:
                m.num_som += 1
                somatic_tsbs_lst.append(
                    (
                        chrom,
                        tpos,
                        ref,
                        alt,
                        "PASS",
                        germ_gq,
                        alt_bq,
                        read_depth,
                        ref_count,
                        alt_count,
                        alt_vaf,
                        ".",
                    )
                )
                
    chrom2tsbs_lst[chrom] = natsort.natsorted(
        list(set(somatic_tsbs_lst + filtered_somatic_tsbs_lst))
    )
    chrom2tsbs_log[chrom] = [
        m.num_ccs,
        m.num_sbs,
        m.num_het_sbs,
        m.num_hetalt_sbs,
        m.num_homalt_sbs,
        m.num_somrev_sbs,
        m.num_homref_sbs,
        m.num_uncallable_sbs,
        m.num_low_gq_sbs,
        m.num_low_bq_sbs,
        m.num_pon_filtered_sbs,
        m.num_pop_filtered_sbs,
        m.num_md_filtered_sbs,
        m.num_ab_filtered_sbs,
        m.num_som,
    ]
    alignments.close()


def call_somatic_substitutions(
    bam_file: str,
    ref_file: str,
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
    create_panel_of_normals: bool,
    version: str,
    out_file: str,
) -> None:

    cpu_start = time.time() / 60
    chrom2ps2hbit_lst = defaultdict(dict)
    chrom2ps2hpos_lst = defaultdict(dict)
    chrom2ps2hetsnp_lst = defaultdict(dict)
    _, tname2tsize = himut.bamlib.get_tname2tsize(bam_file)
    chrom_lst, chrom2chunkloci_lst = himut.util.load_loci(region, region_list, tname2tsize)
    if phase:
        (
            chrom2ps2hbit_lst,
            chrom2ps2hpos_lst,
            chrom2ps2hetsnp_lst,
            chrom2chunkloci_lst,
        ) = himut.vcflib.load_phased_hetsnps(phased_vcf_file, chrom_lst, tname2tsize)
    qlen_lower_limit, qlen_upper_limit, md_threshold = himut.bamlib.get_thresholds(
        bam_file, chrom_lst, tname2tsize
    )
    himut.util.check_caller_input_exists(
        bam_file,
        ref_file,
        vcf_file,
        phased_vcf_file,
        common_snps,
        panel_of_normals,
        chrom_lst,
        tname2tsize,
        phase,
        non_human_sample,
        create_panel_of_normals,
        out_file,
    )
    if create_panel_of_normals:
        print("himut is calling substitutions for panel of normal preparation")
        (
            min_bq,
            min_qv,
            min_mapq,
            min_trim,
            min_hap_count,
            min_sequence_identity,
            phase,
        ) = himut.util.load_pon_params()

    if non_human_sample:
        germline_snv_prior, germline_indel_prior = himut.vcflib.get_germline_priors(
            chrom_lst, ref_file, vcf_file, reference_sample
        )

    vcf_header = himut.vcflib.get_himut_vcf_header(
        bam_file,
        vcf_file,
        phased_vcf_file,
        region,
        region_list,
        tname2tsize,
        common_snps,
        panel_of_normals,
        min_qv,
        min_mapq,
        qlen_lower_limit, 
        qlen_upper_limit, 
        min_sequence_identity,
        min_gq,
        min_bq,
        min_trim,
        mismatch_window,
        max_mismatch_count,
        md_threshold,
        min_ref_count,
        min_alt_count,
        min_hap_count,
        threads,
        somatic_snv_prior,
        germline_snv_prior,
        germline_indel_prior,
        phase,
        non_human_sample,
        reference_sample,
        create_panel_of_normals,
        version,
        out_file,
    )
    if phase:
        print(
            "himut is calling and phasing substitutions with {} threads".format(threads)
        )
    else:
        print("himut is calling substitutions with {} threads".format(threads))

    p = mp.Pool(threads)
    manager = mp.Manager()
    chrom2tsbs_lst = manager.dict()
    chrom2tsbs_log = manager.dict()
    get_somatic_substitutions_arg_lst = [
        (
            chrom,
            bam_file,
            common_snps,
            panel_of_normals,
            chrom2chunkloci_lst[chrom],
            chrom2ps2hbit_lst[chrom],
            chrom2ps2hpos_lst[chrom],
            chrom2ps2hetsnp_lst[chrom], 
            min_qv,
            min_mapq,
            qlen_lower_limit,
            qlen_upper_limit,
            min_sequence_identity,
            min_gq,
            min_bq,
            min_trim,
            mismatch_window,
            max_mismatch_count,
            md_threshold,
            min_ref_count,
            min_alt_count,
            min_hap_count,
            somatic_snv_prior,
            germline_snv_prior,
            germline_indel_prior,
            phase,
            non_human_sample,
            create_panel_of_normals,
            chrom2tsbs_lst,
            chrom2tsbs_log
        )
        for chrom in chrom_lst
    ]
    p.starmap(
        get_somatic_substitutions,
        get_somatic_substitutions_arg_lst,
    )
    p.close()
    p.join()

    if phase:
        himut.vcflib.dump_call_log(chrom_lst, chrom2tsbs_log)
        himut.vcflib.dump_phased_sbs(out_file, vcf_header, chrom_lst, chrom2tsbs_lst)
    else:
        himut.vcflib.dump_call_log(chrom_lst, chrom2tsbs_log)
        himut.vcflib.dump_sbs(out_file, vcf_header, chrom_lst, chrom2tsbs_lst)

    if phase:
        print(
            "himut finished phasing, calling and returning substitutions with {} threads".format(
                threads
            )
        )
    else:
        print(
            "himut finished calling and returning substitutions with {} threads".format(
                threads
            )
        )
    cpu_end = time.time() / 60
    duration = cpu_end - cpu_start
    print(
        "himut single molecule somatic mutation detection took {} minutes".format(
            duration
        )
    )
    himut.util.exit()
