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


def init_allelecounts():
    rpos2allelecounts = defaultdict(lambda: np.zeros(6))
    rpos2allele2bq_lst = defaultdict(lambda: {0: [], 1: [], 2: [], 3: [], 4: [], 5: []})
    rpos2allele2ccs_lst = defaultdict(
        lambda: {0: [], 1: [], 2: [], 3: [], 4: [], 5: []}
    )
    return rpos2allele2bq_lst, rpos2allele2ccs_lst, rpos2allelecounts


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


def is_hap_phased(
    hap2count: Dict[str, int],
    som_hap_count: int,
    min_hap_count: int,
) -> bool:

    h0_count = hap2count["0"]
    h1_count = hap2count["1"]
    if som_hap_count == 1 and (h0_count >= min_hap_count and h1_count >= min_hap_count):
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
    md_threshold: int,
    min_hap_count: int,
    somatic_snv_prior: float,
    germline_snv_prior: float,
    germline_indel_prior: float,
    phase: bool,
    non_human_sample: bool,
    create_panel_of_normals: bool,
    chrom2tsbs_lst: Dict[
        str, List[List[Tuple[str, int, str, str, int, int, int, int, str]]]
    ],
) -> List[Tuple[str, int, str, str, int, int, int, float, float]]:

    seen = set()
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

    somatic_tsbs_lst = []
    filtered_somatic_tsbs_lst = []
    alignments = pysam.AlignmentFile(bam_file, "rb")
    for (chrom, chunk_start, chunk_end) in chunkloci_lst: 
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

        somatic_tsbs_candidate_lst = []
        rpos2allele2bq_lst, rpos2allele2ccs_lst, rpos2allelecounts = init_allelecounts()
        for i in alignments.fetch(chrom, chunk_start, chunk_end):
            ccs = himut.bamlib.BAM(i)
            if not ccs.is_primary:
                continue
            himut.cslib.update_allelecounts(
                ccs, rpos2allelecounts, rpos2allele2bq_lst, rpos2allele2ccs_lst
            )
            if is_low_mapq(ccs.mapq, min_mapq):
                continue
            if ccs.qlen < qlen_lower_limit or ccs.qlen > qlen_upper_limit:
                continue
            if ccs.get_hq_base_proportion() < min_hq_base_proportion:
                continue
            if ccs.get_blast_sequence_identity() < min_sequence_identity:
                continue
            if ccs.get_query_alignment_proportion() < min_alignment_proportion:
                continue

            ccs.cs2mut()
            ccs_somatic_tsbs_candidate_lst = [
                tsbs for tsbs in ccs.tsbs_lst if tsbs[0] not in seen
            ]
            if len(ccs_somatic_tsbs_candidate_lst) == 0:
                continue
            somatic_tsbs_candidate_lst.extend(ccs_somatic_tsbs_candidate_lst)

        for tsbs in set(somatic_tsbs_candidate_lst):
            tpos, ref, alt = tsbs
            if tpos in seen:
                continue
            if not is_chunk(tpos, chunk_start, chunk_end):
                continue

            rpos = tpos - 1 # reference FASTA file position 
            som_gt = "{}{}".format(ref, alt)
            allelecounts = rpos2allelecounts[rpos]
            allele2bq_lst = rpos2allele2bq_lst[rpos]
            germ_gt, _, germ_gt_state, gt2gt_state = himut.gtlib.get_germ_gt(ref, allele2bq_lst)
            if is_germ_gt(som_gt, germ_gt, germ_gt_state, allelecounts):
                continue

            seen.add(tpos)
            germ_gq = himut.gtlib.get_germ_gq(som_gt, gt2gt_state, allele2bq_lst)
            (
                ref_count,
                alt_bq,
                alt_vaf,
                alt_count,
                indel_count,
                read_depth,
            ) = himut.bamlib.get_sbs_allelecounts(ref, alt, allelecounts, allele2bq_lst)
            if germ_gt_state == "het":
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
            if indel_count != 0:
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
            if is_low_gq(germ_gq, min_gq):
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
                filtered_somatic_tsbs_lst.append(
                    (
                        chrom,
                        tpos,
                        ref,
                        alt,
                        "ContRead",
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
            if ref_count < min_ref_count and alt_count < min_alt_count:
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
            elif ref_count >= min_ref_count and alt_count < min_alt_count:
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
            elif ref_count < min_ref_count and alt_count >= min_alt_count:
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

            som_state = 0
            wt_ccs_set = rpos2allele2ccs_lst[rpos][himut.util.base2idx[ref]]
            alt_ccs_set = set(rpos2allele2ccs_lst[rpos][himut.util.base2idx[alt]])
            trimmed_qstart, trimmed_qend = himut.bamlib.get_trimmed_range(ccs.qlen, min_trim)
            for ij in alignments.fetch(chrom, tpos, tpos + 1):
                ccs = himut.bamlib.BAM(ij)
                if ccs.qname not in alt_ccs_set:
                    continue
                ccs.cs2subindel()
                qpos = ccs.qsbs_lst[ccs.tsbs_lst.index(tsbs)][0]
                if himut.bamlib.is_trimmed(qpos, trimmed_qstart, trimmed_qend):
                    som_state = 1
                    filtered_somatic_tsbs_lst.append(
                        (
                            chrom,
                            tpos,
                            ref,
                            alt,
                            "Trimmed",
                            germ_gq,
                            alt_bq,
                            read_depth,
                            ref_count,
                            alt_count,
                            alt_vaf,
                            ".",
                        )
                    )
                    break
                if himut.bamlib.is_mismatch_conflict(
                    ccs, tpos, qpos, mismatch_window, max_mismatch_count
                ):
                    som_state = 1
                    filtered_somatic_tsbs_lst.append(
                        (
                            chrom,
                            tpos,
                            ref,
                            alt,
                            "MismatchConflict",
                            germ_gq,
                            alt_bq,
                            read_depth,
                            ref_count,
                            alt_count,
                            alt_vaf,
                            ".",
                        )
                    )
                    break
            if som_state:
                continue

            if phase:
                som_hap_set = set()
                phase_set = str(chunk_start)
                hap2count = defaultdict(lambda: 0)
                ccs_hap_lst = himut.haplib.get_loci_hap(
                    alignments,
                    (chrom, tpos, tpos + 1),
                    phase_set2hbit_lst[phase_set],
                    phase_set2hpos_lst[phase_set],
                    phase_set2hetsnp_lst[phase_set],
                )
                for qname, hap in ccs_hap_lst:
                    if hap == ".":
                        continue
                    if qname in wt_ccs_set:
                        hap2count[hap] += 1
                    if qname in alt_ccs_set:
                        som_hap_set.add(hap)

                if is_hap_phased(hap2count, len(som_hap_set), min_hap_count):
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
                    continue
            else:
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
    create_panel_of_normals: bool,
    version: str,
    out_file: str,
) -> None:

    cpu_start = time.time() / 60
    _, tname2tsize = himut.bamlib.get_tname2tsize(bam_file)
    chrom2ps2hbit_lst = defaultdict(dict)
    chrom2ps2hpos_lst = defaultdict(dict)
    chrom2ps2hetsnp_lst = defaultdict(dict)
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
            min_mapq,
            min_trim,
            min_sequence_identity,
            min_hq_base_proportion,
            min_alignment_proportion,
            min_bq,
            min_hap_count,
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
        md_threshold,
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
            min_hap_count,
            somatic_snv_prior,
            germline_snv_prior,
            germline_indel_prior,
            phase,
            non_human_sample,
            create_panel_of_normals,
            chrom2tsbs_lst,
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
        himut.vcflib.dump_phased_sbs(out_file, vcf_header, chrom_lst, chrom2tsbs_lst)
    else:
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
