import time
import pysam
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
from typing import Dict, List, Tuple
from himut.mutlib import (
    purine,
    tri_lst,
    purine2pyrimidine,
)


def get_allelecounts(allelecounts: np.ndarray):
    ins_count = allelecounts[4]
    del_count = allelecounts[5]
    indel_count = ins_count + del_count
    read_depth = sum(allelecounts) - ins_count
    return read_depth, indel_count


def update_tricounts(
    pos: int,
    base: str,
    base_bq_count: int,
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
    tri2count[tri] += base_bq_count


def get_callable_tricounts(
    chrom: str,
    chrom_seq: str,
    bam_file: str,
    common_snps: str,
    panel_of_normals: str,
    chunkloci_lst: List[Tuple[str, int, int]],
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

    for (chrom, chunk_start, chunk_end) in chunkloci_lst: 
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

        (
            tpos2allele2bq_lst,
            tpos2allele2ccs_lst,
            tpos2allelecounts,
        ) = himut.util.init_allelecounts()
        for i in alignments.fetch(chrom, chunk_start, chunk_end):
            ccs = himut.bamlib.BAM(i)
            if not ccs.is_primary:
                continue
            himut.cslib.update_allelecounts(
                ccs, tpos2allelecounts, tpos2allele2bq_lst, tpos2allele2ccs_lst
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

        for tpos in range(chunk_start, chunk_end):  # 1-based coordinate
            rpos = tpos - 1  # 0-based coordinate
            ref = chrom_seq[rpos]
            if ref == "N":
                continue
            allelecounts = tpos2allelecounts[tpos]
            allele2bq_lst = tpos2allele2bq_lst[tpos]
            read_depth, indel_count = get_allelecounts(allelecounts)
            germ_gt, germ_gq, germ_gt_state = himut.gtlib.get_germ_gt(
                ref, allele2bq_lst
            )
            if himut.caller.is_low_gq(germ_gq, min_gq):
                continue
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

            gt_base = germ_gt[1]
            gt_base_count = allelecounts[himut.util.base2idx[gt_base]]
            gt_min_bq_count = allele2bq_lst[himut.util.base2idx[gt_base]].count(min_bq)
            if (
                read_depth == gt_base_count
            ):  # homozygous reference without single molecule single-base-substitution
                update_tricounts(rpos, gt_base, gt_min_bq_count, chrom_seq, tri2count)
                continue

            alt_base_lst = himut.util.base_set.difference(gt_base)
            for alt_base in alt_base_lst:
                som_state = 0
                tsbs = (tpos, ref, alt_base)
                alt_ccs_set = set(
                    tpos2allele2ccs_lst[tpos][himut.util.base2idx[alt_base]]
                )
                alt_min_bq_count = allele2bq_lst[himut.util.base2idx[alt_base]].count(
                    min_bq
                )
                if len(alt_ccs_set) == 0:
                    continue
                if alt_min_bq_count == 0:
                    continue
                if not non_human_sample:
                    if tsbs in pon_sbs_set:
                        continue
                    if tsbs in common_snp_set:
                        continue
                for ij in alignments.fetch(chrom, tpos, tpos + 1):
                    ccs = himut.bamlib.BAM(ij)
                    if ccs.qname not in alt_ccs_set:
                        continue
                    ccs.cs2subindel()
                    qpos = ccs.qsbs_lst[ccs.tsbs_lst.index(tsbs)][0]
                    if himut.bamlib.is_trimmed(ccs, qpos, min_trim):
                        som_state = 1
                        break
                    if himut.bamlib.is_mismatch_conflict(
                        ccs, tpos, qpos, mismatch_window, max_mismatch_count
                    ):
                        som_state = 1
                        break
                if som_state:
                    continue
                update_tricounts(rpos, alt_base, alt_min_bq_count, chrom_seq, tri2count)
            update_tricounts(rpos, gt_base, gt_min_bq_count, chrom_seq, tri2count)
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
        chrom2chunkloci_lst = himut.vcflib.get_chrom2hblock_loci(
            phased_vcf_file, chrom_lst, tname2tsize
        )
    else:
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
    himut.mutlib.dump_norm_sbs96_plt(
        out_file, himut.bamlib.get_sample(bam_file), "{}.pdf".format(out_file)
    )
    print("finished returning normcounts")
    cpu_end = time.time() / 60
    duration = cpu_end - cpu_start
    print("SBS96 count normalisation took {} minutes".format(duration))
    himut.util.exit()
