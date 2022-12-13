import time
import math
import pysam
import bisect
import pyfastx
import himut.util
import himut.cslib
import himut.gtlib
import himut.haplib
import himut.bamlib
import himut.reflib
import numpy as np
import multiprocessing as mp
from collections import defaultdict
from typing import Dict, List, Tuple
from himut.mutlib import (
    purine,
    sbs2tri,
    tri_lst,
    purine2pyrimidine,
)


def get_allelecounts(
    allelecounts: np.ndarray
):
    read_depth = sum(allelecounts)  
    indel_count = allelecounts[4] + allelecounts[5] 
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


def get_cumsum_tricounts(chrom2tri2count: Dict[str, Dict[str, int]]):

    tri2count = defaultdict(lambda: 0)
    for chrom in chrom2tri2count:
        for tri in tri_lst:
            tri2count[tri] += chrom2tri2count[chrom][tri]            
    return tri2count 


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
    reference_sample: bool,
    non_human_sample: bool,
    chrom2ccs_tri2count: Dict[str, Dict[str, int]],
) -> Dict[str, int]:


    pon_sbs_set = set()
    common_snp_set = set()
    tri2count = defaultdict(lambda: 0)
    for tri in tri_lst: tri2count[tri] = 0
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
                common_snp_set = himut.vcflib.load_bgz_common_snp((chrom, chunk_start - qlen_upper_limit, chunk_end + qlen_upper_limit), common_snps)
            
            if panel_of_normals.endswith(".bgz"):
                pon_sbs_set = himut.vcflib.load_bgz_pon((chrom, chunk_start - qlen_upper_limit, chunk_end + qlen_upper_limit), panel_of_normals) 
       
        tpos2allele2bq_lst, tpos2allele2ccs_lst, tpos2allelecounts =  himut.util.init_allelecounts()
        for i in alignments.fetch(chrom, chunk_start, chunk_end):
            ccs = himut.bamlib.BAM(i)
            if not ccs.is_primary:
                continue
            
            himut.util.update_allelecounts(ccs, tpos2allelecounts, tpos2allele2bq_lst, tpos2allele2ccs_lst)
            if ccs.mapq < min_mapq:
                continue

            if ccs.qlen < qlen_lower_limit or ccs.qlen > qlen_upper_limit:
                continue

            if ccs.get_hq_base_proportion() < min_hq_base_proportion:
                continue

            if ccs.get_blast_sequence_identity() < min_sequence_identity:
                continue

            if ccs.get_query_alignment_proportion() < min_alignment_proportion:
                continue 

        for tpos in range(chunk_start, chunk_end):
            rpos = tpos - 1
            ref = chrom_seq[rpos]
            if ref == "N": 
                continue
            allelecounts = tpos2allelecounts[tpos]
            allele2bq_lst = tpos2allele2bq_lst[tpos]
            read_depth, indel_count = get_allelecounts(allelecounts) 
            if indel_count != 0:
                continue
            if read_depth < min_read_depth:
                continue
            if read_depth > md_threshold:               
                continue

            germ_gt, germ_gq, germ_gt_state = himut.gtlib.get_germline_gt(ref, allelecounts, allele2bq_lst, germline_snv_prior)
            if germ_gq < min_gq:
                continue
            if germ_gt_state == "het":
                continue
            elif germ_gt_state == "hetalt":
                continue
            elif germ_gt_state == "homalt":
                if reference_sample: # assembly error
                    continue

            gt_base = germ_gt[1]
            gt_base_count = allelecounts[himut.util.base2idx[gt_base]]
            gt_min_bq_count = allele2bq_lst[himut.util.base2idx[gt_base]].count(min_bq)
            if read_depth == gt_base_count:
                update_tricounts(rpos, gt_base, gt_min_bq_count, chrom_seq, tri2count)
                continue
            
            if germ_gt_state == "homref":
                alt_base_lst = himut.util.base_set.difference(gt_base) 
            else:
                alt_base_lst = himut.util.base_set.difference(list([ref, gt_base])) 

            for alt_base in alt_base_lst:
                state = 0
                tsbs = (tpos, ref, alt_base) 
                alt_ccs_set = set(tpos2allele2ccs_lst[tpos][himut.util.base2idx[alt_base]])
                if len(alt_ccs_set) == 0:
                    continue
                
                alt_min_bq_count = allele2bq_lst[himut.util.base2idx[alt_base]].count(min_bq)
                if alt_min_bq_count == 0:
                    continue

                if not non_human_sample:
                    if tsbs in pon_sbs_set:
                        continue
                    if tsbs in common_snp_set:
                        continue

                for j in alignments.fetch(chrom, tpos, tpos+1):
                    qccs = himut.bamlib.BAM(j)
                    if qccs.qname not in alt_ccs_set: 
                        continue
                    qccs.cs2subindel() 
                    qpos = qccs.qsbs_lst[qccs.tsbs_lst.index(tsbs)][0]
                    trimmed_qstart = math.floor(min_trim * qccs.qlen)
                    trimmed_qend = math.ceil((1 - min_trim) * qccs.qlen)
                    if qpos < trimmed_qstart:
                        state = 1
                        break     
                    elif qpos > trimmed_qend:
                        state = 1
                        break     
                    mpos_lst = [mismatch[0] for mismatch in qccs.mismatch_lst]
                    mismatch_start, mismatch_end = himut.util.get_mismatch_range(tpos, qpos, qccs.qlen, mismatch_window)
                    mismatch_count = bisect.bisect_right(mpos_lst, mismatch_end) - bisect.bisect_left(mpos_lst, mismatch_start) - 1
                    if mismatch_count > max_mismatch_count:
                        state = 1
                        break
                if state:
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
    min_mapq: int,
    min_trim: float,
    min_sequence_identity: float,
    min_hq_base_proportion: float,
    min_alignment_proportion: float,
    common_snps: str,
    panel_of_normals: str,
    min_gq: int,
    min_bq: int,
    mismatch_window: int,
    max_mismatch_count: int,
    min_ref_count: int,
    min_alt_count: int,
    somatic_snv_prior: float,
    germline_snv_prior: float,
    germline_indel_prior: float,
    threads: int,
    phase: bool,
    reference_sample: bool,
    non_human_sample: bool,
    out_file: str,
) -> None:

    cpu_start = time.time() / 60
    tname_lst, tname2tsize = himut.bamlib.get_tname2tsize(bam_file)
    chrom_lst, sbs2count = himut.util.check_normcounts_input_exists(
        bam_file,
        ref_file,
        sbs_file,
        vcf_file,
        phased_vcf_file,
        common_snps,
        panel_of_normals, 
        tname_lst,
        tname2tsize,
        phase,
        reference_sample,
        non_human_sample,
        out_file,
    )

    refseq = pyfastx.Fasta(ref_file)
    _, qlen_lower_limit, qlen_upper_limit, md_threshold = himut.bamlib.get_thresholds(
        bam_file, chrom_lst, tname2tsize
    )
    if phase:
        chrom2chunkloci_lst = himut.vcflib.get_chrom2hblock_loci(phased_vcf_file, chrom_lst, tname2tsize)
    else:
        chrom2chunkloci_lst = himut.reflib.load_seq_loci(ref_file, chrom_lst)
            
    if non_human_sample:
        germline_snv_prior, germline_indel_prior = himut.vcflib.get_germline_priors(chrom_lst, ref_file, vcf_file, reference_sample)

    print(
        "starting himut SBS96 count normalisation with {} threads".format(threads)
    )
    p = mp.Pool(threads)
    manager = mp.Manager()
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
            reference_sample,
            non_human_sample,
            chrom2ccs_tri2count,
        )
        for chrom in chrom_lst
    ]
    p.starmap(get_callable_tricounts, tricount_arg_lst)
    p.close()
    p.join()
    print(
        "finished himut SBS96 count normalisation with {} threads".format(threads)
    )
    chrom2ref_tri2count = himut.reflib.get_ref_tricounts(refseq, chrom2chunkloci_lst, threads) 
    ref_tricounts = get_cumsum_tricounts(chrom2ref_tri2count)
    ccs_tricounts = get_cumsum_tricounts(chrom2ccs_tri2count)
    himut.mutlib.dump_normcounts(sbs2count, ref_tricounts, ccs_tricounts, out_file)
    # himut.mutlib.dump_sbs96_plt(out_file, himut.bamlib.get_sample(bam_file), out_file.replace(".tsv", ".pdf")) 
    print("finished returning normcounts")
    cpu_end = time.time() / 60
    duration = cpu_end - cpu_start
    print("SBS96 count normalisation took {} minutes".format(duration))
    himut.util.exit()