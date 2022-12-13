import time
import math
import pysam
import bisect
import natsort
import himut.util
import himut.cslib
import himut.gtlib
import himut.dbslib
import himut.haplib
import himut.mutlib
import himut.bamlib
import himut.vcflib
import multiprocessing as mp
from collections import defaultdict
from typing import Dict, List, Tuple


def is_germ_gt(
    som_gt: str,
    germ_gt: str
):
    if som_gt == germ_gt:
        return True
    else:
        return False 


def get_n_choose_k(
    k: int,
    n: int, 
) -> int:
    n_choose_k = math.factorial(n)/(math.factorial(k)*math.factorial((n-k))) 
    return n_choose_k


def get_somatic_substitutions(
    chrom: str,
    chrom_len: int,
    bam_file: str,
    phased_vcf_file: str,
    common_snps: str,
    panel_of_normals: str,
    loci_lst: List[Tuple[str, int, int]],
    min_mapq: int,
    min_trim: float,
    qlen_mean: int,
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
    reference_sample: bool,
    create_panel_of_normals: bool,
    chrom2tsbs_lst: Dict[str, List[List[Tuple[str, int, str, str, int, int, int, int, str]]]],
) -> List[Tuple[str, int, str, str, int, int, int, float, float]]:

    seen = set()
    pon_sbs_set = set()
    common_snp_set = set()
    if not non_human_sample and common_snps.endswith(".vcf"):
        common_snp_set = himut.vcflib.load_common_snp(chrom, common_snps)

    if not create_panel_of_normals and not non_human_sample and panel_of_normals.endswith(".vcf"):
        pon_sbs_set = himut.vcflib.load_pon(chrom, panel_of_normals)
         
    if phase:
        (
            hpos_lst,
            hpos2phase_set,
            phase_set2hbit_lst,
            phase_set2hpos_lst,
            phase_set2hetsnp_lst
        ) = himut.vcflib.load_phased_hetsnps(phased_vcf_file, chrom, chrom_len)

    if phase:
        chunkloci_lst = himut.haplib.phase_set2chunkloci_lst(chrom, phase_set2hpos_lst)
    else:
        chunkloci_lst = [chunkloci for loci in loci_lst for chunkloci in himut.util.chunkloci(loci)]
  
    somatic_tsbs_lst = []
    filtered_somatic_tsbs_lst = [] 
    alignments = pysam.AlignmentFile(bam_file, "rb")
    for (chrom, chunk_start, chunk_end) in chunkloci_lst:
        if not non_human_sample and common_snps.endswith(".bgz"):
            common_snp_set = himut.vcflib.load_bgz_common_snp((chrom, (chunk_start - qlen_upper_limit), (chunk_end + qlen_upper_limit)), common_snps)
            
        if not non_human_sample and not create_panel_of_normals and panel_of_normals.endswith(".bgz"):
            pon_sbs_set = himut.vcflib.load_bgz_pon((chrom, (chunk_start - qlen_upper_limit), (chunk_end + qlen_upper_limit)), panel_of_normals) 

        somatic_tsbs_candidate_lst = []
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

            ccs.cs2mut()
            ccs_somatic_tsbs_candidate_lst =  [tsbs for tsbs in ccs.tsbs_lst if tsbs not in seen]
            if len(ccs_somatic_tsbs_candidate_lst) == 0:
                continue
            somatic_tsbs_candidate_lst.extend(ccs_somatic_tsbs_candidate_lst)
            
        for tsbs in set(somatic_tsbs_candidate_lst):
            if tsbs in seen:
                continue
           
            seen.add(tsbs)
            tpos, ref, alt = tsbs
            som_gt = "{}{}".format(ref, alt)
            allelecounts = tpos2allelecounts[tpos]
            allele2bq_lst = tpos2allele2bq_lst[tpos]
            bq, vaf, ref_count, alt_count, indel_count, read_depth = himut.bamlib.get_sbs_allelecounts(ref, alt, allelecounts, allele2bq_lst)
            germ_gt, germ_gq, germ_gt_state = himut.gtlib.get_germline_gt(ref, allelecounts, allele2bq_lst, germline_snv_prior)
            if is_germ_gt(som_gt, germ_gt):
                continue 
           
            if germ_gq < min_gq:
                filtered_somatic_tsbs_lst.append((chrom, tpos, ref, alt, "LowGQ", bq, read_depth, ref_count, alt_count, vaf, "."))
                continue

            if germ_gt_state == "het":
                filtered_somatic_tsbs_lst.append((chrom, tpos, ref, alt, "HetSite", bq, read_depth, ref_count, alt_count, vaf, "."))
                continue
            elif germ_gt_state == "hetalt":
                filtered_somatic_tsbs_lst.append((chrom, tpos, ref, alt, "HetAltSite", bq, read_depth, ref_count, alt_count, vaf, "."))
                continue
            elif germ_gt_state == "homref":
                if read_depth == (ref_count + indel_count):
                    continue
            elif germ_gt_state == "homalt":
                if reference_sample: # assembly error
                    continue
                if read_depth == (alt_count + indel_count):
                    continue
               
            if indel_count != 0:
                filtered_somatic_tsbs_lst.append((chrom, tpos, ref, alt, "IndelSite", bq, read_depth, ref_count, alt_count, vaf, "."))
                continue
                
            if not non_human_sample: ## haplotype based estimation of contamiantion prior ## TODO
                if tsbs in common_snp_set:
                    filtered_somatic_tsbs_lst.append((chrom, tpos, ref, alt, "ContRead", bq, read_depth, ref_count, alt_count, vaf, "."))
                    continue 

            alt_bq_lst = allele2bq_lst[himut.util.base2idx[alt]]
            if alt_bq_lst.count(min_bq) == 0:
                filtered_somatic_tsbs_lst.append((chrom, tpos, ref, alt, "LowBQ", bq, read_depth, ref_count, alt_count, vaf, "."))
                continue

            if not non_human_sample and not create_panel_of_normals:
                if tsbs in pon_sbs_set:
                    filtered_somatic_tsbs_lst.append((chrom, tpos, ref, alt, "PanelOfNormal", bq, read_depth, ref_count, alt_count, vaf, "."))
                    continue

            # print(chrom, pos, ref, alt, ref_count, alt_count, read_depth, somatic_pD, germline_pD, somatic_pD/germline_pD)
            # ref_count = allelecounts[himut.util.base2idx[ref]] 
            # alt_count = allelecounts[himut.util.base2idx[alt]]
            # ref_alt_count = ref_count + alt_count
            # nck = get_n_choose_k(alt_count, read_depth)
            # som_likelihood = somatic_snv_prior
            # # som_likelihood = nck*(som_p**ref_count)*(som_q**alt_count)
            # germ_likelihood = nck*(0.5**read_depth)*germline_snv_prior
            # som_pd = som_likelihood/(som_likelihood+germ_likelihood)
            # print((chrom, pos, ref, alt, bq, read_depth, ref_count, alt_count, vaf, som_likelihood, germ_likelihood, som_pd))

            # print(chrom, pos, ref, alt, ref_count, alt_count, indel_count, read_depth)
            if ref_count < min_ref_count and alt_count < min_alt_count:
                filtered_somatic_tsbs_lst.append((chrom, tpos, ref, alt, "LowDepth", bq, read_depth, ref_count, alt_count, vaf, "."))
                continue 
            elif ref_count >= min_ref_count and alt_count < min_alt_count:
                filtered_somatic_tsbs_lst.append((chrom, tpos, ref, alt, "LowDepth", bq, read_depth, ref_count, alt_count, vaf, "."))
                continue 
            elif ref_count < min_ref_count and alt_count >= min_alt_count:
                filtered_somatic_tsbs_lst.append((chrom, tpos, ref, alt, "LowDepth", bq, read_depth, ref_count, alt_count, vaf, "."))
                continue
            if read_depth > md_threshold:
                filtered_somatic_tsbs_lst.append((chrom, tpos, ref, alt, "HighDepth", bq, read_depth, ref_count, alt_count, vaf, ".")) 
                continue
              
            som_state = 0
            alt_ccs_set = set(tpos2allele2ccs_lst[tpos][himut.util.base2idx[alt]])
            for j in alignments.fetch(chrom, tpos, tpos+1):
                qccs = himut.bamlib.BAM(j)
                if qccs.qname not in alt_ccs_set: 
                    continue
                qccs.cs2subindel() 
                qpos = qccs.qsbs_lst[qccs.tsbs_lst.index(tsbs)][0]
                trimmed_qstart = math.floor(min_trim * qccs.qlen)
                trimmed_qend = math.ceil((1 - min_trim) * qccs.qlen)
                if qpos < trimmed_qstart:
                    som_state = 1
                    filtered_somatic_tsbs_lst.append((chrom, tpos, ref, alt, "Trimmed", bq, read_depth, ref_count, alt_count, vaf, "."))
                    break     
                elif qpos > trimmed_qend:
                    som_state = 1
                    filtered_somatic_tsbs_lst.append((chrom, tpos, ref, alt, "Trimmed", bq, read_depth, ref_count, alt_count, vaf, "."))
                    break     
                mpos_lst = [mismatch[0] for mismatch in qccs.mismatch_lst]
                mismatch_start, mismatch_end = himut.util.get_mismatch_range(tpos, qpos, qccs.qlen, mismatch_window)
                mismatch_count = bisect.bisect_right(mpos_lst, mismatch_end) - bisect.bisect_left(mpos_lst, mismatch_start) - 1
                if mismatch_count > max_mismatch_count:
                    som_state = 1
                    filtered_somatic_tsbs_lst.append((chrom, tpos, ref, alt, "MismatchConflict", bq, read_depth, ref_count, alt_count, vaf, "."))
                    break
            if som_state:
                continue 

            if phase:
                som_hap_lst = []
                hap2count = defaultdict(lambda: 0)
                wt_ccs_set = tpos2allele2ccs_lst[tpos][himut.util.base2idx[germ_gt[-1]]]
                phase_set = himut.haplib.get_phase_set(chunk_start, chunk_end, hpos_lst, hpos2phase_set)
                for j in alignments.fetch(chrom, tpos, tpos+1):
                    qccs = himut.bamlib.BAM(j)
                    ccs_hap = himut.haplib.get_ccs_hap(
                        qccs,
                        phase_set2hbit_lst[phase_set],
                        phase_set2hpos_lst[phase_set],
                        phase_set2hetsnp_lst[phase_set],
                    )
                    if qccs.qname in wt_ccs_set:
                        hap2count[ccs_hap] += 1
                    elif qccs.qname in alt_ccs_set:
                        hap2count[ccs_hap] += 1
                        som_hap_lst.append(ccs_hap) 

                som_hap_count = len(set(som_hap_lst))
                if hap2count["0"] >= min_hap_count and hap2count["1"] >= min_hap_count and som_hap_count == 1:
                    somatic_tsbs_lst.append((chrom, tpos, ref, alt, "PASS", bq, read_depth, ref_count, alt_count, vaf, phase_set))
                else:
                    filtered_somatic_tsbs_lst.append((chrom, tpos, ref, alt, "Unphased", bq, read_depth, ref_count, alt_count, vaf, "."))
                    continue                    
            else: 
                somatic_tsbs_lst.append((chrom, tpos, ref, alt, "PASS", bq, read_depth, ref_count, alt_count, vaf, "."))
    chrom2tsbs_lst[chrom] = natsort.natsorted(list(set(somatic_tsbs_lst + filtered_somatic_tsbs_lst)))
    alignments.close()


def call_somatic_substitutions(
    bam_file: str,
    ref_file: str,
    vcf_file: str,
    phased_vcf_file: str,
    region: str,
    region_list: str,
    min_mapq: int,
    min_trim: float,
    min_sequence_identity: float,
    min_hq_base_proportion: float,
    min_alignment_proportion: float,
    common_snps: str,
    panel_of_normals: str,
    min_bq: int,
    min_gq: int,
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
    chrom_lst, chrom2loci_lst = himut.util.load_loci(region, region_list, tname2tsize)
    qlen_mean, qlen_lower_limit, qlen_upper_limit, md_threshold = himut.bamlib.get_thresholds(
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
        germline_snv_prior, germline_indel_prior = himut.vcflib.get_germline_priors(chrom_lst, ref_file, vcf_file, reference_sample)

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
        min_trim,
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
        threads,
        somatic_snv_prior, 
        germline_snv_prior,
        germline_indel_prior,
        phase,
        reference_sample,
        non_human_sample,
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
            tname2tsize[chrom],
            bam_file,
            phased_vcf_file,
            common_snps,
            panel_of_normals,
            chrom2loci_lst[chrom],
            min_mapq,
            min_trim,
            qlen_mean,
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
            reference_sample,
            create_panel_of_normals,
            chrom2tsbs_lst,
        )
        for chrom in chrom_lst
    ]
    p.starmap(
        get_somatic_substitutions, get_somatic_substitutions_arg_lst,
    )
    p.close()
    p.join()
    himut.vcflib.dump_sbs(chrom_lst, chrom2tsbs_lst, phase, vcf_header, out_file)  
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
# 