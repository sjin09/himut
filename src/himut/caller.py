import json
import time
import math
import pysam
import bisect
import natsort
import numpy as np
import himut.util
import himut.cslib
import himut.dbslib
import himut.haplib
import himut.mutlib
import himut.bamlib
import himut.vcflib
import multiprocessing as mp
from dataclasses import dataclass
from typing import Dict, List, Tuple
from collections import defaultdict, Counter

def init_lists():
    somatic_tsbs_candidate_lst = []
    somatic_tdbs_candidate_lst = []
    tsbs_annotation = defaultdict(set)
    filtered_somatic_tsbs_candidate_lst = [] 
    return tsbs_annotation, somatic_tsbs_candidate_lst, somatic_tdbs_candidate_lst, filtered_somatic_tsbs_candidate_lst 


def init_allelecounts():
    tpos2allelecounts = defaultdict(lambda: np.zeros(6)) 
    tpos2allele2bq_lst = defaultdict(lambda: {0: [], 1:[], 2:[], 3:[], 4:[], 5:[]})
    tpos2allele2ccs_lst = defaultdict(lambda: {0: [], 1:[], 2:[], 3:[], 4:[], 5:[]})
    return tpos2allele2bq_lst, tpos2allele2ccs_lst, tpos2allelecounts 


def update_allelecounts(
    ccs,
    ccs2tpos2qbase: Dict[str, Dict[int, Tuple[str, int]]],
    tpos2allelecounts: Dict[int, np.ndarray],
    tpos2allele2bq_lst: Dict[int, Dict[int, List[int]]],
    tpos2allele2ccs_lst: Dict[int, Dict[int, List[str]]],
    phase: bool
) -> None:

    tpos2qbase = {}
    tpos = ccs.tstart
    qpos = ccs.qstart
    if phase:
        for cstuple in ccs.cstuple_lst:
            state, ref, alt, ref_len, alt_len, = cstuple
            if state == 1:  # match
                for i, alt_base in enumerate(alt):
                    tpos2qbase[tpos + i + 1] = (alt_base, ccs.bq_int_lst[qpos + i])
                    tpos2allelecounts[tpos + i + 1][himut.util.base2idx[alt_base]] += 1
                    tpos2allele2ccs_lst[tpos + i + 1][himut.util.base2idx[alt_base]].append(ccs.qname)
                    tpos2allele2bq_lst[tpos + i + 1][himut.util.base2idx[alt_base]].append(ccs.bq_int_lst[qpos + i])
            elif state == 2:  # sub
                tpos2qbase[tpos + 1] = (alt, ccs.bq_int_lst[qpos])
                tpos2allelecounts[tpos + 1][himut.util.base2idx[alt]] += 1
                tpos2allele2ccs_lst[tpos + 1][himut.util.base2idx[alt]].append(ccs.qname)
                tpos2allele2bq_lst[tpos + 1][himut.util.base2idx[alt]].append(ccs.bq_int_lst[qpos])
            elif state == 3:  # insertion
                tpos2allelecounts[tpos + 1][4] += 1
                tpos2allele2ccs_lst[tpos + 1][4].append(ccs.qname)
                tpos2allele2bq_lst[tpos + 1][4].append(ccs.bq_int_lst[qpos])
            elif state == 4:  # deletion
                for j in range(len(ref[1:])):
                    tpos2qbase[tpos + j + 1] = ("-", 0)
                    tpos2allelecounts[tpos + j + 1][5] += 1
                    tpos2allele2bq_lst[tpos + j + 1][5].append(0)
                    tpos2allele2ccs_lst[tpos + j + 1][5].append(ccs.qname)
            tpos += ref_len
            qpos += alt_len
    else:
        for cstuple in ccs.cstuple_lst:
            state, ref, alt, ref_len, alt_len, = cstuple
            if state == 1:  # match
                for i, alt_base in enumerate(alt):
                    tpos2qbase[tpos + i + 1] = (alt_base, ccs.bq_int_lst[qpos + i])
                    tpos2allelecounts[tpos + i + 1][himut.util.base2idx[alt_base]] += 1
                    tpos2allele2bq_lst[tpos + i + 1][himut.util.base2idx[alt_base]].append(ccs.bq_int_lst[qpos + i])
            elif state == 2:  # sub
                tpos2qbase[tpos + 1] = (alt, ccs.bq_int_lst[qpos])
                tpos2allelecounts[tpos + 1][himut.util.base2idx[alt]] += 1
                tpos2allele2bq_lst[tpos + 1][himut.util.base2idx[alt]].append(ccs.bq_int_lst[qpos])
            elif state == 3:  # insertion
                tpos2allelecounts[tpos + 1][4] += 1
                tpos2allele2bq_lst[tpos + 1][4].append(ccs.bq_int_lst[qpos])
            elif state == 4:  # deletion
                for j in range(len(ref[1:])):
                    tpos2qbase[tpos + j + 1] = ("-", 0)
                    tpos2allelecounts[tpos + j + 1][5] += 1
                    tpos2allele2bq_lst[tpos + j + 1][5].append(0)
            tpos += ref_len
            qpos += alt_len
    ccs2tpos2qbase[ccs.qname] = tpos2qbase


def get_germline_gt_prior(b1, b2, ref, germline_snv_prior): 
   
    gt_state = 0 
    if b1 == b2 == ref: # homref
        gt_state = "homref"
        gt_prior = 1 - ((1.5 * germline_snv_prior) + (germline_snv_prior * germline_snv_prior)) 
    elif (b1 == ref and b2 != ref) or (b1 != ref and b2 == ref): # het
        gt_state = "het"
        gt_prior = germline_snv_prior
    elif b1 != ref and b2 != ref and b1 != b2: # hetalt
        gt_state = "hetalt"
        gt_prior = germline_snv_prior * germline_snv_prior * 2 ## i don't understand the normalisation
    elif b1 != ref and b2 != ref and b1 == b2: # homalt
        gt_state = "homalt"
        gt_prior = germline_snv_prior/2
    return gt_state, gt_prior


def get_log10_germline_gt_prior(b1, b2, ref, germline_snv_prior): 
   
    gt_state, gt_prior = get_germline_gt_prior(b1, b2, ref, germline_snv_prior)
    gt_prior = math.log10(gt_prior)
    return gt_state, gt_prior


def get_epsilon(
    bq: int
)-> float:
    epsilon = 10**(-bq/10)
    return epsilon


def get_one_minus_epsilon(
    bq: int
) -> float:
    return 1 - get_epsilon(bq)     


def get_one_half_minus_epsilon(
    bq: int
) -> float:
    return 0.5 - get_epsilon(bq)/2.0


def get_log10_epsilon(
    bq: int
) -> float:
    return math.log10(get_epsilon(bq))


def get_log10_one_minus_epsilon(
    bq: int
) -> float:
    return math.log10(1 - get_epsilon(bq))


def get_log10_one_half_minus_epsilon(
    bq: int
) -> float:
    return math.log10(0.5 - get_epsilon(bq)/2.0)


def get_gt_pD(
    b1: str, 
    b2: str,
    allele2bq_lst: Dict[int, List[int]], 
):
    gt_pD = 0
    for base in himut.util.base_lst:
        base_bq_lst = allele2bq_lst[himut.util.base2idx[base]]
        if (b1 == b2) and (base == b1 or base == b2): ## hom
            gt_pD += sum([get_log10_one_minus_epsilon(base_bq) for base_bq in base_bq_lst])
        elif (b1 != b2) and (base == b1 or base == b2): # het
            gt_pD += sum([get_log10_one_half_minus_epsilon(base_bq) for base_bq in base_bq_lst])
        else: ## error
            gt_pD += sum([get_log10_epsilon(base_bq/3) for base_bq in base_bq_lst])
    return gt_pD
    

def get_germline_gt_pD(
    ref: str,
    allelecounts: Dict[int, int],
    allele2bq_lst: Dict[int, List[int]], 
    germline_snv_prior: float,
):

    gt_lst = []
    gt_pl_lst = []
    gt2gt_state = {}
    for gt in himut.util.gt_lst:
        b1, b2 = list(gt)
        b1b2_count = allelecounts[himut.util.base2idx[b1]] + allelecounts[himut.util.base2idx[b2]]
        if b1b2_count == 0:
            continue
        else:
            gt_pD = get_gt_pD(b1, b2, allele2bq_lst)
            gt_state, gt_prior = get_log10_germline_gt_prior(b1, b2, ref, germline_snv_prior) 
            gt_pD += gt_prior
            gt_pl = -10*gt_pD
            gt_lst.append(gt)
            gt_pl_lst.append(gt_pl) 
            gt2gt_state[gt] = gt_state
    return gt_lst, gt_pl_lst, gt2gt_state


def get_germline_gt(
    ref: str,
    allelecounts: Dict[int, int],
    allele2bq_lst: Dict[int, List[int]], 
    germline_snv_prior: float,
):
  
    gt_lst, gt_pl_lst, gt2gt_state = get_germline_gt_pD(
        ref, 
        allelecounts, 
        allele2bq_lst,
        germline_snv_prior, 
    ) 

    j_lst  = np.argsort(gt_pl_lst) 
    gt_lst = [gt_lst[j] for j in j_lst]
    gt_pl_lst = [gt_pl_lst[j] for j in j_lst]
    norm_pl_lst = np.array(gt_pl_lst) - min(gt_pl_lst)
    gt = gt_lst[0]
    gt_state = gt2gt_state[gt] 
    gq = int(norm_pl_lst[1]) if norm_pl_lst[1] < 99  else 99
    return gt, gq, gt_state


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
    vcf_file: str,
    phased_vcf_file: str,
    common_snps: str,
    panel_of_normals: str,
    loci_lst: List[Tuple[str, int, int]],
    ploidy: str,
    min_mapq: int,
    min_trim: float,
    qlen_mean: int,
    qlen_lower_limit: int,
    qlen_upper_limit: int,
    min_sequence_identity: float,
    min_hq_base_proportion: float,
    min_alignment_proportion: float,
    min_bq: int,
    min_gq: int,
    mismatch_window: int,
    max_mismatch_count: int,
    min_ref_count: int,
    min_alt_count: int,
    md_threshold: int,
    min_hap_count: int,
    contamination_prior: float,
    somatic_snv_prior: float,
    germline_snv_prior: float,
    germline_indel_prior: float,
    phase: bool,
    tree_of_life_sample: bool,
    create_panel_of_normals: bool,
    chrom2tsbs_lst: Dict[str, List[List[Tuple[str, int, str, str, int, int, int, int, str]]]],
) -> List[Tuple[str, int, str, str, int, int, int, float, float]]:

    somatic_tsbs_lst = []
    filtered_somatic_tsbs_lst = [] 
    if vcf_file.endswith(".vcf"):
        sample_snp_set = himut.vcflib.load_snp(chrom, vcf_file)
        
    common_snp_set = set()
    if not tree_of_life_sample and common_snps.endswith(".vcf"):
        common_snp_set = himut.vcflib.load_common_snp(chrom, common_snps)
       
    pon_sbs_set = set()
    pon_dbs_set = set() 
    if not create_panel_of_normals and not tree_of_life_sample and panel_of_normals.endswith(".vcf"):
        pon_sbs_set, pon_dbs_set, = himut.vcflib.load_pon(chrom, panel_of_normals)
         
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
 
    seen = set()
    alignments = pysam.AlignmentFile(bam_file, "rb")
    for loci in loci_lst:
        chunkloci_lst = himut.util.chunkloci(loci)
        for (chrom, chunk_start, chunk_end) in chunkloci_lst:
            if vcf_file.endswith(".bgz"):
                sample_snp_set = himut.vcflib.load_bgz_snp((chrom, (chunk_start - qlen_upper_limit), (chunk_end + qlen_upper_limit)), vcf_file)
                
            if not tree_of_life_sample and common_snps.endswith(".bgz"):
                common_snp_set = himut.vcflib.load_bgz_common_snp((chrom, (chunk_start - qlen_upper_limit), (chunk_end + qlen_upper_limit)), common_snps)
                
            if not create_panel_of_normals and not tree_of_life_sample and panel_of_normals.endswith(".bgz"):
                pon_sbs_set, pon_dbs_set = himut.vcflib.load_bgz_pon((chrom, (chunk_start - qlen_upper_limit), (chunk_end + qlen_upper_limit)), panel_of_normals) 

            ccs2tpos2qbase = defaultdict(dict)
            tpos2allele2bq_lst, tpos2allele2ccs_lst, tpos2allelecounts =  init_allelecounts()
            tsbs_annotation, somatic_tsbs_candidate_lst, somatic_tdbs_candidate_lst, filtered_somatic_tsbs_candidate_lst  = init_lists()
            for line in alignments.fetch(chrom, chunk_start, chunk_end):
                ccs = himut.bamlib.BAM(line)
                ccs_somatic_tsbs_candidate_lst = [] 
                ccs_somatic_qsbs_candidate_lst = [] 
                ccs_somatic_tdbs_candidate_lst = [] 
                ccs_somatic_qsbs_candidate_bq_lst = [] 
                for i, tsbs in enumerate(ccs.tsbs_lst):
                    if tsbs in sample_snp_set:
                        continue
                    ccs_somatic_tsbs_candidate_lst.append(tsbs)
                    ccs_somatic_qsbs_candidate_lst.append(ccs.qsbs_lst[i])
                    ccs_somatic_qsbs_candidate_bq_lst.append(ccs.qsbs_bq_lst[i]) 
                
                if not ccs.is_primary:
                    continue

                update_allelecounts(ccs, ccs2tpos2qbase, tpos2allelecounts, tpos2allele2bq_lst, tpos2allele2ccs_lst, phase)
                if ccs.mapq < min_mapq:
                    continue

                if ccs.qlen < qlen_lower_limit or ccs.qlen > qlen_upper_limit:
                    continue

                if ccs.query_alignment_proportion < min_alignment_proportion:
                    continue 

                if himut.bamlib.get_hq_base_proportion(ccs) < min_hq_base_proportion:
                    continue

                if himut.util.get_blast_sequence_identity(ccs) < min_sequence_identity:
                    continue

                common_cnt = 0
                if not tree_of_life_sample: ## probabilistic?
                    for tsbs in ccs_somatic_tsbs_candidate_lst:
                        if tsbs in common_snp_set:
                            common_cnt += 1
                            tsbs_annotation[tsbs].add("CommonVariant")
                            filtered_somatic_tsbs_candidate_lst.append(tsbs)
                    if common_cnt > 0:
                        continue
                    
                if len(ccs_somatic_tsbs_candidate_lst) == 0 and len(ccs_somatic_tdbs_candidate_lst) == 0:
                    continue
               
                trimmed_qstart = math.floor(min_trim * ccs.qlen)
                trimmed_qend = math.ceil((1 - min_trim) * ccs.qlen)
                mismatch_lst = natsort.natsorted(ccs_somatic_tsbs_candidate_lst + ccs.mismatch_lst)
                mismatch_set = set(mismatch_lst)
                mpos_lst = [mismatch[0] for mismatch in mismatch_lst]
                for j, (tsbs, qsbs) in enumerate(zip(ccs_somatic_tsbs_candidate_lst, ccs_somatic_qsbs_candidate_lst)):
                    if tsbs in seen:
                        continue

                    bq = ccs_somatic_qsbs_candidate_bq_lst[j]
                    if bq < min_bq:
                        tsbs_annotation[tsbs].add("LowBQ")  
                        filtered_somatic_tsbs_candidate_lst.append(tsbs) 
                        continue
                    
                    tpos = tsbs[0]
                    qpos = qsbs[0]
                    if qpos < trimmed_qstart:
                        tsbs_annotation[tsbs].add("Trimmed")
                        filtered_somatic_tsbs_candidate_lst.append(tsbs) 
                        continue
                    elif qpos > trimmed_qend:
                        tsbs_annotation[tsbs].add("Trimmed")
                        filtered_somatic_tsbs_candidate_lst.append(tsbs) 
                        continue
                   
                    if not tree_of_life_sample and not create_panel_of_normals:
                        if tsbs in pon_sbs_set:
                            tsbs_annotation[tsbs].add("PanelOfNormal")  
                            filtered_somatic_tsbs_candidate_lst.append(tsbs) 
                            continue

                    
                    mismatch_start, mismatch_end = himut.util.get_mismatch_range(tpos, qpos, ccs.qlen, mismatch_window)
                    idx = bisect.bisect_left(mpos_lst, mismatch_start)
                    jdx = bisect.bisect_right(mpos_lst, mismatch_end)
                    if tsbs in mismatch_set:
                        mismatch_count = (jdx - idx) - 1 
                    else:
                        mismatch_count = jdx - idx 

                    if mismatch_count > max_mismatch_count:
                        tsbs_annotation[tsbs].add("MismatchConflict")
                        filtered_somatic_tsbs_candidate_lst.append(tsbs) 
                        continue
                    somatic_tsbs_candidate_lst.append(tsbs)
                    seen.add(tsbs)

            for tsbs in set(somatic_tsbs_candidate_lst):
                pos, ref, alt = tsbs
                allelecounts = tpos2allelecounts[pos]
                allele2bq_lst = tpos2allele2bq_lst[pos]
                bq, vaf, ref_count, alt_count, indel_count, read_depth = himut.bamlib.get_sbs_allelecounts(ref, alt, allelecounts, allele2bq_lst)
                if indel_count != 0:
                    filtered_somatic_tsbs_lst.append((chrom, pos, ref, alt, "IndelSite", bq, read_depth, ref_count, alt_count, vaf, "."))
                    continue
               
                germline_gt, germline_gq, germline_gt_state = get_germline_gt(ref, allelecounts, allele2bq_lst, germline_snv_prior)
                if germline_gq < min_gq:
                    filtered_somatic_tsbs_lst.append((chrom, pos, ref, alt, "LowGQ", bq, read_depth, ref_count, alt_count, vaf, "."))
                    continue

                if germline_gt_state == "het":
                    filtered_somatic_tsbs_lst.append((chrom, pos, ref, alt, "HetSite", bq, read_depth, ref_count, alt_count, vaf, "."))
                    continue
                elif germline_gt_state == "hetalt":
                    filtered_somatic_tsbs_lst.append((chrom, pos, ref, alt, "HetAltSite", bq, read_depth, ref_count, alt_count, vaf, "."))
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
                
                if read_depth > md_threshold:
                    filtered_somatic_tsbs_lst.append((chrom, pos, ref, alt, "HighDepth", bq, read_depth, ref_count, alt_count, vaf, "."))
                    continue

                if ref_count >= min_ref_count and alt_count >= min_alt_count:
                    if phase: 
                        bidx = himut.haplib.get_bidx(
                            pos - qlen_mean/2, 
                            pos + qlen_mean/2, 
                            hpos_lst, 
                            hetsnp_lst, 
                            hetsnp2bidx
                        )
                        if bidx == ".":
                            filtered_somatic_tsbs_lst.append((chrom, pos, ref, alt, "Unphased", bq, read_depth, ref_count, alt_count, vaf, "."))
                            continue

                        h0_count, h1_count, som_hap, som_hap_count = himut.haplib.get_hap_counts(
                            ref,
                            alt,
                            hblock_lst[bidx], 
                            hidx2hetsnp,
                            ccs2tpos2qbase,
                            tpos2allele2ccs_lst[pos]
                        )
                        if som_hap == ".": 
                            filtered_somatic_tsbs_lst.append((chrom, pos, ref, alt, "Unphased", bq, read_depth, ref_count, alt_count, vaf, "."))
                            continue

                        if h0_count >= min_hap_count and h1_count >= min_hap_count:
                            phase_set = hidx2hetsnp[hblock_lst[bidx][0][0]][0]
                            # print(chrom, pos, ref, alt, h0_count, h1_count, som_hap, som_hap_count) 
                            somatic_tsbs_lst.append((chrom, tpos, ref, alt, "PASS", bq, read_depth, ref_count, alt_count, vaf, phase_set))
                        else:
                            filtered_somatic_tsbs_lst.append((chrom, pos, ref, alt, "Unphased", bq, read_depth, ref_count, alt_count, vaf, "."))
                            continue                    
                    else:
                        somatic_tsbs_lst.append( 
                            (chrom, pos, ref, alt, "PASS", bq, read_depth, ref_count, alt_count, vaf, ".")
                        )
                        # print(chrom, pos, ref, alt, bq, read_depth, ref_count, alt_count, vaf, ".", germline_gt, germline_gq, germline_gt_state)
                else:
                    filtered_somatic_tsbs_lst.append((chrom, pos, ref, alt, "LowDepth", bq, read_depth, ref_count, alt_count, vaf, "."))
                    continue

            for (pos, ref, alt) in set(filtered_somatic_tsbs_candidate_lst): # filtered single base substitutions
                allelecounts = tpos2allelecounts[pos]
                allele2bq_lst = tpos2allele2bq_lst[pos]
                annot = ";".join(natsort.natsorted(list(tsbs_annotation[(pos, ref, alt)])))
                filtered_somatic_tsbs_lst.append((chrom, pos, ref, alt, annot, bq, read_depth, ref_count, alt_count, vaf, "."))

    chrom2tsbs_lst[chrom] = natsort.natsorted(list(set(somatic_tsbs_lst + filtered_somatic_tsbs_lst)))
    alignments.close()


def call_somatic_substitutions(
    bam_file: str,
    vcf_file: str,
    phased_vcf_file: str,
    region: str,
    region_list: str,
    ploidy: str,
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
    contamination_prior: float,
    somatic_snv_prior: float,
    germline_snv_prior: float,
    germline_indel_prior: float,
    threads: int,
    phase: bool,
    non_human_sample: bool,
    create_panel_of_normals: bool,
    version: str,
    out_file: str,
) -> None:

    cpu_start = time.time() / 60
    _, tname2tsize = himut.bamlib.get_tname2tsize(bam_file)
    chrom_lst, chrom2loci_lst = himut.util.load_loci(region, region_list, tname2tsize)
    qlen_std, qlen_mean, qlen_lower_limit, qlen_upper_limit, md_threshold = himut.bamlib.get_thresholds(
        bam_file, chrom_lst, tname2tsize
    )
    himut.util.check_caller_input_exists(
        bam_file,
        vcf_file,
        phased_vcf_file,
        common_snps,
        panel_of_normals,
        chrom_lst,
        tname2tsize,
        ploidy,
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

    # min_ref_count = get_min_ref_count(somatic_mutation_prior, germline_mutation_prior) 
    vcf_header = himut.vcflib.get_himut_vcf_header(
        bam_file,
        vcf_file,
        phased_vcf_file,
        region,
        region_list,
        ploidy,         
        tname2tsize,
        common_snps,
        panel_of_normals,
        min_mapq,
        min_trim,
        min_sequence_identity,
        min_hq_base_proportion,
        min_alignment_proportion,
        min_bq,
        mismatch_window,
        max_mismatch_count,
        min_ref_count,
        min_alt_count,
        md_threshold,
        min_hap_count, 
        threads,
        phase,
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
            min_gq,
            mismatch_window,
            max_mismatch_count,
            min_ref_count,
            min_alt_count,
            md_threshold,
            min_hap_count,
            contamination_prior,
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
        get_somatic_substitutions, get_somatic_substitutions_arg_lst,
    )
    p.close()
    p.join()
    himut.vcflib.dump_himut_sbs(chrom_lst, chrom2tsbs_lst, phase, vcf_header, out_file)  
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