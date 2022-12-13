import math
import numpy as np
import himut.util
import himut.cslib
import himut.haplib
import himut.mutlib
import himut.bamlib
import himut.vcflib
from typing import Dict, List, Tuple


def get_germline_gt_prior(
    b1: str, 
    b2: str, 
    ref: str, 
    germline_snv_prior: float
) -> Tuple[str, float]: 
   
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


def get_log10_germline_gt_prior(
    b1: str, 
    b2: str, 
    ref: str, 
    germline_snv_prior: float
) -> Tuple[str, float]: 
   
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
) -> float:
    
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
) -> Tuple[List[str], List[float], str]:

    gt_lst = []
    gt_pl_lst = []
    gt2gt_state = {}
    for gt in himut.util.gt_lst:
        b1, b2 = list(gt)
        b1b2_count = allelecounts[himut.util.base2idx[b1]] + allelecounts[himut.util.base2idx[b2]]
        if b1b2_count == 0:
            continue
        else:
            # print(gt, b1, b2, b1b2_count)
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
    allelecounts: np.ndarray,
    allele2bq_lst: Dict[int, List[int]], 
    germline_snv_prior: float,
) -> Tuple[List[str], List[float], str]:
  
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
    if gt[0] != ref and gt.count(ref) == 1:
        gt = gt[::-1] 
    return gt, gq, gt_state
