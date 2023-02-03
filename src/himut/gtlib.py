import math
import numpy as np
import himut.util
import himut.cslib
import himut.haplib
import himut.bamlib
import himut.vcflib
from typing import Dict, List, Tuple
gt_lst = ["AA", "TA", "CA", "GA", "TT", "CT", "GT", "CC", "GC", "GG"]


def init(germline_snv_prior: float):
    global gt_state2gt_prior
    gt_state2gt_prior = {}
    gt_state2gt_prior["het"] = germline_snv_prior
    gt_state2gt_prior["hetalt"] = germline_snv_prior * germline_snv_prior * 2
    gt_state2gt_prior["homref"] = 1 - (
        (1.5 * germline_snv_prior) + (germline_snv_prior * germline_snv_prior)
    )
    gt_state2gt_prior["homalt"] = germline_snv_prior / 2


def get_germ_gt_state(
    b1: str,
    b2: str,
    ref: str,
) -> Tuple[str, float]:

    gt_state = 0
    if b1 == b2 == ref:  # homref
        gt_state = "homref"
    elif (b1 == ref and b2 != ref) or (b1 != ref and b2 == ref):  # het
        gt_state = "het"
    elif b1 != ref and b2 != ref and b1 != b2:  # hetalt
        gt_state = "hetalt"
    elif b1 != ref and b2 != ref and b1 == b2:  # homalt
        gt_state = "homalt"
    return gt_state


def get_log10_germ_gt_prior(gt_state: str) -> Tuple[str, float]:

    gt_prior = math.log10(gt_state2gt_prior[gt_state])
    return gt_prior


def get_epsilon(bq: int) -> float:
    epsilon = 10 ** (-bq / 10)
    return epsilon


def get_one_minus_epsilon(bq: int) -> float:
    return 1 - get_epsilon(bq)


def get_one_half_minus_epsilon(bq: int) -> float:
    return 0.5 - get_epsilon(bq) / 2.0


def get_log10_epsilon(bq: int) -> float:
    return math.log10(get_epsilon(bq))


def get_log10_one_minus_epsilon(bq: int) -> float:
    return math.log10(get_one_minus_epsilon(bq))


def get_log10_one_half_minus_epsilon(bq: int) -> float:
    return math.log10(get_one_half_minus_epsilon(bq))


def get_log10_gt_pD(
    gt: str,
    ref: str,
    allele2bq_lst: Dict[int, List[int]],
) -> float:

    gt_pD = 0
    b1, b2 = list(gt)
    gt_state = get_germ_gt_state(b1, b2, ref)
    gt_prior = get_log10_germ_gt_prior(gt_state)
    for base in himut.util.base_lst:
        base_bq_lst = allele2bq_lst[himut.util.base2idx[base]]
        if (b1 == b2) and (base == b1 or base == b2):  ## hom
            gt_pD += sum(
                [get_log10_one_minus_epsilon(base_bq) for base_bq in base_bq_lst]
            )
        elif (b1 != b2) and (base == b1 or base == b2):  # het
            gt_pD += sum(
                [get_log10_one_half_minus_epsilon(base_bq) for base_bq in base_bq_lst]
            )
        else:  ## error
            gt_pD += sum([get_log10_epsilon(base_bq / 3) for base_bq in base_bq_lst])
    gt_pD += gt_prior
    gt_pl = -10 * gt_pD
    return gt_pl, gt_state


def get_germ_gt_pD(
    ref: str,
    allele2bq_lst: Dict[int, List[int]],
) -> Tuple[List[str], List[float], str]:

    gt_pl_lst = []
    gt2gt_state = {}
    for gt in gt_lst:
        gt_pl, gt_state = get_log10_gt_pD(gt, ref, allele2bq_lst)
        gt_pl_lst.append(gt_pl)
        gt2gt_state[gt] = gt_state
    return gt_pl_lst, gt2gt_state


def get_argmin_gt(gt_pl_lst):
    ilst = np.argsort(gt_pl_lst)
    gt = [gt_lst[i] for i in ilst][0]
    gt_pl_lst = [gt_pl_lst[i] for i in ilst]
    gq = (np.array(gt_pl_lst) - min(gt_pl_lst))[1]
    gq = int(gq) if gq < 99 else 99
    return gt, gq


def get_germ_gt(
    ref: str,
    allele2bq_lst: Dict[int, List[int]],
) -> Tuple[List[str], List[float], str]:

    gt_pl_lst, gt2gt_state = get_germ_gt_pD(
        ref,
        allele2bq_lst,
    )
    gt, gq = get_argmin_gt(gt_pl_lst) 
    gt_state = gt2gt_state[gt]
    if gt[0] != ref and gt.count(ref) == 1:
        gt = gt[::-1]
    return gt, gq, gt_state, gt2gt_state 


def get_germ_gq(
    alt: str,
    gt2gt_state: Dict[str, str],
    allele2bq_lst: Dict[int, List[int]]
):
    germ_gt_lst = []
    germ_gt_pl_lst = []
    for gt in gt_lst:
        gt_pl = 0
        b1, b2 = list(gt)
        gt_state = gt2gt_state[gt] 
        gt_prior = get_log10_germ_gt_prior(gt_state)
        for base in himut.util.base_lst:
            if alt == base:
                continue
            base_bq_lst = allele2bq_lst[himut.util.base2idx[base]]
            if (b1 == b2) and (base == b1 or base == b2):  ## hom
                gt_pl += sum(
                    [get_log10_one_minus_epsilon(base_bq) for base_bq in base_bq_lst]
                )
            elif (b1 != b2) and (base == b1 or base == b2):  # het
                gt_pl += sum(
                    [get_log10_one_half_minus_epsilon(base_bq) for base_bq in base_bq_lst]
                )
            else:  ## error
                gt_pl += sum([get_log10_epsilon(base_bq / 3) for base_bq in base_bq_lst])
        gt_pl += gt_prior
        gt_pl *= -10 
        germ_gt_lst.append(gt)
        germ_gt_pl_lst.append(gt_pl)

    ilst = np.argsort(germ_gt_pl_lst)
    germ_gt_lst = [germ_gt_lst[i] for i in ilst]
    germ_gt_pl_lst = [germ_gt_pl_lst[i] for i in ilst]
    gq = (np.array(germ_gt_pl_lst) - min(germ_gt_pl_lst))[1]
    gq = int(gq) if gq < 99 else 99
    return gq 

