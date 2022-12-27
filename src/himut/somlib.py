import math
import numpy as np
import himut.util
import himut.cslib
import himut.gtlib
import himut.haplib
import himut.mutlib
import himut.bamlib
import himut.vcflib
from typing import Dict, List, Tuple


def init(germline_snv_prior: float):
    global gt_state2gt_prior
    gt_state2gt_prior = {}
    gt_state2gt_prior["het"] = germline_snv_prior
    gt_state2gt_prior["homalt"] = germline_snv_prior / 2
    gt_state2gt_prior["homref"] = 1 - (1.5 * germline_snv_prior)


def get_log2_epsilon(bq: int) -> float:
    return math.log2(himut.gtlib.get_epsilon(bq))


def get_log2_one_minus_epsilon(bq: int) -> float:
    return math.log2(himut.gtlib.get_one_minus_epsilon(bq))


def get_log2_one_half_minus_epsilon(bq: int) -> float:
    return math.log2(himut.gtlib.get_one_half_minus_epsilon(bq))


def get_log10_epsilon(bq: int) -> float:
    return math.log10(himut.gtlib.get_epsilon(bq))


def get_log10_one_minus_epsilon(bq: int) -> float:
    return math.log10(himut.gtlib.get_one_minus_epsilon(bq))


def get_log10_one_half_minus_epsilon(bq: int) -> float:
    return math.log10(himut.gtlib.get_one_half_minus_epsilon(bq))


def get_som_pD(
    ref: str,
    alt: str,
    somatic_snv_prior: float,
    allele2bq_lst: Dict[int, List[int]],
):

    hom_gt_p_lst = []
    het_gt_p_lst = []
    for base in himut.util.base_lst:
        base_bq_lst = allele2bq_lst[himut.util.base2idx[base]]
        if base == ref or base == alt:  ## hom
            hom_gt_p_lst.extend(
                [himut.gtlib.get_one_minus_epsilon(base_bq) for base_bq in base_bq_lst]
            )
            het_gt_p_lst.extend(
                [
                    himut.gtlib.get_one_half_minus_epsilon(base_bq)
                    for base_bq in base_bq_lst
                ]
            )
        else:  ## error
            error_p_lst = [
                himut.gtlib.get_epsilon(base_bq / 3) for base_bq in base_bq_lst
            ]
            hom_gt_p_lst.extend(error_p_lst)
            het_gt_p_lst.extend(error_p_lst)
    het_gt_pD = np.prod(het_gt_p_lst) * gt_state2gt_prior["het"]
    hom_gt_pD = np.prod(hom_gt_p_lst) * gt_state2gt_prior["homref"]
    som_pD = somatic_snv_prior * (hom_gt_pD + het_gt_pD)
    return som_pD


def get_log2_som_pD(
    ref: str,
    alt: str,
    somatic_snv_prior: float,
    allele2bq_lst: Dict[int, List[int]],
):

    hom_gt_pD = 0
    het_gt_pD = 0
    for base in himut.util.base_lst:
        base_bq_lst = allele2bq_lst[himut.util.base2idx[base]]
        if base == ref or base == alt:  ## hom
            hom_gt_pD += sum(
                [get_log2_one_minus_epsilon(base_bq) for base_bq in base_bq_lst]
            )
            het_gt_pD += sum(
                [get_log2_one_half_minus_epsilon(base_bq) for base_bq in base_bq_lst]
            )
        else:  ## error
            error_pD = sum([get_log2_epsilon(base_bq / 3) for base_bq in base_bq_lst])
            hom_gt_pD += error_pD
            het_gt_pD += error_pD
    het_gt_pD += math.log2(gt_state2gt_prior["het"])
    hom_gt_pD += math.log2(gt_state2gt_prior["homref"])
    som_pD = math.log2(somatic_snv_prior) + hom_gt_pD + het_gt_pD
    return som_pD


def get_not_som_pD(
    ref: str,
    alt: str,
    allele2bq_lst: Dict[int, List[int]],
) -> float:

    hom_gt_p_lst = []
    het_gt_p_lst = []
    for base in himut.util.base_lst:
        if alt == base:
            continue
        base_bq_lst = allele2bq_lst[himut.util.base2idx[base]]
        if base == ref or base == alt:  ## hom
            hom_gt_p_lst.extend(
                [himut.gtlib.get_one_minus_epsilon(base_bq) for base_bq in base_bq_lst]
            )
            het_gt_p_lst.extend(
                [
                    himut.gtlib.get_one_half_minus_epsilon(base_bq)
                    for base_bq in base_bq_lst
                ]
            )
        else:  ## error
            error_p_lst = [
                himut.gtlib.get_epsilon(base_bq / 3) for base_bq in base_bq_lst
            ]
            hom_gt_p_lst.extend(error_p_lst)
            het_gt_p_lst.extend(error_p_lst)

    het_gt_pD = np.prod(het_gt_p_lst) * gt_state2gt_prior["het"]
    hom_gt_pD = np.prod(hom_gt_p_lst) * gt_state2gt_prior["homref"]
    not_som_pD = het_gt_pD + hom_gt_pD
    return not_som_pD


def get_log2_not_som_pD(
    ref: str,
    alt: str,
    allele2bq_lst: Dict[int, List[int]],
) -> float:

    hom_gt_pD = 0
    het_gt_pD = 0
    for base in himut.util.base_lst:
        if alt == base:
            continue
        base_bq_lst = allele2bq_lst[himut.util.base2idx[base]]
        if base == ref or base == alt:  ## hom
            hom_gt_pD += sum(
                [get_log2_one_minus_epsilon(base_bq) for base_bq in base_bq_lst]
            )
            het_gt_pD += sum(
                [get_log2_one_half_minus_epsilon(base_bq) for base_bq in base_bq_lst]
            )
        else:  ## error
            error_pD = sum([get_log2_epsilon(base_bq / 3) for base_bq in base_bq_lst])
            hom_gt_pD += error_pD
            het_gt_pD += error_pD

    het_gt_pD += math.log2(gt_state2gt_prior["het"])
    hom_gt_pD += math.log2(gt_state2gt_prior["homref"])
    not_som_pD = hom_gt_pD + het_gt_pD
    return not_som_pD
