import natsort
from collections import defaultdict

purine = set(["A", "G"])
pyrimidine = set(["T", "C"])
BASE_COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}


sbs96_to_sbs52 = {
    "A[C>A]A":"A[C>A]A",
    "A[C>A]C":"A[C>A]C",
    "A[C>A]G":"A[C>A]G",
    "A[C>A]T":"A[C>A]T",
    "C[C>A]A":"C[C>A]A",
    "C[C>A]C":"C[C>A]C",
    "C[C>A]G":"C[C>A]G",
    "C[C>A]T":"C[C>A]T",
    "G[C>A]A":"G[C>A]A",
    "G[C>A]C":"G[C>A]C",
    "G[C>A]G":"G[C>A]G",
    "G[C>A]T":"G[C>A]T",
    "T[C>A]A":"T[C>A]A",
    "T[C>A]C":"T[C>A]C",
    "T[C>A]G":"T[C>A]G",
    "T[C>A]T":"T[C>A]T",
    "T[T>G]T":"A[C>A]A",
    "G[T>G]T":"A[C>A]C",
    "C[T>G]T":"A[C>A]G",
    "A[T>G]T":"A[C>A]T",
    "T[T>G]G":"C[C>A]A",
    "G[T>G]G":"C[C>A]C",
    "C[T>G]G":"C[C>A]G",
    "A[T>G]G":"C[C>A]T",
    "T[T>G]C":"G[C>A]A",
    "G[T>G]C":"G[C>A]C",
    "C[T>G]C":"G[C>A]G",
    "A[T>G]C":"G[C>A]T",
    "T[T>G]A":"T[C>A]A",
    "G[T>G]A":"T[C>A]C",
    "C[T>G]A":"T[C>A]G",
    "A[T>G]A":"T[C>A]T",
    "A[C>T]A":"A[C>T]A",
    "A[C>T]C":"A[C>T]C",
    "A[C>T]G":"A[C>T]G",
    "A[C>T]T":"A[C>T]T",
    "C[C>T]A":"C[C>T]A",
    "C[C>T]C":"C[C>T]C",
    "C[C>T]G":"C[C>T]G",
    "C[C>T]T":"C[C>T]T",
    "G[C>T]A":"G[C>T]A",
    "G[C>T]C":"G[C>T]C",
    "G[C>T]G":"G[C>T]G",
    "G[C>T]T":"G[C>T]T",
    "T[C>T]A":"T[C>T]A",
    "T[C>T]C":"T[C>T]C",
    "T[C>T]G":"T[C>T]G",
    "T[C>T]T":"T[C>T]T",
    "A[T>C]A":"A[C>T]A",
    "A[T>C]C":"A[C>T]C",
    "A[T>C]G":"A[C>T]G",
    "A[T>C]T":"A[C>T]T",
    "C[T>C]A":"C[C>T]A",
    "C[T>C]C":"C[C>T]C",
    "C[T>C]G":"C[C>T]G",
    "C[T>C]T":"C[C>T]T",
    "G[T>C]A":"G[C>T]A",
    "G[T>C]C":"G[C>T]C",
    "G[T>C]G":"G[C>T]G",
    "G[T>C]T":"G[C>T]T",
    "T[T>C]A":"T[C>T]A",
    "T[T>C]C":"T[C>T]C",
    "T[T>C]G":"T[C>T]G",
    "T[T>C]T":"T[C>T]T",
    "A[C>G]A":"A[C>G]A",
    "A[C>G]C":"A[C>G]C",
    "A[C>G]G":"A[C>G]G",
    "A[C>G]T":"A[C>G]T",
    "C[C>G]A":"C[C>G]A",
    "C[C>G]C":"C[C>G]C",
    "C[C>G]G":"C[C>G]G",
    "C[C>G]T":"A[C>G]G",
    "G[C>G]A":"G[C>G]A",
    "G[C>G]C":"G[C>G]C",
    "G[C>G]G":"C[C>G]C",
    "G[C>G]T":"A[C>G]C",
    "T[C>G]A":"T[C>G]A",
    "T[C>G]C":"G[C>G]A",
    "T[C>G]G":"C[C>G]A",
    "T[C>G]T":"A[C>G]A",
    "A[T>A]A":"A[T>A]A",
    "A[T>A]C":"A[T>A]C",
    "A[T>A]G":"A[T>A]G",
    "A[T>A]T":"A[T>A]T",
    "C[T>A]A":"C[T>A]A",
    "C[T>A]C":"C[T>A]C",
    "C[T>A]G":"C[T>A]G",
    "C[T>A]T":"A[T>A]G",
    "G[T>A]A":"G[T>A]A",
    "G[T>A]C":"G[T>A]C",
    "G[T>A]G":"C[T>A]C",
    "G[T>A]T":"A[T>A]C",
    "T[T>A]A":"T[T>A]A",
    "T[T>A]C":"G[T>A]A",
    "T[T>A]G":"C[T>A]A",
    "T[T>A]T":"A[T>A]A"
}

# def get_complementary_sbs96(sbs96):

#     sbs96_complement = ""
#     for s in sbs96:
#         if s.isalpha():
#             sbs96_complement += BASE_COMPLEMENT[s]
#         else:
#             sbs96_complement += s
#     return sbs96_complement


def get_reverse_complementary(tri):
   
    tri_rc = "".join([BASE_COMPLEMENT[nt] for nt in tri[::-1]])
    return tri_rc
   
    
sbs52_to_sbs96_lst = defaultdict(list) 
for (sbs96, sbs52) in sbs96_to_sbs52.items():
    sbs52_to_sbs96_lst[sbs52].append(sbs96)

sbs96_to_sbs52 = {} 
sbs96_to_pyr_count = {}
for (sbs52, sbs96_lst) in sbs52_to_sbs96_lst.items():
    for sbs96 in sbs96_lst:
        counter = 0
        ubase, _, ref, _, alt, _, dbase = list(sbs96)
        tri = f"{ubase}{ref}{dbase}"
        for nt in tri: 
            if nt in pyrimidine:
                counter += 1
        sbs96_to_pyr_count[sbs96] = counter

updated_sbs96_to_sbs52 = {}
updated_sbs52_to_sbs96_lst = defaultdict(list)
for (sbs52, sbs96_lst) in sbs52_to_sbs96_lst.items():
    sbs96_pyr_count_lst = [sbs96_to_pyr_count[sbs96] for sbs96 in sbs96_lst]
    sbs96_pyr_max_count = max(sbs96_pyr_count_lst)
    idx_lst = [i for (i, pyr_count) in enumerate(sbs96_pyr_count_lst) if pyr_count == sbs96_pyr_max_count]
    if len(idx_lst) == 1:
        new_sbs52 = sbs96_lst[idx_lst[0]]
        for sbs96 in sbs96_lst:
            updated_sbs96_to_sbs52[sbs96] = new_sbs52 
    else:
        sorted_sbs96_lst = natsort.natsorted(sbs96_lst)
        if sbs96_lst != sorted_sbs96_lst:
            new_sbs52 = sorted_sbs96_lst[0] 
        new_sbs52 = sorted_sbs96_lst[0] 
        for sbs96 in sbs96_lst:
            updated_sbs96_to_sbs52[sbs96] = new_sbs52 
    
    for sbs96 in sbs96_lst:
        ubase, _, ref, _, alt, _, dbase = list(sbs96)
        fwd_ref = f"{ubase}{ref}{dbase}"
        fwd_alt = f"{ubase}{alt}{dbase}" 
        rev_ref = get_reverse_complementary(fwd_ref) 
        rev_alt = get_reverse_complementary(fwd_alt)
        if ref in pyrimidine:
            sbs96_pyr = f"{fwd_ref}>{fwd_alt}"
            sbs96_pur = f"{rev_alt}>{rev_ref}"
        else:
            sbs96_pur = f"{fwd_ref}>{fwd_alt}"
            sbs96_pyr = f"{rev_alt}>{rev_ref}"
        updated_sbs52_to_sbs96_lst[new_sbs52].append(sbs96_pur)
        updated_sbs52_to_sbs96_lst[new_sbs52].append(sbs96_pyr) 
  
for sbs96, sbs52 in updated_sbs96_to_sbs52.items():
    print('"{}": "{}",'.format(sbs96, sbs52))  
        