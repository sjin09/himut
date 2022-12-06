import re
import sys
import math
import numpy as np
from typing import Dict, List, Tuple


def cs2lst(cs_tag):
    cslst = [cs for cs in re.split("(:[0-9]+|\*[a-z][a-z]|[=\+\-][A-Za-z]+)", cs_tag)]
    cslst = [cs.upper() for cs in cslst if cs != ""]
    return cslst


def cs2tuple(read) -> List[Tuple[int, str, str, int, int]]:

    qpos = read.qstart
    read.cstuple_lst = []
    cs_lst = cs2lst(read.cs_tag)
    for cs in cs_lst:
        m = cs[1:]
        mlen = len(m)
        qstart = qpos
        if cs.startswith("="):  # match # --cs=long
            cs = ":{}".format(mlen)
            t = (1, m, m, mlen, mlen)
        elif cs.startswith(":"):  # match # --cs=short
            mlen = int(m)
            qend = qpos + mlen
            m = read.qseq[qstart:qend]
            t = (1, m, m, mlen, mlen)
        elif cs.startswith("*"):  # snp # target and query
            mlen = 1
            ref, alt = list(m)
            t = (2, ref, alt, 1, 1)
        elif cs.startswith("+"):  # insertion # query
            ref = read.qseq[qpos - 1]
            alt = ref + m
            t = (3, ref, alt, 0, mlen)
        elif cs.startswith("-"):  # deletion # target
            alt = read.qseq[qpos - 1]
            ref = alt + m
            t = (4, ref, alt, mlen, 0)
            mlen = 0
        qpos += mlen
        read.cstuple_lst.append(t)


def cs2subindel(read):

    tpos = read.tstart
    qpos = read.qstart
    read.tsbs_lst = []
    read.qsbs_lst = []
    read.qsbs_bq_lst = []
    read.mismatch_lst = []
    for cstuple in read.cstuple_lst:
        state, ref, alt, ref_len, alt_len, = cstuple
        if state == 2:  # snp 
            if ref == "N": 
                continue
            read.qsbs_bq_lst.append(read.bq_int_lst[qpos])
            read.tsbs_lst.append((tpos + 1, ref, alt))
            read.qsbs_lst.append((qpos, alt, ref))
        elif state == 3 or state == 4:  # insertion
            read.mismatch_lst.append((tpos, ref, alt))
        tpos += ref_len 
        qpos += alt_len


def cs2mut(read):

    state = 0
    counter = 0
    read.tsbs_lst = []
    read.tdbs_lst = []
    read.qsbs_lst = []
    read.qdbs_lst = [] 
    read.qsbs_bq_lst = []
    read.qdbs_bq_lst = []
    read.mismatch_lst = []
    tpos = read.tstart
    qpos = read.qstart
    for cstuple in read.cstuple_lst:
        mstate, ref, alt, ref_len, alt_len, = cstuple
        if state == 0 and mstate == 1:  # init # match
            state = 0
        elif state == 0 and mstate != 1:  # init # mismatch
            counter += 1
            ref_lst = [ref]
            alt_lst = [alt]
            if mstate == 2:  # snp
                state = 1
                tstart = tpos
                qstart = qpos
            elif mstate == 3 or mstate == 4:  # insertion # deletion
                state = 2
                tstart = tpos - 1
                qstart = qpos - 1
        elif state != 0 and mstate == 2:  # snp
            state = 1
            counter += 1
            ref_lst.append(ref)
            alt_lst.append(alt)
        elif state != 0 and mstate == 3:  # insertion
            state = 2
            counter += 1
            ref_lst.append(ref)
            alt_lst.append(alt)
        elif state != 0 and mstate == 4:  # deletion
            state = 2
            counter += 1
            ref_lst.append(ref)
            alt_lst.append(alt)
        elif (
            state != 0 and mstate == 1 and ref_len <= 10
        ):  # match # mnp: condition # snp, match, snp
            counter += 1
            state = state
            ref_lst.append(ref)
            alt_lst.append(ref)
        elif state != 0 and mstate == 1 and ref_len > 10:  # match # return
            state = 4
        tpos += ref_len 
        qpos += alt_len

        # return
        if state == 4:
            ref = "".join(ref_lst)
            alt = "".join(alt_lst)
            if len(ref) == len(alt) == 1: # snv
                if ref == "N": 
                    continue
                read.qsbs_bq_lst.append(read.bq_int_lst[qstart])
                read.tsbs_lst.append((tstart + 1, ref, alt))
                read.qsbs_lst.append((qstart, alt, ref))
            else:
                if len(ref) == len(alt) == 2: # dbs
                    if ref.count("N") > 0: 
                        continue
                    read.tdbs_lst.append((tstart + 1, ref, alt))
                    read.qdbs_lst.append((qstart, alt, ref))
                    read.qdbs_bq_lst.append(read.bq_int_lst[qstart:qstart+2])
                read.mismatch_lst.append((tstart + 1, ref, alt))
            state = 0 
            counter = 0 


def cs2mutation(read):

    state = 0
    tsbs_lst = []
    tdbs_lst = []
    qsbs_lst = []
    qdbs_lst = []
    sbs_bq_lst = []
    dbs_bq_lst = []
    read.tmbs_lst = []
    read.tindel_lst = []
    read.tcomplex_lst = []
    tpos = read.tstart
    qpos = read.qstart
    for cstuple in read.cstuple_lst:
        mstate, ref, alt, ref_len, alt_len, = cstuple
        if state == 0 and mstate == 1:  # init # match
            state = 0
            counter = 0
        elif state == 0 and mstate != 1:  # init # mismatch
            counter += 1
            ref_lst = [ref]
            alt_lst = [alt]
            if mstate == 2:  # snp
                state = 1
                tstart = tpos
                qstart = qpos
            elif mstate == 3 or mstate == 4:  # insertion # deletion
                state = 2
                tstart = tpos - 1
                qstart = qpos - 1
        elif state != 0 and mstate == 2:  # snp
            state = 1
            counter += 1
            ref_lst.append(ref)
            alt_lst.append(alt)
        elif state != 0 and mstate == 3:  # insertion
            state = 2
            counter += 1
            ref_lst.append(ref)
            alt_lst.append(alt)
        elif state != 0 and mstate == 4:  # deletion
            state = 2
            counter += 1
            ref_lst.append(ref)
            alt_lst.append(alt)
        elif (
            state != 0 and mstate == 1 and ref_len <= 10
        ):  # match # mnp: condition # snp, match, snp
            counter += 1
            state = state
            ref_lst.append(ref)
            alt_lst.append(ref)
        elif state != 0 and mstate == 1 and ref_len > 11:  # match # return
            state = 4
        tpos += ref_len 
        qpos += alt_len

        # return
        if state == 4:
            ref = "".join(ref_lst)
            alt = "".join(alt_lst)
            ref_len = len(ref)
            alt_len = len(alt)
            counts = len(ref_lst)
            if counts == 1 and ref_len == 1 and alt_len == 1: # snv
                sbs_bq_lst.append(read.bq_int_lst[qstart])
                tsbs_lst.append((tstart + 1, ref, alt))
                qsbs_lst.append((qstart, alt, ref))
            elif counts == 1 and ref_len < alt_len:  # insertion
                read.tindel_lst.append(tstart + 1, ref, alt)
            elif counts == 1 and ref_len > alt_len:  # deletion
                read.tindel_lst.append(tstart + 1, ref, alt)
            elif counts > 1 and ref_len == alt_len:  # mbs
                if ref_len == alt_len == 2:
                    tdbs_lst.append(tstart + 1, ref, alt)
                    qdbs_lst.append(qstart, alt, ref)
                    dbs_bq_lst.append(read.bq_int_lst[qstart:qstart+2])
                else:
                    read.tmbs_lst.append(tstart + 1, ref, alt)
            elif counts > 1 and ref_len != alt_len:  # complex
                read.tcomplex_lst.append(tstart + 1, ref, alt)
            state = 0  
    return tsbs_lst, qsbs_lst, sbs_bq_lst


def cs2tpos2qbase(read) -> Dict[int, Tuple[str, int]]:
    """
    Converts cstuple to 1-coordinate based read allele
    
    Parameters:
        tpos (int): reference alignemnt start position
        qpos (int): query alignment start position
        qbq_lst: list containing base quality scores for CCS bases

    Returns:
        dictionary mapping reference position to tuple containing reference base, alternative base and base quality score
    """

    tpos2qbase = {}
    tpos = read.tstart
    qpos = read.qstart
    for cstuple in read.cstuple_lst:
        state, ref, alt, ref_len, alt_len, = cstuple
        if state == 1:  # match
            for i, alt_base in enumerate(alt):
                tpos2qbase[tpos + i + 1] = (alt_base, read.bq_int_lst[qpos + i])
        elif state == 2:  # sub
            tpos2qbase[tpos + 1] = (alt, read.bq_int_lst[qpos])
        elif state == 3:  # insertion
            pass
        elif state == 4:  # deletion
            for j in range(len(ref[1:])):
                tpos2qbase[tpos + j + 1] = ("-", 0)
        tpos += ref_len
        qpos += alt_len
    return tpos2qbase

