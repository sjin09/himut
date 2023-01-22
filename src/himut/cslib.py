import re
import himut.util
import numpy as np
from typing import Dict, List, Tuple


def cs2lst(cs_tag):
    cslst = [cs for cs in re.split("(:[0-9]+|\*[a-z][a-z]|[=\+\-][A-Za-z]+)", cs_tag)]
    cslst = [cs.upper() for cs in cslst if cs != ""]
    return cslst


def cs2tuple(ccs, cs_tag) -> List[Tuple[int, str, str, int, int]]:

    qpos = ccs.qstart
    ccs.cstuple_lst = []
    cs_lst = cs2lst(cs_tag)
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
            m = ccs.qseq[qstart:qend]
            t = (1, m, m, mlen, mlen)
        elif cs.startswith("*"):  # snp # target and query
            mlen = 1
            ref, alt = list(m)
            t = (2, ref, alt, 1, 1)
        elif cs.startswith("+"):  # insertion # query
            ref = ccs.qseq[qpos - 1]
            alt = ref + m
            t = (3, ref, alt, 0, mlen)
        elif cs.startswith("-"):  # deletion # target
            alt = ccs.qseq[qpos - 1]
            ref = alt + m
            t = (4, ref, alt, mlen, 0)
            mlen = 0
        qpos += mlen
        ccs.cstuple_lst.append(t)


def cs2subindel(ccs):

    tpos = ccs.tstart # 0-coordinate
    qpos = ccs.qstart
    ccs.tsbs_lst = []
    ccs.qsbs_lst = []
    ccs.qsbs_bq_lst = []
    ccs.mismatch_lst = []
    for (state, ref, alt, ref_len, alt_len) in ccs.cstuple_lst:
        if ref != "N" and state == 2:  # snp
            ccs.qsbs_bq_lst.append(ccs.bq_int_lst[qpos])
            ccs.tsbs_lst.append((tpos + 1, ref, alt)) # 1-coordinate
            ccs.qsbs_lst.append((qpos, alt, ref))
            ccs.mismatch_lst.append((tpos+1, ref, alt))
        elif state == 3 or state == 4:  # insertion or deletion
            ccs.mismatch_lst.append((tpos+1, ref, alt))
        tpos += ref_len
        qpos += alt_len


def cs2mut(ccs):

    state = 0
    counter = 0
    ccs.tsbs_lst = []
    ccs.tdbs_lst = []
    ccs.qsbs_lst = []
    ccs.qdbs_lst = []
    ccs.qsbs_bq_lst = []
    ccs.qdbs_bq_lst = []
    ccs.tmbs_lst = []
    ccs.tindel_lst = []
    ccs.tcomplex_lst = []
    ccs.mismatch_lst = []
    tpos = ccs.tstart # 0-coordinate
    qpos = ccs.qstart
    for (mstate, ref, alt, ref_len, alt_len) in ccs.cstuple_lst:
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
            mcount = len(ref_lst)
            ref = "".join(ref_lst)
            alt = "".join(alt_lst)
            if ref != "N" and (len(ref) == len(alt) == mcount == 1):  # snv
                ccs.qsbs_bq_lst.append(ccs.bq_int_lst[qstart])
                ccs.tsbs_lst.append((tstart + 1, ref, alt)) # 1-coordinate
                ccs.qsbs_lst.append((qstart, alt, ref))
            elif len(ref) < len(alt) and mcount == 1:  # insertion
                ccs.tindel_lst.append((tstart + 1, ref, alt))
            elif len(ref) > len(alt) and mcount == 1:  # deletion
                ccs.tindel_lst.append((tstart + 1, ref, alt))
            elif len(ref) == len(alt) and mcount > 1:  # mbs
                if ref_len == alt_len == 2:
                    ccs.tdbs_lst.append((tstart + 1, ref, alt))
                    ccs.qdbs_lst.append((qstart, alt, ref))
                    ccs.dbs_bq_lst.append(ccs.bq_int_lst[qstart : qstart + 2])
                else:
                    ccs.tmbs_lst.append((tstart + 1, ref, alt))
            elif ref_len != alt_len and mcount > 1:  # complex
                ccs.tcomplex_lst.append((tstart + 1, ref, alt))
            ccs.mismatch_lst.append((tstart + 1, ref, alt))
            counter = 0
            state = 0


def cs2tpos2qbase(ccs):

    tpos = ccs.tstart
    qpos = ccs.qstart
    ccs.tpos2qbase = {}
    for (state, ref, alt, ref_len, alt_len) in ccs.cstuple_lst:
        if state == 1:  # match
            for i, alt_base in enumerate(alt):
                ccs.tpos2qbase[tpos + i + 1] = (alt_base, ccs.bq_int_lst[qpos + i])
        elif state == 2:  # sub
            ccs.tpos2qbase[tpos + 1] = (alt, ccs.bq_int_lst[qpos])
        elif state == 3:  # insertion
            pass
        elif state == 4:  # deletion
            for j in range(len(ref[1:])):
                ccs.tpos2qbase[tpos + j + 1] = ("-", 0)
        tpos += ref_len
        qpos += alt_len


