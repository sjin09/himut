#!/usr/bin/env python3

import re
import sys
import pysam
import tabix
import cyvcf2
import bisect
import natsort
import argparse
import multiprocessing as mp
from collections import defaultdict
from typing import Dict, List, Tuple, Set
bit_complement_hsh = {"0": "1", "1": "0", "-": "-"}

def parse_args(args):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-i",
        "--bam",
        type=str,
        required=True,
        help="BAM file to read",
    )
    parser.add_argument(
        "--sbs",
        type=str,
        required=False,
        help="VCF file to read somatic single-base-substitutions",
    )
    parser.add_argument(
        "--phased_vcf",
        type=str,
        required=False,
        help="VCF file to read phased heterozygous SNPs",
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=1,
        required=False,
        help="number of threads to use",
    )
    parser.add_argument(
        "-o",
        "--out",
        type=str,
        required=True,
        help="file to write",
    )
    args = args[1:]
    return parser.parse_args(args)


class BAM:
    def __init__(self, line):
        # target
        if line.is_secondary:
            self.is_primary = False
        else: 
            self.is_primary = True
            self.tname = line.reference_name
            self.tstart = line.reference_start
            self.tend = line.reference_end
            # query
            self.qname = line.query_name
            self.qstart = line.query_alignment_start
            self.qend = line.query_alignment_end
            self.qseq = line.query_sequence
            self.qlen = len(self.qseq)
            self.mapq = line.mapping_quality
            self.bq_int_lst = line.query_qualities
            self.cs2tuple(line.get_tag("cs"))
            
    def cs2lst(self, cs_tag):
        cslst = [cs for cs in re.split("(:[0-9]+|\*[a-z][a-z]|[=\+\-][A-Za-z]+)", cs_tag)]
        cslst = [cs.upper() for cs in cslst if cs != ""]
        return cslst

    def cs2tuple(self, cs_tag) -> List[Tuple[int, str, str, int, int]]:
        qpos = self.qstart
        self.cstuple_lst = []
        cs_lst = self.cs2lst(cs_tag)
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
                m = self.qseq[qstart:qend]
                t = (1, m, m, mlen, mlen)
            elif cs.startswith("*"):  # snp # target and query
                mlen = 1
                ref, alt = list(m)
                t = (2, ref, alt, 1, 1)
            elif cs.startswith("+"):  # insertion # query
                ref = self.qseq[qpos - 1]
                alt = ref + m
                t = (3, ref, alt, 0, mlen)
            elif cs.startswith("-"):  # deletion # target
                alt = self.qseq[qpos - 1]
                ref = alt + m
                t = (4, ref, alt, mlen, 0)
                mlen = 0
            self.cstuple_lst.append(t)
            qpos += mlen

    def cs2subindel(ccs):

        tpos = ccs.tstart
        qpos = ccs.qstart
        ccs.tsbs_lst = []
        for cstuple in ccs.cstuple_lst:
            state, ref, alt, ref_len, alt_len, = cstuple
            if state == 2 and ref != "N":  # snp 
                ccs.tsbs_lst.append((tpos + 1, ref, alt))
            tpos += ref_len 
            qpos += alt_len

    def cs2tpos2qbase(self):

        tpos = self.tstart
        qpos = self.qstart
        self.tpos2qbase = {}
        for cstuple in self.cstuple_lst:
            state, ref, alt, ref_len, alt_len, = cstuple
            if state == 1:  # match
                for i, alt_base in enumerate(alt):
                    self.tpos2qbase[tpos + i + 1] = (alt_base, self.bq_int_lst[qpos + i])
            elif state == 2:  # sub
                self.tpos2qbase[tpos + 1] = (alt, self.bq_int_lst[qpos])
            elif state == 3:  # insertion
                pass
            elif state == 4:  # deletion
                for j in range(len(ref[1:])):
                    self.tpos2qbase[tpos + j + 1] = ("-", 0)
            tpos += ref_len
            qpos += alt_len


class VCF:
    def __init__(self, line):
        arr = line.strip().split()
        self.chrom = arr[0]
        self.pos = int(arr[1])
        self.id = arr[2]
        self.ref = arr[3]
        self.alt_lst = arr[4].split(",")
        self.qual = arr[5]
        self.qual = float(self.qual) if self.qual != "." else self.qual
        self.is_pass = True if arr[6] == "PASS" else False
        self.info = arr[7]
        self.format_lst = arr[8].split(":")
        self.sample_format_lst = arr[9].split(":")
        hsh = {i:j for i,j in zip(self.format_lst, self.sample_format_lst)}
        if "GT" in hsh: self.sample_gt = hsh["GT"] 
        if "PS" in hsh: self.sample_phase_set = hsh["PS"] 
        if "AD" in hsh: 
            arr = hsh["AD"].split(",")
            self.ref_count = arr[0]
            self.alt_count_arr = arr[1:]

        self.is_snp = False
        self.is_dbs = False
        self.is_indel = False
        if len(self.alt_lst) == 1:  # bi-allelic
            self.is_biallelic = True
            self.alt = self.alt_lst[0]
            if len(self.ref) == 1 and len(self.alt) == 1:  # snp
                self.is_snp = True
            elif len(self.ref) == len(self.alt) == 2:
                self.is_dbs = True
            elif len(self.ref) > len(self.alt): # del
                self.is_indel = True
            elif len(self.ref) < len(self.alt): # ins
                self.is_indel = True
        else:
            self.is_biallelic = False


def get_tname2tsize(bam_file: str) -> Tuple[List[str], Dict[str, int]]:

    tname2tsize = {}
    alignments = pysam.AlignmentFile(bam_file, "rb")
    bam_header_lst = str(alignments.header).strip().split("\n")
    for h in bam_header_lst:
        if h.startswith("@SQ"):
            arr = h.split("\t")
            tname = arr[1].replace("SN:", "")
            tsize = int(arr[2].replace("LN:", ""))
            tname2tsize[tname] = tsize
    alignments.close()
    return tname2tsize


def load_sbs(vcf_file: str) -> Dict[str, Tuple[int, str, str]]:

    chrom2sbs2ps = defaultdict(dict) 
    chrom2sbs_lst = defaultdict(list)
    if vcf_file.endswith(".vcf"):
        for line in open(vcf_file).readlines():
            if line.startswith("#"):
                continue
            v = VCF(line)
            if v.is_pass:
                sbs = (v.pos, v.ref, v.alt)
                chrom2sbs2ps[v.chrom][sbs] = v.sample_phase_set
                chrom2sbs_lst[v.chrom].append((v.pos, v.ref, v.alt))
    elif vcf_file.endswith(".bgz"):
        for i in cyvcf2.VCF(vcf_file):
            v = VCF(str(i))
            if v.is_pass:
                sbs = (v.pos, v.ref, v.alt)
                chrom2sbs2ps[v.chrom][sbs] = v.sample_phase_set
                chrom2sbs_lst[v.chrom].append((v.pos, v.ref, v.alt))
    chrom_lst = natsort.natsorted(chrom2sbs_lst.keys())
    return chrom_lst, chrom2sbs2ps, chrom2sbs_lst


def load_phase_set(
    vcf_file: str, 
    chrom: str,
    chrom_len: int,
) -> Dict[str, List[List[Tuple[int, str]]]]:

    phase_set2h0_lst = defaultdict(list) 
    phase_set2hpos_lst = defaultdict(list) 
    phase_set2hetsnp_lst = defaultdict(list) 
    if vcf_file.endswith(".vcf"):
        for line in open(vcf_file):
            if line.startswith("#"):
                continue
            v = VCF(line)
            if chrom == v.chrom: 
                if v.sample_gt == "0|1" or v.sample_gt == "1|0":
                    hstate = v.sample_gt.split("|")[0]                 
                    phase_set2h0_lst[v.sample_phase_set].append(hstate)
                    phase_set2hpos_lst[v.sample_phase_set].append(v.pos)
                    phase_set2hetsnp_lst[v.sample_phase_set].append((v.pos, v.ref, v.alt))
    elif vcf_file.endswith(".bgz"):
        tb = tabix.open(vcf_file)
        records = tb.query(chrom, 0, chrom_len)
        for record in records:
            v = VCF("\t".join(record))           
            if v.sample_gt == "0|1" or v.sample_gt == "1|0":
                hstate = v.sample_gt.split("|")[0]                 
                phase_set2h0_lst[v.sample_phase_set].append(hstate)
                phase_set2hpos_lst[v.sample_phase_set].append(v.pos)
                phase_set2hetsnp_lst[v.sample_phase_set].append((v.pos, v.ref, v.alt))
    return phase_set2h0_lst, phase_set2hpos_lst, phase_set2hetsnp_lst


def get_bit_complement(bits):
    bit_complement = "".join([bit_complement_hsh[bit] for bit in bits])
    return bit_complement    


def get_ccs_hbit(
    tpos2qbase: Dict[int, Tuple[str, int]], 
    hetsnp_subset_lst: List[Tuple[str, int, str, str]]
) -> List[str]:
    
    hbit = ""
    for hetsnp in hetsnp_subset_lst: # get read haplotype bits
        qbase, _ = tpos2qbase[hetsnp[0]]
        if qbase == hetsnp[1]:
            hbit += "0"
        elif qbase == hetsnp[2]:
            hbit += "1"
        else:
            hbit += "-"
    return hbit


def get_ccs_hap(
    ccs,
    h0_lst,
    hpos_lst,
    hetsnp_lst
) -> str: 

    idx = bisect.bisect_right(hpos_lst, ccs.tstart)
    jdx = bisect.bisect_right(hpos_lst, ccs.tend)
    if (jdx - idx) < 2:
        ccs_hap = "."
    else:
        hetsnp_subset_lst = [hetsnp_lst[kdx] for kdx in range(idx, jdx)]
        h0_hbit = "".join([h0_lst[kdx] for kdx in range(idx, jdx)])
        h1_hbit = get_bit_complement(h0_hbit)
        ccs_hbit = get_ccs_hbit(ccs.tpos2qbase, hetsnp_subset_lst)
        if h0_hbit == ccs_hbit:
            ccs_hap = "0"
        elif h1_hbit == ccs_hbit:
            ccs_hap = "1"
        else:
            ccs_hap = "."
    return ccs_hap


def sbs2hap(
    chrom: str,
    chrom_len: int,
    bam_file: str,
    phased_vcf_file: str,
    sbs_lst: List[Tuple[int, str, str]],
    sbs2ps: Dict[Tuple[int, str, str], int],
    chrom2sbs2hap_count: Dict[str, Dict[str, Tuple[str, str, str, int, int]]]
):

    sbs2hap_count = {}
    alignments = pysam.AlignmentFile(bam_file, "rb")
    ps2h0_lst, ps2hpos_lst, ps2hetsnp_lst = load_phase_set(phased_vcf_file, chrom, chrom_len)
    for sbs in sbs_lst:
        pos = sbs[0]
        ccs2hap = {} 
        ps = sbs2ps[sbs]
        h0_lst = ps2h0_lst[ps]
        hpos_lst = ps2hpos_lst[ps]
        hetsnp_lst = ps2hetsnp_lst[ps]
        hap2count = {"0": 0, "1": 0, ".": 0}
        for i in alignments.fetch(chrom, pos, pos+1):
            ccs = BAM(i)
            if not ccs.is_primary:
                continue
            ccs.cs2subindel() 
            ccs.cs2tpos2qbase()
            ccs_hap = get_ccs_hap(ccs, h0_lst, hpos_lst, hetsnp_lst)
            if sbs in ccs.tsbs_lst:
                ccs2hap[ccs.qname] = ccs_hap
            hap2count[ccs_hap] += 1
        som_hap_lst = list(set(ccs2hap.values()))
        if len(som_hap_lst) != 1:
            continue
        som_hap = som_hap_lst[0]     
        h0_count = hap2count["0"]  
        h1_count = hap2count["1"] 
        sbs2hap_count[sbs] = (som_hap, h0_count, h1_count) 
    chrom2sbs2hap_count[chrom] = sbs2hap_count


def dump_sbs2hap(
    bam_file: str,
    sbs_file: str,
    phased_vcf_file: str,
    threads: int,
    out_file: str
): 
    
    p = mp.Pool(threads)
    manager = mp.Manager()
    tname2tsize = get_tname2tsize(bam_file)
    chrom_lst, chrom2sbs2ps, chrom2sbs_lst = load_sbs(sbs_file) 
    chrom2sbs2hap_count = manager.dict()
    sbs2hap_arg_lst = [
        (
            chrom,
            tname2tsize[chrom],
            bam_file,
            phased_vcf_file,
            chrom2sbs_lst[chrom],
            chrom2sbs2ps[chrom],
            chrom2sbs2hap_count,
        )
        for chrom in chrom_lst
    ]
    p.starmap(
        sbs2hap, sbs2hap_arg_lst,
    )
    p.close()
    p.join()

    o = open(out_file, "w")
    o.write("{}\t{}\t{}\t{}\n".format("SBS", "SBS_HAP", "h0_count", "h1_count"))
    for chrom in chrom_lst:
        sbs2hap_count = chrom2sbs2hap_count[chrom]
        sbs_lst = natsort.natsorted(sbs2hap_count.keys())
        for sbs in sbs_lst:
            pos, ref, alt = sbs
            som_hap, h0_count, h1_count = sbs2hap_count[sbs]
            o.write("{}:{}_{}/{}\t{}\t{}\t{}\n".format(chrom, pos, ref, alt, som_hap, h0_count, h1_count)) 
    o.close()


def main():
    options = parse_args(sys.argv)
    dump_sbs2hap(
        options.bam, 
        options.sbs, 
        options.phased_vcf,
        options.threads, 
        options.out
    )
    sys.exit(0)


if __name__ == "__main__":
    main()
