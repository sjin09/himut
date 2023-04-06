#!/usr/bin/env python

import sys
import cyvcf2
import natsort
import pyfastx
import argparse
import himut.vcflib
from collections import defaultdict 
from typing import Set, Dict, List, Tuple
from himut.mutlib import purine, purine2pyrimidine, sbs48_lst, sbs96_lst, sbs96_to_tri, sbs96_to_sub, sbs96_to_sbs48 


def parse_args(args):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--vcf",
        type=str,
        required=True,
        help="VCF file to read"
    )
    parser.add_argument(
        "--ref",
        type=str,
        required=True,
        help="reference FASTA file to read"
    )
    parser.add_argument(
        "--tgt",
        type=str,
        required=True,
        help="list of chromosomes separated by new line"
    )
    parser.add_argument(
        "--som",
        required=False,
        action="store_true",
        help="VCF file has somatic mutations",
    )
    parser.add_argument(
        "--reference_sample",
        required=False,
        action="store_true",
        help="reads from the sample has been used to create the reference genome",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="file to return SBS48 counts"
    )
    args = args[1:]
    return parser.parse_args(args)


def load_target(tgt_file):
    tgt_lst = [line.strip() for line in open(tgt_file)]
    return tgt_lst


def get_sample(vcf_file):

    if vcf_file.endswith(".vcf"):
        for line in open(vcf_file).readlines():
            if line.startswith("#CHROM"):
                sample = line.strip().split()[-1]
                break
    elif vcf_file.endswith(".bgz"):
        v = cyvcf2.VCF(vcf_file)
        sample = v.samples[0]
    return sample


def get_tname2tsize(vcf_file):

    tname2tsize = {}
    if vcf_file.endswith(".vcf"):
        for line in open(vcf_file).readlines():
            if line.startswith("##"):
                if line.startswith("##contig"):
                    arr = line.strip().replace("##contig=<ID=", "").split(",")
                    tname = arr[0]
                    tsize = int(arr[1].replace("length=", "").replace(">", ""))
                    tname2tsize[tname] = tsize
            elif line.startswith("#CHROM"):
                break
    elif vcf_file.endswith(".bgz"):
        for line in cyvcf2.VCF(vcf_file).raw_header.split("\n"):
            if line.startswith("##"):
                if line.startswith("##contig"):
                    arr = line.strip().replace("##contig=<ID=", "").split(",")
                    tname = arr[0]
                    tsize = int(arr[1].replace("length=", "").replace(">", ""))
                    tname2tsize[tname] = tsize
            elif line.startswith("CHROM"):
                break
    return tname2tsize


def get_sbs96(
    chrom: str,
    pos: str,
    ref: str,
    alt: str,
    refseq,
) -> str:

    if ref in purine:
        upstream = purine2pyrimidine.get(refseq[chrom][pos + 1], "N")
        downstream = purine2pyrimidine.get(refseq[chrom][pos - 1], "N")
        sbs96 = "{}[{}>{}]{}".format(
            upstream,
            purine2pyrimidine.get(ref, "N"),
            purine2pyrimidine.get(alt, "N"),
            downstream,
        )
    else:
        upstream = refseq[chrom][pos - 1]
        downstream = refseq[chrom][pos + 1]
        sbs96 = "{}[{}>{}]{}".format(upstream, ref, alt, downstream)
    return sbs96


def load_somatic_sbs96_counts(
    vcf_file: str,
    ref_file: str,
    tname_lst: List[str]
):

    refseq = pyfastx.Fasta(ref_file)
    tname2sbs2counts = {tname: defaultdict(lambda: 0) for tname in tname_lst}
    if vcf_file.endswith(".vcf"):
        for line in open(vcf_file).readlines():
            if line.startswith("#"):
                continue
            v = himut.vcflib.VCF(line)
            if v.is_snp and v.is_pass and v.is_biallelic:
                sbs96 = get_sbs96(v.chrom, int(v.pos) - 1, v.ref, v.alt, refseq)
                tname2sbs2counts[v.chrom][sbs96] += 1
    elif vcf_file.endswith(".vcf.bgz"):
        for line in cyvcf2.VCF(vcf_file):
            v = himut.vcflib.VCF(str(line))
            if v.is_snp and v.is_pass and v.is_biallelic:
                sbs96 = get_sbs96(v.chrom, int(v.pos) - 1, v.ref, v.alt, refseq)
                tname2sbs2counts[v.chrom][sbs96] += 1
    return tname2sbs2counts


def load_germline_sbs96_counts(
    vcf_file: str,
    ref_file: str,
    tname_lst: List[str],
    sample_state: bool
):

    refseq = pyfastx.Fasta(ref_file)
    tname2sbs2counts = {tname: defaultdict(lambda: 0) for tname in tname_lst}
    if sample_state:
        if vcf_file.endswith(".vcf"):
            for line in open(vcf_file).readlines():
                if line.startswith("#"):
                    continue
                v = himut.vcflib.VCF(line)
                if v.is_snp and v.is_pass and v.is_biallelic and (v.sample_gt == "0/1" or v.sample_gt == "1/0"):
                    sbs96 = get_sbs96(v.chrom, int(v.pos) - 1, v.ref, v.alt, refseq)
                    tname2sbs2counts[v.chrom][sbs96] += 1
        elif vcf_file.endswith(".vcf.bgz"):
            for line in cyvcf2.VCF(vcf_file):
                v = himut.vcflib.VCF(str(line))
                if v.is_snp and v.is_pass and v.is_biallelic and (v.sample_gt == "0/1" or v.sample_gt == "1/0"):
                    sbs96 = get_sbs96(v.chrom, int(v.pos) - 1, v.ref, v.alt, refseq)
                    tname2sbs2counts[v.chrom][sbs96] += 1
    else:
        if vcf_file.endswith(".vcf"):
            for line in open(vcf_file).readlines():
                if line.startswith("#"):
                    continue
                v = himut.vcflib.VCF(line)
                if v.is_snp and v.is_pass and v.is_biallelic: 
                    sbs96 = get_sbs96(v.chrom, int(v.pos) - 1, v.ref, v.alt, refseq)
                    tname2sbs2counts[v.chrom][sbs96] += 1
        elif vcf_file.endswith(".vcf.bgz"):
            for line in cyvcf2.VCF(vcf_file):
                v = himut.vcflib.VCF(str(line))
                if v.is_snp and v.is_pass and v.is_biallelic:
                    sbs96 = get_sbs96(v.chrom, int(v.pos) - 1, v.ref, v.alt, refseq)
                    tname2sbs2counts[v.chrom][sbs96] += 1
    return tname2sbs2counts


def load_sbs96_counts(
    vcf_file: str, 
    ref_file: str,
    tname_lst: List[str],
    som_state: bool,
    sample_state: bool,
) -> Dict[str, int]:


    if som_state: 
        tname2sbs2counts = load_somatic_sbs96_counts(vcf_file, ref_file, tname_lst)
    else:
        tname2sbs2counts = load_germline_sbs96_counts(vcf_file, ref_file, tname_lst, sample_state)
    return tname2sbs2counts


def load_sbs48_counts(chrom_lst, tname2sbs96_counts):

    sbs48_counts = defaultdict(lambda: 0)
    for chrom in chrom_lst:
        for sbs96 in sbs96_lst:
            count = tname2sbs96_counts[chrom][sbs96]
            sbs48 = sbs96_to_sbs48[sbs96]
            sbs48_counts[sbs48] += count 
    return sbs48_counts


def dump_sbs48_counts(vcf_file, ref_file, tgt_file, som_state, sample_state, outfile):

    tname2tsize = get_tname2tsize(vcf_file)
    tname_lst = natsort.natsorted(list(tname2tsize.keys()))
    chrom_lst = [line.strip() for line in open(tgt_file).readlines()]
    tname2sbs96_counts = load_sbs96_counts(vcf_file, ref_file, tname_lst, som_state, sample_state)
    sbs48_counts = load_sbs48_counts(chrom_lst, tname2sbs96_counts)

    o = open(outfile, "w")
    o.write("sub\ttri\tsbs48\tcounts\tnormcounts\tref_tri_ratio\tref_ccs_tri_ratio\tref_tri_count\tref_callable_tri_count\tccs_callable_tri_count\n")
    for sbs48 in sbs48_lst:
        tri = sbs96_to_tri[sbs48]
        sub = sbs96_to_sub[sbs48]
        o.write("{}\t{}\t{}\t.\t{}\t.\t.\t.\t.\t.\n".format(sub, tri, sbs48, sbs48_counts[sbs48])) 
    o.close()

        
def main():
    options = parse_args(sys.argv)
    dump_sbs48_counts(options.vcf, options.ref, options.tgt, options.som, options.reference_sample, options.output)
    sys.exit(0)


if __name__ == "__main__":
    main()
