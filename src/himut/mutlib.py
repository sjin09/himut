import cyvcf2
import natsort
import pyfastx
import itertools
import himut.util
import himut.vcflib
import numpy as np 
import pandas as pd
from plotnine import *
from collections import defaultdict
from typing import Dict, List, Tuple


sbs_lst = []
tri_set = set()
sbs2sub = {}
sbs2tri = {}
purine = set(["A", "G"])
pyrimidine = set(["T", "C"])
sub_lst = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
purine2pyrimidine = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
for sub in sub_lst:
    for upstream, downstream in itertools.product(list("ATGC"), repeat=2):
        ref, alt = sub.split(">")
        tri = "{}{}{}".format(upstream, ref, downstream)
        sbs = "{}[{}]{}".format(upstream, sub, downstream)
        sbs_lst.append(sbs)
        sbs2tri[sbs] = tri
        sbs2sub[sbs] = sub
        tri_set.add(tri)
sbs_lst = natsort.natsorted(sbs_lst)
tri_lst = natsort.natsorted(list(tri_set))


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


def load_sbs96_counts(
    vcf_file: str, 
    ref_file: str, 
    chrom_lst: List[str]
) -> Dict[str, int]:

    tname_lst = []
    tname2sbs2counts = defaultdict()
    refseq = pyfastx.Fasta(ref_file)
    if vcf_file.endswith(".vcf"):
        for line in open(vcf_file).readlines():
            if line.startswith("##"):
                if line.startswith("##contig"):
                    tname = line.strip().replace("##contig=<ID=", "").split(",")[0]
                    tname_lst.append(tname)
                continue
            elif line.startswith("#CHROM"):
                for tname in tname_lst:
                    tname2sbs2counts[tname] = defaultdict(lambda: 0)
                continue
            v = himut.vcflib.VCF(line)
            if v.is_snp and v.is_pass and v.is_biallelic:
                sbs96 = get_sbs96(v.chrom, int(v.pos) - 1, v.ref, v.alt, refseq)
                tname2sbs2counts[v.chrom][sbs96] += int(v.alt_count_arr[0])
    elif vcf_file.endswith(".vcf.bgz"):
        for line in cyvcf2.VCF(vcf_file).raw_header.split("\n"):
            if line.startswith("##"):
                if line.startswith("##contig"):
                    tname = line.replace("##contig=<ID=", "").split(",")[0]
                    tname_lst.append(tname)
                continue
            elif line.startswith("#CHROM"):
                for tname in tname_lst:
                    tname2sbs2counts[tname] = defaultdict(lambda: 0)
        for i in cyvcf2.VCF(vcf_file):
            v = himut.vcflib.VCF(str(i))
            if v.is_snp and v.is_pass and v.is_biallelic:
                sbs96 = get_sbs96(v.chrom, int(v.pos) - 1, v.ref, v.alt, refseq)
                tname2sbs2counts[v.chrom][sbs96] += int(v.alt_count_arr[0])

    sbs2count = defaultdict(lambda: 0)
    for chrom in chrom_lst:
        for sbs, count in tname2sbs2counts[chrom].items():
            sbs2count[sbs] += count
    return sbs2count


def dump_sbs96_counts(
    vcf_file: str, ref_file: str, chrom: str, chrom_fofn: str, out_file: str
) -> None:

    chrom2sbs2counts = load_sbs96_counts(vcf_file, ref_file)
    chrom_lst = himut.util.load_chrom(chrom, chrom_fofn)
    sbs2counts = defaultdict(lambda: 0)
    if len(chrom_lst) == 0:
        for chrom in chrom2sbs2counts:
            for sbs, counts in chrom2sbs2counts[chrom].items():
                sbs2counts[sbs] += counts
    else:
        for chrom in chrom_lst:
            for sbs, counts in chrom2sbs2counts[chrom].items():
                sbs2counts[sbs] += counts

    o = open(out_file, "w")
    o.write("{}\t{}\t{}\t{}\n".format("sub", "tri", "sbs96", "counts"))
    for sbs in sbs_lst:
        o.write(
            "{}\t{}\t{}\t{}\n".format(sbs2sub[sbs], sbs2tri[sbs], sbs, sbs2counts[sbs])
        )
    o.close()


def dump_sbs96_plt(infile: str, sample: str, outfile: str) -> None:

    if sample is None:
        print(
            "Sample cannot be None\nBAM or VCF file might be missing sample information"
        )
        himut.util.exit()

    df = pd.read_csv(infile, sep="\t")
    plot = (
        ggplot(df, aes(x="tri", y="counts", fill="sub"))
        + geom_bar(stat="identity")
        + theme_bw()
        + facet_grid(". ~ sub", scales="free")
        + scale_fill_manual(
            values=("#98D7EC", "#212121", "#FF003A", "#A6A6A6", "#83A603", "#F5ABCC")
        )
        + labs(x="\nTrinucleotide Context\n", y="\nCounts\n")
        + ggtitle(sample)
        + theme(
            legend_title=element_blank(),
            axis_text_x=element_text(
                family="monospace", size=10, angle=90, ha="center"
            ),
            text=element_text(size=12),
        )
    )
    plot.save(outfile, width=22, height=10)


def dump_norm_sbs96_plt(infile: str, sample: str, outfile: str) -> None:

    if sample is None:
        print(
            "Sample cannot be None\nBAM or VCF file might be missing sample information"
        )
        himut.util.exit()

    df = pd.read_csv(infile, sep="\t", comment="#")
    plot = (
        ggplot(df, aes(x="tri", y="normcounts", fill="sub"))
        + geom_bar(stat="identity")
        + theme_bw()
        + facet_grid(". ~ sub", scales="free")
        + scale_fill_manual(
            values=("#98D7EC", "#212121", "#FF003A", "#A6A6A6", "#83A603", "#F5ABCC")
        )
        + labs(x="\nTrinucleotide Context\n", y="\nCounts\n")
        + ggtitle(sample)
        + theme(
            legend_title=element_blank(),
            axis_text_x=element_text(
                family="monospace", size=10, angle=90, ha="center"
            ),
            text=element_text(size=12),
        )
    )
    plot.save(outfile, width=22, height=10)


def get_trifreq(tri2count: Dict[str, int]) -> Dict[str, float]:

    tri_sum = sum(tri2count.values())
    tri2freq = {tri: count / float(tri_sum) for tri, count in tri2count.items()}
    return tri2freq


def get_callable_trifreq_ratio(
    ref_tri2count: Dict[str, int],
    ccs_tri2count: Dict[str, int],
) -> Dict[str, float]:

    ref_tri2freq = get_trifreq(ref_tri2count)
    ccs_tri2freq = get_trifreq(ccs_tri2count)
    trifreq_ratio = {
        tri: ref_tri2freq[tri] / float(ccs_tri2freq[tri]) for tri in ref_tri2freq
    }
    return trifreq_ratio


def get_cumsum_tricounts(chrom2tri2count: Dict[str, Dict[str, int]]) -> Tuple[int, Dict[str, int]]:

    tri2count = defaultdict(lambda: 0)
    for chrom in chrom2tri2count:
        for tri in tri_lst:
            tri_count = chrom2tri2count[chrom][tri]
            tri2count[tri] += tri_count
    return tri2count


def get_normcounts_cmdline(
    bam_file: str,
    ref_file: str,
    sbs_file: str,
    vcf_file: str,
    phased_vcf_file: str,
    min_mapq: int,
    min_sequence_identity: float,
    min_hq_base_proportion: float,
    min_alignment_proportion: float,
    min_gq: int,
    min_bq: int,
    min_trim: float,
    mismatch_window: int,
    max_mismatch_count: int,
    min_ref_count: int,
    min_alt_count: int,
    min_hap_count: int,
    common_snps: str,
    panel_of_normals: str,
    somatic_snv_prior: float,
    germline_snv_prior: float,
    germline_indel_prior: float,
    threads: int,
    phase: bool,
    non_human_sample: bool,
    reference_sample: bool,
    out_file: str,
):

    param = "--min_mapq {} --min_sequence_identity {} --min_hq_base_proportion {} --min_alignment_proportion {} --min_gq {} --min_bq {} --min_ref_count {} --min_alt_count {} --min_hap_count {} --min_trim {} --mismatch_window {} --max_mismatch_count {} --somatic_snv_prior {} --germline_snv_prior {} --germline_indel_prior {} --threads {} -o {}".format(
        min_mapq,
        min_sequence_identity,
        min_hq_base_proportion,
        min_alignment_proportion,
        min_gq,
        min_bq,
        min_ref_count,
        min_alt_count,
        min_hap_count,
        min_trim,
        mismatch_window,
        max_mismatch_count,
        somatic_snv_prior,
        germline_snv_prior,
        germline_indel_prior,
        threads,
        out_file,
    )

    if phase and non_human_sample and reference_sample:
        cmdline = "##himut_command=himut normcounts -i {} --ref {} --sbs {} --vcf {} --phased_vcf {} {} --phase --non_human_sample --reference_sample".format(
            bam_file, ref_file, sbs_file, vcf_file, phased_vcf_file, param
        )
    elif phase and non_human_sample and not reference_sample:
        cmdline = "##himut_command=himut normcounts -i {} --ref {} --sbs {} --vcf {} --phased_vcf {} {} --phase --non_human_sample".format(
            bam_file, ref_file, sbs_file, vcf_file, phased_vcf_file, param
        )
    elif phase and not non_human_sample and not reference_sample:
        cmdline = "##himut_command=himut normcounts -i {} --ref {} --sbs {} --phased_vcf {} {} --common_snps {} --panel_of_normals {} --phase".format(
            bam_file,
            ref_file,
            sbs_file,
            phased_vcf_file,
            param,
            common_snps,
            panel_of_normals,
        )
    elif not phase and non_human_sample and reference_sample:
        cmdline = "##himut_command=himut normcounts -i {} --ref {} --sbs {} --vcf {} {} --non_human_sample --reference_sample".format(
            bam_file, ref_file, sbs_file, vcf_file, param
        )
    elif not phase and non_human_sample and not reference_sample:
        cmdline = "##himut_command=himut normcounts -i {} --ref {} --sbs {} --vcf {} {} --non_human_sample".format(
            bam_file, ref_file, sbs_file, vcf_file, param
        )
    elif not phase and not non_human_sample and not reference_sample:
        cmdline = "##himut_command=himut normcounts -i {} --ref {} --sbs {} --vcf {} {} --common_snps {} --panel_of_normals {}".format(
            bam_file,
            ref_file,
            sbs_file,
            vcf_file,
            param,
            common_snps,
            panel_of_normals,
        )
    return cmdline


def dump_normcounts(
    sbs2count: Dict[str, int],
    ref_tri2count: Dict[str, int],
    chrom2ref_callable_tri2count: Dict[str, int],
    chrom2ccs_callable_tri2count: Dict[str, int],
    cmdline: str,
    out_file: str,
) -> None:

    ref_callable_tri2count = get_cumsum_tricounts(chrom2ref_callable_tri2count)
    ccs_callable_tri2count = get_cumsum_tricounts(chrom2ccs_callable_tri2count)
    ref_tri2freq_ratio = get_callable_trifreq_ratio(ref_callable_tri2count, ref_tri2count)
    callable_tri2freq_ratio = get_callable_trifreq_ratio(ref_callable_tri2count, ccs_callable_tri2count)

    o = open(out_file, "w")
    o.write("{}\n".format(cmdline))
    o.write(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
            "sub",
            "tri",
            "sbs96",
            "counts",
            "normcounts",
            "ref_tri_ratio",
            "ref_ccs_tri_ratio",
            "ref_tri_count",
            "ref_callable_tri_count",
            "ccs_callable_tri_count",
        )
    )

    for sbs in sbs_lst:
        sub = sbs2sub[sbs]
        tri = sbs2tri[sbs]
        count = sbs2count[sbs]
        ref_tri_count = ref_tri2count[tri]
        ref_trifreq_ratio = ref_tri2freq_ratio[tri] 
        ref_callable_tricount = ref_callable_tri2count[tri]
        ccs_callable_tricount = ccs_callable_tri2count[tri]
        ref_ccs_trifreq_ratio = callable_tri2freq_ratio[tri]
        normcount = count * ref_trifreq_ratio * ref_ccs_trifreq_ratio
        o.write(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                sub,
                tri,
                sbs,
                count,
                normcount,
                ref_trifreq_ratio,
                ref_ccs_trifreq_ratio, 
                ref_tri_count,
                ref_callable_tricount,
                ccs_callable_tricount,
            )
        )
    o.close()


def load_ref_tricount(
    ref_file: str,
    tri_file: str,
    region_list: str,
    threads: int,
):
 
    tri2count = 0
    chrom_lst = [line.strip() for line in open(region_list).readlines()]
    if ref_file is None and tri_file is not None:
        tri2count = dict(line.strip().split() for line in open(tri_file).readlines())
    elif ref_file is not None and tri_file is None:
        refseq = pyfastx.Fasta(ref_file)
        tri2count = himut.reflib.get_genome_tricounts(refseq, chrom_lst, threads)
    elif ref_file is not None and tri_file is not None:
        tri2count = dict(line.strip().split() for line in open(tri_file).readlines())
    else:
        print("Please provide either --ref or --tri file")
        himut.util.exit()

    for tri, count in tri2count.items(): tri2count[tri] = int(count)
    return tri2count


def get_burden_per_cell(
    infile,
    ref_file,
    tri_file,
    region_list, 
    threads,
    outfile,
):

    tri2mut_rate = defaultdict(lambda: 0)  
    ref_tri2count = load_ref_tricount(ref_file, tri_file, region_list, threads) 
    for line in open(infile).readlines():
        if line.startswith("#"):
            continue
        elif line.startswith("sub"):
            continue
        _sub, tri, _sbs, _count, normcounts, _ref_tri_ratio, _ref_ccs_tri_ratio, _ref_tri_count, _ref_callable_tri_count, ccs_callable_tri_count = line.strip().split()
        tri2mut_rate[tri] += float(normcounts)/int(ccs_callable_tri_count)
        
    ref_tri2count = load_ref_tricount(ref_file, tri_file, region_list, threads) 
    total_mut_count = sum([tri2mut_rate[tri] * ref_tri2count[tri] for tri in tri_lst])
    burden_per_cell = total_mut_count * 2
    o = open(outfile, "w")
    o.write("{}\n".format(burden_per_cell))
    o.close()


def dump_norm_log(
    chrom_lst: List[str],
    chrom2norm_log: Dict[str, List[int]],
):

    row_names = [
        "num_ccs",
        "num_bases",
        "num_unphased_bases",
        "num_het_bases",
        "num_hetalt_bases",
        "num_homalt_bases",
        "num_homref_bases",
        "num_uncallable_bases",
        "num_md_filtered_bases",
        "num_ab_filtered_bases",
        "num_low_gq_bases",
        "num_pon_filtered_bases",  
        "num_pop_filtered_bases",  
        "num_callable_bases",
    ]
    ncol = len(chrom_lst)
    nrow = len(row_names)
    dt = np.zeros((nrow, ncol))
    for i, chrom in enumerate(chrom_lst):
        for j, count in enumerate(chrom2norm_log[chrom]): 
            dt[j][i] = count

    o = open("norm.log", "w")
    col_lst = chrom_lst + ["total"] 
    o.write("{:30}{}\n".format("", "\t".join(col_lst)))
    for k in range(nrow):
        rsum =  str(int(np.sum(dt[k])))
        rlst = [str(int(r)) for r in dt[k].tolist()] + [rsum]
        o.write("{:30}{}\n".format(row_names[k], "\t".join(rlst)))
    o.close() 

# dbs2dbs78 = {
#     "AA>CC":("TT","GG"),
#     "AA>CG":("TT","CG"),
#     "AA>CT":("TT","AG"),
#     "AA>GC":("TT","GC"),
#     "AA>GG":("TT","CC"),
#     "AA>GT":("TT","AC"),
#     "AA>TC":("TT","GA"),
#     "AA>TG":("TT","CA"),
#     "AA>TT":("TT","AA"),
#     "AC>CA":("AC","CA"),
#     "AC>CG":("AC","CG"),
#     "AC>CT":("AC","CT"),
#     "AC>GA":("AC","GA"),
#     "AC>GG":("AC","GG"),
#     "AC>GT":("AC","GT"),
#     "AC>TA":("AC","TA"),
#     "AC>TG":("AC","TG"),
#     "AC>TT":("AC","TT"),
#     "AG>CA":("CT","TG"),
#     "AG>CC":("CT","GG"),
#     "AG>CT":("CT","AG"),
#     "AG>GA":("CT","TC"),
#     "AG>GC":("CT","GC"),
#     "AG>GT":("CT","AC"),
#     "AG>TA":("CT","TA"),
#     "AG>TC":("CT","GA"),
#     "AG>TT":("CT","AA"),
#     "AT>CA":("AT","CA"),
#     "AT>CC":("AT","CC"),
#     "AT>CG":("AT","CG"),
#     "AT>GA":("AT","GA"),
#     "AT>GC":("AT","GC"),
#     "AT>TA":("AT","TA"),
#     "CA>AC":("TG","GT"),
#     "CA>AG":("TG","CT"),
#     "CA>AT":("TG","AT"),
#     "CA>GC":("TG","GC"),
#     "CA>GG":("TG","CC"),
#     "CA>GT":("TG","AC"),
#     "CA>TC":("TG","GA"),
#     "CA>TG":("TG","CA"),
#     "CA>TT":("TG","AA"),
#     "CC>AA":("CC","AA"),
#     "CC>AG":("CC","AG"),
#     "CC>AT":("CC","AT"),
#     "CC>GA":("CC","GA"),
#     "CC>GG":("CC","GG"),
#     "CC>GT":("CC","GT"),
#     "CC>TA":("CC","TA"),
#     "CC>TG":("CC","TG"),
#     "CC>TT":("CC","TT"),
#     "CG>AA":("CG","TT"),
#     "CG>AC":("CG","GT"),
#     "CG>AT":("CG","AT"),
#     "CG>GA":("CG","TC"),
#     "CG>GC":("CG","GC"),
#     "CG>TA":("CG","TA"),
#     "GA>AC":("TC","GT"),
#     "GA>AG":("TC","CT"),
#     "GA>AT":("TC","AT"),
#     "GA>CC":("TC","GG"),
#     "GA>CG":("TC","CG"),
#     "GA>CT":("TC","AG"),
#     "GA>TC":("TC","GA"),
#     "GA>TG":("TC","CA"),
#     "GA>TT":("TC","AA"),
#     "GC>AA":("GC","AA"),
#     "GC>AG":("GC","AG"),
#     "GC>AT":("GC","AT"),
#     "GC>CA":("GC","CA"),
#     "GC>CG":("GC","CG"),
#     "GC>TA":("GC","TA"),
#     "TA>AC":("TA","GT"),
#     "TA>AG":("TA","CT"),
#     "TA>AT":("TA","AT"),
#     "TA>CC":("TA","GG"),
#     "TA>CG":("TA","CG"),
#     "TA>GC":("TA","GC"),
#     "TT>GG":("TT","GG"),
#     "TT>CG":("TT","CG"),
#     "TT>AG":("TT","AG"),
#     "TT>GC":("TT","GC"),
#     "TT>CC":("TT","CC"),
#     "TT>AC":("TT","AC"),
#     "TT>GA":("TT","GA"),
#     "TT>CA":("TT","CA"),
#     "TT>AA":("TT","AA"),
#     "GT>TG":("AC","CA"),
#     "GT>CG":("AC","CG"),
#     "GT>AG":("AC","CT"),
#     "GT>TC":("AC","GA"),
#     "GT>CC":("AC","GG"),
#     "GT>AC":("AC","GT"),
#     "GT>TA":("AC","TA"),
#     "GT>CA":("AC","TG"),
#     "GT>AA":("AC","TT"),
#     "CT>TG":("CT","TG"),
#     "CT>GG":("CT","GG"),
#     "CT>AG":("CT","AG"),
#     "CT>TC":("CT","TC"),
#     "CT>GC":("CT","GC"),
#     "CT>AC":("CT","AC"),
#     "CT>TA":("CT","TA"),
#     "CT>GA":("CT","GA"),
#     "CT>AA":("CT","AA"),
#     "AT>TG":("AT","CA"),
#     "AT>GG":("AT","CC"),
#     "AT>CG":("AT","CG"),
#     "AT>TC":("AT","GA"),
#     "AT>GC":("AT","GC"),
#     "AT>TA":("AT","TA"),
#     "TG>GT":("TG","GT"),
#     "TG>CT":("TG","CT"),
#     "TG>AT":("TG","AT"),
#     "TG>GC":("TG","GC"),
#     "TG>CC":("TG","CC"),
#     "TG>AC":("TG","AC"),
#     "TG>GA":("TG","GA"),
#     "TG>CA":("TG","CA"),
#     "TG>AA":("TG","AA"),
#     "GG>TT":("CC","AA"),
#     "GG>CT":("CC","AG"),
#     "GG>AT":("CC","AT"),
#     "GG>TC":("CC","GA"),
#     "GG>CC":("CC","GG"),
#     "GG>AC":("CC","GT"),
#     "GG>TA":("CC","TA"),
#     "GG>CA":("CC","TG"),
#     "GG>AA":("CC","TT"),
#     "CG>TT":("CG","TT"),
#     "CG>GT":("CG","GT"),
#     "CG>AT":("CG","AT"),
#     "CG>TC":("CG","TC"),
#     "CG>GC":("CG","GC"),
#     "CG>TA":("CG","TA"),
#     "TC>GT":("TC","GT"),
#     "TC>CT":("TC","CT"),
#     "TC>AT":("TC","AT"),
#     "TC>GG":("TC","GG"),
#     "TC>CG":("TC","CG"),
#     "TC>AG":("TC","AG"),
#     "TC>GA":("TC","GA"),
#     "TC>CA":("TC","CA"),
#     "TC>AA":("TC","AA"),
#     "GC>TT":("GC","AA"),
#     "GC>CT":("GC","AG"),
#     "GC>AT":("GC","AT"),
#     "GC>TG":("GC","CA"),
#     "GC>CG":("GC","CG"),
#     "GC>TA":("GC","TA"),
#     "TA>GT":("TA","GT"),
#     "TA>CT":("TA","CT"),
#     "TA>AT":("TA","AT"),
#     "TA>GG":("TA","GG"),
#     "TA>CG":("TA","CG"),
#     "TA>GC":("TA","GC")
# }


# def load_dbs78_counts(vcf_file):

#     chrom_lst = []
#     chrom2dbs2counts = defaultdict(dict)
#     if vcf_file.endswith(".vcf"):
#         for line in open(vcf_file).readlines():
#             if line.startswith("##"):
#                 if line.startswith("##contig"):
#                     chrom = line.strip().replace("##contig=<ID=", "").split(",")[0]
#                     chrom_lst.append(chrom)
#                 continue
#             elif line.startswith("#CHROM"):
#                 for chrom in chrom_lst:
#                     chrom2dbs2counts[chrom] = defaultdict(lambda: 0)
#                 continue
#             v = himut.vcflib.VCF(line)
#             if v.is_dbs and v.is_pass:
#                 chrom2dbs2counts[v.chrom][dbs2dbs78["{}>{}".format(v.ref, v.alt)]] += 1
#     elif vcf_file.endswith(".vcf.bgz"):
#         for line in cyvcf2.VCF(vcf_file).raw_header.split("\n"):
#             if line.startswith("##"):
#                 if line.startswith("##contig"):
#                     chrom = line.replace("##contig=<ID=", "").split(",")[0]
#                     chrom_lst.append(chrom)
#                 continue
#             elif line.startswith("#CHROM"):
#                 for chrom in chrom_lst:
#                     chrom2dbs2counts[chrom] = defaultdict(lambda: 0)
#         for i in cyvcf2.VCF(vcf_file):
#             v = himut.vcflib.VCF(str(i))
#             if v.is_dbs and v.is_pass:
#                 chrom2dbs2counts[v.chrom][dbs2dbs78["{}>{}".format(v.ref, v.alt)]] += 1
#     return chrom2dbs2counts


# def dump_dbs78_counts(vcf_file: str, chrom: str, chrom_fofn: str, out_file: str) -> None:

#     chrom2dbs2counts = load_dbs78_counts(vcf_file)
#     chrom_lst = himut.util.load_chrom(chrom, chrom_fofn)
#     dbs2counts = {dbs78: 0 for dbs78 in set(dbs2dbs78.values())}
#     if len(chrom_lst) == 0:
#         for chrom in chrom2dbs2counts:
#             for dbs, counts in chrom2dbs2counts[chrom].items():
#                 dbs2counts[dbs] += counts
#     else:
#         for chrom in chrom_lst:
#             for dbs, counts in chrom2dbs2counts[chrom].items():
#                 dbs2counts[dbs] += counts

#     o = open(out_file, "w")
#     o.write("{}\t{}\t{}\n".format("col", "dbs", "counts"))
#     for (ref, alt), counts in dbs2counts.items():
#         o.write("{}>NN\t{}\t{}\n".format(ref, alt, counts))
#     o.close()


# def dump_dbs78_plt(infile: str, sample: str, outfile: str) -> None:

#     if sample is None:
#         print("Sample cannot be None\nBAM or VCF file might be missing sample information")
#         himut.util.exit()

#     df = pd.read_csv(infile, sep="\t")
#     plot = (ggplot(df, aes(x="dbs", y="counts", fill="col")) +
#       geom_bar(stat="identity")  +
#       theme_bw() +
#       facet_grid(". ~ col", scales = "free") +
#       scale_fill_manual(values = ("#57F2F2","#0141FF","#60FF26","#2E8021","#FFCDC7","#FF2403", "#FFB04F", "#FF8D03", "#E2A8FF", "#61048F")) +
#       labs(x = "\nDouble base substitutions\n", y = "\nCounts\n") +
#       ggtitle(sample) +
#       theme(
#           legend_title = element_blank(),
#           axis_text_x = element_text(family = "monospace", size = 10, angle = 90, ha="center"),
#           text = element_text(size=12)
#         )
#     )
#     plot.save(outfile, width = 22, height = 10)
