import cyvcf2
import natsort
import pyfastx
import itertools
import himut.util
import himut.vcflib
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
        tri = upstream + ref + downstream
        sbs = "{}[{}]{}".format(upstream, sub, downstream)
        tri_set.add(tri)
        sbs_lst.append(sbs)
        sbs2sub[sbs] = sub
        sbs2tri[sbs] = tri
sbs_lst = natsort.natsorted(sbs_lst)
tri_lst = natsort.natsorted(list(tri_set))

dbs2dbs78 = {
    "AA>CC":("TT","GG"),
    "AA>CG":("TT","CG"),
    "AA>CT":("TT","AG"),
    "AA>GC":("TT","GC"),
    "AA>GG":("TT","CC"),
    "AA>GT":("TT","AC"),
    "AA>TC":("TT","GA"),
    "AA>TG":("TT","CA"),
    "AA>TT":("TT","AA"),
    "AC>CA":("AC","CA"),
    "AC>CG":("AC","CG"),
    "AC>CT":("AC","CT"),
    "AC>GA":("AC","GA"),
    "AC>GG":("AC","GG"),
    "AC>GT":("AC","GT"),
    "AC>TA":("AC","TA"),
    "AC>TG":("AC","TG"),
    "AC>TT":("AC","TT"),
    "AG>CA":("CT","TG"),
    "AG>CC":("CT","GG"),
    "AG>CT":("CT","AG"),
    "AG>GA":("CT","TC"),
    "AG>GC":("CT","GC"),
    "AG>GT":("CT","AC"),
    "AG>TA":("CT","TA"),
    "AG>TC":("CT","GA"),
    "AG>TT":("CT","AA"),
    "AT>CA":("AT","CA"),
    "AT>CC":("AT","CC"),
    "AT>CG":("AT","CG"),
    "AT>GA":("AT","GA"),
    "AT>GC":("AT","GC"),
    "AT>TA":("AT","TA"),
    "CA>AC":("TG","GT"),
    "CA>AG":("TG","CT"),
    "CA>AT":("TG","AT"),
    "CA>GC":("TG","GC"),
    "CA>GG":("TG","CC"),
    "CA>GT":("TG","AC"),
    "CA>TC":("TG","GA"),
    "CA>TG":("TG","CA"),
    "CA>TT":("TG","AA"),
    "CC>AA":("CC","AA"),
    "CC>AG":("CC","AG"),
    "CC>AT":("CC","AT"),
    "CC>GA":("CC","GA"),
    "CC>GG":("CC","GG"),
    "CC>GT":("CC","GT"),
    "CC>TA":("CC","TA"),
    "CC>TG":("CC","TG"),
    "CC>TT":("CC","TT"),
    "CG>AA":("CG","TT"),
    "CG>AC":("CG","GT"),
    "CG>AT":("CG","AT"),
    "CG>GA":("CG","TC"),
    "CG>GC":("CG","GC"),
    "CG>TA":("CG","TA"),
    "GA>AC":("TC","GT"),
    "GA>AG":("TC","CT"),
    "GA>AT":("TC","AT"),
    "GA>CC":("TC","GG"),
    "GA>CG":("TC","CG"),
    "GA>CT":("TC","AG"),
    "GA>TC":("TC","GA"),
    "GA>TG":("TC","CA"),
    "GA>TT":("TC","AA"),
    "GC>AA":("GC","AA"),
    "GC>AG":("GC","AG"),
    "GC>AT":("GC","AT"),
    "GC>CA":("GC","CA"),
    "GC>CG":("GC","CG"),
    "GC>TA":("GC","TA"),
    "TA>AC":("TA","GT"),
    "TA>AG":("TA","CT"),
    "TA>AT":("TA","AT"),
    "TA>CC":("TA","GG"),
    "TA>CG":("TA","CG"),
    "TA>GC":("TA","GC"),
    "TT>GG":("TT","GG"),
    "TT>CG":("TT","CG"),
    "TT>AG":("TT","AG"),
    "TT>GC":("TT","GC"),
    "TT>CC":("TT","CC"),
    "TT>AC":("TT","AC"),
    "TT>GA":("TT","GA"),
    "TT>CA":("TT","CA"),
    "TT>AA":("TT","AA"),
    "GT>TG":("AC","CA"),
    "GT>CG":("AC","CG"),
    "GT>AG":("AC","CT"),
    "GT>TC":("AC","GA"),
    "GT>CC":("AC","GG"),
    "GT>AC":("AC","GT"),
    "GT>TA":("AC","TA"),
    "GT>CA":("AC","TG"),
    "GT>AA":("AC","TT"),
    "CT>TG":("CT","TG"),
    "CT>GG":("CT","GG"),
    "CT>AG":("CT","AG"),
    "CT>TC":("CT","TC"),
    "CT>GC":("CT","GC"),
    "CT>AC":("CT","AC"),
    "CT>TA":("CT","TA"),
    "CT>GA":("CT","GA"),
    "CT>AA":("CT","AA"),
    "AT>TG":("AT","CA"),
    "AT>GG":("AT","CC"),
    "AT>CG":("AT","CG"),
    "AT>TC":("AT","GA"),
    "AT>GC":("AT","GC"),
    "AT>TA":("AT","TA"),
    "TG>GT":("TG","GT"),
    "TG>CT":("TG","CT"),
    "TG>AT":("TG","AT"),
    "TG>GC":("TG","GC"),
    "TG>CC":("TG","CC"),
    "TG>AC":("TG","AC"),
    "TG>GA":("TG","GA"),
    "TG>CA":("TG","CA"),
    "TG>AA":("TG","AA"),
    "GG>TT":("CC","AA"),
    "GG>CT":("CC","AG"),
    "GG>AT":("CC","AT"),
    "GG>TC":("CC","GA"),
    "GG>CC":("CC","GG"),
    "GG>AC":("CC","GT"),
    "GG>TA":("CC","TA"),
    "GG>CA":("CC","TG"),
    "GG>AA":("CC","TT"),
    "CG>TT":("CG","TT"),
    "CG>GT":("CG","GT"),
    "CG>AT":("CG","AT"),
    "CG>TC":("CG","TC"),
    "CG>GC":("CG","GC"),
    "CG>TA":("CG","TA"),
    "TC>GT":("TC","GT"),
    "TC>CT":("TC","CT"),
    "TC>AT":("TC","AT"),
    "TC>GG":("TC","GG"),
    "TC>CG":("TC","CG"),
    "TC>AG":("TC","AG"),
    "TC>GA":("TC","GA"),
    "TC>CA":("TC","CA"),
    "TC>AA":("TC","AA"),
    "GC>TT":("GC","AA"),
    "GC>CT":("GC","AG"),
    "GC>AT":("GC","AT"),
    "GC>TG":("GC","CA"),
    "GC>CG":("GC","CG"),
    "GC>TA":("GC","TA"),
    "TA>GT":("TA","GT"),
    "TA>CT":("TA","CT"),
    "TA>AT":("TA","AT"),
    "TA>GG":("TA","GG"),
    "TA>CG":("TA","CG"),
    "TA>GC":("TA","GC")
}


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
        sbs96 = "{}[{}>{}]{}".format(upstream, purine2pyrimidine.get(ref, "N"), purine2pyrimidine.get(alt, "N"), downstream)
    else:
        upstream = refseq[chrom][pos -  1]
        downstream = refseq[chrom][pos + 1]
        sbs96 = "{}[{}>{}]{}".format(upstream, ref, alt, downstream)    
    return sbs96


def load_sbs96_counts(vcf_file: str, ref_file: str) -> Dict[str, int]:

    chrom_lst = []
    refseq = pyfastx.Fasta(ref_file)
    chrom2sbs2counts = defaultdict()
    if vcf_file.endswith(".vcf"):
        for line in open(vcf_file).readlines():
            if line.startswith("##"):
                if line.startswith("##contig"):
                    chrom = line.strip().replace("##contig=<ID=", "").split(",")[0]
                    chrom_lst.append(chrom)
                continue
            elif line.startswith("#CHROM"):
                for chrom in chrom_lst:
                    chrom2sbs2counts[chrom] = defaultdict(lambda: 0) 
                continue
            v = himut.vcflib.VCF(line)
            if v.is_snp and v.is_pass:
                chrom2sbs2counts[v.chrom][get_sbs96(v.chrom, int(v.pos) - 1, v.ref, v.alt, refseq)] += 1
    elif vcf_file.endswith(".vcf.bgz"):
        for line in cyvcf2.VCF(vcf_file).raw_header.split("\n"):
            if line.startswith("##"):
                if line.startswith("##contig"):
                    chrom = line.replace("##contig=<ID=", "").split(",")[0]
                    chrom_lst.append(chrom)
                continue
            elif line.startswith("#CHROM"):
                for chrom in chrom_lst:
                    chrom2sbs2counts[chrom] = defaultdict(lambda: 0) 
        for i in cyvcf2.VCF(vcf_file):
            v = himut.vcflib.VCF(str(i))
            if v.is_snp and v.is_pass:
                chrom2sbs2counts[v.chrom][get_sbs96(v.chrom, int(v.pos) - 1, v.ref, v.alt, refseq)] += 1
    return chrom2sbs2counts


def dump_sbs96_counts(vcf_file: str, ref_file: str, chrom: str, chrom_fofn: str, out_file: str) -> None:


    himut.util.check_mutpatterns_input_exists(vcf_file, ref_file, chrom, chrom_fofn, out_file)
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
        o.write("{}\t{}\t{}\t{}\n".format(sbs2sub[sbs], sbs2tri[sbs], sbs, sbs2counts[sbs]))
    o.close()

    
def dump_sbs96_plt(infile: str, sample: str, outfile: str) -> None:

    if sample is None:
        print("Sample cannot be None\nBAM or VCF file might be missing sample information")
        himut.util.exit()
        
    df = pd.read_csv(infile, sep="\t")
    plot = (ggplot(df, aes(x="tri", y="counts", fill="sub")) +
      geom_bar(stat="identity")  +
      theme_bw() +
      facet_grid(". ~ sub", scales = "free") +
      scale_fill_manual(values = ("#98D7EC","#212121","#FF003A","#A6A6A6","#83A603","#F5ABCC")) +
      labs(x = "\nTrinucleotide Context\n", y = "\nCounts\n") +
      ggtitle(sample) +
      theme(
          legend_title = element_blank(),
          axis_text_x = element_text(family = "monospace", size = 10, angle = 90, ha="center"),
          text = element_text(size=12)
        )
    )
    plot.save(outfile, width = 22, height = 10)


def load_dbs78_counts(vcf_file):

    chrom_lst = []
    chrom2dbs2counts = defaultdict(dict)
    if vcf_file.endswith(".vcf"):
        for line in open(vcf_file).readlines():
            if line.startswith("##"):
                if line.startswith("##contig"):
                    chrom = line.strip().replace("##contig=<ID=", "").split(",")[0]
                    chrom_lst.append(chrom)
                continue
            elif line.startswith("#CHROM"):
                for chrom in chrom_lst:
                    chrom2dbs2counts[chrom] = defaultdict(lambda: 0) 
                continue
            v = himut.vcflib.VCF(line)
            if v.is_dbs and v.is_pass:
                chrom2dbs2counts[v.chrom][dbs2dbs78["{}>{}".format(v.ref, v.alt)]] += 1
    elif vcf_file.endswith(".vcf.bgz"):
        for line in cyvcf2.VCF(vcf_file).raw_header.split("\n"):
            if line.startswith("##"):
                if line.startswith("##contig"):
                    chrom = line.replace("##contig=<ID=", "").split(",")[0]
                    chrom_lst.append(chrom)
                continue
            elif line.startswith("#CHROM"):
                for chrom in chrom_lst:
                    chrom2dbs2counts[chrom] = defaultdict(lambda: 0) 
        for i in cyvcf2.VCF(vcf_file):
            v = himut.vcflib.VCF(str(i))
            if v.is_dbs and v.is_pass:
                chrom2dbs2counts[v.chrom][dbs2dbs78["{}>{}".format(v.ref, v.alt)]] += 1
    return chrom2dbs2counts


def dump_dbs78_counts(vcf_file: str, chrom: str, chrom_fofn: str, out_file: str) -> None:

    chrom2dbs2counts = load_dbs78_counts(vcf_file)
    chrom_lst = himut.util.load_chrom(chrom, chrom_fofn)
    dbs2counts = {dbs78: 0 for dbs78 in set(dbs2dbs78.values())}
    if len(chrom_lst) == 0:
        for chrom in chrom2dbs2counts:
            for dbs, counts in chrom2dbs2counts[chrom].items():
                dbs2counts[dbs] += counts
    else:
        for chrom in chrom_lst: 
            for dbs, counts in chrom2dbs2counts[chrom].items():
                dbs2counts[dbs] += counts

    o = open(out_file, "w")
    o.write("{}\t{}\t{}\n".format("col", "dbs", "counts"))
    for (ref, alt), counts in dbs2counts.items():
        o.write("{}>NN\t{}\t{}\n".format(ref, alt, counts))
    o.close()


def dump_dbs78_plt(infile: str, sample: str, outfile: str) -> None:

    if sample is None:
        print("Sample cannot be None\nBAM or VCF file might be missing sample information")
        himut.util.exit()
        
    df = pd.read_csv(infile, sep="\t")
    plot = (ggplot(df, aes(x="dbs", y="counts", fill="col")) +
      geom_bar(stat="identity")  +
      theme_bw() +
      facet_grid(". ~ col", scales = "free") +
      scale_fill_manual(values = ("#57F2F2","#0141FF","#60FF26","#2E8021","#FFCDC7","#FF2403", "#FFB04F", "#FF8D03", "#E2A8FF", "#61048F")) +
      labs(x = "\nDouble base substitutions\n", y = "\nCounts\n") +
      ggtitle(sample) +
      theme(
          legend_title = element_blank(),
          axis_text_x = element_text(family = "monospace", size = 10, angle = 90, ha="center"),
          text = element_text(size=12)
        )
    )
    plot.save(outfile, width = 22, height = 10)


def load_sbs_count(vcf_file: str, ref_file: str) -> Tuple[List[str], Dict[str, Dict[str, int]]]:

    chrom_set = set()
    refseq = pyfastx.Fasta(ref_file)
    sbs2count = defaultdict(lambda: 0)
    if vcf_file.endswith(".vcf"):
        for line in open(vcf_file).readlines():
            if line.startswith("#"):
                continue
            v = himut.vcflib.VCF(line)
            if v.is_pass:
                chrom_set.add(v.chrom)
                sbs2count[get_sbs96(v.chrom, int(v.pos) - 1, v.ref, v.alt, refseq)] += 1
    elif vcf_file.endswith(".bgz"):
        for i in cyvcf2.VCF(vcf_file):
            v = himut.vcflib.VCF(str(i))
            if v.is_pass:
                chrom_set.add(v.chrom)
                sbs2count[get_sbs96(v.chrom, int(v.pos) - 1, v.ref, v.alt, refseq)] += 1
    chrom_lst = natsort.natsorted(list(chrom_set))
    return chrom_lst, sbs2count


def get_trifreq(tri2count: Dict[str, int]) -> Dict[str, float]:

    tri_sum = sum(tri2count.values())
    tri2freq = {tri: count/float(tri_sum) for tri, count in tri2count.items()}
    return tri2freq


def get_trifreq_ratio(
    ref_tri2count: Dict[str, int],
    ccs_tri2count: Dict[str, int],
) -> Dict[str, float]:

    ref_tri2freq = get_trifreq(ref_tri2count)
    ccs_tri2freq = get_trifreq(ccs_tri2count)
    trifreq_ratio = {tri: ref_tri2freq[tri]/float(ccs_tri2freq[tri]) for tri in ref_tri2freq}
    return trifreq_ratio


def dump_normcounts(
    sbs2counts: Dict[str, int], 
    ref_tri2count: Dict[str, int],
    ccs_tri2count: Dict[str, int],
    out_file: str
) -> None:

    print("returning normalised counts")
    o = open(out_file, "w")
    tri2freq_ratio = get_trifreq_ratio(ref_tri2count, ccs_tri2count)
    o.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("sub", "tri", "sbs96", "counts", "tri_ratio", "normcounts", "ref_tri_count", "ccs_tri_count"))
    for sbs in sbs_lst:
        tri = sbs2tri[sbs]
        count = sbs2counts[sbs] 
        ref_tricount = ref_tri2count[tri]
        ccs_tricount = ccs_tri2count[tri]
        trifreq_ratio = tri2freq_ratio[tri]
        normcounts = count * trifreq_ratio
        o.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(sbs2sub[sbs], tri, sbs, count, trifreq_ratio, normcounts, ref_tricount, ccs_tricount))
    print("finished returning normalised counts")


