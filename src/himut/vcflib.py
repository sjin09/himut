import tabix
import cyvcf2
import pyfastx
import natsort
import numpy as np
import himut.util
import himut.bamlib
from datetime import datetime
from collections import defaultdict
from typing import Dict, List, Set, Tuple


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
        hsh = {i: j for i, j in zip(self.format_lst, self.sample_format_lst)}
        if "GT" in hsh:
            self.sample_gt = hsh["GT"]
        if "PS" in hsh:
            self.sample_phase_set = hsh["PS"]
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
            elif len(self.ref) > len(self.alt):  # del
                self.is_indel = True
            elif len(self.ref) < len(self.alt):  # ins
                self.is_indel = True
        else:
            self.is_biallelic = False


def get_sample(vcf_file: str) -> str:
    if vcf_file.endswith(".vcf"):
        for line in open(vcf_file).readlines():
            if line.startswith("#CHROM"):
                sample = line.strip().split()[-1]
                return sample
    elif vcf_file.endswith(".bgz"):
        v = cyvcf2.VCF(vcf_file)
        sample = v.samples[0]
        return sample


def get_phased_vcf_header(
    bam_file: str,
    vcf_file: str,
    region: str,
    region_list: str,
    tname2tsize: Dict[str, int],
    min_bq: int,
    min_mapq: int,
    min_p_value: float,
    min_phase_proportion: float,
    threads: int,
    version: str,
    out_file,
) -> str:

    vcf_header_lst = [
        "##fileformat=VCFv4.2",
        '##FILTER=<ID=PASS,Description="All filters passed">',
        "##fileDate={}".format(datetime.now().strftime("%d%m%Y")),
        "##source=himut",
        "##source_version={}".format(version),
        "##content=himut somatic single base substitutions",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Conditional genotype quality">',
        '##FORMAT=<ID=BQ,Number=1,Type=Float,Description="Average base quality">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Read depth for each allele">',
        '##FORMAT=<ID=VAF,Number=A,Type=Float,Description="Variant allele fractions">',
        '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled genotype likelihoods rounded to the closest integer">',
        '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">',
    ]
    for tname in natsort.natsorted(list(tname2tsize.keys())):
        vcf_header_lst.append(
            "##contig=<ID={},length={}>".format(tname, tname2tsize[tname])
        )

    if region is None and region_list is not None:
        region_param = "--region_list {}".format(region_list)
    elif region is not None and region_list is None:
        region_param = "--region {}".format(region)
    elif region is not None and region_list is not None:
        region_param = "--region_list {}".format(region_list)

    cmdline = "##himut_command=himut phase -i {} --vcf {} {} --min_bq {} --min_mapq {} --min_p_value {} --min_phase_proportion {} --threads {} -o {}".format(
        bam_file,
        vcf_file,
        region_param,
        min_bq,
        min_mapq,
        min_p_value,
        min_phase_proportion,
        threads,
        out_file,
    )
    vcf_header_lst.append(cmdline)
    vcf_header_lst.append(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}".format(
            himut.bamlib.get_sample(bam_file)
        )
    )
    vcf_header = "\n".join(vcf_header_lst)
    return vcf_header


def get_himut_vcf_header(
    bam_file: str,
    vcf_file: str,
    phased_vcf_file: str,
    region: str,
    region_list: str,
    tname2tsize: Dict[str, int],
    common_snps: str,
    panel_of_normals: str,
    min_qv: int,
    min_mapq: int,
    qlen_lower_limit: float, 
    qlen_upper_limit: float,
    min_sequence_identity: float,
    min_gq: int,
    min_bq: int,
    min_trim: float,
    mismatch_window: int,
    max_mismatch_count: int,
    md_threshold: int,
    min_ref_count: int,
    min_alt_count: int,
    min_hap_count: int,
    threads: int,
    somatic_snv_prior: float,
    germline_snv_prior: float,
    germline_indel_prior: float,
    phase: bool,
    non_human_sample: bool,
    reference_sample: bool,
    create_panel_of_normals: bool,
    version: str,
    out_file: str,
) -> str:

    vcf_header_lst = [
        "##fileformat=VCFv4.2",
        "##fileDate={}".format(datetime.now().strftime("%d%m%Y")),
        "##source=himut",
        "##source_version={}".format(version),
        "##content=himut somatic single base substitutions",
        '##FILTER=<ID=PASS,Description="All filters passed">',
        '##FILTER=<ID=LowBQ,Description="Base quality score is below the minimum base quality score of {}">'.format(
            min_bq
        ),
        '##FILTER=<ID=LowGQ,Description="Germline genotype quality score is below the minimum genotype quality score of {}">'.format(
            min_gq
        ),
        '##FILTER=<ID=IndelSite,Description="Somatic substitution at indel site is not considered">',
        '##FILTER=<ID=HetSite,Description="Somatic substitution at heterzygous SNP site is not considered">',
        '##FILTER=<ID=HetAltSite,Description="Somatic substitution at tri-allelic SNP site is not considered">',
        '##FILTER=<ID=HomAltSite,Description="Somatic substitution at homozygous alternative SNP site is not considered">',
        '##FILTER=<ID=ComSnp,Description="Substitution is potentially from genomic DNA contamination">',
        '##FILTER=<ID=PanelOfNormal,Description="Substitution is found within the Panel of Normal VCF file">',
        '##FILTER=<ID=LowDepth,Description="Read depth is below the minimum reference allele and/or alterantive allele depth threshold">',
        '##FILTER=<ID=HighDepth,Description="Read depth is above the maximum depth threshold of {:.1f}">'.format(
            md_threshold
        ),
        '##FILTER=<ID=Trimmed,Description="Substitution is positioned near the end of reads">',
        '##FILTER=<ID=MismatchConflict,Description="Substitution is found next to a mismatch within a given mismatch window">',
        '##FILTER=<ID=Unphased,Description="CCS read is not haplotype phased">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=GQ,Number=1,Type=String,Description="Genotype quality score">',
        '##FORMAT=<ID=BQ,Number=1,Type=Float,Description="Base quality score">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Read depth for each allele">',
        '##FORMAT=<ID=VAF,Number=A,Type=Float,Description="Variant allele fractions">',
        '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">',
    ]
    for tname in natsort.natsorted(list(tname2tsize.keys())):
        vcf_header_lst.append(
            "##contig=<ID={},length={}>".format(tname, tname2tsize[tname])
        )

    if region is None and region_list is not None:
        region_param = "--region_list {}".format(region_list)
    elif region is not None and region_list is None:
        region_param = "--region {}".format(region)
    elif region is not None and region_list is not None:
        region_param = "--region_list {}".format(region_list)

    if phase:
            param = "{} --min_qv {} --min_mapq {} --qlen_lower_limit {} --qlen_upper_limit {} --min_gq {} --min_bq {} --min_trim {} --mismatch_window {} --max_mismatch_count {} --min_ref_count {} --min_alt_count {} --min_hap_count {} --somatic_snv_prior {} --germline_snv_prior {} --germline_indel_prior {} --threads {} -o {}".format(
            region_param,
            min_qv,
            min_mapq,
            qlen_lower_limit,
            qlen_upper_limit,
            min_sequence_identity,
            min_gq,
            min_bq,
            min_trim,
            mismatch_window,
            max_mismatch_count,
            min_ref_count,
            min_alt_count,
            min_hap_count,
            somatic_snv_prior,
            germline_snv_prior,
            germline_indel_prior,
            threads,
            out_file,
    )
    else:
        param = "{} --min_qv {} --min_mapq {} --qlen_lower_limit {} --qlen_upper_limit {} --min_sequence_identity {} --min_gq {} --min_bq {} --min_trim {} --mismatch_window {} --max_mismatch_count {} --min_ref_count {} --min_alt_count {} --somatic_snv_prior {} --germline_snv_prior {} --germline_indel_prior {} --threads {} -o {}".format(
            region_param,
            min_qv,
            min_mapq,
            qlen_lower_limit,
            qlen_upper_limit,
            min_sequence_identity,
            min_gq,
            min_bq,
            min_trim,
            mismatch_window,
            max_mismatch_count,
            min_ref_count,
            min_alt_count,
            somatic_snv_prior,
            germline_snv_prior,
            germline_indel_prior,
            threads,
            out_file,
    )
    if phase and non_human_sample and reference_sample: 
        cmdline = "##himut_command=himut call -i {} --vcf {} --phased_vcf {} {} --phase --non_human_sample --reference_sample".format(
            bam_file, 
            vcf_file, 
            phased_vcf_file, 
            param
        )
    elif phase and non_human_sample and not reference_sample:
        cmdline = "##himut_command=himut call -i {} --vcf {} --phased_vcf {} {} --phase --non_human_sample".format(
            bam_file, 
            vcf_file, 
            phased_vcf_file, 
            param
        )
    elif phase and not non_human_sample and not reference_sample:
        cmdline = "##himut_command=himut call -i {} --phased_vcf {} {} --common_snps {} --panel_of_normals {} --phase".format(
            bam_file,
            phased_vcf_file,
            param,
            common_snps,
            panel_of_normals,
        )
    elif (
        not phase
        and non_human_sample
        and reference_sample
        and not create_panel_of_normals
    ):
        cmdline = "##himut_command=himut call -i {} --vcf {} {} --non_human_sample --reference_sample".format(
            bam_file, 
            vcf_file, 
            param
        )
    elif (
        not phase
        and non_human_sample
        and not reference_sample
        and not create_panel_of_normals
    ):
        cmdline = (
            "##himut_command=himut call -i {} --vcf {} {} --non_human_sample".format(
                bam_file, 
                vcf_file,
                param
            )
        )
    elif (
        not phase
        and not non_human_sample
        and not reference_sample
        and create_panel_of_normals
    ):
        cmdline = "##himut_command=himut call -i {} --vcf {} {} --create_panel_of_normals".format(
            bam_file, 
            vcf_file, 
            param
        )
    elif (
        not phase
        and not non_human_sample
        and not reference_sample
        and not create_panel_of_normals
    ):
        cmdline = "##himut_command=himut call -i {} --vcf {} {} --common_snps {} --panel_of_normals {}".format(
            bam_file, 
            vcf_file, 
            param, 
            common_snps, 
            panel_of_normals
        )
    vcf_header_lst.append(cmdline)
    vcf_header_lst.append(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}".format(
            himut.bamlib.get_sample(bam_file)
        )
    )
    vcf_header = "\n".join(vcf_header_lst)
    return vcf_header


def load_snp(chrom: str, vcf_file: str) -> Set[Tuple[int, str, str]]:

    snp_set = set()
    for line in open(vcf_file):
        if line.startswith("#"):
            continue
        v = VCF(line)
        if v.chrom != chrom:
            continue
        if v.is_pass:
            if v.is_biallelic:
                if v.is_snp:
                    snp_set.add((v.pos, v.ref, v.alt))
            else:
                for alt in v.alt_lst:
                    if len(v.ref) == 1 and len(alt) == 1:
                        snp_set.add((v.pos, v.ref, alt))
    return snp_set


def load_bgz_snp(
    loci: Tuple[str, int, int], vcf_file: str
) -> Set[Tuple[int, str, str]]:

    snp_set = set()
    tb = tabix.open(vcf_file)
    records = tb.query(*loci)
    for record in records:
        v = VCF("\t".join(record))
        if v.is_pass:
            if v.is_biallelic:
                if v.is_snp:
                    snp_set.add((v.pos, v.ref, v.alt))
            else:
                for alt in v.alt_lst:
                    if len(v.ref) == 1 and len(alt) == 1:
                        snp_set.add((v.pos, v.ref, alt))
    return snp_set


def load_pon(
    chrom: str, vcf_file: str
) -> Tuple[Set[Tuple[str, int, str, str]], Set[Tuple[str, int, str, str]]]:

    sbs_set = set()
    for line in open(vcf_file):
        if line.startswith("#"):
            continue
        v = VCF(line)
        if chrom != v.chrom:
            continue
        if v.is_snp and v.is_pass and v.is_biallelic:
            sbs_set.add((v.pos, v.ref, v.alt))
    return sbs_set


def load_bgz_pon(
    loci: Tuple[str, int, int], vcf_file: str
) -> Tuple[Set[Tuple[str, int, str, str]], Set[Tuple[str, int, str, str]]]:

    sbs_set = set()
    tb = tabix.open(vcf_file)
    records = tb.query(*loci)
    for record in records:
        v = VCF("\t".join(record))
        if v.is_snp and v.is_pass and v.is_biallelic:
            sbs_set.add((v.pos, v.ref, v.alt))
    return sbs_set


def load_common_snp(chrom: str, vcf_file: str) -> Set[Tuple[str, int, str, str]]:

    snp_set = set()
    for line in open(vcf_file).readlines():
        if line.startswith("#"):
            continue
        arr = line.strip().split()
        alt_lst = arr[4].split(",")
        if chrom != arr[0] and arr[6] == "PASS" and len(alt_lst) == 1:
            pos = int(arr[1])
            ref = arr[3]
            alt = alt_lst[0]
            if len(ref) == 1 and len(alt) == 1:
                snp_set.add((pos, ref, alt))
    return snp_set


def load_bgz_common_snp(
    loci: Tuple[str, int, int],
    vcf_file: str,
) -> Set[Tuple[str, int, str, str]]:

    snp_set = set()
    tb = tabix.open(vcf_file)
    records = tb.query(*loci)
    for arr in records:
        alt_lst = arr[4].split(",")
        if arr[6] == "PASS" and len(alt_lst) == 1:
            pos = int(arr[1])
            ref = arr[3]
            alt = alt_lst[0]
            if len(ref) == 1 and len(alt) == 1:
                snp_set.add((pos, ref, alt))
    return snp_set


def load_hetsnps(
    vcf_file: str,
    chrom: str,
    chrom_len: int,
):

    hetsnp_lst = []
    hidx2hetsnp = {}
    hetsnp2hidx = {}
    if vcf_file.endswith(".vcf"):
        for line in open(vcf_file):
            if line.startswith("#"):
                continue
            v = VCF(line)
            if (
                chrom == v.chrom
                and v.is_snp
                and v.is_pass
                and v.is_biallelic
                and (v.sample_gt == "0/1" or v.sample_gt == "1/0")
            ):
                hetsnp_lst.append((v.pos, v.ref, v.alt))
    elif vcf_file.endswith(".bgz"):
        tb = tabix.open(vcf_file)
        records = tb.query(chrom, 0, chrom_len)
        for record in records:
            v = VCF("\t".join(record))
            if (
                v.is_snp
                and v.is_pass
                and v.is_biallelic
                and (v.sample_gt == "0/1" or v.sample_gt == "1/0")
            ):
                hetsnp_lst.append((v.pos, v.ref, v.alt))
    for hidx, hetsnp in enumerate(hetsnp_lst):
        hidx2hetsnp[hidx] = hetsnp
        hetsnp2hidx[hetsnp] = hidx
    return hetsnp_lst, hidx2hetsnp, hetsnp2hidx


def load_hblock_hsh(
    vcf_file: str,
    chrom_lst: List[str],
    tname2tsize: Dict[str, int],
    chrom2hblock_lst: Dict[str, List[Tuple[str, str]]],
):

    chrom2hetsnp_lst = defaultdict(list)
    chrom2hidx2hetsnp = defaultdict(dict)
    if vcf_file.endswith(".vcf"):
        for line in open(vcf_file):
            if line.startswith("#"):
                continue
            v = VCF(line)
            if (
                v.is_snp
                and v.is_pass
                and v.is_biallelic
                and (v.sample_gt == "0/1" or v.sample_gt == "1/0")
            ):
                chrom2hetsnp_lst[v.chrom].append((v.pos, v.ref, v.alt))
    elif vcf_file.endswith(".bgz"):
        tb = tabix.open(vcf_file)
        for chrom in chrom_lst:
            records = tb.query(chrom, 0, tname2tsize[chrom])
            for record in records:
                v = VCF("\t".join(record))
                if (
                    v.is_snp
                    and v.is_pass
                    and v.is_biallelic
                    and (v.sample_gt == "0/1" or v.sample_gt == "1/0")
                ):
                    chrom2hetsnp_lst[v.chrom].append((v.pos, v.ref, v.alt))

    hetsnp2hstate = {}
    hetsnp2phase_set = {}
    for chrom, hetsnp_lst in chrom2hetsnp_lst.items():
        for hidx, hetsnp in enumerate(hetsnp_lst):
            chrom2hidx2hetsnp[chrom][hidx] = hetsnp

    for chrom, hblock_lst in chrom2hblock_lst.items():
        for hblock in hblock_lst:
            phase_set = chrom2hidx2hetsnp[chrom][hblock[0][0]][0]
            for (hidx, hstate) in hblock:
                hetsnp = chrom2hidx2hetsnp[chrom][hidx]
                hetsnp2hstate[hetsnp] = hstate
                hetsnp2phase_set[hetsnp] = phase_set
    return hetsnp2hstate, hetsnp2phase_set


def load_germline_counts(vcf_file: str, chrom_lst: List[str]):

    chrom2hetsnp_count = defaultdict(lambda: 0)
    chrom2homsnp_count = defaultdict(lambda: 0)
    chrom2hetindel_count = defaultdict(lambda: 0)
    chrom2homindel_count = defaultdict(lambda: 0)
    if vcf_file.endswith(".vcf"):
        for line in open(vcf_file).readlines():
            if line.startswith("#"):
                continue
            v = himut.vcflib.VCF(line)
            if v.is_pass and v.is_biallelic:
                if v.is_snp:
                    if v.sample_gt == "0/1":
                        chrom2hetsnp_count[v.chrom] += 1
                    elif v.sample_gt == "1/1":
                        chrom2homsnp_count[v.chrom] += 1
                elif v.is_indel:
                    if v.sample_gt == "0/1":
                        chrom2hetindel_count[v.chrom] += 1
                    elif v.sample_gt == "1/1":
                        chrom2homindel_count[v.chrom] += 1
    elif vcf_file.endswith(".vcf.bgz"):
        for i in cyvcf2.VCF(vcf_file):
            v = himut.vcflib.VCF(str(i))
            if v.is_pass and v.is_biallelic:
                if v.is_snp:
                    if v.sample_gt == "0/1":
                        chrom2hetsnp_count[v.chrom] += 1
                    elif v.sample_gt == "1/1":
                        chrom2homsnp_count[v.chrom] += 1
                elif v.is_indel:
                    if v.sample_gt == "0/1":
                        chrom2hetindel_count[v.chrom] += 1
                    elif v.sample_gt == "1/1":
                        chrom2homindel_count[v.chrom] += 1
    hetsnp_count = sum([chrom2hetsnp_count[chrom] for chrom in chrom_lst])
    homsnp_count = sum([chrom2homsnp_count[chrom] for chrom in chrom_lst])
    hetindel_count = sum([chrom2hetindel_count[chrom] for chrom in chrom_lst])
    homindel_count = sum([chrom2homindel_count[chrom] for chrom in chrom_lst])
    return hetsnp_count, homsnp_count, hetindel_count, homindel_count


def get_germline_priors(
    chrom_lst: List[str], ref_file: str, vcf_file: str, reference_sample: bool
):

    refseq = pyfastx.Fasta(ref_file)
    target_sum = sum([len(str(refseq[chrom])) for chrom in chrom_lst])
    hetsnp_count, homsnp_count, hetindel_count, homindel_count = load_germline_counts(
        vcf_file, chrom_lst
    )
    if reference_sample:
        germline_snv_frequency = (2 * hetsnp_count) / target_sum
        germline_indel_frequency = (2 * hetindel_count) / target_sum
    else:
        germline_snv_frequency = (hetsnp_count + homsnp_count) / target_sum
        germline_indel_frequency = (hetindel_count + homindel_count) / target_sum
    germline_snv_prior = himut.util.get_truncated_float(germline_snv_frequency)
    germline_indel_prior = himut.util.get_truncated_float(germline_indel_frequency)
    return germline_snv_prior, germline_indel_prior



def load_phased_hetsnps(
    vcf_file: str,
    chrom_lst: List[str],
    tname2tsize: Dict[str, int]
) -> Dict[str, List[List[Tuple[int, str]]]]:

    chrom2phase_set2hbit_lst = defaultdict()
    chrom2phase_set2hpos_lst = defaultdict()
    chrom2phase_set2hetsnp_lst = defaultdict()
    for tname in tname2tsize:
        chrom2phase_set2hbit_lst[tname] = defaultdict(list)
        chrom2phase_set2hpos_lst[tname] = defaultdict(list) 
        chrom2phase_set2hetsnp_lst[tname] = defaultdict(list)
    if vcf_file.endswith(".vcf"):
        for line in open(vcf_file):
            if line.startswith("#"):
                continue
            v = VCF(line)
            if chrom == v.chrom and (v.sample_gt == "0|1" or v.sample_gt == "1|0"):
                chrom2phase_set2hpos_lst[v.chrom][v.sample_phase_set].append(v.pos)
                chrom2phase_set2hetsnp_lst[v.chrom][v.sample_phase_set].append((v.pos, v.ref, v.alt))
                chrom2phase_set2hbit_lst[v.chrom][v.sample_phase_set].append(v.sample_gt.split("|")[0])
    elif vcf_file.endswith(".bgz"):
        tb = tabix.open(vcf_file)
        for chrom in chrom_lst:
            records = tb.query(chrom, 0, tname2tsize[chrom])
            for record in records:
                v = VCF("\t".join(record))
                if v.sample_gt == "0|1" or v.sample_gt == "1|0":
                    chrom2phase_set2hpos_lst[v.chrom][v.sample_phase_set].append(v.pos)
                    chrom2phase_set2hetsnp_lst[v.chrom][v.sample_phase_set].append((v.pos, v.ref, v.alt))
                    chrom2phase_set2hbit_lst[v.chrom][v.sample_phase_set].append(v.sample_gt.split("|")[0])

    chrom2chunkloci_lst = {} 
    chrom_set = set(chrom_lst)
    for tname, tlen in tname2tsize.items():
        if tname not in chrom_set:
            del chrom2phase_set2hbit_lst[tname] 
            del chrom2phase_set2hpos_lst[tname] 
            del chrom2phase_set2hetsnp_lst[tname] 
            continue
        chrom2chunkloci_lst[tname] = []
        for _ps, hpos_lst in chrom2phase_set2hpos_lst[tname].items():
            start = hpos_lst[0]
            end = hpos_lst[-1]
            chrom2chunkloci_lst[tname].append((tname, start, end))
    return chrom2phase_set2hbit_lst, chrom2phase_set2hpos_lst, chrom2phase_set2hetsnp_lst, chrom2chunkloci_lst, 


def get_thresholds(sbs_file: str) -> Tuple[float, float ,float]:


    hsh = {}
    if sbs_file.endswith(".vcf"):
        for line in open(sbs_file).readlines():
            if line.startswith("##FILTER=<ID=HighDepth"):
                md_threshold = line.strip().split()[-1].replace('">', "")
            elif line.startswith("##himut_command"):
                arr = line.strip().replace("##himut_command=himut call", "").split() 
                for x, y in enumerate(arr):
                    z = x%2
                    if z == 0:
                        k = y
                    elif z == 1 and not y.startswith("--"):
                        hsh[k] = y
            elif line.startswith("#CHROM"):
                break
    elif sbs_file.endswith(".vcf.bgz"):
        for line in cyvcf2.VCF(sbs_file).raw_header.split("\n"):
            if line.startswith("##FILTER=<ID=HighDepth"):
                md_threshold = line.strip().split()[-1].replace('">', "")
            elif line.startswith("##himut_command"):
                arr = line.strip().replace("##himut_command=himut call", "").split() 
                for x, y in enumerate(arr):
                    z = x%2
                    if z == 0:
                        k = y
                    elif z == 1 and not y.startswith("--"):
                        hsh[k] = y
            elif line.startswith("#CHROM"):
                break
    qlen_lower_limit = int(hsh["--qlen_lower_limit"])
    qlen_upper_limit = int(hsh["--qlen_upper_limit"])
    return qlen_lower_limit, qlen_upper_limit, float(md_threshold)



def dump_phased_hetsnps(
    bam_file: str,
    vcf_file: str,
    region: str,
    region_list: str,
    tname2tsize: Dict[str, int],
    min_bq: int,
    min_mapq: int,
    min_p_value: float,
    min_phase_proportion: float,
    threads: int,
    chrom_lst: List[str],
    chrom2hblock_lst: List[List[Tuple[int, int]]],
    version: str,
    out_file: str,
) -> None:

    o = open(out_file, "w")
    o.write(
        "{}\n".format(
            himut.vcflib.get_phased_vcf_header(
                bam_file,
                vcf_file,
                region,
                region_list,
                tname2tsize,
                min_bq,
                min_mapq,
                min_p_value,
                min_phase_proportion,
                threads,
                version,
                out_file,
            )
        )
    )
    hetsnp2hstate, hetsnp2phase_set = load_hblock_hsh(
        vcf_file, chrom_lst, tname2tsize, chrom2hblock_lst
    )
    if vcf_file.endswith(".vcf"):
        for line in open(vcf_file):
            if line.startswith("#"):
                continue
            v = VCF(line)
            if (
                v.is_snp
                and v.is_pass
                and v.is_biallelic
                and (v.sample_gt == "0/1" or v.sample_gt == "1/0")
            ):
                hetsnp = (v.pos, v.ref, v.alt)
                fmt, sample_fmt = line.strip().split()[-2:]
                if hetsnp in hetsnp2hstate:
                    phase_set = hetsnp2phase_set[hetsnp]
                    sample_gt = "0|1" if hetsnp2hstate[hetsnp] == "0" else "1|0"
                    sample_fmt_arr = sample_fmt.split(":")
                    sample_fmt_arr[0] = sample_gt
                    o.write(
                        "{}\t{}\t.\t{}\t{}\t{}\tPASS\t.\t{}:PS\t{}:{}\n".format(
                            v.chrom,
                            v.pos,
                            v.ref,
                            v.alt,
                            v.qual,
                            fmt,
                            ":".join(sample_fmt_arr),
                            phase_set,
                        )
                    )
                else:
                    o.write(
                        "{}\t{}\t.\t{}\t{}\t{}\tPASS\t.\t{}:PS\t{}:.\n".format(
                            v.chrom, v.pos, v.ref, v.alt, v.qual, fmt, sample_fmt
                        )
                    )
    elif vcf_file.endswith(".bgz"):
        tb = tabix.open(vcf_file)
        for chrom in chrom_lst:
            records = tb.query(chrom, 0, tname2tsize[chrom])
            for record in records:
                line = "\t".join(record)
                v = VCF(line)
                fmt, sample_fmt = line.strip().split()[-2:]
                if (
                    v.is_snp
                    and v.is_pass
                    and v.is_biallelic
                    and (v.sample_gt == "0/1" or v.sample_gt == "1/0")
                ):
                    hetsnp = (v.pos, v.ref, v.alt)
                    if hetsnp in hetsnp2hstate:
                        phase_set = hetsnp2phase_set[hetsnp]
                        sample_gt = "0|1" if hetsnp2hstate[hetsnp] == "0" else "1|0"
                        sample_fmt_arr = sample_fmt.split(":")
                        sample_fmt_arr[0] = sample_gt
                        o.write(
                            "{}\t{}\t.\t{}\t{}\t{}\tPASS\t.\t{}:PS\t{}:{}\n".format(
                                v.chrom,
                                v.pos,
                                v.ref,
                                v.alt,
                                v.qual,
                                fmt,
                                ":".join(sample_fmt_arr),
                                phase_set,
                            )
                        )
                    else:
                        o.write(
                            "{}\t{}\t.\t{}\t{}\t{}\tPASS\t.\t{}:PS\t{}:.\n".format(
                                v.chrom, v.pos, v.ref, v.alt, v.qual, fmt, sample_fmt
                            )
                        )
    o.close()


def dump_sbs(
    vcf_file: str,
    vcf_header: str,
    chrom_lst: List[str],
    chrom2tsbs_lst: Dict[
        str, List[Tuple[str, int, str, str, str, int, int, int, float]]
    ],
):

    if not vcf_file.endswith(".vcf"):
        print("VCF file must have .vcf suffix")
        himut.util.exit()

    o = open(vcf_file, "w")
    p = open(vcf_file.replace(".vcf", ".single_molecule_mutations.vcf"), "w")
    o.write("{}\n".format(vcf_header))
    p.write("{}\n".format(vcf_header))
    for chrom in chrom_lst:
        for (
            chrom,
            pos,
            ref,
            alt,
            status,
            gq,
            bq,
            read_depth,
            ref_count,
            alt_count,
            vaf,
            _phase_set,
        ) in chrom2tsbs_lst[chrom]:
            if status == "HetAltSite":
                o.write(
                    "{}\t{}\t.\t{}\t{}\t.\t{}\t.\tGT:GQ:BQ:DP:AD:VAF\t./.:{}:{}:{:0.0f}:{:0.0f},{}:{}\n".format(
                        chrom,
                        pos,
                        ref,
                        alt,
                        status,
                        gq,
                        bq,
                        read_depth,
                        ref_count,
                        alt_count,
                        vaf,
                    )
                )
                if int(ref_count) == 1:
                    p.write(
                        "{}\t{}\t.\t{}\t{}\t.\t{}\t.\tGT:GQ:BQ:DP:AD:VAF\t./.:{}:{}:{:0.0f}:{:0.0f},{}:{}\n".format(
                            chrom,
                            pos,
                            ref,
                            alt,
                            status,
                            gq,
                            bq,
                            read_depth,
                            ref_count,
                            alt_count,
                            vaf,
                        )
                    )
            else:
                o.write(
                    "{}\t{}\t.\t{}\t{}\t.\t{}\t.\tGT:GQ:BQ:DP:AD:VAF\t./.:{}:{:0.1f}:{:0.0f}:{:0.0f},{:0.0f}:{:.2f}\n".format(
                        chrom,
                        pos,
                        ref,
                        alt,
                        status,
                        gq,
                        bq,
                        read_depth,
                        ref_count,
                        alt_count,
                        vaf,
                    )
                )
                if int(alt_count) == 1:
                    p.write(
                        "{}\t{}\t.\t{}\t{}\t.\t{}\t.\tGT:GQ:BQ:DP:AD:VAF\t./.:{}:{:0.1f}:{:0.0f}:{:0.0f},{:0.0f}:{:.2f}\n".format(
                            chrom,
                            pos,
                            ref,
                            alt,
                            status,
                            bq,
                            gq,
                            read_depth,
                            ref_count,
                            alt_count,
                            vaf,
                        )
                    )
    o.close()
    p.close()


def dump_phased_sbs(
    vcf_file: str,
    vcf_header: str,
    chrom_lst: List[str],
    chrom2tsbs_lst: Dict[
        str, List[Tuple[str, int, str, str, str, int, int, int, float]]
    ],
):

    if not vcf_file.endswith(".vcf"):
        print("VCF file must have .vcf suffix")
        himut.util.exit()

    o = open(vcf_file, "w")
    p = open(vcf_file.replace(".vcf", ".single_molecule_mutations.vcf"), "w")
    o.write("{}\n".format(vcf_header))
    p.write("{}\n".format(vcf_header))
    for chrom in chrom_lst:
        for (
            chrom,
            pos,
            ref,
            alt,
            status,
            gq,
            bq,
            read_depth,
            ref_count,
            alt_count,
            vaf,
            phase_set,
        ) in chrom2tsbs_lst[chrom]:
            if status == "HetAltSite":
                o.write(
                    "{}\t{}\t.\t{}\t{}\t.\t{}\t.\tGT:GQ:BQ:DP:AD:VAF\t./.:{}:{}:{:0.0f}:{:0.0f},{}:{}:{}\n".format(
                        chrom,
                        pos,
                        ref,
                        alt,
                        status,
                        gq,
                        bq,
                        read_depth,
                        ref_count,
                        alt_count,
                        vaf,
                        phase_set,
                    )
                )
                if int(ref_count) == 1:
                    p.write(
                        "{}\t{}\t.\t{}\t{}\t.\t{}\t.\tGT:GQ:BQ:DP:AD:VAF:PS\t./.:{}:{}:{:0.0f}:{:0.0f},{}:{}:{}\n".format(
                            chrom,
                            pos,
                            ref,
                            alt,
                            status,
                            gq,
                            bq,
                            read_depth,
                            ref_count,
                            alt_count,
                            vaf,
                            phase_set,
                        )
                    )
            else:
                o.write(
                    "{}\t{}\t.\t{}\t{}\t.\t{}\t.\tGT:GQ:BQ:DP:AD:VAF:PS\t./.:{}:{:0.1f}:{:0.0f}:{:0.0f},{:0.0f}:{:.2f}:{}\n".format(
                        chrom,
                        pos,
                        ref,
                        alt,
                        status,
                        gq,
                        bq,
                        read_depth,
                        ref_count,
                        alt_count,
                        vaf,
                        phase_set,
                    )
                )
                if int(alt_count) == 1:
                    p.write(
                        "{}\t{}\t.\t{}\t{}\t.\t{}\t.\tGT:GQ:BQ:DP:AD:VAF:PS\t./.:{}:{:0.1f}:{:0.0f}:{:0.0f},{:0.0f}:{:.2f}:{}\n".format(
                            chrom,
                            pos,
                            ref,
                            alt,
                            status,
                            gq,
                            bq,
                            read_depth,
                            ref_count,
                            alt_count,
                            vaf,
                            phase_set,
                        )
                    )
    o.close()
    p.close()


def dump_call_log(
    chrom_lst: List[str], 
    chrom2tsbs_log: Dict[str, List[int]], 
) -> None:

    row_names = [
        "num_ccs",
        "num_sbs",
        "num_het_sbs",
        "num_hetalt_sbs",
        "num_homalt_sbs",
        "num_somrev_sbs",
        "num_homref_sbs",
        "num_uncallable_sbs",
        "num_low_gq_sbs",
        "num_low_bq_sbs",
        "num_pon_filtered_sbs",
        "num_pop_filtered_sbs",
        "num_md_filtered_sbs",
        "num_ab_filtered_sbs",
        "num_som",
    ]
    ncol = len(chrom_lst)
    nrow = len(row_names)
    dt = np.zeros((nrow, ncol))
    for i, chrom in enumerate(chrom_lst):
        for j, count in enumerate(chrom2tsbs_log[chrom]): 
            dt[j][i] = count
    
    o = open("himut.log", "w")
    col_lst = chrom_lst + ["total"] 
    o.write("{:30}{}\n".format("", "\t".join(col_lst)))
    for k in range(nrow):
        rsum =  str(int(np.sum(dt[k])))
        rlst = [str(int(r)) for r in dt[k].tolist()] + [rsum]
        o.write("{:30}{}\n".format(row_names[k], "\t".join(rlst)))
    o.close()