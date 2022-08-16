import tabix
import cyvcf2
import natsort
import numpy as np 
import himut.util
import himut.bamlib
from datetime import datetime
from collections import defaultdict
from typing import Dict, List, Set, Tuple


ROW_NAMES = [
    "num_reads", 
    "num_not_primary",  
    "num_low_mapq", 
    "num_abnormal", 
    "num_low_qv", 
    "num_low_identity", 
    "num_contamination", 
    "num_hq_reads", 
    "num_sub_candidates", 
    "num_bq_filtered_sbs", 
    "num_pon_filtered_sbs", 
    "num_trimmed_sbs", 
    "num_mismatch_filtered_sbs", 
    "num_uncallable_sbs", 
    "num_ab_filtered_sbs", 
    "num_md_filtered_sbs", 
    "num_sbs", 
    "num_unphased_sbs", 
    "num_phased_sbs"
]


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


def get_vcf_header(
    bam_file: str,
    version: str,
) -> str:

    vcf_header_lst = [
        "##fileformat=VCFv4.2",
        '##FILTER=<ID=PASS,Description="All filters passed">',
        "##fileDate={}".format(datetime.now().strftime("%d%m%Y")),
        "##source=himut",
        "##source_version={}".format(version),
        '##content=himut somatic single base substitutions',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Conditional genotype quality">',
        '##FORMAT=<ID=BQ,Number=1,Type=Float,Description="Average base quality">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Read depth for each allele">',
        '##FORMAT=<ID=VAF,Number=A,Type=Float,Description="Variant allele fractions">',
        '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled genotype likelihoods rounded to the closest integer">',
        '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">',
    ]

    tname_lst, tname2tsize = himut.bamlib.get_tname2tsize(bam_file)
    for tname in tname_lst:
        vcf_header_lst.append("##contig=<ID={},length={}>".format(tname, tname2tsize[tname]))

    vcf_header_lst.append(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}".format(himut.bamlib.get_sample(bam_file))
    )
    vcf_header = "\n".join(vcf_header_lst)
    return vcf_header


def get_himut_vcf_header(
    bam_file: str,
    vcf_file: str,
    phased_vcf_file: str,
    region: str,
    region_list: str,
    ploidy: str,
    tname2tsize: Dict[str, int],
    common_snps: str,
    panel_of_normals: str,
    min_mapq: int,
    min_trim: float,
    min_sequence_identity: float,
    min_hq_base_proportion: float,
    min_alignment_proportion: float,
    min_bq: int,
    mismatch_window: int,
    max_mismatch_count: int,
    min_ref_count: int,
    min_alt_count: int,
    md_threshold: int,
    min_hap_count: int,
    threads: int,
    phase: bool,
    non_human_sample: bool,
    create_panel_of_normals: bool,
    version: str,
    out_file: str,
) -> str:

    vcf_header_lst = [
        "##fileformat=VCFv4.2",
        "##fileDate={}".format(datetime.now().strftime("%d%m%Y")),
        "##source=himut",
        "##source_version={}".format(version),
        '##content=himut somatic single base substitutions',
        '##FILTER=<ID=PASS,Description="All filters passed">',
        '##FILTER=<ID=LowBQ,Description="Base quality score is below minimum base quality score of {}">'.format(min_bq),
        '##FILTER=<ID=CommonVariant,Description="Substitution is found within the provided VCF file">',
        '##FILTER=<ID=PanelOfNormal,Description="Substitution is found within the Panel of Normal VCF file">',
        '##FILTER=<ID=Trimmed,Description="Substitution are positioned near the end of reads">',
        '##FILTER=<ID=HighDepth,Description="Read depth is above the maximum depth threshold of {}">'.format(md_threshold),
        '##FILTER=<ID=InsufficientDepth,Description="Reference allele depth is below minimum allele balance threshold of {}">'.format(min_ref_count),
        '##FILTER=<ID=MismatchConflict,Description="A mismatch(es) is found near the substitution">',
        '##FILTER=<ID=IndelConflict,Description="CCS reads have insertions or deletions in the same position as the substitution">',
        '##FILTER=<ID=UnphasedRead,Description="CCS read is not assigned to a haplotype block">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=BQ,Number=1,Type=Float,Description="Average base quality">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Read depth for each allele">',
        '##FORMAT=<ID=VAF,Number=A,Type=Float,Description="Variant allele fractions">',
        '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">',
    ]
    tname_lst = natsort.natsorted(list(tname2tsize.keys()))
    for tname in tname_lst:
        vcf_header_lst.append("##contig=<ID={},length={}>".format(tname, tname2tsize[tname]))

    if region is None and region_list is not None:
        region_param = "--region_list {}".format(region_list)
    elif region is not None and region_list is None:
        region_param = "--region {}".format(region)
    elif region is not None and region_list is not None:
        region_param = "--region_list {}".format(region_list)

    cmdline = "{} --ploidy {} --min_mapq {} --min_trim {} --min_sequence_identity {} --min_hq_base_proportion {} --min_alignment_proportion {} --min_bq {} --min_ref_count {} --min_alt_count {} --min_hap_count {} --mismatch_window {} --max_mismatch_count {} --threads {} -o {}".format(
        region_param,
        ploidy,
        min_mapq,
        min_trim,
        min_sequence_identity,
        min_hq_base_proportion,
        min_alignment_proportion,
        min_bq,
        min_ref_count,
        min_alt_count,
        min_hap_count,
        mismatch_window,
        max_mismatch_count,
        threads,
        out_file 
    )
    if phase:
        if non_human_sample and not create_panel_of_normals:
            cmdline = "##himut_command=himut call -i {} --vcf {} --phased_vcf {} {} --phase --tree_life_sample".format(bam_file, vcf_file, phased_vcf_file, cmdline)
        elif not non_human_sample and create_panel_of_normals:
            cmdline = "##himut_command=himut call -i {} --vcf {} --phased_vcf {} {} --phase --create_panel_of_normals".format(bam_file, vcf_file, phased_vcf_file, cmdline)
        elif not non_human_sample and not create_panel_of_normals:
            cmdline = "##himut_command=himut call -i {} --vcf {} --phased_vcf {} --common_snps {} --panel_of_normals {} {} --phase".format(bam_file, vcf_file, phased_vcf_file, common_snps, panel_of_normals, cmdline)
    else:
        if non_human_sample and not create_panel_of_normals:
            cmdline = "##himut_command=himut call -i {} --vcf {} {} --tree_life_sample".format(bam_file, vcf_file, cmdline)
        elif not non_human_sample and create_panel_of_normals:
            cmdline = "##himut_command=himut call -i {} --vcf {} {} --create_panel_of_normals".format(bam_file, vcf_file, cmdline)
        elif not non_human_sample and not create_panel_of_normals:
            cmdline = "##himut_command=himut call -i {} --vcf {} --common_snps {} --panel_of_normals {} {}".format(bam_file, vcf_file, common_snps, panel_of_normals, cmdline)
    vcf_header_lst.append(cmdline)
    vcf_header_lst.append(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}".format(himut.bamlib.get_sample(bam_file))
    )
    vcf_header = "\n".join(vcf_header_lst)
    return vcf_header


def load_snp(
    chrom: str,
    vcf_file: str
) -> Set[Tuple[str, int, str, str]]:

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
                    snp_set.add((v.chrom, v.pos, v.ref, v.alt))
            else:
                for alt in v.alt_lst:
                    if len(v.ref) == 1 and len(alt) == 1:
                        snp_set.add((v.chrom, v.pos, v.ref, alt))
    return snp_set 
      
                     
def load_bgz_snp(
    loci: Tuple[str, int, int],
    vcf_file: str
) -> Set[Tuple[str, int, str, str]]:

    snp_set = set()
    tb = tabix.open(vcf_file)
    records = tb.query(*loci)
    for record in records:
        v = VCF("\t".join(record))            
        if v.is_pass:
            if v.is_biallelic:
                if v.is_snp:
                    snp_set.add((v.chrom, v.pos, v.ref, v.alt))
            else:
                for alt in v.alt_lst:
                    if len(v.ref) == 1 and len(alt) == 1:
                        snp_set.add((v.chrom, v.pos, v.ref, alt))
    return snp_set


def load_pon(
    chrom: str,
    vcf_file: str
) -> Tuple[Set[Tuple[str, int, str, str]], Set[Tuple[str, int, str, str]]]:

    state = 0
    sbs_set = set() 
    dbs_set = set()
    if vcf_file.endswith(".vcf"):
        for line in open(vcf_file):
            if line.startswith("#"):
                continue
            v = VCF(line)
            if chrom != v.chrom: 
                if state == 0:
                    continue
                else:
                    break
            else:
                state = 1
                if v.is_pass and v.is_biallelic:
                    if v.is_snp:
                        sbs_set.add((v.chrom, v.pos, v.ref, v.alt))
                    else:
                        if len(v.ref) == len(v.alt) == 2:
                            dbs_set.add((v.chrom, v.pos, v.ref, v.alt))
    return sbs_set, dbs_set


def load_bgz_pon(
    loci: Tuple[str, int, int],
    vcf_file: str
) -> Tuple[Set[Tuple[str, int, str, str]], Set[Tuple[str, int, str, str]]]:
   
    dbs_set = set() 
    sbs_set = set() 
    tb = tabix.open(vcf_file)
    records = tb.query(*loci)
    for record in records:
        v = VCF("\t".join(record))            
        if v.is_pass and v.is_biallelic:
            if v.is_snp:
                sbs_set.add((v.chrom, v.pos, v.ref, v.alt))
            else:
                if len(v.ref) == len(v.alt) == 2:
                    dbs_set.add((v.chrom, v.pos, v.ref, v.alt))
    return sbs_set, dbs_set


def load_common_snp(
    chrom: str,
    vcf_file: str
) -> Set[Tuple[str, int, str, str]]:
    
    snp_set = set()
    if vcf_file.endswith(".vcf"):
        for line in open(vcf_file).readlines():
            if line.startswith("#"):
                continue
            arr = line.strip().split()
            if chrom != arr[0]:
                continue
            if arr[6] == "PASS":
                alt_lst = arr[4].split(",") 
                if len(alt_lst) == 1:            
                    ref = arr[3]
                    alt = alt_lst[0]
                    if len(ref) == 1 and len(alt) == 1:
                        snp_set.add((chrom, int(arr[1]), ref, alt))
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
        if arr[6] == "PASS":
            if len(alt_lst) == 1:
                ref = arr[3]
                alt = alt_lst[0]
                if len(ref) == 1 and len(alt) == 1:
                    snp_set.add((arr[0], int(arr[1]), ref, alt))
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
            if v.is_snp and v.is_pass and v.is_biallelic and (v.sample_gt == "0/1" or v.sample_gt == "1/0"):
                hetsnp_lst.append((v.chrom, v.pos, v.ref, v.alt))
    elif vcf_file.endswith(".bgz"):
        tb = tabix.open(vcf_file)
        records = tb.query(chrom, 0, chrom_len)
        for record in records:
            v = VCF("\t".join(record))            
            if v.is_snp and v.is_pass and v.is_biallelic and (v.sample_gt == "0/1" or v.sample_gt == "1/0"):
                hetsnp_lst.append((v.chrom, v.pos, v.ref, v.alt))
    for hidx, hetsnp in enumerate(hetsnp_lst):
        hidx2hetsnp[hidx] = hetsnp
        hetsnp2hidx[hetsnp] = hidx
    return hetsnp_lst, hidx2hetsnp, hetsnp2hidx


def load_phased_hetsnp_index(
    vcf_file: str, 
    chrom: str,
    chrom_len: int,
) -> Dict[str, List[List[Tuple[int, str]]]]:

    hidx = 0  
    state = 0
    hidx2hetsnp = {}
    hetsnp2hidx = {}
    hblock_lst = defaultdict(list)
    phase_set2hblock = defaultdict(list) 
    if vcf_file.endswith(".vcf"):
        for line in open(vcf_file):
            if line.startswith("#"):
                continue
            v = VCF(line)
            if chrom != v.chrom:
                if state == 0:
                    continue
                elif state == 1:
                    break
            elif chrom == v.chrom: 
                if state == 0:
                    state = 1
                if v.sample_gt == "0|1" or v.sample_gt == "1|0":
                    hetsnp = (v.chrom, v.pos, v.ref, v.alt)
                    hidx2hetsnp[hidx] = hetsnp
                    hetsnp2hidx[hetsnp] = hidx
                    hstate = v.sample_gt.split("|")[0]                    
                    phase_set2hblock[v.sample_phase_set].append((hidx, hstate))
                hidx += 1 
    elif vcf_file.endswith(".bgz"):
        tb = tabix.open(vcf_file)
        records = tb.query(chrom, 0, chrom_len)
        for hidx, record in enumerate(records):
            v = VCF("\t".join(record))           
            if v.sample_gt == "0|1" or v.sample_gt == "1|0":
                hetsnp = (v.chrom, v.pos, v.ref, v.alt)
                hidx2hetsnp[hidx] = hetsnp
                hetsnp2hidx[hetsnp] = hidx                
                hstate = v.sample_gt.split("|")[0]
                phase_set2hblock[v.sample_phase_set].append((hidx, hstate))
    hblock_lst = [hblock for hblock in phase_set2hblock.values()]
    return hblock_lst, hidx2hetsnp, hetsnp2hidx


def get_phased_hetsnps(
    vcf_file: str, 
    chrom: str,
    chrom_len: int, 
):
    hetsnp2bidx = {}
    hidx2hstate = {} 
    filtered_hblock_lst = []
    hblock_lst, hidx2hetsnp, hetsnp2hidx = load_phased_hetsnp_index(vcf_file, chrom, chrom_len)
    for hblock in hblock_lst:
        ipos = 0
        filtered_hblock = []
        for (hidx, hstate) in hblock:
            jpos = hidx2hetsnp[hidx][1]
            if jpos - ipos > 1:
                filtered_hblock.append((hidx, hstate))
            ipos = jpos
        if len(filtered_hblock) > 1:
            filtered_hblock_lst.append(filtered_hblock)
   
    phased_hetsnp_lst = []
    for bidx, hblock in enumerate(filtered_hblock_lst):
        for (hidx, hstate) in hblock:
            hidx2hstate[hidx] = hstate 
            hetsnp = hidx2hetsnp[hidx]
            hetsnp2bidx[hetsnp] = bidx 
            phased_hetsnp_lst.append(hetsnp)
    phased_hetsnp_lst = natsort.natsorted(phased_hetsnp_lst) 
    hpos_lst = [hetsnp[1] for hetsnp in phased_hetsnp_lst]
    return hpos_lst, filtered_hblock_lst, phased_hetsnp_lst, hidx2hetsnp, hidx2hstate, hetsnp2bidx, hetsnp2hidx


def dump_phased_hetsnps(
    bam_file: str,
    vcf_file: str,
    chrom_lst: List[str],
    tname2tsize: Dict[str, int],
    chrom2hblock_lst: List[List[Tuple[int, int]]],
    version: str,
    out_file: str,
) -> None:

    hetsnp2hstate = {}
    hetsnp2phase_set = {}
    chrom2hetsnp_lst = defaultdict(list)
    chrom2hidx2hetsnp = defaultdict(dict)
    chrom2hetsnp2hidx = defaultdict(dict)
    if vcf_file.endswith(".vcf"):
        for line in open(vcf_file):
            if line.startswith("#"):
                continue
            v = VCF(line)
            if v.is_snp and v.is_pass and v.is_biallelic and (v.sample_gt == "0/1" or v.sample_gt == "1/0"):
                chrom2hetsnp_lst[v.chrom].append((v.chrom, v.pos, v.ref, v.alt))
    elif vcf_file.endswith(".bgz"):
        tb = tabix.open(vcf_file)
        for chrom in chrom_lst:
            records = tb.query(chrom, 0, tname2tsize[chrom])
            for record in records:
                v = VCF("\t".join(record))            
                if v.is_snp and v.is_pass and v.is_biallelic and (v.sample_gt == "0/1" or v.sample_gt == "1/0"):
                    chrom2hetsnp_lst[v.chrom].append((v.chrom, v.pos, v.ref, v.alt))
    for chrom, hetsnp_lst in chrom2hetsnp_lst.items():
        for hidx, hetsnp in enumerate(hetsnp_lst):
            chrom2hidx2hetsnp[chrom][hidx] = hetsnp
            chrom2hetsnp2hidx[chrom][hetsnp] = hidx
    del chrom2hetsnp_lst 
    
    for chrom, hblock_lst in chrom2hblock_lst.items():
        for hblock in hblock_lst:
            phase_set = chrom2hidx2hetsnp[chrom][hblock[0][0]][1]
            for (hidx, hstate) in hblock:
                hetsnp = chrom2hidx2hetsnp[chrom][hidx]
                hetsnp2hstate[hetsnp] = hstate
                hetsnp2phase_set[hetsnp] = phase_set 

    o = open(out_file, "w")
    o.write("{}\n".format(himut.vcflib.get_vcf_header(bam_file, version)))
    if vcf_file.endswith(".vcf"):
        for line in open(vcf_file):
            if line.startswith("#"):
                continue
            v = VCF(line)
            if v.is_snp and v.is_pass and v.is_biallelic and (v.sample_gt == "0/1" or v.sample_gt == "1/0"):
                hetsnp = (v.chrom, v.pos, v.ref, v.alt)
                fmt, sample_fmt = line.strip().split()[-2:]
                if hetsnp in hetsnp2hstate:
                    phase_set = hetsnp2phase_set[hetsnp]
                    sample_gt = "0|1" if hetsnp2hstate[hetsnp] == "0" else "1|0"
                    sample_fmt_arr = sample_fmt.split(":")
                    sample_fmt_arr[0] = sample_gt
                    o.write("{}\t{}\t.\t{}\t{}\t{}\tPASS\t.\t{}:PS\t{}:{}\n".format(v.chrom, v.pos, v.ref, v.alt, v.qual, fmt, ":".join(sample_fmt_arr), phase_set))
                else:
                    o.write("{}\t{}\t.\t{}\t{}\t{}\tPASS\t.\t{}:PS\t{}:.\n".format(v.chrom, v.pos, v.ref, v.alt, v.qual, fmt, sample_fmt))
    elif vcf_file.endswith(".bgz"):
        tb = tabix.open(vcf_file)
        for chrom in chrom_lst:
            records = tb.query(chrom, 0, tname2tsize[chrom])
            for record in records:
                line = "\t".join(record)
                v = VCF(line)            
                fmt, sample_fmt = line.strip().split()[-2:]
                if v.is_snp and v.is_pass and v.is_biallelic and (v.sample_gt == "0/1" or v.sample_gt == "1/0"):
                    hetsnp = (v.chrom, v.pos, v.ref, v.alt)
                    if hetsnp in hetsnp2hstate:
                        phase_set = hetsnp2phase_set[hetsnp]
                        sample_gt = "0|1" if hetsnp2hstate[hetsnp] == "0" else "1|0"
                        sample_fmt_arr = sample_fmt.split(":")
                        sample_fmt_arr[0] = sample_gt
                        o.write("{}\t{}\t.\t{}\t{}\t{}\tPASS\t.\t{}:PS\t{}:{}\n".format(v.chrom, v.pos, v.ref, v.alt, v.qual, fmt, ":".join(sample_fmt_arr), phase_set))
                    else:
                        o.write("{}\t{}\t.\t{}\t{}\t{}\tPASS\t.\t{}:PS\t{}:.\n".format(v.chrom, v.pos, v.ref, v.alt, v.qual, fmt, sample_fmt))
    o.close()


def dump_hblock_statistics(
    chrom_lst: List[str],
    chrom2hblock_statistics: Dict[str, Tuple[int, int, int, int, int]],
    log_file: str, 
) -> None:

    o = open(log_file, "w")
    o.write(
        "chrom\thetsnp_count\tblock_count\tsmallest_block\tlargest_block\tshortest_block (bp)\tlongest_block (bp)\n"
    )
    for chrom in chrom_lst:
        statistics = chrom2hblock_statistics[chrom]
        o.write(
            "{0}\t{1[0]}\t{1[1]}\t{1[2]}\t{1[3]:,}\t{1[4]:,}\t{1[5]}\n".format(
                chrom, statistics
            )
        )
    o.close()


def dump_himut_sbs(
    chrom_lst: List[str],
    chrom2tsbs_lst: Dict[str, List[Tuple[str, int, str, str, str, int, int, int, float]]],
    phase: bool,
    header: str,
    out_file: str
) -> None:

    if not out_file.endswith(".vcf"):
        print("VCF file must have .vcf suffix")
        himut.util.exit()
        
    o = open(out_file, "w")
    p = open(out_file.replace(".vcf", ".single_molecule_mutations.vcf"), "w")
    o.write("{}\n".format(header))
    p.write("{}\n".format(header))
    for chrom in chrom_lst:
        for (chrom, pos, ref, alt, annot, bq, total_count, ref_count, alt_count, vaf, phase_set) in chrom2tsbs_lst[chrom]:
            if phase:
                o.write(
                    "{}\t{}\t.\t{}\t{}\t.\t{}\t.\tGT:BQ:DP:AD:VAF:PS\t./.:{}:{:0.0f}:{:0.0f},{:0.0f}:{:.2f}:{}\n".format(
                        chrom, pos, ref, alt, annot, bq, total_count, ref_count, alt_count, vaf, phase_set
                    )
                )
                if alt_count == 1:
                    p.write(
                        "{}\t{}\t.\t{}\t{}\t.\t{}\t.\tGT:BQ:DP:AD:VAF:PS\t./.:{}:{:0.0f}:{:0.0f},{:0.0f}:{:.2f}:{}\n".format(
                            chrom, pos, ref, alt, annot, bq, total_count, ref_count, alt_count, vaf, phase_set
                        )
                    )
            else:
                o.write(
                    "{}\t{}\t.\t{}\t{}\t.\t{}\t.\tGT:BQ:DP:AD:VAF\t./.:{}:{:0.0f}:{:0.0f},{:0.0f}:{:.2f}\n".format(
                        chrom, pos, ref, alt, annot, bq, total_count, ref_count, alt_count, vaf
                    )
                )
                if alt_count == 1:
                    p.write(
                        "{}\t{}\t.\t{}\t{}\t.\t{}\t.\tGT:BQ:DP:AD:VAF\t./.:{}:{:0.0f}:{:0.0f},{:0.0f}:{:.2f}\n".format(
                            chrom, pos, ref, alt, annot, bq, total_count, ref_count, alt_count, vaf
                        )
                    )
    o.close()
    p.close()
    
    
def dump_himut_dbs(
    chrom_lst: List[str],
    chrom2tdbs_lst: Dict[str, List[Tuple[str, int, str, str, int, int, int, float]]],
    phase: bool,
    header: str,
    out_file: str
) -> None:

    o = open(out_file.replace(".vcf", ".dbs.vcf"), "w")
    o.write("{}\n".format(header))
    for chrom in chrom_lst:
        for (chrom, pos, ref, alt, bq, total_count, ref_count, alt_count, vaf, phase_set) in chrom2tdbs_lst[chrom]:
            if phase: 
                o.write(
                    "{}\t{}\t.\t{}\t{}\t.\tPASS\t.\tGT:BQ:DP:AD:VAF:PS\t./.:{}:{:0.0f}:{:0.0f},{:0.0f}:{:.2f}:{}\n".format(
                        chrom, pos, ref, alt, bq, total_count, ref_count, alt_count, vaf, phase_set
                    )
                )
            else:
                o.write(
                    "{}\t{}\t.\t{}\t{}\t.\tPASS\t.\tGT:BQ:DP:AD:VAF\t./.:{}:{:0.0f}:{:0.0f},{:0.0f}:{:.2f}\n".format(
                        chrom, pos, ref, alt, bq, total_count, ref_count, alt_count, vaf
                    )
                )
    o.close()


def dump_himut_statistics(
    chrom_lst: List[str], 
    genome_stat_hsh: Dict[str, List[int]], 
    stats_file: str
) -> None:
 
    ncol = len(chrom_lst)
    nrow = len(ROW_NAMES)
    dt = np.zeros((nrow, ncol))
    for idx, chrom in enumerate(chrom_lst):
        for jdx, count in enumerate(genome_stat_hsh[chrom]): 
            dt[jdx][idx] = count
    
    o = open(stats_file, "w")
    genome_lst = chrom_lst + ["total"] 
    o.write("{:30}{}\n".format("", "\t".join(genome_lst)))
    for kdx in range(nrow):
        row_sum =  str(int(np.sum(dt[kdx])))
        stats = "\t".join([str(int(stat)) for stat in dt[kdx].tolist()] + [row_sum])
        o.write("{:30}{}\n".format(ROW_NAMES[kdx], stats))
    o.close()
