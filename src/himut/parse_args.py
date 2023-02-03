# modules
import sys
import warnings
import argparse


def make_wide(formatter, w=120, h=36):
    """Return a wider HelpFormatter, if possible."""
    try:
        # https://stackoverflow.com/a/5464440
        # beware: "Only the name of this class is considered a public API."
        kwargs = {"width": w, "max_help_position": h}
        formatter(None, **kwargs)
        return lambda prog: formatter(prog, **kwargs)
    except TypeError:
        warnings.warn("argparse help formatter failed, falling back.")
        return formatter


def parse_args(program_version, arguments=sys.argv[1:]):
    # main_arguments
    parser = argparse.ArgumentParser(
        add_help=True,
        formatter_class=make_wide(argparse.ArgumentDefaultsHelpFormatter),
        description="himut identifies high-confidence single molecule somatic single-base substitutions from PacBio CCS reads",
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="%(prog)s {version}".format(version=program_version),
    )
    # subcommands: init
    subparsers = parser.add_subparsers(dest="sub", metavar="")

    # subcommands: call
    parser_call = subparsers.add_parser(
        "call",
        formatter_class=make_wide(argparse.ArgumentDefaultsHelpFormatter, w=180, h=60),
        help="detects somatic mutations from circular consensus seuqence (CCS) reads",
    )
    parser_call.add_argument(
        "-i",
        "--bam",
        type=str,
        required=True,
        help="minimap2 (parameters: -ax map-hifi --cs=short) aligned SAM/BAM files",
    )
    parser_call.add_argument(
        "--ref",
        type=str,
        required=False,
        help="reference genome FASTA file",
    )
    parser_call.add_argument(
        "--vcf",
        type=str,
        required=False,
        help="deepvariant VCF file with germline mutations",
    )
    parser_call.add_argument(
        "--phased_vcf",
        type=str,
        required=False,
        help="phased deepvariant VCF file",
    )
    parser_call.add_argument(
        "--common_snps",
        type=str,
        required=False,
        help="1000G common SNPs VCF file",
    )
    parser_call.add_argument(
        "--panel_of_normals",
        type=str,
        required=False,
        help="panel of normal VCF file",
    )
    parser_call.add_argument(
        "--region",
        type=str,
        required=False,
        help="target chromosome",
    )
    parser_call.add_argument(
        "--region_list",
        type=str,
        required=False,
        help="list of target chromosomes separated by new line",
    )
    parser_call.add_argument(
        "--min_qv",
        type=int,
        default=30,
        required=False,
        help="minimum average read accuracy",
    )
    parser_call.add_argument(
        "--min_mapq",
        type=int,
        default=60,
        required=False,
        help="minimum mapping quality score",
    )
    parser_call.add_argument(
        "--min_sequence_identity",
        type=float,
        default=0.99,
        required=False,
        help="minimum sequence identity threshold",
    )
    parser_call.add_argument(
        "--min_gq",
        type=int,
        default=40,
        required=False,
        help="minimum germline genotype quality (GQ) score ",
    )
    parser_call.add_argument(
        "--min_bq",
        type=int,
        default=93,
        required=False,
        help="minimum base quality score threshold",
    )
    parser_call.add_argument(
        "--min_ref_count",
        type=int,
        default=3,
        required=False,
        help="minimum reference allele depth at single base substitution site",
    )
    parser_call.add_argument(
        "--min_alt_count",
        type=int,
        default=1,
        required=False,
        help="minimum alternative allele depth at single base substitution site",
    )
    parser_call.add_argument(
        "--min_hap_count",
        type=int,
        default=3,
        required=False,
        help="minimum h0 and h1 haplotype count",
    )
    parser_call.add_argument(
        "--min_trim",
        type=float,
        default=0.01,
        required=False,
        help="minimum proportion of bases to be trimmed from the start and end of the read",
    )
    parser_call.add_argument(
        "--mismatch_window",
        type=int,
        default=20,
        required=False,
        help="mismatch window size",
    )
    parser_call.add_argument(
        "--max_mismatch_count",
        type=int,
        default=0,
        required=False,
        help="maximum number of mismatches within the mismatch window",
    )
    parser_call.add_argument(
        "--somatic_snv_prior",
        type=float,
        default=1/(10**6),
        required=False,
        help="somatic single-base-substitution prior",
    )
    parser_call.add_argument(
        "--germline_snv_prior",
        type=float,
        default=1/(10**3),
        required=False,
        help="germline snv prior",
    )
    parser_call.add_argument(
        "--germline_indel_prior",
        type=float,
        default=1/(10**4),
        required=False,
        help="germline indel prior",
    )
    parser_call.add_argument(
        "-t",
        "--threads",
        type=int,
        default=1,
        required=False,
        help="number of threads to use",
    )
    parser_call.add_argument(
        "--phase",
        required=False,
        action="store_true",
        help="return phased somatic substitutions",
    )
    parser_call.add_argument(
        "--non_human_sample",
        required=False,
        action="store_true",
        help="human (default) or non-human sample",
    )
    parser_call.add_argument(
        "--reference_sample",
        required=False,
        action="store_true",
        help="reads from the sample has been used to create the reference genome",
    )
    parser_call.add_argument(
        "--create_panel_of_normal",
        required=False,
        action="store_true",
        help="call somatic mutations with relaxed parameters for panel of normal preparation",
    )
    parser_call.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="VCF file to write the somatic substitutions",
    )
    # subcommands: sbs96
    parser_sbs96 = subparsers.add_parser(
        "sbs96",
        formatter_class=make_wide(argparse.ArgumentDefaultsHelpFormatter, w=180, h=60),
        help="returns SBS96 counts and barplot",
    )
    parser_sbs96.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="himut VCF file to read somatic single base substitutions",
    )
    parser_sbs96.add_argument(
        "--ref", type=str, required=True, help="reference FASTA file"
    )
    parser_sbs96.add_argument(
        "--region",
        type=str,
        required=False,
        help="target chromosome",
    )
    parser_sbs96.add_argument(
        "--region_list",
        type=str,
        required=False,
        help="list of target chromosomes separated by new line",
    )
    parser_sbs96.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="file to return SBS96 counts (.tsv suffix)",
    )
    # subcommands: phase
    parser_phase = subparsers.add_parser(
        "phase",
        formatter_class=make_wide(argparse.ArgumentDefaultsHelpFormatter, w=180, h=60),
        help="returns phased hetsnps",
    )
    parser_phase.add_argument(
        "-i",
        "--bam",
        type=str,
        required=True,
        help="minimap2 (parameters: -ax map-hifi --cs=short) aligned SAM/BAM files",
    )
    parser_phase.add_argument(
        "--vcf",
        type=str,
        required=True,
        help="deepvariant VCF file with germline mutations",
    )
    parser_phase.add_argument(
        "--region",
        type=str,
        required=False,
        help="target chromosome",
    )
    parser_phase.add_argument(
        "--region_list",
        type=str,
        required=False,
        help="list of target chromosomes separated by new line",
    )
    parser_phase.add_argument(
        "--min_bq",
        type=int,
        default=20,
        required=False,
        help="minimum base quality score threshold",
    )
    parser_phase.add_argument(
        "--min_mapq",
        type=int,
        default=20,
        required=False,
        help="minimum mapping quality score",
    )
    parser_phase.add_argument(
        "--min_p_value",
        type=float,
        default=0.0001,
        required=False,
        help="minimum proportion of phase consistent edges",
    )
    parser_phase.add_argument(
        "--min_phase_proportion",
        type=float,
        default=0.2,
        required=False,
        help="minimum proportion of phase consistent edges",
    )
    parser_phase.add_argument(
        "-t", 
        "--threads", 
        type=int, 
        default=1, 
        required=False, 
        help="number of threads"
    )
    parser_phase.add_argument(
        "-o", 
        "--output", 
        type=str, 
        required=True, 
        help="VCF file to write phased hetsnps"
    )
    # subcommands: burden
    parser_burden = subparsers.add_parser(
        "burden",
        formatter_class=make_wide(argparse.ArgumentDefaultsHelpFormatter, w=180, h=60),
        help="calculates mutation burden per cell",
    )
    parser_burden.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="normalised SBS96 counts",
    )
    parser_burden.add_argument(
        "--ref",
        type=str,
        required=True,
        help="reference FASTA file",
    )
    parser_burden.add_argument(
        "--tri",
        type=str,
        required=True,
        help="reference trinucleotide sequence context",
    )
    parser_burden.add_argument(
        "--region_list",
        type=str,
        required=False,
        help="list of autosomes and sex chromosomes separated by new line",
    )
    parser_burden.add_argument(
        "-t",
        "--threads",
        type=int,
        default=1,
        required=False,
        help="number of threads to use",
    )
    parser_burden.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="file to return normalised ation burden",
    )
    # subcommands: tricount
    parser_tricount = subparsers.add_parser(
        "tricount",
        formatter_class=make_wide(argparse.ArgumentDefaultsHelpFormatter, w=180, h=60),
        help="calculates and returns reference trinucletide context counts",
    )
    parser_tricount.add_argument(
        "-i",
        "--ref",
        type=str,
        required=True,
        help="reference FASTA file",
    )
    parser_tricount.add_argument(
        "--region",
        type=str,
        required=False,
        help="target chromosome",
    )
    parser_tricount.add_argument(
        "--region_list",
        type=str,
        required=False,
        help="list of target chromosomes separated by new line",
    )
    parser_tricount.add_argument(
        "-t",
        "--threads",
        type=int,
        default=1,
        required=False,
        help="number of threads to use",
    )
    parser_tricount.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="file to return normalised mutation burden",
    )
    # subcommands: normcounts
    parser_normcounts = subparsers.add_parser(
        "normcounts",
        formatter_class=make_wide(argparse.ArgumentDefaultsHelpFormatter, w=180, h=60),
        help="normalises SBS96 mutation counts based on genome and read trinucletide context counts",
    )
    parser_normcounts.add_argument(
        "--bam",
        type=str,
        required=True,
        help="minimap2 (parameters: -ax map-hifi --cs=short) aligned SAM/BAM files",
    )
    parser_normcounts.add_argument(
        "--ref",
        type=str,
        required=True,
        help="reference FASTA file",
    )
    parser_normcounts.add_argument(
        "--sbs",
        type=str,
        required=True,
        help="himut VCF file with somatic single-base substitutions",
    )
    parser_normcounts.add_argument(
        "--vcf",
        type=str,
        required=False,
        help="deepvariant VCF file with germline mutations",
    )
    parser_normcounts.add_argument(
        "--phased_vcf",
        type=str,
        required=False,
        help="phased deepvariant VCF file",
    )
    parser_normcounts.add_argument(
        "--common_snps",
        type=str,
        required=False,
        help="1000G common SNPs VCF file",
    )
    parser_normcounts.add_argument(
        "--panel_of_normals",
        type=str,
        required=False,
        help="panel of normal VCF file",
    )
    parser_normcounts.add_argument(
        "--region",
        type=str,
        required=False,
        help="target chromosome",
    )
    parser_normcounts.add_argument(
        "--region_list",
        type=str,
        required=False,
        help="list of target chromosomes separated by new line",
    )
    parser_normcounts.add_argument(
        "--min_qv",
        type=int,
        default=30,
        required=False,
        help="minimum average read accuracy",
    )
    parser_normcounts.add_argument(
        "--min_mapq",
        type=int,
        default=60,
        required=False,
        help="minimum mapping quality score",
    )
    parser_normcounts.add_argument(
        "--min_sequence_identity",
        type=float,
        default=0.99,
        required=False,
        help="minimum sequence identity threshold",
    )
    parser_normcounts.add_argument(
        "--min_gq",
        type=int,
        default=40,
        required=False,
        help="minimum germline genotype quality (GQ) score",
    )
    parser_normcounts.add_argument(
        "--min_bq",
        type=int,
        default=93,
        required=False,
        help="minimum base quality score threshold",
    )
    parser_normcounts.add_argument(
        "--min_ref_count",
        type=int,
        default=3,
        required=False,
        help="minimum reference allele depth",
    )
    parser_normcounts.add_argument(
        "--min_alt_count",
        type=int,
        default=1,
        required=False,
        help="minimum alternative allele depth",
    )
    parser_normcounts.add_argument(
        "--min_hap_count",
        type=int,
        default=3,
        required=False,
        help="minimum h0 and h1 haplotype count",
    )
    parser_normcounts.add_argument(
        "--min_trim",
        type=float,
        default=0.01,
        required=False,
        help="minimum proportion of bases to be trimmed from the start and end of the read",
    )
    parser_normcounts.add_argument(
        "--mismatch_window",
        type=int,
        default=20,
        required=False,
        help="mismatch window size",
    )
    parser_normcounts.add_argument(
        "--max_mismatch_count",
        type=int,
        default=0,
        required=False,
        help="maximum number of mismatches within the mismatch window",
    )
    parser_normcounts.add_argument(
        "--somatic_snv_prior",
        type=float,
        default=1/(10**6),
        required=False,
        help="somatic single-base-substitution prior",
    )
    parser_normcounts.add_argument(
        "--germline_snv_prior",
        type=float,
        default=1/(10**3),
        required=False,
        help="germline snv prior",
    )
    parser_normcounts.add_argument(
        "--germline_indel_prior",
        type=float,
        default=1/(10**4),
        required=False,
        help="germline indel prior",
    )
    parser_normcounts.add_argument(
        "-t",
        "--threads",
        type=int,
        default=1,
        required=False,
        help="number of threads to use",
    )
    parser_normcounts.add_argument(
        "--phase",
        required=False,
        action="store_true",
        help="return phased somatic substitutions",
    )
    parser_normcounts.add_argument(
        "--non_human_sample",
        required=False,
        action="store_true",
        help="human or non_human_sample",
    )
    parser_normcounts.add_argument(
        "--reference_sample",
        required=False,
        action="store_true",
        help="reads from the sample has been used to create the reference genome",
    )
    parser_normcounts.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="file to return normalised SBS96 counts",
    )
    if len(arguments) == 0:
        parser.print_help()
        parser.exit()
    else:
        return parser, parser.parse_args(arguments)
