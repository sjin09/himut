# modules
import argparse
import sys
import warnings
from pathlib import Path


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
        description="""
        himut identifies high-confidence single molecule somatic single-base substitutions from PacBio CCS reads
        """
    )
    parser.add_argument("--force", action="store_true", version="force overwriting of exisiting files")
    parser.add_argument("--region", type=str, required=False, help="target chromosome")
    parser.add_argument(
        "--region-list", type=Path, required=False, help="a file with list of target chromosomes separated by new line"
    )
    parser.add_argument("-t", "--threads", type=int, default=1, required=False, help="number of threads to use")
    parser.add_argument(
        "-v", "--version", action="version", version="%(prog)s {version}".format(version=program_version)
    )
    # subcommands: initialize
    subparsers = parser.add_subparsers(dest="subcommand", metavar="")
    # subcommands: call
    parser_call = subparsers.add_parser(
        "CallSomaticMutations",
        formatter_class=make_wide(argparse.ArgumentDefaultsHelpFormatter, w=180, h=60),
        help="detects somatic mutations from circular consensus seuqence (CCS) reads"
    )
    parser_call.add_argument(
        "-i", "--bam", type=Path, required=True, help="minimap2 BAM file",
    )
    parser_call.add_argument("--ref-fasta", type=Path, required=False, help="reference genome FASTA file")
    parser_call.add_argument("--germline-vcf", type=Path, required=False, help="VCF file with germline mutations")
    parser_call.add_argument(
        "--germline-phased-vcf", type=Path, required=False, help="VCF file with phased germline mutations"
    )
    parser_call.add_argument("--common-snps", type=Path, required=False, help="VCF file with 1000G common SNPs")
    parser_call.add_argument("--panel-of-normals", type=Path, required=False, help="panel of normal VCF file")
    parser_call.add_argument(
        "--min-rq", type=int, default=30, required=False,
        help="minimum read quality (rq is the mean of the read base quality scores)"
    )
    parser_call.add_argument("--min-mq", type=int, default=60, required=False, help="minimum mapping quality score")
    parser_call.add_argument(
        "--min-sequence-identity", type=float, default=0.99, required=False, help="minimum sequence identity"
    )
    parser_call.add_argument(
        "--min-gq", type=int, default=20, required=False, help="minimum germline genotype quality (GQ) score"
    )
    parser_call.add_argument("--min-bq", type=int, default=93, required=False, help="minimum base quality score")
    parser_call.add_argument(
        "--min-ref-count", type=int, default=3, required=False, help="minimum reference allele depth"
    )
    parser_call.add_argument(
        "--min-alt-count", type=int, default=1, required=False, help="minimum alternative allele depth"
    )
    parser_call.add_argument(
        "--min-hap-count", type=int, default=3, required=False, help="minimum h1 and h2 haplotype count"
    )
    parser_call.add_argument(
        "--min-trim", type=float, default=0.01, required=False,
        help="minimum proportion of bases to be trimmed from the start and end of the read"
    )
    parser_call.add_argument(
        "--max-mismatch-count", type=int, default=0, required=False,
        help="maximum number of mismatches within the mismatch window"
    )
    parser_call.add_argument(
        "--mismatch-window-size", type=int, default=20, required=False, help="mismatch window size"
    )
    parser_call.add_argument(
        "--somatic-sbs-prior", type=float, default=1/(10**6), required=False,
        help="human somatic single-base-substitution prior"
    )
    parser_call.add_argument(
        "--germline-snp-prior", type=float, default=1/(10**3), required=False,
        help="human germline single nucleotide polymorphism prior"
    )
    parser_call.add_argument(
        "--phase", required=False, action="store_true", help="return phased somatic substitutions"
    )
    parser_call.add_argument(
        "--non-human-sample", required=False, action="store_true", help="human (default) or non-human sample"
    )
    parser_call.add_argument(
        "--reference-sample", required=False, action="store_true",
        help="sample reads have been used to create the reference genome"
    )
    parser_call.add_argument(
        "--create-panel-of-normal", required=False, action="store_true",
        help="call somatic mutations with relaxed parameters to prepare panel of normal VCF file"
    )
    parser_call.add_argument(
        "-o", "--output", type=Path, required=True, help="VCF file to write the somatic single-base-substitutions"
    )
    # subcommands: sbs52
    parser_sbs52 = subparsers.add_parser(
        "SBS52",
        help="returns SBS52 counts and barplot",
        formatter_class=make_wide(argparse.ArgumentDefaultsHelpFormatter, w=180, h=60)
    )
    parser_sbs52.add_argument(
        "-i", "--vcf", type=Path, required=True, help="himut VCF file to read somatic single-base-substitutions"
    )
    parser_sbs52.add_argument("--ref-fasta", type=Path, required=True, help="reference FASTA file")
    # subcommands: sbs96
    parser_sbs96 = subparsers.add_parser(
        "SBS96",
        help="returns SBS96 counts and barplot",
        formatter_class=make_wide(argparse.ArgumentDefaultsHelpFormatter, w=180, h=60)
    )
    parser_sbs96.add_argument(
        "-i", "--vcf", type=Path, required=True, help="himut VCF file to read somatic single-base-substitutions"
    )
    parser_sbs96.add_argument("--ref-fasta", type=Path, required=True, help="reference FASTA file")
    # subcommands: sbs1536
    parser_sbs1536 = subparsers.add_parser(
        "SBS1536",
        help="returns SBS1536 counts and barplot",
        formatter_class=make_wide(argparse.ArgumentDefaultsHelpFormatter, w=180, h=60)
    )
    parser_sbs1536.add_argument(
        "-i", "--vcf", type=Path, required=True, help="himut VCF file to read somatic single-base-substitutions"
    )
    parser_sbs1536.add_argument("--ref-fasta", type=Path, required=True, help="reference FASTA file")
    # subcommands: phase
    parser_phase = subparsers.add_parser(
        "Phase",
        help="phase heterozygous SNPs and return a phased VCF file",
        formatter_class=make_wide(argparse.ArgumentDefaultsHelpFormatter, w=180, h=60)
    )
    parser_phase.add_argument("-i", "--bam", type=Path, required=True, help="minimap2 BAM files")
    parser_phase.add_argument("--vcf", type=Path, required=True, help="VCF file with germline mutations")
    parser_phase.add_argument("--min-bq", type=int, default=20, required=False, help="minimum base quality score")
    parser_phase.add_argument("--min-mq", type=int, default=20, required=False, help="minimum mapping quality score")
    parser_phase.add_argument(
        "--min-p-value", type=float, default=0.0001, required=False,
        help="minimum proportion of phase consistent edges"
    )
    parser_phase.add_argument(
        "--min-phase-proportion", type=float, default=0.2, required=False,
        help="minimum proportion of phase consistent edges"
    )
    parser_phase.add_argument(
        "-o", "--output", type=Path, required=True, help="VCF file to write phased heterozygous SNPs"
    )
    # subcommands: burden
    parser_burden = subparsers.add_parser(
        "MutationBurden",
        help="calculates mutation burden per cell",
        formatter_class=make_wide(argparse.ArgumentDefaultsHelpFormatter, w=180, h=60)
    )
    parser_burden.add_argument("-i", "--input", type=Path, required=True, help="file to read normalised SBS96 counts")
    parser_burden.add_argument("--ref-fasta", type=Path, required=False, help="reference FASTA file")
    # subcommands: tricount
    parser_tricount = subparsers.add_parser(
        "CountTrinucleotides",
        help="counts and returns reference trinucletide context counts",
        formatter_class=make_wide(argparse.ArgumentDefaultsHelpFormatter, w=180, h=60)
    )
    parser_tricount.add_argument("-i", "--ref-fasta", type=Path, required=True, help="reference FASTA file")
    # subcommands: normcounts
    parser_normcounts = subparsers.add_parser(
        "NormaliseSomaticSubstitutionCounts",
        formatter_class=make_wide(argparse.ArgumentDefaultsHelpFormatter, w=180, h=60),
        help="returns normalised SBS96 mutation counts based on genome and read trinucletide context counts"
    )
    parser_normcounts.add_argument("-i", "--bam", type=str, required=True, help="minimap2 BAM file")
    parser_normcounts.add_argument("--ref-fasta", type=str, required=True, help="reference FASTA file")
    parser_normcounts.add_argument(
        "--somatic-vcf", type=Path, required=True, help="VCF file with somatic single-base substitutions"
    )
    parser_normcounts.add_argument(
        "--germline-vcf", type=Path, required=False, help="VCF file with germline mutations"
    )
    parser_normcounts.add_argument(
        "--germline-phased-vcf",
        type=Path, required=False, help="VCF file with phased germline mutations",
    )
    parser_normcounts.add_argument("--common-snps", type=Path, required=False, help="VCF file with 1000G common SNPs")
    parser_normcounts.add_argument("--panel-of-normals", type=Path, required=False, help="panel of normal VCF file")
    parser_normcounts.add_argument(
        "--min-rq", type=int, default=30, required=False,
        help="minimum read quality (rq is the mean of the read base quality scores)"
    )
    parser_normcounts.add_argument(
        "--min-mq", type=int, default=60, required=False, help="minimum mapping quality score"
    )
    parser_normcounts.add_argument(
        "--min-sequence-identity", type=float, default=0.99, required=False, help="minimum sequence identity"
    )
    parser_normcounts.add_argument(
        "--min-gq", type=int, default=20, required=False, help="minimum germline genotype quality (GQ) score"
    )
    parser_normcounts.add_argument("--min-bq", type=int, default=93, required=False, help="minimum base quality score")
    parser_normcounts.add_argument(
        "--min-ref-count", type=int, default=3, required=False, help="minimum reference allele depth"
    )
    parser_normcounts.add_argument(
        "--min-alt-count", type=int, default=1, required=False, help="minimum alternative allele depth"
    )
    parser_normcounts.add_argument(
        "--min-hap-count", type=int, default=3, required=False, help="minimum h1 and h2 haplotype count"
    )
    parser_normcounts.add_argument(
        "--min-trim", type=float, default=0.01, required=False,
        help="minimum proportion of bases to be trimmed from the start and end of the read"
    )
    parser_normcounts.add_argument(
        "--mismatch-window", type=int, default=20, required=False, help="mismatch window size"
    )
    parser_normcounts.add_argument(
        "--max-mismatch-count", type=int, default=0, required=False,
        help="maximum number of mismatches within the mismatch window"
    )
    parser_normcounts.add_argument(
        "--somatic-sbs-prior", type=float, default=1/(10**6), required=False,
        help="human somatic single-base-substitution prior"
    )
    parser_normcounts.add_argument(
        "--germline-snp-prior", type=float, default=1/(10**3), required=False,
        help="human germline single-nucleotide-polymorphism prior"
    )
    parser_normcounts.add_argument(
        "--phase", required=False, action="store_true",
        help="return normalised SBS96 counts for phased somatic single-base-substitutions",
    )
    parser_normcounts.add_argument(
        "--non-human-sample", required=False, action="store_true", help="human (default) or non_human_sample",
    )
    parser_normcounts.add_argument(
        "--reference-sample", required=False, action="store_true",
        help="sample reads have been used to create the reference genome",
    )
    parser_normcounts.add_argument(
        "-o", "--output", type=str, required=True, help="file to return normalised SBS96 counts",
    )
    if len(arguments) == 0:
        parser.print_help()
        parser.exit()
    else:
        return parser, parser.parse_args(arguments)
