#!/usr/bin/env python3
__version__ = "1.0.0"
__author__ = "Sangjin Lee"

# modules
import himut.util
import himut.caller
import himut.vcflib
import himut.mutlib
import himut.phaselib
import himut.normcounts
from himut.parse_args import parse_args


def main():
    parser, options = parse_args(program_version=__version__)
    if options.sub == "call":  # call somatic substitutions
        himut.util.check_num_threads(options.threads)
        himut.caller.call_somatic_substitutions(
            options.bam,  # BAM file
            options.ref,  # reference FASTA file
            options.vcf,  # VCF file: germline mutations
            options.phased_vcf,  # VCF file: phased germline mutations
            options.region,  # target contigs/scaffolds/chromosomes
            options.region_list,  # target contigs/scaffolds/chromosomes fofn
            options.min_mapq,  # int: 0-60
            options.min_sequence_identity,  # float: blast sequence identity
            options.min_hq_base_proportion,  # float: proportion of BQ=93 bases
            options.min_alignment_proportion,  # float: query alignment length / query read length
            options.min_gq,  # minimum germline genotype quality score: int
            options.min_bq,  # minimum base quality score: int
            options.min_trim,  # float: 0.01 - 0.1
            options.mismatch_window,  # mismatch window size
            options.max_mismatch_count,  # maximum number of mismatches within a window
            options.min_ref_count,  # number of reads supporting the reference base
            options.min_alt_count,  # number of reads supporting the alterantive base
            options.min_hap_count,  # number of reads supporting h0 and h1 haplotype
            options.common_snps,  # common snps
            options.panel_of_normals,  # panel of normals
            options.somatic_snv_prior,
            options.germline_snv_prior,
            options.germline_indel_prior,
            options.threads,  # maxminum number of threads
            options.phase,  # bool
            options.non_human_sample,  # bool
            options.reference_sample,  # bool
            options.create_panel_of_normal,  # bool
            __version__,  # str
            options.out,  # output # himut vcf file
        )
    elif options.sub == "phase":  # returns phased hetsnps
        himut.phaselib.get_chrom_hblock(
            options.bam,
            options.vcf,
            options.region,
            options.region_list,
            options.min_bq,
            options.min_mapq,
            options.min_p_value,
            options.min_phase_proportion,
            options.threads,
            __version__,
            options.out,
        )
    # elif options.sub == "dbs78": # returns dbs78 counts
    #     sample = himut.vcflib.get_sample(options.input)
    #     himut.mutlib.dump_dbs78_counts(options.input, options.region, options.region_list, options.output)
    #     himut.mutlib.dump_dbs78_plt(options.output, sample, options.output.replace(".tsv", ".pdf"))
    elif options.sub == "sbs96":  # returns sbs96 counts
        sample = himut.vcflib.get_sample(options.input)
        himut.util.check_mutpatterns_input_exists(
            options.input,
            options.ref,
            options.region,
            options.region_list,
            options.output,
        )
        himut.mutlib.dump_sbs96_counts(
            options.input,
            options.ref,
            options.region,
            options.region_list,
            options.output,
        )
        himut.mutlib.dump_sbs96_plt(
            options.output, sample, "{}.pdf".format(options.output)
        )
    elif options.sub == "normcounts":  # returns normalised sbs96 counts
        himut.util.check_num_threads(options.threads)
        himut.normcounts.get_normcounts(
            options.bam,
            options.ref,
            options.sbs,
            options.vcf,
            options.phased_vcf,
            options.min_mapq,
            options.min_sequence_identity,
            options.min_hq_base_proportion,
            options.min_alignment_proportion,
            options.min_gq,
            options.min_bq,
            options.min_trim,
            options.mismatch_window,
            options.max_mismatch_count,
            options.min_ref_count,
            options.min_alt_count,
            options.common_snps,
            options.panel_of_normals,
            options.somatic_snv_prior,
            options.germline_snv_prior,
            options.germline_indel_prior,
            options.threads,
            options.phase,
            options.non_human_sample,
            options.reference_sample,
            options.out,
        )
    else:
        print("The subcommand does not exist!\n")
        parser.print_help()
        parser.exit()


if __name__ == "__main__":
    main()
    himut.util.exit()
