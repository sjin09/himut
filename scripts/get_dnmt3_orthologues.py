#!/usr/bin/env python

import argparse
import os
from pathlib import Path
from typing import Dict
import sys

import natsort

DNMT3_GENES = ["DNMT3A", "DNMT3B", "DNMT3L"]
HSAPIENS_ENSEMBLE_TO_GENE_LOOKUP = {
    "ENSG00000119772": "DNMT3A",
    "ENSG00000088305": "DNMT3B",
    "ENSG00000142182": "DNMT3L"
}
HSAPIENS_DNMT3_ENSEMBLE_SET = set(["ENSG00000119772", "ENSG00000088305", "ENSG00000142182"])


def parse_args(args):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        required=True,
        help="sample name fofn",
    )
    parser.add_argument(
        "-d",
        "--dir",
        type=Path,
        required=True,
        help="orthoFinder results directory",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="file to write Homo sapiens DNMT3 orthologues",
    )
    args = args[1:]
    return parser.parse_args(args)


def get_orthogroup_count(orthogroup_gene_count_file_path: Path, sample: str):
    # init lookup
    sample_orthogroup_count = {}
    hsapiens_orthogroup_count = {}
    header = open(orthogroup_gene_count_file_path).readline().rstrip().split("\t")
    sample_idx = header.index(sample)
    hsapiens_idx = header.index("Homo_sapiens")
    for line in open(orthogroup_gene_count_file_path).readlines()[1:]:
        orthogroup, hsapiens_count, sample_count, total_count = line.rstrip().split("\t")
        sample_orthogroup_count[orthogroup] = int(sample_count)
        hsapiens_orthogroup_count[orthogroup] = int(hsapiens_count)
    return sample_orthogroup_count, hsapiens_orthogroup_count


def get_sample_dnmt3_orthologues(
    orthogroup_file_path: Path,
    hsapiens_orthogroup_counts: Dict[str, int],
    sample_orthogroup_counts: Dict[str, int]
):
    sample_dnmt3_orthologues = {dnmt3_gene: "." for dnmt3_gene in DNMT3_GENES}
    for line in open(orthogroup_file_path).readlines()[1:]:
        fields = line.rstrip().split()
        orthogroup = fields[0]
        gene_annotations = fields[1:]
        is_dnmt3_gene = False
        for gene_annotation in gene_annotations:
            gene_annotation_prefix = gene_annotation.split(".")[0]
            if gene_annotation_prefix in HSAPIENS_DNMT3_ENSEMBLE_SET:
                is_dnmt3_gene =True
                break
        if is_dnmt3_gene:
            gene = HSAPIENS_ENSEMBLE_TO_GENE_LOOKUP[gene_annotation_prefix]
            human_orthogroup_count = hsapiens_orthogroup_counts[orthogroup]
            sample_orthogroup_count = sample_orthogroup_counts[orthogroup]
            human_idx_end = 1 + human_orthogroup_count
            sample_idx_end = human_idx_end + sample_orthogroup_count 
            human_ensemble_genes = fields[1:human_idx_end]
            sample_ensemble_genes = fields[human_idx_end:sample_idx_end]
            sample_ensemble_genes = [sample_ensemble_gene.replace(",", "") for sample_ensemble_gene in sample_ensemble_genes]
            sample_dnmt3_orthologues[gene] = ",".join(sample_ensemble_genes)
    return sample_dnmt3_orthologues


def get_dnmt3_orthologues(
    input_file_path: Path, 
    orthofinder_dir_path: Path,
    output_file_path: Path
):

    with open(output_file_path, "w") as outfile:
        print("SAMPLE", "DNMT3A", "DNMT3B", "DNMT3L", sep="\t", file=outfile)
        for line in open(input_file_path).readlines():
            # init path
            sample = line.rstrip()
            sample_orthofinder_dir_path = os.path.join(orthofinder_dir_path, "Hsapiens.{}".format(sample))
            orthogroup_gene_count_file_path = os.path.join(sample_orthofinder_dir_path, "Orthogroups.GeneCount.tsv")
            orthogroup_file_path = os.path.join(sample_orthofinder_dir_path, "Homo_sapiens__v__{}.tsv".format(sample))
            sample_orthogroup_counts, hsapiens_orthogroup_counts = get_orthogroup_count(orthogroup_gene_count_file_path, sample)
            sample_dnmt3_orthologues = get_sample_dnmt3_orthologues(orthogroup_file_path, hsapiens_orthogroup_counts, sample_orthogroup_counts)
            print(
                sample, 
                sample_dnmt3_orthologues["DNMT3A"], 
                sample_dnmt3_orthologues["DNMT3B"], 
                sample_dnmt3_orthologues["DNMT3L"],
                sep="\t",
                file=outfile
            )

def main():
    options = parse_args(sys.argv)
    get_dnmt3_orthologues(options.input, options.dir, options.output)
    sys.exit(0)


if __name__ == "__main__":
    main()


