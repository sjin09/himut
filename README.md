## <a name="started"></a>Getting started

```sh
## clone github repository and install himut
git clone https://github.com/sjin09/himut
cd himut
/bin/bash install.sh

## use miniamp2 and samtools to align, sort (and merge) PacBio CCS read alignments
minimap2 "@RG\tSM:sample" -ax map-hifi --cs ref.fa pacbio.ccs.fastq.gz | samtools sort -o aln.sorted.bam # SM tag must be provided to retrieve sample ID
samtools view -bh -F 0x900 aln.sorted.bam > aln.primary_alignments.sorted.bam # select primary alignments
samtools index aln.primary_alignments.sorted.bam
samtools merge *.primary_alignments.sorted.bam | samtools sort -o aln.primary_alignments.mergeSorted.bam -
samtools index aln.primary_alignments.mergeSorted.bam 

## use deepvariant to call germline mutations
deepvariant.simg /opt/deepvariant/bin/run_deepvariant --model_type=PACBIO --ref ref.fa --reads=aln.primary_alignments.sorted.bam --output_vcf=germline.vcf

## use himut to call somatic mutations 
himut call -i aln.primary_alignments.sorted.bam --vcf germline.vcf -o somatic.vcf --non_human_sample
himut call -i aln.primary_alignments.sorted.bam --vcf germline.vcf --common_snps ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.decomposed.normalised.common_snps.vcf.bgz --panel_of_normals b37.himut.pon.sbs.dbs.vcf.bgz -o somatic.vcf 
```

## Table of Contents

- [Getting Started](#started)
- [Users' Guide](#uguide)
  - [Installation](#install)
  - [General usage](#general)
  - [Citing himut](#cite)
  - [Limitations](#limits)

## <a name="uguide"></a>Users' Guide

himut (high-fidelity mutation) is a somatic mutation caller that leverages the base accuracy and read length of Pacific Biosciences (PacBio) CCS reads to call single molecule somatic single-base substitutions (SBS). 

A typical somatic mutation caller requires support from multiple reads to determine that a mismatch is a somatic mutation and not a sequencing error. If the read, however, perfectly represents the original DNA molecule and if a mismatch between a read alignment and reference is detected, the mismatch should theoretically reflect the biological process generating the mutation. In addition, the mutational pattern generated from the mismatches should also be concordant with the expected mutational signatures associated with the somatic mutational processes active in the sample. 

We show that PacBio CCS base accuracy is sufficiently accurate to call single molecule somatic SBS and that himut can be used to call these mutations with high confidence.

### <a name="install"></a>Installation

Download and install the latest release:

```sh
wget https://github.com/sjin09/himut/archive/refs/tags/v1.0.0.tar.gz
tar -zxvf v1.0.0.tar.gz
cd himut-1.0.0
bash install.sh
```

You can download and install the source code with the latest developments that might contain code still under development.

```sh
git clone https://github.com/sjin09/himut 
cd himut
/bin/bash install.sh 
````

Installation through Bioconda will be supported in the future.

### <a name="general"></a>General Usage

himut accepts as input BAM file with CCS read alignments and VCF file with germline mutations and returns a VCF file with somatic mutations (.vcf suffix is required). himut will also return a VCF file with single molecule somatic single base substitutions (.single_molecule_mutations.vcf) and another VCF file with single molecule double base substitutions (.dbs.vcf), which is experimental and is not for production uses. himut is also able to accept BGZIP compressed VCF files as a more memory-efficient input.

himut, currently, only accepts minimap2 generated BAM files as `--cs` tag is required for somatic mutation calling. The read alignments must be sorted, compressed and indexed using SAMtools. In addition, secondary and supplementary alignment needs to be removed from the BAM file as deepvariant can only process primary alignments.

```sh
minimap2 "@RG\tSM:sample" -ax map-hifi --cs ref.fa pacbio.ccs.fastq.gz | samtools sort -o aln.sorted.bam # SM tag must be provided to retrieve sample ID
samtools view -bh -F 0x900 aln.sorted.bam > aln.primary_alignments.sorted.bam # select primary alignments
samtools merge *.primary_alignments.sorted.bam | samtools sort -o aln.primary_alignments.mergeSorted.bam -
samtools index aln.primary_alignments.sorted.bam
samtools index aln.primary_alignments.mergeSorted.bam 
```

#### Create a Panel of Normal VCF file

Panel of Normal VCF file is used to remove substitutions arising from systematic errors. If you wish to generate a separate Panel of Normal VCF file than the one provided, please follow the following step. When the `--create_panel_of_normal` parameter is provided, himut relaxes the hard filter thresholds to maximise the number of mismatches called against the reference genome. 

```sh
himut call -i aln.primary_alignments.sorted.bam --vcf germline.vcf --common_snps ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.decomposed.normalised.common_snps.vcf.bgz --panel_of_normals b37.himut.pon.sbs.dbs.vcf.bgz -o pon.vcf --create_panel_of_normal
```

#### Haplotype phase heterozygous single nucleotide polymorphisms 

CCS reads have an average read length of ~10-20kb and can span a number of heterozygous single nucleotide polymorphisms (hetSNPs). To haplotype phase the hetSNPs, we treat each hetSNP as a node in a graph and as CCS reads are able to span multiple hetSNPs, we are able to count the number of edges between two nodes and determine whether the two nodes belong to the same haplotype. Two or more haplotype consistent nodes are connected to construct a haplotype block across the chromosome and we are able to assign CCS reads to a haplotype block based on their hetSNP composition. If the CCS read belongs to multiple haplotype blocks or if the CCS read does not have a hetSNP, CCS read is determined to be not phased. We would like to highlight that the length of haplotype blocks and the number of phased hetSNPs will be dependent on SNP density, heterozygosity and CCS read length.

```sh
## use himut to phase CCS reads and to call phased somatic mutations
himut phase --bam aln.primary_alignments.sorted.bam --vcf germline.vcf.bgz -o germline.phased.vcf 
```

hetSNPs belonging to the same phase set (PS) in the VCF file belongs to the same haplotype block.

#### call somatic mutations 

The Darwin Tree of Life Project often sequences and assembles one sample per species and hence, we are not able to prepare a VCF file with common SNPs and a Panel of Normal VCF file to remove erroneous substitutions resulting from DNA contamination or systematic errors, respectively. If CCS dataset from multiple samples are not available for the species in question, please use the `--non_human_sample` parameter to indicate the absence of a VCF file with common SNPs and a Panel of Normal VCF file.

```sh
himut call -i aln.primary_alignments.sorted.bam --vcf germline.vcf -o somatic.vcf --non_human_sample
```

We have selected common SNPs (>1% minor allele frequency) from the 1000 Genomes Project VCF file and prepared a Panel of Normal VCF file against the b37 and the GRCh38 reference genomes using publicly available CCS reads, and the VCF files are available in the Google Drive (https://drive.google.com/drive/folders/14-cSTBocwdyIemXTExuhW-mnTJb87SLT?usp=sharing). 

```sh
himut call -i aln.primary_alignments.sorted.bam --vcf germline.vcf --common_snps ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.decomposed.normalised.common_snps.vcf.bgz --panel_of_normals b37.himut.pon.sbs.dbs.vcf.bgz -o somatic.vcf 
```

#### call haplotype phased somatic mutations

If the `--phase` parameter is used, VCF file with haplotype phased hetSNPs must be provided to assign CCS reads to a haplotype block. himut will return somatic SBS only from phased CCS reads.

```sh
himut call -i aln.primary_alignments.sorted.bam --vcf germline.vcf --phased_vcf germline.phased.vcf -o somatic.phased.vcf --phase --non_human_sample
himut call -i aln.primary_alignments.sorted.bam --vcf germline.vcf --phased_vcf germline.phased.vcf --common_snps ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.decomposed.normalised.common_snps.vcf.bgz --panel_of_normals b37.himut.pon.sbs.dbs.vcf.bgz -o somatic.vcf --phase 
```

#### SBS96 counts

himut returns SBS96 counts and a bar-plot of SBS96 counts in the pyrimidine context. 

```sh
himut sbs96 -i somatic.vcf --ref ref.fa -o somatic.sbs96.tsv 
```

#### SBS96 count normalisation

To compare SBS96 counts and to extract mutational signatures from multiple samples, SBS96 counts must be normalised against the number of reference genome trinucleotide sequence contexts and the number of callable CCS trinucleotide sequence contexts from which the SBS could have been potentially called. To count the number of callable trinucleotide sequence contexts, the same set of parameters as the one that was used for somatic mutation calling must be used.

```sh
## non-human sample
himut normcounts --bam aln.primary_alignments.sorted.bam --ref ref.fa --vcf germline.vcf --sample_sbs somatic.vcf -o somatic.normcounts.tsv
himut normcounts --bam aln.primary_alignments.sorted.bam --ref ref.fa --vcf germline.vcf --phased_vcf germline.phased.vcf --sample_sbs somatic.phased.vcf -o somatic.phased.normcounts.tsv --phase
## human sample
himut normcounts --bam aln.primary_alignments.sorted.bam --ref ref.fa --vcf germline.vcf --common_snps ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.decomposed.normalised.common_snps.vcf.bgz --panel_of_normals b37.himut.pon.sbs.dbs.vcf.bgz --sample_sbs somatic.vcf -o somatic.normcounts.tsv
himut normcounts --bam aln.primary_alignments.sorted.bam --ref ref.fa --vcf germline.vcf --phased_vcf germline.phased.vcf --common_snps ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.decomposed.normalised.common_snps.vcf.bgz --panel_of_normals b37.himut.pon.sbs.dbs.vcf.bgz --sample_sbs somatic.phased.vcf -o somatic.phased.normcounts.tsv --phase
```

#### DBS78 counts (experimental)

himut also returns DBS78 counts and a bar-plot of DBS78 counts in pyrimidine context.

```sh
himut dbs78 -i somatic.dbs78.vcf -o somatic.dbs78.tsv
```

### Limitations

PacBio CCS base accuracy is the limiting factor for somatic mutation calling using himut. PacBio circular consensus sequence (pbccs) algorithm calculates the base quality score depending on the dinucleotide sequence context and the number of supporting subread bases. Our research suggests that the base quality score does not correctly reflect the true base accuracy and that the consensus sequence caller does not account for the different substitution error rates in different trinucleotide sequence contexts. We discuss how this issue can be mitigated in our manuscript.

