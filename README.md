## <a name="started"></a>Getting started

```sh
## pip installation
pip install himut

## download and install the latest release
wget https://github.com/sjin09/himut/archive/refs/tags/v1.0.2.tar.gz
tar -zxvf v1.0.3.tar.gz
cd himut-1.0.3
bash install.sh

## use miniamp2 and samtools to align, sort (and merge) PacBio CCS read alignments
minimap2 "@RG\tSM:sample" -ax map-hifi --cs ref.fa pacbio.ccs.fastq.gz | samtools sort -o aln.sorted.bam # SM tag must be provided to retrieve sample ID
samtools merge *.sorted.bam | samtools sort -o aln.mergeSorted.bam - ## if there are multiple BAM files, merge and sort the BAM files 
samtools view -bh -F 0x900 aln.sorted.bam > aln.primary_alignments.sorted.bam # select primary alignments
samtools index aln.sorted.bam
samtools index aln.primary_alignments.sorted.bam 

## germline mutation detection and haplotype phasing
deepvariant.simg /opt/deepvariant/bin/run_deepvariant --model_type=PACBIO --ref ref.fa --reads=aln.primary_alignments.sorted.bam --output_vcf=germline.vcf ## use deepvariant to call germline mutations 
bgzip -c germline.vcf > germline.vcf.bgz ## bgzip compress
tabix -p vcf germline.vcf.bgz ## tabix index
himut phase --bam aln.primary_alignments.sorted.bam --vcf germline.vcf.bgz -o germline.phased.vcf ## phsae
bgzip -c germline.phased.vcf  > germline.phased.vcf.bgz ## bgzip compress
tabix -p vcf germline.phased.vcf.bgz ## tabix index

## somatic mutation detection in human samples

### use himut to call somatic mutations in human samples
himut call -i aln.primary_alignments.sorted.bam -o somatic.vcf  # somatic mutation detection
himut call -i aln.primary_alignments.sorted.bam --phased_vcf germline.phased.vcf.bgz -o somatic.vcf --phase # haplotype phased somatic mutation detection
himut call -i aln.primary_alignments.sorted.bam --common_snps ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.decomposed.normalised.common_snps.vcf.bgz --panel_of_normals b37.himut.pon.sbs.dbs.vcf.bgz -o somatic.vcf # provide a VCF file with common SNPs or a panel of normal VCF file
himut call -i aln.primary_alignments.sorted.bam --phased_vcf germline.phased.vcf.bgz --common_snps ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.decomposed.normalised.common_snps.vcf.bgz --panel_of_normals b37.himut.pon.sbs.dbs.vcf.bgz -o somatic.vcf --phase # somatic mutation detection with the highest sensitivity 
bgzip -c somatic.vcf > somatic.vcf.bgz
tabix -p vcf somatic.vcf.bgz

### normalise SBS96 counts ## use the same parameters as that used for somatic mutation detection 
himut normcounts --bam aln.primary_alignments.sorted.bam --ref ref.fa --sbs somatic.vcf.bgz -o somatic.normcounts.tsv
himut normcounts --bam aln.primary_alignments.sorted.bam --ref ref.fa --sbs somatic.vcf.bgz --phased_vcf germline.phased.vcf.bgz -o somatic.normcounts.tsv --phase
himut normcounts --bam aln.primary_alignments.sorted.bam --ref ref.fa --sbs somatic.vcf.bgz --common_snps ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.decomposed.normalised.common_snps.vcf.bgz --panel_of_normals b37.himut.pon.sbs.dbs.vcf.bgz -o somatic.normcounts.tsv
himut normcounts --bam aln.primary_alignments.sorted.bam --ref ref.fa --sbs somatic.vcf.bgz --common_snps ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.decomposed.normalised.common_snps.vcf.bgz --panel_of_normals b37.himut.pon.sbs.dbs.vcf.bgz --phased_vcf germline.phased.vcf.bgz -o somatic.normcounts.tsv --phase

## somatic mutation detection in non-human samples

### use himut to call somatic mutations in non-human samples
himut call -i aln.primary_alignments.sorted.bam --vcf germline.vcf.bgz -o somatic.vcf --non_human_sample 
himut call -i aln.primary_alignments.sorted.bam --vcf germline.vcf.bgz --phased_vcf germline.phased.vcf.bgz -o somatic.vcf --phase --non_human_sample 
bgzip -c somatic.vcf > somatic.vcf.bgz
tabix -p vcf somatic.vcf.bgz

## normalise SBS96 counts
himut normcounts --bam aln.primary_alignments.sorted.bam --ref ref.fa --vcf germline.vcf.bgz --sbs somatic.vcf.bgz -o somatic.normcounts.tsv
himut normcounts --bam aln.primary_alignments.sorted.bam --ref ref.fa --vcf germline.vcf.bgz --phased_vcf germline.phased.vcf.bgz --sbs somatic.vcf.bgz -o somatic.normcounts.tsv
```

## Table of Contents

- [Getting Started](#started)
- [Users' Guide](#uguide)
  - [Installation](#install)
  - [minimap2 CCS read alignment](#alignment)
  - [General usage](#general)
    - [Somatic mutation detection](#som)
    - [Haplotype phased somatic mutation detection](#phasedsom)
    - [SBS96 count calculation](#SBS96) 
    - [SBS96 counts normalisation](#normcounts)
    - [Somatic mutation detection across the Tree of Life](#tol)
  - [Panel of Normal VCF file](#pon)
  - [VCF file with common SNPs](#pop)
  - [Citing himut](#cite)
  - [Limitations](#limits)

## <a name="uguide"></a>Users' Guide

himut (high-fidelity mutation) is a somatic mutation caller that leverages the base accuracy and read length of Pacific Biosciences (PacBio) CCS reads to call single molecule somatic single-base substitutions (SBS). 

A single-read mismatch bewteen a sample and the reference genome is derived from library errors, sequencing errors or somatic mutations. If sequencing error rate is lower than the somatic mutation rate, single-molecule somatic mutation detection is possible. On the other hand, if sequencing error rate is higher than the somatic mutation rate, but if the number of accumulated somatic mutations is higher than the number of sequencing errors, which is a constant, somatic mutational process associated with the sample can be directly observed from the mutational spectrum.

### <a name="install"></a>Installation

Use pip to install himut:

```sh
pip install himut
```

Download and install the latest release:

```sh
wget https://github.com/sjin09/himut/archive/refs/tags/v1.0.3.tar.gz
tar -zxvf v1.0.2.tar.gz
cd himut-1.0.2
bash install.sh
```

You can download and install the source code with the latest developments that might contain code still under development.

```sh
git clone https://github.com/sjin09/himut 
cd himut
/bin/bash install.sh 
````

Installation through Bioconda will be supported in the future.

#### <a name="alignment"></a>minimap2 CCS read alignment

Currently, himut only accepts minimap2 CCS read alignment SAM/BAM file with a `--cs` tag, like the `cigar` string provides information about the matches and mismatches between the read and the reference genome. pbmm2 and ngmlr does not provide the `--cs` tag and as a result, pbmm2 and ngmlr SAM/BAM files are incompatible with himut. The `--cs` tag is favored over the `cigar` string because it offers a more elegant representation of the differences between the read and the reference genome (source: https://lh3.github.io/2018/03/27/the-history-the-cigar-x-operator-and-the-md-tag).

```sh
minimap2 "@RG\tSM:sample" -ax map-hifi --cs ref.fa pacbio.ccs.fastq.gz | samtools sort -o aln.sorted.bam # SM tag must be provided to retrieve sample ID
samtools merge *.sorted.bam | samtools sort -o aln.mergeSorted.bam - # merge read alignments
samtools view -bh -F 0x900 aln.mergeSorted.bam > aln.primary_alignments.mergeSorted.bam # select primary alignments
samtools index aln.primary_alignments.sorted.bam
samtools index aln.primary_alignments.mergeSorted.bam 
```

### <a name="general"></a>General Usage

Currently, himut only accepts minimap2 generated BAM files, as `--cs` tag is required for somatic mutation detection. The CCS read alignments must be sorted, \(merged\), compressed and indexed using SAMtools.

#### <a name="som"></a>Somatic mutation detection

Himut accepts as input a BAM file with CCS read alignments and returns a VCF file with somatic mutations. 

```sh
himut call -i sample.primary_alignments.sorted.bam -o sample.somatic.vcf  
```

himut can optionally accept a VCF file with common SNPs and/or a panel of normal VCF file to filter substitutions derived from genomic DNA contamination and/or from non-canonical error modes, respectively.

```sh
himut call -i aln.primary_alignments.sorted.bam --common_snps ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.decomposed.normalised.common_snps.vcf.bgz --panel_of_normals b37.himut.pon.sbs.dbs.vcf.bgz -o somatic.vcf # provide a VCF file with common SNPs or a panel of normal VCF file
```

#### <a name="phasedsom"></a>Haplotype phased somatic mutation detection 

Himut leverages CCS read length and base accuracy to haplotype phase germline mutations. CCS reads have an average read length between 10kb and 20kb and consequently, a CCS read can span mulitple heterozaygous single nucleotide polymorphisms (hetSNPs). Himut identifies haplotype consistent hetSNPs that are connected to construct haplotype blocks. In the VCF file, haplotype phased hetSNPs in the same haplotype block have the same phase set (PS), which is the position of the first hetSNP in the haplotype block. 

```sh
himut phase --bam aln.primary_alignments.sorted.bam --vcf germline.vcf.bgz -o germline.phased.vcf ## phase
bgzip -c germline.phased.vcf  > germline.phased.vcf.bgz ## bgzip compress
tabix -p vcf germline.phased.vcf.bgz ## tabix index
```

Himut can take the haplotype phased VCF file as an additional input to haplotype phase CCS reads and limit somatic mutation detection to haplotype phased CCS reads. 

```sh
himut call -i aln.primary_alignments.sorted.bam --phased_vcf germline.phased.vcf.bgz -o somatic.vcf --phase # haplotype phased somatic mutation detection
himut call -i aln.primary_alignments.sorted.bam --phased_vcf germline.phased.vcf.bgz --common_snps ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.decomposed.normalised.common_snps.vcf.bgz --panel_of_normals b37.himut.pon.sbs.dbs.vcf.bgz -o somatic.vcf --phase # somatic mutation detection with the highest sensitivity 
bgzip -c somatic.vcf > somatic.vcf.bgz
tabix -p vcf somatic.vcf.bgz
```

#### <a name="SBS96"></a>SBS96 count calculation

There are 6 possible substitutions in the pyrimidine sequence context (C>A, C>G, C>T, T>A, T> and T>G) and 16 possible trinucleotide sequence contexts for each substitution, creating the SBS96 classification system (6 * 16 = 96). 

Himut accepts as input a VCF file with somatic mutations and a reference genome FASTA file, classifies each substitution and returns a file with SBS96 counts.

```sh
himut sbs96 -i somatic.vcf --ref ref.fa -o somatic.sbs96.tsv 
```

#### <a name="normcounts"></a>SBS96 count normalisation 

Before mutational signature extraction and analysis, SBS96 count needs to be normalised based on the trinucleotide frequencines in the reference genome FASTA file and the callable trinucleotide frequencies in CCS reads. This is because number of callable bases, the CCS bases from which somatic mutation detection could have been performed, is unique to each sample and each sequencing run. To calculate the number of callable CCS bases, SBS96 count normalisation parameters needs to be identical to the somatic mutation detection parameters, as the parameters changes the number of callable CCS bases. 

```sh
himut normcounts --bam aln.primary_alignments.sorted.bam --ref ref.fa --common_snps ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.decomposed.normalised.common_snps.vcf.bgz --panel_of_normals b37.himut.pon.sbs.dbs.vcf.bgz --sbs somatic.vcf.bgz -o somatic.normcounts.tsv ## unphased normalisation
himut normcounts --bam aln.primary_alignments.sorted.bam --ref ref.fa --phased_vcf germline.phased.vcf.bgz --common_snps ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.decomposed.normalised.common_snps.vcf.bgz --panel_of_normals b37.himut.pon.sbs.dbs.vcf.bgz --sbs somatic.vcf.bgz -o somatic.normcounts.tsv --phase ## phased normalisation
```

#### <a name="tol"></a>Somatic mutation detection across the Tree of Life

Himut enable somatic mutation detection from bulk normal tissue, agnostic of species and clonality, and can be used to detect somatic mutations across the Tree of Life. We used himut to call somatic mutations from Darwin Tree of Life (DToL) project samples. Himut requires germline SNP prior to calculate the genotype quality score for each germline muation and because the germline SNP prior is unknown, a deepvariant VCF file needs to be provided to calculate the germline SNP prior. How the germline SNP prior is calculated depends on whether CCS reads and reference genome is derived from the same sample or from a different sample. 

To indicate to himut that germline SNP prior needs to be calculated, provide `--non_human_sample` as a parameter. In addition, if the CCS reads and reference genome are derived from the same sample, provide `--reference_sample` as a parameter. In addition, because the DToL project sequences and assembles one sample per species, a VCF file with common SNPs and a Panel of Normal VCF file is not available. 

```sh
## unphased somatic mutation detection and SBS96 count normalisation where CCS reads and reference genome is derived from the same sample
himut call -i aln.primary_alignments.sorted.bam --vcf germline.vcf.bgz -o somatic.vcf --non_human_sample --reference_sample 
himut normcounts --bam aln.primary_alignments.sorted.bam --ref ref.fa --vcf germline.vcf --phased_vcf germline.phased.vcf.bgz --non_human_sample --reference_sample

## phased somatic mutation detection and SBS96 count normalisation where CCS reads and reference genome is derived from the same sample
himut call -i aln.primary_alignments.sorted.bam --vcf germline.vcf.bgz --phased_vcf germline.phased.vcf.bgz -o somatic.vcf --phase --non_human_sample --reference_sample
himut normcounts --bam aln.primary_alignments.sorted.bam --ref ref.fa --vcf germline.vcf --phased_vcf germline.phased.vcf.bgz --non_human_sample --reference_sample
```

If the CCS reads and reference genome are derived from a different sample, use the following commands.

```sh
## unphased somatic mutation detection and SBS96 count normalisation where CCS reads and reference genome is derived from a different sample
himut call -i aln.primary_alignments.sorted.bam --vcf germline.vcf.bgz -o somatic.vcf --non_human_sample
himut normcounts --bam aln.primary_alignments.sorted.bam --ref ref.fa --vcf germline.vcf --phased_vcf germline.phased.vcf.bgz --non_human_sample

## phased somatic mutation detection and SBS96 count normalisation where CCS reads and reference genome is derived from a different sample
himut call -i aln.primary_alignments.sorted.bam --vcf germline.vcf.bgz --phased_vcf germline.phased.vcf.bgz -o somatic.vcf --phase --non_human_sample
himut normcounts --bam aln.primary_alignments.sorted.bam --ref ref.fa --vcf germline.vcf --phased_vcf germline.phased.vcf.bgz --non_human_sample
```

#### <a name="pon"></a>Panel of Normal VCF file

A Panel of Normal VCF file is used to remove false positive substitutions resulting from non-canonical errors. To create a Panel of Normal VCF file, please provide himut with a minimap2 BAM file from an unrelated individual from the same species and the `--create_panel_of_normal` parameter. Himut will use less stringent thresholds to maximise the number of substitutions called from the sample and if the same found in a sample of interest, the substitution is not considered as a somatic mutation. 

```sh
himut call -i sample.aln.primary_alignments.sorted.bam -o sample.pon.vcf --create_panel_of_normal
```

We prepared a Panel of Normal VCF file for public use Google Drive ((https://drive.google.com/drive/folders/14-cSTBocwdyIemXTExuhW-mnTJb87SLT?usp=sharing).

#### <a name="pop"></a>VCF file with common SNPs

We selected common SNPs (>1% minor allele frequency) from the 1000 Genomes Project Phase 3 VCF file (https://drive.google.com/drive/folders/14-cSTBocwdyIemXTExuhW-mnTJb87SLT?usp=sharing). 


### <a name="limits"></a>Limitations

- Himut does not support Revio CCS reads, which has BQ score that ranges from Q1 to Q40. 
- PacBio CCS base accuracy is the limiting factor for somatic mutation calling using himut. PacBio circular consensus sequence (pbccs) algorithm calculates the base quality score depending on the dinucleotide sequence context and the number of supporting subread bases. Our research suggests that the base quality score does not correctly reflect the true base accuracy and that the consensus sequence caller does not account for the different substitution error rates in different trinucleotide sequence contexts. We discuss how this issue can be mitigated in our manuscript.
