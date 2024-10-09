#!/bin/bash
#$ -cwd
#$ -V

# Script to extract only the selected chrZ scaffolds from the vcf file

# 1 Define files and pathways

VCF_file=/users-d3/jferrer/pmau_popgen/SNP_calling/vcfs/Puffinus_intersect_raw.vcf.gz
guillem=/users-d3/jferrer/gizquierdo/TFG/chr_split/vcfs

bed_file=$guillem/auto/auto.bed

# 2 Filter the desired scaffolds using vcftools

vcftools --gzvcf $VCF_file --bed $bed_file --recode --out $guillem/Puffinus_intersect_raw.auto
