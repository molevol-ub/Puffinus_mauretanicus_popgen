#!/bin/bash
#$ -cwd
#$ -V

# Script to extract only the selected chrZ or autosomal scaffolds from the vcf file

# 1 Define files and pathways

VCF_file=/set/your/path/Puffinus_intersect_raw.vcf.gz
guillem=/set/your/path

bed_file=$guillem/auto/auto.bed
bed_file2=$guillem/auto/chrZ.bed

# 2 Filter the desired scaffolds using vcftools

vcftools --gzvcf $VCF_file --bed $bed_file --recode --out $guillem/Puffinus_intersect_raw.auto
vcftools --gzvcf $VCF_file --bed $bed_file2 --recode --out $guillem/Puffinus_intersect_raw.chrZ
