#!/bin/bash
#$ -cwd
#$ -R y
#$ -e fixed_diff.err
#$ -o fixed_diff.out
#$ -q h13.q
#$ -pe ompi255h13 4
#$ -V                    #export environment var
#$ -N fixed_diff            #name Job
#$ -M 000izquierdoguillem@gmail.com
#$ -m be

# Script to prepare the input vcf for SNAPP

# 1. Make a vcf file with only representatative individuals of the populations/groups you want to include

vcftools --gzvcf /path/to/your/main_file.vcf.gz --indv ind1 --indv ind2 --indv ind3 --recode --our dataset_downsized

mv dataset_downsized.recode.vcf dataset_downsized.vcf
bgzip dataset_downsized.vcf
tabix dataset_downsized.vcf.gz

# 2. Make a 5 random subsamples of your file with 1k datasets

# -r --> ratio of the nº of snps you want to keep vs total. For ex.: if you want to keep 1k SNPs from a 1M SNP original dataset, then -r 0.001

for i in {1..5}
do

bcftools view dataset_downsized.vcf.gz | vcfrandomsample -r 0.001 > dataset_downsized.vcf.gz.$i.vcf 

done
