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
#$ -m e

# Script to find the fixed differences between 2 populations (whose .list files you must have prepared beforehand)

cd /users-d3/jferrer/gizquierdo/TFM/fixed_differences

dataset=Pmed_Ppuf
pop1=Pmed
pop2=Ppuf

# 1 Generate VCF file with only inds from the 2 populations that you wanna compare

vcftools --gzvcf ~/gizquierdo/TFG/chr_split/vcfs/chrZ/Puffinus_SNP.maxmiss65.filtered.wPP.merged.nomono.masked.vcf.gz --keep ./lists/$dataset.list --recode --out $dataset

mv $dataset.recode.vcf $dataset.vcf
bgzip $dataset.vcf
tabix $dataset.vcf.gz

# 2 Use BCFtools to remove non-monomorphic sites

bcftools view -Oz -o $dataset.nomono.vcf.gz -e 'COUNT(GT="AA")=N_SAMPLES || COUNT(GT="RR")=N_SAMPLES' $dataset.vcf.gz
tabix $dataset.nomono.vcf.gz

# 3 Use BCFtools -x (--private) twice to keep positions where only 1 population retains the derived alleles

bcftools view -Oz -o $dataset.$pop1.unique_derived.vcf.gz -x -S ./lists/$pop1.list $dataset.nomono.vcf.gz
bcftools view -Oz -o $dataset.$pop2.unique_derived.vcf.gz -x -S ./lists/$pop2.list $dataset.nomono.vcf.gz

tabix $dataset.$pop1.unique_derived.vcf.gz
tabix $dataset.$pop2.unique_derived.vcf.gz

# 4 Now eliminate ell HET and HOM.REF inds so you obtain only fixed alleles

bcftools view -Oz -o $dataset.$pop1.unique_derived.fixed.vcf.gz -e 'COUNT(GT="AR")!=0 || COUNT(GT="RR")!=0' $dataset.$pop1.unique_derived.vcf.gz
bcftools view -Oz -o $dataset.$pop2.unique_derived.fixed.vcf.gz -e 'COUNT(GT="AR")!=0 || COUNT(GT="RR")!=0' $dataset.$pop2.unique_derived.vcf.gz

tabix $dataset.$pop1.unique_derived.fixed.vcf.gz
tabix $dataset.$pop1.unique_derived.fixed.vcf.gz

