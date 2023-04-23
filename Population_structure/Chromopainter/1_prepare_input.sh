#!/bin/bash
#$ -cwd
#$ -R y
#$ -e gatk_reference.err
#$ -o gatk_reference.out
#$ -q h12.q
#$ -pe ompi128h12 4
#$ -V                    #export environment var
#$ -N gatk_reference             #name Job
#$ -M 000izquierdoguillem@gmail.com
#$ -m be


workdir=/users-d3/jferrer/gizquierdo/TFM/chromopainter
cd $workdir

# 1. Transform VCF to chromopainter input (if it works) using the script provided in: https://github.com/sahwa/vcf_to_chromopainter
# You must keep the INFO field of shapeit4-phased VCFs to allow the inference of recombination rates

invcf=/users-d3/jferrer/pmau_popgen/SNP_calling/vcfs/phased_vcfs_def/autosomes/Puffinus_raw.auto.shapeit4_whatshap_phased.vcf.gz
outvcf=/users-d3/jferrer/pmau_popgen/SNP_calling/vcfs/phased_vcfs_def/autosomes/Puffinus_def.auto.shapeit4_whatshap_phased.nPP.maxmiss80
outvcf_nomono=/users-d3/jferrer/pmau_popgen/SNP_calling/vcfs/phased_vcfs_def/autosomes/Puffinus_def.auto.shapeit4_whatshap_phased.nPP.maxmiss80.vcf

vcftools --gzvcf $invcf --remove-indv M22 --remove-indv sacella --remove-indv M7 --remove-indv COP1 --remove-indv LT2 --recode-INFO-all --recode --out $outvcf
mv $outvcf.recode.vcf $outvcf.vcf

grep "^#" $outvcf.vcf > prova
bedtools intersect -a $outvcf.vcf -b /users-d3/jferrer/gizquierdo/TFG/chr_split/vcfs/auto/Puffinus_SNP.maxmiss80.filtered.noPP.merged.nomono.masked.vcf.gz* | grep -v "^#" >> prova
mv prova $outvcf.vcf

bcftools view -Ov -o $outvcf_nomono --threads 4 -e 'COUNT(GT="AA")=N_SAMPLES || COUNT(GT="RR")=N_SAMPLES' $outvcf.vcf



vcf_to_chromopainter_main.R -g $outvcf_nomono -o Puffinus.auto.nPP.maxmiss80
