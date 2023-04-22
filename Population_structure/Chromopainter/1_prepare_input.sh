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

vcftools --gzvcf $invcf --max-missing 0.8 --min-alleles 2 --max-alleles 2 --minDP 4 --maxDP 50 --minQ 30 --remove-indv M22 --remove-indv sacella
--remove-indv M7 --remove-indv COP1 --remove-indv LT2 --remove-filtered-all --recode-INFO-all --recode --out 
mv $outvcf.recode.vcf $outvcf.vcf

bcftools view -Ov -o $outvcf_nomono --threads 4 -e 'COUNT(GT="AA")=N_SAMPLES || COUNT(GT="RR")=N_SAMPLES' $outvcf.vcf

vcf_to_chromopainter_main.R -g $outvcf_nomono -o Puffinus.auto.nPP.maxmiss80
