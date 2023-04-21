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

# 1. Transform VCF to chromopainter input (if it works)

invcf=/users-d3/jferrer/pmau_popgen/SNP_calling/vcfs/phased_vcfs_def/autosomes/Puffinus_def.auto.shapeit4_whatshap_phased.nPP.maxmiss80.vcf.gz
invcf_nomono=/users-d3/jferrer/pmau_popgen/SNP_calling/vcfs/phased_vcfs_def/autosomes/Puffinus_def.auto.shapeit4_whatshap_phased.nPP.maxmiss80.vcf

bcftools view -Ov -o $invcf_nomono --threads 4 -e 'COUNT(GT="AA")=N_SAMPLES || COUNT(GT="RR")=N_SAMPLES' $invcf

vcf_to_chromopainter_main.R -g $invcf_nomono -o Puffinus.auto.nPP.maxmiss80
