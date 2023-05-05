#!/bin/bash
#$ -cwd
#$ -R y
#$ -e divide_vcf.err
#$ -o divide_vcf.out
#$ -q h12.q
#$ -pe ompi128h12 2
#$ -V                    #export environment var
#$ -N divide_vcf             #name Job
#$ -M 000izquierdoguillem@gmail.com
#$ -m be

workdir=/users-d3/jferrer/gizquierdo/TFM/chromopainter/scaffolds
mkdir $workdir
cd $workdir

# 1. Divide your VCF by scaffolds

outvcf_nomono=/users-d3/jferrer/pmau_popgen/SNP_calling/vcfs/phased_vcfs_def/autosomes/Puffinus_def.auto.shapeit4_whatshap_phased.nPP.maxmiss80.nomono.vcf
scf_list=/users-d3/jferrer/pmau_popgen/SNP_calling/vcfs/phased_vcfs_def/autosomes/scf_list.txt

while read line; do
  vcftools --vcf $outvcf_nomono --chr $line --recode-INFO-all --recode --out $line
  mv $line.recode.vcf $line.vcf
done < $scf_list
