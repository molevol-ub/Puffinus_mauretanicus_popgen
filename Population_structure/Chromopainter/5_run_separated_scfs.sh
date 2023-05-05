#!/bin/bash
#$ -cwd
#$ -R y
#$ -e divide_vcf.err
#$ -o divide_vcf.out
#$ -q h11.q
#$ -pe ompi64h11 1
#$ -V                    #export environment var
#$ -N divide_vcf             #name Job
#$ -t 1-876
#$ -tc 24

workdir=/users-d3/jferrer/gizquierdo/TFM/chromopainter/scaffolds
#mkdir $workdir
cd $workdir

# 0. Make your list of scaffolds

outvcf_nomono=/users-d3/jferrer/pmau_popgen/SNP_calling/vcfs/phased_vcfs_def/autosomes/Puffinus_def.auto.shapeit4_whatshap_phased.nPP.maxmiss80.nomono.vcf
cat $outvcf_nomono | grep -v "^#" | cut -f1 | sort | uniq > scaffolds.txt
scf=($(cat scaffolds.txt))

# 1. Divide your VCF by scaffolds

#while read line; do
#  vcftools --vcf $outvcf_nomono --chr $line --recode-INFO-all --recode --out $line
#  mv $line.recode.vcf $line.vcf
#done < $scf_list

vcftools --vcf $outvcf_nomono --chr ${scf[(($SGE_TASK_ID))]} --recode-INFO-all --recode --out ${scf[(($SGE_TASK_ID))]}
mv ${scf[(($SGE_TASK_ID))]}\.recode.vcf ${scf[(($SGE_TASK_ID))]}\.vcf
