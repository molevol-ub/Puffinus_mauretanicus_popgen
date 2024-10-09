#!/bin/bash
#$ -cwd
#$ -R y
#$ -t 1-876
#$ -tc 7

# Script to perform statistical phasing using Shapeit4 over the Whatshap-phased VCF

# 0. Activate the conda environment where Shapeit4 was installed

source activate whatshap

outdir=/path/to/phased_vcfs/autosomes
ref_genome=/path/to/genome.fa

scfdir=$outdir/scaffolds_hom

input_vcf=$outdir/whatshap_phased/Puffinus_raw.auto.whatshap_phased.vcf.gz
output_vcf=$outdir/Puffinus_raw.auto.shapeit4_whatshap_phased.vcf.gz

# 1. Retrieve the names of all scaffolds

zcat $input_vcf | grep -v "^#" | cut -f1 | sort | uniq > $outdir/scf_list.txt

scf=($(cat $outdir/scf_list.txt))

# 2. Perform phasing using whatshap - use the non-masked reference genome

shapeit4 --input $input_vcf --region ${scf[(($SGE_TASK_ID))]} --use-PS 0.0001 --thread 6 --output $scfdir/${scf[(($SGE_TASK_ID))]}\_whatshap_shapeit4_phased.vcf.gz

# 3. Index the output files

tabix $scfdir/${scf[(($SGE_TASK_ID))]}\_whatshap_shapeit4_phased.vcf.gz

# 4. Merge individual vcfs - previouly write a file (merge.txt) with all the phased vcf names

ls $scfdir/*_whatshap_shapeit4_phased.vcf.gz > $outdir/merge_scf.txt

bcftools concat -f $outdir/merge_scf.txt -Oz -o $output_vcf --threads 6
tabix $output_vcf
