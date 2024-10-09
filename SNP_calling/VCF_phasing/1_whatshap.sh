#!/bin/bash
#$ -cwd
#$ -R y
#$ -t 1-37
#$ -tc 6

# Script to perform read-based phasing of the raw vcfs using WhatsHap

# 0. Activate the conda environment where whatshap was installed

source activate whatshap

bamdir=/path/to/bams/clean_bams
vcf_file=/path/to/vcfs/auto/Puffinus_intersect_raw.auto.recode.vcf.gz
outdir=/path/to/phased_vcfs/autosomes
ref_genome=/path/to/genome.fa

final_vcf=$outdir/Puffinus_raw.auto.whatshap_phased.vcf.gz

# 1. Retrieve the names of all the individuals

bcftools query -l $vcf_file > $outdir/ind_list.txt
inds=($(cat $outdir/ind_list.txt))

# 2. Should we make sure that no phasing is included in the original vcf? Not really necessary, as it was done only to very obvious positions using freebayes and the result should be the same

# 3. Divide vcf into different individual vcfs. Don't remove monomorphic sites! Obtained from: https://www.biostars.org/p/224702/

bcftools view -Oz -s ${inds[(($SGE_TASK_ID))]} -o $outdir/ind_vcfs/${inds[(($SGE_TASK_ID))]}\.vcf.gz $vcf_file

tabix $outdir/ind_vcfs/${inds[(($SGE_TASK_ID))]}\.vcf.gz

# 4. Perform phasing using whatshap - use the non-masked reference genome

whatshap phase -o $outdir/${inds[(($SGE_TASK_ID))]}\.whatshap_phased.vcf.gz --reference=$ref_genome $outdir/ind_vcfs/${inds[(($SGE_TASK_ID))]}\.vcf.gz $bamdir/${inds[(($SGE_TASK_ID))]}\.sorted.dups.bam

# 5. Index the output files

tabix $outdir/${inds[(($SGE_TASK_ID))]}\.whatshap_phased.vcf.gz

# 6. Merge individual vcfs - previouly write a file (merge.txt) with all the phased vcf names

ls $outdir/*.whatshap_phased.vcf.gz > $outdir/merge.txt

bcftools merge -l $outdir/merge.txt -Oz -o $final_vcf --threads 6
tabix $final_vcf
