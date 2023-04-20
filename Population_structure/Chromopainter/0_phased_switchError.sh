#!/bin/bash
#$ -cwd
#$ -R y
#$ -e logs/2_shapeit4.err
#$ -o logs/2_shapeit4.out
#$ -q h13.q
#$ -pe ompi255h13 6
#$ -V                    #export environment var
#$ -N shapeit4             #name Job
#$ -t 1-876
#$ -tc 7

# Information about switchError is found in pull requests in: https://github.com/SPG-group/switchError

outdir=/users-d3/jferrer/pmau_popgen/SNP_calling/vcfs/phased_vcfs_def/autosomes
ref_genome=/users/ccuevas/Functional_Genome_Annotation/scripts/braker3/3pmaureta/genome.fa
vcfdir=$outdir

cd $outdir

# 1. Activate the conda enviroment were switchError is installed

source activate whatshap
#source activate switchError

# 2. Run shapeit with all but 1 individual; just for trial we'll iterate through the removal of 5 individuals: CZA11, G15, M19, G11, M17

for ind in (CZA11 G15 M19 G11 M17)
do
vcftools --gzvcf $vcfdir/whatshap_phased/Puffinus_raw.auto.whatshap_phased.vcf.gz --remove-indv $ind --recode-INFO-all --recode --out $vcfdir/whatshap_phased/Puffinus_raw.$ind.whatshap_phased
mv $vcfdir/whatshap_phased/Puffinus_raw.$ind.whatshap_phased.recode.vcf $vcfdir/whatshap_phased/Puffinus_raw.$ind.whatshap_phased.vcf
bgzip $vcfdir/whatshap_phased/Puffinus_raw.$ind.whatshap_phased.vcf
tabix $vcfdir/whatshap_phased/Puffinus_raw.$ind.whatshap_phased.vcf.gz
done

for ind in CZA11 G15 M19 G11 M17
do

# Start of the script of shapeit4---------------------------

scfdir=$outdir/scaffolds_$ind
mkdir $scfdir
input_vcf=$outdir/whatshap_phased/Puffinus_raw.$ind.whatshap_phased.vcf.gz
output_vcf=$outdir/Puffinus_raw.$ind.shapeit4_whatshap_phased.vcf.gz

# 2.1. Retrieve the names of all scaffolds

zcat $input_vcf | grep -v "^#" | cut -f1 | sort | uniq > $outdir/scf_list.txt
scf=($(cat $outdir/scf_list.txt))

# 2.2. Perform phasing using whatshap - use the non-masked reference genome

shapeit4 --input $input_vcf --region ${scf[(($SGE_TASK_ID))]} --use-PS 0.0001 --thread 6 --output $scfdir/${scf[(($SGE_TASK_ID))]}\_whatshap_shapeit4_phased.vcf.gz

# 2.3. Index the output files

tabix $scfdir/${scf[(($SGE_TASK_ID))]}\_whatshap_shapeit4_phased.vcf.gz

# 2.4. Merge individual vcfs - previouly write a file (merge.txt) with all the phased vcf names

ls $scfdir/*_whatshap_shapeit4_phased.vcf.gz > $outdir/merge_scf.txt
bcftools concat -f $outdir/merge_scf.txt -Oz -o $output_vcf --threads 6
tabix $output_vcf

# Stop of the script of shapeit4----------------------------

done

# 3. Generate input bcf files for switchError (1st shapeit, 2nd whatshap-phased)

cd 
vcfdir=/users-d3/jferrer/pmau_popgen/SNP_calling/vcfs/phased_vcfs_def/autosomes

for ind in (CZA11 G15 M19 G11 M17)
do
bcftools --view --Ou -o Puffinus_shapeit.$ind.bcf $vcfdir/Puffinus_raw.$ind.shapeit4_whatshap_phased.vcf.gz
done

bcftools --view --Ou -o Puffinus_whatshap.bcf $vcfdir/whatshap_phased/Puffinus_raw.auto.whatshap_phased.vcf.gz
