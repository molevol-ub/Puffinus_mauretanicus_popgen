#!/bin/bash
#$ -cwd
#$ -R y
#$ -e logs/4_switchError.err
#$ -o logs/4_switchError.out
#$ -q h13.q
#$ -pe ompi255h13 6
#$ -V                    #export environment var
#$ -N switchError             #name Job
#$ -t 1-876
#$ -tc 6

# Information about switchError is found in pull requests in: https://github.com/SPG-group/switchError

outdir=/users-d3/jferrer/pmau_popgen/SNP_calling/vcfs/phased_vcfs_def/autosomes
ref_genome=/users/ccuevas/Functional_Genome_Annotation/scripts/braker3/3pmaureta/genome.fa
vcfdir=$outdir

cd $outdir

# 1. Activate the conda enviroment were switchError is installed

source activate whatshap

# 2. Run shapeit with all but 1 individual; just for trial we'll iterate through the removal of 5 individuals: CZA11, G15, M19, G11, M17


#vcftools --gzvcf $vcfdir/whatshap_phased/Puffinus_raw.auto.whatshap_phased.vcf.gz --remove-indv M19 --recode-INFO-all --recode --out $vcfdir/whatshap_phased/Puffinus_raw.M19.whatshap_phased
#mv $vcfdir/whatshap_phased/Puffinus_raw.M19.whatshap_phased.recode.vcf $vcfdir/whatshap_phased/Puffinus_raw.M19.whatshap_phased.vcf
#bgzip $vcfdir/whatshap_phased/Puffinus_raw.M19.whatshap_phased.vcf
#tabix $vcfdir/whatshap_phased/Puffinus_raw.M19.whatshap_phased.vcf.gz


vcftools --gzvcf $vcfdir/whatshap_phased/Puffinus_raw.auto.whatshap_phased.vcf.gz --indv M19 --recode-INFO-all --recode --out $vcfdir/whatshap_phased/Puffinus_raw.only_M19.whatshap_phased
mv $vcfdir/whatshap_phased/Puffinus_raw.only_M19.whatshap_phased.recode.vcf $vcfdir/whatshap_phased/Puffinus_raw.only_M19.whatshap_phased.vcf

cat $vcfdir/whatshap_phased/Puffinus_raw.only_M19.whatshap_phased.vcf | sed 's/1|1/1\/1/g' > prova
mv prova $vcfdir/whatshap_phased/Puffinus_raw.only_M19.whatshap_phased.vcf
cat $vcfdir/whatshap_phased/Puffinus_raw.only_M19.whatshap_phased.vcf | sed 's/0|0/0\/0/g' > prova
mv prova $vcfdir/whatshap_phased/Puffinus_raw.only_M19.whatshap_phased.vcf
cat $vcfdir/whatshap_phased/Puffinus_raw.only_M19.whatshap_phased.vcf | sed 's/1|0/1\/0/g' > prova
mv prova $vcfdir/whatshap_phased/Puffinus_raw.only_M19.whatshap_phased.vcf
cat $vcfdir/whatshap_phased/Puffinus_raw.only_M19.whatshap_phased.vcf | sed 's/0|1/0\/1/g' > prova
mv prova $vcfdir/whatshap_phased/Puffinus_raw.only_M19.whatshap_phased.vcf

bgzip $vcfdir/whatshap_phased/Puffinus_raw.only_M19.whatshap_phased.vcf
tabix $vcfdir/whatshap_phased/Puffinus_raw.only_M19.whatshap_phased.vcf.gz

bcftools merge -Oz -o $vcfdir/whatshap_phased/Puffinus_raw.M19.def.whatshap_phased.vcf.gz $vcfdir/whatshap_phased/Puffinus_raw.only_M19.whatshap_phased.vcf.gz $vcfdir/whatshap_phased/Puffinus_raw.M19.whatshap_phased.vcf.gz
tabix Puffinus_raw.M19.def.whatshap_phased.vcf.gz

# Start of the script of shapeit4---------------------------

scfdir=$outdir/scaffolds_M19
input_vcf=$outdir/whatshap_phased/Puffinus_raw.M19.def.whatshap_phased.vcf.gz
output_vcf=$outdir/Puffinus_raw.M19.shapeit4_whatshap_phased.vcf.gz

# 2.1. Retrieve the names of all scaffolds

zcat $input_vcf | grep -v "^#" | cut -f1 | sort | uniq > $outdir/scf_list.txt
scf=($(cat $outdir/scf_list.txt))

# 2.2. Perform phasing using whatshap - use the non-masked reference genome

shapeit4 --input $input_vcf --region ${scf[(($SGE_TASK_ID))]} --use-PS 0.0001 --thread 5 --output $scfdir/${scf[(($SGE_TASK_ID))]}\_whatshap_shapeit4_phased.vcf.gz

# 2.3. Index the output files

tabix $scfdir/${scf[(($SGE_TASK_ID))]}\_whatshap_shapeit4_phased.vcf.gz

# 2.4. Merge individual vcfs - previouly write a file (merge.txt) with all the phased vcf names

ls $scfdir/*_whatshap_shapeit4_phased.vcf.gz > $outdir/merge_scf.txt
bcftools concat -f $outdir/merge_scf.txt -Oz -o $output_vcf --threads 5
tabix $output_vcf

# Stop of the script of shapeit4----------------------------

# 3. Generate input bcf files for switchError (1st shapeit, 2nd whatshap-phased)

cd /users-d3/jferrer/pmau_popgen/SNP_calling/vcfs/phased_vcfs_def/autosomes/switchError
vcfdir=/users-d3/jferrer/pmau_popgen/SNP_calling/vcfs/phased_vcfs_def/autosomes

bcftools view -Ou -o Puffinus_shapeit.M19.bcf -s M19 $vcfdir/Puffinus_raw.M19.shapeit4_whatshap_phased.vcf.gz
bcftools view -Ou -o Puffinus_whatshap.bcf -s M19 $vcfdir/whatshap_phased/Puffinus_raw.auto.whatshap_phased.vcf.gz

# 4. Run switchError were, if I understand it correctly (https://github.com/SPG-group/switchError/pull/3/files)
#       --gen includes the phased independently phased individuals (whatshap)
#       --hap includes the statistically phased individuals (shapeit4)
#       --reg includes the regions to be evaluated

source activate switchError    #Activate the environment

switchError --gen Puffinus_whatshap.bcf --hap Puffinus_shapeit.M19.bcf --reg /users-d3/jferrer/gizquierdo/TFM/genome_scans/PopGenome/auto_def_wM8/25_windows.bed --out Puffinus.M19
