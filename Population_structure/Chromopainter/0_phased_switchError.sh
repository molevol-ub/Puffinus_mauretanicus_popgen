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

cd /users-d3/jferrer/pmau_popgen/SNP_calling/vcfs/phased_vcfs_def/autosomes/switchError

# 3. Use bedtools intersecte to keep only positions in common between both files (?)
# 4. Sort both vcfs in the same order

grep "^#" Puffinus_shapeit.M19.def.vcf > prova
grep -v "^#" Puffinus_shapeit.M19.def.vcf | sort -k1,1V -k2,2g >> prova 
mv prova Puffinus_shapeit.M19.def.vcf 

grep "^#" Puffinus_whatshap.M19.def.vcf > prova
grep -v "^#" Puffinus_whatshap.M19.def.vcf | sort -k1,1V -k2,2g >> prova 
mv prova Puffinus_whatshap.M19.def.vcf 

# 5. Calculate switch error

vcftools --vcf Puffinus_whatshap.M19.def.vcf --diff Puffinus_shapeit.M19.def.vcf --diff-switch-error --out Puffinus.M19

# We can try with just the phased positions for the whatshap file: 1st keep those positions from the whatshap file (and sort them) and then intersect with the shapeit file

grep "^#" Puffinus_whatshap.M19.def.vcf > Puffinus_whatshap.M19.only_phased.vcf
grep "1|1" Puffinus_whatshap.M19.def.vcf > prova
grep "0|1" Puffinus_whatshap.M19.def.vcf >> prova
grep "1|0" Puffinus_whatshap.M19.def.vcf >> prova
grep "0|0" Puffinus_whatshap.M19.def.vcf >> prova
grep -v "^#" prova | sort -k1,1V -k2,2g >> Puffinus_whatshap.M19.only_phased.vcf 

grep "^#" Puffinus_whatshap.M19.def.vcf > Puffinus_shapeit.M19.only_phased.vcf
bedtools intersect -a Puffinus_shapeit.M19.def.vcf -b Puffinus_whatshap.M19.only_phased.vcf >> Puffinus_shapeit.M19.only_phased.vcf

vcftools --vcf Puffinus_whatshap.M19.only_phased.vcf --diff Puffinus_shapeit.M19.only_phased.vcf --diff-switch-error --out Puffinus.M19.only_phased

#----------------------------------------------------------------------------------------------------

# 6. Test whether switch error decreases when eliminating windows with low SNP density, which should be more prone to bad phasing

# 6.1. Obtain the SNP density in 25 kb windows (the non-modified whatshap-phased VCF shoudl work)

vcftools --vcf Puffinus_whatshap.M19.def.vcf --SNPdensity 25000 --out Puffinus.M19

# 6.2 Calculate the threshold (90% above)

thres=$(sort -g -r -k4 Puffinus.M19.snpden | awk '{all[NR] = $4} END{print all[int(NR*0.9 - 0.5)]}')

# 6.3. Keep windows above threshold

awk -v thres=$thres '$4<thres{next}1' Puffinus.M19.snpden > prova
cat prova | cut -f 1,2 > Puffinus.M19.snpden
awk '$1 == "CHROM" {print $1, $2, $3 = "BIN_END"}' Puffinus.M19.snpden > prova
awk '$1 != "CHROM" {print $1, $2, $3 = $2+24999}' Puffinus.M19.snpden >> prova
# Change separators to tabs
sed 's/ /\t/g' prova > Puffinus.M19.snpden

# 6.4. Include only these positions in the VCFs and later calculate the switch error

vcftools --vcf Puffinus_shapeit.M19.only_phased.vcf --bed Puffinus.M19.snpden --recode-INFO-all --recode --out prova
mv prova.recode.vcf Puffinus_shapeit.M19.only_phased.vcf

vcftools --vcf Puffinus_whatshap.M19.only_phased.vcf --bed Puffinus.M19.snpden --recode-INFO-all --recode --out prova
mv prova.recode.vcf Puffinus_whatshap.M19.only_phased.vcf
