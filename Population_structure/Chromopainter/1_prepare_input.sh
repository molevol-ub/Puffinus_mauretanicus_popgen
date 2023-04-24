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

# 1. Prepare your VCFs with the desired SNPs
# You must keep the INFO field of shapeit4-phased VCFs to allow the inference of recombination rates

invcf=/users-d3/jferrer/pmau_popgen/SNP_calling/vcfs/phased_vcfs_def/autosomes/Puffinus_raw.auto.shapeit4_whatshap_phased.vcf.gz
outvcf=/users-d3/jferrer/pmau_popgen/SNP_calling/vcfs/phased_vcfs_def/autosomes/Puffinus_def.auto.shapeit4_whatshap_phased.nPP.maxmiss80
outvcf_nomono=/users-d3/jferrer/pmau_popgen/SNP_calling/vcfs/phased_vcfs_def/autosomes/Puffinus_def.auto.shapeit4_whatshap_phased.nPP.maxmiss80.nomono.vcf

vcftools --gzvcf $invcf --remove-indv M22 --remove-indv sacella --remove-indv M7 --remove-indv COP1 --remove-indv LT2 --recode-INFO-all --recode --out $outvcf
mv $outvcf.recode.vcf $outvcf.vcf

grep "^#" $outvcf.vcf > prova
bedtools intersect -a $outvcf.vcf -b /users-d3/jferrer/gizquierdo/TFG/chr_split/vcfs/auto/Puffinus_SNP.maxmiss80.filtered.noPP.merged.nomono.masked.vcf.gz* | grep -v "^#" >> prova
mv prova $outvcf.vcf

bcftools view -Ov -o $outvcf_nomono --threads 4 -e 'COUNT(GT="AA")=N_SAMPLES || COUNT(GT="RR")=N_SAMPLES' $outvcf.vcf

# 2. Transform VCF to chromopainter input (if it works) using the script provided in: https://github.com/sahwa/vcf_to_chromopainter

vcf_to_chromopainter_main.R -g $outvcf_nomono -o Puffinus.auto.nPP.maxmiss80

# 3. Remove spaces from the "haplotype" lines in the haplotype input file

head Puffinus.auto.nPP.maxmiss80.chromopainter.inp -n 3 > prova
tail -n +4 Puffinus.auto.nPP.maxmiss80.chromopainter.inp | sed 's/ //g' >> prova
mv prova Puffinus.auto.nPP.maxmiss80.chromopainter.inp

#----------------------------------------------------------------------------------------

# 4. Repeat these two last steps with 1/2 of the indvs and 1/10 of the scaffolds to generate a reduced dataset to run the parameter estimation step

#4.1. Count scaffolds from a list (scaffold.txt) and divide into 10; then write into a new file
nline=$(wc -l scaffolds.txt | sed 's/ scaffolds.txt//g')
nline=$(echo $nline/10 | bc)
shuf -n $nline scaffolds.txt > scaffolds.reduced.txt 

#4.2. Format the file as a bed file (with large "chrom_end" values so that it includes all snps)

awk -F'\t' 'BEGIN {OFS=FS} {print $1, $2=0, $3 = 50000000}' scaffolds.reduced.txt > prova.txt
mv prova.txt scaffolds_reduced.txt

#4.3. Filter the vcf

vcftools --vcf $outvcf_nomono --keep Puffinus.txt --bed scaffolds_reduced.txt --recode-INFO-all --recode --out Puffinus_subset
mv Puffinus_subset.recode.vcf Puffinus_subset.vcf

# If the number of SNPs is too small (due to the fact that the randomly-chosen scaffolds are too few) repeat all three steps

#4.4. Transform VCF to chromopainter input (if it works) using the script provided in: https://github.com/sahwa/vcf_to_chromopainter

vcf_to_chromopainter_main.R -g Puffinus_subset.vcf -o Puffinus_subset

#4.5. Remove spaces from the "haplotype" lines in the haplotype input file

head Puffinus_subset.chromopainter.inp -n 3 > prova
tail -n +4 Puffinus_subset.chromopainter.inp | sed 's/ //g' >> prova
mv prova Puffinus_subset.chromopainter.inp

#-------------------------------------------------------------------------------------------

# 5. Modify the input recombination file so that the last recombination of every gene is substituted by "-9" to accomodate the requirements of chromopainter

scripts_dir=/users-d3/jferrer/gizquierdo/TFM/chromopainter/scripts
python $scripts_dir/correct_recomb_file.py Puffinus_subset.recomrates.txt prova.txt

# Copy the last line, that is not included in the script

tail -n 1 Puffinus_subset.recomrates.txt >> prova.txt

mv prova.txt Puffinus_subset.recomrates.txt
