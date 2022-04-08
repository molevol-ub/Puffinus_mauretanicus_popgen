#!/bin/bash
#$ -R y
#$ -N refine_filter
#$ -e logs/refine_filter_mtDNA.err.txt
#$ -o logs/refine_filter_mtDNA.out.txt
#$ -cwd
#$ -q h13.q
#$ -pe ompi255h13 16

source $HOME/.profile

### Script to filter and refine vcfs ###
MASTER=/users-d3/jferrer/pmau_popgen
reference=/$MASTER/SNP_calling/reference/Pmau_mitogenome_19885.fasta
raw_vcf=$MASTER/SNP_calling/mt_freebayes/tmp_mtDNA/Puffinus_mitogenome.vcf
DATASET=Puffinus_mtDNA
WORKING_DIR=$MASTER/SNP_calling/mt_vcfs

## Popmaps ## - These have been filtered for low coverage individuals, make sure all are in one directory
SNP_filtered=$WORKING_DIR/${DATASET}_SNP_filter.vcf
allele_filtered=$WORKING_DIR/${DATASET}_SNP.minmax2.mindp9.filtered
maxmiss_50_filtered=$WORKING_DIR/${DATASET}_SNP.maxmiss50.filtered
maxmiss_80_filtered=$WORKING_DIR/${DATASET}_SNP.maxmiss80.filtered
maxmiss_100_filtered=$WORKING_DIR/${DATASET}_SNP.maxmiss100.filtered

## Processing ...

### Use vcftools to filter out indels and SNPs with quality below 750 (and above 7500)? ##

cat $raw_vcf | vcffilter -f "TYPE = snp & QUAL > 750"> $SNP_filtered

# Other way:
# vcftools --vcf $raw_vcf --remove-indels --minQ 750 --recode --remove-filtered-all --out $SNP_filtered

# Remember to remove related inds!

## Retain only biallelic SNPs with a min depth of 9 and min&max_mean_depth between 41-55 - similarly, filter out the "repeated inds" and P. puffinus (if that's the case)
vcftools --vcf $SNP_filtered --min-alleles 2 --max-alleles 2 --min-meanDP 41 --max-meanDP 55 --minDP 9 --remove-indv M22 --remove-indv sacella --remove-indv M7 --recode --remove-filtered-all --out $allele_filtered

#Intuitively I would have chosen 1 allele, but I guess not (erroneous results)

## Filter by missing data (should we separate pops?) and mac both 1 & 2
vcftools --vcf $allele_filtered.recode.vcf --max-missing 0.50 --recode --remove-filtered-all --out ${maxmiss_50_filtered}
vcftools --vcf $allele_filtered.recode.vcf --max-missing 0.80 --recode --remove-filtered-all --out ${maxmiss_80_filtered}
vcftools --vcf $allele_filtered.recode.vcf --max-missing 1.00 --recode --remove-filtered-all --out ${maxmiss_100_filtered}

vcftools --vcf $allele_filtered.recode.vcf --max-missing 0.50 --mac 2 --recode --remove-filtered-all --out ${maxmiss_50_filtered}.mac2
vcftools --vcf $allele_filtered.recode.vcf --max-missing 0.80 --mac 2 --recode --remove-filtered-all --out ${maxmiss_80_filtered}.mac2
vcftools --vcf $allele_filtered.recode.vcf --max-missing 1.00 --mac 2 --recode --remove-filtered-all --out ${maxmiss_100_filtered}.mac2

# Also: remove monomorphic sites using bcftools with the gzipped files

bgzip $maxmiss_50_filtered.recode.vcf
tabix $maxmiss_50_filtered.recode.vcf.gz
bgzip $maxmiss_80_filtered.recode.vcf
tabix $maxmiss_80_filtered.recode.vcf.gz
bgzip $maxmiss_100_filtered.recode.vcf
tabix $maxmiss_100_filtered.recode.vcf.gz

bgzip $maxmiss_50_filtered.mac2.recode.vcf
tabix $maxmiss_50_filtered.mac2.recode.vcf.gz
bgzip $maxmiss_80_filtered.mac2.recode.vcf
tabix $maxmiss_80_filtered.mac2.recode.vcf.gz
bgzip $maxmiss_100_filtered.mac2.recode.vcf
tabix $maxmiss_100_filtered.mac2.recode.vcf.gz

bcftools view -Oz -o $maxmiss_50_filtered.nomono.vcf.gz -e 'COUNT(GT="AA")=N_SAMPLES || COUNT(GT="RR")=N_SAMPLES' $maxmiss_50_filtered.recode.vcf.gz
bcftools view -Oz -o $maxmiss_80_filtered.nomono.vcf.gz -e 'COUNT(GT="AA")=N_SAMPLES || COUNT(GT="RR")=N_SAMPLES' $maxmiss_80_filtered.recode.vcf.gz
bcftools view -Oz -o $maxmiss_100_filtered.nomono.vcf.gz -e 'COUNT(GT="AA")=N_SAMPLES || COUNT(GT="RR")=N_SAMPLES' $maxmiss_100_filtered.recode.vcf.gz

bcftools view -Oz -o $maxmiss_50_filtered.mac2.nomono.vcf.gz -e 'COUNT(GT="AA")=N_SAMPLES || COUNT(GT="RR")=N_SAMPLES' $maxmiss_50_filtered.mac2.recode.vcf.gz
bcftools view -Oz -o $maxmiss_80_filtered.mac2.nomono.vcf.gz -e 'COUNT(GT="AA")=N_SAMPLES || COUNT(GT="RR")=N_SAMPLES' $maxmiss_80_filtered.mac2.recode.vcf.gz
bcftools view -Oz -o $maxmiss_100_filtered.mac2.nomono.vcf.gz -e 'COUNT(GT="AA")=N_SAMPLES || COUNT(GT="RR")=N_SAMPLES' $maxmiss_100_filtered.mac2.recode.vcf.gz

## Prepare a msa in fasta format to build the tree
#bgzip ${maxmiss_100_filtered}.recode.vcf
#tabix ${maxmiss_100_filtered}.recode.vcf.gz
#
#samples=($(zcat ${maxmiss_100_filtered}.recode.vcf.gz | grep -v "^##" | grep "^#" | cut -f10- | tr '\t' '\n'))
#for sample in ${samples[@]}
# do cat $reference | bcftools consensus -s $sample ${maxmiss_100_filtered}.recode.vcf.gz | sed -E "s/^>.*/>$sample/" >> $WORKING_DIR/mtDNA_fasta/mtDNA_msa.fa
# done
