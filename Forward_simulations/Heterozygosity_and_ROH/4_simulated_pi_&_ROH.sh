#!/bin/bash
#$ -cwd
#$ -R y

# Script to iterate over every VCF file and individual generated in the forward simulations and calculate heterozygosity and FROH

# To do that you need to previously generate 2 lists and iterate through them: the list of VCFs ("vcf_list3.txt") and that of idividuals ("indv.list")

# 1. First calculate the observed heterozygosity (pi) using vcftools

while read -r line; do while read -r line2; do bgzip $line.vcf; tabix $line.vcf.gz; vcftools --gzvcf $line.vcf.gz --indv $line2 --window-pi 25000000 --out $line2.$line; tail -n +2 $line2.$line.windowed.pi >> $line.windowed.pi ; done < indv.list; done < vcf_list3.txt

# 2. Then infer the length of ROH using plink (for more information see the ROH section in the Conservation_genomics folder)

while read -r line; do plink --vcf $line.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# -indep-pairwise 50 10 0.1 --homozyg-window-snp 100 --homozyg-window-het 1 --homozyg-window-threshold 0.05 --homozyg-snp 25 --homozyg-kb 100 --homozyg-density 50 --homozyg-gap 1000 --homozyg-het 750 --out $line; done < vcf_list3.txt
