#!/bin/bash
#$ -cwd
#$ -R y
#$ -e roh_plink.err
#$ -o roh_plink.out
#$ -q h13.q
#$ -pe ompi255h13 12
#$ -V                    #export environment var
#$ -N ROH_plink             #name Job

# Script to calculate ROHs using the parameters in the kakapo paper: https://www.sciencedirect.com/science/article/pii/S2666979X21000021

# Create and acess the directory were you'll be working

cd /users-d3/jferrer/gizquierdo/TFM/conservation_genomics/ROHs/plink

VCF_DIR=/users-d3/jferrer/gizquierdo/TFM/conservation_genomics/het/PopGenome/DEF_HET		

# 1. First merge all individual vcfs - the 1st one is manual

cp $VCF_DIR/ALT78.masked.vcf.gz .
mv ALT78.masked.vcf.gz Puffinus_ROHs.vcf.gz
tabix Puffinus_ROHs.vcf.gz

list=(CZA11 ILA13 ILA2 PORQ TZE1 M8 M1 M5 M20 M4 G3 M2 M3 M18 M13 M19 M14 M12 M11 M17 M10 M21 G12 G9 G10 G11 G14 G15 M16 G4 M6 COP1 LT2)

for ind in ${list[*]}
do

VCF=$VCF_DIR/$ind.masked.vcf.gz
bcftools merge -o prova.vcf --threads 12 Puffinus_ROHs.vcf.gz $VCF

bgzip prova.vcf
tabix prova.vcf.gz
mv prova.vcf.gz Puffinus_ROHs.vcf.gz

rm Puffinus_ROHs.vcf.gz.tbi
tabix Puffinus_ROHs.vcf.gz

done

name=Puffinus_std
	
# 2. Run ROH calculation
	
plink --vcf Puffinus_ROHs.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --homozyg --out $name
