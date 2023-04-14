#!/bin/bash
#$ -cwd
#$ -R y
#$ -e gatk_reference.err
#$ -o gatk_reference.out
#$ -q h13.q
#$ -pe ompi255h13 6
#$ -V                    #export environment var
#$ -N gatk_reference             #name Job
#$ -M 000izquierdoguillem@gmail.com
#$ -m be

# Before anything else, enter the conda environment where you have installed Cactus (in this case with Python version 3.9)

source activate py39
cd /users-d3/jferrer/gizquierdo/TFM/Cactus

#---------------------------------------------------------------------
#---------------CREATE  A  P.puffinus  "REFERENCE"  GENOME------------
#---------------------------------------------------------------------

# Technically you can do this after the cactus allignment as well

cd /users-d3/jferrer/gizquierdo/TFM/Cactus/Ppuf

# 1. Create VCF with only P.puffinus puffinus (you can remove missing data as it is only 1 indv) and remove HOM REF files (but keep HOM ALT) to run gatk faster
# Repeat the step twice, for the autosomal and chrZ scaffolds

#Autosomal

bcftools view -Oz -o COP1.auto.vcf.gz --threads 6 -e 'COUNT(GT="RR")==N_SAMPLES' /users-d3/jferrer/gizquierdo/TFM/conservation_genomics/het/PopGenome/DEF_HET/COP1.masked.vcf.gz
tabix COP1.auto.vcf.gz

#chrZ

vcftools --vcf /users-d3/jferrer/gizquierdo/TFG/chr_split/vcfs/chrZ/Puffinus_SNP.minmax2.mindp4maxdp50.filtered.PP.recode.vcf --indv COP1 --recode-INFO-all --recode --out COP1
mv COP1.recode.vcf COP1.vcf
vcftools --vcf COP1.vcf --max-missing 1.0 --recode-INFO-all --recode --out COP1
mv COP1.recode.vcf COP1.vcf
bgzip COP1.vcf
tabix COP1.vcf.gz

bcftools view -Oz -o COP1.chrZ.vcf.gz --threads 6 -e 'COUNT(GT="RR")==N_SAMPLES' COP1.vcf.gz
tabix COP1.chrZ.vcf.gz

# 2. Use GATK to generate the P.puffinus  

gatk FastaAlternateReferenceMaker -R /users-d3/jferrer/pmau_popgen/genome/genome.fa -O Ppuf.chrZ.fa -V COP1.chrZ.vcf.gz
gatk FastaAlternateReferenceMaker -R /users-d3/jferrer/pmau_popgen/genome/genome.fa -O Ppuf.fa -V COP1.auto.vcf.gz

#---------------------------------------------------------------------
#---------------PREPARE  THE  CACTUS  363-WAY  ALLIGNMENT------------
#---------------------------------------------------------------------

# 1. Download the 363-way allignment

wget -P /media/guillem/BC90A1CD90A18F08/Guillem/TFM_GIA/Cactus https://cgl.gi.ucsc.edu/data/cactus/363-avian-2020-hub/Gallus_gallus/Gallus_gallus.maf.gz




# Remember to run with Singularity instead of Docker

--binariesMode singularity
