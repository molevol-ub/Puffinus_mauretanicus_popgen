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

# 3. Keep only sex-linked scaffolds in Ppuf.chrZ.fa, while remove them from Ppuf.fa - afterwards, combine them into one single file: Ppuf.def.fa
# To do so, 1st you have to list scaffold names in the chrZ VCF

vcftools --gzvcf COP1.chrZ.vcf.gz --kept-sites --out COP1.chrZ
cat COP1.chrZ.kept.sites | cut -f 1 | uniq > prova
mv prova COP1.chrZ.kept.sites

#Include these from chrZ fasta

awk '{ if ((NR>1)&&($0~/^>/)) { printf("\n%s", $0); } else if (NR==1) { printf("%s", $0); } else { printf("\t%s", $0); } }' Ppuf.chrZ.fa | grep -Ff COP1.chrZ.kept.sites - | tr "\t" "\n" > prova
mv prova Ppuf.def2.fa

#Exclude from autosomal fasta

awk '{ if ((NR>1)&&($0~/^>/)) { printf("\n%s", $0); } else if (NR==1) { printf("%s", $0); } else { printf("\t%s", $0); } }' Ppuf.fa | grep -v -Ff COP1.chrZ.kept.sites - | tr "\t" "\n" > prova
mv prova Ppuf.def.fa

# Concatenate files

cat Ppuf.def2.fa >> Ppuf.def.fa
rm Ppuf.def2.fa

#Eliminate any white spaces that might have appeared by error - later remove estra characters to end up with just the scaffold

sed 's/ /_/g' Ppuf.def.fa > prova
mv prova Ppuf.def.fa

cat Ppuf.def.fa | sed 's/:[1234567890-]\+//g' | sed 's/[1234567890]\+_//g' > prova
mv prova Ppuf.def.fa

# 4. Proceed to soft mask these genomes
# 1st obtain the soft masked positions in the reference genom using a custom script

soft_pos_script=/users-d3/jferrer/gizquierdo/TFM/Cactus/scripts/generate_masked_bed.py
python $soft_pos_script /users-d3/jferrer/pmau_popgen/genome/genome.fa soft_masked.bed

# Now soft mask

bedtools maskfasta -fi Ppuf.def.fa -bed soft_masked.bed -soft > Ppuf.masked.fa

#---------------------------------------------------------------------
#---------------PREPARE  THE  CACTUS  363-WAY  ALLIGNMENT------------
#---------------------------------------------------------------------

# 1. Download the 363-way allignment

wget -P /media/guillem/BC90A1CD90A18F08/Guillem/TFM_GIA/Cactus https://cgl.gi.ucsc.edu/data/cactus/363-avian-2020-hub/Gallus_gallus/Gallus_gallus.maf.gz
cd /your/dir

# 2. Convert MAF.gz to FASTA files for each genome (the script maf2fasta.pl is required)

gunzip Gallus_gallus.maf.gz
perl /users-d3/jferrer/gizquierdo/TFM/Cactus/scripts/maf2fasta.pl < Gallus_gallus.maf > Gallus_gallus.fasta
#rm Gallus_gallus.maf  if gunzipped file still present


# Remember to run with Singularity instead of Docker

--binariesMode singularity
