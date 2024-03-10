#!/bin/bash
#$ -cwd
#$ -R y
#$ -e wig2bed.err
#$ -o wig2bed.out
#$ -q h14.q
#$ -pe ompi511h14 12
#$ -V                    #export environment var
#$ -N SFS_polarized             #name Job
#$ -M 000izquierdoguillem@gmail.com
#$ -m be

date

export PATH=$PATH:~/programari/cactus-bin-v2.6.8/bin
export PYTHONPATH=~/programari/cactus-bin-v2.6.8/lib:\$PYTHONPATH
source activate py390
source /users-d3/jferrer/programari/cactus-bin-v2.6.8/venv-cactus-v2.6.8/bin/activate

VCF=~/gizquierdo/TFM/dadi/regions/Puffinus.gene_masked.macro_masked.vcf.gz

# 1. Make bedfile with each SNP to be considered and the run halBranchMutations to get which ones have to be repolarized according to the ancestor with Calonectris

zcat $VCF | grep -v "#" | cut -f 1,2 | awk -F '\t' '{print $1, $2 - 1, $2}' | sed "s/ /\t/g" > SFS.SNPs.bed
halBranchMutations ~/gizquierdo/TFM/Cactus/Cactus/363-avian-2020.wPmau.hal Pmau --refTargets SFS.SNPs.bed --snpFile SFS.polarized_SNPs.csv

cut -f 1,3 SFS.polarized_SNPs.csv > SFS.polarized_SNPs.pos

# 2. Subset VCFs with REFdel and ALTdel poisitions

vcftools --gzvcf $VCF --exclude-positions SFS.polarized_SNPs.pos --recode-INFO-all --recode --out Puffinus.SNP.ALTdel
mv Puffinus.SNP.ALTdel.recode.vcf Puffinus.SNP.ALTdel.vcf
bgzip Puffinus.SNP.ALTdel.vcf
tabix Puffinus.SNP.ALTdel.vcf.gz

vcftools --gzvcf $VCF --positions SFS.polarized_SNPs.pos --recode-INFO-all --recode --out Puffinus.SNP.REFdel
mv Puffinus.SNP.REFdel.recode.vcf Puffinus.SNP.REFdel.vcf
bgzip Puffinus.SNP.REFdel.vcf
tabix Puffinus.SNP.REFdel.vcf.gz

# 3. Repolarize the REFdel dataset and join both in a single file

gunzip Puffinus.SNP.REFdel.vcf.gz

cat Puffinus.SNP.REFdel.vcf | sed 's/0\/0/0|0/g' > prova
mv prova Puffinus.SNP.REFdel.vcf
cat Puffinus.SNP.REFdel.vcf | sed 's/1\/1/1|1/g' > prova
mv prova Puffinus.SNP.REFdel.vcf
cat Puffinus.SNP.REFdel.vcf | sed 's/1|1/0\/0/g' > prova
mv prova Puffinus.SNP.REFdel.vcf
cat Puffinus.SNP.REFdel.vcf | sed 's/0|0/1\/1/g' > prova
mv prova Puffinus.SNP.REFdel.vcf

bgzip Puffinus.SNP.REFdel.vcf
tabix Puffinus.SNP.REFdel.vcf.gz

bcftools concat -Oz -o Puffinus.SFS.SNP.polarized.vcf.gz -a --threads 12 Puffinus.SNP.REFdel.vcf.gz Puffinus.SNP.ALTdel.vcf.gz
tabix Puffinus.SFS.SNP.polarized.vcf.gz

# Add Ancestral allele (AA) information:

zgrep "#" Puffinus.SFS.SNP.polarized.vcf.gz > Puffinus.SFS.polarized.def.vcf
echo '##INFO=<ID=AA,Number=1,Type=Character,Description="Ancestral allele">' > hdr.txt

zgrep -v "#" Puffinus.SFS.SNP.polarized.vcf.gz | awk -v OFS="\t" -F "\t" '{$8="AA="$4; print}' >> Puffinus.SFS.polarized.def.vcf
bgzip Puffinus.SFS.polarized.def.vcf
tabix Puffinus.SFS.polarized.def.vcf.gz

# 4. Generate SFSs

dadi-cli GenerateFs --vcf Puffinus.SFS.polarized.def.vcf.gz --pop-info ../LDpruned/Puffinus.dadi.popfile.txt --pop-ids Pyel Pmau --projections 24 40 --polarized --output Puffinus.SFS.dadi.fs
dadi-cli GenerateFs --vcf Puffinus.SFS.polarized.def.vcf.gz --pop-info ../LDpruned/Pmed.dadi.popfile.txt --pop-ids Pmau --projections 64 --polarized --output Pmed.SFS.dadi.fs
