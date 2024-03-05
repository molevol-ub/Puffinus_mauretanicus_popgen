#!/bin/bash
#$ -cwd
#$ -R y
#$ -e wig2bed.err
#$ -o wig2bed.out
#$ -q h14.q
#$ -pe ompi511h14 2
#$ -V                    #export environment var
#$ -N phylop             #name Job
#$ -M 000izquierdoguillem@gmail.com
#$ -m be

date

export PATH=$PATH:~/programari/cactus-bin-v2.6.8/bin
export PYTHONPATH=~/programari/cactus-bin-v2.6.8/lib:\$PYTHONPATH
source activate py390
source /users-d3/jferrer/programari/cactus-bin-v2.6.8/venv-cactus-v2.6.8/bin/activate
cd /users-d3/jferrer/gizquierdo/TFM/Cactus/Cactus/PhyloP

# 1. Transform each wig file to bed and append to the same file

ls *.wig > prova.list

while read -r line
do
cat $line | wig2bed >> Puffinus.phyloP.CON.bed
done < prova.list

rm prova.list

# 2. Keep only macro chromosomes and unmasked regions

bedtools intersect -a Puffinus.phyloP.CON.bed -b ~/gizquierdo/TFM/dadi/regions/macro_chr.bed > Puffinus.phyloP.CON.masked.bed
cut -f 1,2,3 /users-d3/jferrer/gizquierdo/TFG/genome_scans/masking/final_masked.bed > final_masked.bed
bedtools intersect -a Puffinus.phyloP.CON.masked.bed -b final_masked.bed > prova
mv prova Puffinus.phyloP.CON.masked.bed
rm final_masked.bed

# 3. Lines have become duplicated, keep only odd numbered ones (or even, idc)

cat Puffinus.phyloP.CON.masked.bed | sed -n 'n;p' > prova.bed
mv prova.bed Puffinus.phyloP.CON.masked.bed

# 4. Split files to lower computational intensity and run script to calculate p.values and q-values
split -l59900453 Puffinus.phyloP.CON.masked.bed Puffinus.phyloP.CON.masked.subsets
for i in aa ab ac ad ae af ag ah ai aj ak al
do
Rscript ../8_qval.r Puffinus.phyloP.CON.masked.subsets$i
cat Puffinus.phyloP.CON.masked.subsets$i.csv >> Puffinus.phyloP.CON.masked.qval.csv
rm *subsets*
done

# 5. Make style modifications
cut -d " " -f 2,3,4,5,6,7,8 Puffinus.phyloP.CON.masked.qval.csv | sed "s/ /\t/g" | grep -v "start" > Puffinus.phyloP.CON.masked.def.csv
mv Puffinus.phyloP.CON.masked.def.csv Puffinus.phyloP.CON.masked.qval.csv

# 6. Select SNPs from the VCF used for ROHs
# And do intersect with the previous files

zless ~/gizquierdo/TFM/conservation_genomics/ROHs/plink/proves/Puffinus_ROHs.vcf.gz | grep -v "#" | cut -f 1,2 | awk -F '\t' '{print $1, $2 - 1, $2}' | sed "s/ /\t/g" > prova.bed
bedtools intersect -a ../../PhyloP/Puffinus.phyloP.CON.masked.qval.csv -b prova.bed > Puffinus.phyloP.CON.masked.qval.SNPs.csv

# 7. Keep only SNPs keep only SNPs below the threshold of your choice (in our case 0.00000001)

#cat Puffinus.phyloP.CON.masked.qval.SNPs.csv | awk -F"\t" '$7<0.00000001' > Puffinus.phyloP.CON.masked.qval.SNPs.threshold1e-8.csv
#cut -f 1,3 Puffinus.phyloP.CON.masked.qval.SNPs.threshold1e-8.csv > conserved_pos.csv

# 8. Make VCF from the ROH VCF with each of the conserved positions

vcftools --gzvcf ~/gizquierdo/TFM/conservation_genomics/ROHs/plink/proves/Puffinus_ROHs.vcf.gz --positions conserved_pos.csv --recode-INFO-all --recode --out Puffinus.SNP.CONS
mv Puffinus.SNP.CONS.recode.vcf Puffinus.SNP.CONS.vcf
bgzip Puffinus.SNP.CONS.vcf
tabix Puffinus.SNP.CONS.vcf.gz

# 9. Get list of those SNPs (CONS.SNPs.bed) and also write a file with the SNP and the reference allele (CONS.SNPs.csv)

zcat Puffinus.SNP.CONS.vcf.gz | grep -v "#" | cut -f 1,2 | awk -F '\t' '{print $1, $2 - 1, $2}' | sed "s/ /\t/g" > CONS.SNPs.bed
zcat Puffinus.SNP.CONS.vcf.gz | grep -v "#" | cut -f 1,2,4 > CONS.SNPs.csv

# 10. For each SNP, convert to maf and count the number of each base pair; write down the major alelel

touch CONS.SNPs.refallele.csv

while read -r line
do

echo $line | sed "s/ /\t/g" > prova.bed

j=$(cut -f1 prova.bed)
i=$(cut -f2 prova.bed)
