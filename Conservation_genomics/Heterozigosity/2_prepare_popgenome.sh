#!/bin/bash
#$ -cwd
#$ -R y
#$ -e het.err
#$ -o het.out
#$ -q h13.q
#$ -pe ompi255h13 2
#$ -V                    #export environment var
#$ -N vcftools.het             #name Job

cd /users-d3/jferrer/gizquierdo/TFM/conservation_genomics/het/PopGenome

VCF=/users-d3/jferrer/gizquierdo/TFG/chr_split/vcfs/auto/Puffinus_SNP.maxmiss80.filtered.wPP.merged.nomono.masked.vcf.gz

list=(ALT78 CZA11 ILA13 ILA2 PORQ TZE1 M8 M1 M5 M20 M4 G3 M2 M3 M18 M13 M19 M14 M12 M11 M17 M10 M21 G12 G9 G10 G11 G14 G15 M16 G4 M6 COP1 LT2)

# To correctly run the PopGenome analyses, we have to repeat the samples twice

for ind in ${list[*]}
do

vcftools --gzvcf $VCF --indv $ind --recode --out $ind.maxmiss80.masked
mv $ind.maxmiss80.masked.recode.vcf $ind.maxmiss80.masked.vcf

bgzip $ind.maxmiss80.masked.vcf
tabix $ind.maxmiss80.masked.vcf.gz

bcftools view -Oz -o $ind.maxmiss80.masked.nomono.vcf.gz -e 'COUNT(GT="AA")=N_SAMPLES || COUNT(GT="RR")=N_SAMPLES' $ind.maxmiss80.masked.vcf.gz
tabix $ind.maxmiss80.masked.nomono.vcf.gz

vcftools --gzvcf $ind.maxmiss80.masked.nomono.vcf.gz --window-pi 25000 --window-pi-step 25000 --out $ind

done
