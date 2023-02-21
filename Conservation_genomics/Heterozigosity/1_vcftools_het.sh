#!/bin/bash
#$ -cwd
#$ -R y
#$ -e het.err
#$ -o het.out
#$ -q h11.q
#$ -pe ompi64h11 2
#$ -V                    #export environment var
#$ -N vcftools.het             #name Job

cd /users-d3/jferrer/gizquierdo/TFM/conservation_genomics/het

VCF=/users-d3/jferrer/gizquierdo/TFG/chr_split/vcfs/auto/Puffinus_SNP.maxmiss100.filtered.noPP.merged.nomono.masked.vcf.gz
VCF2=/users-d3/jferrer/gizquierdo/TFG/chr_split/vcfs/auto/Puffinus_SNP.maxmiss100.filtered.wPP.merged.nomono.masked.vcf.gz
VCF3=/users-d3/jferrer/gizquierdo/TFG/chr_split/vcfs/auto/Puffinus_SNP.maxmiss100.filtered.noPP.merged.nomono.minmac2.masked.vcf.gz
VCF4=/users-d3/jferrer/gizquierdo/TFG/chr_split/vcfs/auto/Puffinus_SNP.maxmiss80.filtered.noPP.merged.nomono.masked.vcf.gz
VCF5=/users-d3/jferrer/gizquierdo/TFG/chr_split/vcfs/auto/Puffinus_SNP.maxmiss80.filtered.noPP.merged.nomono.minmac2.masked.vcf.gz
VCF6=/users-d3/jferrer/gizquierdo/TFG/chr_split/vcfs/chrZ/Puffinus_SNP.maxmiss65.filtered.noPP.merged.nomono.masked.vcf.gz
VCF7=/users-d3/jferrer/gizquierdo/TFG/chr_split/vcfs/chrZ/Puffinus_SNP.maxmiss65.filtered.noPP.merged.nomono.minmac2.masked.vcf.gz

vcftools --gzvcf $VCF --hardy --out vcftools_hardy.noPP.maxmiss100
vcftools --gzvcf $VCF2 --hardy --out vcftools_hardy.wPP.maxmiss100 
vcftools --gzvcf $VCF3 --hardy --out vcftools_hardy.noPP.maxmiss100.minmac2
vcftools --gzvcf $VCF4 --hardy --out vcftools_hardy.noPP.maxmiss80
vcftools --gzvcf $VCF5 --hardy --out vcftools_hardy.noPP.maxmiss80.minmac2
vcftools --gzvcf $VCF6 --hardy --out vcftools_hardy.chrZ.noPP.maxmiss65
vcftools --gzvcf $VCF7 --hardy --out vcftools_hardy.chrZ.noPP.maxmiss65.minmac2
