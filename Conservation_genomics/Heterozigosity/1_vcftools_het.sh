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
outname=vcftools_het.noPP.maxmiss100

vcftools --gzvcf $VCF --het --out $outname
