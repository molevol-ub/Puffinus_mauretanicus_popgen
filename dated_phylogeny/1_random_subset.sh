#!/bin/bash
#$ -cwd
#$ -R y
#$ -e fixed_diff.err
#$ -o fixed_diff.out
#$ -q h13.q
#$ -pe ompi255h13 4
#$ -V                    #export environment var
#$ -N fixed_diff            #name Job
#$ -M 000izquierdoguillem@gmail.com
#$ -m be

#bcftools view Puffinus_downsized.wPP.maxmiss80.nomono.vcf.gz | vcfrandomsample -r 0.00012645 > Puffinus_downsized.wPP.maxmiss80.subset1.vcf
#bcftools view Puffinus_downsized.wPP.maxmiss80.nomono.vcf.gz | vcfrandomsample -r 0.00012645 > Puffinus_downsized.wPP.maxmiss80.subset2.vcf
#bcftools view Puffinus_downsized.wPP.maxmiss80.nomono.vcf.gz | vcfrandomsample -r 0.00012645 > Puffinus_downsized.wPP.maxmiss80.subset3.vcf
#bcftools view Puffinus_downsized.wPP.maxmiss80.nomono.vcf.gz | vcfrandomsample -r 0.00012645 > Puffinus_downsized.wPP.maxmiss80.subset4.vcf
#bcftools view Puffinus_downsized.wPP.maxmiss80.nomono.vcf.gz | vcfrandomsample -r 0.00012645 > Puffinus_downsized.wPP.maxmiss80.subset5.vcf

#vcftools --gzvcf Puffinus_downsized.wPP.maxmiss80.nomono.vcf.gz --max-missing 1 --recode --out Puffinus_downsized.wPP.maxmiss100
#mv Puffinus_downsized.wPP.maxmiss100.recode.vcf Puffinus_downsized.wPP.maxmiss100.nomono.vcf

#bgzip Puffinus_downsized.wPP.maxmiss100.nomono.vcf
#tabix Puffinus_downsized.wPP.maxmiss100.nomono.vcf.gz

bcftools view Puffinus_downsized.wPP.maxmiss100.nomono.vcf.gz | vcfrandomsample -r 0.00093691 > Puffinus_downsized.wPP.maxmiss100.subset6.vcf 
bcftools view Puffinus_downsized.wPP.maxmiss100.nomono.vcf.gz | vcfrandomsample -r 0.00093691 > Puffinus_downsized.wPP.maxmiss100.subset7.vcf
bcftools view Puffinus_downsized.wPP.maxmiss100.nomono.vcf.gz | vcfrandomsample -r 0.00093691 > Puffinus_downsized.wPP.maxmiss100.subset8.vcf
bcftools view Puffinus_downsized.wPP.maxmiss100.nomono.vcf.gz | vcfrandomsample -r 0.00093691 > Puffinus_downsized.wPP.maxmiss100.subset9.vcf
bcftools view Puffinus_downsized.wPP.maxmiss100.nomono.vcf.gz | vcfrandomsample -r 0.00093691 > Puffinus_downsized.wPP.maxmiss100.subset10.vcf

