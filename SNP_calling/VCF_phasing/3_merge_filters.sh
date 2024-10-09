#!/bin/bash
#$ -cwd
#$ -R y

# Script to leave only those SNPs included in the filtered VCF of interest (e.g. completeness=80% & noPP)

# 0. Activate the conda environment where whatshap was installed

outdir=/path/to/phased_vcfs/autosomes

scfdir=$outdir/temp_vcfs
rm -r $scfdir
mkdir $scfdir

phased_vcf=$outdir/Puffinus_raw.auto.shapeit4_whatshap_phased.vcf
unphased_vcf=/path/to/vcfs/auto/Puffinus_SNP.maxmiss80.filtered.noPP.merged.nomono.masked.vcf.gz
output_vcf=$outdir/Puffinus_def.auto.shapeit4_whatshap_phased.nPP.maxmiss80.masked

kept_sites=$scfdir/kept_sites.100.noPP

# 1. Use bedtools intersect to keep only those sites present in both vcfs (they should be all those present in th unphased vcf) -previouls uncompress gzvcfs; later compress them

gunzip $phased_vcf.gz
gunzip $unphased_vcf.gz

bedtools intersect -a $unphased_vcf -b $phased_vcf -bed > $kept_sites

bgzip $phased_vcf
bgzip $unphased_vcf

# 2. Keep only the 2 first columns of the output

cat $kept_sites | cut -f 1,2 > $kept_sites.prov
mv $kept_sites.prov $kept_sites

# 3. Use vcftools to keep only these SNPs on the final vcf; also, remove the individuals that shouldn't be present in that vcf (as they were removed from the unphased one)

vcftools --gzvcf $phased_vcf.gz --positions-overlap $kept_sites --remove-indv sacella --remove-indv M22 --remove-indv M7 --remove-indv COP1 --remove-indv LT2 --recode --out $output_vcf

mv $output_vcf.recode.vcf $output_vcf.vcf
bgzip $output_vcf.vcf
tabix $output_vcf.vcf.gz
