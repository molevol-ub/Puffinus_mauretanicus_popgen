#!/usr/bin/bash

# 1. Use bcftools to intersect the VCF files resulting from the GATK and Freebayes SNP calling pipelines

bcftools isec -p out_dir -n=2 vcf_gatk.vcf.gz vcf_freebayes.vcf.gz

# 2. Rename the output file

mv out_dir/0000.vcf.gz Puffinus_intersect.raw.vcf.gz
tabix Puffinus_intersect.raw.vcf.gz
