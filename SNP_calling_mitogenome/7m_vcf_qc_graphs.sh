#!/bin/bash

# Script to calculate quality statistics regarding our data (before filtering) - this will allow you to decide which filters to apply when filtering the VCF file

# 1 Create subset VCF

work_dir=/users-d3/jferrer/pmau_popgen/SNP_calling/mt_vcfs/VCF_QC

# 2 Create calculus of allele frequency, mean depth (COVERAGE) per ind & site, site quality, and missing data per ind & site

IN_VCF=/users-d3/jferrer/pmau_popgen/SNP_calling/mt_freebayes/tmp_mtDNA/Puffinus_mitogenome.vcf
OUT=pmau_mitogenome
vcftools --vcf $IN_VCF --freq2 --out $work_dir/$OUT --max-alleles 2
vcftools --vcf $IN_VCF --depth --out $work_dir/$OUT
vcftools --vcf $IN_VCF --site-mean-depth --out $work_dir/$OUT
vcftools --vcf $IN_VCF --site-quality --out $work_dir/$OUT
vcftools --vcf $IN_VCF --missing-indv --out $work_dir/$OUT
vcftools --vcf $IN_VCF --missing-site --out $work_dir/$OUT
