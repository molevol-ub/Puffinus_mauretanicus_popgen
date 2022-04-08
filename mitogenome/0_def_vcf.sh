#!/bin/bash

# Bash script to remove the ending and starting sites of the mitogenome vcf - we only use sites between 200 and 15600 bp
# These sites were decided according to the mapping quality profiles visualized using IGV v.2.12.2

MASTER=/users-d3/jferrer/pmau_popgen
reference=/$MASTER/SNP_calling/reference/Pmau_mitogenome_200_15600.fasta
vcf=$MASTER/SNP_calling/mt_vcfs/wPP/Puffinus_mtDNA_SNP.maxmiss80.filtered.nomono.vcf.gz
work_dir=/users-d3/jferrer/gizquierdo/TFG/mitogenome/def_vcfs

# Filter out the positions

vcftools --gzvcf $vcf --out $work_dir/Puffinus_mtDNA_SNP.maxmiss80.filtered.nomono.def --chr Contig01+0404545312 --from-bp 200 --to-bp 15600 

bgzip $work_dir/Puffinus_mtDNA_SNP.maxmiss80.filtered.nomono.def.vcf
tabix $work_dir/Puffinus_mtDNA_SNP.maxmiss80.filtered.nomono.def.vcf.gz 

#---------------------------------------------------------------------------------------------------------------------------------------

# Also use this to run a PCA inside your plink directory

mkdir ../plink_directory
cd ../plink_directory
VCF=/users-d3/jferrer/pmau_popgen/SNP_calling/mt_vcfs/noPP/Puffinus_mtDNA_SNP.maxmiss80.filtered.mac2.nomono.vcf.gz

plink --vcf $VCF --double-id --allow-extra-chr --set-missing-var-ids @:# --make-bed --pca --out output_name
