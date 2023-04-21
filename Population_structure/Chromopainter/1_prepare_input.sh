
workdir=/users-d3/jferrer/gizquierdo/TFM/chromopainter
cd $workdir

# 1. First create IMPUTE file from phased VCFs

invcf=/users-d3/jferrer/pmau_popgen/SNP_calling/vcfs/phased_vcfs_def/autosomes/Puffinus_def.auto.shapeit4_whatshap_phased.nPP.maxmiss80.vcf.gz

vcftools --gzvcf $invcf --IMPUTE --out Puffinus.auto.nPP.maxmiss80
