#!/bin/bash

#Script to generate the input files needed to run treemix

# Use vcf without missing data and with P.puffinus

vcf_file=/users-d3/jferrer/gizquierdo/TFG/chr_split/vcfs/auto/Puffinus_SNP.maxmiss100.filtered.wPP.merged.nomono.minmac2.vcf.gz
work_dir=/users-d3/jferrer/gizquierdo/TFG/treemix
output=pmau_treemix

# 1 Prune for linkage desequilibrium

plink --vcf $vcf_file --double-id --allow-extra-chr --indep-pairwise 50 10 0.1 --recode vcf --out $work_dir/$output.LDpruned

# 2 Generate the input format needed for running treemix.
# Previously generate the clust file (with the correct individual IDs!) used in this command. It must includes the desired populations in the format: "Sample_ID\tSample_ID\tPop_ID"

vcf2treemix.sh $work_dir/$output.LDpruned.vcf $work_dir/$output.clust
