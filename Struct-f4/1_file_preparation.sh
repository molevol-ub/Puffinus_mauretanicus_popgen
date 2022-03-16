# Script to generate the input "treemix-like" file needed to run Calc-f4

# 1 Define directories and files

VCF=/users-d3/jferrer/pmau_popgen/SNP_calling/mt_vcfs/noPP/Puffinus_mtDNA_SNP.maxmiss80.filtered.nomono.vcf.gz
work_dir=/users-d3/jferrer/gizquierdo/TFG/struct-f4/mitogenome
output=pmau.80.mito

# 2 Use PLINK to exclude linkage-desequilibrium sites and generate tped and tfam files

#plink --vcf $VCF --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 10 0.1 --recode --transpose --out $work_dir/$output

# Just in this case for mitogenome, we don't filter for linkage desequilibrium

plink --vcf $VCF --double-id --allow-extra-chr --set-missing-var-ids @:# --recode --transpose --out $work_dir/$output

# 3 Run Pablo's Perl script to generate the new treemix-like files

Tped2Structf4.pl $work_dir/$output.tped $work_dir/$output.tfam 5000000 $work_dir/treemix_files/$output 0
