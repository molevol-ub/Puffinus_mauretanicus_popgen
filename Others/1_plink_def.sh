#!/bin/bash
#$ -cwd
#$ -R y
#$ -e plink.err
#$ -o plink.out
#$ -q h13.q
#$ -pe ompi255h13 4
#$ -V                    #export environment var
#$ -N plink             #name Job

# 1 Make "plink" directory and define your in_file

#mkdir /users-d3/jferrer/gizquierdo/TFG/chr_split/unmasked/vcfs/auto/noLD_plink
#cd /users-d3/jferrer/gizquierdo/calonectris_Joan

VCF=./indels_scfD.recode.vcf

# 2 Perform linkage pruning (identify prune sites)

#plink --vcf $VCF --double-id --allow-extra-chr --set-missing-var-ids @:# -indep-pairwise 50 10 1 --vcf-half-call m --out indels_scfE

# 3 Prune and create PCA

plink --vcf $VCF --double-id --allow-extra-chr --set-missing-var-ids @:# --vcf-half-call m --make-bed --pca --out indels_scfD

