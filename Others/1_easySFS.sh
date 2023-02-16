#!/bin/bash
#$ -cwd
#$ -R y
#$ -e logs/easySFS.err
#$ -o logs/easySFS.out
#$ -q h10.q
#$ -pe ompi64h10 2
#$ -V                    #export environment var
#$ -N easySFS             #name Job

source activate whatshap
export PATH=/users-d3/jferrer/gizquierdo/TFM/SFS/scripts:$PATH

VCF=/users-d3/jferrer/gizquierdo/TFG/chr_split/vcfs/auto/Puffinus_SNP.maxmiss100.filtered.noPP.merged.nomono.masked.vcf.gz
popfile=/users-d3/jferrer/gizquierdo/TFM/SFS/pops_file.txt

#easySFS.py -a -i $VCF -p $popfile --preview

easySFS.py -a -i $VCF -p $popfile --proj 12,12 -o /users-d3/jferrer/gizquierdo/TFM/SFS/Pmau_auto_100.4 --prefix Pmau_Pyel
#easySFS.py -a -i $VCF -p $popfile --proj 22,8 -o /users-d3/jferrer/gizquierdo/TFM/SFS/Z_65 --prefix Puffinus_SFS.Z.65

