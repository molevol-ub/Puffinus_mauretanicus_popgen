#!/bin/bash
#$ -cwd
#$ -R y
#$ -e logs/easySFS.err
#$ -o logs/easySFS.out
#$ -q h10.q
#$ -pe ompi64h10 2
#$ -V                    #export environment var
#$ -N easySFS             #name Job

# Use easySFS.py to create folded SFS using a the projection into the x,y desired haplotype indvs (--proj).

source activate whatshap
export PATH=/users-d3/jferrer/gizquierdo/TFM/SFS/scripts:$PATH

VCF=/users-d3/jferrer/gizquierdo/TFG/chr_split/vcfs/auto/Puffinus_SNP.maxmiss100.filtered.noPP.merged.nomono.masked.vcf.gz
popfile=/users-d3/jferrer/gizquierdo/TFM/SFS/pops_file.txt         # Popfile containing the indvs of each population in the vcf

#easySFS.py -a -i $VCF -p $popfile --preview               # Use to choose best projection to use when dealing with missing data

easySFS.py -a -i $VCF -p $popfile --proj 12,12 -o /users-d3/jferrer/gizquierdo/TFM/SFS/Pmau_auto_100.4 --prefix Pmau_Pyel

