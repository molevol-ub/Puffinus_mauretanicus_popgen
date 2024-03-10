#!/bin/bash
#$ -cwd
#$ -R y
#$ -e SFS.err
#$ -o SFS.out
#$ -q h12.q
#$ -pe ompi128h12 2
#$ -V                    #export environment var
#$ -N dadi_SFS             #name Job
#$ -M 000izquierdoguillem@gmail.com
#$ -m be

WORKDIR=/users-d3/jferrer/gizquierdo/TFM/dadi/chrZ
cd $WORKDIR
VCF=$WORKDIR/Puffinus_chrZ.gene_masked.vcf.gz

source activate py39

# 1. Genetate the folded SFS for 2 pops

dadi-cli GenerateFs --vcf $VCF --pop-info ../Puffinus.dadi.popfile.txt --pop-ids Pyel Pmau --projections 19 33 --output Puffinus.SFS.1933.fs

# And then assuming all of them belong to the same population

dadi-cli GenerateFs --vcf $VCF --pop-info Pmed.dadi.popfile.txt --pop-ids Pmau --projections 64 --output Pmed.SFS.dadi.fs
