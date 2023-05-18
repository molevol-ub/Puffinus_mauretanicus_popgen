#!/bin/bash
#$ -cwd
#$ -R y
#$ -e divide_vcf.err
#$ -o divide_vcf.out
#$ -q h11.q
#$ -pe ompi64h11 1
#$ -V                    #export environment var
#$ -N own_pops             #name Job

# Script to modify chromopainter input files in order paint samples according to their own population

workdir=/users-d3/jferrer/gizquierdo/TFM/chromopainter/double
#mkdir $workdir
cd $workdir

# 1. Modifiy the input file by repeating the phased information twice

head Puffinus.auto.nPP.maxmiss80.chromopainter.inp -n 3 > Puffinus_double.chromopainter.inp 
tail Puffinus.auto.nPP.maxmiss80.chromopainter.inp -n +4 >> Puffinus_double.chromopainter.inp
tail Puffinus.auto.nPP.maxmiss80.chromopainter.inp -n +4 >> Puffinus_double.chromopainter.inp

#In the idfile repeat rows at the end, 1st set with ind names, 2nd set with popnames; in the popfile, use the pops as donors, the inds as targets
#nano Puffinus.auto.nPP.maxmiss80.idfile.txt
#nano Puffinus.auto.nPP.maxmiss80.popfile.txt
