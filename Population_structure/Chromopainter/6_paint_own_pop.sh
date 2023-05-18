#!/bin/bash
#$ -cwd
#$ -R y
#$ -e divide_vcf.err
#$ -o divide_vcf.out
#$ -q h11.q
#$ -pe ompi64h11 1
#$ -V                    #export environment var
#$ -N divide_vcf             #name Job
#$ -t 1-876
#$ -tc 36

# Script to modify chromopainter input files in order paint samples according to their own population

workdir=/users-d3/jferrer/gizquierdo/TFM/chromopainter/scaffolds_double
#mkdir $workdir
cd $workdir

# 0. Make your list of scaffolds

scf=($(cat scaffolds.txt))

tail scf7180000013763.chromopainter.inp -n +4
