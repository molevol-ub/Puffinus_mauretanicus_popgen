#!/bin/bash
#$ -cwd
#$ -R y
#$ -e SNAPP.80.normal.err
#$ -o SNAPP.80.normal.out
#$ -q h13.q
#$ -pe ompi255h13 24
#$ -V                    #export environment var
#$ -N SNAPP.80.normal            #name Job
#$ -M 000izquierdoguillem@gmail.com
#$ -m be

WORKDIR=/users-d3/jferrer/gizquierdo/TFM/dated_phylo/SNAPP/maxmiss80/scfE

/users-d3/jferrer/gizquierdo/beast/bin/beast -threads 24 $WORKDIR/scfE_downsized.maxmiss80.phased.xml

