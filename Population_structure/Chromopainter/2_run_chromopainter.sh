#!/bin/bash
#$ -cwd
#$ -R y
#$ -e gatk_reference.err
#$ -o gatk_reference.out
#$ -q h12.q
#$ -pe ompi128h12 4
#$ -V                    #export environment var
#$ -N gatk_reference             #name Job
#$ -M 000izquierdoguillem@gmail.com
#$ -m be

# 1. Estimate the required parameters fro chromopainter using a SUBSAMPLE OF INDIVIDUALS AND SITES; use options:
#      · “-s 0” to specify that no paintingsamples are needed in this step
#      · “-i 10 -in -iM” specifies to use 10 iterations of Expectation-Maximisation to infer the switch (“-in”) and emission (“-iM”) parameters. 

ChromoPainterv2 -g prova.txt -r Puffinus.auto.nPP.maxmiss80.recomrates.txt -a 0 0 -s 0 -i 10 -in -iM -o Puffinus.auto.nPP.maxmiss80
