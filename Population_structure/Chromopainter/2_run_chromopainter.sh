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

cd /users-d3/jferrer/gizquierdo/TFM/chromopainter

# 1. Estimate the required parameters fro chromopainter using a SUBSAMPLE OF INDIVIDUALS AND SITES; use options:
#      · “-s 0” to specify that no paintingsamples are needed in this step
#      · “-i 10 -in -iM” specifies to use 10 iterations of Expectation-Maximisation to infer the switch (“-in”) and emission (“-iM”) parameters. 

ChromoPainterv2 -g Puffinus_subset.chromopainter.inp -r Puffinus_subset.recomrates.txt -a 0 0 -s 0 -i 10 -in -iM -o Puffinus_subset

# 2. Calculate the average estimated values for n and M to use in the final chromopainter run. The script might need some polishing to use the appropiate input (parameters "$infilePREFIX", "$infileSUFFIX", "@chromovec" and "@chromolengths")
# For more info see: https://github.com/hellenthal-group-UCL/fastGLOBETROTTER/blob/master/tutorial/ChromoPainterv2EstimatedNeMutExtractEM.pl

perl scripts/ChromoPainterv2EstimatedNeMutExtractEM.pl 

# 3. Using the output values of the previous step, run chromopainter
