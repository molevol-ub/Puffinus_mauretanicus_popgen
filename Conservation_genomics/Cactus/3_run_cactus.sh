#!/bin/bash
#$ -cwd
#$ -R y
#$ -e gatk_reference.err
#$ -o gatk_reference.out
#$ -q h13.q
#$ -pe ompi255h13 36
#$ -V                    #export environment var
#$ -N gatk_reference             #name Job
#$ -M 000izquierdoguillem@gmail.com
#$ -m be

# We add our genomes to the 363-way allignment by solving 2 subproblems
# More info here

source activate py39
cd /users-d3/jferrer/gizquierdo/TFM/Cactus/Cactus

# Remember to run with Singularity instead of Docker

cactus jobstore Puffinus_run1.txt Puffinus_run1.hal --root birdAnc361 --binariesMode singularity --consCores 36 --batchSystem slurm --logFile Cactus.log --rotatingLogging

