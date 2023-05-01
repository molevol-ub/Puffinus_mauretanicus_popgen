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

# 1. Define the population with which you are going to work and enter its directory

POP=Menorca_noMenCab
cd /users-d3/jferrer/gizquierdo/TFM/chromopainter/$POP/fastGLOBETROTTER

