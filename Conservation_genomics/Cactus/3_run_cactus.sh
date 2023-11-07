#!/bin/bash
#$ -M 000izquierdoguillem@gmail.com
#$ -m be

# We add our genomes to the 363-way allignment by solving 2 subproblems
# More info here

date

export PATH=$PATH:~/programari/cactus-bin-v2.5.1/bin
source activate py39
cd /users-d3/jferrer/gizquierdo/TFM/Cactus/Cactus

# Remember to run with Singularity instead of Docker

mkdir LOGS
mkdir WDIR

cactus jobstore Puffinus_run1.txt Puffinus_run1.hal --root birdAnc361 \
        --maxCores 10 \
        --logFile Cactus.log \
        --binariesMode local \
        --writeLogs LOGS/ \
        --workDir WDIR \
        --disableCaching \
        --defaultDisk 8G \
        --maxDisk 50G \
        --nodeStorage 20 \
        --consCores 10 \
        --defaultMemory 8G \
        --maxLocalJobs 10 \
        --maxServiceJobs 10 

date
