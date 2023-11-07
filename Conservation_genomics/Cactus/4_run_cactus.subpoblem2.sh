#!/bin/bash
#$ -M 000izquierdoguillem@gmail.com
#$ -m be

# We add our genomes to the 363-way allignment by solving 2 subproblems
# More info here

export PATH=$PATH:~/programari/cactus-bin-v2.6.8/bin
export PYTHONPATH=~/programari/cactus-bin-v2.6.8/lib:\$PYTHONPATH
source activate py390
source /users-d3/jferrer/programari/cactus-bin-v2.6.8/venv-cactus-v2.6.8/bin/activate
#cd /users-d3/jferrer/gizquierdo/TFM/Cactus/Cactus

mkdir LOGS2
mkdir WDIR2

cactus jobstore Puffinus_run2.txt Puffinus_run2.hal \
        --maxCores 10 \
        --logFile Cactus.log \
        --binariesMode local \
        --writeLogs LOGS2/ \
        --workDir WDIR2 \
        --disableCaching \
        --defaultDisk 8G \
        --maxDisk 50G \
        --nodeStorage 20 \
        --consCores 10 \
        --defaultMemory 8G \
        --maxLocalJobs 10 \
        --maxServiceJobs 10 

date
