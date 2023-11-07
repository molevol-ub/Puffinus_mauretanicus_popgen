#!/bin/bash
#$ -cwd
#$ -R y
#$ -e SFS.err
#$ -o SFS.out
#$ -q h13.q
#$ -pe ompi255h13 2
#$ -V                    #export environment var
#$ -N halAddToBranch            #name Job
#$ -M 000izquierdoguillem@gmail.com
#$ -m be

# We add our genomes to the 363-way allignment by solving 2 subproblems
# More info here

date

export PATH=$PATH:~/programari/cactus-bin-v2.6.8/bin
export PYTHONPATH=~/programari/cactus-bin-v2.6.8/lib:\$PYTHONPATH
source activate py390
source /users-d3/jferrer/programari/cactus-bin-v2.6.8/venv-cactus-v2.6.8/bin/activate
#cd /users-d3/jferrer/gizquierdo/TFM/Cactus/Cactus

cp 363-avian-2020.hal 363-avian-2020.wPmau.hal

# 1. Make list of renamings and rename using halRenameSequence (repeat this step for every genome that gives problems when running halAddToBranch)
genome=birdAnc274
while read -r genome
do
grep ">" ../genomes/$genome.fa | sed 's/>birdAnc274\.//g' > $genome.names
grep ">" ../genomes/$genome.fa | sed 's/>//g' > $genome.old.names
paste $genome.old.names $genome.names > prova    
mv prova $genome.names 
rm $genome.old.names 
halRenameSequences Puffinus_run1.hal $genome $genome.names
#halRenameSequences Puffinus_run2.hal $genome $genome.names
done < genomes.list

rm *.names

halAddToBranch 363-avian-2020.wPmau.hal Puffinus_run1.hal Puffinus_run2.hal birdAnc274 birdAnc362 Calonectris_borealis Pmau 0.0101 0.0157
date
