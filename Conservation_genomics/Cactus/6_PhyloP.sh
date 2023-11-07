#!/bin/bash
#$ -cwd
#$ -R y
#$ -e phylop.err
#$ -o phylop.out
#$ -q h13.q
#$ -pe ompi255h13 1
#$ -V                    #export environment var
#$ -N phylop             #name Job
#$ -t 1-25
#$ -tc 25

date

export PATH=$PATH:~/programari/cactus-bin-v2.6.8/bin
export PYTHONPATH=~/programari/cactus-bin-v2.6.8/lib:\$PYTHONPATH
source activate py390
source /users-d3/jferrer/programari/cactus-bin-v2.6.8/venv-cactus-v2.6.8/bin/activate
cd /users-d3/jferrer/gizquierdo/TFM/Cactus/Cactus/PhyloP

# 1. Create a list of all the scaffolds and their length, ordered by length from bigger to smaller

cat ~/pmau_popgen/genome/genome.fa.fai | cut -f 1,2 | sort -gr -k 2 > Pmau.scf.length

# 2. Split into 24 files of equal length

for i in {0..23}; do awk -v var="$i" 'NR%24==var' Pmau.scf.length > Puffinus_scf.$i.list; done

# 3. For half of them invert the order of the scaffolds: in this way, we won't decompress all the giant mafs at the same time and we will occupy less teras

numlist=("0" "2" "4" "6" "8" "10" "12" "14" "16" "18" "20" "22")
for i in ${numlist[*]}
do
sort -g -k 2 Puffinus_scf.$i.list > prova
mv prova Puffinus_scf.$i.list
done

# And remove second columns

for i in {0..23}; do cut -f 1 Puffinus_scf.$i.list > prova; mv prova Puffinus_scf.$i.list; done
mv Puffinus_scf.0.list Puffinus_scf.24.list

# 4. Then make a loop where you iterate through all scaffold lists:

# Activate the conda environment where whatshap was installed

scf=(0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24)

while read -r scaffold
do

# a) Transform the scaffold to maf

hal2maf ../363-avian-2020.wPmau.hal $scaffold.maf --noAncestors --noDupes --refGenome Pmau --refSequence $scaffold

# b) Run phylop

phyloP --msa-format MAF --mode CON --wig-scores --method LRT ../b10k_model_363_macros.mod $scaffold.maf > $scaffold.wig

# c) Remove the maf because it takes up so much space

rm $scaffold.maf

done < Puffinus_scf.${scf[(($SGE_TASK_ID))]}.list
