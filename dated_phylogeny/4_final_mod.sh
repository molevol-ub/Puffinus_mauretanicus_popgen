#!/bin/bash
#$ -cwd
#$ -R y
#$ -e fixed_diff.err
#$ -o fixed_diff.out
#$ -q h13.q
#$ -pe ompi255h13 4
#$ -V                    #export environment var
#$ -N fixed_diff            #name Job
#$ -M 000izquierdoguillem@gmail.com
#$ -m be

# Depending on your beast version, you might have to change the following utility:

for file in ./*.xml
do
sed 's/beast.util.TreeParser/beast.base.evolution.tree.TreeParser/g' $file > $file.def
mv $file.def $file
done
