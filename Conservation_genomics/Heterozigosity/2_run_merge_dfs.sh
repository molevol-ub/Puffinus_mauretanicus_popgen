#!/bin/bash
#$ -cwd
#$ -R y
#$ -e covgenome.err
#$ -o covgenome.out
#$ -q h12.q
#$ -pe ompi128h12 2
#$ -V                    #export environment var
#$ -N covgenome             #name Job

# Script must be run in hercules12 to use the correct R version

cd /users-d3/jferrer/gizquierdo/TFM/conservation_genomics/het/PopGenome

list=(ALT78 CZA11 ILA13 ILA2 PORQ TZE1 M8 M1 M5 M20 M4 G3 M2 M3 M18 M13 M19 M14 M12 M11 M17 M10 M21 G12 G9 G10 G11 G14 G15 M16 G4 M6 COP1 LT2)

for ind in ${list[*]}
do

/soft/R-4.1.1/bin/Rscript /users-d3/jferrer/gizquierdo/TFM/conservation_genomics/het/scripts/merge_dfs.R $ind

done

