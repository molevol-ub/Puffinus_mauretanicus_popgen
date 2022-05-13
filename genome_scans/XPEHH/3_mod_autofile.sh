#!/bin/bash
#$ -cwd
#$ -R y
#$ -e selscan.err
#$ -o selscan.out
#$ -q h12.q
#$ -pe ompi128h12 2
#$ -V                    #export environment var
#$ -N xpehh_avg             #name Job

# Script to create mapfile - later we can run selscan

window_size=25

WORKDIR=/users-d3/jferrer/gizquierdo/TFG/genome_scans/selscan/output

infile=/users-d3/jferrer/gizquierdo/TFG/genome_scans/PopGenome/Z/$window_size\_windows/Z_statistics.gc.unmasked.csv
window_file=/users-d3/jferrer/gizquierdo/TFG/genome_scans/PopGenome/Z/$window_size\_windows/$window_size\_windows.Z.bed			# Take care which window file to use!!
selscan_file=$WORKDIR/Puffinus_selscan_Z.out

prov_file=$WORKDIR/provisonal.bed

# Keep only position and XP_EHH columns and run bedtools intersect with window_file

cat $selscan_file | cut -f 1 > col1
cat $selscan_file | cut -f 2 > col2
cat $selscan_file | cut -f 8 > col8

paste -d "\t" col1 col2 col2 col8 > $prov_file
rm col*

bedtools intersect -wa -wb -a $window_file -b $prov_file > temp
mv temp $prov_file

# Run python script to make the average

python3 average_xpehh.py $prov_file $infile

mv temp $WORKDIR/Z_selscan.$window_size.csv
rm xpehh
