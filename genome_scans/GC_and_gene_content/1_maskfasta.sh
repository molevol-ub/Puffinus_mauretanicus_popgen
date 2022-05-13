#!/bin/bash
#$ -cwd
#$ -R y
#$ -e mask_bedtools.err
#$ -o mask_bedtools.out
#$ -q h12.q
#$ -pe ompi128h12 2
#$ -V                    #export environment var
#$ -N bedtools_mask             #name Job

# Bash script to calculate GC content using bedtools - don't use redux file (that for autosomes?)

WORK_DIR=/users-d3/jferrer/gizquierdo/TFG/genome_scans/masking
mask_file=$WORK_DIR/final_masked.redux.bed
outfile=$WORK_DIR/genome_masked.def.fa
genome_file=/users-d3/jferrer/pmau_popgen/genome/genome.fa

# 1 Use bedtools intersect to obtain a bedfile of the MASKED regions (instead of unmasked) which you need for bedtools maskfasta

true_mask=$WORK_DIR/masked_regions_genome.bed
genome_bed=/users-d3/jferrer/gizquierdo/TFG/chr_split/scf_info/scf_info.txt

cat $genome_bed | cut -f 1,2 | sed '1d' > $WORK_DIR/temp	#Prepares genome_bed
python3 insert_zeros.py $WORK_DIR/temp

cat $mask_file | cut -f 1,2,3 >  $WORK_DIR/temp2			#Prepares mask_file

bedtools subtract -a $WORK_DIR/temp -b $WORK_DIR/temp2 -bed > $true_mask

# 2 Run bedtools maskfasta

bedtools maskfasta -fi $genome_file -bed $true_mask -fo $outfile

rm $WORK_DIR/temp*
