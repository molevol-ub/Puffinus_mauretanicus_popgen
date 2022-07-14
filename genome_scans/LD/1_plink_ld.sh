#!/bin/bash
#$ -cwd
#$ -R y
#$ -e plink_LD.err
#$ -o plink_LD.out
#$ -q h13.q
#$ -pe ompi255h13 2
#$ -V                    #export environment var
#$ -N plink_LD             #name Job

# Script to obtain linkage desequilibrium table (r2) using PLINK - for both groups both coupled and separated

# 1 Set up your working directories

scaffold=your_scf

WORKDIR=/your/workdir
mkdir $WORKDIR

indir=/your/indir
infile=$indir/only_$scaffold.recode.vcf.gz 	# Make sure it's already filtered for your scaffold of interest

output=$WORKDIR/only_$scaffold

# 2 Run PLINK - options:	--ld-window [INT] --> ignore comparisons with more than [INT] SNPs between them (set as HIGH as possible)
#				--ld-window-kb [INT] --> ignore comparisons with SNPs more than [INT] kb apart (set as HIGH as possible)
#				--ld-window-r2 [FLOAT] --> ignore comparisons with r2 lower than [FLOAT] (set at 0)
#				square --> output as matrix format

plink --vcf $infile --double-id --allow-extra-chr --r2 square --out $output

# The output file has the following columns: chr1, pos1, snp_id1, chr2, pos2, snp_id2, r2_value

# 3 Cretate file (derived from bim) to use as genome map in the LDhetamap rscript

plink --vcf $infile --double-id --allow-extra-chr --make-just-bim --out $WORKDIR/prova
cat $WORKDIR/prova.bim | cut -f 4 > $output.mapfile
rm $WORKDIR/prova*
