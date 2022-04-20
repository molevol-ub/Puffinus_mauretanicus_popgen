#!/bin/bash
#$ -cwd
#$ -R y
#$ -e logs/8_bamcaller.err.txt
#$ -o logs/8_bamcaller.out.txt
#$ -q h13.q
#$ -pe ompi255h13 24
#$ -V                    #export environment var
#$ -N bamcaller
#$ -t 1-2
#$ -tc 20

# Script to run bamCaller (masks under a certain coverage) for all scaffolds above 50kbp - separation of these scaffolds done previously with a python script ("scf_lengths.csv")

## Set your master path
MASTER=/users-d3/jferrer/MSMC2_pmau
REF=$MASTER/input/genome.fa
bamcaller=/users-d3/jferrer/programari/msmc-tools
guillem=/users-d3/jferrer/gizquierdo/TFG/genome_scans

## Read bams and chromosomes ##
bam=$MASTER/bam/allgenome_bwa_sorted.bam

## Estimate the average coverage from chromosome SUPER_1 (?)

depth=$(samtools depth -r scf7180000012602 $bam | awk '{sum += $3} END {print sum / NR}')

## Call SNP and generate mask and vcf files

while read line
	do
		chroms=($(echo $line| cut -f 1))
		samtools mpileup -B -q 20 -Q 20 -C 50 -g -r $chroms -f $REF $bam | bcftools call -c -V indels |\
		$bamcaller/bamCaller.py $depth $guillem/masking/$chroms.mask.bed.gz | gzip -c > $guillem/masking/$chroms.vcf.gz
	done < $guillem/vcfs/scf_length.csv

## Unite all mask files under a single 

touch $guillem/masking/scaff_coverage.mask.bed.gz

while read line
	do
		chroms=($(echo $line| cut -f 1))
		zcat $guillem/masking/$chroms.mask.bed.gz >> $guillem/masking/scaff_coverage.mask.bed.gz
	done < $guillem/vcfs/scf_length.csv

