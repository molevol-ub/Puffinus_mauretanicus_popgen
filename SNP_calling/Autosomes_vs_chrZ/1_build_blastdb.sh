#!/bin/bash
#$ -V
#$ -cwd

# Script to build a blast database for both, the Spheniscus humboldtii genome and the Puffinus genome using makeblastdb according to https://www.ncbi.nlm.nih.gov/books/NBK569841/

# 1 # Set up main directories and files

guillem=/your/path/sphen_genome
programari=/path/to/bin
outfiles=$guillem/input_files_db

# 2 # Define all files

sphen_chr=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chr23 chr24 chr25)

# 3 # Run makeblast db for all

for chr in "${sphen_chr[@]}"
	do
		sphen_genome=$guillem/${chr}_sphen.fasta
		echo "Fasta file opened successfully"
		$programari/makeblastdb -in $sphen_genome -dbtype nucl -out $outfiles/${chr}_sphen
	done
