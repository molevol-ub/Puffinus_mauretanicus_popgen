#!/bin/bash
#$ -cwd
#$ -V

# Script to extract only those scaffolds corresponding to the Z chromosome of Spheniscus

# 1 # Define general pathways and files

guillem=/set/your/path
programari=/path/to/bin
pmau_genome=$guillem/genome_masked.fa

# 2 # Index the genome
samtools faidx $pmau_genome

# 3 # Define all files

sphen_chr=(chrz chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chr23 chr24 chr25)

# 4 # Run blast db for all

for chr in "${sphen_chr[@]}"
	do
		sphen_genome=$guillem/chr_split/sphen_genome/input_files_db/${chr}_sphen
		
		$programari/blastn -query $pmau_genome -db $sphen_genome -evalue 1e-20 -max_target_seqs 20 -max_hsps 1000 -out $guillem/chr_split/sphen_genome/${chr}_sphen.1000.txt -outfmt 6 -num_threads 12
	done
# Search meaning of options: -evalue. For the moment we leave -evalue 1e-3

