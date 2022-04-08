#!/bin/bash
###$ -l h_vmem=50G
#$ -R y
#$ -N master_read_groups_freebayes
#$ -e logs/master_read_groups_freebayes.err.txt
#$ -o logs/master_read_groups_freebayes.out.txt
#$ -cwd
#$ -t 1-60
#$ -q h0809.q
#$ -pe ompi24 4

### Script to add readgroups to bams

## Load your system modules
# Required modules are picard tools

## Load your system modules
# Required modules are: picard tools

## Set your master path
MASTER=/users-d3/jferrer/pmau_popgen/SNP_calling
PICARD=/users-d3/jferrer/programari/anaconda3/pkgs/picard-2.20.4-0/share/picard-2.20.4-0

## Fill in directories if different from the workspace setup
bam_in=$MASTER/mt_bams/clean_bams
bam_out=$MASTER/mt_bams/freebayes_bams

## Fill in path for population specific metadata
metadata=$MASTER/metadata/pmau_popgen_metadata.tsv

# Metadata file for each population
while read line
	do
		ID=$(echo $line | cut -d " " -f 2 )
		
		LB=$(echo $line | cut -d " " -f 2 )
		
		PL="ILLUMINA"
		
		SM=$(echo $line | cut -d " " -f 2 )
		
		## In array
		insampleID_array=$(echo $line | cut -d " " -f 2 )
		insampleID=$bam_in/$insampleID_array
		
		## Out array
		outsampleID_array=$(echo $line | cut -d " " -f 2 )
		outsampleID=$bam_out/$outsampleID_array
		
		## Run picard tools AddreplaceRGs
		java -jar $PICARD/picard.jar AddOrReplaceReadGroups \
		I=${insampleID}.sorted.dups.bam \
		O=${outsampleID}.merged.sorted.dups.freebayesrg.bam \
		RGSM=${SM} \
		RGLB=${LB} \
		RGID=${ID} \
		RGPU=${ID} \
		RGPL=${PL} \
		
		## Index the readgroup bam files
		java -jar $PICARD/picard.jar BuildBamIndex \
		I=${outsampleID}.merged.sorted.dups.freebayesrg.bam VALIDATION_STRINGENCY=LENIENT
	done < $metadata


