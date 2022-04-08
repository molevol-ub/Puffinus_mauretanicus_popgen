#!/bin/bash
###$ -l h_vmem=50G
#$ -R y
#$ -N master_dedup_mtDNA
#$ -e logs/master_dedup_mtDNA.err.txt
#$ -o logs/master_dedup_mtDNA.out.txt
#$ -cwd
#$ -t 1-908
#$ -q h13.q
#$ -pe ompi255h13 2

### Script to mark duplicates and remove duplicates

## Load your system modules
# Required modules are: picard tools

## Set your master path
MASTER=/users-d3/jferrer/pmau_popgen/SNP_calling
PICARD=/users-d3/jferrer/programari/anaconda3/pkgs/picard-2.20.4-0/share/picard-2.20.4-0

## Fill in directories if different from the workspace setup
bam_in=$MASTER/mt_bams/interim_bams
bam_out=$MASTER/mt_bams/clean_bams

## Fill in path for population specific metadata
metadata=$MASTER/metadata/pmau_popgen_metadata.tsv

while read line
	do
		
		## In array ##
		insampleID_array=$(echo $line | cut -d " " -f 2 )
		insampleID=$bam_in/$insampleID_array
		
		## Out array
		outsampleID_array=$(echo $line | cut -d " " -f 2 )
		outsampleID=$bam_out/$outsampleID_array
		
		## Run picardtools MarkDuplicates
		java -jar $PICARD/picard.jar MarkDuplicates \
		I=$insampleID.sorted.rg.bam \
		O=$outsampleID.sorted.dups.bam \
		METRICS_FILE=$outsampleID.metrics.txt \
		REMOVE_DUPLICATES=true \
		VALIDATION_STRINGENCY=LENIENT AS=true \
		MAX_FILE_HANDLES=1000
		
		## Index the deduped bam files
		java -jar $PICARD/picard.jar BuildBamIndex \
		I=$outsampleID.sorted.dups.bam VALIDATION_STRINGENCY=LENIENT
	done < $metadata
