#!/bin/bash
###$ -l h_vmem=50G
#$ -R y
#$ -cwd
#$ -t 1-908
#$ -q h0809.q
#$ -pe ompi24 4

### Script to add readgroups to bams

## Load your system modules
# Required modules are picard tools

## Set your master path
PICARD=/users-d3/jferrer/programari/anaconda3/pkgs/picard-2.20.4-0/share/picard-2.20.4-0

## Fill in path for population specific metadata ##
metadata=/users-d3/jferrer/pmau_popgen/SNP_calling/metadata/pmau_popgen_metadata.tsv
bam_in=/users-d3/jferrer/pmau_popgen/SNP_calling/mt_bams/raw_bams
bam_out=/users-d3/jferrer/pmau_popgen/SNP_calling/mt_bams/interim_bams

# Metadata file for each population
# 1 sample_name
# 2 pop_ID
# 3 simple_ID
# 4 sample_well
# 5 ring_ID
# 6 read_1
# 7 read_2
# 8 flowcell
# 9 lane
# 10 library
# 11 barcode
# 12 instrument

while read line
	do
		simpleID=$(echo $line | cut -d " " -f 1 )
		
		instrument=$(echo $line | cut -d " " -f 5 )
		
		seqnum=$(echo $line | cut -d " " -f 9 )
		
		flowcell=$(echo $line | cut -d " " -f 6 )
		
		lane=$(echo $line | cut -d " " -f 7 )
		
		barcode=$(echo $line | cut -d " " -f 8 )
		
		## In array
		insampleID_array=$(echo $line | cut -d " " -f 2 )
		insampleID=$bam_in/$insampleID_array
		
		## Out array
		outsampleID_array=$(echo $line | cut -d " " -f 2 )
		outsampleID=$bam_out/$outsampleID_array
		
		## Run picard tools AddreplaceRGs
		java -jar $PICARD/picard.jar AddOrReplaceReadGroups \
		I=${insampleID}.sorted.raw.bam \
		O=${outsampleID}.sorted.rg.bam \
		RGSM=${simpleID} \
		RGLB=${simpleID}.${seqnum} \
		RGID=${flowcell}.${lane} \
		RGPU=${flowcell}.${lane}.${barcode} \
		RGPL=${instrument} \
		
		## Index the readgroup bam files
		java -jar $PICARD/picard.jar BuildBamIndex \
		I=${outsampleID}.sorted.rg.bam VALIDATION_STRINGENCY=LENIENT
	done < $metadata
