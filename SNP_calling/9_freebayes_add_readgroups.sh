#!/bin/bash
#$ -t 1-38
#$ -q h13.q
#$ -pe ompi255h13 2

### Script to add readgroups to bams before running freebayes

## Load your system modules
# Required modules are picard tools

## 1. Set your master path
MASTER=/users/jferrer/pmau_popgen
PICARD=/users/jferrer/programari/anaconda3/pkgs/picard-2.20.4-0/share/picard-2.20.4-0

## Fill in path for population specific metadata ##
metadata=$MASTER/SNP_calling/pmau_popgen_metadata.tsv
bam_in=$MASTER/SNP_calling/bams/clean_bams
bam_out=$MASTER/SNP_calling/bams/interim_bams

# Metadata file for each population
# 1 simple_ID
# 2 sample_ID
# 3 read1
# 4 read2
# 5 instrument
# 6 flow_cell
# 7 lane
# 8 barcode
# 9 seq_num
# 10 run_num
# 11 sex
# 12 species
# 13 pop
# 14 colour

simpleID_array=( `cat $metadata | cut -f 1` )
simpleID=${simpleID_array[(($SGE_TASK_ID))]}

instrument_array=( `cat $metadata | cut -f 5` )
instrument=${instrument_array[(($SGE_TASK_ID))]}

seqnum_array=( `cat $metadata | cut -f 9` )
seqnum=${seqnum_array[(($SGE_TASK_ID))]}

flowcell_array=( `cat $metadata | cut -f 6` )
flowcell=${flowcell_array[(($SGE_TASK_ID))]}

lane_array=( `cat $metadata | cut -f 7` )
lane=${lane_array[(($SGE_TASK_ID))]}

barcode_array=( `cat $metadata | cut -f 8` )
barcode=${barcode_array[(($SGE_TASK_ID))]}

## In array
insampleID_array=( `cat $metadata | cut -f 2` )
insampleID=$bam_in/${insampleID_array[(($SGE_TASK_ID))]}

## Out array
outsampleID_array=( `cat $metadata | cut -f 2` )
outsampleID=$bam_out/${outsampleID_array[(($SGE_TASK_ID))]}

## 2. Run picard tools AddreplaceRGs
java -jar $PICARD/picard.jar AddOrReplaceReadGroups \
I=${insampleID}.sorted.dups.bam \
O=${outsampleID}.sorted.freebayes.bam \
RGSM=${simpleID} \
RGLB=${simpleID}.${seqnum} \
RGID=${simpleID} \
RGPU=${flowcell}.${lane}.${barcode} \
RGPL=${instrument} \

## 3. Index the readgroup bam files
java -jar $PICARD/picard.jar BuildBamIndex \
I=${outsampleID}.sorted.freebayes.bam VALIDATION_STRINGENCY=LENIENT

