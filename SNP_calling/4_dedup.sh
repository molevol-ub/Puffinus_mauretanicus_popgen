#!/bin/bash
#$ -t 1-38
#$ -q h13.q
#$ -pe ompi255h13 2

### Script to mark duplicates and remove duplicates

## Load your system modules
# Required modules are: picard tools

## 1. Set your master path
MASTER=~/set/your/path
PICARD=/users/jferrer/programari/anaconda3/pkgs/picard-2.20.4-0/share/picard-2.20.4-0

## 2. Fill in directories if different from the workspace setup
bam_in=$MASTER/SNP_calling/bams/interim_bams
bam_out=$MASTER/SNP_calling/bams/interim_bams

## Fill in path for population specific metadata
metadata=$MASTER/SNP_calling/pmau_popgen_metadata.tsv

## 3. Define "IN" array ##
insampleID_array=( `cat $metadata | cut -f 2` )
insampleID=$bam_in/${insampleID_array[(($SGE_TASK_ID))]}

## Define "OUT" array
outsampleID_array=( `cat $metadata | cut -f 2` )
outsampleID=$bam_out/${outsampleID_array[(($SGE_TASK_ID))]}

## 4. Run picardtools MarkDuplicates
java -jar $PICARD/picard.jar MarkDuplicates \
I=$insampleID.sorted.rg.bam \
O=$outsampleID.sorted.dups.bam \
METRICS_FILE=$outsampleID.metrics.txt \
REMOVE_DUPLICATES=true \
VALIDATION_STRINGENCY=LENIENT AS=true

## 5. Index the deduped bam files
java -jar $PICARD/picard.jar BuildBamIndex \
I=$outsampleID.sorted.dups.bam VALIDATION_STRINGENCY=LENIENT
