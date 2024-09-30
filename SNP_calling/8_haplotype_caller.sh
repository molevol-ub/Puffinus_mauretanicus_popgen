#!/bin/bash
#$ -t 9-38
#$ -tc 12
#$ -q h11.q
#$ -pe ompi64h11 4

### Script to run haplotype caller

## 1. Set your master path
MASTER=/set/your/wd

## Fill in directories if different from the workspace setup
## Also add path to reference
bam_in=$MASTER/SNP_calling/bams/clean_bams
gvcfs=$MASTER/SNP_calling/gvcfs
reference=$MASTER/SNP_calling/reference/genome.fa

## Fill in path for population specific metadata
metadata=$MASTER/SNP_calling/pmau_popgen_metadata.tsv

## 2. Define "IN" array ##
insampleID_array=( `cat $metadata | cut -f 1` )
insampleID=$bam_in/${insampleID_array[(($SGE_TASK_ID))]}

## Define "OUT" array
outsampleID_array=( `cat $metadata | cut -f 1` )
outsampleID=$gvcfs/${outsampleID_array[(($SGE_TASK_ID))]}

## 3. Run haplotype caller
/users/jferrer/programari/gatk-4.1.9.0/gatk --java-options "-Xmx10g" HaplotypeCaller \
-R $reference \
-I ${insampleID}.sorted.dups.bam \
-O ${outsampleID}.gvcf.gz -ERC GVCF
