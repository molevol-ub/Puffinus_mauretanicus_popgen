#!/bin/bash
#$ -t 1-2
#$ -q h13.q
#$ -pe ompi255h13 4

source $HOME/.profile

### Script to run a bwa mem, bam and sort bam
# Required modules are: bwa, samtools, bcftools

## 1. Set your master path
MASTER=~/set/your/path

## Fill in directories if different from the workspace setup
## Also add path to indexed reference
reference=/users/ccuevas/Functional_Genome_Annotation/scripts/braker3/3pmaureta/genome.fa
input_reads=$MASTER/SNP_calling/reads/clean_reads
bam_dir=$MASTER/SNP_calling/bams/raw_bams

echo "input directory" $input_reads
echo "bam directory" $bam_dir

## 2. Fill in path for population specific metadata
metadata=$MASTER/SNP_calling/pmau_popgen_failed.tsv
echo "metadata" $metadata

## 3. Read 1 array ##
read1_array=(`cat $metadata | cut -f 3`)
read1=$input_reads/${read1_array[(($SGE_TASK_ID))]}

## Read 2 array ##
read2_array=( `cat $metadata | cut -f 4` )
read2=$input_reads/${read2_array[(($SGE_TASK_ID))]}

## 4. Output array ##
out_array=( `cat $metadata | cut -f 2` )
bam_out=$bam_dir/${out_array[(($SGE_TASK_ID))]}

echo "reference" $reference
echo "read1" $read1
echo "read2" $read2
echo "alignment" ${bam_out}.unsorted.raw.sam

## 5. Align with bwa mem using 24 cores. Again, make sure the read name prefixes match!
bwa mem -t 4 $reference ${read1} ${read2} > ${bam_out}.unsorted.raw.sam

## Convert bam to sam, sort bam, index, flagstat
samtools view -bS ${bam_out}.unsorted.raw.sam > ${bam_out}.unsorted.raw.bam
samtools sort ${bam_out}.unsorted.raw.bam -o ${bam_out}.sorted.raw.bam
samtools index ${bam_out}.sorted.raw.bam
samtools flagstat ${bam_out}.sorted.raw.bam > ${bam_out}.mappingstats.txt

# ## Remove the sam and unsorted bam files
#rm ${bam_dir}/*.sam
#rm ${bam_dir}/*.unsorted.raw.bam

## 6. To check that the bams are not corrupted, run (in the directory where the bams are):

samtools quickcheck -v *.sorted.raw.bam > bad_bams.fofn   && echo 'all ok' || echo 'some files failed check, see bad_bams.fofn'
