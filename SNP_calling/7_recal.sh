#!/bin/bash
#SBATCH -D . 
#SBATCH -p pq  
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH -A Research_Project-T110748
#SBATCH --job-name=master_base_recal 
#SBATCH --error=master_base_recal.err.txt 
#SBATCH --output=master_base_recal.out.txt 
#SBATCH --export=All
#SBATCH --array=1-25%5

### Script to perform recalibration on the bam files

# Required modules are: GATK4
module load GATK/4.0.5.1-foss-2018a-Python-3.6.4

## 1. Set your master path
MASTER=/set/your/path

## Fill in directories if different from the workspace setup
## Also add path to indexed reference and recalibration file
reference=/gpfs/ts0/home/jrp228/startup/STAR/STAR.chromosomes.release.fasta
recal=/gpfs/ts0/home/jrp228/startup/STAR/STAR_intersect_recal.vcf
bam_in=$MASTER/SNP_calling/bams/interim_bams
bam_out=$MASTER/SNP_calling/bams/clean_bams

## This just catches the array in case it's running for a value with no individual (this screws with the outputs)
IND_N=$(cut -f1 $metadata | tail -n+2 | uniq | awk 'NF > 0' | wc -l)

if [ $SLURM_ARRAY_TASK_ID -le $IND_N ]
then

## 2. Define "IN" array ##
insampleID_array=( `cat $samples | cut -f 1` )
insampleID=$bam_path/${insampleID_array[(($SLURM_ARRAY_TASK_ID))]}

## Define "OUT" array
outsampleID_array=( `cat $samples | cut -f 1` )
outsampleID=$bam_path/${outsampleID_array[(($SLURM_ARRAY_TASK_ID))]}

## 3. Apply base quality score recalibration

## BaseRecal
gatk --java-options "-Xmx8g" BaseRecalibrator \
-I $insampleID.merged.sorted.dups.bam \
-R $reference --known-sites $recal \
-O $outsampleID.table

## 4. ApplyBQSR
gatk --java-options "-Xmx8g" ApplyBQSR \
$insampleID.merged.sorted.dups.bam \
-R $reference --bqsr-recal-file $outsampleID.table \
-O $outsampleID.recal.bam

## 5. Index the recalibrated bam files
java -Xmx10g -jar $EBROOTPICARD/picard.jar BuildBamIndex \
I=$outsampleID.recal.bam VALIDATION_STRINGENCY=LENIENT

fi
