#!/bin/bash
#SBATCH -D .
#SBATCH -p <queue_name>
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -A <research_project>
#SBATCH --job-name=qc_clean
#SBATCH --error=./logs/qc_clean.err.txt
#SBATCH --output=./logs/qc_clean.out.txt
#SBATCH --export=All
#SBATCH --array=1-20

## Load your system modules
# Required modules are: FastQC, Python2, Trim_Galore

## 1. Set your master path
MASTER=<path>

## Fill in path for population specific metadata ##
metadata=$MASTER/SNP_calling/metadata.tsv

## Change these directories if different from the workspace setup
raw_reads=$MASTER/SNP_calling/reads/raw_reads
clean_reads=$MASTER/SNP_calling/reads/clean_reads
fastqc_raw=$MASTER/SNP_calling/reads/raw_reads/fastqc
fastqc_clean=$MASTER/SNP_calling/reads/clean_reads/fastqc
                                                                                                                                                   
## 2. Create an array to hold all of the files within the raw reads location
read1=( `cat $metadata | cut -f 2` )
read1_array=$raw_reads/${read1[(($SLURM_ARRAY_TASK_ID))]}

read2=( `cat $metadata | cut -f 2` )
read2_array=$raw_reads/${read2[(($SLURM_ARRAY_TASK_ID))]}


## 3. Testing that all variables are working correctly                                                                                                                                                                  
echo "read1" $read1
echo "read2" $read2
echo "output directory" $clean_reads

## This just catches the array in case it's running for a value with no individual (this screws with the outputs)
IND_N=$(cut -f1 $metadata | tail -n+2 | uniq | awk 'NF > 0' | wc -l)

if [ $SLURM_ARRAY_TASK_ID -le $IND_N ]
then

## 4. Run fastqc on raw reads

### Change the suffix of the file names to reflect your read 1 and read 2 names :)
fastqc ${read1_array}_R1.fastq.gz ${read2_array}_R2_001.fastq.gz -o $fastqc_raw

## Run trim_galore default settings (but adjust adaptors if needed) - change your phred score if old Illumina style with the option --phred64
trim_galore -q 20 --path_to_cutadapt cutadapt -o $clean_reads --paired ${read1_array}_R1.fastq.gz ${read2_array}_R2.fastq.gz

## 5. Run fastqc on clean reads
fastqc ${read1_array}_R1_001_val_1.fq.gz ${read2_array}_R1_001_val_2.fq.gz -o $fastqc_raw

### DONE! ###
echo "1_qc_clean.sh is complete"
