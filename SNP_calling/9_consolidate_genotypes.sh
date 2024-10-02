#!/bin/bash
#$ -t 1-4171
#$ -tc 25
#$ -q h13.q
#$ -pe ompi255h13 2

### Script to combine gVCFs and then genotype across the cohort
source $HOME/.profile

## 1. Set your master path
MASTER=/users/jferrer/pmau_popgen

## Fill in directories if different from the workspace setup
## Also add path to reference
gvcfs=$MASTER/SNP_calling/gvcfs
output_vcfs=$MASTER/SNP_calling/vcfs/intervals
reference=$MASTER/SNP_calling/reference/genome.fa

## Name your dataset
DATASET=Puffinus

# 2. Make the command for variants
infiles=(`ls -1 ${gvcfs}/*.gvcf.gz`)
let len=${#infiles[@]}-1
cmd=""
for i in `seq 0 $len`
do
    cmd+="--variant ${infiles[$i]} "
done

# 3. Fetch the interval of interest
interval=$(awk '{print $1}' ${reference}.fai | sed "${SGE_TASK_ID}q;d")
echo $interval
# Remove if done already
#rm -rf $gvcfs/INTERVAL_${interval}_db

# 4. Make database
gatk --java-options "-Xmx16g -Xms4g" GenomicsDBImport \
    $cmd \
    --genomicsdb-workspace-path $output_vcfs/INTERVAL_${SGE_TASK_ID}_db \
    --reader-threads 2 \
    --intervals $interval

# 5. Genotype
#cd $output_vcfs
cd $gvcfs

gatk --java-options "-Xmx16g -Xms4g" GenotypeGVCFs \
    -R $reference \
    -V gendb://INTERVAL_${SGE_TASK_ID}_db \
    -O ${interval}.vcf.gz

# Tidy 
rm -rf *INTERVAL_${SGE_TASK_ID}_db*

### 6. Run this when the above is all good
ls $gvcfs/*vcf.gz >> batch_filter.txt
grep -Fxf batch_filter.txt batch_inputs.txt > batch_inputs_2.txt
bcftools concat -o $gvcfs/${DATASET}_cohort_batch_genotyped.g.vcf -f $gvcfs/intervals/batch_inputs_2.txt

