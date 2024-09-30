#!/bin/bash
#$ -R y
#$ -N freebayes_parallel
#$ -e logs/freebayes_parallel.err.txt
#$ -o logs/freebayes_parallel.out.txt
#$ -cwd
#$ -q h12.q
#$ -pe ompi128h12 64


## This script submits the freebayes parallel script (which has been edited from the original, by removal of tke -k option)
source $HOME/.profile

## 1. Set paths

export TMPDIR=/users/jferrer/scratch

MASTER=/users/jferrer/pmau_popgen/SNP_calling
freebayes=/users/jferrer/programari/anaconda3/bin/freebayes
ref=$MASTER/reference/genome.fa
regions=$MASTER/reference/targets.region
bams=$MASTER/bams/clean_bams
bam1=$bams/PORQ.sorted.dups.bam
bam_list=$MASTER/freebayes/all_bams2.fofn
out=$MASTER/freebayes/tmp

## 2. Run freebayes-parallel

$MASTER/scripts/freebayes-parallel_no_sort $regions 64 -f $ref -L $bam_list -v $out/{}.vcf -g 500 --use-best-n-alleles 4
