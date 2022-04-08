#!/bin/bash
#$ -R y
#$ -N freebayes_parallel_mtDNA
#$ -cwd
#$ -q h13.q
#$ -pe ompi255h13 12

## This script submits the freebayes script (but not in parallel!) for mitochondrial DNA (don't use with nuclear DNA)
source $HOME/.profile

## Set paths

export TMPDIR=/users-d3/jferrer/scratch

MASTER=/users-d3/jferrer/pmau_popgen/SNP_calling
freebayes=/users-d3/jferrer/programari/anaconda3/bin/freebayes
reference=/$MASTER/reference/Pmau_mitogenome_200_15600.fasta
bams=$MASTER/mt_bams/freebayes_bams

# The two following files must have been prepared beforehand with the regions of mtDNA to include (whole mitogenome - generated with fasta_genome_regions.py) and the bam file list

regions=$MASTER/mt_freebayes/pmau.mtDNA.regions
bam_list=$MASTER/mt_freebayes/mtDNA_bams.fofn

out=$MASTER/mt_freebayes/tmp_mtDNA

# Option -p 1 is used when treating with haploid data

freebayes -r Contig01+0404545312:0..15399 -f $reference -p 1 -L $bam_list -v $out/Puffinus_mitogenome.vcf

# Check indexed genome if you get a "unable to find certain scaffold" error
