#!/bin/bash
#$ -cwd
#$ -R y
#$ -e plink_map.err
#$ -o plink_map.out
#$ -q h12.q
#$ -pe ompi128h12 2
#$ -V                    #export environment var
#$ -N plink_mapfile             #name Job

# Script to create mapfile - later we can run selscan

WORKDIR=/users-d3/jferrer/gizquierdo/TFG/genome_scans/selscan/input
infile=/users-d3/jferrer/gizquierdo/TFG/genome_scans/vcfs/Z/Puffinus_SNP.maxmiss100.filtered.noPP.merged.nomono.masked.vcf.gz

# Before anything, we shall also use vcftools to generate vcfs with each population (reference --> yelkouan)

# Check inds always!

vcftools --gzvcf $infile --indv M21 --indv M4 --indv M2 --indv M18 --indv M12 --indv M13 --indv M17 --indv M5 --remove-filtered-all --recode --out $WORKDIR/Puffinus_selscan_Z.mauretanicus
vcftools --gzvcf $infile --remove-indv M21 --remove-indv M4 --remove-indv M2 --remove-indv M18 --remove-indv M12 --remove-indv M13 --remove-indv M17 --remove-indv M5 --remove-filtered-all --recode --out $WORKDIR/Puffinus_selscan_Z.yelkouan

bgzip $WORKDIR/Puffinus_selscan_Z.mauretanicus.recode.vcf
tabix $WORKDIR/Puffinus_selscan_Z.mauretanicus.recode.vcf.gz

bgzip $WORKDIR/Puffinus_selscan_Z.yelkouan.recode.vcf
tabix $WORKDIR/Puffinus_selscan_Z.yelkouan.recode.vcf.gz

#-----------------------------------------------------------------------------------

# Then first make list with all scaffolds:

zcat $infile | grep -v "^#" | cut -f 1 | sort | uniq > $WORKDIR/scaff_names

# And make vcf file with each one of them: also make mapfile using plink and format it correctly

while read scaff;
	do
		echo $scaff
		bcftools view $WORKDIR/Puffinus_selscan_Z.mauretanicus.recode.vcf.gz $scaff > $WORKDIR/$scaff.mauretanicus.vcf;
		bcftools view $WORKDIR/Puffinus_selscan_Z.yelkouan.recode.vcf.gz $scaff > $WORKDIR/$scaff.yelkouan.vcf;
		plink --vcf $WORKDIR/$scaff.mauretanicus.vcf --double-id --allow-extra-chr --recode --out $WORKDIR/$scaff.selscan;
		
		cat $WORKDIR/$scaff.selscan.map | cut -f 1 > $WORKDIR/col1;
		cat $WORKDIR/$scaff.selscan.map | cut -f 4 > $WORKDIR/col4;
		
		paste -d "\t" $WORKDIR/col1 $WORKDIR/col1 $WORKDIR/col4 $WORKDIR/col4 > $WORKDIR/$scaff.selscan.map;
		
		rm $WORKDIR/col*;
		
	done < $WORKDIR/scaff_names


rm $WORKDIR/*.ped
rm $WORKDIR/*.nosex
rm $WORKDIR/*.log

#----------------------------------------------------------------------------------
