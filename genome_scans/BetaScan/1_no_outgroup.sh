#!/bin/bash
#$ -cwd
#$ -R y
#$ -e betascan.def.err
#$ -o betascan.def.out
#$ -q h10.q
#$ -pe ompi64h10 2
#$ -V                    #export environment var
#$ -N betascan.def            #name Job
#$ -M 000izquierdoguillem@gmail.com
#$ -m be

workdir=/users-d3/jferrer/gizquierdo/TFM/genome_scans/BetaScan/Mallorca/scfs
cd $workdir
# 1. Convert vcf to gacf using glactools - in this case do not include outgroup

# Use datasets with missing data

# Options:	--onlyGT --> ignores GT and GLs
#		--fold --> when you are not using anc/derived data

vcf_file=$workdir/Mallorca.maxmiss80.masked.nomono.vcf.gz
genome_file=/users-d3/ccuevas/Functional_Genome_Annotation/scripts/braker3/3pmaureta/genome.fa

# First, for the fastas, divide the main one using:

#csplit $genome_file /\>scf.*/ {*};					# Use csplit to extract each scaffold
#for a in x*; do echo $a; mv $a $(head -1 $a| cut -c 2-).fa; done;	# Rename each fasta file

# Then, for the vcfs, first make list with all scaffolds:

#zcat $vcf_file | grep -v "^#" | cut -f 1 | sort | uniq > $workdir/scaff_names

# And make vcf file with each one of them; also index the fastas

while read scaff;
	do
		echo $scaff;
		bcftools view $vcf_file $scaff > $workdir/$scaff.vcf;
		
		samtools faidx $scaff.fa;				# Index the fasta files
		
	done < $workdir/scaff_names

# Now run glactools to obtain the betascan input file

while read scaff;
	do
		vcf_scf=$workdir/$scaff.vcf;
		fasta_scf=$workdir/$scaff.fa.fai;
		
		~/programari/anaconda3/bin/glactools/glactools vcfm2acf --onlyGT --fai $fasta_scf $vcf_scf > $workdir/$scaff.acf.gz;
		~/programari/anaconda3/bin/glactools/glactools acf2betascan --fold $workdir/$scaff.acf.gz | gzip > $workdir/$scaff.no_out.beta.txt.gz;
		
	done < $workdir/scaff_names

# ----------------------------------------------------------------------------------------------------------------------

# 2. Run BetaScan for each scaffold.
# Options:	-w --> window size (start by 25000, to compare with PopGenome?)
#		-m --> maf filter (choose 0.12 (mac=4) to avoid low-freq variants and sequencing errors)
#		-fold --> use the folded version of beta to take into account the absence of outgroup

#while read scaff;
#	do
		
#		python /users-d3/jferrer/gizquierdo/TFM/genome_scans/scripts/BetaScan/BetaScan.py -i $workdir/$scaff.no_out.beta.txt.gz -w 25000 -m 0.25 -fold -o $workdir/$scaff.25k.m25.txt
#		python /users-d3/jferrer/gizquierdo/TFM/genome_scans/scripts/BetaScan/BetaScan.py -i $workdir/$scaff.no_out.beta.txt.gz -w 25000 -m 0.25 -p 10 -fold -o $workdir/$scaff.25k.m25.p10.txt
#		python /users-d3/jferrer/gizquierdo/TFM/genome_scans/scripts/BetaScan/BetaScan.py -i $workdir/$scaff.no_out.beta.txt.gz -w 25000 -m 0.25 -p 20 -fold -o $workdir/$scaff.25k.m25.p20.txt

#	done < $workdir/scaff_names

#---------------------------------------------------------------------------------------------

# 3. Join all files under one, with the name of the scaffold in the first column

#while read scaff;
#	do
#
#	while read line;
#		do
			
#			echo "$scaff	$line" >> $workdir/../BetaScan.no_out.25k.m25.txt

#		done < $workdir/$scaff.25k.m25.txt
#
#	while read line;
#		do
#
#			echo "$scaff	$line" >> $workdir/../BetaScan.no_out.25k.m25.p10.txt
#
#		done < $workdir/$scaff.25k.m25.p10.txt
#
#	while read line;
#		do
#
#			echo "$scaff	$line" >> $workdir/../BetaScan.no_out.25k.m25.p20.txt
#
#		done < $workdir/$scaff.25k.m25.p20.txt
#
#	while read line;
#		do
#
#			echo "$scaff	$line" >> $workdir/../BetaScan.no_out.25k.m06.txt
#
#		done < $workdir/$scaff.25k.m06txt
#
#	while read line;
#		do
#
#			echo "$scaff	$line" >> $workdir/../BetaScan.no_out.25k.m25.txt
#
#		done < $workdir/$scaff.25k.m25.txt
#
#	done < $workdir/scaff_names

#---------------------------------------------------------------------------------------------------------------------

# 4. Prepare file for plotting:

# a) First create alternative *.csv files for each summary *.txt file created before; they should only contain a header: "scf	Postition	Beta1"

# b) Then append the contents of your *.txt files int these *.csv while changing all their "headers" into NAs

# E.g.: sed 's/Position\tBeta1\*/NA\tNA/g' BetaScan.no_out.25k.m0.txt >> BetaScan.no_out.25k.m0.csv

# c ) Finally, remove all your old *.txt files if the transfer has been correct

#----------------------------------------------------------------------------------------------------------------------

# 5. Calculate the average and write it down:

#cd $workdir/..

#scf_file=/users-d3/jferrer/gizquierdo/TFM/genome_scans/BetaScan/no_outgroup/BetaScan.no_out.25k.m25.p20.def.csv
#beta_file=/users-d3/jferrer/gizquierdo/TFM/genome_scans/BetaScan/no_outgroup/BetaScan.no_out.25k.m25.p20.csv

# a) Extract the position and scaffold columns:

#cat ~/gizquierdo/TFM/genome_scans/selscan/auto_def_wM8/output_unphased/auto_selscan.25.csv | cut -f 9,10,11 > $scf_file

# b) Run python script for average

#script_dir=/users-d3/jferrer/gizquierdo/TFM/genome_scans/scripts/BetaScan

#python $script_dir/average_beta.py $beta_file $scf_file

