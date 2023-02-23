#!/bin/bash
#$ -cwd
#$ -R y
#$ -e ballermix.err
#$ -o ballermix.out
#$ -q h13.q
#$ -pe ompi255h13 2
#$ -V                    #export environment var
#$ -N ballermix.def            #name Job
#$ -M 000izquierdoguillem@gmail.com
#$ -m be

# This first script is intended to generate the input files and run BalLerMix+ to calculate the B_0_MAF statistic for your scaffold:

workdir=/users-d3/jferrer/gizquierdo/TFM/genome_scans/BalLerMix+/Mallorca/scfs
cd $workdir

# 1. Generating input files: you will aid yourself with the already prepared files for BetaScan and SFS

#	a) Generate the input file from the BetaScan input - this command eliminates any sites with MAF=0 (caused by missing data) and adds a genetic distance column (which we will not use):

#while read scaff;
#        do


#		beta_file=/users-d3/jferrer/gizquierdo/TFM/genome_scans/BetaScan/Mallorca/scfs/$scaff.no_out.beta.txt.gz
#		input_file=$workdir/$scaff.beta.txt

#		zcat $beta_file | cut -f 1 > $input_file
#		zcat $beta_file > temp

#		paste $input_file temp > temp2 | mv temp2 $input_file				# Up to here to generate the 4 columns

#		awk '$3 == 0 {next} {print}' $input_file > temp

#		mv temp $input_file

		# Also, all allelle counts have to be based on the total number of indvs, irrespective of the amount of missing data. Therefore, if your number of indvs is, for ex, 34, you have to scale it to that:

#		awk '{temp=($3*24)/$4; printf"%.0f\t24\n", temp}' $input_file > temp
#		cat $input_file | cut -f 1,2 > temp2
#		paste temp2 temp > $input_file
#		rm temp*


#	b) Generate helper file (whole-genome SFS) manually from the file you obtained using easySFS.py
#
# It's structure should be:	1	34	0.4345
#				2	34	0.0679
#				...	...	...

# Up until the halfway point (maximum MAC) - the 2nd column is the total haploid indv count - the 3rd column should sum 1!

#		help_file=$workdir/auto_SFS.txt

#--------------------------------------------------------------------------------------------------

# 2. Run BalLerMix+
#
# Options:	--noSub --MAF --> use both options to indicate we want to calculate B_0_MAF
#		-i --> inputfile
#		--spect --> golbal sfs file
#		--usePhysPos --> use bp instead of genetic_distance (cM) - assume uniform recombnation rate of 1e-6!

#		--fixX --> fix the equilibrium frequency to search at; as we don't want to consider sites that could be considered neutral/positive selection, make iterations only for x= {0.2,0.3,0.4,0.5} - also reduces computaionla intensiveness
#		--listA <min,max,step> --> range for the parameter A; as we only wan to detect recent sweeps, low values of A (long regions of linkage) should be ideal; set the step to 200 to reduce computational intensity

# window size options

#		ballermix=/users-d3/jferrer/gizquierdo/TFM/genome_scans/scripts/BalLerMix+/BalLeRMix+_v1.py

#		python $ballermix -i $input_file --noSub --MAF --usePhysPos --fixX 0.5 --listA 200,400,600,800,1000,1200,1400,1600,1800,2000 --fixWinSize -w 12500 --fixAlpha 2 --spect $help_file -o $workdir/$scaff.x05.a2.w25k.ballermix.out
#		python $ballermix -i $input_file --noSub --MAF --usePhysPos --fixX 0.5 --listA 200,400,600,800,1000,1200,1400,1600,1800,2000 --fixWinSize -w 12500 --fixAlpha 5 --spect $help_file -o $workdir/$scaff.x05.a5.w25k.ballermix.out
#		python $ballermix -i $input_file --noSub --MAF --usePhysPos --fixX 0.5 --listA 200,400,600,800,1000,1200,1400,1600,1800,2000 --fixWinSize -w 25000 --fixAlpha 10 --spect $help_file -o $workdir/$scaff.x05.a10.w50k.ballermix.out
#		python $ballermix -i $input_file --noSub --MAF --usePhysPos --fixX 0.5 --listA 200,400,600,800,1000,1200,1400,1600,1800,2000 --fixWinSize -w 12500 --fixAlpha 25 --spect $help_file -o $workdir/$scaff.x05.a25.w25k.ballermix.out

#	done < $workdir/scaff_names

#---------------------------------------------------------------------------------------------

# 3. Join all files under one, with the name of the scaffold in the first column

#while read scaff;
#	do

#	while read line;
#		do
                        
#			echo "$scaff	$line" >> $workdir/../Pmau.x05.a10.w50k.ballermix.txt

#		done < $workdir/$scaff.x05.a10.w50k.ballermix.out

#	done < $workdir/scaff_names

#---------------------------------------------------------------------------------------------

# 4. Calculate the average and write it down:

cd $workdir/..

scf_file=/users-d3/jferrer/gizquierdo/TFM/genome_scans/BalLerMix+/Mallorca/Pmau.x05.a10.w50k.ballermix.def.csv
beta_file=/users-d3/jferrer/gizquierdo/TFM/genome_scans/BalLerMix+/Mallorca/Pmau.x05.a10.w50k.ballermix.txt

# a) Extract the position and scaffold columns:

cat ~/gizquierdo/TFM/genome_scans/selscan/auto_def_wM8/output_unphased/auto_selscan.25.csv | cut -f 9,10,11 > $scf_file

# b) Run python script for average

script_dir=/users-d3/jferrer/gizquierdo/TFM/genome_scans/scripts/BalLerMix+

python $script_dir/average_BalLerMix.py $beta_file $scf_file

