#!/bin/bash
#$ -cwd
#$ -R y
#$ -e TreeMix.err
#$ -o TreeMix.out
#$ -q h12.q
#$ -pe ompi128h12 2
#$ -V                    #export environment var
#$ -N TreeMix             #name Job

# Script to run TreeMix to infer:	a) Phylogenetic tree with 100 bootstrap pseudoreplicates (0 migration edges) using your ind_input
#					b) Population tree with migration edges with 100 bootstrap pseudoreplicates (0-10 migration edges) using your pop_input

MASTER=/your/working/directory
input_file_inds=$MASTER/your_ind_file.treemix.frq.gz
input_file_pops=$MASTER/your_pop_file.treemix.frq.gz
output=output_name

# Options: -m "INT" --> specifies the number of admixture events
#          -o "STRING" --> output
#	   -k --> SNP blocks to account for linkage disequilibrium
#          -bootstrap -k "INT" --> the bootstrap blocks
#          -noss --> to avoid making correction for few inds (not applied)
#	   -root --> the outgroup pop/inds; change to your convenience

#------------------------------------------------------------------------------------------------------------------------------------------------------

# 1 Generate the main phylogenetic tree without bootstraping; k in this case used as to take into account linkage desequilibrium (blocks of 1000 SNPs considered together?)

treemix -i $input_file_inds -m 0 -o $MASTER/$output.inds.main -root COP1,LT2 -k 1000 > $MASTER/$output.inds.main.log

# 2 Now generate a 100 bootstrap replicates using k=1000 SNP blocks

for i in {0..100}
	do
		treemix -i $input_file_inds -m 0 -o $MASTER/$output.inds.$i.replicate -root COP1,LT2 -bootstrap -k 1000 > $MASTER/$output.inds.$i.replicate.log
	done

#-------------------------------------------------------------------------------------------------------------------------------------------------------

# Now, infer population trees; first infer the 100 bootstrap replicates

# 3 To use OptM afterwards, run 5 iterations (j=[1,10]) for each tree

for i in {0..10}
	do
	for j in {1..100}
		do
		treemix -i $input_file_pops -m $i -o $MASTER/$output.pops.$j.$i -root P.puffinus -bootstrap -k 1000 >> $MASTER/$output.pops.$j.$i.log
	done
done

rm $MASTER/*edg* | rm $MASTER/*vert* | rm $MASTER/*covse* | rm $MASTER/*log*

# 4 And now calculate the main tree (as to not remove any of its results with the previous step)

for i in {0..10}
	do
	treemix -i $input_file_pops -m $i -o $MASTER/$output.pops.$i.main -root P.puffinus -k 1000 >> $MASTER/$output.pops.main.$j.$i.log
done

