#!/bin/bash

# Script to run treemix

MASTER=/users-d3/jferrer/gizquierdo/TFG/treemix
input_file=$MASTER/input_files/pmau_treemix.LDpruned.treemix.frq.gz
output=pmau_treemix

# Options: -m "INT" --> specifies the number of admixture events
#          -o "STRING" --> output
#          -bootstrap -k "INT" --> the number of bootstraps
#          -noss --> to avoid making correction for few inds (not applied)

for i in {0..10}
	do
		treemix -i $input_file -m $i -o $MASTER/bootstrap_1000_output/$output.$i -root P.puffinus -bootstrap -k 1000 > $MASTER/bootstrap_1000_output/$output.$i.log
	done

# Now we repeat the process with -noss (which shouldn't change much the output)

for i in {0..10}
        do
                treemix -i $input_file -m $i -o $MASTER/uncorrected_1000_output/$output.$i.uncorrected -root P.puffinus -bootstrap -k 1000 > $MASTER/uncorrected_1000_output/$output.$i.uncorrected.log
        done

