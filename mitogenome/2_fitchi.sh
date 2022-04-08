#!/bin/bash

# Short script to run fitchi (it generate the haplotype genealogy graph)
# Check out the requirements to install fitchi (and its manual) in: https://evoinformatics.group/fitchi.html

# The output of raxml-ng (besttree(?)) has to be transformed manually into a nexus file. Also, the individual names must include their population (eg. M16_Menorca)

# Options:	-e "INT" --> minimum edge length - this will cluster together individuals separated by les than "INT" substitutions
#			-p "list" --> a list of all the populations you want to colour the graph by

python3 /path/to/fitchi.py mtDNA.nexus mtDNA_pmau.html -e 5 --haploid -p puffinus yelkouan Menorca Malllorca_Malgrats Mallorca_SaCella Cabrera Pitiuses
