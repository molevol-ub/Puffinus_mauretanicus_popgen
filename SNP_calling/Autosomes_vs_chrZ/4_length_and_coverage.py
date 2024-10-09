#!/usr/bin/env python

# Script to filter results according to length and add coverage information to the filtered_hit_data

from datetime import datetime
import re
import os
print (datetime.now())

def whole_process(ch):
	
	new_file=open("/set/your/path/sphen_genome/" + ch + "_final_results.1000.txt", "w")
	
	#1 Take those scaffolds that you will use
	
	scaffold_list=[]
	length_list=[]
	
	scaffold_file=open("/set/your/path/pmau_scaffolds_75k.txt")
	
	for line in scaffold_file:
		
		expression=re.search(r"(scf.+)\t(.+)\n", line)
		
		scf=expression.group(1)
		length=expression.group(2)
		
		if scf not in scaffold_list:
			scaffold_list.append(scf)
			length_list.append(length)
	
	if len(scaffold_list)>1:
		print ("Scaffold_list_complete\n" + str(datetime.now()))
		
	# 2 Filter out those scaffolds not included in the scaffold list and add scaffold_length to file
	
	hit_file=open("/set/your/path/sphen_genome/" + ch + "_hit_length.1000.txt")
	
	for line in hit_file:
		
		expression=re.search(r"(scf.+)\t(.+)\n", line)
		scf=expression.group(1)
		
		if scf in scaffold_list:
			index=scaffold_list.index(scf)
			new_line = line.rstrip()
			new_file.write(new_line + "\t" + length_list[index] + "\n")
			
	new_file.close()
	hit_file.close()
	
	# 3 Calculation of % of scaffold length
	
	info_file=open("/set/your/path/sphen_genome/" + ch + "_final_results.1000.txt")
	new_file=open("/set/your/path/sphen_genome/" + ch + "_final_hit_ratio.1000.txt", "w")
	
	for line in info_file:
	
		expression=re.search(r"(scf.+)\t(.+)\t(.+)\n", line)
		
		hit_length=float(expression.group(2))
		scf_length=float(expression.group(3))
		
		hit_ratio = (100*hit_length)/scf_length
		new_line=line.rstrip()
		
		new_file.write(new_line + "\t" + str(hit_ratio) + "\n")
	
	info_file.close()
	new_file.close()
	
	os.remove("/set/your/path/sphen_genome/" + ch + "_final_results.1000.txt")
	
	# 4 Add coverage information
	
	coverage_file=open("/set/your/path/scf_info/coverage_ratios/scaffolds_def_75k.txt")
	info_file=open("/set/your/path/sphen_genome/" + ch + "_final_hit_ratio.1000.txt")
	new_file=open("/set/your/path/sphen_genome/" + ch + "_final_results.1000.txt", "w")
	
	scaffold_list_2=[]
	coverage_list=[]
	
	for line in coverage_file:
		
		expression=re.search(r"(scf.+)\t(.+)\n", line)
		
		scf=expression.group(1)
		coverage=expression.group(2)
		
		if scf in scaffold_list:
			scaffold_list_2.append(scf)
			coverage_list.append(coverage)
			
	for line in info_file:
	
		expression=re.search(r"(scf.+)\t(.+)\t(.+)\t(.+)\n", line)
		scf=expression.group(1)
		
		if scf in scaffold_list:
			index=scaffold_list.index(scf)
			new_line = line.rstrip()
			new_file.write(new_line + "\t" + coverage_list[index] + "\n")
	
	coverage_file.close()
	info_file.close()
	new_file.close()
	
	os.remove("/set/your/path/sphen_genome/" + ch + "_final_hit_ratio.1000.txt")
	
	print ("Procedure complete\n" + str(datetime.now()))


# 2 Define list of chromosomes and loop through them

chr_list=("chrz", "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chr23", "chr24", "chr25")

for chrom in chr_list:
	
	whole_process(chrom)
