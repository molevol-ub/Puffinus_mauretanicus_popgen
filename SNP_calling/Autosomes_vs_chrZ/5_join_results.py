#!/usr/bin/env python

# Script to join result summaries of all chromosomes into one and compare with chrZ

from datetime import datetime
import re
import os
print (datetime.now())

# 1 Prepare the chrz file

Z_file=open("/set/your/path/sphen_genome/chrz_final_results.1000.txt")
new_file=open("/set/your/path/sphen_genome/provisional.txt", "w")

for line in Z_file:

	expression=re.search(r"(scf.+)\t(.+)\t(.+)\t(.+)\t(.+)\n", line)
	scf_id=expression.group(1)
	hit_ratio=expression.group(4)
	coverage=expression.group(5)
	
	new_file.write(scf_id + "\t" + coverage + "\t" + hit_ratio + "bar\n")

Z_file.close()
new_file.close()

# 2 Pass all info of hit_ratio to 1 file and compare along the way; if the conclusion is that a scaffold does NOT belong to chrZ we will introduce "signal" to the line in file

def merge(ch):
	
	Z_file=open("/set/your/path/sphen_genome/provisional.txt")
	new_file=open("/set/your/path/sphen_genome/provisional2.txt", "w")
	
	for line in Z_file:
		
		chr_file=open("/set/your/path/sphen_genome/" + ch + "_final_results.1000.txt")
		
		expression=re.search(r"(scf[0123456789]{13})\t(.+)\t(.+)bar", line)
		scf_id=expression.group(1)
		hit_ratio=float(expression.group(3))
		
		break_parameter=0
		
		for line2 in chr_file:
			
			expression2=re.search(r"(scf.+)\t(.+)\t(.+)\t(.+)\t(.+)\n", line2)
			scf_id2=expression2.group(1)
			hit_ratio2=float(expression2.group(4))
			
			if break_parameter!= 0:
				
				chr_file.close()
				break
				
			
			elif scf_id2 == scf_id:
				
				if hit_ratio > hit_ratio2:
				
					new_line=line.rstrip()
					new_file.write(new_line + "\t" + str(hit_ratio2) + "\n")
					break_parameter+=1
					
				else:
					new_line=line.rstrip()
					new_file.write(new_line + "\t" + str(hit_ratio2) + "signal\n")
					break_parameter+=1
	
	chr_file.close()
	Z_file.close()
	new_file.close()
	
	os.rename("/set/your/path/sphen_genome/provisional2.txt", "/set/your/path/sphen_genome/provisional.txt")

chr_list=("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chr23", "chr24", "chr25")

for chrom in chr_list:
	
	merge(chrom)

# 3 Finnally write two output files: one with the bad scaffolds ("signal") and one without them. In both output files we will only include the coverage and hit-ratio for chrZ.

info_file=open("/set/your/path/sphen_genome/provisional.txt")
good_file=open("/set/your/path/sphen_genome/chr_Z_summary.txt", "w")
bad_file=open("/set/your/path/sphen_genome/chr_auto_summary.txt", "w")

for line in info_file:

	expression=re.search(r"(scf[0123456789]{13})\t(.+)\t(.+)bar", line)
	chr_id=expression.group(1)
	coverage=expression.group(2)
	hit_ratio=expression.group(3)
	
	bad_parameter="signal"
	
	if bad_parameter in line:
	
		bad_file.write(chr_id + "\t" + coverage + "\t" + hit_ratio + "\n")
	
	else:
		
		good_file.write(chr_id + "\t" + coverage + "\t" + hit_ratio + "\n")
		
info_file.close()
good_file.close()
bad_file.close()

os.remove("/set/your/path/sphen_genome/provisional.txt")

print("Procedure completed: " + str(datetime.now()))
