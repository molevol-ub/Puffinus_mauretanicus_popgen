#!/bin/python

# Script dedicated to making a list of scaffolds, its length and its coverage ratio

from datetime import datetime
import re
import os

print (datetime.now())

# 1 Prepare scaffold list

scaffold_file=open("/set/your/path/pmau_scaffolds.txt")
list_file=open("/set/your/path/coverage_ratios/scaffolds_pre.txt", "w")

scaffold_list = []

for line in scaffold_file:
	
	expression=re.search(r"(scf.+)\t(.+)\n", line)
	
	scf=expression.group(1)
	
	if scf not in scaffold_list:
		scaffold_list.append(scf)
		list_file.write(scf + "\n")

if len(scaffold_list)>1000:
	print ("Scaffold_list_complete\n" + str(datetime.now()))
	
scaffold_file.close()
list_file.close()

# 2 Define function: identify scaffold and coverage and write it into output

def data_extraction(ind):

	pmau_file=open("/users-d3/jferrer/pmau_popgen/SNP_calling/bams/qualimap/" + ind + ".sorted.dups.bam_qualimap/genome_results.txt")
	
	list_file=open("/set/your/path/coverage_ratios/scaffolds_pre.txt")
	
	new_file=open("/set/your/path/coverage_ratios/scaffolds_def.txt", "w")
	
	for line in pmau_file:
		
		if re.search(r"(scf.+)\t.+\t.+\t(.+)\t.+\n", line):
			
			ind_line=re.search(r"(scf.+)\t.+\t.+\t(.+)\t.+\n", line)
			scaffold_id=ind_line.group(1)
			
			if scaffold_id in scaffold_list:
			
				scaffold_coverage=ind_line.group(2)
				
				for line2 in list_file:
					
					if re.search(r"(scf[0123456789]{13})", line2):
						
						ind_line2=re.search(r"(scf[0123456789]{13})", line2)
						scaffold_id2=ind_line.group(1)
						
						if scaffold_id2 == scaffold_id and float(scaffold_coverage) > 0:
							
							new_file.write(str(scaffold_coverage) + "\t" + line2)
						
						break
						
	pmau_file.close()
	list_file.close()
	new_file.close()		
	
	os.rename("/set/your/path/coverage_ratios/scaffolds_def.txt", "/set/your/path/coverage_ratios/scaffolds_pre.txt")

# 3 Repeat this function for all individuals (and previously identify all individuals)

ind_file=open("/users-d3/jferrer/gizquierdo/TFG/pmau_popgen_metadata.tsv")

ind_list=[]

for line in ind_file:
	if re.search(r"(.+)\t.+\t.+\.fastq\.gz\t.+\.fastq", line):
		ind_name=re.search(r"(.+)\t.+\t.+\.fastq\.gz\t.+\.fastq", line)
		ind_list.append(ind_name.group(1))
		
if len(ind_list) > 20:

	print ("Ind. list correctly assembled\n") 
	

for ind in ind_list:

	data_extraction(ind)
	
	print ("Data extraction on ind. " + ind + " has been completed\n" + str(datetime.now()))
	

ind_file.close()

# 4 Sex ratio calcul

list_file=open("/set/your/path/coverage_ratios/scaffolds_pre.txt")
new_file=open("/set/your/path/coverage_ratios/scaffolds_def.txt", "w")

for line in list_file:
	
	value=line.split("\t")
	
	scf_id=re.search(r"(scf[0123456789]{13})", line)
	
	male_mean=(float(value[33])+float(value[30])+float(value[29])+float(value[28])+float(value[26])+float(value[24])+float(value[21])+float(value[18])+float(value[17])+float(value[15])+float(value[14])+float(value[13])+float(value[12])+float(value[11])+float(value[10])+float(value[8])+float(value[7])+float(value[2])+float(value[1])+float(value[0]))/20
	
	female_mean=(float(value[36])+float(value[34])+float(value[32])+float(value[31])+float(value[27])+float(value[25])+float(value[22])+float(value[20])+float(value[19])+float(value[16])+float(value[9])+float(value[6])+float(value[5])+float(value[4])+float(value[3]))/15
	
	sex_ratio=female_mean/male_mean
	
	new_file.write(scf_id.group(1) + "\t" + str(sex_ratio) + "\n")

list_file.close()
new_file.close()

print ("Coverage-ratio calculated: " + str(datetime.now()))

print(datetime.now())

	
