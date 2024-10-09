#!/usr/bin/env python

# Python script to assess the scaffold length occupied by hits (not filtered by length yet)

from datetime import datetime
import re
print (datetime.now())

# 1 Define function

def extraction (ch):
	
	sphen_file=open("/set/your/path/sphen_genome/" + ch + "_sphen.1000.txt")
	new_file=open("/set/your/path/sphen_genome/" + ch + "_hit_length.1000.txt", "w")
	
	scaffolds_expression=["0-0"]
	previous_scaffold="scf7180000012825"
	
	for line in sphen_file:
	
		expression=re.search(r"(scf.+)\t(.+)\t(.+)\t(.+)\t(.+)\t(.+)\t(.+)\t(.+)\t(.+)\t(.+)\t(.+)\t(.+)", line)
		scf_id=expression.group(1)
			
		if scf_id == previous_scaffold:
			
			initial_pos=int(expression.group(7))
			final_pos=int(expression.group(8))
					
			counter=0
			
			for expression in scaffolds_expression:
				
				scf_exp=re.search(r'(.+)-(.+)', expression)
				
				scf_in=int(scf_exp.group(1))
				scf_out=int(scf_exp.group(2))
				
				resta_1 = scf_in - initial_pos
				resta_2 = final_pos - scf_out
				
				if initial_pos in range(scf_in-1, scf_out+1):
					
					if final_pos in range (scf_in-1, scf_out+1):
						
						counter+=1
						continue
					
					else:
					
						scaffolds_expression.remove (expression)
						exp_input= str(scf_in) + "-" + str(final_pos)
						scaffolds_expression.append (exp_input)
						counter+=1
				
				elif final_pos in range (scf_in-1, scf_out+1):
				
					scaffolds_expression.remove (expression)
					exp_input= str(initial_pos) + "-" + str(scf_out)
					scaffolds_expression.append (exp_input)
					counter+=1
						
				
				elif resta_1 >=0 and resta_2 >=0:
				
					scaffolds_expression.remove (expression)
					exp_input= str(initial_pos) + "-" + str(final_pos)
					scaffolds_expression.append (exp_input)
					counter+=1
					
			if counter == 0:
			
				exp_input= str(initial_pos) + "-" + str(final_pos)
				scaffolds_expression.append (exp_input)
		
		else:
			
			scaffold_length=0
			
			for expression in scaffolds_expression:
	
				scf_exp=re.search(r'(.+)-(.+)', expression)
				
				scf_in=int(scf_exp.group(1))
				scf_out=int(scf_exp.group(2))
				
				calc=scf_out - scf_in
				
				scaffold_length += calc
	
			print("Cumulative hit length for scaffold " + previous_scaffold + " is: " + str(scaffold_length))
	
			print (datetime.now())
			new_file.write(previous_scaffold + "\t" + str(scaffold_length) + "\n")
			
			scaffolds_expression=["0-0"]
			previous_scaffold=scf_id
			
	
	sphen_file.close()
	new_file.close()
	print ("Procedure for " + ch + " finished: " + str(datetime.now()))
	
# 2 Define list of chromosomes and loop through them

chr_list=("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chr23", "chr24", "chr25", "chrz")

for chrom in chr_list:
	
	extraction(chrom)
