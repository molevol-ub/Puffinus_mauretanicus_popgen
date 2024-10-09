#!/bin/python

# Script to add information to each scaffold: coverage_ratio, hit_ratios, selection and nº genes

# We will use the following selection code: A - chromosome Z; B - coverage ratio ~ 75% (not confidently assigned to either); C - coverage ratio ~ 100% but not syntenic (not confidently assigned to either category); D - Autosomes; n - not assigned to any category

import re
import os

# 1 Adding coverage_ratio data:

scf_file=open("/set/your/path/scf_info/scf_info.txt")
coverage_file=open("/set/your/path/scf_info/coverage_ratios/scaffolds_def.txt")
new_file=open("/set/your/path/scf_info/scf_info_def.txt", "w")

coverage_dict={}
scaffold_list=[]

for line in coverage_file:
	if re.search(r"(scf[1234567890]{13})", line):
		exp=re.search(r"(scf[1234567890]{13})\t(.+)\n", line)
		scf_id=exp.group(1)
		scf_cov=exp.group(2)
		
		coverage_dict[scf_id]=scf_cov
		scaffold_list.append(scf_id)

for line in scf_file:

	if re.search(r"(scf[1234567890]{13})", line):

		newline=line.rstrip()
		
		exp=re.search(r"(scf[1234567890]{13})", line)
		scf_id=exp.group(1)
		
		if scf_id in scaffold_list:
			
			new_file.write(newline + "\t" + coverage_dict[scf_id] + "\n")
		else:
			
			new_file.write(newline + "\t0\n")

new_file.close()
coverage_file.close()
scf_file.close()

# 2 Adding hit_ratio data:

def hit_extraction(file):
	
	scf_file=open("/set/your/path/scf_info/scf_info_def.txt")
	new_file=open("/set/your/path/scf_info/scf_info_def2.txt", "w")
	in_file=open("/set/your/path/"+file+"_genome/chrz_final_results.1000.txt")

	hit_dict={}
	scaffold_list=[]

	for line in in_file:
		if re.search(r"(scf[1234576890]{13})", line):
			exp=re.search(r"(scf[1234576890]{13})\t.+\t.+\t(.+)\t.+\n", line)
			scf_id=exp.group(1)
			scf_hit=exp.group(2)

			hit_dict[scf_id]=scf_hit
			scaffold_list.append(scf_id)

	for line in scf_file:
		if re.search(r"(scf[1234567890]{13})", line):

			newline=line.rstrip()
			
			exp=re.search(r"(scf[1234567890]{13})", line)
			scf_id=exp.group(1)

			if scf_id in scaffold_list:
			
				new_file.write(newline + "\t" + hit_dict[scf_id] + "\n")
			else:
				
				new_file.write(newline + "\tna\n")

	scf_file.close()
	new_file.close()
	in_file.close()

	os.rename("/set/your/path/scf_info/scf_info_def2.txt", "/set/your/path/scf_info/scf_info_def.txt")

hit_extraction("cic")
hit_extraction("sphen")

# 3 Add selection info

def selection(work_file):
	
	A_list=[]
	B_list=[]
	C_list=[]
	D_list=[]

	A_select_file=open("/set/your/path/" + work_file + "_consensus.A.csv")
	for line in A_select_file:
		scf_id=re.search(r"(scf[1234567890]{13})", line).group(1)
		A_list.append(scf_id)
	
	B_select_file=open("/set/your/path/" + work_file + "_consensus.B.csv")
	for line in B_select_file:
		scf_id=re.search(r"(scf[1234567890]{13})", line).group(1)
		B_list.append(scf_id)
	
	C_select_file=open("/set/your/path/" + work_file + "_consensus.C.csv")
	for line in C_select_file:
		scf_id=re.search(r"(scf[1234567890]{13})", line).group(1)
		C_list.append(scf_id)
	
	D_select_file=open("/set/your/path/" + work_file + "_consensus.D.csv")
	for line in D_select_file:
		scf_id=re.search(r"(scf[1234567890]{13})", line).group(1)
		D_list.append(scf_id)
	
	scf_file=open("/set/your/path/scf_info/scf_info_def.txt")
	new_file=open("/set/your/path/scf_info/scf_info_def2.txt", "w")
	
	
	for line in scf_file:
		if re.search(r"(scf[1234567890]{13})", line):

			newline=line.rstrip()
			
			exp=re.search(r"(scf[1234567890]{13})", line)
			scf_id=exp.group(1)

			if scf_id in A_list:

				new_file.write(newline + "\tA\n")

			elif scf_id in B_list:

				new_file.write(newline + "\tB\n")

			elif scf_id in C_list:

				new_file.write(newline + "\tC\n")

			elif scf_id in D_list:

				new_file.write(newline + "\tD\n")

			else:
				
				new_file.write(newline + "\tn\n")

	os.rename("/set/your/path/scf_info/scf_info_def2.txt", "/set/your/path/scf_info/scf_info_def.txt")


selection("sphen")
