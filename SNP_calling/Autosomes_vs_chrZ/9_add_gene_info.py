#!/usr/bin/python

# Script to add information to each scaffold: nº genes (Z and auto)

import re
import os

# 1 First, make dictionaries of scaffolds-gene_length for Z and auto

scaffold_list_z=[]

z_file=open("/set/your/path/annotation_summary.txt")

z_count={}
z_dict={}

for line in z_file:
	if re.search(r"(scf[1234567890]{13})", line):
		exp=re.search(r"(scf[1234567890]{13})\t(.+)\t(.+)\n", line)
		scf_id=exp.group(1)
		scf_count=exp.group(2)
		scf_length=exp.group(3)
		
		z_count[scf_id]=scf_count
		z_dict[scf_id]=scf_length
		scaffold_list_z.append(scf_id)

z_file.close()

# 1.2 Now for autosomes

scaffold_list_auto=[]

auto_file=open("/set/your/path/annotation_summary.auto.txt")

auto_count={}
auto_dict={}

for line in auto_file:
	if re.search(r"(scf[1234567890]{13})", line):
		exp=re.search(r"(scf[1234567890]{13})\t(.+)\t(.+)\n", line)
		scf_id=exp.group(1)
		scf_count=exp.group(2)
		scf_length=exp.group(3)
		
		auto_count[scf_id]=scf_count
		auto_dict[scf_id]=scf_length
		scaffold_list_auto.append(scf_id)

auto_file.close()

# 2 Add to the info file

new_file=open("/set/your/path/scf_info/scf_info_def2.txt", "w")
scf_file=open("/set/your/path/scf_info/scf_info_def.txt")

for line in scf_file:
	if re.search(r"(scf[1234567890]{13})", line):
		newline=line.rstrip()
		exp=re.search(r"(scf[1234567890]{13})", line)
		scf_id=exp.group(1)

		if scf_id in scaffold_list_z:
			new_file.write(newline + "\t" + z_count[scf_id] + "\t" + z_dict[scf_id] + "\n")

		else:
			new_file.write(newline + "\t0\t0\n")
			

scf_file.close()
new_file.close()
os.rename("/set/your/path/scf_info/scf_info_def2.txt", "/set/your/path/scf_info/scf_info_def.txt")

# 2.2 Now for autosomes

new_file=open("/set/your/path/scf_info/scf_info_def2.txt", "w")
scf_file=open("/set/your/path/scf_info/scf_info_def.txt")

new_file.write("ID\tLength\tMasked_ratio\tCoverage_ratio\tHit_ratio_c\tHit_ratio_t\tSelection_c\tSelection_t\tConclusion\tGene_count_z\tGene_length_z\tGene_count_auto\tGene_length_auto\n")

for line in scf_file:
	if re.search(r"(scf[1234567890]{13})", line):
		newline=line.rstrip()
		exp=re.search(r"(scf[1234567890]{13})", line)
		scf_id=exp.group(1)

		if scf_id in scaffold_list_auto:
			new_file.write(newline + "\t" + auto_count[scf_id] + "\t" + auto_dict[scf_id] + "\n")

		else:
			new_file.write(newline + "\t0\t0\n")

scf_file.close()
new_file.close()

os.rename("/set/your/path/scf_info/scf_info_def2.txt", "/set/your/path/scf_info/scf_info_def.txt")
