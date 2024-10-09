#!/usr/bin/python

# Calculate cumulative scaffold length of those scaffolds belonging to chrZ: we considered those with over 65% coverage_ratio and 60% hit_ratio

# In the case of autosomes these thresholds are: 80% - 120% coverage ratio & 20% hit-ratio

import re

# 1 Choose the scaffolds

scaffolds_file=open("/set/your/path/vcfs/auto/selected_scaffolds.auto.csv")

scaffold_list=[]

for line in scaffolds_file:

	if re.search(r"(scf[0123456789]{13}),", line):
		expression=re.search(r"(scf[0123456789]{13}),", line)
		scf_id=expression.group(1)
		
		scaffold_list.append(scf_id)

scaffolds_file.close()

# 2 Compare with the penguin genome

scaffolds_file=open("/set/your/path/vcfs/auto/selected_scaffolds_sphen.auto.csv")

scaffold_def=[]

for line in scaffolds_file:

	if re.search(r"(scf[0123456789]{13}),", line):
		expression=re.search(r"(scf[0123456789]{13}),", line)
		scf_id=expression.group(1)
		
		if scf_id in scaffold_list:
			
			scaffold_def.append(scf_id)

scaffolds_file.close()

# 3 Extract and sum scaffold length

info_file=open("/set/your/path/sphen_genome/final_results.1000.txt")
bed_file=open("/set/your/path/vcfs/auto/auto.bed", "w")

for line in info_file:
	
	expression=re.search(r"(scf[0123456789]{13})\t(.+)\t(.+)\t(.+)\t(.+)\n", line)
	
	scf_id=expression.group(1)
	scf_length=int(expression.group(3))
	
	if scf_id in scaffold_list:
	
		bed_file.write(scf_id + "\t0\t" + str(scf_length) + "\n")
	
info_file.close()
bed_file.close()

print("The total number of scaffolds identified as chrZ is: " + str(len(scaffold_def)))
