#!/bin/python

# Script to calculate the scaffold length and the percentage of masked regions (%)

import re

# 1 Locate files

genome_file=open("/set/your/path/genome_masked.fa")
scf_file=open("//set/your/path/scf_info/scf_info.txt", "w")

scf_file.write("Scf_id\tScf_length\tMask_ratio(%)\n")

char_count=0
mask_count=0
char_list=["A","C","T","G","N"]

# 2 Count both all and masked length

for line in genome_file:

	if re.search(r"(scf[1234567890]{13})", line):

		expression = re.search(r"(scf[1234567890]{13})", line)
		scf_id=expression.group(1)

		if char_count != 0:

			mask_ratio=(100*mask_count/char_count)

			scf_file.write(previous_scf+"\t" + str(char_count) +"\t" + str(mask_ratio) + "\n")

		char_count=0
		mask_count=0

		previous_scf=scf_id


	else:

		for char in line:

			if char in char_list:
				
				char_count+=1

			if char == "N":

				mask_count +=1

scf_file.close()
genome_file.close()
