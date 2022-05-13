# Python script to calculate gene proportion

import os
import re
import sys

input_file=sys.argv[1]
window_size=int(sys.argv[2])

if len(sys.argv) != 3:
		print("Error: the python script requires just 2 input values")
		sys.exit()

infile=open(input_file)
outfile=open("temp", "w")

first_scf=""
first_opening=0
hit_amount=0

for line in infile:
	first_scf=re.search(r"(scf[1234567890]{13})", line).group(1)
	first_opening=int(re.search(r"(.+)\t(.+)\t.+\t.+\t.+\t.+\t(.+)\n", line).group(2))
	break

infile.close()
infile=open(input_file)

for line in infile:

	scf=re.search(r"(.+)\t(.+)\t.+\t.+\t.+\t.+\t(.+)\n", line).group(1)
	start=int(re.search(r"(.+)\t(.+)\t.+\t.+\t.+\t.+\t(.+)\n", line).group(2))
	local_hit=int(re.search(r"(.+)\t(.+)\t.+\t.+\t.+\t.+\t(.+)\n", line).group(3))

	if first_scf == scf and start==first_opening:
		hit_amount+=local_hit

	else:
		if hit_amount > (window_size*1000):
			hit_amount = window_size*1000			# Just in case there are overlapping genes
		
		outfile.write(str(float(hit_amount)/float(window_size*1000)) + "\n")

		first_scf=re.search(r"(scf[1234567890]{13})", line).group(1)
		first_opening=start
		hit_amount=local_hit

infile.close()
outfile.close()

os.rename ("temp", input_file)
