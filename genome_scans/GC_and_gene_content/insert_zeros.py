# Python script to calculate gene proportion

import os
import re
import sys

input_file=sys.argv[1]

if len(sys.argv) != 2:
		print("Error: the python script requires just one input value")
		sys.exit()

infile=open(input_file)
outfile=open("temp", "w")

for line in infile:
	name=re.search(r"(.+)\t(.+)\n", line).group(1)
	length=re.search(r"(.+)\t(.+)\n", line).group(2)
	
	outfile.write(name + "\t0\t" + str(length) + "\n")
	
os.rename ("temp", input_file)
