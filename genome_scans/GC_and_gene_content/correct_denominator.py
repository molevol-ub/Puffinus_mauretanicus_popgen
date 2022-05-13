# Python script to calculate gene proportion

import os
import re
import sys

input_file=sys.argv[1]
window_size=int(sys.argv[2])

if len(sys.argv) != 3:
                print("Error: the python script requires just one input value")
                sys.exit()

# Now correct the denominator

infile=open(input_file)
outfile=open("temp", "w")

for line in infile:
	if re.search (r"\"scf[1234567890]{13}\"\t.+\t.+\t.+\t.+\t(.+)\t(.+)\t(.+)\n", line):
		exp=re.search (r"(.+\t\"scf[1234567890]{13}\"\t.+\t.+\t.+\t.+\t)(.+)\t(.+)\t(.+)\n", line)
		rest=exp.group(1)
		unmask=float(exp.group(4))
		gc=float(exp.group(2))
		gene=float(exp.group(3))
	
		if unmask != 0:
			gc_correct=gc/unmask
			gene_correct=gene/unmask
			
			if gene_correct > 1:
				gene_correct = 1		# Those few cases of overlap (check beforehand)
			
			outfile.write(rest+str(gc_correct) + "\t" + str(gene_correct) + "\n")
	else:
		outfile.write(line)

infile.close()
outfile.close()

os.rename("temp", input_file)
