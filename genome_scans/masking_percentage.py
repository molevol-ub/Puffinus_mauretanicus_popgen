# Python script to calculate the masking percentages for each window

import sys
window_size=sys.argv[1]

if len(sys.argv) != 2:
	print("Error: the python script requires just one input value")
	sys.exit()

import os
import re

infile=open("/users-d3/jferrer/gizquierdo/TFG/genome_scans/PopGenome/"+str(window_size)+"_windows/masking_amount.bed")
outfile=open("/users-d3/jferrer/gizquierdo/TFG/genome_scans/PopGenome/"+str(window_size)+"_windows/masking_prop.csv", "w")

outfile.write("Unmasked proportion(%)\tScaffold\tStart\n")

count=0

first_start=0
first_scf="scf7180000012825"

for line in infile:
	exp=re.search(r"(.+)\t(.+)\t(.+)\t.+\t.+\t.+\t.+\t(.+)\n", line)
	scf=exp.group(1)
	start=int(exp.group(2))
	stop=int(exp.group(3))
	mask=int(exp.group(4))

	if start == first_start and scf==first_scf:
		count+= mask

	elif start != first_start and scf == first_scf:
		outfile.write(str(count/(1000*float(window_size))) + "\t" + str(first_scf) + "\t" + str(first_start) + "\n")
		first_start=start
		count=0
		count+=mask

	elif scf != first_scf:
		outfile.write(str(count/(0.001*float(window_size))) + "\t" + str(first_scf) + "\t" + str(first_start) + "\n")
		first_scf=scf
		first_start=0
		count=0
		count+=mask

infile.close()
outfile.close()

# Now add to the autosome statistics file

stat_file=open("/users-d3/jferrer/gizquierdo/TFG/genome_scans/PopGenome/"+str(window_size)+"_windows/autosome_statistics.csv")
infile=open("/users-d3/jferrer/gizquierdo/TFG/genome_scans/PopGenome/"+str(window_size)+"_windows/masking_prop.csv")
outfile=open("/users-d3/jferrer/gizquierdo/TFG/genome_scans/PopGenome/"+str(window_size)+"_windows/autosome_statistics_def.csv", "w")

info_list=[]
info_dict={}

for line in infile:
	exp=re.search(r"(.+)\t(.+\t.+)\n", line)
	prop=exp.group(1)
	scf=exp.group(2)
	info_dict[scf]=prop
	info_list.append(scf)

for line in stat_file:
	if re.search(r"scf[1234567890]{13}", line):
		exp1=re.search(r".+\t.+\t.+\t.+\t.+\t.+\t.+\t.+\t\"(.+)\"(\t.+)\t.+\t.+\n", line).group(1)
		exp2=re.search(r".+\t.+\t.+\t.+\t.+\t.+\t.+\t.+\t\"(.+)\"(\t.+)\t.+\t.+\n", line).group(2)
		exp=exp1+exp2
		if exp in info_list:
			outfile.write(line.rstrip() + "\t" + str(info_dict[exp]) + "\n")
		else:
			outfile.write(line.rstrip() + "\tNA\n")		
	else:
		outfile.write(line.rstrip() + "\tUnmasked_proportion\n")
		
stat_file.close()
infile.close()
outfile.close()

os.rename ("/users-d3/jferrer/gizquierdo/TFG/genome_scans/PopGenome/"+str(window_size)+"_windows/autosome_statistics_def.csv", "/users-d3/jferrer/gizquierdo/TFG/genome_scans/PopGenome/"+str(window_size)+"_windows/autosome_statistics.csv")
