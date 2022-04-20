# Python script to calculat pi, dxy & deltapi using the correct amount of valid positions (multiply by the inversed masked percentage)

window_size=1

import re
import os

infile=open("/users-d3/jferrer/gizquierdo/TFG/genome_scans/PopGenome/"+str(window_size)+"_windows/autosome_statistics.csv")
outfile=open("/users-d3/jferrer/gizquierdo/TFG/genome_scans/PopGenome/"+str(window_size)+"_windows/autosome_statistics_prov.csv", "w")

for line in infile:
	if re.search(r"[1234567890.]\t.+\t.+\t.+\t.+\t.+\t.+\t.+\t.+\t.+\t.+\t.+\t.+\n", line):
		exp=re.search(r"(.+)\t(.+)\t(.+)\t(.+\t.+)\t(.+)\t.+\t(.+\t.+\t.+\t.+\t.+)\t(.+)\n", line)
		rest1=exp.group(1)
		pi1=float(exp.group(2))
		pi2=float(exp.group(3))
		rest2=exp.group(4)
		dxy=float(exp.group(5))
		rest3=exp.group(6)
		prop=exp.group(7)
		if prop == "NA":
			prop = 100
		else: 
			prop = float(exp.group(7))
		outfile.write(rest1 + "\t" +str(pi1/prop)+"\t"+str(pi2/prop)+"\t"+rest2+"\t"+str(dxy/prop)+"\t"+str(abs(pi1-pi2))+"\t"+rest3+"\t"+str(prop)+"\n")
	else:
		outfile.write(line)

infile.close()
outfile.close()

os.rename ("/users-d3/jferrer/gizquierdo/TFG/genome_scans/PopGenome/"+str(window_size)+"_windows/autosome_statistics_prov.csv", "/users-d3/jferrer/gizquierdo/TFG/genome_scans/PopGenome/"+str(window_size)+"_windows/autosome_statistics.csv")
