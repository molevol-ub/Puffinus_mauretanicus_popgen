# Python script to make modifications to the PopGenome output

import re
import os

# First script to modify Fst (no 0s), pi & dxy (remember to use correct window size!!!) & calculate deltapi

window_size=25

infile=open("/users-d3/jferrer/gizquierdo/TFG/genome_scans/PopGenome/"+str(window_size)+"_windows/Fst_pairwise_auto.csv")
outfile=open("/users-d3/jferrer/gizquierdo/TFG/genome_scans/PopGenome/"+str(window_size)+"_windows/autosome_statistics.csv", "w")

outfile.write("Fst\tpi_yelkouan\tpi_mauretanicus\tTajimaD_yelkouan\tTajimaD_mauretanicus\tdxy\tdeltapi\tSNP_number\tscf_name\n")

for line in infile:
	if re.search(r"[1234567890.]\t.+\t.+\t.+\t.+\t.+\t.+\t.+\n", line):
		exp=re.search(r"(.+)\t(.+)\t(.+)\t(.+\t.+)\t(.+)\t(.+\t.+\n)", line)
		Fst=exp.group(1)
		pi1=float(exp.group(2))
		pi2=float(exp.group(3))
		rest1=exp.group(4)
		dxy=float(exp.group(5))
		rest2=exp.group(6)
		deltapi=abs(pi1-pi2)
		if "-" in Fst:
			outfile.write("0\t"+str(pi1/(window_size*1000))+"\t"+str(pi2/(window_size*1000))+"\t"+rest1+"\t"+str(dxy/(window_size*1000))+"\t"+str(deltapi/(window_size*1000))+"\t"+rest2)
		else:
			outfile.write(Fst + "\t"+str(pi1/(window_size*1000))+"\t"+str(pi2/(window_size*1000))+"\t"+rest1+"\t"+str(dxy/(window_size*1000))+"\t"+str(deltapi/(window_size*1000))+"\t"+rest2)

infile.close()
outfile.close()

# Then you can add the positions of the windows and the total length of the scaffold (use this last one to sort the file by scaffold size in excel previously to plotting)

infile=open("/users-d3/jferrer/gizquierdo/TFG/genome_scans/PopGenome/"+str(window_size)+"_windows/autosome_statistics.csv")
lenfile=open("/users-d3/jferrer/gizquierdo/TFG/genome_scans/vcfs/scf_length.csv")
outfile=open("/users-d3/jferrer/gizquierdo/TFG/genome_scans/PopGenome/"+str(window_size)+"_windows/autosome_statistics_prov.csv", "w")

incontents=infile.read()
first_scf=re.search(r"name\n.+\t.+\t.+\t.+\t.+\t.+\t.+\t.+\t\"(.+)\"\n",incontents).group(1)
window_count=0

len_dict={}

for line in lenfile:
	scf=re.search(r"(.+)\t(.+)\n", line).group(1)
	length=re.search(r"(.+)\t(.+)\n", line).group(2)
	len_dict[scf]=length

outfile.write("Fst\tpi_yelkouan\tpi_mauretanicus\tTajimaD_yelkouan\tTajimaD_mauretanicus\tdxy\tdeltapi\tSNP_number\tscf_name\tstart\tstop\tscf_length\n")

infile.close()
infile=open("/users-d3/jferrer/gizquierdo/TFG/genome_scans/PopGenome/"+str(window_size)+"_windows/autosome_statistics.csv")

for line in infile:
	newline=line.rstrip()
	if re.search(r"[1234567890.]\t.+\t.+\t.+\t.+\t.+\t.+\t.+\t.+\n", line):
		scf=re.search(r"[1234567890.]\t.+\t.+\t.+\t.+\t.+\t.+\t.+\t\"(.+)\"\n", line).group(1)
		length=len_dict[scf]
		if scf==first_scf:
			window_count+=1
		else:
			window_count=1
			first_scf=scf
		outfile.write(newline + "\t" + str((window_size*1000)*(window_count-1)) + "\t" + str((window_size*1000)*window_count) + "\t"+ str(length) + "\n")


os.rename("/users-d3/jferrer/gizquierdo/TFG/genome_scans/PopGenome/"+str(window_size)+"_windows/autosome_statistics_prov.csv", "/users-d3/jferrer/gizquierdo/TFG/genome_scans/PopGenome/"+str(window_size)+"_windows/autosome_statistics.csv")
