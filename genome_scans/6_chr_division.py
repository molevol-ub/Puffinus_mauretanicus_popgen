# Python script to extract autosomal and chrZ windows

window_size=25

import re

infile=open("/users-d3/jferrer/gizquierdo/TFG/genome_scans/PopGenome/"+str(window_size)+"_windows/Fst_paiwise.csv")
auto_file=open("/users-d3/jferrer/gizquierdo/TFG/genome_annotation/consensus_scfs/common_consensus.D.csv")
auto_contents=auto_file.read()
outfile=open("/users-d3/jferrer/gizquierdo/TFG/genome_scans/PopGenome/"+str(window_size)+"_windows/Fst_pairwise_auto.csv", "w")

for line in infile:
	if re.search(r".+\t.+\t.+\t.+\t.+\t.+\t.+\n", line):
		scf=re.search(r".+\t.+\t.+\t.+\t.+\t.+\t\"(.+)\"\n", line).group(1)
		if scf in auto_contents:
			outfile.write(line)

infile.close()
auto_file.close()
outfile.close()


infile=open("/users-d3/jferrer/gizquierdo/TFG/genome_scans/PopGenome/"+str(window_size)+"_windows/Fst_paiwise.csv")
chrZ_file=open("/users-d3/jferrer/gizquierdo/TFG/genome_annotation/consensus_scfs/common_consensus.A.csv")
Z_contents=chrZ_file.read()
outfile=open("/users-d3/jferrer/gizquierdo/TFG/genome_scans/PopGenome/"+str(window_size)+"_windows/Fst_pairwise_chrZ.csv", "w")

for line in infile:
	if re.search(r".+\t.+\t.+\t.+\t.+\t.+\t.+\n", line):
		scf=re.search(r".+\t.+\t.+\t.+\t.+\t.+\t\"(.+)\"\n", line).group(1)
		if scf in Z_contents:
			outfile.write(line)

infile.close()
auto_file.close()
outfile.close()
