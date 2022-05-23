#Python script to generate the input csv from where you can modify the order of individuals and generate a heatmap

import os

infile=open("pmau.100.minmac2.auto.f4.out")
outfile=open("prova.csv", "w")

char1="Initial K"
char2="lik"

for line in infile:
	if char1 in line:
		break
	else:
		continue

for line in infile:
	if char2 in line:
		break
	else:
		outfile.write(line)
		
outfile.close()
infile.close()

os.rename("prova2.csv", "prova.csv")

# And then edit it in "excel" and replace all distances between the same individual with 0.0001
