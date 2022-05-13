# Python script to calculate average XP-EHH

import os
import re
import sys

input_file=sys.argv[1]
stat_file=sys.argv[2]

if len(sys.argv) != 3:
                print("Error: the python script requires just one input value")
                sys.exit()

workfile=open(input_file)
outfile=open("temp", "w")
outfile2=open("xpehh", "w")

first_scf=""
first_start=0
first_stop=0
count=0
xpehh_count=0

# First choose the first scf

for line in workfile:
	scf=re.search(r"(.+)\t(.+)\t(.+)\t(.+)\t(.+)\t(.+)\t(.+)\n", line).group(1)
	start=float(re.search(r"(.+)\t(.+)\t(.+)\t(.+)\t(.+)\t(.+)\t(.+)\n", line).group(2))
	stop=float(re.search(r"(.+)\t(.+)\t(.+)\t(.+)\t(.+)\t(.+)\t(.+)\n", line).group(3))
	first_scf=scf
	first_start=start
	first_stop=stop
	break

workfile.close()
workfile=open(input_file)


for line in workfile:
	scf=re.search(r"(.+)\t(.+)\t(.+)\t(.+)\t(.+)\t(.+)\t(.+)\n", line).group(1)
	start=float(re.search(r"(.+)\t(.+)\t(.+)\t(.+)\t(.+)\t(.+)\t(.+)\n", line).group(2))
	stop=float(re.search(r"(.+)\t(.+)\t(.+)\t(.+)\t(.+)\t(.+)\t(.+)\n", line).group(3))
	
	xpehh=float(re.search(r"(.+)\t(.+)\t(.+)\t(.+)\t(.+)\t(.+)\t(.+)\n", line).group(7))

	if scf == first_scf and first_start == start and first_stop == stop:

		count+=1
		xpehh_count+=xpehh

	else:

		xpehh_value=xpehh_count/float(count)
		outfile.write("\"" + first_scf + "\"\t" + str(int(first_start)) + "\t" + str(int(first_stop)) + "\n")
		outfile2.write(str(xpehh_value) + "\n")

		count=0
		xpehh_count=0

		first_scf=scf
		first_stop=stop
		first_start=start

		count+=1
		xpehh_count+=xpehh

# Repeat the later condition when the file ends

xpehh_value=xpehh_count/float(count)
outfile.write("\"" + first_scf + "\"\t" + str(int(first_start)) + "\t" + str(int(first_stop)) + "\n")
outfile2.write(str(xpehh_value) + "\n")

# Rename file

workfile.close()
outfile.close()
outfile2.close()

os.rename ("temp", input_file)

# Finally add to the general statistics file

outfile=open("temp", "w")
final_input=open(stat_file)

infile=open(input_file)
incontents=infile.read()
inlist=incontents.split("\n")

xpehh_file=open("xpehh")
xpehh_cont=xpehh_file.read()
xpehh_list=xpehh_cont.split("\n")

dict={}

for num in range(0,len(inlist)):
	dict[inlist[num]]=xpehh_list[num]

for line in final_input:
	if re.search(r"scf[1234567890]{13}", line):
		exp=re.search(r"(\"scf[1234567890]{13}\"\t.+\t.+)\t.+\t.+\t.+\t.+\n", line).group(1)
		if exp in inlist:
			newline=line.rstrip()
			outfile.write(newline + "\t" + dict[exp] + "\n")
		else:
                        newline=line.rstrip()
                        outfile.write(newline + "\tNA\n")

	else:
		newline=line.rstrip()
		outfile.write(newline + "\tavg_XP-EHH\n")

outfile.close()
infile.close()
xpehh_file.close()
final_input.close()

os.rename("temp", input_file)
