# Python script to calculate average XP-EHH

import os
import re
import sys

beta_file=sys.argv[1]
scf_file=sys.argv[2]

if len(sys.argv) != 3:
                print("Error: the python script requires just one input value")
                sys.exit()

scffile=open(scf_file)														# scf_file must contain only the scf, start and stop columns
												# beta_file consists of the 3 columns used as output in Beta
outfile=open("temp_2","w")
outfile.write("scf\tstart\tstop\tavg_beta\n")											# Write the header of the output_file

# First choose the first scf

first_scf=""
first_start=0
first_stop=first_start+25000

for line in scffile:
	
	if re.search(r".+(scf[1234567890]{13}).+\t(.+)\t(.+)\n",line):							# Exclude header
		
		exp=re.search(r".+(scf[1234567890]{13}).+\t(.+)\t(.+)\n",line)						# Remember to exclude commas (if there are any) from your scf_name
		
		scf_name=exp.group(1)
		start=float(exp.group(2))
		stop=float(exp.group(3))

		first_scf=scf_name
		break

scffile.close()

# And now calulate and write down averages simultaneously:			
		
betafile=open(beta_file)

value_sum=0
count_sum=0
		
for line2 in betafile:
	
	if re.search(r"(.+)\t(.+)\t(.+)\n",line2):
		
		exp=re.search(r"(.+)\t(.+)\t(.+)\n",line2)
				
		scf=exp.group(1)
		position=exp.group(2)
				
		if position != "NA" and position != "Position" and exp.group(3) != "Beta1":		# Exclude header and use only the desired positions
	
			value=float(exp.group(3))
			
			if scf==first_scf and float(position) > first_start and float(position) <= first_stop:
	
				value_sum+=value
				count_sum+=1

			else:
	
				if count_sum == 0:
					outfile.write(first_scf + "\t" + str(first_start) + "\t" + str(first_stop) + "\tNA\n")
				else:
					outfile.write(first_scf + "\t" + str(first_start) + "\t" + str(first_stop)  + "\t" + str(value_sum/count_sum) + "\n")
	
				count_sum=1
				value_sum=value
				
				if scf==first_scf:
				
					first_start=first_stop
					first_stop=first_start+25000
				
				else:
					first_scf=scf
					first_start=0
					first_stop=first_start+25000


# Repeat the later condition when the file ends

if count_sum == 0:
	outfile.write(line.rstrip() + "\tNA\n")
else:
	outfile.write(line.rstrip() + "\t" + str(value_sum/count_sum) + "\n")

betafile.close()
outfile.close()
		
os.rename("temp_2", scf_file)
