# Python script to add separation between chromosomes in the recombination input file (add -9 to the last position of each chromosome)

beta_file=sys.argv[1]
scf_file=sys.argv[2]

import re
infile=open(beta_file)
outfile=open(scf_file, "w")

count=0
first_line=""
first_recomb="0.00000001"

for line in infile:
  if "start" not in line:
    exp=re.search(r"(.+) (.+)\n", line)
    pos=int(exp.group(1))
    recomb=exp.group(2)
    
    if pos > count:                                                  # If the current position is smaller than the preivous one...
      if count != 0:
        outfile.write(str(count) + " " + first_recomb + "\n")
    else:
      outfile.write(str(count) + " -9\n")                             # ... add the -9 to the previous position
    
    count=pos
    first_line=line
    first_recomb=recomb
    
  else:
    outfile.write(line)           # Write the header
    

infile.close()
outfile.close()

# This script does not consider the last line of the recombination file, which must be copied manually
