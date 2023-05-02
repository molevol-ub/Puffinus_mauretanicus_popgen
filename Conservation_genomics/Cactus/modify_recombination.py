import re

infile=open("Puffinus.auto.nPP.maxmiss80.recomrates.txt")
outfile=open("Puffinus.auto.nPP.maxmiss80.recomrates_positive.txt","w")

last_pos_chr=1
first_pos_chr=0
current_pos=1

for line in infile:
  if "start.pos" in line:
    outfile.write(line)
  else:
    exp=re.search(r"(.+) (.+)\n", line)
    pos=int(exp.group(1))
    rec=exp.group(2)
    
    
    if "-9" in line:
      current_pos=first_pos_chr+pos
      outfile.write(str(current_pos) + " 0.00000001\n")
      last_pos_chr=current_pos
    
    else:
      if last_pos_chr > first_pos_chr:
        first_pos_chr=last_pos_chr
        current_pos=first_pos_chr+pos                                           # This remedy is to avoid altering positions of the 1st chromosome
        outfile.write(str(current_pos) + " " + rec + "\n")

      else:
        first_pos_chr=last_pos_chr+1000000
        current_pos=first_pos_chr+pos                                            # After a change of chromosome, add 1 million bp to the previous position (last SNP of the previous chromosome)
        outfile.write(str(current_pos) + " " + rec + "\n")

      
infile.close()
outfile.close()
