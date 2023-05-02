import re

infile=open("Puffinus.auto.nPP.maxmiss80.recomrates.txt")
outfile=open("Puffinus.auto.nPP.maxmiss80.recomrates_positive.txt","w")

last_pos=0
last_pos_chr=0
count_position=0

for line in infile:
  if "start.pos" in line:
    outfile.write(line)
  else:
    exp=re.search(r"(.+)\t(.+)\n", line)
    pos=int(exp.group(1))
    rec=int(exp.group(2))
    
    
    if "-9" in line:
      outfile.write(str(pos) + "\t0.00000001\n")
      last_pos_chr=pos
      last_pos=pos
    
    else:
      if last_pos==last_pos_chr:
        last_pos=pos+1000000                                            # After a change of chromosome, add 1 million bp to the previous position (last SNP of the previous chromosome)
        outfile.write(str(last_pos) + "\t" + str(rec) + "\n")
      
      else:
        last_pos+=pos                                                   # In the rest of lines sum to the previous one, which already has the milion bp incorporated
        outfile.write(str(last_pos) + "\t" + str(rec) + "\n")
      
infile.close()
outfile.close()
