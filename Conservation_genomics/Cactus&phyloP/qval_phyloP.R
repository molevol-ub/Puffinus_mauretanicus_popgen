#!/usr/bin/env Rscript

setwd("/users-d3/jferrer/gizquierdo/TFM/Cactus/Cactus/PhyloP")

# Pass infile name as argument
args = commandArgs(trailingOnly=TRUE)
print(args)
infile<-read.table(args[1])

#Transform -log(pval) to pval
infile$V6<-10^(-infile$V5)

#Transform pval to qval
infile$V7<-p.adjust (infile$V6, method="fdr")

#Write_output
colnames(infile)<-c("scf", "start", "stop", "id", "-log(pval)", "pval", "qval")
write.table(infile, paste0(args[1],".csv"), quote=FALSE)
