# R script to summarise the results of snpEff

library(dplyr)
library(tidyr)

setwd("/home/guillem/Documentos/TFM/conservation_genomics/snpEff")
input <- as.data.frame(read.table("outfile.no_MODIFIER.txt", sep="\t", header=TRUE))

good_cols <- colnames(input)[!colnames(input) %in% c("CHROM","POS","INFO","RISK","GENE")]
hola <- pivot_longer(input, cols = good_cols, names_to="IND", values_to="GENOTYPE")

# And now count: 

count_mutation <- hola[!hola$GENOTYPE %in% c("0/0","./."),] %>% count(RISK,GENOTYPE,IND)

HIGH <- count_mutation[count_mutation$RISK=="HIGH",]
LOW <- count_mutation[count_mutation$RISK=="LOW",]
MODERATE <- count_mutation[count_mutation$RISK=="MODERATE",]
