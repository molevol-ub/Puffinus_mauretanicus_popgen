# R script to summarise the results of snpEff

library(dplyr)
library(tidyr)

setwd("/home/guillem/Documentos/TFM/conservation_genomics/snpEff")
input <- as.data.frame(read.table("outfile.no_MODIFIER.txt", sep="\t", header=TRUE))

good_cols <- colnames(input)[!colnames(input) %in% c("CHROM","POS","INFO","RISK","GENE")]
hola <- pivot_longer(input, cols = good_cols, names_to="IND", values_to="GENOTYPE")

# And now count: 

number_LOF <- input %>% count(V4)

# In function of genotype:

number_LOF <- input %>% count(RISK,IND,GENOTYPE)

# Or characterizing only HIGH risk variants

number_LOF <- input[input$RISK=="HIGH",] %>% count(IND,GENOTYPE)
