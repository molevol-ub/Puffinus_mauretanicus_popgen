# R script to summarise the results of snpEff

library(dplyr)
library(tidyr)
library(ggplot2)

setwd("/home/guillem/Documentos/TFM/conservation_genomics/snpEff/def")
input <- as.data.frame(read.table("outfile.no_MODIFIER.txt", sep="\t", header=TRUE))

good_cols <- colnames(input)[!colnames(input) %in% c("CHROM","POS","INFO","RISK","GENE")]
hola <- pivot_longer(input, cols = good_cols, names_to="IND", values_to="GENOTYPE")

# And now count: 

count_mutation <- hola[!hola$GENOTYPE %in% c("0/0","./."),] %>% count(RISK,GENOTYPE,IND)

HIGH <- count_mutation[count_mutation$RISK=="HIGH",]
LOW <- count_mutation[count_mutation$RISK=="LOW",]
MODERATE <- count_mutation[count_mutation$RISK=="MODERATE",]

HIGH$POP<- c("Pyel", "Pyel", "Cab", "Cab", "Pit", "Men", "Men", "Men", "Mal", "Cab", "Pyel", "Pyel", "Mal", "Pit", "Pit", "Mal", "Mal", "Mal", "Men", "Pit", "Mal", "Mal", "Mal", "Mal", "Pit", "Mal", "Mal", "Mal", "Men", "Men", "Pyel", "Pyel")

# And plot:

png("snpEff.mac2.HIGH.noPP.png", width=2500, height=1000, res=150)
violin_plot<-ggplot(HIGH, aes(x=IND, y=n, fill=GENOTYPE)) + geom_bar(stat="identity") +
	scale_fill_brewer (palette="Paired") + theme_classic() +
	labs(x="Individual", y= "# HIGH risk mutations") + 
	scale_x_discrete(limits = c("ALT78","CZA11","ILA13","ILA2","PORQ","TZE1","M6","M8","G14","G15","M16","G4","G3","M12","M13","M14","M18","M19","M20","M1","M2","M3","M4","M5","G9","G10","G11","G12","M10","M11","M17","M21"))
violin_plot
dev.off()


png("snpEff.mac2.MODERATE.noPP.png", width=2500, height=1000, res=150)
violin_plot<-ggplot(MODERATE, aes(x=IND, y=n, fill=GENOTYPE)) + geom_bar(stat="identity") +
	scale_fill_brewer (palette="Paired") + theme_classic() +
	labs(x="Individual", y= "# MODERATE risk mutations") + 
	scale_x_discrete(limits = c("ALT78","CZA11","ILA13","ILA2","PORQ","TZE1","M6","M8","G14","G15","M16","G4","G3","M12","M13","M14","M18","M19","M20","M1","M2","M3","M4","M5","G9","G10","G11","G12","M10","M11","M17","M21"))
violin_plot
dev.off()


png("snpEff.mac2.LOW.noPP.png", width=2500, height=1000, res=150)
violin_plot<-ggplot(LOW, aes(x=IND, y=n, fill=GENOTYPE)) + geom_bar(stat="identity") +
	scale_fill_brewer (palette="Paired") + theme_classic() +
	labs(x="Individual", y= "# LOW risk mutations") + 
	scale_x_discrete(limits = c("ALT78","CZA11","ILA13","ILA2","PORQ","TZE1","M6","M8","G14","G15","M16","G4","G3","M12","M13","M14","M18","M19","M20","M1","M2","M3","M4","M5","G9","G10","G11","G12","M10","M11","M17","M21"))
violin_plot
dev.off()
