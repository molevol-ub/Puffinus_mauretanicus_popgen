#!/soft/R-4.1.1/bin/Rscript

setwd("/home/guillem/Documentos/TFM/conservation_genomics/phyloP")

library(tidyverse)
library(ggplot2)
library(grid)

results_file<-as.data.frame(read.table("Puffinus.load.count.def.csv", sep=",", h=T))
pop_file <- as.data.frame(read.table("Puffinus_list.csv", sep=",", h=T))
results_file <- merge(results_file, pop_file, by="IND")

pdf("phyloP.wPP.def.pdf", width=22, height=10)
manhplot <- ggplot(results_file, aes(x = IND, y=COUNT, fill=GENOTYPE, col=GENOTYPE)) + 
	geom_bar(stat="identity") +
	theme_classic() +
	theme (legend.title = element_text(size=15), legend.text = element_text(size=12)) +
	scale_colour_brewer (palette="Blues", name = "GENOTYPE") + 
	scale_fill_brewer (palette="Blues", name = "GENOTYPE") +
	labs(x="Individual", y="# deleterious mutations") +
	theme (axis.text=element_text(size=12.5), axis.title=element_text(size=16)) +
	scale_x_discrete(limits = c("COP1","LT2","ALT78","CZA11","ILA13","ILA2","PORQ","TZE1","M6","M8","G14","G15","M16","G4","G3","M12","M13","M14","M18","M19","M20","M1","M2","M3","M4","M5","G9","G10","G11","G12","M10","M11","M17","M21"))
	
manhplot
dev.off()
#--------------------------------------------------------------------------------

results_file<-results_file[results_file$POP!="P.p.puffinus",]
results_file<-results_file[results_file$POP!="P.p.canariensis",]
results_file<-results_file[results_file$IND!="M19",]

# Get sum of HOM and HET for the jitter plot
sum_IND <- aggregate (COUNT ~ IND + POP, data=results_file, FUN=sum)
colnames(sum_IND)<-c("IND","POP","SUM")
results_file<-merge(results_file, sum_IND)
results_file[results_file$GENOTYPE=="HOM",]$SUM<-results_file[results_file$GENOTYPE=="HOM",]$COUNT

#Step for dots of >1Mb category to be above the average 100kb-1Mb always

provisional<-aggregate(SUM ~ POP, data=results_file[results_file$GENOTYPE=="HOM",], FUN=mean)
colnames(provisional)<-c("POP","SUM2")
results_file<-merge(results_file, provisional)
results_file[results_file$GENOTYPE=="HET",]$SUM<-results_file[results_file$GENOTYPE=="HET",]$COUNT+results_file[results_file$GENOTYPE=="HET",]$SUM2

mean_pop<-aggregate (SUM ~ POP + GENOTYPE, data=results_file, FUN=mean)
sd_pop<-aggregate (SUM ~ POP + GENOTYPE, data=results_file, FUN=sd)
colnames(mean_pop)<-c("POP", "GENOTYPE","AVG")
colnames(sd_pop)<-c("POP", "GENOTYPE" ,"SD")
results_file<-merge(merge(results_file, mean_pop),sd_pop)

# And plot:

pdf("phyloP.def.pdf", width=5, height=6.5)
manhplot <- ggplot(results_file, aes(x = POP, y=COUNT, fill=POP, alpha=GENOTYPE)) + 
	geom_bar(stat="summary") +
	theme_classic() +
	scale_fill_manual(values = c("#C34632", "#8B0001", "#789FF2", "#0052A2", "#FF8532"))+
	labs(x="Population", y="# deleterious mutations") +
	theme (axis.text=element_text(size=11.5), axis.title=element_text(size=17)) +
	scale_x_discrete(limits=c("P.yelkouan","Menorca","Cabrera","Mallorca","Pitiuses")) +
	theme(legend.position="none")+
	geom_jitter(aes(y=SUM), shape=16, position=position_jitter(0.1))+
	scale_alpha_discrete(range=c(0.45,1))+
	geom_errorbar(aes(ymin=AVG-SD, ymax=AVG+SD, alpha=GENOTYPE), width=.2)
	
manhplot
dev.off()

#--------------------------------------------------------------------------------

# To test for significant differences using Welch's 2-sided t test)

# 8. Create dataframe with sum of totals (important to include the population)

sum_IND <- aggregate (P_ROH ~ FID + POP, data=ROHs, FUN=sum)

# 9. Test for differences

Mallorca <- results_file[results_file$POP=="Mallorca",]
Menorca <- results_file[results_file$POP=="Menorca",]
Pitiuses <- results_file[results_file$POP=="Pitiuses",]
Cabrera <- results_file[results_file$POP=="Cabrera",]
Pyel <- results_file[results_file$POP=="P.yelkouan",]

t.test(Mallorca$P_ROH, Menorca$P_ROH)

# And if you want to calculate the average length per pop and category:

aggregate (P_ROH ~ POP, data=sum_IND, FUN=mean)
