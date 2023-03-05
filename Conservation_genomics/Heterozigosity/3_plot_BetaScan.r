#!/soft/R-4.1.1/bin/Rscript
#$ -cwd
#$ -R y
#$ -e plot_betascan.err
#$ -o plot_betascan.out
#$ -q h12.q
#$ -pe ompi128h12 2
#$ -V                    #export environment var
#$ -N plot_betascan

setwd("/home/guillem/Documentos/TFM/gencons/het/final")

library(tidyverse)
library(ggplot2)
library(grid)

results_file<-as.data.frame(read.table("het_popgenome.def.csv", sep="\t", h=T))

# Use pivot longer to put all inds in a column (check your columns)

results_file<-pivot_longer(results_file,cols=colnames(results_file)[4:37],
	names_to="IND", values_to="HET")

# Add population info using merge with a populations file

popfile <- as.data.frame(read.table("Puffinus_pops.csv", sep=",", h=T))

results_file <- merge(results_file, popfile, by="IND")

# Set the populations

Mallorca<-results_file %>% filter(POP == "Mallorca")
Pyel<-results_file %>% filter(POP == "P.yelkouan")
Menorca<-results_file %>% filter(POP == "Menorca")
Cabrera<-results_file %>% filter(POP == "Cabrera")
Pitiuses<-results_file %>% filter(POP == "Pitiuses")
Puffinus<-results_file %>% filter(POP == "P.p.puffinus")
Canariensis<-results_file %>% filter(POP == "P.p.canariensis")

#png("het_80.noPP.png", width=1000, height=1400, res=150)
#Tajima_box<-boxplot(Pyel$HET, Menorca$HET, Mallorca$HET, Cabrera$HET, Pitiuses$HET, Puffinus$HET, Canariensis$HET, las=2, outline=FALSE, 
#	col = c("#0052A2", "#789FF2", "#8B0001", "#C34632", "#FF8532"), names=c("P.yelkouan","Menorca","Mallorca","Cabrera","Pitiüses"))
#dev.off()

png("violin.trimmed.01.png", width=1500, height=1000, res=150)
violin_plot<-ggplot(results_file, aes(x=results_file$POP, y=results_file$HET, fill=results_file$POP)) + geom_violin(trim=TRUE) + geom_boxplot(width = 0.1, fill="white", outlier.shape=NA) +
	scale_fill_manual (values = c("#C34632", "#8B0001", "#789FF2", "#A3CD67", "#3E8E28", "#0052A2", "#FF8532")) + theme_classic() +
	scale_x_discrete(limits = c("P.yelkouan","Menorca","Mallorca","Cabrera","Pitiuses", "P.p.canariensis", "P.p.puffinus")) + theme(legend.position="none") +
	labs(x="Population", y= "Heterozygosity by individual") + ylim(0,0.009)
violin_plot
dev.off()

#-------------------------------------------------------------------------------------------------------------------

# For individuals

png("violin_80.inds.trimmed.png", width=2500, height=1000, res=150)
violin_plot<-ggplot(results_file, aes(x=results_file$IND, y=results_file$HET, fill=results_file$POP)) + geom_violin(trim=TRUE) + geom_boxplot(width = 0.1, fill="white", outlier.shape=NA) +
	scale_fill_manual (values = c("#C34632", "#8B0001", "#789FF2", "#A3CD67", "#3E8E28", "#0052A2", "#FF8532")) + theme_classic() + theme(legend.position="none") +
	labs(x="Population", y= "Heterozygosity by individual") + #ylim(0,0.009) +
	scale_x_discrete(limits = c("COP1","LT2","ALT78","CZA11","ILA13","ILA2","PORQ","TZE1","M6","M8","G14","G15","M16","G4","G3","M12","M13","M14","M18","M19","M20","M1","M2","M3","M4","M5","G9","G10","G11","G12","M10","M11","M17","M21"))
violin_plot
dev.off()

#--------------------------------------------------------------------------------------------------------------------

# Barplots of each individual (average)

barplot_df <- aggregate (HET ~ IND + POP,
				data=results_file,
				FUN = mean)

png("average_barplot.png", width=2500, height=1000, res=150)
violin_plot<-ggplot(barplot_df, aes(x=IND, y=HET, fill=POP)) + geom_bar(stat="identity") +
	scale_fill_manual (values = c("#C34632", "#8B0001", "#789FF2", "#A3CD67", "#3E8E28", "#0052A2", "#FF8532")) + theme_classic() + theme(legend.position="none") +
	labs(x="Population", y= "Heterozygosity by individual") + 
	geom_text(aes(label=round(HET, digits = 5)), vjust=-1.1, hjust=0.15, size=3, angle=25) +
	scale_x_discrete(limits = c("COP1","LT2","ALT78","CZA11","ILA13","ILA2","PORQ","TZE1","M6","M8","G14","G15","M16","G4","G3","M12","M13","M14","M18","M19","M20","M1","M2","M3","M4","M5","G9","G10","G11","G12","M10","M11","M17","M21"))
violin_plot
dev.off()

