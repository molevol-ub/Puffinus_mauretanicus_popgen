#!/soft/R-4.1.1/bin/Rscript

setwd("/home/guillem/Documentos/TFM/gencons/ROH")

library(tidyverse)
library(ggplot2)
library(grid)

results_file<-as.data.frame(read.table("Puffinus.window_het_5.thres_01.100kb.csv", sep=",", h=T))

# 1. Create column with ROH category

results_file$ROH_category <- ifelse(results_file$KB > 1000, ">1000kb",
				ifelse (results_file$KB <= 1000 & results_file$KB > 500, "500kb-1000kb",
				ifelse (results_file$KB <= 500 & results_file$KB > 250, "250kb-500kb", "100kb-250kb")))

# 2. Delete all columns with PHET < 0.03 (arbitrary cutoff based on observations on scaffold E, ...)

results_file <- results_file[results_file$PHET <= 0.03,]

# 3. Calculate proportion of genome of type of ROH by indv

ROHs <- aggregate(KB ~ ROH_category + FID,
			data=results_file,
			FUN=sum)  
			
# 4. Merge with a popmap that also has the coverage/mappability mask for each individual

pop_file <- as.data.frame(read.table("Puffinus_list.csv", sep="\t", h=T))

ROHs <- merge(ROHs, pop_file, by="FID")

# 5. Calculate the proportion of the unmasked genome covered by ROHs

ROHs$P_ROH <- 100*(ROHs$KB/ROHs$UNMASKED)

# 6. Convert category to factor and reorder

ROHs$ROH_category <- factor(ROHs$ROH_category, levels=c(">1000kb", "500kb-1000kb", "250kb-500kb", "100kb-250kb"))

# 6. Plot

pdf("ROHs.wPP.pdf", width=22, height=10)
manhplot <- ggplot(ROHs, aes(x = FID, y=P_ROH, fill=ROH_category, col=ROH_category)) + 
	geom_bar(stat="identity") +
	theme_classic() +
	theme (legend.title = element_text(size=15), legend.text = element_text(size=12)) +
	scale_colour_brewer (palette="YlOrRd", name = "ROH length") + 
	scale_fill_brewer (palette="YlOrRd", name = "ROH length") +
	labs(x="Individual", y="% of genome in ROHs") +
	theme (axis.text=element_text(size=12.5), axis.title=element_text(size=16)) +
	scale_x_discrete(limits = c("COP1","LT2","ALT78","CZA11","ILA13","ILA2","PORQ","TZE1","M6","M8","G14","G15","M16","G4","G3","M12","M13","M14","M18","M19","M20","M1","M2","M3","M4","M5","G9","G10","G11","G12","M10","M11","M17","M21"))
	
manhplot
dev.off()


# And if you want to calculate the average length per pop and category:

# aggregate (P_ROH ~ ROH_category + POP, data=ROHs, FUN=mean)
#--------------------------------------------------

# And if you want, make averages per population
