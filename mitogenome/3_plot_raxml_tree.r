#!/soft/R-4.1.1/bin/Rscript

library(vcfR)
library(ape)
#library(adegenet)
#library(phangorn)

setwd("/home/guillem/Documents/TFG/phylogenies/mitochondrial")

library(ggplot2)
library(ggtree)
library(treeio)
library (tibble)

# 1 The "support" ouput or RaxML-ng can be read with read.tree. Then, bootstrap values can be found as "label" (label=label)

tree <- read.tree("Puffinus_mtDNA_modified.raxml.support")

# 2 You could try to highlight nodes according to its bootstrap value (according to the approximate labels in a modified text file)
# It is better to group the bootstrap values into cathegories (eg. 100-90, 90-70, ...) in your nexus file

ggtree(tree) + geom_tiplab(size=3) + geom_nodepoint(aes(col=label), size = 2.5) + geom_treescale()

# 3 as well as locations 

#loc<- data.frame(node=c(1,2,8,7,11,10,4,6,9,25,28,29,5,12,32,17,20,21,22,23,3,13,14,15,16,18,19,24,26,27,30,31), type=c("yelkouan","yelkouan","yelkouan","yelkouan","yelkouan","yelkouan","Menorca","Menorca","Menorca","Menorca","Menorca","Menorca","Cabrera","Cabrera","Cabrera", "Pitiüses","Pitiüses","Pitiüses","Pitiüses","Pitiüses", "Mallorca", "Mallorca", "Mallorca", "Mallorca", "Mallorca", "Mallorca", "Mallorca", "Mallorca", "Mallorca", "Mallorca", "Mallorca", "Mallorca"))

#ggtree(tree) + geom_tiplab(size=3) + geom_nodepoint(aes(col=label), size = 2.5) + geom_hilight(data=loc, aes(node=node, fill=type), type = "roundrect")

# 4 You can try to apply node dating; root set at 2.82 (or 2.04-3.59) and P.puffinus split ser at 0.05 (0.03-0.07 Ma); then use a vector (written manually in a .csv file) where you specify the tip dates and the root date (while the rest are set to NA) - in our case, you can choose between Mean, Max and Min.

dates<-read.csv("datation.csv")
dates_def<-dates[["Mean"]]
mutation_rate <- estimate.mu(tree, dates_def)
node.dates.mean <- estimate.dates(tree, dates_def, mutation_rate, nsteps = 100)

dates<-read.csv("datation.csv")
dates_def<-dates[["Max"]]
mutation_rate <- estimate.mu(tree, dates_def)
node.dates.max <- estimate.dates(tree, dates_def, mutation_rate, nsteps = 100)

dates<-read.csv("datation.csv")
dates_def<-dates[["Min"]]
mutation_rate <- estimate.mu(tree, dates_def)
node.dates.min <- estimate.dates(tree, dates_def, mutation_rate, nsteps = 100)

dates_df<-data.frame(node.dates.mean, node.dates.max, node.dates.min, stringsAsFactors=FALSE)

write.csv(dates_df, "./node_dates.csv", row.names=TRUE)

# However, in our case, datation differs to much whether we use one calibration point or 2