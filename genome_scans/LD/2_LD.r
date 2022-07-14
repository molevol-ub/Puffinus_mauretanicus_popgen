#!/soft/R-4.1.1

# Rscript to created heatmap of the given scaffold using the LDheatmap library

# Dependencies: LDheatmap, grid

library(LDheatmap)

setwd("/your/dir")

results_file<-as.matrix(read.table("input.ld", sep="\t", h=F))
mapfile<-as.list(read.table("input.mapfile", h=F))

# We need to coerce the list into a numeric object

mapfile<-as.numeric(unlist(mapfile))

# Now you can plot: if you don't want to include the scale, don't specify genetic.distances

png("your_name.png", width=1500, height=1500, res=150)
LD_hmap<- LDheatmap(results_file, genetic.distances=mapfile, flip=TRUE, color=heat.colors(20), title = "LD heatmap")
dev.off()
