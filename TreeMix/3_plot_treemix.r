# R script to plot the results of the treemix analysis

library(RColorBrewer)
library(R.utils)

setwd("/home/guillem/Documents/TFG/treemix/bootstrap_1000_output")
prefix="pmau_treemix"

#Also set path and source the required script

source("../plotting_funcs.R")

# Now you can run the script for every "edge" you want (

par(mfrow=c(2,3))

for(edge in 0:10){
	plot_tree(cex=0.8,paste0(prefix,".",edge))
	title(paste(edge,"edges"))
	file_name = paste("pmau_treemix_plot",edge,".png", sep="")
	png(file_name, width=1500, height=700, res=150)
}

# Remember: plot_tree(cex=0.8,paste0(prefix,".",edge,".uncorrected")) if using the uncorrected files

# This last bit is to check which parts of the tree are not well modelled for the different runs

for(edge in 0:10){
	plot_resid(stem=paste0(prefix,".",edge),pop_order="../pmau_pops.list")
	file_name = paste("pmau_residuals_plot",edge,".png", sep="")
	png(file_name, width=1500, height=700, res=150)
}
