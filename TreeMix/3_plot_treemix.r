#!/soft/R-4.1.1

# R script to plot the results of the TreeMix analysis, as well as to tun OptM
# Required packages are: RColorBrewer & R.utils
# It also requires a script ("plotting_funcs.R") whose path you should specify

library(RColorBrewer)
library(R.utils)
source("../plotting_funcs.R")

setwd("/your/dir")
prefix="your_file_prefix"

# 1. Now you can run the script to plot for every "edge" you want; you also check the % variance explained (csv output)

for(edge in 0:10){
	file_name = paste("pmau_treemix.",edge,".main.png", sep="")
	png(file_name, width=1300, height=650, res=150)
	plot_tree(cex=1,paste0(prefix,".",edge, ".main"), mbar = F, arrow=0.1, ybar=0.55, plus=0.03)
	title(paste(edge,"edges"))
	dev.off()

	var_val<-get_f(stem=paste0(prefix,".",edge, ".main"))
	output<-as.data.frame(var_val)
	table_name=paste("explained_variance.",edge,".csv")
	write.table(output, table_name)
	
}

# 2. This bit is to check which parts of the tree are not well modelled for the different runs (plot the residuals)

for(edge in 0:10){

	file_name = paste("pmau_residuals.main.",edge,".png", sep="")
	png(file_name, width=1500, height=1100, res=150)
	plot_resid(stem=paste0(prefix,".",edge, ".main"),pop_order="../pmau_pops.list")
	dev.off()
}


#-----------------------------------------------------------------------------------------------------------

# You can use OptM to check which TreeMix graph is most viable using the 3 algorithms proposed ("Evanno", "linear" and "SiZer")

setwd("/home/guillem/Documents/TFG/treemix_v2")

library("OptM")

# 3.1. "Evanno"-like method

evanno.optM=optM("/home/guillem/Documents/TFG/treemix_v2/subgroups_output_3/OptM")

plot_optM(evanno.optM, method = "Evanno", pdf="./subgroups_output_3/Evanno-OptM.pdf")
dev.off()

# 3.2. Linear modeling estimates method

linear.optM=optM("/home/guillem/Documents/TFG/treemix_v2/subgroups_output_3/OptM", method = "linear")

file_name = paste("./subgroups_output_3/Linear-OptM.png", sep="")
png(file_name, width=1500, height=950, res=150)
plot_optM(linear.optM, method = "linear")
dev.off()

# 3.3. SiZer method

SiZer.optM=optM("/home/guillem/Documents/TFG/treemix_v2/subgroups_output_3/OptM", method = "SiZer")

file_name = paste("./subgroups_output_3/SiZer-OptM.png", sep="")
png(file_name, width=1500, height=950, res=150)
plot_optM(SiZer.optM, method = "SiZer")
dev.off()
