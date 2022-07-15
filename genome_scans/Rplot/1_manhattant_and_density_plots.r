#!/soft/R-4.1.1/bin/R

library(tidyverse)
library(ggplot2)
library(grid)

setwd("/home/guillem/Documentos/TFG/genome_scans/25_auto_auto/def")

results_file<-as.data.frame(read.table("../auto_selscan_100.25.csv", sep="\t", h=T))

results_file<-results_file[order(-results_file$scf_length),]
head(results_file)

results_file$WindowNumber <- 1:nrow(results_file)

# Look! Up to top 0.05%

Fst.thresh=as.numeric(quantile(results_file$Fst, probs=0.995,na.rm=T))
Fst.mean=as.numeric(mean(results_file$Fst, na.rm=TRUE))
Fst_outliers <- rbind(subset(results_file, results_file$Fst > Fst.thresh))
write.table(x = Fst_outliers, 
	file = paste("PopGen_Fst_outliers.tsv", sep=""),
	quote=FALSE,  col.names = T, row.names = FALSE, sep= "\t")

XP.thresh=as.numeric(quantile(results_file$avg_XP.EHH, probs=0.9975,na.rm=T))
XP.negative.thresh=as.numeric(quantile(results_file$avg_XP.EHH,probs=0.0025,na.rm=T)) # For Tajima's D
XP.mean=as.numeric(mean(results_file$avg_XP.EHH, na.rm=TRUE))

XP_outliers <- rbind(subset(results_file, results_file$avg_XP.EHH > XP.thresh))
XP_outliers2 <- rbind(subset(results_file, results_file$avg_XP.EHH < XP.negative.thresh))

manhplot <- ggplot(results_file, aes(x = WindowNumber, y = Fst, colour=avg_XP.EHH)) + 
	geom_point(alpha = 0.75, stat = "identity", size=0.65, color="#3c81ba") +
	#scale_color_distiller(palette="Blues") +
	geom_hline(yintercept=Fst.thresh,linetype="dashed", size=0.3, colour="dimgrey") + 
	#geom_hline(yintercept=XP.negative.thresh,linetype="dashed", size=0.3, colour="dimgrey") + 
	geom_hline(yintercept=Fst.mean,linetype="dashed", size=0.2, colour="dimgrey") + 
	geom_point (data=XP_outliers, aes(x = WindowNumber, y = Fst), color="red", size=0.9) +
	geom_point (data=XP_outliers2, aes(x = WindowNumber, y = Fst), color="red", size=0.9) +
	theme_classic() +
	
#	annotate("rect",xmin=33490, xmax=33700, ymin=-Inf, ymax=Inf, fill="red", alpha=0.1) +
#	annotate("rect",xmin=28450, xmax=28700, ymin=-Inf, ymax=Inf, fill="red", alpha=0.1) +
#	annotate("rect",xmin=24200, xmax=24450, ymin=-Inf, ymax=Inf, fill="red", alpha=0.1) +
#	annotate("rect",xmin=23000, xmax=23250, ymin=-Inf, ymax=Inf, fill="red", alpha=0.1) +
#	annotate("rect",xmin=9350, xmax=9600, ymin=-Inf, ymax=Inf, fill="red", alpha=0.1) +
	
	theme (axis.text=element_text(size=12.5), axis.title=element_text(size=16)) + labs(x="Window Number", y="Fst")
	
	# + TITLE # LEGEND # X AXIS SCAFFOLD NAMES BELOW
	

png(paste("manhplot_Fst.25.png", sep=""), width=2000, height=600, res=150)
print(manhplot)
dev.off()


#----------------------------------------------------------------------------------------------------------------------------------------------------------------

# Also print a density plot (with Tajima's D & pi use both pops per graph)

Fst_density<-ggplot(results_file, aes(x=results_file$Fst)) + geom_density() + scale_x_continuous(name="Fst")
png("density_Fst_non0.png", width=1000, height=500, res=150)
Fst_density
dev.off()

dxy_density<-ggplot(results_file, aes(x=results_file$dxy)) + geom_density()  + scale_x_continuous(name="dxy")
png("density_dxy.png", width=1000, height=500, res=150)
dxy_density
dev.off()

pi_density<-ggplot(results_file) + geom_density(aes (x=results_file$pi_yelkouan, col="P. yelkouan")) + geom_density(aes (x=results_file$pi_mauretanicus, col="P. mauretanicus")) + scale_x_continuous(name="pi")
png("density_pi.png", width=1000, height=500, res=150)
pi_density
dev.off()

Tajima_density<-ggplot(results_file) + geom_density(aes (x=results_file$TajimaD_yelkouan, col="P. yelkouan")) + geom_density(aes (x=results_file$TajimaD_mauretanicus, col="P. mauretanicus"))  + scale_x_continuous(name="Tajima's D")
png("density_TajimaD.png", width=1000, height=500, res=150)
Tajima_density
dev.off()

deltapi_density<-ggplot(results_file, aes(x=results_file$deltapi)) + geom_density()  + scale_x_continuous(name="deltapi")
png("density_deltapi.png", width=1000, height=500, res=150)
deltapi_density
dev.off()


#---------------------------------------------------------------------------------------------------------------------------------------------------------------------

# I fer Manhattan plots comparatius amb els outliers de pi i Tajima's D (and we add 2 thresholds: negative and positive for the oppoiste and current population respectively)

outlier_file<-as.data.frame(read.table("PopGen_Fst_outliers.tsv", h=T))
main_file<-as.data.frame(read.table("Z_statistics.gc.unmasked.csv", sep="\t", h=T))

upper.thresh=as.numeric(quantile(main_file$TajimaD_yelkouan,probs=0.999,na.rm=T))
#lower.thresh=as.numeric(quantile(main_file$TajimaD_yelkouan,probs=0.001,na.rm=T))

upper.m.thresh=as.numeric(quantile(main_file$TajimaD_mauretanicus,probs=0.999,na.rm=T))
#lower.m.thresh=as.numeric(quantile(main_file$TajimaD_mauretanicus,probs=0.001,na.rm=T))

outlierplot <- ggplot(outlier_file) + geom_point(aes(x = WindowNumber, y = TajimaD_yelkouan, col="P. yelkouan")) +
	geom_point(aes(x = WindowNumber, y = TajimaD_mauretanicus, col="P. mauretanicus")) +
	geom_hline(yintercept=upper.thresh,linetype="dashed", size=0.5, colour="blue") +
	geom_hline(yintercept=upper.m.thresh,linetype="dashed", size=0.5, colour="red") + scale_y_continuous(name="Tajima's D") +
	ggtitle ("Fst outliers") #+ xlim (0, 250) +
	#geom_hline(yintercept=lower.m.thresh,linetype="dashed", size=0.5, colour="red") +  geom_hline(yintercept=lower.thresh,linetype="dashed", size=0.5, colour="blue")

png("outliers_Fst-TajimaD.25.png", width=2000, height=750, res=150)
outlierplot
dev.off()

# Also for dxy

outlier_file<-as.data.frame(read.table("PopGen_Fst_outliers.tsv", h=T))
main_file<-as.data.frame(read.table("Z_statistics.gc.unmasked.csv", sep="\t", h=T))

upper.thresh=as.numeric(quantile(main_file$deltapi,probs=0.99,na.rm=T))
results_mean=as.numeric(mean(main_file$deltapi, na.rm=TRUE))
lower.thresh=as.numeric(quantile(main_file$deltapi,probs=0.01,na.rm=T))

outlierplot <- ggplot(outlier_file) + geom_point(aes(x = WindowNumber, y = deltapi)) +
	geom_hline(yintercept=upper.thresh,linetype="dashed", size=0.5, colour="red") +
	geom_hline(yintercept=results_mean,linetype="dashed", size=0.25, colour="red") +
	geom_hline(yintercept=lower.thresh,linetype="dashed", size=0.5, colour="red") + scale_y_continuous(name="deltapi") +
	ggtitle ("Fst outliers") #+ xlim (0,1700)

png("outlier
