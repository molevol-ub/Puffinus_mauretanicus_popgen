#!/soft/R-4.1.1/bin/R
#$ -cwd
#$ -R y
#$ -e popgenome.err
#$ -o popgenome.out
#$ -q h12.q
#$ -pe ompi128h12 2
#$ -V                    #export environment var
#$ -N PopGenome
#$ -M 000izquierdoguillem@gmail.com
#$ -m be

# Manual found in: https://github.com/tonig-evo/workshop-popgenome/blob/tut2019/Whole_genome_analyses_using_VCF_files.pdf

# Previously get a list of unique scaffolds with "zcat Puffinus_subset.recode.vcf.gz | grep -v "^#" | cut -f1 | sort | uniq > scaff_names" and make table with its first and last positions from scf_info_def.txt
library(PopGenome)
library(tidyverse)
library(ggplot2)

setwd("/users-d3/jferrer/gizquierdo/TFM/genome_scans/PopGenome")

#First we take out the scaffold names and lengths

scf_list <- read.table("../../vcfs/scf_length.csv", sep="\t", h=F)

# Now we iterate through all of them and add them to a common dataframe; trycatch is used to detect and skip those scaffolds with length inferior to the length of the windows

# REMEMBER TO USE "include.unknown=TRUE" TO HANDLE MISSING DATA!

# First we calculate Fst between populations (nucleotide) & pi; then dxy (where, as it is the last one, we add the scf_ID & position informations)

full_F_ST<-data.frame()
full_PI<-data.frame()
full_dxy<-data.frame()
full_Tajima_D<-data.frame()
full_SNP<-data.frame()

for (n in 1:(length(scf_list$V1))){

	message(scf_list$V1[n])

	vcf_file<-readVCF("../../vcfs/Puffinus_SNP.maxmiss80.filtered.noPP.merged.nomono.masked.vcf.gz",10000, scf_list$V1[n], 1, scf_list$V2[n], include.unknown=TRUE)
	
	vcf_file <- set.populations (vcf_file, list (c("ALT78", "CZA11", "ILA13", "ILA2", "PORQ", "TZE1", "G4", "M8"), c("M1", "M11", "M14", "M13", "M18", "M2", "M21", "M4")))
	
	tryCatch({vcf_file.slide <- sliding.window.transform(vcf_file,width=25000, jump=25000, type=2, whole.data=TRUE)
		
		vcf_file.slide<-F_ST.stats(vcf_file.slide)
				
		F_ST<-as.data.frame(vcf_file.slide@nucleotide.F_ST)
		full_F_ST <- data.frame(rbind(full_F_ST, data.frame(F_ST)))
		
		PI<-as.data.frame(vcf_file.slide@Pi)
		full_PI <- data.frame(rbind(full_PI, data.frame(PI)))

		vcf_file.slide<-neutrality.stats(vcf_file.slide)
		
		Tajima_D<-as.data.frame(vcf_file.slide@Tajima.D)
		full_Tajima_D <- data.frame(rbind(full_Tajima_D, data.frame(Tajima_D)))		

		vcf_file.slide<-diversity.stats.between(vcf_file.slide)

		dxy<-as.data.frame(vcf_file.slide@nuc.diversity.between)
		dxy$scf<-scf_list$V1[n]

		full_dxy<-data.frame(rbind(full_dxy, data.frame(dxy)))

		SNP<-lengths(vcf_file.slide@region.data@biallelic.sites)
		full_SNP<-data.frame(rbind(full_SNP, data.frame(SNP)))},error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

}

all_data <- data.frame(full_F_ST$FST..Nucleotide., full_PI$pop.1, full_PI$pop.2, full_Tajima_D$pop.1, full_Tajima_D$pop.2, full_dxy$pop1.pop2, full_SNP$SNP, full_dxy$scf)

write.table(all_data, file="Fst_paiwise.csv", sep="\t", row.names=FALSE)

#------------------------------------------------------------------------------------------------------------------------------------------------------------------

stop() #If you want to stop after the first part is done

# Now create the Manhattan plot (first adding the window number and window length) - previously turn negative Fsts to 0
# If using autosomes/chrZ, use: V1-Fst; V2&3 - Pi; V4&5 - Tajima's D; V6 - dxy; V7 - scf_names

results_file<-as.data.frame(read.table("autosome_statistics.unmasked.csv", sep="\t", h=T))

results_file<-results_file[order(-results_file$scf_length),]
head(results_file)

results_file$WindowNumber <- 1:nrow(results_file)

for (var in colnames(results_file[, c(1,2,3,4,5,6,7)])){

results_file.thresh=as.numeric(quantile(results_file[[var]], probs=0.999,na.rm=T))
#results_file.negative.thresh=as.numeric(quantile(results_file$V4,probs=0.001,na.rm=T)) # For Tajima's D
results_mean=as.numeric(mean(results_file[[var]], na.rm=TRUE))

dxy_outliers <- subset(results_file, Fst > results_file.thresh)

write.table(x = dxy_outliers, 
	file = paste("PopGen_", var, "_outliers.tsv", sep=""),
	quote=FALSE,  col.names = T, row.names = FALSE, sep= "\t")

manhplot <- ggplot(results_file, aes(x = WindowNumber, y = results_file[[var]])) +
	geom_point(alpha = 0.75, stat = "identity", size=0.8) +
	scale_color_manual(values= c("Black","darkgray")) + # TITLE # LEGEND # X AXIS SCAFFOLD NAMES BELOW
	geom_hline(yintercept=results_file.thresh,linetype="dashed", size=0.5, colour="red") + 
	#geom_hline(yintercept=results_file.negative.thresh,linetype="dashed", size=0.5, colour="red") + 
	geom_hline(yintercept=results_mean,linetype="dashed", size=0.3, colour="red") + 
	scale_y_continuous(name=var)

png(paste("manhplot_", var, ".10.png", sep=""), width=2000, height=750, res=150)
print(manhplot)
dev.off()

}

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

# Headers: full_F_ST.FST..Nucleotide. ; full_dxy.pop1.pop2 ; full_PI.pop.1 (yelkouan) ; full_PI.pop.2 (mauretanicus), full_Tajima_D.pop.1 (or 2)

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------

# I fer Manhattan plots comparatius amb els outliers de pi i Tajima's D (and we add 2 thresholds: negative and positive for the oppoiste and current population respectively)

outlier_file<-as.data.frame(read.table("PopGen_Fst_outliers.tsv", h=T))
main_file<-as.data.frame(read.table("autosome_statistics.unmasked.csv", sep="\t", h=T))

upper.thresh=as.numeric(quantile(main_file$pi_yelkouan,probs=0.999,na.rm=T))
#lower.thresh=as.numeric(quantile(main_file$TajimaD_yelkouan,probs=0.001,na.rm=T))

upper.m.thresh=as.numeric(quantile(main_file$pi_mauretanicus,probs=0.999,na.rm=T))
#lower.m.thresh=as.numeric(quantile(main_file$TajimaD_mauretanicus,probs=0.001,na.rm=T))

outlierplot <- ggplot(outlier_file) + geom_point(aes(x = WindowNumber, y = pi_yelkouan, col="P. yelkouan")) + geom_point(aes(x = WindowNumber, y = pi_mauretanicus, col="P. mauretanicus")) + geom_hline(yintercept=upper.thresh,linetype="dashed", size=0.5, colour="blue") + geom_hline(yintercept=upper.m.thresh,linetype="dashed", size=0.5, colour="red") + scale_y_continuous(name="Tajima's D") + ggtitle ("Fst outliers") + #geom_hline(yintercept=lower.m.thresh,linetype="dashed", size=0.5, colour="red") +  geom_hline(yintercept=lower.thresh,linetype="dashed", size=0.5, colour="blue")

png("outliers_Fst-pi.10.png", width=2000, height=750, res=150)
outlierplot
dev.off()

# Also for dxy

outlier_file<-as.data.frame(read.table("PopGen_Fst_outliers.tsv", h=T))
main_file<-as.data.frame(read.table("autosome_statistics.unmasked.csv", sep="\t", h=T))

upper.thresh=as.numeric(quantile(main_file$deltapi,probs=0.999,na.rm=T))

outlierplot <- ggplot(outlier_file) + geom_point(aes(x = WindowNumber, y = deltapi)) + geom_hline(yintercept=upper.thresh,linetype="dashed", size=0.5, colour="blue") + scale_y_continuous(name="dxy") + ggtitle ("Fst outliers")

png("outliers_Fst-deltapi.10.png", width=2000, height=750, res=150)
outlierplot
dev.off()
