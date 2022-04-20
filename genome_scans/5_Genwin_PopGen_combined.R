#!/soft/R-4.1.1/bin/Rscript

library(PopGenome)

setwd("/users-d3/jferrer/gizquierdo/TFG/genome_scans/combined")

#Using the windows from the splineAnalyze file (should we add the faulty scaffolds that weren't included there?)

window_list <- as.data.frame(read.table("../splineAnalyze/Genwin_windows_output.csv", sep="\t", h=T))

# Iterate through list to define the regions of the 

full_F_ST<-data.frame()
full_PI<-data.frame()
full_dxy<-data.frame()
full_Tajima_D<-data.frame()

for (n in 1:(length(window_list$WindowStart))){
	positions<-list(c(window_list$WindowStart[n]:window_list$WindowStop[n]))
	tryCatch({vcf_file<-readVCF("../vcfs/Puffinus_SNP.maxmiss80.filtered.noPP.merged.nomono.vcf.gz",10000, window_list$scaffold[n], window_list$WindowStart[n], window_list$WindowStop[n])
	vcf_file <- set.populations (vcf_file, list (c("ALT78", "CZA11", "ILA13", "ILA2", "PORQ", "TZE1", "G4", "M8"), c("M1", "M11", "M12", "M13", "M18", "M2", "M21", "M4")))

	vcf_file<-F_ST.stats(vcf_file)
				
	F_ST<-as.data.frame(vcf_file@nucleotide.F_ST)
	full_F_ST <- data.frame(rbind(full_F_ST, data.frame(F_ST)))
	
	PI<-as.data.frame(vcf_file@Pi)
	full_PI <- data.frame(rbind(full_PI, data.frame(PI)))

	vcf_file<-neutrality.stats(vcf_file)
	
	Tajima_D<-as.data.frame(vcf_file@Tajima.D)
	full_Tajima_D <- data.frame(rbind(full_Tajima_D, data.frame(Tajima_D)))		

	vcf_file<-diversity.stats.between(vcf_file)

	dxy<-as.data.frame(vcf_file@nuc.diversity.between)
	dxy$scf<- window_list$scaffold[n]
	dxy$start <- window_list$WindowStart[n]
	dxy$stop <- window_list$WindowStop[n]

	full_dxy<-data.frame(rbind(full_dxy, data.frame(dxy)))},error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

}

all_data <- data.frame(full_F_ST$FST..Nucleotide., full_PI$pop.1, full_PI$pop.2, full_Tajima_D$pop.1, full_Tajima_D$pop.2, full_dxy$pop1.pop2, full_dxy$scf, full_dxy$start, full_dxy$stop)

write.table(all_data, file="Genwin_PopGen_combined.csv", sep="\t", row.names=FALSE)
