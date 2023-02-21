#!/soft/R-4.1.1/bin/Rscript
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

# Previously get a list of unique scaffolds with "zcat Puffinus_subset.recode.vcf.gz | grep -v "^#" | cut -f1 | sort | uniq > scaff_names" and make table with its first and last positions from scf_info_def.txt
library(PopGenome)

setwd("/users-d3/jferrer/gizquierdo/TFM/conservation_genomics/het/PopGenome")

# We reuse the scaffold ID & length table from the PopGenome scans of the outlier windows

scf_list <- read.table("/users-d3/jferrer/gizquierdo/TFM/genome_scans/vcfs/scf_length.csv", sep="\t", h=F)

# Now we iterate through all of them and add them to a common dataframe; trycatch is used to detect and skip those scaffolds with length inferior to 10kbp
# First we calculate Fst between populations (nucleotide) & pi; then dxy (where, as it is the last one, we add the scf_ID & position informations)

full_PI<-data.frame()

for (n in 1:(length(scf_list$V1))){

	message(scf_list$V1[n])

	vcf_file<-readVCF("/users-d3/jferrer/gizquierdo/TFG/chr_split/vcfs/auto/Puffinus_SNP.maxmiss100.filtered.noPP.merged.nomono.masked.vcf.gz",10000, scf_list$V1[n], 1, scf_list$V2[n], include.unknown=TRUE)
	
	vcf_file <- set.populations (vcf_file, list ("ALT78", "CZA11", "ILA13", "ILA2", "PORQ", "TZE1", "M8", "M1", "M5", "M20", "M4", "G3", "M2", "M3", "M18", "M13", "M19", "M14", "M12", "M11", "M17", "M10", "M21", "G12", "G9", "G10", "G11", "G14", "G15", "M16", "G4", "M6", "COP1", "LT2"))
	
	tryCatch({vcf_file.slide <- sliding.window.transform(vcf_file,width=25000, jump=25000, type=2, whole.data=TRUE)
		
		vcf_file.slide<-neutrality.stats(vcf_file.slide)
		
		PI<-as.data.frame(vcf_file.slide@theta_Watterson)
		PI$scf<-scf_list$V1[n]

		full_PI <- data.frame(rbind(full_PI, data.frame(PI)))},error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

}

write.table(full_PI, file="het_popgenome.raw.csv", sep="\t", row.names=FALSE)

#------------------------------------------------------------------------------------------------------------------------------------------------------------------