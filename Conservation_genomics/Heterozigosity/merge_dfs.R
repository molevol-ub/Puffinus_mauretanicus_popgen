#!/soft/R-4.1.1/bin/Rscript

setwd("/users-d3/jferrer/gizquierdo/TFM/conservation_genomics/het/PopGenome")

# Script to merge the windowed estimates of heterozigosity for all indvs and correct them

# First use the specified arguments (individual name)

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=1) {
  stop("An ind name must be specified")
} 


#-------------------------------------------------------------------------------------------------------------------------------

# Open files

def <- read.table("het_popgenome.def.csv", sep="\t", header=TRUE)
M12 <- read.table(paste0(args[1],".windowed.pi"), sep="\t", header=TRUE)

# Add +1 to the start of the PopGenome Table so that its ID matches with the vcftools one

M12$BIN_START <- M12$BIN_START -1

# Define columns to merge to get the "position" column in each dataframe

cols_1 <- c("scf_name","start", "stop")
cols_2 <- c("CHROM","BIN_START","BIN_END")

# Create "position" columns for each dataframe

def$position <- apply( def[, cols_1] , 1 , paste , collapse = "-" )
M12$position <- apply( M12[, cols_2] , 1 , paste , collapse = "-" )

#-------------------------------------------------------------------------------------------------------------------------------

# Open file and create position column

maskfile <- read.table(paste0(args[1],".mask_overlap.bed"), sep="\t", header=FALSE)
cols_3 <- c("V1", "V2", "V3")
maskfile$position <- apply( maskfile[, cols_3 ] , 1 , paste , collapse = "-" )

# Use aggregate to know the sum for position

mask_total <- aggregate(V7 ~ position,
	data=maskfile,
	function(x) sum(x)/50000)                # 50k instead of 25k because for some reason every entry is repeated in the mask file


#-------------------------------------------------------------------------------------------------------------------------------

# Combine the 3 dataframes

all <- merge(merge(def, M12, by="position", all=TRUE), mask_total, by="position", all=TRUE)

# Remove the columns where there's no information for the PopGenome file (as they represent deleted windows with too much repeat content)

all <- all[!is.na(all$scf_name),]

# Change all NA values for zeroes

all[is.na(all)] <- 0

# Remove the columns you are not interested in, as well as the V7 column (unmasked proportion)

all <- subset(all, select = -c(position,CHROM,BIN_START,BIN_END,N_VARIANTS))
colnames(all)[which(names(all) == "V7")] <- "Unmasked_proportion"


#-------------------------------------------------------------------------------------------------------------------------------

# Correct pi for the masked window proportion and remove the "Unmasked proportion" column

all$PI <- all$PI/all$Unmasked_proportion
all <- subset(all, select = -c(Unmasked_proportion))

colnames(all)[which(names(all) == "PI")] <- as.character(args[1])      # and change PI column name,

# Save with the PopGenome name

write.table(x=all, file="het_popgenome.def.csv", sep="\t", quote=FALSE, row.names=FALSE)
