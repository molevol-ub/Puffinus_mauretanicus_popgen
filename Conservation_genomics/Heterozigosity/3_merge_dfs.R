# Open files

def <- read.table("het_popgenome.def.csv", sep="\t", header=TRUE)
TZE1 <- read.table("TZE1.windowed.pi", sep="\t", header=TRUE)

# Add +1 to the start of the PopGenome Table so that its ID matches with the vcftools one

#def$start <- def$start + 1 # Only the first time!!!

# Define columns to merge to get the "position" column in each dataframe

cols_1 <- c("scf_name","start", "stop")
cols_2 <- c("CHROM","BIN_START","BIN_END")

# Create "position" columns for each dataframe

def$position <- apply( def[ , cols_1 ] , 1 , paste , collapse = "-" )
TZE1$position <- apply( TZE1[ , cols_2 ] , 1 , paste , collapse = "-" )

# Combine the 2 dataframes

all <- merge(def, TZE1, by="position", all=TRUE)

# Remove the columns where there's no information for the PopGenome file (as they represent deleted windows with too much repeat content)

all <- all[!is.na(all$scf_name),]

# Change all NA values for zeroes

all[is.na(all)] <- 0

# Remove the columns you are not interested in, and change PI column name

all <- subset(all, select = -c(position,CHROM,BIN_START,BIN_END,N_VARIANTS))
colnames(all)[which(names(all) == "PI")] <- "PI_TZE1"

# Correct pi for the masked window proportion

all$PI_TZE1 <- all$PI_TZE1/all$Unmasked_proportion

# Save with the PopGenome name

write.table(x=all, file="het_popgenome.def.csv", sep="\t", quote=FALSE, row.names=FALSE)