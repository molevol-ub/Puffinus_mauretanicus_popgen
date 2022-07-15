#!/soft/R-4.1.1

# Rscript to generate 2 genotype plots: one with P.puffinus, one without; based on https://github.com/JimWhiting91/genotype_plot
# You might want to consider: a) Varying missingness to include the max # SNPs and b) Set a high mac filter - this will reduce low-freq variants, which might be better for haplotype visualization purposes, but you might loose information.

library(GenotypePlot)
library(vcfR)
library(ggplot2)

# 1. Preparing input files: aside from a population map with 2 columns, a VCF file of the desired region has to be compressed (bgzip) and indexed (tabix) before reading it

setwd("/your/directory")

# 2. Read in the subsetted VCF as well as the "population map"

our_popmap<-as.data.frame(read.table("../../Puffinus.list", h=F, sep=",", stringsAsFactors = FALSE))
our_vcf <- read.vcfR("your_file.vcf.gz")

# 3. Generate plot and save as png; you can only do 1 scaffold at a time!

# Options:	cluster --> whether or not to cluster haplotypes for similarity
#		snp_label_size --> interval between position labels
#		polarise_genotypes --> IMPORTANT! You may polarize for: a) The ancestral population (outgroup) - but this might be less efficient if that one is heterozygous - or b) Either of the most dissimilar populations (better for visualization purposes)

plot_scf1 <- genotype_plot(vcf_object=our_vcf, popmap=our_popmap, cluster =FALSE, snp_label_size=100000, 
	colour_scheme=c("#7A0001","#d56d5d", "#789FF2"), invariant_filter=TRUE, polarise_genotypes="Pitiuses")

pdf("scf2_map.def.pdf", width=10, height=5)
combine_genotype_plot(plot_scf1)
dev.off()

# 4. Try also clustering individuals

plot_scf1 <- genotype_plot(vcf_object=our_vcf, popmap=our_popmap, cluster =TRUE, snp_label_size=50000, 
	colour_scheme=c("#7A0001","#d56d5d", "#789FF2"), invariant_filter=TRUE, polarise_genotypes="Pitiuses")

dendro_tips<-plot_scf1$dendrogram + geom_text(aes(x=1:length(plot_scf1$dendro_labels), y=-0.7, label=plot_scf1$dendro_labels), size=2.4)

plot_scf1$dendrogram<-dendro_tips

pdf("scf12.def.pdf", width=14, height=7)
combine_genotype_plot(plot_scf1)
dev.off()

# 5. Finally, also do allele frequencies (plot_allele_frequency=TRUE)

plot_scf1 <- genotype_plot(vcf_object=our_vcf, popmap=our_popmap, cluster =FALSE, snp_label_size=50000, 
	colour_scheme=c("#7A0001","#d56d5d", "#789FF2"), invariant_filter=TRUE, plot_allele_frequency=TRUE, polarise_genotypes="Pitiuses")

pdf("scf12.alfreq.def.pdf", width=14, height=7)
combine_genotype_plot(plot_scf1)
dev.off()
