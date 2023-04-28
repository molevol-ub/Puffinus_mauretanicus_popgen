library(tidyverse)
library(ggplot2)

setwd("/home/guillem/Documentos/TFM/chromopainter")
chunklengths <- as.data.frame(read.table("std.chunklengths.csv", header=T, sep=","))
chunklengths <- chunklengths %>% remove_rownames %>% column_to_rownames(var="Recipient")        # Anomenar columnes
order_inds <- c("PORQ","CZA11","ALT78","ILA13","ILA2","TZE1","M8","G14","G15","G4","M16","M6","G9","G10","G11","M11","G12","M17","M10","M21","M19","M12","M13","M14","M2","M3","M5","M20","M18","M4","G3","M1")
chunklengths <- chunklengths[order_inds, order_inds]

# Script from heatmap in struct-f4

suppressMessages(library(Rcpp))
suppressMessages(library(RcppParallel))
suppressMessages(library(RcppZiggurat))
suppressMessages(library(data.table))
suppressMessages(library(ape))
suppressMessages(library(optparse))
suppressMessages(library(this.path))
suppressMessages(library(R.utils))
suppressMessages(library(gplots))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))

# Rowv & Colv set to FALSE if you want to mantain the order of individuals set on your file

plotDRIFT<-function(drift){heatmap.2(drift,dendrogram="row",Rowv=FALSE, Colv=FALSE, symm=TRUE,scale="none",sepwidth=c(0,0),trace="none")}

chunklengths<-data.matrix(chunklengths, rownames.force=NA)

pdf(file="chunklengths.pdf", width=1000, height=1000)
plotDRIFT(abs(chunklengths))
dev.off()
