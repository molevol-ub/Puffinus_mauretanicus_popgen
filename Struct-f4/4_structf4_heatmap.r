# Scripts to repeat barplot and heatplot (using individuals)

# 1 Heatplot (repeat)

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

# Input files must consist only of the distance table included in the struct-f4 output between the lines "Initial K ...) and "lik". 

setwd ("/home/guillem/Documentos/programari/plibradosanz-structf4-ac34f3d2d2a4/minmac_100_k2_auto")
tab_file<-read.table("prova2.csv")

plotDRIFT<-function(drift){heatmap.2(drift,dendrogram="row",Rowv=FALSE, Colv=FALSE, symm=TRUE,scale="none",sepwidth=c(0,0),trace="none")}

table_def<-data.matrix(tab_file, rownames.force=NA)

png(file="./heatmap_pmau.100.minmac2.auto.k2_prova.png", width=1000, height=1000)
plotDRIFT(abs(table_def))
dev.off()

# 2 "Barplot" - modifying aesthetics
