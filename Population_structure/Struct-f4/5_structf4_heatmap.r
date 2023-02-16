#!/soft/R-4.1.1

# R script to plot the results of generate heatmap with the modified individual-order; required packages - see below

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

setwd ("/path/to/your/directory")
tab_file<-read.table("your_file.csv")

# Rowv & Colv set to FALSE if you want to mantain the order of individuals set on your file

plotDRIFT<-function(drift){heatmap.2(drift,dendrogram="row",Rowv=FALSE, Colv=FALSE, symm=TRUE,scale="none",sepwidth=c(0,0),trace="none")}

table_def<-data.matrix(tab_file, rownames.force=NA)

png(file="./your_output.png", width=1000, height=1000)
plotDRIFT(abs(table_def))
dev.off()
