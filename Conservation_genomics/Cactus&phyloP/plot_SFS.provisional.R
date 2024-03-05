library(ggplot2)
infile_plot<-read.table("Puffinus.SFS.qval.plot.csv", header=T)
ggplot(infile_plot, aes(x=SFS, weight=FREQ, color=QVAL_CAT, fill=QVAL_CAT)) + geom_bar(position="dodge") + theme_classic() + scale_fill_brewer(palette="YlOrRd") + scale_color_brewer(palette="YlOrRd") + scale_x_continuous(breaks=seq(1,17,1))

ggsave("SFS_per_qval_cat.NUEDEL.pdf", width = 10, height = 6)

infile_plot[infile_plot$QVAL_CAT=="HIGHLY DELETERIOUS",]$QVAL_CAT<-"HIGH phyloP"
infile_plot[infile_plot$QVAL_CAT=="HIGH snpEff_phyloP",]$QVAL_CAT<-"A. HIGH snpEff_phyloP"
infile_plot[infile_plot$QVAL_CAT=="HIGH phyloP",]$QVAL_CAT<-"B. HIGH phyloP"
infile_plot[infile_plot$QVAL_CAT=="HIGH snpEff",]$QVAL_CAT<-"C. HIGH snpEff"
infile_plot[infile_plot$QVAL_CAT=="NEUTRAL",]$QVAL_CAT<-"D. NEUTRAL"
ggplot(infile_plot[(infile_plot$QVAL_CAT=="C. HIGH snpEff")|(infile_plot$QVAL_CAT=="A. HIGH snpEff_phyloP")|(infile_plot$QVAL_CAT=="B. HIGH phyloP")|(infile_plot$QVAL_CAT=="D. NEUTRAL"),], aes(x=SFS, weight=FREQ, color=QVAL_CAT, fill=QVAL_CAT)) + geom_bar(position="dodge") + theme_classic() + scale_fill_brewer(palette="YlOrRd") + scale_color_brewer(palette="YlOrRd") + scale_x_continuous(breaks=seq(1,17,1)) + ylab("FREQ")

ggsave("SFS_per_RISK.pdf", width = 10, height = 6)
