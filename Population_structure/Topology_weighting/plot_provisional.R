pdf("topology_weighting.def.pdf", width=9, height=8.85)
ggplot(infile, aes(x=Topology, y=Weights, fill=Topology)) + geom_violin(trim=FALSE) +geom_boxplot(width = 0.15, fill="white", outlier.shape=NA) + theme_classic() + scale_fill_brewer(palette="Dark2") + theme(legend.position="none") + ylab ("Fraction of the genome") + theme(axis.text = element_text(size = 15), axis.title=element_text(size=18))
dev.off()
