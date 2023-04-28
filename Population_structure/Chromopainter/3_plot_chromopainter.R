library(tidyverse)
library(ggplot2)

setwd("/home/guillem/Documentos/TFM/chromopainter")
chunklengths <- as.data.frame(read.table("std.chunklengths.csv", header=T, sep="\t"))
chunklengths <- chunklengths %>% remove_rownames %>% column_to_rownames(var="Recipient")        # Anomenar columnes
order_inds <- c("PORQ","CZA11","ALT78","ILA13","ILA2","TZE1","M8","G14","G15","G4","M16","M6","G9","G10","G11","M11","G12","M17","M10","M21","M19","M12","M13","M14","M2","M3","M5","M20","M18","M4","G3","M1")
chunklengths <- chunklengths[order_inds, order_inds]
