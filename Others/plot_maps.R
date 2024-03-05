library(rworldmap)
library(rworldxtra)
library(sp)
library(rgeos)
library(sf)
library(ggplot2)
library("rnaturalearth")
library("rnaturalearthdata")

setwd("/home/guillem/Documentos/TFG/maps")
df <- read.delim("map.csv", header=TRUE)
coords <- df[,c(4,5)]
sp <- SpatialPoints(coords)
data <- df[,c(1,2,3)]
spdf = SpatialPointsDataFrame(coords, data)
spdf = SpatialPointsDataFrame(sp, data)

world <- ne_coastline(scale = "medium", returnclass = "sf")

#world <- ne_countries(scale="large", returnclass = "sf")


data_def<-data.frame(spdf)

ggplot(data=world) +
	theme(panel.grid.major=element_blank()) +
	geom_sf(col="gray50", fill="#e9e9e9") +
	#geom_point(data = df, aes(x = long, y = lat, col=colour), size = df$log*5, alpha=0.85) +
	#scale_color_manual(values=c("#C34632", "#8B0001", "#789FF2","#F26B38", "#0052A2", "#FF8532")) + #scale_shape_manual(values=c(17,16,18, 15)) +
	coord_sf(xlim = c(-21, 81), ylim = c(29, 56), expand = FALSE) +
	labs(x="", y="") + theme(axis.text=element_text(size=25))

ggsave("map_mad.3.pdf",plot=last_plot(),device="pdf", path=getwd(), height = 10 , width = 17)

# Global: xlim = c(-21, 46), ylim = c(29, 56)
# Mediterrani: xlim = c(-10, 20), ylim = c(32, 46)
# Balears: xlim = c(1, 5), ylim = c(38, 41)
