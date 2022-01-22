################################################
# CH1 telomere length vs genomic vulnerability #
################################################
library(dplyr)
library(ggplot2)
library(usmap)
library(cowplot)
library(sp)
library(rgdal)
library(viridis)
library(maps)
library(RColorBrewer)
library(mapdata)
library(plotly)
library(ggspatial)
library(sf)
library(rgeos)
library(ggmap)
library(devtools)



setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Research/PHD_Bio_CSU/CH1_YEWA_telomeres/Analysis")

## Mapping points and TL ####

#read a file containing YEWA TL and metdata
YEWATL_csv=read.csv("Data/YEWA_TL_data/YEWA_TL_GV_NoHY_2021.csv")

class(YEWATL_csv$TS_Mean)
class(YEWATL_csv$Lat)
class(YEWATL_csv$Long)


#Explore some simple statistics 
ggplot(data=YEWATL_csv)+
  geom_histogram(mapping=aes(TS_Mean), bins = 20, fill="darkgreen")+
  theme_classic()+
  labs(y = "Frequency", x = "Average Normalized TL")+
  theme(text = element_text(size=15))


range(YEWATL_csv$TS_Mean)
mean(YEWATL_csv$TS_Mean)


#ggplot has some map data too
states_map <- map_data("state")

#use ggplot to map a state map
TL_map=ggplot() +
  geom_polygon(data = states_map, aes(x=long, y = lat, group = group), fill = NA, color ="black", lwd=.1) +
  geom_point(data = YEWATL_csv, position=position_jitter(h=0.4,w=0.4), aes(x=Long, y = Lat, colour=TS_Mean), alpha =.6) +
  coord_map("bonne", lat0 = 40)+
  theme_bw()+
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank())+
  theme( panel.border=element_blank())+
  scale_color_viridis(breaks=seq(0, 6, by=1), name = "Average Normalized TL")+
  guides(color= guide_legend(), size=guide_legend())+
  labs(y = "Latitude", x = "Longitude")
TL_map

ggplotly(TL_map)


#### Mapping genomic vulnerability ######

#read a file containing Rachael Bay's genomic vulnerabilty estimates
YEWA_GV_csv=read.csv("EuclideanGenomicVulnerability_02.09.17.csv")



class(YEWATL_csv$TF85_1)
class(YEWATL_csv$Lat)
class(YEWATL_csv$Long)


ggplot(data=YEWATL_csv)+
  geom_histogram(mapping=aes(TF85_1), bins = 30, fill="darkgreen")+
  theme_classic()+
  labs(y = "Frequency", x = "GV TF85")+
  theme(text = element_text(size=15))


#use ggplot to map a state map
YEWA_GV=ggplot() +
  geom_polygon(data = states_map, aes(x=long, y = lat, group = group), fill = "white", color ="black", lwd=.1) +
  geom_point(data = YEWATL_csv, aes(x=Long, y = Lat, colour=TF85_1), size=1, alpha = 20) +
  coord_map()+
  scale_colour_gradientn(colours = c("deepskyblue3", "yellow", "orange","orange","red", "red", "red", "red", "red"),breaks = c(0, 0.001,0.003, 0.004,0.005, 0.006, 0.01, 0.085), limits=c(0,0.085), name = "Genomic Vulnerability")+
  labs(y = "Latitude", x = "Longitude")
YEWA_GV


#read in breeding population shapefile
YEWA_pops=st_read("YWAR_Pops_shapefile/E.shp")
YEWA_pops_df=as.data.frame(YEWA_pops)
plot(YEWA_pops)

#create dataframe with population numbers and names
ID <- c("0", "1", "2", "3", "4")
YWARpops_names <- c("Arctic", "Central","Western", "Eastern", "Northwest" )
YWARpops_df <- data.frame(ID, YWARpops_names, row.names = ID)

#read in population shapefile
YWARpops <- readOGR("YWAR_Pops_shapefile/E.shp")
summary(YWARpops)

#change projection
YWARpops1 <- spTransform(YWARpops, CRS("+init=epsg:4326"))

#create spatial polygons data frame by combing shapefile and data frame
YWARpops2 <- SpatialPolygonsDataFrame(YWARpops1, YWARpops_df, match.ID = FALSE)
plot(YWARpops2)


#### Map on ggogle maps

YEWA_GV_map=get_googlemap(center= c(lon=-97, lat = 40), zoom = 4, maptype = "hybrid") %>%
  ggmap() +
  geom_point(data = YEWATL_csv, aes(x=Long, y = Lat, colour=TF85_1), size=0.8, alpha = 1) +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank())+
  theme(panel.border=element_blank())+
  scale_colour_gradientn(colours = c("deepskyblue3", "yellow", "red","red","red", "red", "red", "red", "red"),breaks = c(0, 0.001,0.003, 0.004,0.005, 0.006, 0.01, 0.085), limits=c(0,0.085), name = "Genomic Vulnerability")+
  labs(y = "Latitude", x = "Longitude")

YEWA_GV_map



