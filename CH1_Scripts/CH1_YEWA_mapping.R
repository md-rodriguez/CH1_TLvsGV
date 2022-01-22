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
YEWATL_csv=read.csv("YEWA_TL_data.csv")

YEWATL_csv <- YEWATL_csv[-c(199,230:273), ]

class(YEWATL_csv$AveNormTASQ)

class(YEWATL_csv$Lat)
class(YEWATL_csv$Long)

YEWATL_csv$Lat <- as.numeric(as.character(YEWATL_csv$Lat))
YEWATL_csv$Long <- as.numeric(as.character(YEWATL_csv$Long))

class(YEWATL_csv$Lat)
class(YEWATL_csv$Long)

###Just messing with data ###
#subset to CO
YEWA_PA <- subset(YEWATL_csv, State=="Pennsylvania")
YEWA_MI <- subset(YEWATL_csv, State=="Michigan")

#subset to TL > 1
Long_YEWATL=subset(YEWATL_csv, AveNormTASQ>1)
head(Long_YEWATL)


#how many samples? 
nrow(YEWATL_csv)

#Explore some simple statistics 
ggplot(data=YEWATL_csv)+
  geom_histogram(mapping=aes(AveNormTASQ), bins = 20, fill="darkgreen")+
  theme_classic()+
  labs(y = "Frequency", x = "Average Normalized TL")+
  theme(text = element_text(size=15))


range(YEWATL_csv$AveNormTASQ)
mean(YEWATL_csv$AveNormTASQ)

#Let's make a shapefile of the locations
YEWATL_shp=YEWATL_csv
coordinates(YEWATL_shp)=~Long+Lat
proj4string(YEWATL_shp)= CRS("+proj=longlat +datum=WGS84")

YEWATL_csv=as.data.frame(YEWATL_csv)

#let's determine the range and save the valeus. We use these below for plotting
min_TL=min(YEWATL_csv$AveNormTASQ)
max_TL=max(YEWATL_csv$AveNormTASQ)

#Within R, there are some datasets we can manipulate and plot
states<-as.data.frame(state.x77)
states$region <- tolower(rownames(states))


dev.off()



#ggplot has some map data too
states_map <- map_data("state")

#use ggplot to map a state map
TL_map=ggplot() +
  geom_polygon(data = states_map, aes(x=long, y = lat, group = group), fill = NA, color ="black", lwd=.1) +
  geom_point(data = YEWATL_csv, position=position_jitter(h=0.4,w=0.4), aes(x=Long, y = Lat, colour=AveNormTASQ), alpha =.6) +
  coord_map("bonne", lat0 = 40)+
  theme_bw()+
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank())+
  theme( panel.border=element_blank())+
  scale_color_viridis(limits=c(min_TL, max_TL), breaks=seq(0, 6, by=1), name = "Average Normalized TL")+
  guides(color= guide_legend(), size=guide_legend())+
  labs(y = "Latitude", x = "Longitude")
TL_map

ggplotly(TL_map)


#### Mapping genomic vulnerability ######

#read a file containing Rachael Bay's genomic vulnerabilty estimates
YEWA_GV_csv=read.csv("EuclideanGenomicVulnerability_02.09.17.csv")



class(YEWA_GV_csv$tf5045)

class(YEWATL_csv$Lat)
class(YEWATL_csv$Long)

###Just messing with data ###
head(YEWA_GV_csv)

#how many samples? 
nrow(YEWA_GV_csv)

HighGV=subset(YEWA_GV_csv, tf5045>0.01)
nrow(HighGV)

ggplot(data=YEWA_GV_csv)+
  geom_histogram(mapping=aes(tf5045), bins = 30, fill="darkgreen")+
  theme_classic()+
  labs(y = "Frequency", x = "Average Normalized TL")+
  theme(text = element_text(size=15))


#Let's make a shapefile of the locations
YEWA_GV_shp=YEWA_GV_csv
coordinates(YEWA_GV_shp)=~Long+Lat
proj4string(YEWA_GV_shp)= CRS("+init=epsg:4326")



YEWA_GV_csv=as.data.frame(YEWA_GV_csv)

#let's determine the range and save the valeus. We use these below for plotting
min_GV_5045=min(YEWA_GV_csv$tf5045)
max_GV_5045=max(YEWA_GV_csv$tf5045)
Mean_GV_5045=mean(YEWA_GV_csv$tf5045)


#use ggplot to map a state map
YEWA_GV=ggplot() +
  geom_polygon(data = states_map, aes(x=long, y = lat, group = group), fill = NA, color ="black", lwd=.1) +
  geom_point(data = YEWA_GV_csv, aes(x=Long, y = Lat, colour=tf5045), size=0.01, alpha = .9) +
  coord_map()+
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank())+
  theme(panel.border=element_blank())+
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

#make sure CRS is the same for both shapefiles
st_crs(YWARpops2)==st_crs(YEWA_GV_shp)

out <- gIntersection(YEWA_GV_shp, YWARpops2, byid=TRUE)

YEWA_clipped <- as.data.frame(out)

names(YEWA_clipped)[names(YEWA_clipped) == "x"] <- "Long"
names(YEWA_clipped)[names(YEWA_clipped) == "y"] <- "Lat"

YEWA_GV2 <- merge(YEWA_clipped, YEWA_GV_csv, by = c("Lat","Long"),all.x=T)

saveRDS(YEWA_GV2, file="YEWA_GV2_map.Rda")
YEWA_GV2 <- readRDS(file="YEWA_GV2_map.Rda")

#### Map on ggogle maps

YEWA_GV_map=get_googlemap(center= c(lon=-110, lat = 50), zoom = 3, maptype = "hybrid") %>%
  ggmap() +
  geom_point(data = YEWA_GV_csv, aes(x=Long, y = Lat, colour=tf5045), size=0.8, alpha = 1) +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank())+
  theme(panel.border=element_blank())+
  scale_colour_gradientn(colours = c("deepskyblue3", "yellow", "red","red","red", "red", "red", "red", "red"),breaks = c(0, 0.001,0.003, 0.004,0.005, 0.006, 0.01, 0.085), limits=c(0,0.085), name = "Genomic Vulnerability")+
  labs(y = "Latitude", x = "Longitude")

YEWA_GV_map



