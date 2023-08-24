# Sumava analysis

This was intended as a practice use case the credits goes to: R.Hladký et al.

Main idea was to perform spatial change detection using remote sensing spectral indeces
Area of interest was Sumava National Park, with timeline selected between 2013 and 2018.
![image](https://github.com/StanislavHerber/Sumava_analysis/assets/134272440/e471573b-4b4c-4115-afed-37a6b3ecf7ca)
credit: https://de.wikipedia.org/wiki/Benutzer:Oligoplectrum

As a software I used either R Studio or QGIS

Area of interest: NP Sumava
![Map1](https://github.com/StanislavHerber/Sumava_analysis/assets/134272440/4babf9cc-e6e2-427d-80ec-e796c81d0efb)

1) Input data were Landsat 8 images (2013,2015,2018)

Workflow description:

a) raster clip
![Map2](https://github.com/StanislavHerber/Sumava_analysis/assets/134272440/d1326ff4-ec06-44df-bd52-17ca342836e6)


b) coordinate system selection

c) make smaller areas for more detail
- this was achieved by firstly creating three points with provided coordinates
- then i made 60 m buffer
- lastly from this buffer i made bounding box, which was used as area of interest

Here is a snapshot of the areas
3 points LAT/LON(13.45196,49.02482,13.52290,48.98473,13.47758,49.04568)

![Map3](https://github.com/StanislavHerber/Sumava_analysis/assets/134272440/881b87cb-8d08-4822-b26e-a3c0a2d359b4)


d) calculate indeces: NDVI, NDMI, TCG, TCW

![Map4](https://github.com/StanislavHerber/Sumava_analysis/assets/134272440/f832e642-12ac-430c-9706-7cc157e29ea9)


Normalized Difference Vegetation Index 

![Map5](https://github.com/StanislavHerber/Sumava_analysis/assets/134272440/9a4dbaaa-416d-4bc8-bbbc-5e8b1b74649b)

Normalized Difference Moisture Index

![Map6](https://github.com/StanislavHerber/Sumava_analysis/assets/134272440/85557305-a319-4b66-8c82-3908a41543be)

Tasseled cap greenness index

![Map7](https://github.com/StanislavHerber/Sumava_analysis/assets/134272440/34595cd5-fc78-4fbe-8a7c-5d635ca3cdbf)

Tasseled cap wetness index

![Map8](https://github.com/StanislavHerber/Sumava_analysis/assets/134272440/69089cc8-5dc2-467e-9f37-4ad5666f3a6a)

Lastly I calculated mean from all indeces inside my areas of interest to determine the change in indeces
using Zonal statistics



Whole R script:
#-----------------------------------------#
#Sumava NP project
#timeline: 
#https://www.mdpi.com/2072-4292/12/12/1914

#Packages
library(raster)
library(sf)
library (tmap)
library (ggplot2)
library (lidR)

#Load input data
setwd("C:/Users/stand/OneDrive/Plocha/Sumava_project")

#Create area of interests (60m)
SA1_dist <- st_point(c(13.45196, 49.02482))
SA2_recov <- st_point(c(13.52290, 48.98473))
SA3_non_dist <- st_point(c(13.47758, 49.04568))

#Convert points to an sf object
SA1_dist_sf <- st_sf(geometry = st_sfc(SA1_dist))
SA2_recov_sf <- st_sf(geometry = st_sfc(SA2_recov))
SA3_non_dist_sf <- st_sf(geometry = st_sfc(SA3_non_dist))
#Set tmap mode to view for interactive maps
tmap_mode("view")
#Plot the points on an interactive map
tm_shape(SA1_dist_sf) +
  tm_dots(col="red")+
  tm_shape(SA2_recov_sf)+
  tm_dots(col="red")+
  tm_shape(SA3_non_dist_sf)+
  tm_dots(col="red")

st_transform(SA1_dist_sf,5514)
st_transform(SA2_recov_sf,5514)
st_transform(SA3_non_dist_sf,5514)

#Make rectangle buffer 60 m
#st_buffer doesnt support rectangular so i made it in qgis
#first i had to export sf files to shapefile
st_write(SA1_dist_sf,"results/SA1_dist.shp")
st_write(SA2_recov_sf,"results/SA2_recov.shp")
st_write(SA3_non_dist_sf,"results/SA3_non_dist.shp")

#Set the coordinate system of the sf object
st_crs(SA1_dist_sf) <- 4326
st_crs(SA2_recov_sf) <- 4326
st_crs(SA3_non_dist_sf) <- 4326

#Transform to a projected coordinate system
point1_sjtsk <- st_transform(SA1_dist_sf, 5514)
point2_sjtsk <- st_transform(SA2_recov_sf, 5514)
point3_sjtsk <- st_transform(SA3_non_dist_sf, 5514)

#Create a 60m buffer around the point
buffer1 <- st_buffer(point1_sjtsk, dist = 60)
buffer2 <- st_buffer(point2_sjtsk, dist = 60)
buffer3 <- st_buffer(point3_sjtsk, dist = 60)

#Create a bounding box around the buffer
bbox1 <- st_as_sfc(st_bbox(buffer1))
bbox2 <- st_as_sfc(st_bbox(buffer2))
bbox3 <- st_as_sfc(st_bbox(buffer3))

#Transform back to the original coordinate system
#bbox1 <- st_transform(bbox1, 4326)

#Export the bounding box as a shapefile
st_write(bbox1, "results/SA1_dist_60m.shp")
st_write(bbox1, "results/SA2_recov_60m.shp")
st_write(bbox1, "results/SA3_non_dist_60m.shp")

#load imagery 
#border of NP Sumava from protected planet
border <- read_sf("data/Hranice_NP/Hranice_NP_sjtsk.shp")
raster2013 <- raster("data/Landsat8_L2/2013/LC08_201306_merge.tif")
raster2015 <- raster("data/Landsat8_L2/2015/LC08_201504_merge.tif")
raster2018 <- raster("data/Landsat8_L2/2018/LC08_201804_merge.tif")

new_crs <- CRS("+init=epsg:5514")

raster2013_proj <- projectRaster(raster2013,crs = new_crs)
raster2015_proj <- projectRaster(raster2015,crs = new_crs)
raster2018_proj <- projectRaster(raster2018,crs = new_crs)

#crop data around border of NP
NPboundary <- st_read("data/hranice_NP/Hranice_NP_sjtsk.shp")

raster2013_proj_clip <- crop(x = raster2013_proj, y = NPboundary)
raster2015_proj_clip <- crop(x = raster2015_proj, y = NPboundary)
raster2018_proj_clip <- crop(x = raster2018_proj, y = NPboundary)

writeRaster(raster2013_proj_clip, format = "GTiff")

#raster or brick
raster2013_proj_clip <- brick("C:/Users/stand/OneDrive/Plocha/Sumava_project/raster_analysis/LC08_201305_merge_sjtsk_clip.tif")
raster2015_proj_clip <- brick("C:/Users/stand/OneDrive/Plocha/Sumava_project/raster_analysis/LC08_201504_merge_sjtsk_clip.tif")
raster2018_proj_clip <- brick("C:/Users/stand/OneDrive/Plocha/Sumava_project/raster_analysis/LC08_201804_merge_sjtsk_clip.tif")
plot(raster2013_proj_clip)

#index analysis
NDVI_2013 <- (raster2013_proj_clip[[5]] - raster2013_proj_clip[[4]]) / (raster2013_proj_clip[[5]] + raster2013_proj_clip[[4]])
NDVI_2015 <- (raster2015_proj_clip[[5]] - raster2015_proj_clip[[4]]) / (raster2015_proj_clip[[5]] + raster2015_proj_clip[[4]])
NDVI_2018 <- (raster2018_proj_clip[[5]] - raster2018_proj_clip[[4]]) / (raster2018_proj_clip[[5]] + raster2018_proj_clip[[4]])

#different way to make NDVI index
ggplot(NDVI_2013)+
  geom_raster()
#
NDMI_2013 <- (raster2013_proj_clip[[5]] - raster2013_proj_clip[[6]]) / (raster2013_proj_clip[[5]] + raster2013_proj_clip[[6]])
NDMI_2015 <- (raster2015_proj_clip[[5]] - raster2015_proj_clip[[6]]) / (raster2015_proj_clip[[5]] + raster2015_proj_clip[[6]])
NDMI_2018 <- (raster2018_proj_clip[[5]] - raster2018_proj_clip[[6]]) / (raster2018_proj_clip[[5]] + raster2018_proj_clip[[6]])

#tasseled cap_greenes,wetness
TCG_2013 <- raster2013_proj_clip[[1]]*(-0.2941)+raster2013_proj_clip[[2]]*(-0.243)+raster2013_proj_clip[[3]]*(-0.5424)+raster2013_proj_clip[[4]]*0.7276+raster2013_proj_clip[[5]]*0.0713+raster2013_proj_clip[[6]]*(-0.1608)
TCG_2015 <- raster2015_proj_clip[[1]]*(-0.2941)+raster2015_proj_clip[[2]]*(-0.243)+raster2015_proj_clip[[3]]*(-0.5424)+raster2015_proj_clip[[4]]*0.7276+raster2015_proj_clip[[5]]*0.0713+raster2015_proj_clip[[6]]*(-0.1608)
TCG_2018 <- raster2018_proj_clip[[1]]*(-0.2941)+raster2018_proj_clip[[2]]*(-0.243)+raster2018_proj_clip[[3]]*(-0.5424)+raster2018_proj_clip[[4]]*0.7276+raster2018_proj_clip[[5]]*0.0713+raster2018_proj_clip[[6]]*(-0.1608)
#set same extent

TCW_2013 <- raster2013_proj_clip[[1]]*0.1511+raster2013_proj_clip[[2]]*0.1973+raster2013_proj_clip[[3]]*0.3283+raster2013_proj_clip[[4]]*0.3407+raster2013_proj_clip[[5]]*(-0.7117)+raster2013_proj_clip[[6]]*(-0.4559)
TCW_2015 <- raster2015_proj_clip[[1]]*0.1511+raster2015_proj_clip[[2]]*0.1973+raster2015_proj_clip[[3]]*0.3283+raster2015_proj_clip[[4]]*0.3407+raster2015_proj_clip[[5]]*(-0.7117)+raster2015_proj_clip[[6]]*(-0.4559)
TCW_2018 <- raster2018_proj_clip[[1]]*0.1511+raster2018_proj_clip[[2]]*0.1973+raster2018_proj_clip[[3]]*0.3283+raster2018_proj_clip[[4]]*0.3407+raster2018_proj_clip[[5]]*(-0.7117)+raster2018_proj_clip[[6]]*(-0.4559)


plot(NDVI_2013, main="NDVI")
plot(NDMI_2013,main="NDMI")
plot(TCG_2013,main="TCG")
plot(TCW_2013,main="TCW")

ndvi_stack <- stack(NDVI_2013,NDVI_2015,NDVI_2018)
ndmi_stack <- stack(NDMI_2013,NDMI_2015,NDMI_2018)
tcg_stack <- stack(TCG_2013,TCG_2015,TCG_2018)
tcw_stack <- stack(TCW_2013,TCW_2015,TCW_2018)

library("RColorBrewer")
library("Lattice")
cols = colorRampPalette(brewer.pal(4,"BuGn"))
levelplot(ndvi_stack,layout=c(3,3),col.regions=cols)
                    
#mean indeces around test sites
area1 <- read_sf("C:/Users/stand/OneDrive/Plocha/Sumava_project/results/SA1_dist_60m.shp")
area2 <- read_sf("C:/Users/stand/OneDrive/Plocha/Sumava_project/results/SA2_recov_60m.shp")
area3 <- read_sf("C:/Users/stand/OneDrive/Plocha/Sumava_project/results/SA3_non_dist_60m.shp")

mean_NDVI_2013_a1 <- extract(NDVI_2013, area1, fun='mean', na.rm=TRUE, df=TRUE, weights = TRUE)
mean_NDVI_2013_a2 <- extract(NDVI_2013, area2, fun='mean', na.rm=TRUE, df=TRUE, weights = TRUE)
mean_NDVI_2013_a3 <- extract(NDVI_2013, area3, fun='mean', na.rm=TRUE, df=TRUE, weights = TRUE)

mean_NDMI_2013_a1 <- extract(NDMI_2013, area1, fun='mean', na.rm=TRUE, df=TRUE, weights = TRUE)
mean_NDMI_2013_a1 <- extract(NDMI_2013, area2, fun='mean', na.rm=TRUE, df=TRUE, weights = TRUE)
mean_NDMI_2013_a1 <- extract(NDMI_2013, area3, fun='mean', na.rm=TRUE, df=TRUE, weights = TRUE)

#-----------------------------------------#
