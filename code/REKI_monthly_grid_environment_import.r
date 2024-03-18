##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
########## LOAD GRIDDED REKI DATA WITH ENVIRONMENTAL LAYERS  #############
##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
### created by steffen.oppel@vogelwarte.ch on 18 March 2024
### loads pre-pepared dataset based on Red Kite tracking data and MODIS Terra NDVI and CORINE Land Cover habitat data
### this script loads and visualises the gridded data which contain the following columns

# -	year, month, date: the year, and month for which the grid cell shows red kite and NDVI data (date is just a composite of year and month)
# -	N_ind: the number of individual red kites (from tracking data) that were within that grid cell in that month
# -	tot_points: the total number of GPS locations (across all individual red kites that were tracked) that were within that grid cell in that month – locations were downsampled to 2 hrs so each point represents 2 ‘red kite hours’
# -	NDVI /EVI: the normalised difference vegetation index and the Landsat Enhanced vegetation index averaged per grid cell for the month and year
# -	propSuitableArea: the proportion of the grid cell that has land cover classes (from Corine Land Cover) that are ‘suitable’ for red kites (i.e. excluding cities, water bodies, forests, glaciers etc.)
# -	totalArea, suitableArea: the absolute values of area (in sqm) and suitable habitat area


####### LIBRARIES REQUIRED
library(tidyverse)
library(dplyr, warn.conflicts = FALSE)
library(lubridate)
library(tmap)
library(sf)





##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
########## READ IN NDVI AND HABITAT DATA    #############
##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
### reading in the rds file from GitHub may not work - clone repo to local drive and read from there
ENVgrid<- readRDS(gzcon(url("https://github.com/Vogelwarte/REKI_NDVI/blob/main/data/REKI_ENV_grid_2016_2023.rds")))
ENVgrid<- readRDS("C:/STEFFEN/OneDrive - Vogelwarte/General/ANALYSES/REKI_NDVI/data/REKI_ENV_grid_2016_2023.rds")

str(ENVgrid)
dim(ENVgrid)
summary(ENVgrid)



##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
########## USE TMAP TO VISUALISE THE DATA    #############
##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################

tmap_mode("view")
tm_basemap(server="OpenStreetMap") +
  tm_shape(ENVgrid)  +
  tm_polygons(col = "tot_points", alpha=0.5) +
  tm_facets(by=c("year"), ncol  =3, showNA = FALSE)


