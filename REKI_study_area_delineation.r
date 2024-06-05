##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
########## CREATE MULTIPLE GRIDS ACROSS EUROPE TO BE USED AS 'bands' IN DISTRIBUTION SIMULATION  #############
##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
### created by steffen.oppel@vogelwarte.ch on 18 March 2024
### loads pre-pepared dataset based on Red Kite tracking data and MODIS Terra NDVI and CORINE Land Cover habitat data
### this script loads and visualises the gridded data which contain the following columns

####### LIBRARIES REQUIRED
library(tidyverse)
library(dplyr, warn.conflicts = FALSE)
library(lubridate)
library(tmap)
library(sf)
library(terra)
library("rnaturalearth")
library("rnaturalearthdata")
library(gridExtra)

## set root folder for project
#setwd("C:/Users/sop/OneDrive - Vogelwarte/General/ANALYSES/REKI_NDVI")


#### ###########
###  CREATE GRID FOR EUROPE AT 50 KM RESOLUTION  -----------------------------------------------------------------
#### ###########
# create bounding box
bbox <- st_sfc(st_point(c(-10.2, 36)), st_point(c(20, 55)), crs = 4326) %>% st_bbox()
europe <- ne_countries(scale = "medium", returnclass = "sf") %>%
  #mutate(subregion=ifelse(admin %in% c("Czechia","Poland","Slovakia","Hungary"),"Western Europe",subregion)) %>%
  mutate(subregion=ifelse(admin %in% c("Spain","Portugal","Andorra","Italy"),"Western Europe",subregion)) %>%
  filter(continent=="Europe") %>%
  filter(type!="Dependency") %>%
  filter(subregion %in% c("Western Europe")) %>%
  st_crop(y=bbox) %>%
  st_transform(3035)
plot(europe %>% select(tlc))


##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
########## LOOP OVER DIFFERENT VALUES TO GENERATE STUDY MAP SOLUTIONS---------------    #############
##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################

### REDUCED on 5 June to curtail to sensible solution
SIMS<-expand.grid(latstart=c(1400000,1600000),
                  longstart=c(2600000,2800000),
                  flattop=c(TRUE),
                  gridsize=c(400000,500000)) %>%
  mutate(SIM=seq_along(latstart))
plot_list <- list()
for (m in SIMS$SIM) {
  
  grid_EU <- europe %>%
    st_make_grid(cellsize = SIMS$gridsize[m], what = "polygons",
                 offset=c(SIMS$longstart[m],SIMS$latstart[m]),
                 flat_topped=SIMS$flattop[m],
                 square = FALSE) # This statements leads to hexagons
  tab <- st_intersects(grid_EU, europe)
  grid_EU <- st_sf(n_countries = lengths(tab), geometry = st_cast(grid_EU, "MULTIPOLYGON")) %>%
    dplyr::filter(n_countries>0)
  
  plot_list[[m]]<-ggplot(data = europe) +
    geom_sf(fill = "antiquewhite1") +
    geom_sf(data = grid_EU, fill = "firebrick", alpha=0.3) +
    ggtitle(paste0("long: ",SIMS$longstart[m]/100000," | lat: ",SIMS$latstart[m]/100000," | grid: ", SIMS$gridsize[m]/1000))
  
  
}


out<-grid.arrange(grobs=plot_list,ncol=2)
ggsave("output/REKI_Europe_grid_cells.jpg", out, width=10, height=16)



##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
########## USE TMAP TO VISUALISE THE DATA-------------------------------------------    #############
##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################

### since June 2024 tmap_save no longer works with any combination of device or file path

tmap_mode("plot")
map1<-tm_basemap(server="OpenStreetMap") +
  #tm_basemap(server="http://tile.openstreetmap.org/{z}/{x}/{y}.png") +
  tm_shape(grid_EU)  +
  tm_polygons(col = "firebrick", alpha=0.3) +
  tm_layout(title =paste0("longrange: ",bbox[1]," - ",bbox[3]," | latrange: ",bbox[2]," - ", bbox[4]," | gridsize: ", "500"),
          title.size=14, title.color="grey17",title.position=c('right', 'top'))

# tmap_save(map1,filename=paste0("output/MAP_long",bbox[1],"_",bbox[3],"_lat",bbox[2],"_", bbox[4],"_grid", "500",".jpg"),
#           width=9, height=8)
# 
# saveWidget(mapl,file.path(normalizePath(dirname(paste0("C:/Users/sop/OneDrive - Vogelwarte/General/ANALYSES/REKI_NDVI/output/MAP_long",bbox[1],"_",bbox[3],"_lat",bbox[2],"_", bbox[4],"_grid", "500",".html"))),basename(paste0("C:/Users/sop/OneDrive - Vogelwarte/General/ANALYSES/REKI_NDVI/output/MAP_long",bbox[1],"_",bbox[3],"_lat",bbox[2],"_", bbox[4],"_grid", "500",".html"))))
# 
# tmap_save(map1,
#           filename=file.path(normalizePath(dirname(paste0("C:/Users/sop/OneDrive - Vogelwarte/General/ANALYSES/REKI_NDVI/output/MAP_long",bbox[1],"_",bbox[3],"_lat",bbox[2],"_", bbox[4],"_grid", "500",".html")))),
#           width=9, height=8)





########## PUTTING MAPS INTO A QUARTO REPORT---------------    #############

## https://stackoverflow.com/questions/76878132/how-can-i-use-a-loop-in-rmarkdown-to-programmatically-plot-figures-and-text-on-i
library(quarto)
quarto_render("code/StudyAreaMap_compilation.qmd") # all formats
