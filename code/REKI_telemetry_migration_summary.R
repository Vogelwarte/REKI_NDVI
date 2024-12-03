##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
########## SUMMARISE RED KITE MOVEMENTS BETWEEN HEXAGONS ACROSS EUROPE  #############
##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
### this script uses the latest tracking data from Movebank manipulated in REKI_monthly_migration_index.r
### goal is to summarise data for comparison with NDVI-based predictions by Marius Somveille

### for each 10 day period per year one map with 500 km hexagon grids across Europe

library(tidyverse)
library(sf)
library(terra)
library(ncdf4)
library(tidyterra)
library(ggpubr)


### LOAD TRACKING DATA 
REKI<-readRDS("C:/Users/sop/OneDrive - Vogelwarte/General/DATA/REKI_filtered_15min_tracking_data_20Nov2024.rds")
head(REKI)

###  Create large hexagon grid over Europe ### 
bbox <- st_sfc(st_point(c(-10.2, 35)), st_point(c(30, 58)), crs = "+proj=eqearth +datum=WGS84 +lon_0=-10 +lat_0=36") %>% st_bbox()
d = 450
# Hexagon grid covering the migration system
migration_system_map <- vect("C:/Users/sop/Documents/Steffen/mapdata/ne_10m_land/ne_10m_land.shp") %>% 
  crop(bbox) %>% 
  project("+proj=eqearth +datum=WGS84 +lon_0=-10") %>% st_as_sf()
grid_EU <- migration_system_map %>% 
  st_make_grid(cellsize = d * 1000, what = "polygons",
               offset=c(-10000, 4500000),
               flat_topped=TRUE,
               square = FALSE) %>%
  st_transform("+proj=eqearth +datum=WGS84 +lon_0=-10")
tab <- st_intersects(migration_system_map, grid_EU)
grid_EU <- grid_EU[unique(unlist(tab)),] 
selected_hexs_distrib <- c(10,16,19,20,25,30,35,39,40,44,45,49,52)
grid_EU_REKI <- grid_EU[selected_hexs_distrib,]
grid_EU_REKI <- st_sf(ID = seq_along(selected_hexs_distrib), geometry = st_cast(grid_EU_REKI, "MULTIPOLYGON"))

# Plot the grid over Europe
ggplot() +
  geom_sf(data=migration_system_map) +
  geom_sf(data = grid_EU[selected_hexs_distrib,], fill = "firebrick", alpha=0.3)




###  MANIPULATE TRACKING DATA AND BIN TEMPORALLY AND SPATIALLY ### 
year(REKI$timestamp)<-2020  ## assign all dat to a single year
mois <- c("01","02","03","04","05","06","07","08","09","10","11","12")
jours <- c("01","11","21")

timebins<-expand.grid(y=2020,d=jours,m=mois) %>%
  mutate(start=ymd(paste(y,m,d,sep="-"))) %>%
  mutate(end=dplyr::lead(start)) %>%
  mutate(end=if_else(is.na(end), ymd("2020-12-31"),end)) %>%
  mutate(period_id=seq_along(start)) %>%
  mutate(period=interval(start,end)) %>%
  select(period_id,period,start,end)

REKI_summary<-REKI %>% 
  mutate(mo=month(timestamp), day=day(timestamp)) %>%
  mutate(start=if_else(day<11,ymd(paste("2020",mo,"01", sep="-")),
                      if_else(day<21,ymd(paste("2020",mo,"11", sep="-")),ymd(paste("2020",mo,"21", sep="-"))))) %>%
  left_join(timebins, by="start") %>%
  st_transform("+proj=eqearth +datum=WGS84 +lon_0=-10") %>%
  st_join(grid_EU_REKI,join = st_within) %>%
  st_drop_geometry() %>%
  group_by(ID, period_id) %>%
  summarise(n=length(unique(id)))
write.csv(REKI_summary, "REKI_grid_summary.csv")


# Plot the grid over Europe
grid_EU_plot <- grid_EU_REKI %>% full_join(REKI_summary, by='ID') %>%
  left_join(timebins, by='period_id') %>%
  filter(!is.na(period_id))
saveRDS(grid_EU_plot, "REKI_grid_summary.rds")

ggplot() +
  geom_sf(data=migration_system_map) +
  geom_sf(data = grid_EU_plot, aes(fill = n), alpha=0.6) +
  scale_fill_continuous(type="viridis") +
  facet_wrap(~start)

ggsave("REKI_annual_distribution_telemetry.jpg", width=9, height=8) 
  
  
  
##  Potential spatiotemporal paths across hexagons  ##

# Algorithm to create the spatiotemporal network

# Set parameter values
hexgrid <- grid_EU_REKI

# list of neighbours (including itself)
hexgrid_neighbours <- st_touches(hexgrid)
for(i in 1:length(hexgrid_neighbours)){
  hexgrid_neighbours[[i]] <- c(i, hexgrid_neighbours[[i]])
}
# Make sure the hexagons are ordered along the gradient (here: latitude — from south-west to north-east)
st_coordinates(st_centroid(hexgrid))


# Run the algorithm
spatiotemporal_paths <- list()
for(k in 1:Nreg){
  s = list() 
  s[[1]] <- k
  for(i in 2:timesteps){
    # which hexagon is touching the last occupied hexagons
    #s[[i]] <- unlist(hexgrid_neighbours[s[[i-1]]])
    aaa <- hexgrid_neighbours[s[[i-1]]]
    for(j in 1:(i-1)){
      s[[j]] <- rep(s[[j]], times=unlist(lapply(aaa, length)))
    }
    s[[i]] <- unlist(aaa)
    tokeep <- which(apply(do.call(cbind,s), 1, filtering_fct) == 1)
    for(j in 1:i){
      s[[j]] <- s[[j]][tokeep]
    }
    print(i)
  }
  spat_paths <- do.call(cbind,s)
  # remove paths for which the last hexagon is not a neighbour of the first hexagon (or itself)
  tokeep <- spat_paths[,timesteps] %in% hexgrid_neighbours[[k]]
  #tokeep <- vector()
  #for(i in 1:nrow(spat_paths)){
  #  tokeep[i] <- spat_paths[,timesteps][i] %in% unlist(hexgrid_neighbours[spat_paths[,1][i]])
  #}
  spat_paths <- spat_paths[which(tokeep == T),]
  # remove paths that are not staying at least Smin periods at each seasonal destination
  tokeep <- which(apply(spat_paths, 1, filtering_fct_2) == 1)
  spatiotemporal_paths[[k]] <- spat_paths[tokeep,]
  print(c(k, nrow(spatiotemporal_paths[[k]])))
}
spatiotemporal_paths <- do.call(rbind, spatiotemporal_paths)
write.csv(spatiotemporal_paths, "resources/spatiotemporal_paths_hex3.csv")


## Graph model [deleted code except the code that specifies structure of output data frame]

graph_paths <- as.matrix(read.csv("resources/spatiotemporal_paths_hex3.csv", header = T)[,-1])
ee_path_ID <- which.max(ee_score)
simulated_migr_paths <- c(simulated_migr_paths, ee_path_ID)
simulated_migration_paths <- graph_paths[simulated_migr_paths,]


## Plot results from STDS model ##

# Calculate migratory links for every time steps
simulated_migration_paths_2 <- cbind(simulated_migration_paths, simulated_migration_paths[,1])
migra_links <- vector()
for(t in 1:timesteps){
  links <- which(simulated_migration_paths_2[,t] != simulated_migration_paths_2[,t+1])
  if(length(links) > 1){
    ml = simulated_migration_paths_2[links,][,t:(t+1)] #%>% as_data_frame()
    ml$link <- apply(ml, 1, function(x) paste(x, collapse = " "))
    migra_links <- rbind(migra_links, cbind(t, as.data.frame(table(ml$link))))
  }else if(length(links) == 1){
    ml = simulated_migration_paths_2[links,][t:(t+1)]
    migra_links <- rbind(migra_links, c(t, paste(ml, collapse = " "), 1))
    colnames(migra_links) <- c("t", "Var1", "Freq")
  }
}

# Total weekly migration activity 
total_migration_activity <- vector()
for(i in 1:timesteps){
  total_migration_activity[i] <- sum(as.numeric(migra_links$Freq[which(migra_links$t==i)]))
}
data_to_plot <- data.frame(date = as.Date(seq(1, 360, 10), origin="2015-01-01"),
                           migra_activity = total_migration_activity)
theme_set(theme_bw())
ggplot(data_to_plot, aes(x=date, y=migra_activity)) + 
  geom_line() + 
  geom_point() + xlab("Date") + ylab("Migration activity")

# Seasonal difference in richness
weekly_density <- vector()
for(i in 1:timesteps){
  weekly_dens <- vector()
  for(j in 1:Nreg){
    weekly_dens[j] <- length(which(simulated_migration_paths[,i] == j))
  }
  weekly_density <- rbind(weekly_density, weekly_dens)
}
winter_timestep <- 1
summer_timestep <- 22
density_winter <- weekly_density[winter_timestep,]
density_summer <- weekly_density[summer_timestep,]
data_plot <- data.frame(band = 1:13,
                        winter_density = density_winter,
                        summer_density = density_summer,
                        seasonal_difference = density_summer - density_winter)
theme_set(theme_classic())
ggplot(data_plot) +
  geom_smooth(aes(x=band, y=seasonal_difference),span = 0.5) +
  geom_point(aes(x=band, y=seasonal_difference)) +
  xlab("Latitude") + ylab("Seasonal difference in bird density")

# Plot weekly links on map
migra_links2 <- cbind(migra_links, t(apply(migra_links, 1, function(x) as.numeric(unlist(strsplit(as.character(x[2]), " "))))))
colnames(migra_links2) <- c("week", "unique_ID", "inds", "origin", "destination")
migra_links2 <- migra_links2 %>% dplyr::select(-unique_ID)
hex_centroids <- st_coordinates(st_centroid(grid_EU[selected_hexs_distrib,])) %>% as.data.frame()
myTheme <- theme(plot.title = element_text(size=30),
                 legend.text = element_text(size = 30), 
                 legend.title = element_text(size = 30), 
                 legend.key.size = unit(2, 'cm'))
g_flow_bins <- list()
for(t in 1:timesteps){
  migra_links_w <- migra_links2 %>% filter(week==t)
  migra_links_w_plot <- data.frame(
    x1 = as.numeric(hex_centroids$X[migra_links_w$origin]),
    x2 = as.numeric(hex_centroids$X[migra_links_w$destination]),
    y1 = as.numeric(hex_centroids$Y[migra_links_w$origin]),
    y2 = as.numeric(hex_centroids$Y[migra_links_w$destination]),
    flow = as.numeric(migra_links_w$inds)
  )
  g_flow_bins[[t]] <- ggplot() +
    geom_sf(data=migration_system_map) +
    geom_sf(data = grid_EU[selected_hexs_distrib,], fill = "firebrick", alpha=0.3) + theme_void() + xlim(c(-15500, 3000000)) +
    ggtitle(format(ymd(paste0("15", mois[floor(t/3)+1], jours[round(((t/3) - floor(t/3)) *3) + 1])), "%B %d")) +
    geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2, linewidth=flow),
                 data = migra_links_w_plot[order(migra_links_w_plot$flow),],
                 #lineend = "round", linejoin = "round",
                 arrow = ggplot2::arrow(length = unit(0.04, "npc"))) + myTheme +
    scale_linewidth(range=c(0,ifelse(nrow(migra_links_w_plot) > 0, max(migra_links_w_plot$flow)/3, 1)), name = "Migration flow")
}

pdf(file=paste0("results/flow_results.pdf"), width=45, height=45)
ggarrange(g_flow_bins[[1]], g_flow_bins[[2]], g_flow_bins[[3]], g_flow_bins[[4]], g_flow_bins[[5]], g_flow_bins[[6]],
          g_flow_bins[[7]], g_flow_bins[[8]], g_flow_bins[[9]], g_flow_bins[[10]], g_flow_bins[[11]], g_flow_bins[[12]],
          g_flow_bins[[13]], g_flow_bins[[14]], g_flow_bins[[15]], g_flow_bins[[16]], g_flow_bins[[17]], g_flow_bins[[18]],
          g_flow_bins[[19]], g_flow_bins[[20]], g_flow_bins[[21]], g_flow_bins[[22]], g_flow_bins[[23]], g_flow_bins[[24]],
          g_flow_bins[[25]], g_flow_bins[[26]], g_flow_bins[[27]], g_flow_bins[[28]], g_flow_bins[[29]], g_flow_bins[[30]],
          g_flow_bins[[31]], g_flow_bins[[32]], g_flow_bins[[33]], g_flow_bins[[34]], g_flow_bins[[35]], g_flow_bins[[36]],
          nrow=6, ncol=6, common.legend = T, legend="none")
dev.off()





