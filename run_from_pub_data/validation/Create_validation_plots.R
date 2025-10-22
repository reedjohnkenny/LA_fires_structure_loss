setwd("~/Desktop/Urban_tree_fire/landscape_analysis/")

library(dplyr)
library(sf)
library(tmap)

eaton_crowns <- st_read("~/OneDrive - Cal Poly/Advanced_GIS_data/Tree_crowns_files/Eaton_crowns_filt.dbf")

eaton_grid <- st_as_sf(st_make_grid(eaton_crowns, cellsize = 250))
eaton_grid$UID <- row_number(eaton_grid)

eaton_perim_minus50m <- st_read("~/OneDrive - Cal Poly/Data/WFIGS_Interagency_Perimeters_YearToDate_-5395415287356828930/") %>% 
  filter(poly_Incid == "Eaton") %>% 
  st_transform(crs = 32611) %>% 
  select(poly_Incid, poly_Featu) %>% 
  st_buffer(-50)

grid_in_poly <- st_join(eaton_grid,eaton_crowns, left = F, join = st_intersects) %>% 
  group_by(UID) %>% 
  summarize(num_trees = n()) %>% 
  st_join(eaton_perim_minus50m, left = F, join = st_within)

tmap_mode("view")

tm_shape(eaton_crowns) +
  tm_polygons(fill = NULL, col = "forestgreen") +
  tm_shape(grid_in_poly) +
  tm_polygons(fill = "num_trees") +
  tm_basemap("Esri.WorldImagery")

rand_plots <- grid_in_poly %>% 
  filter(num_trees > 50) %>% 
  sample_n(9)

st_write(rand_plots, "~/OneDrive - Cal Poly/Advanced_GIS_data/Random_plots/eaton_rand_plots.shp")


tm_shape(rand_plots) +
  tm_polygons(fill = NULL, col = "num_trees") +
tm_shape(eaton_crowns) +
  tm_polygons(fill = NULL, col = "forestgreen") +
  tm_basemap("Esri.WorldImagery")



pal_crowns <- st_read("~/OneDrive - Cal Poly/Advanced_GIS_data/Tree_crowns_files/palisades_Crowns_nofilt.shp")

pal_grid <- st_as_sf(st_make_grid(pal_crowns, cellsize = 250))
pal_grid$UID <- row_number(pal_grid)

palisades_perim_minus50m <- st_read("~/OneDrive - Cal Poly/Data/WFIGS_Interagency_Perimeters_YearToDate_-5395415287356828930/") %>% 
  st_transform(crs = 32611) %>% 
  filter(poly_Incid == "PALISADES") %>%
  select(poly_Incid, poly_Featu) %>% 
  st_buffer(-50)

grid_in_poly <- st_join(pal_grid,pal_crowns, left = F, join = st_intersects) %>% 
  group_by(UID.x) %>% 
  summarize(num_trees = n()) %>% 
  st_join(palisades_perim_minus50m, left = F, join = st_within)

tmap_mode("view")

tm_shape(pal_crowns) +
  tm_polygons(fill = NULL, col = "forestgreen") +
  tm_shape(grid_in_poly) +
  tm_polygons(fill = "num_trees") +
  tm_basemap("Esri.WorldImagery")

rand_plots <- grid_in_poly %>% 
  filter(num_trees > 100) %>% 
  sample_n(4)

tm_shape(pal_crowns) +
  tm_polygons(fill = NULL, col = "forestgreen") +
  tm_shape(rand_plots) +
  tm_polygons(fill = NULL, col = "num_trees") +
  tm_basemap("Esri.WorldImagery")


st_write(rand_plots, "~/OneDrive - Cal Poly/Advanced_GIS_data/Random_plots/palisades_rand_plots.shp")



