setwd("~/Desktop/Urban_tree_fire/structure_analysis/")


library(tidyverse)
library(sf)  
library(terra)
library(mapview)

# Palisades stratified random validation points

pal_crowns <- st_read("tmp_data/palisades_burned_tree_crowns_final.shp")

pal_bound <- st_read("tmp_data/palisades_burned_fprints.shp") %>% 
  st_union() %>% 
  st_buffer(50) %>% 
  st_cast("POLYGON") %>% 
  st_as_sf()

pal_crowns_crop <- st_crop(st_union(pal_crowns), pal_bound) %>%
  st_union()

pal_non_crowns <- st_difference(pal_bound, pal_crowns_crop) %>% 
  st_union()

canopy_points <- st_sample(pal_crowns_crop, size = 250, type = "random") %>% 
  st_as_sf()

canopy_points$mapped_class <- "canopy"

non_canopy_points <- st_sample(pal_non_crowns, size = 250, type = "random") %>% 
  st_as_sf()

non_canopy_points$mapped_class <- "non_canopy"


validation_points <- rbind(canopy_points, non_canopy_points)

st_write(validation_points, "Scripts/validation/validation_points.shp")


# Eaton stratified random validation points

eaton_crowns <- st_read("tmp_data/eaton_burned_tree_crowns_final.shp")


eaton_bound <- st_read("tmp_data/eaton_burned_fprints.shp") %>% 
  st_union() %>% 
  st_buffer(50) %>% 
  st_cast("POLYGON") %>% 
  st_as_sf()

eaton_crowns_crop <- st_crop(st_union(eaton_crowns), eaton_bound) %>% 
  st_union

eaton_non_crowns <- st_difference(eaton_bound, eaton_crowns_crop) %>% 
  st_union()

eaton_canopy_points <- st_sample(eaton_crowns_crop, size = 250, type = "random") %>% 
  st_as_sf()

eaton_canopy_points$mapped_class <- "canopy"

eaton_non_canopy_points <- st_sample(eaton_non_crowns, size = 250, type = "random") %>% 
  st_as_sf()

eaton_non_canopy_points$mapped_class <- "non_canopy"


eaton_validation_points <- rbind(eaton_canopy_points, eaton_non_canopy_points)

st_write(eaton_validation_points, "Scripts/validation/eaton_validation_points.shp")


