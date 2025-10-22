library(terra)
library(dplyr)
library(sf)

setwd("~/OneDrive - Cal Poly/data/LiDAR/Eaton/Registered_DEMs/")

# Eaton

tile1 <- terra::rast("all_pan_tile_1_phase.tif")

tile2 <- terra::rast("all_pan_tile_2_phase.tif")

tile3 <- terra::rast("all_pan_tile_3_phase.tif")

tile4 <- terra::rast("all_pan_tile_4_phase.tif")


all_lid <- terra::sprc(tile1, tile2, tile3, tile4)

eaton_burned_bound <- st_read("~/OneDrive - Cal Poly/data/town_polygons/Eaton_clip.shp")[1,] %>% st_transform(32611)

eaton_burned <- terra::crop(all_lid, eaton_burned_bound) %>% 
  terra::merge()

terra::writeRaster(eaton_burned, "eaton_prefire_CHM_burned_reg.tif", overwrite = T)

eaton_unburned_bound <- st_read("~/OneDrive - Cal Poly/data/town_polygons/Eaton_clip.shp")[2,] %>% st_transform(32611)

eaton_unburned <- terra::crop(all_lid, eaton_unburned_bound) %>% 
  terra::merge()

terra::writeRaster(eaton_unburned, "eaton_prefire_CHM_unburned_reg.tif", overwrite = T)


# Palisades


pal_all <- terra::rast("/Users/reedkenny/OneDrive - Cal Poly/data/LiDAR/Palisades/Registered_DEMS/palisades_prefire_CHM_reg.tif")

pal_burn_bound <- st_read("~/OneDrive - Cal Poly/data/Maxar/Palisades/Palisades_pre/GIS_FILES/016569561010_01_ORDER_SHAPE.shp") %>% st_transform(32611)

pal_burned <- terra::crop(pal_all, pal_burn_bound, mask = T)

writeRaster(pal_burned, "/Users/reedkenny/OneDrive - Cal Poly/data/LiDAR/Palisades/Registered_DEMS/palisades_prefire_CHM_burned_reg.tif", overwrite = T)

pal_unburn_bound <- st_read("~/OneDrive - Cal Poly/data/town_polygons/pal_external.shp") %>% st_transform(32611)

pal_unburned <- terra::crop(pal_all, pal_unburn_bound, mask = T)

writeRaster(pal_unburned, "/Users/reedkenny/OneDrive - Cal Poly/data/LiDAR/Palisades/Registered_DEMS/palisades_prefire_CHM_unburned_reg.tif", overwrite = T)
