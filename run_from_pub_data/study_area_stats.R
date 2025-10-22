library(tidyverse)
library(terra)
library(sf)
library(tmap)
library(mapview)
library(geosphere)
library(spatialEco)
library(sp)
library(nngeo)

## read in fire perimeter, pull it back by 50 m to avoid edge effects
eaton_perim_minus50m <- st_read("~/OneDrive - Cal Poly/Data/WFIGS_Interagency_Perimeters_YearToDate_-5395415287356828930/") %>% 
  filter(poly_Incid == "Eaton") %>% 
  st_transform(crs=32611) %>% 
  select(poly_Incid, poly_Featu) %>% 
  st_buffer(-50)

## read in tree detector points, transform to UTM zone 11 
tree_detector_pts <- st_read("~/OneDrive - Cal Poly/data/TreeCounter/Altadena_tree_detector_Z.shp") %>% 
  st_transform(crs=32611)

## read in DINS points
eaton_dins <- st_read("~/OneDrive - Cal Poly/data/DINS/DINS_2025_Eaton_Public_View.geojson") %>%
  st_transform(crs=32611)

## read in building footprints, transform to UTM zone 11, filter to 50 m inside fire perim, and only those with dins point

study_bound <- st_read("~/OneDrive - Cal Poly/data/town_polygons/Eaton_clip.shp") %>% 
  st_transform(crs = 32611)

fprints <- st_read("tmp_data/eaton_burned_fprints.shp")

fprints$UID <- seq(1:nrow(fprints))

fprints$centroid <- st_centroid(fprints$geometry)

fprints$footprint_area <- st_area(fprints$geometry)

## create layer of just destroyed building footprints
destroyed <- fprints %>% filter(DAMAGE == "Destroyed (>50%)")

tmap_mode("view")
tm_shape(study_bound) +
  tm_polygons(fill = NULL, col = "blue") +
  tm_shape(eaton_perim_minus50m) + 
  tm_polygons()

study_bound <- st_intersection(study_bound[1,], eaton_perim_minus50m)

## Study area

st_area(study_bound)

## Building kernel density mean 

Eaton <- read.csv("~/Desktop/Urban_tree_fire/landscape_analysis/model_inputs/eaton_burned_trees_model_inputs.csv")

Eaton_ext <- read.csv("~/Desktop/Urban_tree_fire/landscape_analysis/model_inputs/eaton_unburned_trees_model_inputs.csv")

mean(Eaton$build_dens, na.rm = T)

Eaton_build_dens_rast <- sf.kde(eaton_dins, bw = 200, res = 30) %>% terra::crop(study_bound, mask = T)

## Tree Kernel Density Mean

mean(Eaton$tree_dens, na.rm = T)

## Number of destroyed structures

Eaton_struc <- read.csv("~/Desktop/Urban_tree_fire/landscape_analysis/model_inputs/eaton_burned_struc_model_inputs.csv")

sum(Eaton_struc$DAMAGE == "Destroyed (>50%)")

## Number of undamaged structures

sum(Eaton_struc$DAMAGE %in% c("No Damage", "Affected (1-9%)"))

## Tree within 2 m of structure

mean(Eaton_struc$area_of_trees_2m)

## Dist to nearest tree from struc

mean(Eaton_struc$distance_to_nearest_tree)


## Dist to other struc

mean(Eaton_struc$distance_to_nearest_building, na.rm = T)



# Palisades


## read in fire perimeter, pull it back by 50 m to avoid edge effects
palisades_perim_minus50m <- st_read("~/OneDrive - Cal Poly/Data/WFIGS_Interagency_Perimeters_YearToDate_-5395415287356828930/") %>% 
  st_transform(crs = 32611) %>% 
  filter(poly_Incid == "PALISADES") %>%
  select(poly_Incid, poly_Featu) %>% 
  st_buffer(-50)

palisades_perim <- st_read("~/OneDrive - Cal Poly/Data/WFIGS_Interagency_Perimeters_YearToDate_-5395415287356828930/") %>% 
  st_transform(crs = 32611) %>% 
  filter(poly_Incid == "PALISADES") %>%
  select(poly_Incid, poly_Featu)

study_bound <- st_read("~/OneDrive - Cal Poly/Advanced_GIS_data/Maxar/Palisades/Palisades_pre/GIS_FILES/016569561010_01_ORDER_SHAPE.shp") %>% st_transform(crs = 32611)


# ## read in tree detector points, transform to UTM zone 11 
# tree_detector_pts <- st_read("~/OneDrive - Cal Poly/data/TreeCounter/Palisades_tree_detector_Z.shp") %>% 
#   st_transform(crs = 32611)

## read in DINS points
dins_pts <- st_read("~/OneDrive - Cal Poly/data/DINS/DINS_2025_Palisades_Public_View.geojson") %>%
  st_transform(crs = 32611)

## read in building footprints, transform to UTM zone 11
fprints <- st_read("tmp_data/palisades_burned_fprints.shp")

fprints$UID <- seq(1:nrow(fprints))

fprints$footprint_area <- st_area(fprints$geometry)


## create layer of just destroyed building footprints
destroyed <- fprints %>% 
  filter(DAMAGE == "Destroyed (>50%)")

study_bound <- st_intersection(study_bound, palisades_perim_minus50m)

## study area

st_area(study_bound)


## Building kernel density mean 

Palisades <- read.csv("~/Desktop/Urban_tree_fire/landscape_analysis/model_inputs/palisades_burned_trees_model_inputs.csv")

Palisades_ext <- read.csv("~/Desktop/Urban_tree_fire/landscape_analysis/model_inputs/palisades_unburned_trees_model_inputs.csv")

mean(Palisades_struc$build_dens, na.rm = T)

Pal_build_dens_rast <- sf.kde(dins_pts, bw = 200, res = 30) %>% terra::crop(study_bound, mask = T)

## Tree Kernel Density Mean

mean(Palisades_struc$tree_dens, na.rm = T)

## Number of destroyed structures

Palisades_struc <- read.csv("~/Desktop/Urban_tree_fire/landscape_analysis/model_inputs/palisades_burned_struc_model_inputs.csv")

sum(Palisades_struc$DAMAGE == "Destroyed (>50%)")

## Number of undamaged structures

sum(Palisades_struc$DAMAGE %in% c("No Damage", "Affected (1-9%)"))

## Tree within 10 m of structure

mean(Palisades_struc$area_of_trees_2m)

## Dist to nearest tree from struc

mean(Palisades_struc$distance_to_nearest_tree)


## Dist to other struc

mean(Palisades_struc$distance_to_nearest_building, na.rm = T)



