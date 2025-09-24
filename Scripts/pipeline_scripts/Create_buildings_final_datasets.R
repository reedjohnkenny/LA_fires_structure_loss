library(tidyverse)
library(terra)
library(sf)
library(tmap)
library(mapview)
library(geosphere)
library(spatialEco)
library(sp)
library(nngeo)

setwd("~/Desktop/Urban_tree_fire/structure_analysis/Scripts/pipeline_scripts/")

##########################

# This script reads in the raw buildings shapefiles, filters them and outputs the final layers to be used for calculating proximity


##########################

# Eaton burned

## read in fire perimeter, pull it back by 50 m to avoid edge effects
eaton_perim_minus50m <- st_read("~/OneDrive - Cal Poly/Data/WFIGS_Interagency_Perimeters_YearToDate_-5395415287356828930/") %>% 
  filter(poly_Incid == "Eaton") %>% 
  st_transform(crs=32611) %>% 
  select(poly_Incid, poly_Featu) %>% 
  st_buffer(-50)

## read in DINS points
eaton_dins <- st_read("~/OneDrive - Cal Poly/data/DINS/DINS_2025_Eaton_Public_View.geojson") %>%
  st_transform(crs=32611)

eaton_bound <- st_read("~/OneDrive - Cal Poly/data/town_polygons/Eaton_clip.shp") %>% 
  st_transform(crs = 32611)

eaton_burned_bound <- eaton_bound[1,]

eaton_burned_fprints <- st_read("~/OneDrive - Cal Poly/data/Buildings/LAIAC6_Buildings_2020_Eaton.shp") %>% st_transform(crs = 32611) %>% 
  st_join(eaton_perim_minus50m, left = FALSE) %>% 
  st_join(eaton_burned_bound, left = FALSE) %>% 
  st_join(eaton_dins, left = FALSE)

eaton_burned_fprints$UID <- seq_len(nrow(eaton_burned_fprints))

st_write(eaton_burned_fprints, "../../tmp_data/eaton_burned_fprints.shp", append=FALSE)



# Eaton unburned 

eaton_perim <- st_read("~/OneDrive - Cal Poly/Data/WFIGS_Interagency_Perimeters_YearToDate_-5395415287356828930/") %>% 
  filter(poly_Incid == "Eaton") %>% 
  st_transform(crs=32611) %>% 
  select(poly_Incid, poly_Featu)

eaton_unburned_bound <- eaton_bound[2,]

eaton_unburned_footprints <- st_read("~/OneDrive - Cal Poly/data/Buildings/LAIAC6_Buildings_2020_Eaton.shp") %>% st_transform(crs=32611) %>% 
  st_join(eaton_unburned_bound, left = FALSE) %>% 
  st_difference(eaton_perim)

eaton_unburned_footprints$UID <- seq_len(nrow(eaton_unburned_footprints))

st_write(eaton_unburned_footprints, "../../tmp_data/eaton_unburned_fprints.shp", append=FALSE)

# Palisades

## read in fire perimeter, pull it back by 50 m to avoid edge effects
palisades_perim_minus50m <- st_read("~/OneDrive - Cal Poly/Data/WFIGS_Interagency_Perimeters_YearToDate_-5395415287356828930/") %>% 
  st_transform(crs = 32611) %>% 
  filter(poly_Incid == "PALISADES") %>%
  select(poly_Incid, poly_Featu) %>% 
  st_buffer(-50)

pal_burned_bound <- st_read("~/OneDrive - Cal Poly/Advanced_GIS_data/Maxar/Palisades/Palisades_pre/GIS_FILES/016569561010_01_ORDER_SHAPE.shp") %>% st_transform(crs = 32611)


## read in DINS points
pal_dins <- st_read("~/OneDrive - Cal Poly/data/DINS/DINS_2025_Palisades_Public_View.geojson") %>%
  st_transform(crs = 32611)

## read in building footprints, transform to UTM zone 11
pal_burned_fprints <- st_read("~/OneDrive - Cal Poly/data/Buildings/LARIAC6_Buildings_2020_Palisades.shp") %>% st_transform(crs = 32611) %>%
  st_join(pal_burned_bound, left = FALSE) %>% 
  st_join(palisades_perim_minus50m, left = FALSE) %>% 
  st_join(pal_dins, left = FALSE)

pal_burned_fprints$UID <- seq_len(nrow(pal_burned_fprints))

st_write(pal_burned_fprints, "../../tmp_data/palisades_burned_fprints.shp", append=FALSE)

# Palisades unburned 


palisades_perim <- st_read("~/OneDrive - Cal Poly/Data/WFIGS_Interagency_Perimeters_YearToDate_-5395415287356828930/") %>% 
  st_transform(crs = 32611) %>% 
  filter(poly_Incid == "PALISADES") %>%
  select(poly_Incid, poly_Featu)

pal_unburned_bound <- st_read("~/OneDrive - Cal Poly/data/town_polygons/pal_external.shp") %>% st_transform(32611)

pal_unburned_fprints <- st_read("~/OneDrive - Cal Poly/data/Buildings/LARIAC6_Buildings_2020_Palisades.shp") %>% st_transform(crs = 32611) %>% 
  st_join(pal_unburned_bound, left = FALSE) %>% 
  st_difference(palisades_perim)

pal_unburned_fprints$UID <- seq_len(nrow(pal_unburned_fprints))

st_write(pal_unburned_fprints, "../../tmp_data/palisades_unburned_fprints.shp", append=FALSE)

filter_buildings <- function(fp_path, bounds, perim, dins = NULL, burned = TRUE) {
  fps <- st_read(fp_path) %>% st_transform(32611)
  fps <- st_join(fps, bounds, left = FALSE)
  if (burned) {
    fps <- st_join(fps, st_buffer(perim, -50), left = FALSE)
    if (!is.null(dins)) fps <- st_join(fps, dins, left = FALSE)
  } else {
    fps <- st_filter(fps, predicate = st_disjoint, y = perim)
  }
  fps$UID <- seq_len(nrow(fps))
  return(fps)
}

