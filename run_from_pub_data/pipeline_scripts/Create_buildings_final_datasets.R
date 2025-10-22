library(tidyverse)
library(terra)
library(sf)
library(tmap)
library(mapview)
library(geosphere)
library(spatialEco)
library(sp)
library(nngeo)

setwd("your working directory")

##########################

# This script reads in the raw buildings shapefiles, filters them and outputs the final layers to be used for calculating proximity


##########################

# Eaton burned

## read in fire perimeter, pull it back by 50 m to avoid edge effects
eaton_perim_minus50m <- st_read("your_download/Perimeters.shp") %>% 
  filter(poly_Incid == "Eaton") %>% 
  st_transform(crs=32611) %>% 
  select(poly_Incid, poly_Featu) %>% 
  st_buffer(-50)

## read in DINS points
eaton_dins <- st_read("your_download/DINS/DINS_2025_Eaton_Public_View.geojson") %>%
  st_transform(crs=32611)

eaton_bound <- st_read("your_download/eaton_burned_bound.shp") %>% 
  st_transform(crs = 32611)

eaton_burned_bound <- eaton_bound

eaton_burned_fprints <- st_read("your_download/LARIAC6_Buildings_2020_eaton.shp") %>% st_transform(crs = 32611) %>% 
  st_join(eaton_perim_minus50m, left = FALSE) %>% 
  st_join(eaton_burned_bound, left = FALSE) %>% 
  st_join(eaton_dins, left = FALSE)

eaton_burned_fprints$UID <- seq_len(nrow(eaton_burned_fprints))

st_write(eaton_burned_fprints, "your_tmp_data/eaton_burned_fprints.shp", append=FALSE)



# Palisades

## read in fire perimeter, pull it back by 50 m to avoid edge effects
palisades_perim_minus50m <- st_read("your_download/Perimeters.shp") %>% 
  st_transform(crs = 32611) %>% 
  filter(poly_Incid == "PALISADES") %>%
  select(poly_Incid, poly_Featu) %>% 
  st_buffer(-50)

pal_burned_bound <- st_read("your_download/palisades_burned_bound.shp") %>% st_transform(crs = 32611)


## read in DINS points
pal_dins <- st_read("your_download/DINS_2025_Palisades_Public_View.geojson") %>%
  st_transform(crs = 32611)

## read in building footprints, transform to UTM zone 11
pal_burned_fprints <- st_read("your_download/LARIAC6_Buildings_2020_Palisades.shp") %>% st_transform(crs = 32611) %>%
  st_join(pal_burned_bound, left = FALSE) %>% 
  st_join(palisades_perim_minus50m, left = FALSE) %>% 
  st_join(pal_dins, left = FALSE)

pal_burned_fprints$UID <- seq_len(nrow(pal_burned_fprints))

st_write(pal_burned_fprints, "your_tmp_data/palisades_burned_fprints.shp", append=FALSE)

