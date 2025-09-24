library(dplyr)
library(sf)

setwd("~/Desktop/Urban_tree_fire/landscape_analysis/")

eaton <- st_read("tmp_data/eaton_burned_fprints.shp") %>% 
  st_transform(4326) %>% 
  mutate(longitude = st_coordinates(st_centroid(geometry))[,1],
         latitude = st_coordinates(st_centroid(geometry))[,2],
         start_time = "2025-01-07 18:49:00",
         end_time = "2025-01-08 10:36:00", 
         fire = "eaton") %>% 
  as_tibble() %>% 
  select(fire, UID, longitude, latitude, start_time, end_time)

pal <- st_read("tmp_data/palisades_burned_fprints.shp") %>% 
  st_transform(4326) %>% 
  st_make_valid() %>% 
  mutate(longitude = st_coordinates(st_centroid(geometry))[,1],
         latitude = st_coordinates(st_centroid(geometry))[,2],
         start_time = "2025-01-07 11:00:00",
         end_time = "2025-01-08 20:00:00", 
         fire = "palisades") %>% 
  as_tibble() %>% 
  select(fire, UID, longitude, latitude, start_time, end_time)

all_wind <- rbind(eaton, pal)

write.csv(all_wind, "tmp_data/all_burned_wind_sites.csv")

# run jupyter notebook wind_data

