setwd("~/Desktop/Urban_tree_fire/structure_analysis/")

#install.packages("lidR")
library(lidR)
library(tidyverse)
library(sf)
library(tmap)
library(lasR)
library(terra)


# Function to segment tree canopies

segment_tree_canopies <- function(chm, boundary, footprints, ndvi_rast) {
  chm_agg <- terra::aggregate(chm, fact = 4)
  chm_crop <- terra::crop(chm_agg, terra::vect(boundary), mask = T)
  chm_crop_mask <- terra::mask(chm_crop, terra::vect(footprints), inverse = T)
  # Post-processing median filter
  kernel <- matrix(1,3,3)
  chm_sm <- terra::focal(chm_crop_mask, w = kernel, fun = median, na.rm = TRUE)
  
  f <- function(x) {x * 0.15 + 5}
  heights <- seq(0,30,5)
  ws <- f(heights)
  ttops_chm <- locate_trees(chm_sm, lmf(f))
  message("running segmentation algorithm")
  algo <- dalponte2016(chm_sm, ttops_chm, 
                       th_tree = 1, th_seed = 0.45,
                       th_cr = 0.55,
                       max_cr = 30)
  
  crowns_all <- algo()
  
  crowns_polys <- terra::as.polygons(crowns_all, round = F)
  
  crowns_polys_ndvi <- terra::extract(ndvi_rast, crowns_polys, fun = mean, bind = T) %>%  
    st_as_sf() %>%  
    select(focal_median, maxar_post_ndvi = ms8) 
  
  crowns_polys_filt <- crowns_polys_ndvi %>% filter(maxar_post_ndvi > 0.5)
  return(list(crowns_polys_filt, crowns_all, ttops_chm)) 
}

### Palisades

palisades_pre_chm <- terra::rast("~/OneDrive - Cal Poly/data/LiDAR/Palisades/Registered_DEMS/24OCT20184323-P3DS-016569561010_01_P001_phase.tif")

pal_buildings <- st_read("~/OneDrive - Cal Poly/data/Buildings/LARIAC6_Buildings_2020_Palisades.shp") %>%
  st_transform(crs="epsg:32611") %>%
  st_make_valid() %>%
  st_union()

pal_study_bound <- st_read("~/OneDrive - Cal Poly/Advanced_GIS_data/Maxar/Palisades/Palisades_pre/GIS_FILES/016569561010_01_ORDER_SHAPE.shp") %>% st_transform(crs = 32611)

maxpre_rast1 <- terra::rast("~/OneDrive - Cal Poly/Advanced_GIS_data/Maxar/Palisades/Palisades_pre/Multi_spec/24OCT20184323-M3DS-016569561010_01_P001.TIF") 

maxpre_rast1 <- terra::project(maxpre_rast1, "epsg:32611")

names(maxpre_rast1) <-  c("ms1", "ms2", "ms3", "ms4", "ms5", "ms6", "ms7", "ms8")

maxpre_ndvi_rast <- ((maxpre_rast1$ms8 - maxpre_rast1$ms5)/(maxpre_rast1$ms8 + maxpre_rast1$ms5)) 

pal_segs <- segment_tree_canopies(chm = palisades_pre_chm, boundary = pal_study_bound, footprints = pal_buildings, ndvi_rast = maxpre_ndvi_rast)

#pal_segs

st_write(st_as_sf(pal_segs[[1]]), "tmp_data/palisades_crowns_burned_polys_v2.shp", append = FALSE)

writeRaster(pal_segs[[2]], filename = "tmp_data/palisades_crown_tile_all_v2.tif", overwrite=TRUE)

st_write(pal_segs[[3]], "tmp_data/palisades_ttops_all.csv", layer_options = "GEOMETRY=AS_XY", append = FALSE)

## Palisades external

palisades_ext_pre_chm <- terra::rast("/Users/reedkenny/OneDrive - Cal Poly/data/LiDAR/Palisades/Registered_DEMS/palisades_prefire_CHM_unburned_reg.tif") %>% 
  project("epsg:32611")


Palisades_ext_bound <- st_read("~/OneDrive - Cal Poly/data/town_polygons/pal_external.shp") %>% st_transform(crs = 32611)


Palisades_ext_buildings <- st_read("~/OneDrive - Cal Poly/data/Buildings/LARIAC6_Buildings_2020_Palisades.shp") %>% st_transform(crs = 32611)

pal_ext_mulispec <- terra::rast("~/OneDrive - Cal Poly/Advanced_GIS_data/Maxar/Palisades/Palisades_pre/Multi_spec/24OCT20184322-M2AS-200008990764_01_P001_TOA_cor.TIF") %>%  terra::project("epsg:32611")

names(pal_ext_mulispec) <-  c("ms1", "ms2", "ms3", "ms4", "ms5", "ms6", "ms7", "ms8")

pal_ext_ndvi_rast <- ((pal_ext_mulispec$ms8 - pal_ext_mulispec$ms5)/(pal_ext_mulispec$ms8 + pal_ext_mulispec$ms5)) 

# tm_basemap("Esri.WorldImagery") +
# tm_shape(palisades_ext_pre_chm) +
#   tm_raster() +
#   tm_shape(pal_ext_ndvi_rast) +
#   tm_raster() +
#   tm_shape(Palisades_ext_buildings) + 
#   tm_polygons(fill = NULL) +
#   tm_shape() +
#   tm_polygons(fill = NULL, col = "green")
# 
# tm_basemap("Esri.WorldImagery") + 
#   #tm_shape(pal_ext_ndvi_rast) +
#   #tm_raster() +
#   tm_shape(crowns_polys_ndvi) +
#   tm_polygons(fill = NULL, col = "maxar_post_ndvi")

pal_ext_segs <- segment_tree_canopies(chm = palisades_ext_pre_chm, boundary = Palisades_ext_bound, footprints = Palisades_ext_buildings, ndvi_rast = pal_ext_ndvi_rast)

st_write(st_as_sf(pal_ext_segs[[1]]), "tmp_data/palisades_crowns_unburned_polys_v2.shp", append = FALSE)

writeRaster(pal_ext_segs[[2]], filename = "tmp_data/palisades_crown_tile_all_ext_v2.tif", overwrite=TRUE)

st_write(pal_ext_segs[[3]], "tmp_data/palisades_ttops_all_ext.csv", layer_options = "GEOMETRY=AS_XY", append = FALSE)

chm_agg <- terra::aggregate(palisades_ext_pre_chm, fact = 4)

chm_crop <- terra::crop(chm_agg, terra::vect(Palisades_ext_bound), mask = T)
chm_crop_mask <- terra::mask(chm_crop, terra::vect(Palisades_ext_buildings), inverse = T)
# Post-processing median filter
kernel <- matrix(1,3,3)
chm_sm <- terra::focal(chm_crop_mask, w = kernel, fun = median, na.rm = TRUE)

f <- function(x) {x * 0.15 + 5}
heights <- seq(0,30,5)
ws <- f(heights)
ttops_chm <- locate_trees(chm_sm, lmf(f))
message("running segmentation algorithm")
algo <- dalponte2016(chm_sm, ttops_chm, 
                     th_tree = 1, th_seed = 0.45,
                     th_cr = 0.55,
                     max_cr = 30)

crowns_all <- algo()

crowns_polys <- terra::as.polygons(crowns_all, round = F)

crowns_polys_ndvi <- terra::extract(pal_ext_ndvi_rast, crowns_polys, fun = mean, bind = T) %>%  
  st_as_sf() %>%  
  select(focal_median, maxar_post_ndvi = ms8) 

crowns_polys_filt <- crowns_polys_ndvi %>% filter(maxar_post_ndvi > 0.5)



## eaton


eaton_pre_chm <- terra::rast("/Users/reedkenny/OneDrive - Cal Poly/data/LiDAR/Eaton/Registered_DEMS/eaton_prefire_CHM_burned_reg.tif")

eaton_study_bound <- st_read("~/OneDrive - Cal Poly/data/town_polygons/Eaton_clip.shp") %>% st_transform(crs = 32611)

eaton_burn_bound <- eaton_study_bound[1,]

eaton_buildings <- st_read("~/OneDrive - Cal Poly/data/Buildings/LAIAC6_Buildings_2020_Eaton.shp") %>%
  st_transform(crs=32611) %>%
  st_make_valid() %>%
  st_union()

eaton_pre_ndvi_rast <- rast("~/OneDrive - Cal Poly/Advanced_GIS_data/Maxar/Eaton/PRE/eaton_maxar_pre_all_ndvi.tif")

eaton_segs <- segment_tree_canopies(chm = eaton_pre_chm, boundary = eaton_burn_bound, footprints = eaton_buildings, ndvi_rast = eaton_pre_ndvi_rast)


st_write(st_as_sf(eaton_segs[[1]]), "tmp_data/eaton_crowns_burned_polys_v2.shp", append = FALSE)

writeRaster(eaton_segs[[2]], filename = "tmp_data/Eaton_crown_tile_all_v2.tif", overwrite=TRUE)

st_write(eaton_segs[[3]], "tmp_data/eaton_ttops_all.csv", layer_options = "GEOMETRY=AS_XY", append = FALSE)


## Eaton external

eaton_ext_pre_chm <- terra::rast("~/OneDrive - Cal Poly/data/LiDAR/Eaton/Registered_DEMS/eaton_prefire_CHM_unburned_reg.tif") %>% project("epsg:32611")

Eaton_unburn_bound <- eaton_study_bound[2,]

eaton_ext_buildings <- st_read("~/OneDrive - Cal Poly/data/Buildings/LAIAC6_Buildings_2020_Eaton.shp") %>%
  st_transform(crs=32611) %>%
  st_make_valid() %>%
  st_union()


eaton_ext_pre_ndvi_rast <- rast("~/OneDrive - Cal Poly/Advanced_GIS_data/Maxar/Eaton/PRE/eaton_maxar_pre_all_ndvi.tif") %>% project("epsg:32611")

eaton_ext_segs <- segment_tree_canopies(eaton_ext_pre_chm, Eaton_unburn_bound, eaton_ext_buildings, eaton_ext_pre_ndvi_rast)

st_write(st_as_sf(eaton_ext_segs[[1]]), "tmp_data/eaton_crowns_unburned_polys_v2.shp", append = FALSE)

writeRaster(eaton_ext_segs[[2]], filename = "tmp_data/Eaton_crown_tile_all_ext.tif", overwrite=TRUE)

st_write(eaton_ext_segs[[3]], "tmp_data/eaton_ttops_all_ext.csv", layer_options = "GEOMETRY=AS_XY", append = F)
