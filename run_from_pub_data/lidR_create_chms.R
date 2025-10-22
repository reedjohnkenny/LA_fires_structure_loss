setwd("/Users/Reed/OneDrive - Cal Poly/data/LiDAR/")

library(lidR)
library(lasR)
library(ggplot2)
library(terra)
library(future)
library(sf)
library(tmap)



las_eaton <- readLAScatalog("LARIAC_2023/eaton_LAZ/")

las_check(las_eaton)

# 3. (Optional) If ground points are not already classified, classify them.
#    Here we use the Progressive TIN algorithm to identify ground points.
#    If your LAS already has a ground classification (Classification == 2),
#    you can skip this step.
opt_output_files(las_eaton) <- "/Users/Reed/Desktop/Lidar_tmp/2023_{ID}"
plan(multisession, workers = 14L)
ws <- seq(3,12, 3)
th <- seq(0.1, 1.5, length.out = length(ws))


# 4. Generate a Digital Terrain Model (DTM).
#    We interpolate only the ground-classified points.  
#    Here, resolution = 3 (3ft per pixel); adjust to your needs.
plan(multisession, workers = 15L)
dtm_raster <- rasterize_terrain(las_eaton, res = 1, 
                                algorithm = tin(),
                                pkg = "terra")


# Normalize point cloud, filter points greater than 150 ft

nlas_eaton <- normalize_height(las_eaton, dtm_raster)

las_noisy <- classify_noise(nlas_eaton, algorithm = sor(10, 2))


# filter points below 0 and above 50 m (150 ft)

opt_filter(las_noisy) <- "-drop_z_below 0 -drop_z_above 50 -keep_class 1 2"

opt_output_files(las_noisy) <- "/Users/Reed/Desktop/Lidar_tmp/2023_chm3_{ID}"
dsm_raster <- rasterize_canopy(las_noisy, res = 1, 
                               algorithm = p2r(0.2, na.fill = tin())) 

dsm_raster <- terra::project(dsm_raster, "EPSG:32611")


terra::writeRaster(dsm_raster, "eaton_chm_2023_epsg32611.tif", overwrite = T)



# Palisades 


las_pal <- readLAScatalog("LARIAC_2023/palisades_LAZ/")

las_check(las_pal)

# 3. (Optional) If ground points are not already classified, classify them.
#    Here we use the Progressive TIN algorithm to identify ground points.
#    If your LAS already has a ground classification (Classification == 2),
#    you can skip this step.

plan(multisession, workers = 14L)
ws <- seq(3,12, 3)
th <- seq(0.1, 1.5, length.out = length(ws))


# 4. Generate a Digital Terrain Model (DTM).
#    We interpolate only the ground-classified points.  
#    Here, resolution = 3 (3ft per pixel); adjust to your needs.
plan(multisession, workers = 15L)
dtm_pal <- rasterize_terrain(las_pal, res = 1, 
                             algorithm = tin(),
                             pkg = "terra")


# Normalize point cloud, filter points greater than 50 m (150 ft)
opt_output_files(las_pal) <- "/Users/Reed/Desktop/Lidar_tmp/pal_{ID}"
nlas_pal <- normalize_height(las_pal, dtm_pal)

pal_noisy <- classify_noise(nlas_pal, algorithm = sor(10, 2))


# filter points below 0 and above 50 m (150 ft)

opt_filter(pal_noisy) <- "-drop_z_below 0 -drop_z_above 50 -keep_class 1 2"

opt_output_files(pal_noisy) <- "/Users/Reed/Desktop/Lidar_tmp/palnoisy2_{ID}"

dsm_raster <- rasterize_canopy(pal_noisy, res = 1, 
                               algorithm = p2r(0.2, na.fill = tin())) 

dsm_raster <- terra::project(dsm_raster, "EPSG:32611")


terra::writeRaster(dsm_raster, "pal_chm_2023_epsg32611.tif", overwrite = T)
