
library(dplyr)
library(terra)
library(sf)
library(tmap)

## Eaton Pre

setwd("~/OneDrive - Cal Poly/data/Maxar/Eaton/PRE/")


rast1_crop <- st_read("~/OneDrive - Cal Poly/data/Maxar/Eaton/PRE/ms1_crop.shp") %>% 
  st_transform(crs = 32611)

maxar_pre_rast1 <- terra::rast("~/OneDrive - Cal Poly/data/Maxar/Eaton/PRE/10400100A17E8600-ms.tif") %>% terra::crop(rast1_crop[1,], mask = T)

maxar_pre_rast2 <- terra::rast("~/OneDrive - Cal Poly/data/Maxar/Eaton/PRE/10400100A17E8600-ms (1).tif")

maxar_pre_rast3 <- terra::rast("~/OneDrive - Cal Poly/data/Maxar/Eaton/PRE/10400100A17E8600-ms (2).tif")

maxar_pre_rast4 <- terra::rast("016569561020_01/016569561020_01_P001_MUL/25JAN01182958-M3DS-016569561020_01_P001.TIF") %>% terra::crop(rast1_crop[2,], mask = T)

maxar_pre_rast5 <- terra::rast("016569561030_01/016569561030_01_P001_MUL/25JAN01182958-M3DS-016569561030_01_P001.TIF")


#normalize bands to between 0 and 1

library(terra)

library(terra)

normalize_raster <- function(raster) {
  if (!inherits(raster, "SpatRaster")) {
    stop("Input must be a SpatRaster object.")
  }
  
  # Initialize list to store normalized layers
  norm_layers <- list()
  
  # Normalize each band
  for (i in 1:nlyr(raster)) {
    band <- raster[[i]]
    min_val <- global(band, fun = "min", na.rm = TRUE)[1,1]
    max_val <- global(band, fun = "max", na.rm = TRUE)[1,1]
    norm_band <- (band - min_val) / (max_val - min_val)
    norm_layers[[i]] <- norm_band
  }
  
  # Combine layers using c() which is safe for SpatRaster objects
  normalized <- do.call(c, norm_layers)
  names(normalized) <- c("ms1", "ms2", "ms3", "ms4", "ms5", "ms6", "ms7", "ms8")
  return(normalized)
}


maxar_pre_rast1_norm <- normalize_raster(maxar_pre_rast1)

maxar_pre_rast2_norm <- normalize_raster(maxar_pre_rast2)

maxar_pre_rast3_norm <- normalize_raster(maxar_pre_rast3)

maxar_pre_rast4_norm <- normalize_raster(maxar_pre_rast4)

maxar_pre_rast5_norm <- normalize_raster(maxar_pre_rast5)


maxar_pre_all_norm <- terra::merge(maxar_pre_rast2_norm, maxar_pre_rast3_norm, maxar_pre_rast1_norm, maxar_pre_rast4_norm, maxar_pre_rast5_norm)


## read in fire perimeter
eaton_perim <- st_read("~/OneDrive - Cal Poly/Data/WFIGS_Interagency_Perimeters_YearToDate_-5395415287356828930/") %>% 
  filter(poly_Incid == "Eaton") %>% 
  st_transform(crs=32611) %>% 
  select(poly_Incid, poly_Featu)

tmap_mode("plot")

tm_shape(maxar_pre_all_norm$ms8) +
  tm_raster() +
  tm_shape(eaton_perim) +
  tm_polygons(fill = NULL, col = "red")

maxpre_all_ndvi_rast <- ((maxar_pre_all_norm$ms8 - maxar_pre_all_norm$ms5)/(maxar_pre_all_norm$ms8 + maxar_pre_all_norm$ms5))


tm_shape(maxpre_all_ndvi_rast) +
  tm_raster(style = "cont") +
  tm_shape(eaton_perim) +
  tm_polygons(fill = NULL, col = "red")

terra::writeRaster(maxar_pre_all_norm, "~/OneDrive - Cal Poly/Advanced_GIS_data/Maxar/Eaton/PRE/eaton_maxar_pre_all.TIF", overwrite=TRUE)

terra::writeRaster(maxar_pre_all_norm, "~/OneDrive - Cal Poly/data/Maxar/Eaton/PRE/eaton_maxar_pre_all.tif")

terra::writeRaster(maxpre_all_ndvi_rast, "~/OneDrive - Cal Poly/Advanced_GIS_data/Maxar/Eaton/PRE/eaton_maxar_pre_all_ndvi.tif", overwrite=TRUE)


## Eaton post

crop1 <- st_read("~/OneDrive - Cal Poly/data/Maxar/Eaton/PRE/ms1_crop.shp") %>% st_transform(crs = 32611)

study_bound <- st_read("~/OneDrive - Cal Poly/data/town_polygons/Eaton_clip.shp") %>% 
  st_transform(crs = 32611)


files <- c("../POST/103001010C360000-ms.tif", 
           "../POST/10400100A26E9900-ms (2).tif",
           "../POST/10400100A26E9900-ms (1).tif",
           "../POST/10400100A26E9900-ms.tif",
           "../POST/200007937065_01_200007937065_01_P001_MUL_25JAN14183036-M2AS-200007937065_01_P001_TOA_cor.TIF")

# 2. Pre-allocate a list to hold normalized rasters

norm <- vector("list", length(files))

# 3. Loop: read, normalize, and store

for (i in seq_along(files)) {
  # read the i-th raster
  r <- rast(files[i])
  
  # apply your normalization function
  norm[[i]] <- normalize_raster(r)
}

norm[[6]] <- terra::crop(norm[[5]], crop1[3,], mask = T)

muls <- terra::merge(norm[[6]], norm[[1]], norm[[2]], norm[[3]], norm[[4]])

eaton_ndvi_post <- ((muls$ms8 - muls$ms5)/(muls$ms8 + muls$ms5))

tmap_mode("view")
tm_shape(eaton_ndvi_post) +
  tm_raster() +
  tm_shape(crop1[3,]) +
  #tm_polygons(fill = NULL, col = "red") +
  tm_shape(study_bound) +
  tm_polygons(fill = NULL, col = "yellow")



terra::writeRaster(muls, "Eaton_multispec_post_all.tif")

terra::writeRaster(eaton_ndvi_post, "Eaton_ndvi_post_all.tif")
