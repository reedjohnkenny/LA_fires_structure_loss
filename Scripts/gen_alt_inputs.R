
library(sf)
library(terra)
library(geosphere)
library(nngeo)
library(dplyr)

setwd("~/Desktop/Urban_tree_fire/structure_analysis/")


#bbox <- c(xmin = 396300, xmax = 396450,ymin = 3783800, ymax = 3783990)


all_fprints <- st_read("~/OneDrive - Cal Poly/data/Buildings/LARIAC6_Buildings_2020_eaton.shp") %>% 
  st_transform(32611) 

burned_fprints <- st_read("tmp_data/eaton_burned_fprints.shp")

all_tree_crowns <- st_read("tmp_data/eaton_crowns_burned_polys_v2.shp") 

burned_tree_crowns <- st_read("tmp_data/eaton_burned_tree_crowns_final.shp")


# Calc tree canopy area within 200 m radius of fprint centroid

burned_fprints$cents <- st_centroid(burned_fprints$geometry)

fprints_cents <- st_sf(burned_fprints, sf_column_name = "cents")                        

tcv <- terra::vect(all_tree_crowns)
fpv <- vect(fprints_cents)
cu   <- aggregate(tcv, cores = 4)        
fb2  <- buffer(fpv, width = 200) 
iv2 <- terra::intersect(fb2, cu)
area200 <- terra::expanse(iv2)
area_tree_df_200m  <- data.frame(UID = iv2$UID, area_tree_200m = area200)

# Calc tree canopy area within 100 m radius of fprint centroid

burned_fprints$cents <- st_centroid(burned_fprints$geometry)

fprints_cents <- st_sf(burned_fprints, sf_column_name = "cents")                        

tcv <- terra::vect(all_tree_crowns)
fpv <- vect(fprints_cents)
cu   <- aggregate(tcv, cores = 4)        
fb2  <- buffer(fpv, width = 100) 
iv2 <- terra::intersect(fb2, cu)
area100 <- terra::expanse(iv2)
area_tree_df_100m  <- data.frame(UID = iv2$UID, area_tree_100m = area100)

# Calc tree canopy area within 50 m radius of fprint centroid

burned_fprints$cents <- st_centroid(burned_fprints$geometry)

fprints_cents <- st_sf(burned_fprints, sf_column_name = "cents")                        

tcv <- terra::vect(all_tree_crowns)
fpv <- vect(fprints_cents)
cu   <- aggregate(tcv, cores = 4)        
fb2  <- buffer(fpv, width = 50) 
iv2 <- terra::intersect(fb2, cu)
area50 <- terra::expanse(iv2)
area_tree_df_50m  <- data.frame(UID = iv2$UID, area_tree_50m = area50)



# Calc building area within 200 m radius of fprint centroid

fpv <- vect(fprints_cents)
  fpv_all <- terra::vect(all_fprints)
  fu  <- terra::aggregate(fpv_all, cores = 4)
  fb  <- terra::buffer(fpv, width = 200)
  iv  <- terra::intersect(fb, fu)
  areas <- terra::expanse(iv)
  col_name <- paste0("total_build_200_m")
  df_200 <- data.frame(UID = iv$UID)
  df_200[[col_name]] <- areas
  
  
fpv <- vect(fprints_cents)
  fpv_all <- terra::vect(all_fprints)
  fu  <- terra::aggregate(fpv_all, cores = 4)
  fb  <- terra::buffer(fpv, width = 100)
  iv  <- terra::intersect(fb, fu)
  areas <- terra::expanse(iv)
  col_name <- paste0("total_build_100_m")
  df_100 <- data.frame(UID = iv$UID)
  df_100[[col_name]] <- areas 
  
  
fpv <- vect(fprints_cents)
  fpv_all <- terra::vect(all_fprints)
  fu  <- terra::aggregate(fpv_all, cores = 4)
  fb  <- terra::buffer(fpv, width = 50)
  iv  <- terra::intersect(fb, fu)
  areas <- terra::expanse(iv)
  col_name <- paste0("total_build_50_m")
  df_50 <- data.frame(UID = iv$UID)
  df_50[[col_name]] <- areas 
  

  
  fprints_struc_area <- burned_fprints %>% select(UID) %>% 
    mutate(fprint_area = st_area(geometry)) %>% 
    as.data.frame() %>% 
    select(UID, fprint_area)
  


# Calc number of trees and structures within 200m of focal stucture

burned_fprints$num_trees_200m  <- fprints_cents %>% st_buffer(200) %>% 
  st_intersects(all_tree_crowns) %>% 
  lengths()

burned_fprints$num_builds_200m  <- fprints_cents %>% st_buffer(200) %>% 
  st_intersects(all_fprints) %>% 
  lengths() -1

nums_200 <- burned_fprints %>% as.data.frame() %>% 
  select(UID, num_trees_200m, num_builds_200m)


# Calc Kernel density at diff bandwidths

eaton_tree_ref <- rast(ext(vect(st_centroid(all_tree_crowns))), resolution = 30, crs = st_crs(all_tree_crowns)$wkt)


eaton_tree_dens_500_pdf <- sf.kde(st_centroid(all_tree_crowns), bw = 500, res = 30, ref = eaton_tree_ref, standardize = FALSE, scale.factor = 1, mask = FALSE) 

eaton_tree_dens_500_tph <- eaton_tree_dens_500_pdf * nrow(all_tree_crowns) * 10000

eaton_tree_dens_500 <- terra::extract(eaton_tree_dens_500_tph, burned_fprints, fun  = "mean",
  bind = TRUE
  ) %>%
  st_as_sf() %>%
  select(UID, tree_dens_500 = z)

eaton_tree_dens_300_pdf <- sf.kde(st_centroid(all_tree_crowns), bw = 300, res = 30, ref = eaton_tree_ref, standardize = FALSE, scale.factor = 1, mask = FALSE) 

eaton_tree_dens_300_tph <- eaton_tree_dens_300_pdf * nrow(all_tree_crowns) * 10000

eaton_tree_dens_300 <- terra::extract(eaton_tree_dens_300_tph, burned_fprints, fun  = "mean",
                                      bind = TRUE
) %>%
  st_as_sf() %>%
  select(UID, tree_dens_300 = z)


eaton_build_ref <- rast(ext(vect(st_centroid(all_fprints))), resolution = 30, crs = st_crs(all_fprints)$wkt)


eaton_build_dens_500_pdf <- sf.kde(st_centroid(all_fprints), bw = 500, res = 30, ref = eaton_tree_ref, standardize = FALSE, scale.factor = 1, mask = FALSE) 

eaton_build_dens_500_tph <- eaton_build_dens_500_pdf * nrow(all_fprints) * 10000


eaton_build_dens_500 <- terra::extract(eaton_build_dens_500_tph, burned_fprints, fun  = "mean",
                 bind = TRUE
  ) %>%
  st_as_sf() %>%
  select(UID, build_dens_500 = z)


eaton_build_dens_300_pdf <- sf.kde(st_centroid(all_fprints), bw = 300, res = 30, ref = eaton_tree_ref, standardize = FALSE, scale.factor = 1, mask = FALSE) 

eaton_build_dens_300_tph <- eaton_build_dens_300_pdf * nrow(all_fprints) * 10000


eaton_build_dens_300 <- terra::extract(eaton_build_dens_300_tph, burned_fprints, fun  = "mean",
                                       bind = TRUE
) %>%
  st_as_sf() %>%
  select(UID, build_dens_300 = z)

areas_df <- left_join(area_tree_df_200m, df_200, by = "UID") %>% 
  left_join(df_100, by = "UID") %>% 
  left_join(df_50, by = "UID") %>% 
  left_join(area_tree_df_100m, by = "UID") %>% 
  left_join(area_tree_df_50m, by = "UID") %>% 
  left_join(fprints_struc_area, by = "UID") %>% 
  left_join(nums_200, by = "UID") %>% 
  left_join(as.data.frame(eaton_tree_dens_500), by = "UID") %>% 
  left_join(as.data.frame(eaton_tree_dens_300), by = "UID") %>% 
  left_join(as.data.frame(eaton_build_dens_500), by = "UID") %>% 
  left_join(as.data.frame(eaton_build_dens_300), by = "UID") %>% 
  select(-geometry.y, -geometry.x, -geometry.y.y, -geometry.x.x)

areas_df <- areas_df %>% mutate(mean_perc_struc = total_build_200_m/(pi*200^2),
                                mean_tree_density = num_trees_200m/(pi*200^2)*10000,
                                mean_struc_density = num_builds_200m/(pi*200^2)*10000)

write.csv(as.data.frame(areas_df), "model_inputs/eaton_burned_alt_inputs.csv")




# Palisades



pal_all_fprints <- st_read("~/OneDrive - Cal Poly/data/Buildings/LARIAC6_Buildings_2020_palisades.shp") %>% 
  st_transform(32611)

pal_burned_fprints <- st_read("tmp_data/palisades_burned_fprints.shp")

pal_all_tree_crowns <- st_read("tmp_data/palisades_crowns_burned_polys_v2.shp")

pal_burned_tree_crowns <- st_read("tmp_data/palisades_burned_tree_crowns_final.shp")


# Calc tree canopy area within 200 m radius of fprint centroid

pal_burned_fprints$cents <- st_centroid(pal_burned_fprints$geometry)

pal_fprints_cents <- st_sf(pal_burned_fprints, sf_column_name = "cents")                        

tcv <- terra::vect(pal_all_tree_crowns)
fpv <- vect(pal_fprints_cents)
cu   <- aggregate(tcv, cores = 4)          # union
fb2  <- buffer(fpv, width = 200)  # 2 m buffer
iv2 <- terra::intersect(fb2, cu)
area2 <- terra::expanse(iv2)
area_tree_df_200m  <- data.frame(UID = iv2$UID, area_tree_200m = area2)

# Calc tree canopy area within 100 m radius of fprint centroid

tcv <- terra::vect(pal_all_tree_crowns)
fpv <- vect(pal_fprints_cents)
cu   <- aggregate(tcv, cores = 4)        
fb2  <- buffer(fpv, width = 100) 
iv2 <- terra::intersect(fb2, cu)
area100 <- terra::expanse(iv2)
area_tree_df_100m  <- data.frame(UID = iv2$UID, area_tree_100m = area100)

# Calc tree canopy area within 50 m radius of fprint centroid

tcv <- terra::vect(pal_all_tree_crowns)
fpv <- vect(pal_fprints_cents)
cu   <- aggregate(tcv, cores = 4)        
fb2  <- buffer(fpv, width = 50) 
iv2 <- terra::intersect(fb2, cu)
area50 <- terra::expanse(iv2)
area_tree_df_50m  <- data.frame(UID = iv2$UID, area_tree_50m = area50)



# Calc building area within 200 m radius of fprint centroid

fpv <- vect(pal_fprints_cents)
fpv_all <- terra::vect(pal_all_fprints)
fu  <- terra::aggregate(fpv_all, cores = 4)
fb  <- terra::buffer(fpv, width = 200)
iv  <- terra::intersect(fb, fu)
areas <- terra::expanse(iv)
col_name <- paste0("total_build_200_m")
df_200 <- data.frame(UID = iv$UID)
df_200[[col_name]] <- areas


# Calc building area within 100 m radius of fprint centroid

fpv <- vect(pal_fprints_cents)
fpv_all <- terra::vect(pal_all_fprints)
fu  <- terra::aggregate(fpv_all, cores = 4)
fb  <- terra::buffer(fpv, width = 100)
iv  <- terra::intersect(fb, fu)
areas <- terra::expanse(iv)
col_name <- paste0("total_build_100_m")
df_100 <- data.frame(UID = iv$UID)
df_100[[col_name]] <- areas


# Calc building area within 50 m radius of fprint centroid

fpv <- vect(pal_fprints_cents)
fpv_all <- terra::vect(pal_all_fprints)
fu  <- terra::aggregate(fpv_all, cores = 4)
fb  <- terra::buffer(fpv, width = 50)
iv  <- terra::intersect(fb, fu)
areas <- terra::expanse(iv)
col_name <- paste0("total_build_50_m")
df_50 <- data.frame(UID = iv$UID)
df_50[[col_name]] <- areas




pal_burned_fprints$num_trees_200m  <- pal_fprints_cents %>% st_buffer(200) %>% 
  st_intersects(pal_all_tree_crowns) %>% 
  lengths()

pal_burned_fprints$num_builds_200m  <- pal_fprints_cents %>% st_buffer(200) %>% 
  st_intersects(pal_all_fprints) %>% 
  lengths() -1

nums_200 <- pal_burned_fprints %>% as.data.frame() %>% 
  select(UID, num_trees_200m, num_builds_200m)

# Kernel dens with diff bandwidths

pal_tree_ref <- rast(ext(vect(st_centroid(pal_all_tree_crowns))), resolution = 30, crs = st_crs(pal_all_tree_crowns)$wkt)


pal_tree_dens_500_pdf <- sf.kde(x = st_centroid(pal_all_tree_crowns), bw = 500, ref = pal_tree_ref, standardize = FALSE, scale.factor = 1, mask = FALSE, res = 30)

pal_tree_dens_500_tpha <- pal_tree_dens_500_pdf * nrow(pal_all_tree_crowns) * 10000

pal_tree_dens_500 <- terra::extract(pal_tree_dens_500_tpha, pal_burned_fprints,  fun  = "mean",
                 bind = TRUE
  ) %>%
  st_as_sf() %>%
  select(UID, tree_dens_500 = z)

pal_tree_dens_300_pdf <- sf.kde(st_centroid(pal_all_tree_crowns), ref = pal_tree_ref, bw = 300, res = 30, standardize = FALSE, scale.factor = 1, mask = FALSE) 

pal_tree_dens_300_tpha <- pal_tree_dens_300_pdf * nrow(pal_all_tree_crowns) * 10000

pal_tree_dens_300 <- terra::extract(pal_tree_dens_300_tpha, pal_burned_fprints, fun  = "mean",
                 bind = TRUE
  ) %>%
  st_as_sf() %>%
  select(UID, tree_dens_300 = z)



pal_build_ref <- rast(ext(vect(st_centroid(pal_all_fprints))), resolution = 30, crs = st_crs(pal_all_fprints)$wkt)


pal_build_dens_500_pdf <- sf.kde(st_centroid(pal_all_fprints), ref= pal_build_ref, bw = 500, res = 30, standardize = FALSE, scale.factor = 1, mask = FALSE)
                 
pal_build_dens_500_bph <- pal_build_dens_500_pdf * nrow(pal_all_fprints) * 10000

pal_build_dens_500 <- terra::extract(pal_build_dens_500_bph, pal_burned_fprints, fun = "mean", bind = TRUE) %>%
  st_as_sf() %>%
  select(UID, build_dens_500 = z)


pal_build_dens_300_pdf <- sf.kde(st_centroid(pal_all_fprints), ref = pal_build_ref, bw = 300, res = 30, standardize = FALSE, scale.factor = 1, mask = FALSE)

pal_build_dens_300_bph <- pal_build_dens_300_pdf * nrow(pal_all_fprints) * 10000

pal_build_dens_300 <- terra::extract(pal_build_dens_300_bph, pal_burned_fprints, fun = "mean", bind = TRUE) %>%
  st_as_sf() %>%
  select(UID, build_dens_300 = z)


pal_areas_df <- left_join(area_tree_df_200m, df_200, by = "UID") %>% 
  left_join(df_100, by = "UID") %>% 
  left_join(df_50, by = "UID") %>% 
  left_join(area_tree_df_100m, by = "UID") %>% 
  left_join(area_tree_df_50m, by = "UID") %>% 
  left_join(nums_200, by = "UID") %>% 
  left_join(as.data.frame(pal_tree_dens_500), by = "UID") %>% 
  left_join(as.data.frame(pal_tree_dens_300), by = "UID") %>% 
  left_join(as.data.frame(pal_build_dens_500), by = "UID") %>% 
  left_join(as.data.frame(pal_build_dens_300), by = "UID") %>% 
  select(-geometry.y, -geometry.x, -geometry.y.y, -geometry.x.x)

pal_areas_df <- pal_areas_df %>% mutate(mean_perc_struc = total_build_200_m/(pi*200^2),
                                mean_tree_density = num_trees_200m/(pi*200^2)*10000,
                                mean_struc_density = num_builds_200m/(pi*200^2)*10000)

write.csv(as.data.frame(pal_areas_df), "model_inputs/palisades_burned_alt_inputs.csv")



