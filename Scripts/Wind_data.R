setwd("~/Desktop/Urban_tree_fire/structure_analysis/")


library(terra)
library(dplyr)
library(sf)


# Eaton wind angle files

eaton_ang_files <- list.files("wind_data/windninja_output/eaton", pattern = "_ang\\.asc$", full.names = TRUE)
eaton_ang <- terra::rast(eaton_ang_files)

eaton_fprints <- terra::vect("tmp_data/eaton_burned_fprints.shp")

# extract mean wind angle at each footprint
eaton_ang_vals <- terra::extract(eaton_ang, eaton_fprints, fun = "mean", bind = TRUE)
eaton_ang_df <- as.data.frame(eaton_ang_vals)
ang_cols <- names(eaton_ang_df)[grepl("ang", names(eaton_ang_df))]
eaton_ang_df$mean_wind_ang <- rowMeans(eaton_ang_df[, ang_cols], na.rm = TRUE)
eaton_ang_df <- eaton_ang_df %>% select(UID, mean_wind_ang)

# save csv including UID
write.csv(eaton_ang_df, "wind_data/eaton_windninja_ang.csv", row.names = FALSE)

# Eaton wind velocity files

eaton_vel_files <- list.files("wind_data/windninja_output/eaton", pattern = "_vel\\.asc$", full.names = TRUE)
eaton_vel <- terra::rast(eaton_vel_files)

# Extract mean wind velocity at each footprint
eaton_vel_vals <- terra::extract(eaton_vel, eaton_fprints, fun = "mean", bind = TRUE)
eaton_vel_df <- as.data.frame(eaton_vel_vals)
vel_cols <- names(eaton_vel_df)[grepl("vel", names(eaton_vel_df))]
eaton_vel_df$mean_wind_vel <- rowMeans(eaton_vel_df[, vel_cols], na.rm = TRUE)
eaton_vel_df <- eaton_vel_df %>% select(UID, mean_wind_vel)

# save csv including UID
write.csv(eaton_vel_df, "wind_data/eaton_windninja_vel.csv", row.names = FALSE)


# Palisades wind angle files

pal_ang_files <- list.files("wind_data/windninja_output/palisades", pattern = "_ang\\.asc$", full.names = TRUE)
pal_ang <- terra::rast(pal_ang_files)

pal_fprints <- terra::vect("tmp_data/palisades_burned_fprints.shp")

# extract mean wind angle at each footprint
pal_ang_vals <- terra::extract(pal_ang, pal_fprints, fun = "mean", bind = TRUE)
pal_ang_df <- as.data.frame(pal_ang_vals)
ang_cols <- names(pal_ang_df)[grepl("ang", names(pal_ang_df))]
pal_ang_df$mean_wind_ang <- rowMeans(pal_ang_df[, ang_cols], na.rm = TRUE)
pal_ang_df <- pal_ang_df %>% select(UID, mean_wind_ang)

# save csv including UID
write.csv(pal_ang_df, "wind_data/palisades_windninja_ang.csv", row.names = FALSE)

# Palisades wind velocity files

pal_vel_files <- list.files("wind_data/windninja_output/palisades", pattern = "_vel\\.asc$", full.names = TRUE)
pal_vel <- terra::rast(pal_vel_files)

# Extract mean wind velocity at each footprint
pal_vel_vals <- terra::extract(pal_vel, pal_fprints, fun = "mean", bind = TRUE)
pal_vel_df <- as.data.frame(pal_vel_vals)
vel_cols <- names(pal_vel_df)[grepl("vel", names(pal_vel_df))]
pal_vel_df$mean_wind_vel <- rowMeans(pal_vel_df[, vel_cols], na.rm = TRUE)
pal_vel_df <- pal_vel_df %>% select(UID, mean_wind_vel)

# save csv including UID
write.csv(pal_vel_df, "wind_data/palisades_windninja_vel.csv", row.names = FALSE)
