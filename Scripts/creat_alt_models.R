library(spaMM)
library(dplyr)
library(INLA)
library(tidyr)

source("~/Desktop/Urban_tree_fire/structure_analysis/Scripts/functions/compare_AIC.R")

## Eaton struc

setwd("~/Desktop/Urban_tree_fire/structure_analysis/")

eaton_wind <- read.csv("wind_data/all_hrrr.csv") %>% 
  filter(point_fire == "eaton") %>% group_by(point_UID) %>% 
  summarize(mean_wind_sp = mean(wdsp), 
            mean_wind_dir = mean(wdir))

eaton_areas <- read.csv("model_inputs/eaton_burned_alt_inputs.csv")


Eaton_struc <- read.csv("model_inputs/eaton_burned_struc_model_inputs.csv") %>% 
  select(-X) %>% 
  filter(DAMAGE %in% c("Destroyed (>50%)", "No Damage", "Affected (1-9%)")) %>% 
  na.omit() %>% 
  left_join(eaton_wind, by = c("UID" = "point_UID")) %>% 
  left_join(eaton_areas, by = "UID")

Eaton_struc <- Eaton_struc %>% mutate(angular_diff_build = pmin(abs(bearing_to_nearest_building - mean_wind_dir), 360 - abs(bearing_to_nearest_building - mean_wind_dir)), angular_diff_tree = pmin(abs(bearing_to_nearest_tree - mean_wind_dir), 360 - abs(bearing_to_nearest_tree - mean_wind_dir)), destroyed = if_else(DAMAGE %in% c("No Damage", "Affected (1-9%)"),0,1))

cols <- c("distance_to_nearest_tree", "distance_to_nearest_building", "angular_diff_tree", "number_of_trees_2m", "number_of_trees_5m", "number_of_trees_10m",  "area_of_trees_2m",  "area_of_trees_5m",  "area_of_trees_10m", "area_ext_build_2", "area_ext_build_10", "area_ext_build_25", "area_ext_build_50", "area_ext_build_100", "angular_diff_build", "build_dens", "tree_dens", "mean_perc_tree", "mean_perc_struc", "area_tree_200m", "total_build_200_m", "fprint_area")

names(cols) <- cols

Eaton_struc_scale <- mutate(Eaton_struc, across(all_of(cols), ~ scale(.x, center = F)))

spd <- sp::SpatialPointsDataFrame(
  coords = Eaton_struc[, c("utm_x", "utm_y")],
  data = Eaton_struc)

bound <- inla.nonconvex.hull(spd, convex = -.01, resolution = c(105,80))

mesh <- inla.mesh.2d(boundary = bound, max.edge = 200)

# Create SPDE model
spde <- INLA::inla.spde2.matern(mesh)

spaMM.options(separation_max = 100) # or even higher if suggested


eaton_struc_m_scale <- fitme(destroyed ~ distance_to_nearest_tree * angular_diff_tree + area_of_trees_2m + 
                               area_ext_build_2 + distance_to_nearest_building + angular_diff_build +   mean_wind_sp + 
                               total_build_200_m + area_tree_200m  + fprint_area + IMRF(1 | utm_x + utm_y, model = spde), data = Eaton_struc_scale)

saveRDS(eaton_struc_m_scale, "models/eaton_struc_area_model_scale.rda")

eaton_struc_m <- fitme(destroyed ~ distance_to_nearest_tree * angular_diff_tree + area_of_trees_2m + 
                         area_ext_build_2 + distance_to_nearest_building + angular_diff_build +   mean_wind_sp + 
                         total_build_200_m + area_tree_200m  + fprint_area + IMRF(1 | utm_x + utm_y, model = spde), data = Eaton_struc)

saveRDS(eaton_struc_m, "models/eaton_struc_area_model_unscaled.rda")



# Palisades

# Structures


# Palisades Struc

pal_wind <- read.csv("wind_data/all_hrrr.csv") %>% 
  filter(point_fire == "palisades") %>% group_by(point_UID) %>% 
  summarize(mean_wind_sp = mean(wdsp), 
            mean_wind_dir = mean(wdir))

pal_areas <- read.csv("model_inputs/palisades_burned_alt_inputs.csv")


Palisades_struc <- read.csv("model_inputs/palisades_burned_struc_model_inputs.csv") %>% 
  filter(DAMAGE %in% c("Destroyed (>50%)", "No Damage", "Affected (1-9%)")) %>% 
  left_join(pal_wind, by = c("UID" = "point_UID")) %>% 
  left_join(pal_areas, by = c("UID"))
# make new columns

Palisades_struc <- Palisades_struc %>% mutate(angular_diff_build = pmin(abs(bearing_to_nearest_building - mean_wind_dir), 360 - abs(bearing_to_nearest_building - mean_wind_dir)), angular_diff_tree = pmin(abs(bearing_to_nearest_tree - mean_wind_dir), 360 - abs(bearing_to_nearest_tree - mean_wind_dir)), destroyed = if_else(DAMAGE %in% c("No Damage", "Affected (1-9%)"),0,1)) %>% na.omit()


cols <- c("distance_to_nearest_tree", "distance_to_nearest_building", "angular_diff_tree", "number_of_trees_2m", "number_of_trees_5m", "number_of_trees_10m",  "area_of_trees_2m",  "area_of_trees_5m",  "area_of_trees_10m", "area_ext_build_2", "area_ext_build_10", "area_ext_build_25", "area_ext_build_50", "area_ext_build_100", "angular_diff_build", "build_dens", "tree_dens", "mean_perc_tree", "mean_perc_struc", "area_tree_200m", "total_build_200_m", "fprint_area")

names(cols) <- cols

Palisades_struc_scale <- mutate(Palisades_struc, across(all_of(cols), ~ scale(.x, center = F)))




spd <- sp::SpatialPointsDataFrame(
  coords = Palisades_struc[, c("utm_x", "utm_y")],
  data = Palisades_struc)

bound <- inla.nonconvex.hull(spd, convex = -.05)

mesh <- inla.mesh.2d(boundary = bound, max.edge = 250)

# Create SPDE model
spde <- INLA::inla.spde2.matern(mesh)

spaMM.options(separation_max = 100)

Palisades_struc_m_scaled <- fitme(destroyed ~ distance_to_nearest_tree + angular_diff_tree + area_of_trees_2m + area_ext_build_2 + distance_to_nearest_building + angular_diff_build + total_build_200_m + area_tree_200m + mean_wind_sp + fprint_area + IMRF(1 | utm_x + utm_y, model = spde), data = Palisades_struc_scale)

saveRDS(Palisades_struc_m_scaled, "models/Palisades_struc_area_model_scaled.rda")

Palisades_struc_m <- fitme(destroyed ~ distance_to_nearest_tree + angular_diff_tree + area_of_trees_2m + area_ext_build_2 + distance_to_nearest_building + angular_diff_build + total_build_200_m * area_tree_200m + mean_wind_sp + fprint_area + IMRF(1 | utm_x + utm_y, model = spde), data = Palisades_struc)

saveRDS(Palisades_struc_m, "models/Palisades_struc_area_model.rda")




