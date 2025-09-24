library(spaMM)
library(dplyr)
library(INLA)
library(tidyr)

source("~/Desktop/Urban_tree_fire/landscape_analysis/Scripts/functions/compare_AIC.R")

## Eaton struc

setwd("~/Desktop/Urban_tree_fire/landscape_analysis/")

eaton_wind <- read.csv("wind_data/all_hrrr.csv") %>% 
  filter(point_fire == "eaton") %>% group_by(point_UID) %>% 
  summarize(mean_wind_sp = mean(wdsp), 
            mean_wind_dir = mean(wdir))


Eaton_struc <- read.csv("~/Desktop/Urban_tree_fire/landscape_analysis/model_inputs/eaton_burned_struc_model_inputs.csv") %>% 
  select(-X) %>% 
  filter(DAMAGE %in% c("Destroyed (>50%)", "No Damage", "Affected (1-9%)")) %>% 
  na.omit() %>% 
  left_join(eaton_wind, by = c("UID" = "point_UID"))

Eaton_struc <- Eaton_struc %>% mutate(angular_diff_build = pmin(abs(bearing_to_nearest_building - mean_wind_dir), 360 - abs(bearing_to_nearest_building - mean_wind_dir)), angular_diff_tree = pmin(abs(bearing_to_nearest_tree - mean_wind_dir), 360 - abs(bearing_to_nearest_tree - mean_wind_dir)), destroyed = if_else(DAMAGE %in% c("No Damage", "Affected (1-9%)"),0,1))

cols <- c("distance_to_nearest_tree", "distance_to_nearest_building", "angular_diff_tree", "number_of_trees_2m", "number_of_trees_5m", "number_of_trees_10m",  "area_of_trees_2m",  "area_of_trees_5m",  "area_of_trees_10m", "area_ext_build_2", "area_ext_build_10", "area_ext_build_25", "area_ext_build_50", "area_ext_build_100", "angular_diff_build", "build_dens", "tree_dens", "bl___10", "bl___25", "bl___50", "b___100")

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
                               build_dens + tree_dens  + IMRF(1 | utm_x + utm_y, model = spde), data = Eaton_struc_scale)

saveRDS(eaton_struc_m_scale, "~/Desktop/Urban_tree_fire/landscape_analysis/simple_models/eaton_struc_m1_scale.rda")

eaton_struc_m <- fitme(destroyed ~ distance_to_nearest_tree * angular_diff_tree + area_of_trees_2m + 
                         area_ext_build_2 + distance_to_nearest_building + angular_diff_build +   mean_wind_sp + 
                         build_dens + tree_dens  + IMRF(1 | utm_x + utm_y, model = spde), data = Eaton_struc)

saveRDS(eaton_struc_m, "~/Desktop/Urban_tree_fire/landscape_analysis/simple_models/eaton_struc_m1_unscaled.rda")


# Eaton Tree height



Eaton_tree <- read.csv("~/Desktop/Urban_tree_fire/landscape_analysis/model_inputs/eaton_burned_trees_model_inputs.csv") %>%
  replace_na(list(total_build_2m = 0, total_build_5m = 0, total_build_10m = 0)) %>% 
  na.omit()

Eaton_tree <- Eaton_tree %>% mutate(angular_diff_build = pmin(abs(bearing_to_closest_building - 45), 360 - abs(bearing_to_closest_building - 45)), angular_diff_tree = pmin(abs(bearing_to_nearest_tree - 45), 360 - abs(bearing_to_nearest_tree - 45)), burned = if_else(height_diff < 0 & ndvi_diff < -0.2, 1, 0), perc_diff = (height_diff/pre_height)*100) %>% 
  na.omit()

cols <- c("dist_to_nearest_tree", "dist_to_nearest_building", "angular_diff_tree", "number_of_trees_2m", "number_of_trees_5m", "number_of_trees_10m", "area_ext_tree2",  "area_ext_tree_5", "area_ext_tree_10", "total_build_2m", "total_build_5m", "total_build_10m", "angular_diff_build", "build_dens", "tree_dens")

names(cols) <- cols

Eaton_tree_scale <- mutate(Eaton_tree, across(all_of(cols), ~ scale(.x, center = F)))


spd <- sp::SpatialPointsDataFrame(
  coords = Eaton_tree[, c("utm_x", "utm_y")],
  data = Eaton_tree)

bound <- inla.nonconvex.hull(spd, convex = -.01, resolution = c(105,81))

mesh <- inla.mesh.2d(boundary = bound, max.edge = 200)

# Create SPDE model
spde <- INLA::inla.spde2.matern(mesh)


eaton_height_m_scaled <- fitme(perc_diff ~ dist_to_nearest_tree * angular_diff_tree + area_ext_tree2 + 
                                 total_build_2m + dist_to_nearest_building * angular_diff_build + 
                                 build_dens + tree_dens + IMRF(1 | utm_x + utm_y, model = spde), data = Eaton_tree_scale)

saveRDS(eaton_height_m_scaled, "~/Desktop/Urban_tree_fire/landscape_analysis/simple_models/eaton_tree_height_scaled.rda")


eaton_height_m <- fitme(perc_diff ~ dist_to_nearest_tree * angular_diff_tree + area_ext_tree2 + 
                                 total_build_2m + dist_to_nearest_building * angular_diff_build + 
                                 build_dens + tree_dens + IMRF(1 | utm_x + utm_y, model = spde), data = Eaton_tree)

saveRDS(eaton_height_m, "~/Desktop/Urban_tree_fire/landscape_analysis/simple_models/eaton_tree_height_unscaled.rda")


## Eaton Tree NDVI


Eaton_ndvi_m_scaled <- fitme(ndvi_diff ~ dist_to_nearest_tree * angular_diff_tree + area_ext_tree2 + 
                               total_build_2m + dist_to_nearest_building * angular_diff_build + 
                               build_dens + tree_dens + IMRF(1 | utm_x + utm_y, model = spde), data = Eaton_tree_scale)

saveRDS(Eaton_ndvi_m_scaled, "~/Desktop/Urban_tree_fire/landscape_analysis/simple_models/eaton_m2_scaled.rda")

Eaton_ndvi_m_unscaled <- fitme(ndvi_diff ~ dist_to_nearest_tree * angular_diff_tree + area_ext_tree2 + 
                               total_build_2m + dist_to_nearest_building * angular_diff_build + 
                               build_dens + tree_dens + IMRF(1 | utm_x + utm_y, model = spde), data = Eaton_tree)

saveRDS(Eaton_ndvi_m_unscaled, "~/Desktop/Urban_tree_fire/landscape_analysis/simple_models/eaton_m2_unscaled.rda")



# Eaton Tree unburned

Eaton_unburn <- read.csv("~/Desktop/Urban_tree_fire/landscape_analysis/model_inputs/eaton_unburned_trees_model_inputs.csv") %>% 
  na.omit()

Eaton_unburn <- Eaton_unburn %>% mutate(angular_diff_build = pmin(abs(bearing_to_closest_building - 45), 360 - abs(bearing_to_closest_building - 45)), angular_diff_tree = pmin(abs(bearing_to_nearest_tree - 45), 360 - abs(bearing_to_nearest_tree - 45)), burned = if_else(height_diff < 0 & ndvi_diff < -0.2, 1, 0), in_fire = "no", perc_diff = ((height_diff/pre_height)*100)) %>% 
  na.omit()

Eaton_unburn_clean <- Eaton_unburn %>%
  mutate(
    # if perc_diff came from a ratio, scrub bad values explicitly
    perc_diff = ifelse(is.finite(perc_diff), perc_diff, NA_real_)) %>% 
      na.omit()


cols <- c("dist_to_nearest_tree", "dist_to_nearest_building", "angular_diff_tree", "number_trees_2m", "number_of_trees_5m", "number_of_trees_10m", "area_ext_tree2",  "area_ext_tree_5", "area_ext_tree10", "total_build_2m", "total_build_5m", "total_build_10m", "angular_diff_build", "build_dens", "tree_dens")

names(cols) <- cols

Eaton_unburn_tree_scale <- mutate(Eaton_unburn_clean, across(all_of(cols), ~ scale(.x, center = F)))


spd <- sp::SpatialPointsDataFrame(
  coords = Eaton_unburn_clean[, c("utm_x", "utm_y")],
  data = Eaton_unburn_clean)

bound <- inla.nonconvex.hull(spd, convex = -.05, resolution = c(105,99))

mesh <- inla.mesh.2d(boundary = bound, max.edge = 250)

# Create SPDE model
spde <- INLA::inla.spde2.matern(mesh)

eaton_unburn_m <- fitme(perc_diff ~ 
             dist_to_nearest_tree * angular_diff_tree + 
             area_ext_tree2 + 
             total_build_2m + 
             dist_to_nearest_building * angular_diff_build + 
             build_dens + 
             tree_dens + 
             IMRF(1 | utm_x + utm_y, model = spde), 
           data = Eaton_unburn_clean)


saveRDS(eaton_unburn_m, "~/Desktop/Urban_tree_fire/landscape_analysis/simple_models/Eaton_tree_ext_m1_unscaled.rda")

eaton_unburn_m_scale <- fitme(perc_diff ~ 
                   dist_to_nearest_tree * angular_diff_tree + 
                   area_ext_tree2 + 
                   total_build_2m + 
                   dist_to_nearest_building * angular_diff_build + 
                   build_dens + 
                   tree_dens + 
                   IMRF(1 | utm_x + utm_y, model = spde), 
                 data = Eaton_unburn_tree_scale)


saveRDS(eaton_unburn_m_scale, "~/Desktop/Urban_tree_fire/landscape_analysis/simple_models/Eaton_tree_ext_m1_scaled.rda")

## Eaton tree unburned NDVI

eaton_unburn_m2 <- fitme(ndvi_diff ~ 
              dist_to_nearest_tree * angular_diff_tree + 
              area_ext_tree2 + 
              total_build_2m + 
              dist_to_nearest_building * angular_diff_build + 
              build_dens + 
              tree_dens + 
              IMRF(1 | utm_x + utm_y, model = spde), 
            data = Eaton_unburn)

saveRDS(eaton_unburn_m2, "~/Desktop/Urban_tree_fire/landscape_analysis/simple_models/Eaton_tree_ext_m2_unscaled.rda")

eaton_unburn_m2_scale <- fitme(ndvi_diff ~ 
                           dist_to_nearest_tree * angular_diff_tree + 
                           area_ext_tree2 + 
                           total_build_2m + 
                           dist_to_nearest_building * angular_diff_build + 
                           build_dens + 
                           tree_dens + 
                           IMRF(1 | utm_x + utm_y, model = spde), 
                         data = Eaton_unburn_tree_scale)

saveRDS(eaton_unburn_m2_scale, "~/Desktop/Urban_tree_fire/landscape_analysis/simple_models/Eaton_tree_ext_m2_scaled.rda")




# Palisades

# Structures


# Palisades Struc

pal_wind <- read.csv("wind_data/all_hrrr.csv") %>% 
  filter(point_fire == "palisades") %>% group_by(point_UID) %>% 
  summarize(mean_wind_sp = mean(wdsp), 
            mean_wind_dir = mean(wdir))



Palisades_struc <- read.csv("~/Desktop/Urban_tree_fire/landscape_analysis/model_inputs/palisades_burned_struc_model_inputs.csv") %>% 
  filter(DAMAGE %in% c("Destroyed (>50%)", "No Damage", "Affected (1-9%)")) %>% 
  left_join(pal_wind, by = c("UID" = "point_UID"))
# make new columns

Palisades_struc <- Palisades_struc %>% mutate(angular_diff_build = pmin(abs(bearing_to_nearest_building - mean_wind_dir), 360 - abs(bearing_to_nearest_building - mean_wind_dir)), angular_diff_tree = pmin(abs(bearing_to_nearest_tree - mean_wind_dir), 360 - abs(bearing_to_nearest_tree - mean_wind_dir)), destroyed = if_else(DAMAGE %in% c("No Damage", "Affected (1-9%)"),0,1)) %>% na.omit()


cols <- c("distance_to_nearest_tree", "distance_to_nearest_building", "angular_diff_tree", "number_of_trees_2m", "number_of_trees_5m", "number_of_trees_10m",  "area_of_trees_2m",  "area_of_trees_5m",  "area_of_trees_10m", "area_ext_build_2", "area_ext_build_10", "area_ext_build_25", "area_ext_build_50", "area_ext_build_100", "angular_diff_build", "build_dens", "tree_dens", "bl___10", "bl___25", "bl___50", "b___100")

names(cols) <- cols

Palisades_struc_scale <- mutate(Palisades_struc, across(all_of(cols), ~ scale(.x, center = F)))




spd <- sp::SpatialPointsDataFrame(
  coords = Palisades_struc[, c("utm_x", "utm_y")],
  data = Palisades_struc)

bound <- inla.nonconvex.hull(spd, convex = -.05)

mesh <- inla.mesh.2d(boundary = bound, max.edge = 200)

# Create SPDE model
spde <- INLA::inla.spde2.matern(mesh)

spaMM.options(separation_max = 100)

Palisades_struc_m_scaled <- fitme(destroyed ~ distance_to_nearest_tree + angular_diff_tree + area_of_trees_2m + area_ext_build_2 + distance_to_nearest_building + angular_diff_build + build_dens * tree_dens + mean_wind_sp + IMRF(1 | utm_x + utm_y, model = spde), data = Palisades_struc_scale)

saveRDS(Palisades_struc_m_scaled, "~/Desktop/Urban_tree_fire/landscape_analysis/simple_models/Palisades_struc_m1_scaled.rda")

Palisades_struc_m <- fitme(destroyed ~ distance_to_nearest_tree + angular_diff_tree + area_of_trees_2m + area_ext_build_2 + distance_to_nearest_building + angular_diff_build + build_dens * tree_dens + mean_wind_sp + IMRF(1 | utm_x + utm_y, model = spde), data = Palisades_struc)

saveRDS(Palisades_struc_m, "~/Desktop/Urban_tree_fire/landscape_analysis/simple_models/Palisades_struc_m1.rda")


# Palisades Tree height


Palisades_tree <- read.csv("~/Desktop/Urban_tree_fire/landscape_analysis/model_inputs/palisades_burned_trees_model_inputs.csv") %>%
  replace_na(list(total_build_2m = 0, total_build_5m = 0, total_build_10m = 0)) %>% 
  na.omit()

Palisades_tree <- Palisades_tree %>% mutate(angular_diff_build = pmin(abs(bearing_to_closest_building - 45), 360 - abs(bearing_to_closest_building - 45)), angular_diff_tree = pmin(abs(bearing_to_nearest_tree - 45), 360 - abs(bearing_to_nearest_tree - 45)), burned = if_else(height_diff < 0 & ndvi_diff < -0.2, 1, 0), perc_diff = (height_diff/pre_height)*100) %>% 
  na.omit()

cols <- c("dist_to_nearest_tree", "angular_diff_tree", "number_of_trees_2m", "number_of_trees_5m",   "number_of_trees_10m", "total_build_2m", "total_build_5m", "total_build_10m", "area_ext_tree2",  "area_ext_tree_5",  "area_ext_tree_10", "dist_to_nearest_building", "angular_diff_build", "build_dens", "tree_dens")

names(cols) <- cols

Palisades_tree_scale <- mutate(Palisades_tree, across(all_of(cols), scale))


spd <- sp::SpatialPointsDataFrame(
  coords = Palisades_tree[, c("utm_x", "utm_y")],
  data = Palisades_tree)

bound <- inla.nonconvex.hull(spd, convex = -.01, resolution = c(105,86))

mesh <- inla.mesh.2d(boundary = bound, max.edge = 200)

# Create SPDE model
spde <- INLA::inla.spde2.matern(mesh)

pal_height_m_scaled <- fitme(perc_diff ~ dist_to_nearest_tree * angular_diff_tree + area_ext_tree2 + 
                               total_build_2m + dist_to_nearest_building * angular_diff_build + 
                               build_dens + tree_dens + IMRF(1 | utm_x + utm_y, model = spde), data = Palisades_tree_scale)

saveRDS(pal_height_m_scaled, "~/Desktop/Urban_tree_fire/landscape_analysis/simple_models/Palisades_tree_m1_scaled.rda")


pal_height_m <- fitme(perc_diff ~ dist_to_nearest_tree * angular_diff_tree + area_ext_tree2 + 
                               total_build_2m + dist_to_nearest_building * angular_diff_build + 
                               build_dens + tree_dens + IMRF(1 | utm_x + utm_y, model = spde), data = Palisades_tree)

saveRDS(pal_height_m, "~/Desktop/Urban_tree_fire/landscape_analysis/simple_models/Palisades_tree_m1_unscaled.rda")


# Palisades Tree NDVI


Palisades_ndvi_m_scaled <- fitme(ndvi_diff ~ dist_to_nearest_tree * angular_diff_tree + area_ext_tree2 + 
                                   total_build_2m + dist_to_nearest_building * angular_diff_build + 
                                   build_dens + tree_dens + IMRF(1 | utm_x + utm_y, model = spde), data = Palisades_tree_scale)

saveRDS(Palisades_ndvi_m_scaled, "~/Desktop/Urban_tree_fire/landscape_analysis/simple_models/Palisades_tree_m2_scaled.rda")


Palisades_ndvi_m <- fitme(ndvi_diff ~ dist_to_nearest_tree * angular_diff_tree + area_ext_tree2 + 
                                   total_build_2m + dist_to_nearest_building * angular_diff_build + 
                                   build_dens + tree_dens + IMRF(1 | utm_x + utm_y, model = spde), data = Palisades_tree)

saveRDS(Palisades_ndvi_m_scaled, "~/Desktop/Urban_tree_fire/landscape_analysis/simple_models/Palisades_tree_m2.rda")



 ## Palisades External


Palisades_unburn <- read.csv("~/Desktop/Urban_tree_fire/landscape_analysis/model_inputs/palisades_unburned_trees_model_inputs.csv") %>% 
  mutate(angular_diff_build = pmin(abs(bearing_to_closest_building - 45), 360 - abs(bearing_to_closest_building - 45)), angular_diff_tree = pmin(abs(bearing_to_nearest_tree - 45), 360 - abs(bearing_to_nearest_tree - 45)), burned = if_else(height_diff < 0 & ndvi_diff < -0.2, 1, 0), perc_diff = (height_diff/pre_height)*100) %>% na.omit()

cols <- c("dist_to_nearest_tree", "angular_diff_tree", "number_trees_2m", "number_of_trees_5m",   "number_of_trees_10m", "total_build_2m", "total_build_5m", "total_build_10m", "area_ext_tree2",  "area_ext_tree_5",  "area_ext_tree10", "dist_to_nearest_building", "angular_diff_build", "build_dens", "tree_dens")

names(cols) <- cols

Palisades_unburn_tree_scale <- mutate(Palisades_unburn, across(all_of(cols), ~ scale(.x, center = F)))

# generate mesh

spd <- sp::SpatialPointsDataFrame(
  coords = Palisades_unburn[, c("utm_x", "utm_y")],
  data = Palisades_unburn)

bound <- inla.nonconvex.hull(spd, convex = -.05)

mesh <- inla.mesh.2d(boundary = bound, max.edge = 200)

# Create SPDE model
spde <- INLA::inla.spde2.matern(mesh)




pal_m_scale <- fitme(perc_diff ~ 
                       dist_to_nearest_tree * angular_diff_tree + 
                       area_ext_tree2 + 
                       total_build_2m + 
                       dist_to_nearest_building * angular_diff_build + 
                       build_dens + 
                       tree_dens + 
                       IMRF(1 | utm_x + utm_y, model = spde), 
                     data = Palisades_unburn_tree_scale)





saveRDS(pal_m_scale, "~/Desktop/Urban_tree_fire/landscape_analysis/simple_models/Palisades_tree_m1_ext_scaled.rda")

pal_m <- fitme(perc_diff ~ 
                 dist_to_nearest_tree * angular_diff_tree + 
                 area_ext_tree2 + 
                 total_build_2m + 
                 dist_to_nearest_building * angular_diff_build + 
                 build_dens + 
                 tree_dens + 
                 IMRF(1 | utm_x + utm_y, model = spde), 
               data = Palisades_unburn)

saveRDS(pal_m, "~/Desktop/Urban_tree_fire/landscape_analysis/simple_models/Palisades_tree_m1_ext_unscaled.rda")


pal_m2_scale <- fitme(ndvi_diff ~ dist_to_nearest_tree * angular_diff_tree + 
                        area_ext_tree2 + 
                        total_build_2m + 
                        dist_to_nearest_building * 
                        angular_diff_build + 
                        build_dens + 
                        tree_dens + 
                        IMRF(1 | utm_x + utm_y, model = spde), 
                      data = Palisades_unburn_tree_scale)



saveRDS(pal_m2_scale, "~/Desktop/Urban_tree_fire/landscape_analysis/simple_models/Palisades_tree_m2_ext_scaled.rda")

pal_m2 <- fitme(ndvi_diff ~ dist_to_nearest_tree * angular_diff_tree + 
                  area_ext_tree2 + 
                  total_build_2m + 
                  dist_to_nearest_building * 
                  angular_diff_build + 
                  build_dens + 
                  tree_dens + 
                  IMRF(1 | utm_x + utm_y, model = spde), 
                data = Palisades_unburn)

saveRDS(pal_m2, "~/Desktop/Urban_tree_fire/landscape_analysis/simple_models/Palisades_tree_m2_ext_unscaled.rda")





