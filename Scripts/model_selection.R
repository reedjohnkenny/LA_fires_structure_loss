

library(spaMM)
library(dplyr)
library(INLA)
library(tidyr)

source("~/Desktop/Urban_tree_fire/structure_analysis/Scripts/functions/compare_AIC.R")
source("functions/plot_alt_models.R")

## Select best model for Eaton struc

## Eaton struc

setwd("~/Desktop/Urban_tree_fire/structure_analysis/Scripts/")

eaton_wind <- read.csv("../wind_data/all_hrrr.csv") %>% 
  filter(point_fire == "eaton") %>% group_by(point_UID) %>% 
  summarize(mean_wind_sp = mean(wdsp), 
            mean_wind_dir = mean(wdir))

eaton_areas <- read.csv("../model_inputs/eaton_burned_alt_inputs.csv")


Eaton_struc <- read.csv("../model_inputs/eaton_burned_struc_model_inputs.csv") %>% 
  select(-X) %>% 
  filter(DAMAGE %in% c("Destroyed (>50%)", "No Damage", "Affected (1-9%)")) %>% 
  na.omit() %>% 
  left_join(eaton_wind, by = c("UID" = "point_UID")) %>% 
  left_join(eaton_areas, by = "UID")

Eaton_struc <- Eaton_struc %>% mutate(angular_diff_build = pmin(abs(bearing_to_nearest_building - mean_wind_dir), 360 - abs(bearing_to_nearest_building - mean_wind_dir)), angular_diff_tree = pmin(abs(bearing_to_nearest_tree - mean_wind_dir), 360 - abs(bearing_to_nearest_tree - mean_wind_dir)), destroyed = if_else(DAMAGE %in% c("No Damage", "Affected (1-9%)"),0,1))

cols <- c("distance_to_nearest_tree", "distance_to_nearest_building", "angular_diff_tree", "number_of_trees_2m", "number_of_trees_5m", "number_of_trees_10m",  "area_of_trees_2m",  "area_of_trees_5m",  "area_of_trees_10m", "area_ext_build_2", "area_ext_build_10", "area_ext_build_25", "area_ext_build_50", "area_ext_build_100", "angular_diff_build", "build_dens", "tree_dens", "mean_perc_struc", "area_tree_200m", "total_build_200_m", "fprint_area", "num_trees_200m", "mean_tree_density", "mean_struc_density", "num_builds_200m", "tree_dens_500", "tree_dens_300", "area_tree_100m", "area_tree_50m", "total_build_100_m", "build_dens_500", "build_dens_300", "total_build_50_m")

names(cols) <- cols

Eaton_struc_scale <- mutate(Eaton_struc, across(all_of(cols), ~ scale(.x, center = T)))

spd <- sp::SpatialPointsDataFrame(
  coords = Eaton_struc[, c("utm_x", "utm_y")],
  data = Eaton_struc)

bound <- inla.nonconvex.hull(spd, convex = -.01, resolution = c(105,82))

mesh <- inla.mesh.2d(boundary = bound, max.edge = 200)

# Create SPDE model
spde <- INLA::inla.spde2.matern(mesh)

spaMM.options(separation_max = 100) # or even higher if suggested



# 1. Define the sets of predictors
tree_dens_vars <- c("area_tree_50m", "area_tree_100m", "area_tree_200m")
build_dens_vars  <- c("area_ext_build_25","area_ext_build_50", "area_ext_build_100", "total_build_50_m", "total_build_100_m", "total_build_200_m")

# 4. Make every combination of (tree_var × build_var)
combos <- expand.grid(
  tree_var  = tree_dens_vars,
  build_var = build_dens_vars,
  stringsAsFactors = FALSE
)

# 5. Pre‐allocate a named list to hold each fit
Eaton_model_list <- vector("list", nrow(combos))
names(Eaton_model_list) <- paste0("m", seq_len(nrow(combos)))

# 6. Loop over each row of combos, build the formula, and call fitme()
for(i in seq_len(nrow(combos))) {
  fmla <- as.formula(
    paste0(
      "destroyed ~ ",
      "distance_to_nearest_tree + angular_diff_tree + ",
      combos$tree_var[i], " + ",
      combos$build_var[i], " + ",
      "distance_to_nearest_building + angular_diff_build + ",
      "area_ext_build_2 + ",
      "area_of_trees_2m + ",
      "IMRF(1 | utm_x + utm_y, model = spde)"
    )
  )
  Eaton_model_list[[i]] <- fitme(fmla, data = Eaton_struc_scale, family = binomial(link = "logit"))
}

# 5.1 Pre‐allocate a named list to hold each fit
Eaton_int_model_list <- vector("list", nrow(combos))
names(Eaton_int_model_list) <- paste0("m_int", seq_len(nrow(combos)))

# 6. Loop over each row of combos, build the formula, and call fitme()
for(i in seq_len(nrow(combos))) {
  fmla <- as.formula(
    paste0(
      "destroyed ~ ",
      "distance_to_nearest_tree * angular_diff_tree + ",
      combos$tree_var[i], " * ",
      combos$build_var[i], " + ",
      "distance_to_nearest_building * angular_diff_build + ",
      "area_ext_build_2 + ",
      "area_of_trees_2m + ",
      "IMRF(1 | utm_x + utm_y, model = spde)"
    )
  )
  Eaton_int_model_list[[i]] <- fitme(fmla, data = Eaton_struc_scale, family = binomial(link = "logit"))
}

eaton_int_1.5 <- fitme(destroyed ~ distance_to_nearest_tree + angular_diff_tree + tree_dens * 
                         build_dens + distance_to_nearest_building + angular_diff_build + 
                         area_ext_build_2 + area_of_trees_2m +IMRF(1 | utm_x + utm_y, 
                                                                    model = spde), 
                       data = Eaton_struc_scale, family = binomial(link = "logit"))

Eaton_int_model_list$int_1.5 <- eaton_int_1.5


Eaton_full_model_list <- c(Eaton_model_list, Eaton_int_model_list)

eaton_struc_AICs <- compare_AIC(models = Eaton_full_model_list)


saveRDS(eaton_int_1.5, "~/Desktop/Urban_tree_fire/structure_analysis/models/eaton_struc_m1_scale.rda")

eaton_struc_m_unscale <- fitme(destroyed ~ distance_to_nearest_tree + angular_diff_tree +  tree_dens * 
                                 build_dens + distance_to_nearest_building + angular_diff_build + 
                                 area_ext_build_2 + area_of_trees_2m + IMRF(1 | utm_x + utm_y, 
                                                                            model = spde), 
                               data = Eaton_struc, family = binomial(link = "logit"))

saveRDS(eaton_struc_m_unscale, "../models/eaton_struc_m1_unscaled.rda")


## Palisades Models


# Palisades Struc

pal_wind <- read.csv("../wind_data/all_hrrr.csv") %>% 
  filter(point_fire == "palisades") %>% group_by(point_UID) %>% 
  summarize(mean_wind_sp = mean(wdsp), 
            mean_wind_dir = mean(wdir))

pal_areas <- read.csv("../model_inputs/palisades_burned_alt_inputs.csv")


Palisades_struc <- read.csv("../model_inputs/palisades_burned_struc_model_inputs.csv") %>% 
  filter(DAMAGE %in% c("Destroyed (>50%)", "No Damage", "Affected (1-9%)")) %>% 
  left_join(pal_wind, by = c("UID" = "point_UID")) %>% 
  left_join(pal_areas, by = "UID")
# make new columns

Palisades_struc <- Palisades_struc %>% mutate(angular_diff_build = pmin(abs(bearing_to_nearest_building - mean_wind_dir), 360 - abs(bearing_to_nearest_building - mean_wind_dir)), angular_diff_tree = pmin(abs(bearing_to_nearest_tree - mean_wind_dir), 360 - abs(bearing_to_nearest_tree - mean_wind_dir)), destroyed = if_else(DAMAGE %in% c("No Damage", "Affected (1-9%)"),0,1))


cols <- c("distance_to_nearest_tree", "distance_to_nearest_building", "angular_diff_tree", "number_of_trees_2m", "number_of_trees_5m", "number_of_trees_10m",  "area_of_trees_2m",  "area_of_trees_5m",  "area_of_trees_10m", "area_ext_build_2", "area_ext_build_10", "area_ext_build_25", "area_ext_build_50", "area_ext_build_100", "angular_diff_build", "build_dens", "tree_dens", "mean_perc_struc", "area_tree_200m", "total_build_200_m", "num_trees_200m", "mean_tree_density", "mean_struc_density", "num_builds_200m", "tree_dens_500", "tree_dens_300", "area_tree_100m", "area_tree_50m", "total_build_50_m")

names(cols) <- cols

Palisades_struc_scale <- mutate(Palisades_struc, across(all_of(cols), ~ scale(.x, center = T)))




spd <- sp::SpatialPointsDataFrame(
  coords = Palisades_struc[, c("utm_x", "utm_y")],
  data = Palisades_struc)

bound <- inla.nonconvex.hull(spd, convex = -.05)

mesh <- inla.mesh.2d(boundary = bound, max.edge = 250)

# Create SPDE model
spde <- INLA::inla.spde2.matern(mesh)

spaMM.options(separation_max = 100)

# 1. Define the sets of predictors
tree_dens_vars <- c("area_tree_50m")
build_dens_vars  <- c("area_ext_build_25","area_ext_build_50", "area_ext_build_100", "total_build_50_m", "total_build_100_m", "total_build_200_m")


# 4. Make every combination of (tree_var × build_var)
combos <- expand.grid(
  tree_var  = tree_dens_vars,
  build_var = build_dens_vars,
  stringsAsFactors = FALSE
)

# 5. Pre‐allocate a named list to hold each fit
pal_struc_model_list <- vector("list", nrow(combos))
names(pal_struc_model_list) <- paste0("m", seq_len(nrow(combos)))

# 6. Loop over each row of combos, build the formula, and call fitme()
for(i in seq_len(nrow(combos))) {
  fmla <- as.formula(
    paste0(
      "destroyed ~ ",
      "distance_to_nearest_tree + angular_diff_tree + ",
      combos$tree_var[i], " + ",
      combos$build_var[i], " + ",
      "distance_to_nearest_building + angular_diff_build + ",
      "area_ext_build_2 + ",
      "area_of_trees_2m + ",
      "IMRF(1 | utm_x + utm_y, model = spde)"
    )
  )
  pal_struc_model_list[[i]] <- fitme(fmla, data = Palisades_struc_scale, family = binomial(link = "logit"))
}

# 5. Pre‐allocate a named list to hold each fit
pal_struc_model_int_list <- vector("list", nrow(combos))
names(pal_struc_model_int_list) <- paste0("m_int", seq_len(nrow(combos)))

# 6. Loop over each row of combos, build the formula, and call fitme()
for(i in seq_len(nrow(combos))) {
  fmla <- as.formula(
    paste0(
      "destroyed ~ ",
      "distance_to_nearest_tree * angular_diff_tree + ",
      combos$tree_var[i], " * ",
      combos$build_var[i], " + ",
      "distance_to_nearest_building * angular_diff_build + ",
      "area_ext_build_2 + ",
      "area_of_trees_2m + ",
      "IMRF(1 | utm_x + utm_y, model = spde)"
    )
  )
  pal_struc_model_int_list[[i]] <- fitme(fmla, data = Palisades_struc_scale, family = binomial(link = "logit"))
}

pal_struc_model_int_4.5 <- fitme(destroyed ~ distance_to_nearest_tree + angular_diff_tree + tree_dens * 
                                   build_dens + distance_to_nearest_building + angular_diff_build + 
                                   area_ext_build_2 + area_of_trees_2m +IMRF(1 | utm_x + utm_y, 
                                                                              model = spde), 
                                 data = Palisades_struc_scale, family = binomial(link = "logit"))


pal_struc_model_list$int_4.5 <- pal_struc_model_int_4.5

pal_struc_full_model_list <- c(pal_struc_model_list, pal_struc_model_int_list)



Palisades_struc_AICs <- compare_AIC(models = pal_struc_full_model_list)

saveRDS(pal_struc_full_model_list$int_4.5, "../models/Palisades_struc_m1_scaled.rda")

Palisades_struc_m_unscaled <- fitme(destroyed ~ distance_to_nearest_tree + angular_diff_tree + tree_dens * 
                                      build_dens + distance_to_nearest_building + angular_diff_build + 
                                      area_ext_build_2 + area_of_trees_2m + IMRF(1 | utm_x + utm_y, 
                                                                                 model = spde), 
                                    data = Palisades_struc, family = binomial(link = "logit"))

saveRDS(Palisades_struc_m_unscaled, "../models/Palisades_struc_m1_unscaled.rda")

