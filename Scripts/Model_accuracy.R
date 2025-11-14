library(spaMM)
library(dplyr)
library(INLA)
library(tidyr)
library(pROC)

source("~/Desktop/Urban_tree_fire/structure_analysis/Scripts/functions/compare_AIC.R")
source("functions/plot_alt_models.R")



## Generate training and testing data

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

Eaton_struc_train <- sample_frac(Eaton_struc, size = 0.5)

cols <- c("distance_to_nearest_tree", "distance_to_nearest_building", "angular_diff_tree", "number_of_trees_2m", "number_of_trees_5m", "number_of_trees_10m",  "area_of_trees_2m",  "area_of_trees_5m",  "area_of_trees_10m", "area_ext_build_2", "area_ext_build_10", "area_ext_build_25", "area_ext_build_50", "area_ext_build_100", "angular_diff_build", "build_dens", "tree_dens", "mean_perc_struc", "area_tree_200m", "total_build_200_m", "fprint_area", "num_trees_200m", "mean_tree_density", "mean_struc_density", "num_builds_200m", "tree_dens_500", "tree_dens_1000", "area_tree_100m", "area_tree_50m", "total_build_100_m", "build_dens_500", "build_dens_1000", "total_build_50_m")

names(cols) <- cols

Eaton_struc_train_scale <- mutate(Eaton_struc_train, across(all_of(cols), ~ scale(.x, center = F)))

spd <- sp::SpatialPointsDataFrame(
  coords = Eaton_struc_train[, c("utm_x", "utm_y")],
  data = Eaton_struc_train)

bound <- inla.nonconvex.hull(spd, convex = -.01, resolution = c(105,82))

mesh <- inla.mesh.2d(boundary = bound, max.edge = 200)

# Create SPDE model
spde <- INLA::inla.spde2.matern(mesh)

spaMM.options(separation_max = 100) # or even higher if suggested


eaton_train_m <- fitme(destroyed ~ distance_to_nearest_tree + angular_diff_tree + tree_dens * 
                         build_dens + distance_to_nearest_building + angular_diff_build + 
                         area_ext_build_2 + area_of_trees_2m +IMRF(1 | utm_x + utm_y, 
                                                                   model = spde), 
                       data = Eaton_struc_train_scale, family = binomial(link = "logit"))


## ----- Holdout/test split -----------------------------------------------------
Eaton_test <- Eaton_struc %>% filter(!UID %in% Eaton_struc_train$UID)

## Build a scaled train frame (already done above) and a scaled test frame
Eaton_test_scale <- mutate(Eaton_test, across(all_of(cols), ~ scale(.x, center = F)))

## ----- Predict “destroyed” on the test set -----------------------------------
## For out-of-sample, use population-level prediction (exclude spatial RE):
## re.form = ~0 avoids borrowing spatial random effects estimated on training set.
test_prob <- predict(eaton_train_m,
                     newdata = Eaton_test_scale,
                     type = "response",
                     re.form = ~0)

Eaton_test_pred <- Eaton_test_scale %>%
  mutate(p_hat = as.numeric(test_prob),
         pred_0.5 = as.integer(p_hat >= 0.5))

## ----- Confusion matrix at 0.5 threshold -------------------------------------
cm_05 <- with(Eaton_test_pred, table(
  truth = destroyed,
  pred  = pred_0.5
))

accuracy_05   <- sum(diag(cm_05)) / sum(cm_05)
sensitivity_05 <- ifelse(sum(cm_05[2, ]) > 0, cm_05["1","1"] / sum(cm_05[2, ]), NA) # TPR
specificity_05 <- ifelse(sum(cm_05[1, ]) > 0, cm_05["0","0"] / sum(cm_05[1, ]), NA) # TNR

cat("Confusion matrix (threshold = 0.5):\n"); print(cm_05)
cat(sprintf("Accuracy: %.3f  |  Sensitivity (TPR): %.3f  |  Specificity (TNR): %.3f\n",
            accuracy_05, sensitivity_05, specificity_05))

## ----- ROC / AUC and an ROC-optimized threshold ------------------------------
roc_obj <- roc(response = Eaton_test_pred$destroyed, predictor = Eaton_test_pred$p_hat,
               quiet = TRUE, direction = "<")  # destroyed=1 should correspond to higher p_hat
auc_val <- as.numeric(auc(roc_obj))

## Youden-optimal threshold
opt_coords <- coords(roc_obj, x = "best", best.method = "youden", transpose = TRUE)
thr_best <- unname(opt_coords["threshold"])

Eaton_test_pred <- Eaton_test_pred %>%
  mutate(pred_best = as.integer(p_hat >= thr_best))

cm_best <- with(Eaton_test_pred, table(
  truth = destroyed,
  pred  = pred_best
))

accuracy_best    <- sum(diag(cm_best)) / sum(cm_best)
sensitivity_best <- ifelse(sum(cm_best[2, ]) > 0, cm_best["1","1"] / sum(cm_best[2, ]), NA)
specificity_best <- ifelse(sum(cm_best[1, ]) > 0, cm_best["0","0"] / sum(cm_best[1, ]), NA)

cat(sprintf("\nAUC: %.3f\n", auc_val))
cat(sprintf("Youden-optimal threshold: %.3f\n", thr_best))
cat("Confusion matrix (optimal threshold):\n"); print(cm_best)
cat(sprintf("Accuracy: %.3f  |  Sensitivity (TPR): %.3f  |  Specificity (TNR): %.3f\n",
            accuracy_best, sensitivity_best, specificity_best))

## ----- (Optional) compact results tibble -------------------------------------
perf_summary <- tibble(
  threshold = c(0.5, thr_best),
  accuracy  = c(accuracy_05, accuracy_best),
  sensitivity = c(sensitivity_05, sensitivity_best),
  specificity = c(specificity_05, specificity_best),
  AUC = auc_val
)
print(perf_summary)


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

Palisades_struc_train <- sample_frac(Palisades_struc, size = 0.5)

cols <- c("distance_to_nearest_tree", "distance_to_nearest_building", "angular_diff_tree", "number_of_trees_2m", "number_of_trees_5m", "number_of_trees_10m",  "area_of_trees_2m",  "area_of_trees_5m",  "area_of_trees_10m", "area_ext_build_2", "area_ext_build_10", "area_ext_build_25", "area_ext_build_50", "area_ext_build_100", "angular_diff_build", "build_dens", "tree_dens", "mean_perc_struc", "area_tree_200m", "total_build_200_m", "num_trees_200m", "mean_tree_density", "mean_struc_density", "num_builds_200m", "tree_dens_500", "tree_dens_1000", "area_tree_100m", "area_tree_50m", "total_build_50_m")

names(cols) <- cols

Palisades_struc_train_scale <- mutate(Palisades_struc_train, across(all_of(cols), ~ scale(.x, center = F)))



spd <- sp::SpatialPointsDataFrame(
  coords = Palisades_struc_train[, c("utm_x", "utm_y")],
  data = Palisades_struc_train)

bound <- inla.nonconvex.hull(spd, convex = -.05)

mesh <- inla.mesh.2d(boundary = bound, max.edge = 250)

# Create SPDE model
spde <- INLA::inla.spde2.matern(mesh)

spaMM.options(separation_max = 100)

pal_train_m <- fitme(destroyed ~ distance_to_nearest_tree + angular_diff_tree + tree_dens * 
                                   build_dens + distance_to_nearest_building + angular_diff_build + 
                                   area_ext_build_2 + area_of_trees_2m +IMRF(1 | utm_x + utm_y, 
                                                                             model = spde), 
                                 data = Palisades_struc_train_scale, family = binomial(link = "logit"))


# ----- Holdout/test split -----------------------------------------------------
Palisades_test <- Palisades_struc %>% filter(!UID %in% Palisades_struc_train$UID)

## Use the SAME scaling as training (training used scale(center=FALSE))
## Get training SDs for each covariate (no centering was applied)
train_sds <- sapply(Palisades_struc_train[ , cols], sd, na.rm = TRUE)
# Guard against zeros
train_sds[train_sds == 0 | is.na(train_sds)] <- 1

## Build a scaled train frame (already done above) and a scaled test frame
Palisades_test_scale <- mutate(Palisades_test, across(all_of(cols), ~ scale(.x, center = F))) %>% 
  na.omit

## ----- Predict “destroyed” on the test set -----------------------------------
## For out-of-sample, use population-level prediction (exclude spatial RE):
## re.form = ~0 avoids borrowing spatial random effects estimated on training set.
test_prob <- predict(pal_train_m,
                     newdata = Palisades_test_scale,
                     type = "response",
                     re.form = ~0)


Palisades_test_pred <- Palisades_test_scale %>%
  mutate(p_hat = as.numeric(test_prob),
         pred_0.5 = as.integer(p_hat >= 0.5))

## ----- Confusion matrix at 0.5 threshold -------------------------------------
cm_05 <- with(Palisades_test_pred, table(
  truth = destroyed,
  pred  = pred_0.5
))

accuracy_05   <- sum(diag(cm_05)) / sum(cm_05)
sensitivity_05 <- ifelse(sum(cm_05[2, ]) > 0, cm_05["1","1"] / sum(cm_05[2, ]), NA) # TPR
specificity_05 <- ifelse(sum(cm_05[1, ]) > 0, cm_05["0","0"] / sum(cm_05[1, ]), NA) # TNR

cat("Confusion matrix (threshold = 0.5):\n"); print(cm_05)
cat(sprintf("Accuracy: %.3f  |  Sensitivity (TPR): %.3f  |  Specificity (TNR): %.3f\n",
            accuracy_05, sensitivity_05, specificity_05))

## ----- ROC / AUC and an ROC-optimized threshold ------------------------------
roc_obj <- roc(response = Palisades_test_pred$destroyed, predictor = Palisades_test_pred$p_hat,
               quiet = TRUE, direction = "<")  # destroyed=1 should correspond to higher p_hat
auc_val <- as.numeric(auc(roc_obj))

## Youden-optimal threshold
opt_coords <- coords(roc_obj, x = "best", best.method = "youden", transpose = TRUE)
thr_best <- unname(opt_coords["threshold"])

Palisades_test_pred <- Palisades_test_pred %>%
  mutate(pred_best = as.integer(p_hat >= thr_best))

cm_best <- with(Palisades_test_pred, table(
  truth = destroyed,
  pred  = pred_best
))

accuracy_best    <- sum(diag(cm_best)) / sum(cm_best)
sensitivity_best <- ifelse(sum(cm_best[2, ]) > 0, cm_best["1","1"] / sum(cm_best[2, ]), NA)
specificity_best <- ifelse(sum(cm_best[1, ]) > 0, cm_best["0","0"] / sum(cm_best[1, ]), NA)

cat(sprintf("\nAUC: %.3f\n", auc_val))
cat(sprintf("Youden-optimal threshold: %.3f\n", thr_best))
cat("Confusion matrix (optimal threshold):\n"); print(cm_best)
cat(sprintf("Accuracy: %.3f  |  Sensitivity (TPR): %.3f  |  Specificity (TNR): %.3f\n",
            accuracy_best, sensitivity_best, specificity_best))

## ----- (Optional) compact results tibble -------------------------------------
perf_summary <- tibble(
  threshold = c(0.5, thr_best),
  accuracy  = c(accuracy_05, accuracy_best),
  sensitivity = c(sensitivity_05, sensitivity_best),
  specificity = c(specificity_05, specificity_best),
  AUC = auc_val
)
print(perf_summary)


