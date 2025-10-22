# 1. Load libraries
library(sf)       # for vector spatial operations
library(dplyr)    # for data manipulation
library(tidyr)    # for pivoting to a matrix
library(terra)
library(caret)



## Eaton


study_area_eaton <- st_read("~/OneDrive - Cal Poly/Advanced_GIS_data/Random_plots/eaton_rand_plots_1.shp")


# 1. read your tree‐canopy layers, rename for clarity
truth_e <- st_read("~/OneDrive - Cal Poly/Analysis_Results/ArcExports/Eaton_training_final.shp") %>%
  st_union() %>%                           # merge all into one multipolygon
  st_sf(truth = TRUE, geometry = .)

truth_e_rast_blank <- rast(truth_e, ncols = st_bbox(truth_e)[3] - st_bbox(truth_e)[1], nrows = st_bbox(truth_e)[4] - st_bbox(truth_e)[2])

truth_e_rast <- terra::rasterize(truth_e, truth_e_rast_blank) 

truth_e_rast[is.na(truth_e_rast[])] <- 0 

truth_e_rast <- terra::crop(truth_e_rast, study_area_eaton, mask = T)

pred_e <- st_read("~/Desktop/Urban_tree_fire/landscape_analysis/tmp_data/eaton_burned_tree_crowns_final.shp")

pred_e_plots <- st_join(pred_e, study_area_eaton, left = F) %>%
  st_union() %>%                           # merge all into one multipolygon
  st_sf(pred_e = TRUE, geometry = .)


pred_e_plots_rast <- terra::rasterize(pred_e_plots, truth_e_rast_blank)

pred_e_plots_rast <- terra::crop(pred_e_plots_rast, truth_e_rast)

#pred_e_plots_rast <- terra::resample(pred_e_plots_rast, truth_e_rast)

pred_e_plots_rast[is.na(pred_e_plots_rast[])] <- 0 

pred_e_plots_rast <- terra::crop(pred_e_plots_rast, study_area_eaton, mask = T)




v_truth_e <- values(truth_e_rast)
v_pred_e  <- values(pred_e_plots_rast)

# 4. drop any cells where either is NA
ok <- !is.na(v_truth_e) & !is.na(v_pred_e)
v_truth_e <- v_truth_e[ok]
v_pred_e  <- v_pred_e[ok]


y_true_e <- factor(v_truth_e, levels = c(0,1), labels = c("non-tree","tree"))
y_pred_e <- factor(v_pred_e,  levels = c(0,1), labels = c("non-tree","tree"))

cm_e <- confusionMatrix(y_pred_e, y_true_e, positive = "tree")
print(cm_e)



## Palisades

study_area_pal <- st_read("~/OneDrive - Cal Poly/Advanced_GIS_data/Random_plots/palisades_rand_plots.shp")

ptruth <- st_read("~/OneDrive - Cal Poly/Analysis_Results/ArcExports/Palisades_training_final.shp") %>%
  st_union() %>%                           # merge all into one multipolygon
  st_sf(truth = TRUE, geometry = .) %>% 
  st_transform(32611)

ptruth_rast_blank <- rast(ptruth, ncols = st_bbox(ptruth)[3] - st_bbox(ptruth)[1], nrows = st_bbox(ptruth)[4] - st_bbox(ptruth)[2])

ptruth_rast <- terra::rasterize(ptruth, ptruth_rast_blank)

ptruth_rast[is.na(ptruth_rast[])] <- 0 

ptruth_rast <- terra::crop(ptruth_rast, study_area_pal, mask = T)

ppred  <- st_read("~/Desktop/Urban_tree_fire/landscape_analysis/tmp_data/palisades_burned_tree_crowns_final.shp")

ppred_plots <- st_join(ppred, study_area_pal, left = F) %>% 
  st_union() %>%                           # merge all into one multipolygon
  st_sf(ppred = TRUE, geometry = .)

ppred_plots_rast <- terra::rasterize(ppred_plots, ptruth_rast_blank)

ppred_plots_rast[is.na(ppred_plots_rast[])] <- 0 

ppred_plots_rast <- terra::crop(ppred_plots_rast, study_area_pal, mask = T)

v_truth_p <- values(ptruth_rast)
v_pred_p  <- values(ppred_plots_rast)

# 4. drop any cells where either is NA
ok <- !is.na(v_truth_p) & !is.na(v_pred_p)
v_truth_p <- v_truth_p[ok]
v_pred_p  <- v_pred_p[ok]


y_true_p <- factor(v_truth_p, levels = c(0,1), labels = c("non-tree","tree"))
y_pred_p <- factor(v_pred_p,  levels = c(0,1), labels = c("non-tree","tree"))

cm_p <- confusionMatrix(y_pred_p, y_true_p, positive = "tree")
print(cm_p)




