library(caret)
library(tidyverse)

setwd("~/Desktop/Urban_tree_fire/structure_analysis/validation_data/")


# Palisades


val_points <- read.csv("palisades_validation_points.csv") %>% 
  drop_na() %>% 
  mutate(map_bin = if_else(mppd_cl == "canopy", 1, 0))

y_true_p <- factor(val_points$map_bin, levels = c(0,1), labels = c("non_canopy","canopy"))
y_pred_p <- factor(val_points$ref_class,  levels = c(0,1), labels = c("non_canopy","canopy"))

cm_p <- confusionMatrix(y_pred_p, y_true_p, positive = "canopy")
print(cm_p)


# Eaton

val_points <- read.csv("eaton_validation_points.csv") %>% 
  drop_na() %>% 
  mutate(map_bin = if_else(mppd_cl == "canopy", 1, 0))

y_true_p <- factor(val_points$map_bin, levels = c(0,1), labels = c("non_canopy","canopy"))
y_pred_p <- factor(val_points$ref_class,  levels = c(0,1), labels = c("non_canopy","canopy"))

cm_p <- confusionMatrix(y_pred_p, y_true_p, positive = "canopy")
print(cm_p)
