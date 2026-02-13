library(tidyverse)


setwd("~/Desktop/Urban_tree_fire/structure_analysis/validation_data/")

# Eaton


checks <- read.csv("rmse_checkpoints_eaton.csv")

ref <- checks %>% filter(image == "ref")

shifted <- checks %>% filter(image == "shifted")

rmse <- sqrt(mean((shifted$x - ref$x)^2 + (shifted$y - ref$y)^2))


# Palisades 

checks <- read.csv("rmse_checkpoints_palisades.csv")

ref <- checks %>% filter(image == "ref")

shifted <- checks %>% filter(image == "shifted")

rmse <- sqrt(mean((shifted$x - ref$x)^2 + (shifted$y - ref$y)^2))
