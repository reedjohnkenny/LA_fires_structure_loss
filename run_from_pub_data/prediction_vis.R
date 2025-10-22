library(spaMM)
library(dplyr)

source("~/Desktop/Urban_tree_fire/landscape_analysis/Scripts/functions/prob_change.R")

Eaton_struc_spamm_1 <- readRDS("~/Desktop/Urban_tree_fire/landscape_analysis/simple_models/Eaton_struc_m1_unscaled.rda")

pred_vars <- names(model.frame(Eaton_struc_spamm_1))

pred_vars <- pred_vars[!pred_vars %in% c("(Intercept)", "destroyed")]

changes <- as.data.frame(matrix(ncol = 3, nrow = 0))
colnames(changes) <- c("diff_prob", "lo95", "hi95")

for(i in 1:length(pred_vars)){
  change = predict_prob(fit = Eaton_struc_spamm_1, pred_var = pred_vars[i])[[2]]
  rownames(change) <- pred_vars[i]
  changes = rbind(changes, change)
}



predict_prob(fit = Palisades_struc_spamm_1, pred_var = "build_dens:tree_dens")


# Palisades

Palisades_struc_spamm_1 <- readRDS("~/Desktop/Urban_tree_fire/landscape_analysis/simple_models/Palisades_struc_m1.rda")


pred_vars <- names(model.frame(Palisades_struc_spamm_1))

pred_vars <- pred_vars[!pred_vars %in% c("(Intercept)", "destroyed")]

pal_changes <- as.data.frame(matrix(ncol = 3, nrow = 0))
colnames(pal_changes) <- c("diff_prob", "lo95", "hi95")

for(i in 1:length(pred_vars)){
  change = predict_prob(fit = Palisades_struc_spamm_1, pred_var = pred_vars[i])[[2]]
  rownames(change) <- pred_vars[i]
  pal_changes = rbind(pal_changes, change)
}

