library(spaMM)
library(dplyr)

# Total changes between point estimates

source("~/Desktop/Urban_tree_fire/structure_analysis/Scripts/functions/prob_change.R")

Eaton_struc_spamm_1 <- readRDS("~/Desktop/Urban_tree_fire/structure_analysis/models/Eaton_struc_m1_unscaled.rda")

pred_vars <- names(model.frame(Eaton_struc_spamm_1))

pred_vars <- pred_vars[!pred_vars %in% c("(Intercept)", "destroyed")]

changes <- as.data.frame(matrix(ncol = 9, nrow = 0))
colnames(changes) <- c("prob_mean", "prob_mean_lo95", "prob_mean_hi95", "prob_min", "prob_min_lo95", "prob_min_hi95", "p_diff", "p_diff_high", "p_diff_low")

for(i in 1:length(pred_vars)){
  change = predict_prob(fit = Eaton_struc_spamm_1, pred_var = pred_vars[i])
  rownames(change) <- pred_vars[i]
  changes = rbind(changes, change)
}



# Palisades

Palisades_struc_spamm_1 <- readRDS("~/Desktop/Urban_tree_fire/structure_analysis/models/Palisades_struc_m1.rda")


pred_vars <- names(model.frame(Palisades_struc_spamm_1))

pred_vars <- pred_vars[!pred_vars %in% c("(Intercept)", "destroyed")]

pal_changes <- as.data.frame(matrix(ncol = 9, nrow = 0))
colnames(pal_changes) <- c("prob_mean", "prob_mean_lo95", "prob_mean_hi95", "prob_min", "prob_min_lo95", "prob_min_hi95", "p_diff", "p_diff_high", "p_diff_low")

for(i in 1:length(pred_vars)){
  change = predict_prob(fit = Palisades_struc_spamm_1, pred_var = pred_vars[i])
  rownames(change) <- pred_vars[i]
  pal_changes = rbind(pal_changes, change)
}


# Continuos prob curves

source("Scripts/functions/plot_pred_change.R")

eaton_terms <- names(fixef(Eaton_struc_spamm_scaled))[2:11]

eaton_mean_plots <- list()

for(i in 1:length(eaton_terms)){

eaton_mean_plots[[i]] <- plot_pred_change(Eaton_struc_spamm_1, eaton_terms[i])

}

eaton_mean_plots[[9]]


plot_pred_change(Eaton_struc_spamm_1, eaton_terms[3])


# Palisades

pal_terms <- names(fixef(Palisades_struc_spamm_scaled))[2:11]

pal_mean_plots <- list()

for(i in 1:length(pal_terms)){
  
  pal_mean_plots[[i]] <- plot_pred_change(Palisades_struc_spamm_1, pal_terms[i])
  
}

pal_mean_plots[[3]]


plot_pred_change(Eaton_struc_spamm_1, eaton_terms[3])

