

source("~/Desktop/Urban_tree_fire/landscape_analysis/Scripts/functions/human_labs.R")



library(dplyr)
library(ggplot2)


plot_fire_effects <- function(burned_model,
                              unburned_model,
                              plot_title,
                              label_map = label_map) {
  # 1. Extract fixed‐effects coefficients and standard errors
  coefs_burned   <- fixef(burned_model)
  se_burned      <- sqrt(diag(vcov(burned_model)))
  coefs_unburned <- fixef(unburned_model)
  se_unburned    <- sqrt(diag(vcov(unburned_model)))
  
  # 2. Build data frames for each model
  effect_df_b <- data.frame(
    term      = names(coefs_burned),
    estimate  = coefs_burned,
    std.error = se_burned,
    lower     = coefs_burned - 1.96 * se_burned,
    upper     = coefs_burned + 1.96 * se_burned,
    burn_area = "yes",
    stringsAsFactors = FALSE
  )
  effect_df_u <- data.frame(
    term      = names(coefs_unburned),
    estimate  = coefs_unburned,
    std.error = se_unburned,
    lower     = coefs_unburned - 1.96 * se_unburned,
    upper     = coefs_unburned + 1.96 * se_unburned,
    burn_area = "no",
    stringsAsFactors = FALSE
  )
  
  # 3. Combine, drop intercept, add significance flag
  df <- bind_rows(effect_df_u, effect_df_b) %>%
    filter(term != "(Intercept)") %>%
    mutate(
      significant = (lower > 0) | (upper < 0),
      label       = ifelse(significant, "*", "")
    )
  
  # 4. Determine ordering by the burned‐model estimates (descending)
  ordered_terms <- effect_df_b %>%
    filter(term != "(Intercept)") %>%
    arrange(abs(estimate)) %>%
    pull(term)
  
  # 5. Apply factor levels (and optional relabeling)
  if (!is.null(label_map)) {
    # sanity check
    missing <- setdiff(ordered_terms, names(label_map))
    if (length(missing)) {
      stop("Missing labels for: ", paste(missing, collapse = ", "))
    }
    ordered_labels <- label_map[ordered_terms]
    df <- df %>%
      mutate(term = factor(term,
                           levels = ordered_terms,
                           labels = ordered_labels))
  } else {
    df <- df %>%
      mutate(term = factor(term,
                           levels = ordered_terms))
  }
  
  df <- df %>%
    mutate(burn_area = factor(burn_area, levels = c("no", "yes")))
  
  # 6. Plot
  ggplot(df, aes(x = term, y = estimate, color = burn_area)) +
    geom_point(alpha = 0.7) +
    geom_errorbar(aes(ymin = lower, ymax = upper),
                  width = 0.2, alpha = 0.7) +
    geom_text(aes(label = label),
              nudge_y = -0.0001, nudge_x = 0.1, size = 5) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    coord_flip() +
    labs(title = plot_title,
         x     = "Predictor",
         y     = "Estimate (± 95% CI)",
         color = "Burn area") +
    scale_color_manual(values = c("no" = "black", "yes" = "red")) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
}



