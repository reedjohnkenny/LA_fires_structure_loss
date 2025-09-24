library(ggplot2)
library(dplyr)
library(spaMM)

plot_model_effects <- function(model, plot_title = "Effect Sizes with 95% CI") {
  # Extract fixed effects and their standard errors
  coefs <- fixef(model)
  se <- sqrt(diag(vcov(model)))  # standard errors
  
  # Create data frame of coefficients
  effect_df <- data.frame(
    term = names(coefs),
    estimate = coefs,
    std.error = se,
    lower = coefs - 1.96 * se,
    upper = coefs + 1.96 * se
  )
  
  # Remove the intercept for plotting
  effect_df <- effect_df[!effect_df$term %in% "(Intercept)", ]
  
  # Add significance label
  effect_df <- effect_df %>%
    mutate(
      significant = lower > 0 | upper < 0,
      label = ifelse(significant, "*", "")
    )
  
  # Ensure terms are ordered alphabetically
  effect_df <- effect_df %>%
    mutate(term = factor(term, levels = sort(unique(term))))
  
  # Plot
  ggplot(effect_df, aes(x = term, y = estimate)) +
    geom_point() +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) + 
    geom_text(aes(label = label), nudge_y = -0.0001, nudge_x = 0.1, size = 5) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    coord_flip() +
    labs(
      title = plot_title,
      x = "Predictor",
      y = "Estimate (± 95% CI)"
    ) +
    theme_minimal()
}



