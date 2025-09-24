library(dplyr)
library(forcats)
library(ggplot2)

plot_effects_mag <- function(model, plot_title = NULL, bar_fill = 'red') {
  # extract coefficients & SEs
  coefs <- fixef(model)
  ses   <- sqrt(diag(vcov(model)))
  
  # build summary data.frame
  df <- data.frame(
    term     = names(coefs),
    estimate = coefs,
    std.error= ses,
    lower    = coefs - 1.96 * ses,
    upper    = coefs + 1.96 * ses,
    stringsAsFactors = FALSE
  )
  
  ordered_terms <- df %>%
    filter(term != "(Intercept)") %>%
    arrange(abs(estimate)) %>%
    pull(term)
  
    # drop intercept, compute magnitude & significance
    df <- df %>% filter(term != "(Intercept)") %>%
    mutate(
      mag_est     = abs(estimate),
      significant = (lower > 0) | (upper < 0),
      label       = ifelse(significant, "*", "")
    ) %>% 
    # Ensure terms are ordered alphabetically
    mutate(term = factor(term, levels = ordered_terms))
  
  # plot
  ggplot(df, aes(x = mag_est, y = term)) +
    geom_col(fill = bar_fill) +
    labs(
      title = plot_title,
      x     = "Magnitude of Estimate",
      y     = NULL
    ) +
    theme_minimal(base_size = 14)
}
