library(dplyr)
library(forcats)
library(ggplot2)



plot_effects_mag <- function(model, plot_title, bar_fill, label_map = label_map) {
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
  

    ordered_labels <- label_map[ordered_terms]
    df <- df %>%
      mutate(term = factor(term,
                           levels = ordered_terms,
                           labels = ordered_labels))

  
    # drop intercept, compute magnitude & significance
    df <- df %>% filter(term != "(Intercept)") %>%
    mutate(
      mag_est     = round(abs(estimate),2),
      significant = (lower > 0) | (upper < 0),
      label       = ifelse(significant, paste0(mag_est,"*"), paste0(mag_est,""))
    )
  
  # plot
  ggplot(df, aes(x = mag_est, y = term)) +
    geom_col(fill = bar_fill) +
    geom_text(aes(label = label),  # round values if desired
              hjust = -0.1,                    # adjust horizontal position
              size = 2) +      
    labs(
      title = plot_title,
      x     = "Absolute Magnitude of Standardized Estimate",
      y     = NULL
    ) +
    theme_minimal(base_size = 6) +
    xlim(0, max(df$mag_est) * 1.1) 
}
