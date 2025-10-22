
plot_alt_models <- function(model){


coefs_E <- fixef(model)
se_E    <- sqrt(diag(vcov(model)))


# Build data frames for each model
effect_df_E<- data.frame(
  term      = names(coefs_E),
  estimate  = coefs_E,
  std.error = se_E,
  lower     = coefs_E - 1.96 * se_E,
  upper     = coefs_E + 1.96 * se_E,
  fire      = "Eaton",
  stringsAsFactors = FALSE
)

coefs_scale <- fixef(model)
ses   <- sqrt(diag(vcov(model)))

# build summary data.frame
df <- data.frame(
  term     = names(coefs_scale),
  estimate = coefs_scale,
  std.error= ses,
  lower    = coefs_scale - 1.96 * ses,
  upper    = coefs_scale + 1.96 * ses,
  stringsAsFactors = FALSE
)

ordered_terms <- df %>%
  filter(term != "(Intercept)") %>%
  arrange(abs(estimate)) %>%
  pull(term)

ordered_labels <- ordered_terms

effect_df_E <- effect_df_E %>%
  mutate(
    significant = (lower > 0) | (upper < 0),
    label       = ifelse(significant, "*", "")
  ) %>%
  filter(!effect_df_E$term %in% "(Intercept)") %>% 
  # Ensure consistent ordering of terms (alphabetical)
  mutate(term = factor(term, levels = ordered_terms, 
                       labels = ordered_labels))

# Create the point + error‐bar plot
p <- ggplot(effect_df_E, aes(x = term, y = estimate, color = fire)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  geom_text(aes(label = label),
            nudge_y = -0.0001,
            nudge_x = 0.1,
            size = 5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  labs(
    title = "Standardized effect sizes predicting stucture loss Eaton",
    x     = "Predictor",
    y     = "Estimate (± 95% CI)"
  ) +
  scale_color_manual(values = c("#440154FF" , "#21908CFF")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position = "none")



p2 <- plot_effects_mag(
  model,
  plot_title = "", 
  bar_fill = "#440154FF"
) +
  theme(
    axis.title.y    = element_blank(),
    axis.text.y     = element_blank(),
    axis.ticks.y    = element_blank(),
    axis.title.x    = element_blank(),
    axis.text.x     = element_blank(),
    axis.ticks.x    = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# 3) Combine with aligned y‐axes, 6:1 width ratio
plot_grid(
  p, p2,
  ncol       = 2,
  rel_widths = c(6, 1),
  align      = "h",
  axis       = "y"
)
}
