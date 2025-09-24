library(ggplot2)
library(cowplot)

make_fire_effects_panel <- function(burned_model,
                                    unburned_model,
                                    scaled_model,
                                    plot_title) {
  
  # 1) Left: effect‐size plot
  p1 <- plot_fire_effects(
    burned_model   = burned_model,
    unburned_model = unburned_model,
    plot_title     = plot_title,
    label_map     = label_map
  ) +
    theme(legend.position = "bottom")
  
  # 2) Right: absolute magnitudes (stripped axes)
  p2 <- plot_effects_mag(
    scaled_model,
    plot_title = ""
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
  combined <- plot_grid(
    p1, p2,
    ncol       = 2,
    rel_widths = c(6, 1),
    align      = "h",
    axis       = "y"
  )
  
  return(combined)
}
