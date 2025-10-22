library(ggplot2)
library(dplyr)

simple_slopes_plot <- function(fit, x, modx,
                               modx_vals = c("Q10","Q50","Q90"),
                               xlim,
                               n_x = 100, controls = list()){
  d <- model.frame(fit)
  
  # choose moderator values by quantiles or numeric constants
  pick_vals <- function(z, specs){
    qmap <- c(Q10=0.10, Q25=0.25, Q50=0.50, Q75=0.75, Q90=0.90)
    sapply(specs, function(s){
      if (grepl("^Q", s)) unname(quantile(z, qmap[[s]], na.rm=TRUE)) else as.numeric(s)
    })
  }
  m_vals <- pick_vals(d[[modx]], modx_vals)
  
  x_seq  <- seq(min(d[[x]], na.rm=TRUE), max(d[[x]], na.rm=TRUE), length.out = n_x)
  
  new_terms <- do.call(expand.grid, setNames(list(x_seq, m_vals), c(x, modx)))
  
  median_vals <- fit$data %>%
    summarise(across(everything(), ~ if (is.numeric(.)) median(., na.rm = TRUE) else .[1]))
  
  newdat <- median_vals[rep(1, length(modx_vals)*n_x), ]
  
  newdat[[x]] <- new_terms[,1]
  newdat[[modx]] <- new_terms[,2]
  
  # Predict fixed effects only (remove spatial/random effects)
  newdat$predicted <- as.numeric(predict(fit, newdata = newdat, re.form = NA))
  
  # Get 95% confidence intervals for fixed effect predictions
  intervals <- get_intervals(fit, newdata = newdat, intervals = "fixefVar", re.form = NA)
  newdat <- cbind(newdat, intervals)
  
  ggplot(newdat, aes(x = .data[[x]], y = predicted, color = as.factor(.data[[modx]]))) +
    geom_ribbon(aes(ymin = fixefVar_0.025, ymax = fixefVar_0.975, fill = as.factor(.data[[modx]])), alpha = 0.15) +
    geom_line(linewidth = 1) +
    xlim(0,xlim) +
    labs(y = "Predicted Y",
         title = sprintf("Simple slopes of %s at fixed values of %s", x, modx),
         subtitle = "Lines = fitted values, bands = 95% CI") +
    theme_classic()
}

p1 <- simple_slopes_plot(eaton_struc_m, x="distance_to_nearest_tree", modx="angular_diff_tree", modx_vals=c("Q10","Q50","Q90"), xlim = 10)
p1

p2 <- simple_slopes_plot(eaton_struc_m, x = "angular_diff_tree", modx="distance_to_nearest_tree", modx_vals=c(0.1,5,10), xlim = 200)
p2
