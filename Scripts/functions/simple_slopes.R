

simple_slopes_plot <- function(fit, x, modx,
                               modx_vals = c("Q10","Q50","Q90"),
                               xlim,
                               n_x = 100, controls = list()){
  d <- model.frame(fit)
  
  # choose moderator values by quantiles or numeric constants
  pick_vals <- function(z, specs){
    qmap <- c(Q10=0.10, Q25=0.25, Q50=0.50, Q75=0.75, Q90=0.90)
    sapply(specs, function(s){
      if (grepl("^Q", s)) unname(signif(quantile(z, qmap[[s]], na.rm=TRUE),2)) else as.numeric(s)
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
  newdat$predicted <- as.numeric(predict(fit, newdata = newdat, re.form = NA, type = "response"))
  
  # Get 95% confidence intervals for fixed effect predictions
  intervals <- get_intervals(fit, newdata = newdat, intervals = "predVar", re.form = NA)
  newdat <- cbind(newdat, intervals)
  
  ggplot(newdat, aes(x = .data[[x]], y = predicted, color = as.factor(.data[[modx]]), fill = as.factor(.data[[modx]]))) +
    geom_ribbon(aes(ymin = predVar_0.025, ymax = predVar_0.975), alpha = 0.15) +
    geom_line(linewidth = 0.5) +
    xlim(0,xlim) +
    labs(y = "Predicted Y",
         title = sprintf("Simple slopes of %s at fixed values of %s", x, modx)) +
    theme_classic()
}