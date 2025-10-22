library(ggplot2)
library(dplyr)
library(spaMM)

plot_fixed_effect_relationship <- function(model, predictor, grid_length = 20, pred_min = 0, pred_max = 20, line_col = "red") {
  # Extract model data
  model_data <- model$data
  
  # Confirm predictor exists in data
  if (!predictor %in% colnames(model_data)) {
    stop("Predictor not found in model data.")
  }
  
  # Get response variable name
  response_var <- as.character(model$predictor[2])
  
  # Generate a sequence of values for the chosen predictor
  predictor_seq <- seq(
    pred_min, pred_max,
    length.out = grid_length
  )
  
  # Create a data frame where all predictors are set to their median (or first level if factor)
  median_vals <- model_data %>%
    summarise(across(everything(), ~ if (is.numeric(.)) median(., na.rm = TRUE) else .[1]))
  
  newdat <- median_vals[rep(1, grid_length), ]
  newdat[[predictor]] <- predictor_seq
  
  # Predict fixed effects only (remove spatial/random effects)
  newdat$predicted <- as.numeric(predict(model, newdata = newdat, re.form = NA))
  
  # Get 95% confidence intervals for fixed effect predictions
  intervals <- get_intervals(model, newdata = newdat, intervals = "fixefVar", re.form = NA)
  newdat <- cbind(newdat, intervals)
  
  # Plot observed vs predicted
  ggplot(model_data, aes_string(x = predictor, y = response_var)) +
    geom_point(color = "blue", alpha = 0.1) +
    geom_path(data = newdat, aes_string(x = predictor, y = "predicted"), color = line_col) +
    geom_ribbon(data = newdat, aes(ymin = fixefVar_0.025, ymax = fixefVar_0.975), alpha = 0.2) +
    theme_bw() +
    xlab(predictor) +
    ylab(response_var) +
    xlim(c(pred_min, pred_max))
}
