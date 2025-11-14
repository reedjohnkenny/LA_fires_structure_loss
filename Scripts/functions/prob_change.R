library(spaMM)


# --- Inputs ---
fit      <- Eaton_struc_spamm_1                 # spaMM HLfit, binomial/logit
pred_var <- "tree_dens"             # predictor to vary

predict_prob <- function(fit, pred_var){
  stopifnot(inherits(fit, "HLfit"))
  df <- model.frame(fit)
  stopifnot(pred_var %in% names(df))
  if (!is.numeric(df[[pred_var]])) stop("Requested predictor must be numeric.")
  
  # create new data at mean of pred_var and min pred_var
  
  new_dat <- fit$data
  
  new_dat_mean <- fit$data
  
  new_dat_mean[[pred_var]] <- mean(new_dat[[pred_var]])
  
  new_dat_min <- fit$data
  
  new_dat_min[[pred_var]] <- min(new_dat[[pred_var]])
  

  
  # Assume: fit (spaMM HLfit, binomial/logit), newdat with rows "min","max"
  # Request prediction variances; try to also get covariance between rows
  
  pr_mean <- predict(
    fit, newdata = new_dat_mean, type = "link",
    variances = list(predVar = TRUE, cov = TRUE),
    re.form = NULL
  )
  
  eta_mean   <- as.numeric(pr_mean)
  Veta_mean  <- attr(pr_mean, "predVar")
  
  pr_min <- predict(
    fit, newdata = new_dat_min, type = "link",
    variances = list(predVar = TRUE, cov = TRUE),
    re.form = NULL
  )
  
  eta_min   <- as.numeric(pr_min)
  Veta_min  <- attr(pr_min, "predVar")
  
  
  # Standard errors for each row
  se_eta_mean <- if (is.matrix(Veta_mean)) sqrt(diag(Veta_mean)) else sqrt(as.numeric(Veta_mean))
  
  se_eta_min <- if (is.matrix(Veta_min)) sqrt(diag(Veta_min)) else sqrt(as.numeric(Veta_min))
  
  # 95% Wald CI on link scale, then transform to probability
  z <- 1.96
  eta_lo_mean <- eta_mean - z * se_eta_mean
  eta_hi_mean <- eta_mean + z * se_eta_mean
  
  eta_lo_min <- eta_min - z * se_eta_min
  eta_hi_min <- eta_min + z * se_eta_min
  
  
  p_mean    <- plogis(eta_mean)
  p_lo_mean <- plogis(eta_lo_mean)
  p_hi_mean <- plogis(eta_hi_mean)
  
  p_min    <- plogis(eta_min)
  p_lo_min <- plogis(eta_lo_min)
  p_hi_min <- plogis(eta_hi_min)
  
  pred_tab <- data.frame(
    prob_mean      = p_mean,
    prob_mean_lo95 = p_lo_mean,
    prob_mean_hi95 = p_hi_mean,
    prob_min      = p_min,
    prob_min_lo95 = p_lo_min,
    prob_min_hi95 = p_hi_min,
    p_diff = p_mean - p_min,
    p_diff_high = p_hi_mean - p_lo_min,
    p_diff_low = p_lo_mean - p_hi_min,
    row.names = NULL
  )
  
  mean_pred_tab <- as.data.frame(as.list(colMeans(pred_tab)))
  
  return(mean_pred_tab)
  
  # Contrast: (max - min) on probability scale
  # gradient wrt (eta_min, eta_max) for p_max - p_min
  # grad <- c(-p[1] * (1 - p[1]),  p[2] * (1 - p[2]))
  # 
  # # Use full covariance if available; otherwise assume independence
  # if (is.matrix(Veta)) {
  #   var_delta <- as.numeric(t(grad) %*% Veta %*% grad)
  # } else {
  #   var_delta <- sum((grad * se_eta)^2)  # independence approximation
  # }
  # se_delta <- sqrt(var_delta)
  # 
  # delta_p  <- p[2] - p[1]
  # delta_ci <- c(delta_p - z * se_delta, delta_p + z * se_delta)
  # 
  # list(
  #   point_estimates = pred_tab,
  #   difference_max_minus_min = data.frame(
  #     diff_prob = delta_p,
  #     lo95 = delta_ci[1],
  #     hi95 = delta_ci[2],
  #     min = signif(x_min,4),
  #     max = signif(x_max,4)
  #   )
  #)
}

