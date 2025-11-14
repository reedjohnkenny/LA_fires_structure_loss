suppressPackageStartupMessages({
  library(spaMM)
  library(dplyr)
  library(ggplot2)
  library(purrr)
  library(rlang)
})

plot_pred_change <- function(model1, model2, var,
                             labels = c("Model 1", "Model 2"),
                             delta = 1,
                             n_points = 200) {
  
  stopifnot(length(labels) == 2)
  
  # --- helpers ---------------------------------------------------------------
  typical_row <- function(df) {
    summarise(df, across(everything(), \(x) {
      if (is.numeric(x)) mean(x, na.rm = TRUE)
      else if (is.factor(x)) names(sort(table(x), decreasing = TRUE))[1]
      else if (is.logical(x)) FALSE
      else if (is.character(x)) names(sort(table(x), decreasing = TRUE))[1]
      else NA
    }))
  }
  
  pred_curve <- function(model, var, var_seq, label) {
    df <- model$data
    if (!var %in% names(df))
      stop(sprintf("'%s' not found in model$data.", var))
    
    # grid at typical covariates
    typ <- typical_row(df)
    grid <- typ[rep(1, length(var_seq)), , drop = FALSE]
    grid[[var]] <- var_seq
    
    # population-level curve & SEs on link; then logit^{-1}
    pr <- predict(model, newdata = grid, type = "link",
                  variances = list(predVar = TRUE),
                  re.form = ~0)
    
    eta     <- as.numeric(pr)
    se_fix  <- sqrt(attr(pr, "fixefVar"))
    se_pred <- sqrt(attr(pr, "predVar"))
    
    tibble(
      .var   = var_seq,
      fit    = plogis(eta),
      lwr_CI = plogis(eta - 1.96 * se_fix),
      upr_CI = plogis(eta + 1.96 * se_fix),
      lwr_PI = plogis(eta - 1.96 * se_pred),
      upr_PI = plogis(eta + 1.96 * se_pred),
      model  = label
    )
  }
  
  ame_plus_delta <- function(model, var, delta) {
    df <- model$data
    if (!is.numeric(df[[var]]))
      stop("AME defined here for numeric predictors only.")
    # conditional predictions incl. REs to match user’s original
    eta0 <- predict(model, type = "link", re.form = NULL)
    p0   <- plogis(eta0)
    nd1  <- df
    nd1[[var]] <- nd1[[var]] + delta
    eta1 <- predict(model, newdata = nd1, type = "link", re.form = NULL)
    p1   <- plogis(eta1)
    mean(p1 - p0, na.rm = TRUE)
  }
  
  # --- range logic (use overlap to avoid extrapolation) ----------------------
  x1 <- model1$data[[var]]
  x2 <- model2$data[[var]]
  if (!is.numeric(x1) || !is.numeric(x2))
    stop("Both models must have a numeric 'var' for this plot.")
  
  rng1 <- range(x1, na.rm = TRUE)
  rng2 <- range(x2, na.rm = TRUE)
  x_min <- max(rng1[1], rng2[1])
  x_max <- min(rng1[2], rng2[2])
  
  if (!is.finite(x_min) || !is.finite(x_max) || x_min >= x_max) {
    stop("No overlapping support for 'var' between the two models; cannot plot without extrapolation.")
  }
  
  var_seq <- seq(x_min, x_max, length.out = n_points)
  
  # --- curves & AMEs ---------------------------------------------------------
  df1  <- pred_curve(model1, var, var_seq, labels[1])
  df2  <- pred_curve(model2, var, var_seq, labels[2])
  dff  <- bind_rows(df1, df2)
  
  ame1 <- ame_plus_delta(model1, var, delta)
  ame2 <- ame_plus_delta(model2, var, delta)
  
  # annotate near the right edge at fixed y offsets
  x_annot <- x_max
  ann <- tibble(
    x = c(x_annot, x_annot),
    y = c(0.10, 0.20),
    lab = c(
      sprintf("%s: AME (Δ=%g) = %.4f", labels[1], delta, ame1),
      sprintf("%s: AME (Δ=%g) = %.4f", labels[2], delta, ame2)
    ),
    model = labels
  )
  
  # --- plot ------------------------------------------------------------------
  ggplot(dff, aes(x = .var, y = fit, color = model, fill = model)) +
    # prediction band (wider)
    geom_ribbon(aes(ymin = lwr_PI, ymax = upr_PI), alpha = 0.10, linewidth = 0) +
    # CI for mean
    geom_ribbon(aes(ymin = lwr_CI, ymax = upr_CI), alpha = 0.25, linewidth = 0) +
    geom_line(size = 0.9) +
    coord_cartesian(ylim = c(0, 1)) +
    geom_text(data = ann, aes(x = x, y = y, label = lab, color = model),
              inherit.aes = FALSE, hjust = 1, vjust = 1, size = 3) +
    labs(x = var,
         y = "Predicted probability of structure loss",
         color = "Model", fill = "Model") +
    theme_minimal() +
    theme(text = element_text(size = 9))
}
