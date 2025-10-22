compare_slope_welch <- function(model1, model2,
                                coef1,
                                coef2 = coef1,
                                df1   = NULL,
                                df2   = NULL) {
  # Helper to pull est & SE from a summary matrix
  get_est_se <- function(model, coef_name) {
    # Extract fixed effects and their standard errors
    coefs <- fixef(model)
    se <- sqrt(diag(vcov(model)))  # standard errors
    
    if (!coef_name %in% names(coefs)) {
      stop("Coefficient '", coef_name, "' not found in model.")
    }
    
    # Create data frame of coefficients
    effect_df <- data.frame(
      term = names(coefs),
      estimate = coefs,
      std.error = se,
      lower = coefs - 1.96 * se,
      upper = coefs + 1.96 * se
    )
    
    est <- effect_df[coef_name, "estimate"]
    se  <- effect_df[coef_name, "std.error"]
    list(est = est, se = se)
    
  }
  
  # Extract estimates and SEs
  a1 <- get_est_se(model1, coef1)
  a2 <- get_est_se(model2, coef2)
  
  # Infer dfs if not supplied
  infer_df <- function(mod, supplied_df) {
    if (!is.null(supplied_df)) return(supplied_df)
    if (!is.null(mod$data)) {
      n  <- nrow(mod$data)
      p  <- length(stats::coef(mod))
      return(n - p)
    }
    stop("Cannot infer df for model; please supply df1/df2 manually.")
  }
  df1 <- infer_df(model1, df1)
  df2 <- infer_df(model2, df2)
  
  # Welch t
  diff_est <- a1$est - a2$est
  diff_se  <- sqrt(a1$se^2 + a2$se^2)
  t_val    <- diff_est / diff_se
  
  # Satterthwaite df
  df_welch <- (a1$se^2 + a2$se^2)^2 /
    ((a1$se^4 / df1) + (a2$se^4 / df2))
  
  p_val <- 2 * stats::pt(-abs(t_val), df = df_welch)
  
  # Return a list akin to stats::t.test
  structure(
    list(
      statistic   = t_val,
      parameter   = df_welch,
      p.value     = p_val,
      estimate    = c(model1 = a1$est, model2 = a2$est),
      stderr      = c(se1 = a1$se, se2 = a2$se),
      null.value  = 0,
      conf.int    = diff_est + c(-1,1) * stats::qt(0.975, df_welch) * diff_se,
      estimate.diff = diff_est,
      stderr.diff   = diff_se,
      method        = "Welch t-test for difference in spaMM slopes",
      data.name     = paste0(deparse(substitute(model1)), " vs ", deparse(substitute(model2)))
    ),
    class = "htest"
  )
}
