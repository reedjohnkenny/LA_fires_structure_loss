library(spaMM)


# --- Inputs ---
fit      <- Eaton_struc_spamm_1                 # spaMM HLfit, binomial/logit
pred_var <- "tree_dens"             # predictor to vary

predict_prob <- function(fit, pred_var){
stopifnot(inherits(fit, "HLfit"))
df <- model.frame(fit)
stopifnot(pred_var %in% names(df))
if (!is.numeric(df[[pred_var]])) stop("Requested predictor must be numeric.")

# Reference row (medians for numeric; first level for factors)
ref_row <- as.data.frame(lapply(df, function(x) {
  if (is.numeric(x)) median(x, na.rm = TRUE)
  else if (is.factor(x)) levels(x)[1L]
  else if (is.logical(x)) FALSE
  else if (is.character(x)) unique(x)[1L]
  else x[1]
}), stringsAsFactors = FALSE)

# Values to predict at
x_min <- min(df[[pred_var]], na.rm = TRUE)
x_max <- max(df[[pred_var]], na.rm = TRUE)
newdat <- rbind(ref_row, ref_row)
rownames(newdat) <- c("min","max")
newdat[[pred_var]] <- c(x_min, x_max)

# Assume: fit (spaMM HLfit, binomial/logit), newdat with rows "min","max"
# Request prediction variances; try to also get covariance between rows
pr <- predict(
  fit, newdata = newdat, type = "link",
  variances = list(predVar = TRUE, cov = TRUE),
  re.form = ~0
)

eta   <- as.numeric(pr)
Veta  <- attr(pr, "predVar")

# Standard errors for each row
se_eta <- if (is.matrix(Veta)) sqrt(diag(Veta)) else sqrt(as.numeric(Veta))

# 95% Wald CI on link scale, then transform to probability
z <- 1.96
eta_lo <- eta - z * se_eta
eta_hi <- eta + z * se_eta

p    <- plogis(eta)
p_lo <- plogis(eta_lo)
p_hi <- plogis(eta_hi)

pred_tab <- data.frame(
  setting   = rownames(newdat),
  prob      = p,
  prob_lo95 = p_lo,
  prob_hi95 = p_hi,
  logit     = eta,
  se_logit  = se_eta,
  row.names = NULL
)

# Contrast: (max - min) on probability scale
# gradient wrt (eta_min, eta_max) for p_max - p_min
grad <- c(-p[1] * (1 - p[1]),  p[2] * (1 - p[2]))

# Use full covariance if available; otherwise assume independence
if (is.matrix(Veta)) {
  var_delta <- as.numeric(t(grad) %*% Veta %*% grad)
} else {
  var_delta <- sum((grad * se_eta)^2)  # independence approximation
}
se_delta <- sqrt(var_delta)

delta_p  <- p[2] - p[1]
delta_ci <- c(delta_p - z * se_delta, delta_p + z * se_delta)

list(
  point_estimates = pred_tab,
  difference_max_minus_min = data.frame(
    diff_prob = delta_p,
    lo95 = delta_ci[1],
    hi95 = delta_ci[2],
    min = signif(x_min,4),
    max = signif(x_max,4)
  )
)
}
