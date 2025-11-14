suppressPackageStartupMessages({
  library(spaMM)
  library(dplyr)
  library(ggplot2)
  library(purrr)
})

plot_pred_change <- function(model, var){


# ---------------------------
# A) "Representative-covariate" effect with 95% CI
#    (keep all other covariates fixed at typical values)
# ---------------------------

# Helper: typical (mean for numeric, mode for factor)
typical_row <- model$data %>%
  summarise(across(everything(), \(x) {
    if (is.numeric(x)) mean(x, na.rm = TRUE)
    else if (is.factor(x)) names(sort(table(x), decreasing = TRUE))[1]
    else if (is.logical(x)) FALSE
    else if (is.character(x)) names(sort(table(x), decreasing = TRUE))[1]
    else NA
  }))

var_seq <- seq(min(model$data[[var]], na.rm = TRUE),
              max(model$data[[var]], na.rm = TRUE),
              length.out = 200)

grid_rep <- typical_row[rep(1, length(var_seq)), , drop = FALSE]
grid_rep[[var]] <- var_seq

# Predict on link scale to get standard errors; transform to probability
pred_rep <- predict(model, newdata = grid_rep, type = "link",
                    variances = list(predVar = TRUE),
                    re.form = ~0)

eta <- as.numeric(pred_rep)  # 200x1 vector of link-scale predictions

# Pointwise SEs:
se_fix  <- sqrt(attr(pred_rep, "fixefVar"))  # CI for the MEAN (fixed-effect uncertainty)
se_pred <- sqrt(attr(pred_rep, "predVar"))   # prediction SE (typically wider)

# Transform to probability (logit^{-1})

df <- tibble(
  !!var := var_seq,
  fit  = plogis(eta),
  
  lwr_CI = plogis(eta - 1.96*se_fix),
  upr_CI = plogis(eta + 1.96*se_fix),
  
  lwr_PI = plogis(eta - 1.96*se_pred),
  upr_PI = plogis(eta + 1.96*se_pred)
)

# Predicted probabilities at observed values

eta0 <- predict(model, type = "link", re.form = NA)
p0   <- plogis(eta0)

nd1 <- model$data
nd1[[var]] <- nd1[[var]] + 1  # x is your focal predictor
eta1 <- predict(model, newdata = nd1, type = "link", re.form = NA)
p1   <- plogis(eta1)

AME_exact <- mean(p1 - p0, na.rm = TRUE)   # average change in probability for +1 in x
AME_exact

xpos <- max(nd1[[var]])

ypos <- 0.1


ggplot(df, aes(.data[[var]], fit)) +
  # prediction band (usually wider)
  geom_ribbon(aes(ymin = lwr_PI, ymax = upr_PI), alpha = 0.12) +
  # confidence band for the mean
  geom_ribbon(aes(ymin = lwr_CI, ymax = upr_CI), alpha = 0.25) +
  geom_line(size = 0.9) +
  ylim(c(0,1)) +
  annotate("text",
           x = xpos, y = ypos,
           label = sprintf("Marginal Effect = %.6f", signif(AME_exact, 3)),
           hjust = 1, vjust = 1, size = 3) +
  labs(x = paste(var),
       y = "Predicted probability of structure loss") +
  theme_minimal() +
  theme(text = element_text(size = 8))


}
