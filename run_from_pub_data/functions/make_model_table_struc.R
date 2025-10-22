make_model_table_struc <- function(burned_model, scaled_model, label_map = label_map) {
  # make sure dplyr verbs are available
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Please install the dplyr package")
  }
  
  # 1) extract burned‐model fixed effects + SE + CIs
  coefs_burned <- fixef(scaled_model)
  se_burned    <- sqrt(diag(vcov(scaled_model)))
  
  effect_df <- data.frame(
    term            = names(coefs_burned),
    estimate = round(coefs_burned, 5),
    std.error       = se_burned,
    lower           = coefs_burned - 1.96 * se_burned,
    upper           = coefs_burned + 1.96 * se_burned,
    stringsAsFactors = FALSE
  )
  
 # 3) drop intercept
  effect_df <- effect_df[effect_df$term != "(Intercept)", ]
  
  # 5) variable importance from the scaled model
  coefs_scaled <- abs(fixef(scaled_model))
  # drop intercept
  coefs_scaled <- coefs_scaled[names(coefs_scaled) != "(Intercept)"]
  
  # align order to effect_df$term, then rank
  scaled_for_terms <- coefs_scaled[effect_df$term]
  effect_df$var_importance <- rank(
    -scaled_for_terms,
    ties.method = "first"
  )
  
  # 6) add significance flag, then drop unneeded columns
  effect_df <- effect_df %>%
    dplyr::mutate(
      sig_brn_mod = ifelse(lower > 0 | upper < 0, "yes", "no")
    ) %>%
    dplyr::select(
      term,
      estimate,
      var_importance,
      sig_brn_mod
    )
  
  # 5) variables from the unscaled model
  coefs_raw <- fixef(burned_model)
  # drop intercept
  coefs_raw <- coefs_raw[names(coefs_raw) != "(Intercept)"]
  
  # align order to effect_df$term, then rank
  raw_for_terms <- coefs_raw[effect_df$term]
  
  raw_df <- data.frame(
    raw_term            = names(coefs_raw),
    raw_estimate = round(coefs_raw, 5),
    odds_ratio = round(exp(coefs_raw),5))
  
  
  effect_df <- effect_df %>% left_join(raw_df, by = c("term" = "raw_term"))
  
  
  
  # 4. Determine ordering by the burned‐model estimates (descending)
  ordered_terms <- effect_df %>%
    filter(term != "(Intercept)") %>%
    arrange(var_importance) %>%
    pull(term)
  
  # 5. Apply factor levels (and optional relabeling)
  if (!is.null(label_map)) {
    # sanity check
    missing <- setdiff(ordered_terms, names(label_map))
    if (length(missing)) {
      stop("Missing labels for: ", paste(missing, collapse = ", "))
    }
    ordered_labels <- label_map[ordered_terms]
    effect_df <- effect_df %>%
      mutate(term = factor(term,
                           levels = ordered_terms,
                           labels = ordered_labels))
  } else {
    effect_df <- effect_df %>%
      mutate(term = factor(term,
                           levels = ordered_terms))
  }
  
  effect_df <- effect_df %>% arrange(var_importance)
  

  return(effect_df)
}
