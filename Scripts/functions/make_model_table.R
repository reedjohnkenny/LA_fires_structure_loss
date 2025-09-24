make_model_table <- function(burned_model, unburned_model, scaled_model, label_map = NULL) {
  # make sure dplyr verbs are available
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Please install the dplyr package")
  }
  
  # 1) extract burned‐model fixed effects + SE + CIs
  coefs_burned <- fixef(burned_model)
  se_burned    <- sqrt(diag(vcov(burned_model)))
  
  effect_df <- data.frame(
    term            = names(coefs_burned),
    estimate_burned = round(coefs_burned, 6),
    std.error       = se_burned,
    lower           = coefs_burned - 1.96 * se_burned,
    upper           = coefs_burned + 1.96 * se_burned,
    stringsAsFactors = FALSE
  )
  
  # 2) add unburned estimates
  coefs_unburned <- fixef(unburned_model)
  effect_df$estimate_unburned <- round(coefs_unburned[effect_df$term], 6)
  
  # 3) drop intercept
  effect_df <- effect_df[effect_df$term != "(Intercept)", ]
  
  # 4) compute Welch p‐values for each term
  effect_df$Welch_p_val_brn_unbrn <- vapply(
    effect_df$term,
    function(t) compare_slope_welch(
      model1 = burned_model,
      model2 = unburned_model,
      coef1  = t
    )$p.value,
    numeric(1)
  )
  effect_df$Welch_p_val_brn_unbrn <- round(effect_df$Welch_p_val_brn_unbrn, 3)
  
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
      estimate_burned,
      estimate_unburned,
      Welch_p_val_brn_unbrn,
      var_importance,
      sig_brn_mod
    )
  
  # 4. Determine ordering by 
  ordered_terms <- effect_df %>%
    arrange(var_importance) %>% 
    pull(term)
  
  
  #Apply factor levels (and optional relabeling)
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
