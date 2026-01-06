make_model_table_struc <- function(scaled_model, unscaled_model, label_map = label_map) {
  # make sure dplyr verbs are available
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Please install the dplyr package")
  }
  
  # 1) extract burned‐model fixed effects + SE + CIs
  coefs_burned <- fixef(scaled_model)
  se_burned    <- sqrt(diag(vcov(scaled_model)))
  
  tvals <- summary(scaled_model)$beta_table[, "t-value"]
  
  effect_df <- data.frame(
    term            = names(coefs_burned),
    standardized_estimate = round(coefs_burned, 5),
    std.error       = se_burned,
    lower           = coefs_burned - 1.96 * se_burned,
    upper           = coefs_burned + 1.96 * se_burned,
    tval           = tvals,
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
      significant = dplyr::case_when(
        abs(tval) > 3.3  ~ "***",
        abs(tval) > 2.6  ~ "**",
        abs(tval) > 1.96 ~ "*",
        TRUE             ~ ""
      )
    ) %>%
    dplyr::select(
      term,
      standardized_estimate,
      var_importance,
      significant
    )
  
 
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
  
  ames <- c()
  
  ame_terms <- as.data.frame(effect_df$term[effect_df$term != "tree_dens:build_dens"])
  
  names(ame_terms) <- "term"
   
  for(i in 1:nrow(ame_terms)){
    ame <- ame_plus_delta(unscaled_model, ame_terms$term[i], 1)
    ame_terms$AME[i] <- signif(ame, 3)
  }

  effect_df <- effect_df %>% left_join(ame_terms, by = "term") %>% 
    mutate(standardized_estimate_sig = paste0(standardized_estimate,significant)) %>% 
    select(term, var_importance, standardized_estimate_sig, AME)
  
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
