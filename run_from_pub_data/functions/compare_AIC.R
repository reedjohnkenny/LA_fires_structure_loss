
compare_AIC <- function(..., models = NULL) {
  # 1. Figure out if user passed a list or individual models
  if (!is.null(models)) {
    if (!is.list(models)) 
      stop("`models` must be a list of fitted model objects")
    model_list <- models
    # if the list is unnamed, give it default names
    if (is.null(names(model_list))) {
      names(model_list) <- paste0("Model", seq_along(model_list))
    }
  } else {
    # grab varargs as list
    model_list <- list(...)
    call_args  <- match.call(expand.dots = FALSE)$...
    names(model_list) <- sapply(call_args, deparse)
  }
  
  # 2. Prepare storage
  n <- length(model_list)
  mAICs <- cAICs <- dAICs <- GoFdfs <- numeric(n)
  
  # 3. Loop through each model, extract its AIC table
  for (i in seq_along(model_list)) {
    aic_df <- as.data.frame(AIC(model_list[[i]]))
    # rename rows
    row.names(aic_df) <- c("mAIC", "cAIC", "dAIC", "GoFdf")
    # rename the column to values
    names(aic_df) <- c("values")
    
    # now extract by short.name
    mAICs[i] <- aic_df$value[rownames(aic_df) == "mAIC"]
    cAICs[i] <- aic_df$value[rownames(aic_df) == "cAIC"]
    dAICs[i] <- aic_df$value[rownames(aic_df) == "dAIC"]
    GoFdfs[i] <- aic_df$value[rownames(aic_df) == "GoFdf"]
  }
  
  # 4. Assemble into a data.frame and sort
  result <- data.frame(
    Model = names(model_list),
    mAIC  = mAICs,
    cAIC  = cAICs,
    dAIC  = dAICs,
    GoFdf = GoFdfs,
    stringsAsFactors = FALSE
  )
  result[order(result$mAIC), , drop = FALSE]
}


