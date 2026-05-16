# Extract Risk Score Formula from StepCox[forward] + GBM Model - Version 2
# ============================================================

cat("Loading results object...\n")
load("res_no_GSE14520.Rdata")
load("QWEN0208_no_GSE14520.Rdata")

# Extract the StepCox[forward] + GBM model
model_name <- "StepCox[forward] + GBM"

cat("\n=== Extracting Risk Score Information ===\n")
cat("Model:", model_name, "\n\n")

# Check riskscore component
if(!is.null(res$riskscore)) {
  cat("=== Risk Score Component Structure ===\n")
  print(names(res$riskscore))

  if(model_name %in% names(res$riskscore)) {
    cat("\n=== Risk Scores for", model_name, "===\n")
    rs_data <- res$riskscore[[model_name]]
    print(str(rs_data, max.level = 2))

    # Check if there are risk scores for each dataset
    cat("\n=== Available Datasets ===\n")
    print(names(rs_data))
  }
}

# Check Sig.genes - these are the genes selected by StepCox
if(!is.null(res$Sig.genes)) {
  cat("\n=== Significant Genes Selected ===\n")

  if(model_name %in% names(res$Sig.genes)) {
    selected_genes <- res$Sig.genes[[model_name]]
    cat("Number of genes:", length(selected_genes), "\n")
    cat("Genes:\n")
    print(selected_genes)

    # Save gene list
    write.table(selected_genes,
                file = "StepCox_GBM_selected_genes.txt",
                row.names = FALSE, col.names = FALSE, quote = FALSE)
    cat("\nGene list saved to: StepCox_GBM_selected_genes.txt\n")
  }
}

# Extract the actual model object
if(!is.null(res$ml.res) && model_name %in% names(res$ml.res)) {
  model_obj <- res$ml.res[[model_name]]

  cat("\n=== Model Object Details ===\n")
  cat("Model type:", class(model_obj$fit), "\n")
  cat("Number of trees:", model_obj$fit$n.trees, "\n")
  cat("Best iteration:", model_obj$best, "\n")

  # Get variable importance
  if(!is.null(model_obj$fit$var.names)) {
    cat("\n=== Variables in Model ===\n")
    cat("Number of variables:", length(model_obj$fit$var.names), "\n")
    cat("Variables:\n")
    print(model_obj$fit$var.names)

    # Get relative influence (variable importance)
    library(gbm)
    var_imp <- summary(model_obj$fit, plotit = FALSE)
    cat("\n=== Variable Importance (Top 20) ===\n")
    print(head(var_imp, 20))

    # Save variable importance
    write.csv(var_imp,
              file = "StepCox_GBM_variable_importance.csv",
              row.names = FALSE)
    cat("\nVariable importance saved to: StepCox_GBM_variable_importance.csv\n")
  }

  # Create a comprehensive output file
  sink("StepCox_GBM_model_details.txt")
  cat("============================================================\n")
  cat("StepCox[forward] + GBM Model Details\n")
  cat("============================================================\n\n")
  cat("Model Type: Gradient Boosting Machine (GBM) with Cox PH loss\n")
  cat("Feature Selection: Stepwise Cox Regression (Forward)\n\n")

  cat("=== Model Parameters ===\n")
  cat("Number of trees:", model_obj$fit$n.trees, "\n")
  cat("Best iteration:", model_obj$best, "\n")
  cat("Interaction depth:", model_obj$fit$interaction.depth, "\n")
  cat("Shrinkage:", model_obj$fit$shrinkage, "\n")
  cat("Bag fraction:", model_obj$fit$bag.fraction, "\n")
  cat("Min observations in node:", model_obj$fit$n.minobsinnode, "\n\n")

  cat("=== Selected Features ===\n")
  cat("Number of features:", length(model_obj$fit$var.names), "\n")
  cat("Features:\n")
  for(i in 1:length(model_obj$fit$var.names)) {
    cat(sprintf("  %d. %s\n", i, model_obj$fit$var.names[i]))
  }

  cat("\n=== Variable Importance (Top 30) ===\n")
  print(head(var_imp, 30))

  cat("\n\n=== Risk Score Calculation ===\n")
  cat("For GBM models, risk scores are calculated using:\n")
  cat("  Risk Score = predict(gbm_model, newdata, n.trees = best_iteration)\n\n")
  cat("This is NOT a simple linear formula but a complex ensemble of decision trees.\n")
  cat("The model makes predictions by:\n")
  cat("  1. Passing input features through", model_obj$best, "decision trees\n")
  cat("  2. Each tree contributes a partial prediction\n")
  cat("  3. Final risk score = sum of all tree predictions\n\n")

  cat("To calculate risk scores for new samples:\n")
  cat("  1. Prepare a data frame with the", length(model_obj$fit$var.names), "selected features\n")
  cat("  2. Use: predict(model$fit, newdata, n.trees = model$best)\n")
  cat("  3. Higher scores indicate higher risk\n\n")

  cat("============================================================\n")
  cat("Generated on:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("============================================================\n")
  sink()

  # Save the model object for future use
  gbm_model <- list(
    model_name = model_name,
    gbm_fit = model_obj$fit,
    best_iteration = model_obj$best,
    selected_genes = model_obj$fit$var.names,
    variable_importance = var_imp,
    n_features = length(model_obj$fit$var.names)
  )

  save(gbm_model, file = "StepCox_GBM_model_object.Rdata")

  cat("\n=== Files Created ===\n")
  cat("  1. StepCox_GBM_model_details.txt - Comprehensive model description\n")
  cat("  2. StepCox_GBM_selected_genes.txt - List of selected genes\n")
  cat("  3. StepCox_GBM_variable_importance.csv - Variable importance scores\n")
  cat("  4. StepCox_GBM_model_object.Rdata - Model object for predictions\n")
}

cat("\n=== Script Completed ===\n")
