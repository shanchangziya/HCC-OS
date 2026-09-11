#!/usr/bin/env Rscript

args_all <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args_all[grepl("^--file=", args_all)])
script_path <- normalizePath(script_arg[[1]], mustWork = TRUE)
analysis_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
source(file.path(analysis_root, "scripts", "common.R"))

selection_file <- file.path(processed_root, "TCGA_CV_model_selection.rds")
if (!file.exists(selection_file)) stop("Run 03_aggregate_select_and_lock.R first")
selection <- readRDS(selection_file)
configs <- generate_configurations()
config <- configs[configs$configuration_id == selection$selected_id, , drop = FALSE]
if (nrow(config) != 1L) stop("Selected configuration is not unique")

expression_file <- file.path(raw_root, "TCGA_unclipped_448_expression.csv")
clinical_file <- file.path(raw_root, "TCGA-LIHC_frozen_risks.csv")
expr <- read.csv(expression_file, check.names = FALSE)
clin <- read.csv(clinical_file, check.names = FALSE)
clin <- clin[match(expr$ID, clin$ID), , drop = FALSE]
if (anyNA(clin$ID) || !identical(as.character(expr$ID), as.character(clin$ID))) stop("TCGA ID alignment failed")
candidate_genes <- names(expr)[-1]

prep <- fit_preprocessor(expr[, candidate_genes, drop = FALSE])
x_full <- apply_preprocessor(expr[, candidate_genes, drop = FALSE], prep)
screen <- univariable_cox_screen(x_full, clin$OS.time, clin$OS, p_cutoff = univariable_p_cutoff)
if (length(screen$selected) < 2L) stop("Full TCGA screen retained fewer than two genes")

full_data <- data.frame(
  ID = as.character(clin$ID),
  OS.time = as.numeric(clin$OS.time),
  OS = as.integer(clin$OS),
  x_full[, screen$selected, drop = FALSE],
  check.names = FALSE
)

fit_result <- run_mime_configuration(
  config = config,
  train_data = full_data,
  heldout_data = full_data,
  candidate_genes = screen$selected,
  log_file = file.path(logs_root, "04_final_TCGA_fit.log"),
  seed = master_seed,
  gbm_cores = 1L,
  keep_model = TRUE
)

model_data <- full_data[, -1, drop = FALSE]
names(model_data) <- gsub("-", ".", names(model_data), fixed = TRUE)
model_input_features <- fit_result$input_features
if (!all(model_input_features %in% names(model_data))) stop("Unable to reconstruct final model matrix")
training_model_data <- model_data[, c("OS.time", "OS", model_input_features), drop = FALSE]

bundle <- list(
  model_name = "OSARS",
  model_status = "locked_after_TCGA_only_model_selection",
  configuration_id = config$configuration_id,
  configuration = config$configuration,
  configuration_row = config,
  master_seed = master_seed,
  outer_cv_folds = outer_folds,
  univariable_p_cutoff = univariable_p_cutoff,
  preprocessing = "training-derived 1st/99th percentile winsorization, mean imputation, and near-zero-variance removal",
  preprocessor = prep,
  frozen_candidate_genes = candidate_genes,
  univariable_selected_genes = screen$selected,
  model_input_features = model_input_features,
  orientation = fit_result$orientation,
  fit = fit_result$model_result$ml.res[[1]],
  upstream_model_name = fit_result$upstream_model_name,
  parameter_string = fit_result$parameter_string,
  training_model_data = training_model_data,
  TCGA_expression_sha256 = sha256_file(expression_file),
  TCGA_clinical_sha256 = sha256_file(clinical_file),
  selected_mean_5fold_Cindex = selection$selected_mean_Cindex,
  selection_rule = selection$selection_rule,
  package_versions = vapply(required_model_packages, function(p) as.character(utils::packageVersion(p)), character(1)),
  locked_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
)

locked_file <- file.path(processed_root, "final_OSARS_model_locked.rds")
saveRDS(bundle, locked_file, compress = "xz")

reproduced_risk <- predict_locked_bundle(bundle, expr[, candidate_genes, drop = FALSE])
fit_risk <- fit_result$heldout$risk_score[match(expr$ID, fit_result$heldout$patient_id)]
max_prediction_delta <- max(abs(reproduced_risk - fit_risk), na.rm = TRUE)
if (!is.finite(max_prediction_delta) || max_prediction_delta > 1e-7) {
  stop("Locked prediction method failed reproduction check; max delta = ", max_prediction_delta)
}

tcga_scores <- data.frame(
  patient_id = expr$ID,
  OS.time = clin$OS.time,
  OS = clin$OS,
  locked_risk_score = reproduced_risk,
  evaluation_type = "apparent performance after refitting on complete TCGA; not CV performance",
  stringsAsFactors = FALSE
)
write.csv(tcga_scores, file.path(processed_root, "final_model_TCGA_apparent_scores.csv"), row.names = FALSE)

coef_vector <- extract_linear_coefficients(bundle)
if (!is.null(coef_vector)) {
  write.csv(data.frame(gene = names(coef_vector), coefficient = as.numeric(coef_vector)),
            file.path(processed_root, "final_model_coefficients.csv"), row.names = FALSE)
}

write.csv(screen$table, file.path(processed_root, "final_full_TCGA_univariable_cox.csv"), row.names = FALSE)
write.csv(data.frame(
  cohort = "complete TCGA development cohort",
  patients = nrow(expr),
  events = sum(clin$OS == 1),
  candidates_before_preprocessing = length(candidate_genes),
  removed_near_zero_variance = length(prep$removed),
  genes_after_univariable_cox_p_lt_0_01 = length(screen$selected),
  final_model_input_features = length(model_input_features),
  stringsAsFactors = FALSE
), file.path(processed_root, "final_full_TCGA_feature_selection.csv"), row.names = FALSE)

locked_hash <- sha256_file(locked_file)
manifest_file <- file.path(processed_root, "MODEL_LOCKED.sha256")
writeLines(paste(locked_hash, basename(locked_file)), manifest_file)
write.csv(data.frame(
  item = c("locked_model", "selected_configuration", "selection_basis", "external_data_access_before_lock", "prediction_reproduction_max_abs_delta"),
  value = c(
    locked_hash,
    config$configuration,
    "mean 5-fold cross-validation within TCGA",
    "none",
    format(max_prediction_delta, scientific = TRUE)
  ),
  stringsAsFactors = FALSE
), file.path(processed_root, "final_model_lock_manifest.csv"), row.names = FALSE)

perf <- read.csv(file.path(analysis_root, "all_117_models_5fold_performance.csv"), check.names = FALSE)
selected_perf <- perf[perf$configuration_id == config$configuration_id, , drop = FALSE]
formula_text <- if (is.null(coef_vector)) {
  paste0(
    "This is a non-linear or component-based predictor. The exact score is generated by ",
    "predict_locked_bundle() using the serialized model, the locked TCGA preprocessing vectors, ",
    "the recorded model-input feature order, and orientation=", bundle$orientation, "."
  )
} else {
  kind <- if (config$mode == "single") config$single_ml else config$double_ml2
  transform <- if (kind %in% c("Lasso", "Ridge", "StepCox")) "exp(sum(beta_g * x_g))" else "sum(beta_g * x_g)"
  paste0(
    "risk = orientation * ", transform, ", where orientation=", bundle$orientation,
    "; x_g is processed with the locked TCGA q01/q99/mean values; coefficients are in final_model_coefficients.csv."
  )
}

summary_lines <- c(
  "OSARS FINAL MODEL (LOCKED AFTER TCGA-ONLY SELECTION)",
  paste0("Final algorithm/configuration: ", config$configuration),
  paste0("Selection metric: mean 5-fold cross-validation within TCGA = ", sprintf("%.6f", selected_perf$mean_Cindex)),
  paste0("Fold C-indices: ", paste(sprintf("%.6f", unlist(selected_perf[paste0("fold", 1:5, "_Cindex")])), collapse = ", ")),
  paste0("Fold SD: ", sprintf("%.6f", selected_perf$SD)),
  paste0("Fold median: ", sprintf("%.6f", selected_perf$median_Cindex)),
  paste0("95% CI across fold estimates: ", sprintf("%.6f to %.6f", selected_perf$CI95_low, selected_perf$CI95_high)),
  paste0("Selection rule: ", selection$selection_rule),
  paste0("Full-TCGA candidate genes: ", length(candidate_genes)),
  paste0("Full-TCGA genes after univariable Cox P < 0.01: ", length(screen$selected)),
  paste0("Final model input features: ", length(model_input_features)),
  paste0("Parameters: ", fit_result$parameter_string),
  paste0("Risk-score calculation: ", formula_text),
  paste0("Locked model: ", normalizePath(locked_file)),
  paste0("Locked model SHA256: ", locked_hash),
  paste0("Prediction reproduction max absolute delta: ", format(max_prediction_delta, scientific = TRUE)),
  paste0("Complete-TCGA apparent C-index (descriptive only): ", sprintf("%.6f", harrell_c(clin$OS.time, clin$OS, reproduced_risk))),
  "No external cohort was opened by this script."
)
writeLines(summary_lines, file.path(analysis_root, "final_model_summary.txt"))

cat(paste(summary_lines, collapse = "\n"), "\n")
