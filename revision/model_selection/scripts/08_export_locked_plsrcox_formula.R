#!/usr/bin/env Rscript

args_all <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args_all[grepl("^--file=", args_all)])
script_path <- normalizePath(script_arg[[1]], mustWork = TRUE)
analysis_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
source(file.path(analysis_root, "scripts", "common.R"))

locked_file <- file.path(processed_root, "final_OSARS_model_locked.rds")
lock_file <- file.path(processed_root, "MODEL_LOCKED.sha256")
if (!file.exists(locked_file) || !file.exists(lock_file)) stop("Locked model files are missing")
expected_hash <- strsplit(readLines(lock_file, n = 1L), "[[:space:]]+")[[1]][1]
if (!identical(expected_hash, sha256_file(locked_file))) stop("Locked model hash verification failed")

bundle <- readRDS(locked_file)
config <- bundle$configuration_row
final_kind <- if (config$mode == "single") config$single_ml else config$double_ml2
if (!identical(final_kind, "plsRcox") || !inherits(bundle$fit, "plsRcoxmodel")) {
  stop("This exporter is only applicable when the locked final predictor is plsRcox")
}

fit <- bundle$fit
features <- bundle$model_input_features
k <- as.integer(fit$computed_nt)
weights <- as.matrix(fit$wwetoile[, seq_len(k), drop = FALSE])
if (nrow(weights) != length(features)) stop("PLS weight rows do not match locked features")
rownames(weights) <- features
colnames(weights) <- paste0("component_", seq_len(k))

pls_center <- attr(fit$ExpliX, "scaled:center")
pls_scale <- attr(fit$ExpliX, "scaled:scale")
if (length(pls_center) != length(features) || length(pls_scale) != length(features)) {
  stop("PLS centering/scaling vectors do not match locked features")
}
names(pls_center) <- names(pls_scale) <- features

original_names <- names(bundle$preprocessor$q01)
normalized_names <- gsub("_", ".", gsub("-", ".", original_names, fixed = TRUE), fixed = TRUE)
feature_index <- match(features, normalized_names)
if (anyNA(feature_index)) stop("Unable to map model features to locked preprocessing vectors")

feature_table <- data.frame(
  gene = features,
  winsor_q01 = unname(bundle$preprocessor$q01[feature_index]),
  winsor_q99 = unname(bundle$preprocessor$q99[feature_index]),
  missing_value_imputation_mean = unname(bundle$preprocessor$means[feature_index]),
  pls_center = unname(pls_center[features]),
  pls_scale = unname(pls_scale[features]),
  stringsAsFactors = FALSE
)
write.csv(feature_table, file.path(processed_root, "final_model_feature_preprocessing.csv"), row.names = FALSE)
write.csv(
  data.frame(gene = rownames(weights), weights, check.names = FALSE),
  file.path(processed_root, "final_model_pls_component_weights.csv"),
  row.names = FALSE
)

beta <- as.numeric(fit$FinalModel$coefficients[seq_len(k)])
reference_mean <- as.numeric(fit$FinalModel$means[seq_len(k)])
component_table <- data.frame(
  component = paste0("component_", seq_len(k)),
  cox_coefficient = beta,
  cox_reference_mean = reference_mean,
  stringsAsFactors = FALSE
)
write.csv(component_table, file.path(processed_root, "final_model_pls_cox_coefficients.csv"), row.names = FALSE)

manual_score <- function(raw_expression) {
  x <- prepare_locked_model_matrix(bundle, raw_expression)
  z <- sweep(sweep(as.matrix(x), 2, pls_center), 2, pls_scale, "/")
  component_scores <- z %*% weights
  as.numeric(bundle$orientation * (component_scores %*% beta - sum(reference_mean * beta)))
}

validation_rows <- list()
tcga_file <- file.path(raw_root, "TCGA_unclipped_448_expression.csv")
icgc_file <- file.path(raw_root, "ICGC_original_448_expression.csv")
cohort_files <- c(TCGA = tcga_file, ICGC = icgc_file)
for (cohort in names(cohort_files)) {
  path <- unname(cohort_files[[cohort]])
  if (!file.exists(path)) next
  dat <- read.csv(path, check.names = FALSE)
  raw_expression <- dat[, bundle$frozen_candidate_genes, drop = FALSE]
  serialized_score <- predict_locked_bundle(bundle, raw_expression)
  explicit_score <- manual_score(raw_expression)
  validation_rows[[cohort]] <- data.frame(
    cohort = cohort,
    n = length(serialized_score),
    max_absolute_difference = max(abs(serialized_score - explicit_score)),
    pearson_correlation = stats::cor(serialized_score, explicit_score),
    stringsAsFactors = FALSE
  )
}
validation_table <- do.call(rbind, validation_rows)
write.csv(validation_table, file.path(processed_root, "final_model_formula_validation.csv"), row.names = FALSE)

formula_lines <- c(
  "EXPLICIT LOCKED PLS-COX RISK-SCORE CALCULATION",
  paste0("K = ", k, " latent components; orientation = ", bundle$orientation, "."),
  "For each locked model gene g, cap expression x_g at the stored TCGA q01/q99 bounds and replace missing values with the stored TCGA mean.",
  "Standardize each processed gene: z_g = (x_g - pls_center_g) / pls_scale_g.",
  "Calculate latent components: t_k = sum_g(z_g * W_gk).",
  "Calculate risk score: risk = orientation * sum_k[beta_k * (t_k - cox_reference_mean_k)].",
  "No cohort-specific re-estimation, re-scaling, or exponential transformation is applied.",
  "Parameters: data/processed/final_model_feature_preprocessing.csv; final_model_pls_component_weights.csv; final_model_pls_cox_coefficients.csv.",
  paste0(
    "Formula-versus-serialized prediction max absolute differences: ",
    paste(validation_table$cohort, format(validation_table$max_absolute_difference, scientific = TRUE), sep = "=", collapse = "; "), "."
  ),
  paste0("Locked model SHA256 (unchanged): ", expected_hash)
)
writeLines(formula_lines, file.path(analysis_root, "final_model_risk_formula.txt"))

summary_file <- file.path(analysis_root, "final_model_summary.txt")
summary_lines <- readLines(summary_file, warn = FALSE)
summary_lines <- summary_lines[!grepl("^Explicit PLS-Cox calculation:", summary_lines)]
summary_lines <- c(
  summary_lines,
  paste0(
    "Explicit PLS-Cox calculation: K=", k,
    "; z_g=(processed x_g-center_g)/scale_g; t_k=sum_g(z_g*W_gk); ",
    "risk=orientation*sum_k[beta_k*(t_k-reference_mean_k)]. See final_model_risk_formula.txt and the three parameter CSV files."
  )
)
writeLines(summary_lines, summary_file)

if (!identical(expected_hash, sha256_file(locked_file))) stop("Locked model changed during formula export")
cat(paste(formula_lines, collapse = "\n"), "\n")
