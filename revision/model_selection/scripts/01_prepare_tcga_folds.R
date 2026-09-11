#!/usr/bin/env Rscript

args_all <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args_all[grepl("^--file=", args_all)])
script_path <- normalizePath(script_arg[[1]], mustWork = TRUE)
analysis_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
source(file.path(analysis_root, "scripts", "common.R"))

dir.create(processed_root, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_root, recursive = TRUE, showWarnings = FALSE)

expression_file <- file.path(raw_root, "TCGA_unclipped_448_expression.csv")
clinical_file <- file.path(raw_root, "TCGA-LIHC_frozen_risks.csv")
if (!file.exists(expression_file) || !file.exists(clinical_file)) stop("Frozen TCGA inputs are missing")

expr <- read.csv(expression_file, check.names = FALSE)
clin <- read.csv(clinical_file, check.names = FALSE)
stopifnot(names(expr)[1] == "ID", all(c("ID", "OS.time", "OS") %in% names(clin)))
if (anyDuplicated(expr$ID) || anyDuplicated(clin$ID)) stop("Duplicate TCGA patient IDs")
if (!setequal(expr$ID, clin$ID)) stop("Expression and survival IDs differ")
clin <- clin[match(expr$ID, clin$ID), , drop = FALSE]
if (!identical(as.character(expr$ID), as.character(clin$ID))) stop("TCGA alignment failed")
if (any(!is.finite(clin$OS.time)) || any(clin$OS.time <= 0)) stop("Invalid survival times")
if (!all(clin$OS %in% c(0, 1))) stop("Invalid survival status")

candidate_genes <- names(expr)[-1]
if (length(candidate_genes) != 448L || anyDuplicated(candidate_genes)) {
  stop("Expected exactly 448 unique frozen candidate genes")
}

configs <- generate_configurations()
write.csv(configs, file.path(processed_root, "all_117_model_configurations.csv"), row.names = FALSE)

fold_id <- stratified_folds(clin$OS, k = outer_folds, seed = master_seed)
fold_assignment <- data.frame(
  patient_id = as.character(clin$ID),
  fold = fold_id,
  OS.time = as.numeric(clin$OS.time),
  OS = as.integer(clin$OS),
  stringsAsFactors = FALSE
)
write.csv(fold_assignment, file.path(processed_root, "TCGA_fold_assignments.csv"), row.names = FALSE)

feature_rows <- list()
balance_rows <- list()

for (fold in seq_len(outer_folds)) {
  held_idx <- which(fold_id == fold)
  train_idx <- which(fold_id != fold)
  x_train_raw <- expr[train_idx, candidate_genes, drop = FALSE]
  x_held_raw <- expr[held_idx, candidate_genes, drop = FALSE]

  prep <- fit_preprocessor(x_train_raw)
  x_train <- apply_preprocessor(x_train_raw, prep)
  x_held <- apply_preprocessor(x_held_raw, prep)
  screen <- univariable_cox_screen(
    x_train,
    time = clin$OS.time[train_idx],
    status = clin$OS[train_idx],
    p_cutoff = univariable_p_cutoff
  )
  if (length(screen$selected) < 2L) stop("Fold ", fold, " retained fewer than two genes")

  train_data <- data.frame(
    ID = as.character(clin$ID[train_idx]),
    OS.time = as.numeric(clin$OS.time[train_idx]),
    OS = as.integer(clin$OS[train_idx]),
    x_train[, screen$selected, drop = FALSE],
    check.names = FALSE
  )
  heldout_data <- data.frame(
    ID = as.character(clin$ID[held_idx]),
    OS.time = as.numeric(clin$OS.time[held_idx]),
    OS = as.integer(clin$OS[held_idx]),
    x_held[, screen$selected, drop = FALSE],
    check.names = FALSE
  )

  fold_bundle <- list(
    fold = fold,
    master_seed = master_seed,
    expression_sha256 = sha256_file(expression_file),
    clinical_sha256 = sha256_file(clinical_file),
    candidate_genes = candidate_genes,
    preprocessor = prep,
    univariable_results = screen$table,
    selected_genes = screen$selected,
    train_data = train_data,
    heldout_data = heldout_data
  )
  saveRDS(fold_bundle, file.path(processed_root, sprintf("TCGA_fold%d_input.rds", fold)), compress = "xz")
  write.csv(screen$table, file.path(processed_root, sprintf("TCGA_fold%d_univariable_cox.csv", fold)), row.names = FALSE)

  feature_rows[[fold]] <- data.frame(
    fold = fold,
    training_patients = length(train_idx),
    heldout_patients = length(held_idx),
    candidates_before_preprocessing = length(candidate_genes),
    removed_near_zero_variance = length(prep$removed),
    candidates_entering_univariable_cox = length(prep$keep),
    genes_after_univariable_cox_p_lt_0_01 = length(screen$selected),
    selected_genes = paste(screen$selected, collapse = ";"),
    stringsAsFactors = FALSE
  )
  balance_rows[[fold]] <- data.frame(
    fold = fold,
    partition = c("training", "heldout"),
    n = c(length(train_idx), length(held_idx)),
    events = c(sum(clin$OS[train_idx] == 1), sum(clin$OS[held_idx] == 1)),
    censored = c(sum(clin$OS[train_idx] == 0), sum(clin$OS[held_idx] == 0)),
    event_fraction = c(mean(clin$OS[train_idx] == 1), mean(clin$OS[held_idx] == 1)),
    stringsAsFactors = FALSE
  )
}

feature_table <- do.call(rbind, feature_rows)
balance_table <- do.call(rbind, balance_rows)
write.csv(feature_table, file.path(analysis_root, "feature_selection_by_fold.csv"), row.names = FALSE)
write.csv(balance_table, file.path(processed_root, "TCGA_fold_balance.csv"), row.names = FALSE)

input_manifest <- data.frame(
  role = c("TCGA expression", "TCGA survival/outcome"),
  source_id = c(basename(expression_file), basename(clinical_file)),
  sha256 = c(sha256_file(expression_file), sha256_file(clinical_file)),
  n_patients = c(nrow(expr), nrow(clin)),
  stringsAsFactors = FALSE
)
write.csv(input_manifest, file.path(processed_root, "TCGA_input_manifest.csv"), row.names = FALSE)

session <- capture.output(sessionInfo())
writeLines(c(
  paste0("Prepared at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("Master seed: ", master_seed),
  paste0("TCGA patients: ", nrow(expr)),
  paste0("TCGA events: ", sum(clin$OS == 1)),
  paste0("Frozen candidate genes: ", length(candidate_genes)),
  "No ICGC file is referenced or opened by this script.",
  "",
  session
), file.path(logs_root, "01_prepare_tcga_folds.log"))

cat("Prepared 5 stratified TCGA folds with train-only preprocessing and univariable Cox screening.\n")
print(feature_table[, setdiff(names(feature_table), "selected_genes")])
print(balance_table)
