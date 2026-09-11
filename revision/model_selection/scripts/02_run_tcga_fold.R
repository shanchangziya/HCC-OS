#!/usr/bin/env Rscript

args_all <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args_all[grepl("^--file=", args_all)])
script_path <- normalizePath(script_arg[[1]], mustWork = TRUE)
analysis_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
source(file.path(analysis_root, "scripts", "common.R"))

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("Usage: 02_run_tcga_fold.R <fold 1-5>")
fold <- as.integer(args[[1]])
if (!is.finite(fold) || !fold %in% seq_len(outer_folds)) stop("Invalid fold")

input_file <- file.path(processed_root, sprintf("TCGA_fold%d_input.rds", fold))
if (!file.exists(input_file)) stop("Run 01_prepare_tcga_folds.R first")
bundle <- readRDS(input_file)
configs <- generate_configurations()

requested <- Sys.getenv("OSARS_CONFIG_IDS", unset = "")
if (nzchar(requested)) {
  ids <- as.integer(strsplit(requested, ",", fixed = TRUE)[[1]])
  ids <- ids[is.finite(ids) & ids %in% configs$configuration_id]
  if (!length(ids)) stop("OSARS_CONFIG_IDS did not contain valid IDs")
  configs <- configs[configs$configuration_id %in% ids, , drop = FALSE]
}

checkpoint_dir <- file.path(processed_root, sprintf("fold%d_checkpoints", fold))
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(logs_root, sprintf("02_fold%d_models.log", fold))

for (row_index in seq_len(nrow(configs))) {
  config <- configs[row_index, , drop = FALSE]
  checkpoint <- file.path(checkpoint_dir, sprintf("config_%03d.rds", config$configuration_id))
  if (file.exists(checkpoint)) {
    cat(sprintf("Fold %d | %03d/117 | cached | %s\n", fold, config$configuration_id, config$configuration))
    next
  }

  cat(sprintf("Fold %d | %03d/117 | start  | %s\n", fold, config$configuration_id, config$configuration))
  started <- Sys.time()
  result <- tryCatch(
    {
      out <- run_mime_configuration(
        config = config,
        train_data = bundle$train_data,
        heldout_data = bundle$heldout_data,
        candidate_genes = bundle$selected_genes,
        log_file = log_file,
        seed = master_seed,
        gbm_cores = 1L
      )
      out$success <- TRUE
      out$error <- NA_character_
      out
    },
    error = function(e) {
      list(
        configuration_id = config$configuration_id,
        configuration = config$configuration,
        success = FALSE,
        error = conditionMessage(e),
        fold_cindex = NA_real_,
        training_cindex = NA_real_,
        orientation = NA_real_,
        heldout = data.frame(
          patient_id = bundle$heldout_data$ID,
          OS.time = bundle$heldout_data$OS.time,
          OS = bundle$heldout_data$OS,
          risk_score = NA_real_,
          raw_risk_score = NA_real_,
          stringsAsFactors = FALSE
        ),
        input_features = character(),
        input_feature_count = NA_integer_,
        univariable_feature_count = length(bundle$selected_genes),
        parameter_string = NA_character_,
        warnings = character(),
        elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
        upstream_model_name = NA_character_
      )
    }
  )
  result$fold <- fold
  result$completed_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  saveRDS(result, checkpoint, compress = "xz")
  cat(sprintf(
    "Fold %d | %03d/117 | %-6s | C=%s | %.1fs\n",
    fold,
    config$configuration_id,
    ifelse(result$success, "done", "failed"),
    ifelse(is.finite(result$fold_cindex), sprintf("%.4f", result$fold_cindex), "NA"),
    result$elapsed_seconds
  ))
  if (!result$success) cat("  Error: ", result$error, "\n", sep = "")
  rm(result)
  invisible(gc())
}

checkpoint_files <- list.files(checkpoint_dir, pattern = "^config_[0-9]{3}\\.rds$", full.names = TRUE)
records <- lapply(checkpoint_files, readRDS)
summary_rows <- lapply(records, function(x) data.frame(
  fold = x$fold,
  configuration_id = x$configuration_id,
  configuration = x$configuration,
  success = x$success,
  Cindex = x$fold_cindex,
  training_Cindex = x$training_cindex,
  orientation = x$orientation,
  univariable_feature_count = x$univariable_feature_count,
  final_input_feature_count = x$input_feature_count,
  elapsed_seconds = x$elapsed_seconds,
  parameter_string = x$parameter_string,
  error = x$error,
  warning_count = length(x$warnings),
  stringsAsFactors = FALSE
))
fold_summary <- do.call(rbind, summary_rows)
fold_summary <- fold_summary[order(fold_summary$configuration_id), , drop = FALSE]
write.csv(fold_summary, file.path(processed_root, sprintf("TCGA_fold%d_model_summary.csv", fold)), row.names = FALSE)

cat(sprintf("Fold %d checkpointed configurations: %d; successful: %d.\n", fold, nrow(fold_summary), sum(fold_summary$success)))
