#!/usr/bin/env Rscript

args_all <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args_all[grepl("^--file=", args_all)])
script_path <- normalizePath(script_arg[[1]], mustWork = TRUE)
analysis_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
source(file.path(analysis_root, "scripts", "common.R"))

configs <- generate_configurations()
all_records <- list()
record_index <- 0L
missing_checkpoints <- character()
for (fold in seq_len(outer_folds)) {
  for (config_id in configs$configuration_id) {
    path <- file.path(processed_root, sprintf("fold%d_checkpoints/config_%03d.rds", fold, config_id))
    if (!file.exists(path)) {
      missing_checkpoints <- c(missing_checkpoints, sprintf("fold%d/config_%03d", fold, config_id))
      next
    }
    record_index <- record_index + 1L
    all_records[[record_index]] <- readRDS(path)
  }
}
if (length(missing_checkpoints)) {
  stop("Missing ", length(missing_checkpoints), " model checkpoints; first missing: ", paste(head(missing_checkpoints, 10), collapse = ", "))
}
stopifnot(length(all_records) == 117L * outer_folds)

fold_rows <- lapply(all_records, function(x) data.frame(
  fold = x$fold,
  configuration_id = x$configuration_id,
  configuration = x$configuration,
  success = isTRUE(x$success),
  Cindex = x$fold_cindex,
  training_Cindex = x$training_cindex,
  orientation = x$orientation,
  final_input_feature_count = x$input_feature_count,
  univariable_feature_count = x$univariable_feature_count,
  elapsed_seconds = x$elapsed_seconds,
  parameter_string = x$parameter_string,
  error = x$error,
  warning_count = length(x$warnings),
  stringsAsFactors = FALSE
))
fold_table <- do.call(rbind, fold_rows)
fold_table <- fold_table[order(fold_table$configuration_id, fold_table$fold), , drop = FALSE]
write.csv(fold_table, file.path(processed_root, "all_model_fold_level_results.csv"), row.names = FALSE)

prediction_rows <- lapply(all_records, function(x) {
  tab <- x$heldout
  data.frame(
    patient_id = as.character(tab$patient_id),
    fold = x$fold,
    OS.time = as.numeric(tab$OS.time),
    OS = as.integer(tab$OS),
    configuration_id = x$configuration_id,
    configuration = x$configuration,
    risk_score = as.numeric(tab$risk_score),
    raw_risk_score = as.numeric(tab$raw_risk_score),
    model_success = isTRUE(x$success),
    stringsAsFactors = FALSE
  )
})
predictions <- do.call(rbind, prediction_rows)
predictions <- predictions[order(predictions$configuration_id, predictions$fold, predictions$patient_id), , drop = FALSE]
write.csv(predictions, file.path(analysis_root, "cv_patient_predictions.csv"), row.names = FALSE)

expected_patients <- length(unique(predictions$patient_id))
metrics <- lapply(seq_len(nrow(configs)), function(i) {
  cfg <- configs[i, ]
  d <- fold_table[fold_table$configuration_id == cfg$configuration_id, , drop = FALSE]
  cvals <- d$Cindex[order(d$fold)]
  n_ok <- sum(d$success & is.finite(d$Cindex))
  complete_predictions <- predictions[
    predictions$configuration_id == cfg$configuration_id & is.finite(predictions$risk_score),
    , drop = FALSE
  ]
  is_complete <- n_ok == outer_folds && length(unique(complete_predictions$patient_id)) == expected_patients &&
    nrow(complete_predictions) == expected_patients
  m <- if (n_ok) mean(cvals, na.rm = TRUE) else NA_real_
  s <- if (n_ok > 1L) stats::sd(cvals, na.rm = TRUE) else NA_real_
  med <- if (n_ok) stats::median(cvals, na.rm = TRUE) else NA_real_
  half <- if (n_ok > 1L) stats::qt(0.975, df = n_ok - 1L) * s / sqrt(n_ok) else NA_real_
  fcount <- d$final_input_feature_count[order(d$fold)]
  data.frame(
    configuration_id = cfg$configuration_id,
    configuration = cfg$configuration,
    stage_count = cfg$stage_count,
    fold1_Cindex = cvals[1],
    fold2_Cindex = cvals[2],
    fold3_Cindex = cvals[3],
    fold4_Cindex = cvals[4],
    fold5_Cindex = cvals[5],
    mean_Cindex = m,
    SD = s,
    median_Cindex = med,
    CI95_low = if (is.finite(half)) max(0, m - half) else NA_real_,
    CI95_high = if (is.finite(half)) min(1, m + half) else NA_real_,
    CI_method = "two-sided t interval across five held-out-fold C-indices",
    successful_folds = n_ok,
    five_fold_predictions_complete = is_complete,
    fold1_n_features = fcount[1],
    fold2_n_features = fcount[2],
    fold3_n_features = fcount[3],
    fold4_n_features = fcount[4],
    fold5_n_features = fcount[5],
    mean_n_features = mean(fcount, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
})
performance <- do.call(rbind, metrics)

eligible <- which(performance$five_fold_predictions_complete & is.finite(performance$mean_Cindex))
if (!length(eligible)) stop("No configuration completed all five folds")
selection_order <- eligible[order(
  -performance$mean_Cindex[eligible],
  performance$mean_n_features[eligible],
  performance$stage_count[eligible],
  performance$configuration_id[eligible]
)]
selected_index <- selection_order[1]
selected_id <- performance$configuration_id[selected_index]
selected_configuration <- performance$configuration[selected_index]
best_mean <- performance$mean_Cindex[selected_index]

performance$rank <- NA_integer_
performance$rank[eligible] <- rank(-performance$mean_Cindex[eligible], ties.method = "min")
performance$delta_from_best <- best_mean - performance$mean_Cindex
performance$within_0_01_of_best <- performance$five_fold_predictions_complete &
  is.finite(performance$delta_from_best) & performance$delta_from_best < 0.01
performance$selected_final_configuration <- performance$configuration_id == selected_id

# Audit exact and monotonic duplicate out-of-fold predictors within every fold.
pair_rows <- list()
pair_index <- 0L
complete_ids <- performance$configuration_id[performance$five_fold_predictions_complete]
for (ii in seq_len(length(complete_ids) - 1L)) {
  for (jj in (ii + 1L):length(complete_ids)) {
    id1 <- complete_ids[ii]
    id2 <- complete_ids[jj]
    exact_by_fold <- logical(outer_folds)
    rank_by_fold <- logical(outer_folds)
    max_delta <- numeric(outer_folds)
    min_rho <- numeric(outer_folds)
    for (fold in seq_len(outer_folds)) {
      a <- predictions[predictions$configuration_id == id1 & predictions$fold == fold, c("patient_id", "risk_score")]
      b <- predictions[predictions$configuration_id == id2 & predictions$fold == fold, c("patient_id", "risk_score")]
      b <- b[match(a$patient_id, b$patient_id), , drop = FALSE]
      av <- a$risk_score
      bv <- b$risk_score
      exact_by_fold[fold] <- isTRUE(all.equal(av, bv, tolerance = 1e-10, check.attributes = FALSE))
      max_delta[fold] <- max(abs(av - bv), na.rm = TRUE)
      rho <- suppressWarnings(stats::cor(av, bv, method = "spearman", use = "complete.obs"))
      min_rho[fold] <- rho
      rank_by_fold[fold] <- is.finite(rho) && abs(rho - 1) < 1e-10
    }
    if (all(exact_by_fold) || all(rank_by_fold)) {
      pair_index <- pair_index + 1L
      pair_rows[[pair_index]] <- data.frame(
        configuration_1 = performance$configuration[match(id1, performance$configuration_id)],
        configuration_2 = performance$configuration[match(id2, performance$configuration_id)],
        exact_predictions_in_all_folds = all(exact_by_fold),
        rank_equivalent_in_all_folds = all(rank_by_fold),
        max_absolute_difference = max(max_delta),
        minimum_fold_spearman_rho = min(min_rho),
        stringsAsFactors = FALSE
      )
    }
  }
}
duplicate_pairs <- if (length(pair_rows)) do.call(rbind, pair_rows) else data.frame(
  configuration_1 = character(), configuration_2 = character(),
  exact_predictions_in_all_folds = logical(), rank_equivalent_in_all_folds = logical(),
  max_absolute_difference = numeric(), minimum_fold_spearman_rho = numeric()
)

assign_groups <- function(ids, pairs, exact_only = TRUE) {
  parent <- stats::setNames(ids, ids)
  find_root <- function(x) {
    while (parent[[as.character(x)]] != x) x <- parent[[as.character(x)]]
    x
  }
  union_ids <- function(a, b) {
    ra <- find_root(a); rb <- find_root(b)
    if (ra != rb) parent[[as.character(rb)]] <<- ra
  }
  if (nrow(pairs)) {
    use <- if (exact_only) pairs$exact_predictions_in_all_folds else pairs$rank_equivalent_in_all_folds
    for (i in which(use)) {
      a <- performance$configuration_id[match(pairs$configuration_1[i], performance$configuration)]
      b <- performance$configuration_id[match(pairs$configuration_2[i], performance$configuration)]
      union_ids(a, b)
    }
  }
  roots <- vapply(ids, find_root, integer(1))
  counts <- table(roots)
  group <- rep(NA_character_, length(ids))
  duplicated_roots <- as.integer(names(counts[counts > 1L]))
  for (i in seq_along(duplicated_roots)) group[roots == duplicated_roots[i]] <- sprintf("D%02d", i)
  stats::setNames(group, ids)
}

exact_groups <- assign_groups(performance$configuration_id, duplicate_pairs, exact_only = TRUE)
rank_groups <- assign_groups(performance$configuration_id, duplicate_pairs, exact_only = FALSE)
performance$exact_duplicate_group <- unname(exact_groups[as.character(performance$configuration_id)])
performance$rank_equivalent_group <- unname(rank_groups[as.character(performance$configuration_id)])

performance <- performance[order(is.na(performance$rank), performance$rank, performance$configuration_id), , drop = FALSE]
write.csv(performance, file.path(analysis_root, "all_117_models_5fold_performance.csv"), row.names = FALSE)
write.csv(duplicate_pairs, file.path(processed_root, "duplicate_predictor_audit.csv"), row.names = FALSE)
write.csv(performance[performance$within_0_01_of_best, ], file.path(processed_root, "models_within_0.01_of_best.csv"), row.names = FALSE)

feature_count_long <- fold_table[, c(
  "fold", "configuration_id", "configuration", "univariable_feature_count", "final_input_feature_count"
)]
write.csv(feature_count_long, file.path(processed_root, "model_feature_counts_by_fold.csv"), row.names = FALSE)

selection <- list(
  selected_id = selected_id,
  selected_configuration = selected_configuration,
  selected_mean_Cindex = best_mean,
  selected_row = performance[performance$configuration_id == selected_id, , drop = FALSE],
  selection_rule = paste(
    "Highest mean held-out-fold C-index among configurations with complete five-fold predictions;",
    "exact numerical ties are resolved by lower mean input-feature count, then single-stage over two-stage,",
    "then the frozen configuration order. ICGC is not consulted."
  ),
  completed_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
)
saveRDS(selection, file.path(processed_root, "TCGA_CV_model_selection.rds"), compress = "xz")

cat("Selected configuration: ", selected_configuration, "\n", sep = "")
cat("Mean 5-fold C-index: ", sprintf("%.6f", best_mean), "\n", sep = "")
cat("Complete models: ", sum(performance$five_fold_predictions_complete), "/117\n", sep = "")
cat("Models within <0.01: ", sum(performance$within_0_01_of_best), "\n", sep = "")
cat("Duplicate/rank-equivalent pairs: ", nrow(duplicate_pairs), "\n", sep = "")
