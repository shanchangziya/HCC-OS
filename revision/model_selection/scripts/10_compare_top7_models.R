#!/usr/bin/env Rscript

args_all <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args_all[grepl("^--file=", args_all)])
script_path <- normalizePath(script_arg[[1]], mustWork = TRUE)
analysis_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
processed_root <- file.path(analysis_root, "data", "processed")
dir.create(processed_root, recursive = TRUE, showWarnings = FALSE)

performance <- read.csv(
  file.path(analysis_root, "all_117_models_5fold_performance.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
fold_columns <- paste0("fold", seq_len(5L), "_Cindex")
required <- c(
  "configuration_id", "configuration", "rank", fold_columns,
  "mean_Cindex", "SD", "median_Cindex", "CI95_low", "CI95_high",
  "successful_folds", "five_fold_predictions_complete", "mean_n_features",
  "delta_from_best", "within_0_01_of_best", "exact_duplicate_group"
)
if (!all(required %in% names(performance))) {
  stop("Missing required columns: ", paste(setdiff(required, names(performance)), collapse = ", "))
}

osars_configuration <- "StepCox[forward] + GBM"
performance$is_OSARS <- performance$configuration == osars_configuration
if (sum(performance$is_OSARS) != 1L) stop("OSARS configuration is not unique")

# Publication-ready source table for the 117-configuration panel.
supplementary_table <- performance[, c(
  "configuration_id", "configuration", "rank", fold_columns,
  "mean_Cindex", "SD", "median_Cindex", "CI95_low", "CI95_high",
  "successful_folds", "five_fold_predictions_complete",
  paste0("fold", seq_len(5L), "_n_features"), "mean_n_features",
  "delta_from_best", "within_0_01_of_best", "exact_duplicate_group",
  "rank_equivalent_group", "is_OSARS"
)]
names(supplementary_table)[names(supplementary_table) == "rank"] <- "CV_rank"
names(supplementary_table)[names(supplementary_table) == "CI95_low"] <- "mean_Cindex_CI95_low"
names(supplementary_table)[names(supplementary_table) == "CI95_high"] <- "mean_Cindex_CI95_high"
write.csv(
  supplementary_table,
  file.path(analysis_root, "Supplementary_Table_117_models_TCGA_5fold_Cindex.csv"),
  row.names = FALSE
)

top7 <- performance[performance$within_0_01_of_best, , drop = FALSE]
top7 <- top7[order(-top7$mean_Cindex, top7$configuration_id), , drop = FALSE]
if (nrow(top7) != 7L) stop("Expected seven configurations within 0.01 of the best")
if (!osars_configuration %in% top7$configuration) stop("OSARS is not among the top seven")

fold_matrix <- t(as.matrix(top7[, fold_columns, drop = FALSE]))
storage.mode(fold_matrix) <- "double"
colnames(fold_matrix) <- top7$configuration
rownames(fold_matrix) <- paste0("Fold ", seq_len(5L))

top7_long <- do.call(rbind, lapply(seq_len(ncol(fold_matrix)), function(j) {
  data.frame(
    configuration = colnames(fold_matrix)[j],
    CV_rank = top7$rank[j],
    fold = seq_len(nrow(fold_matrix)),
    Cindex = fold_matrix[, j],
    is_OSARS = colnames(fold_matrix)[j] == osars_configuration,
    stringsAsFactors = FALSE
  )
}))
write.csv(top7_long, file.path(processed_root, "Top7_fold_Cindex_long.csv"), row.names = FALSE)

friedman_result <- stats::friedman.test(fold_matrix)
kendall_w <- unname(as.numeric(friedman_result$statistic)) /
  (nrow(fold_matrix) * (ncol(fold_matrix) - 1L))

paired_difference <- function(x, y) {
  difference <- as.numeric(x - y)
  n <- length(difference)
  mean_difference <- mean(difference)
  sd_difference <- stats::sd(difference)

  if (!is.finite(sd_difference) || sd_difference == 0) {
    ci <- c(mean_difference, mean_difference)
    t_statistic <- if (mean_difference == 0) 0 else sign(mean_difference) * Inf
    paired_t_p <- if (mean_difference == 0) 1 else 0
  } else {
    standard_error <- sd_difference / sqrt(n)
    half_width <- stats::qt(0.975, df = n - 1L) * standard_error
    ci <- mean_difference + c(-1, 1) * half_width
    t_statistic <- mean_difference / standard_error
    paired_t_p <- 2 * stats::pt(-abs(t_statistic), df = n - 1L)
  }

  # Exact paired sign-flip randomization test over the five held-out folds.
  sign_grid <- as.matrix(expand.grid(rep(list(c(-1, 1)), n)))
  permuted_statistics <- abs(as.numeric(sign_grid %*% difference) / n)
  observed_statistic <- abs(mean_difference)
  sign_flip_p <- mean(permuted_statistics >= observed_statistic - 1e-15)

  c(
    mean_difference = mean_difference,
    CI95_low = ci[1],
    CI95_high = ci[2],
    SD_paired_difference = sd_difference,
    paired_t_statistic = t_statistic,
    paired_t_pvalue = paired_t_p,
    exact_sign_flip_pvalue = sign_flip_p
  )
}

pairwise_rows <- list()
row_index <- 0L
for (i in seq_len(ncol(fold_matrix) - 1L)) {
  for (j in (i + 1L):ncol(fold_matrix)) {
    row_index <- row_index + 1L
    stats_row <- paired_difference(fold_matrix[, i], fold_matrix[, j])
    pairwise_rows[[row_index]] <- data.frame(
      configuration_1 = colnames(fold_matrix)[i],
      configuration_2 = colnames(fold_matrix)[j],
      difference_definition = "mean 5-fold C-index: configuration_1 minus configuration_2",
      n_paired_folds = nrow(fold_matrix),
      as.list(stats_row),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }
}
pairwise <- do.call(rbind, pairwise_rows)
pairwise$Holm_adjusted_sign_flip_pvalue_all_21 <- stats::p.adjust(
  pairwise$exact_sign_flip_pvalue,
  method = "holm"
)
pairwise$significant_after_Holm_0_05 <-
  pairwise$Holm_adjusted_sign_flip_pvalue_all_21 < 0.05
write.csv(
  pairwise,
  file.path(processed_root, "Top7_all_pairwise_Cindex_comparisons.csv"),
  row.names = FALSE
)

reference_scores <- fold_matrix[, osars_configuration]
versus_rows <- lapply(seq_len(ncol(fold_matrix)), function(j) {
  configuration <- colnames(fold_matrix)[j]
  stats_row <- paired_difference(fold_matrix[, j], reference_scores)
  data.frame(
    configuration = configuration,
    CV_rank = top7$rank[match(configuration, top7$configuration)],
    mean_5fold_Cindex = top7$mean_Cindex[match(configuration, top7$configuration)],
    comparison_reference = osars_configuration,
    difference_definition = "mean 5-fold C-index: model minus OSARS",
    n_paired_folds = nrow(fold_matrix),
    as.list(stats_row),
    is_OSARS = configuration == osars_configuration,
    exact_duplicate_of_OSARS = all(fold_matrix[, j] == reference_scores),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
})
versus_osars <- do.call(rbind, versus_rows)
nonreference <- !versus_osars$is_OSARS
versus_osars$Holm_adjusted_sign_flip_pvalue_vs_OSARS <- 1
versus_osars$Holm_adjusted_sign_flip_pvalue_vs_OSARS[nonreference] <- stats::p.adjust(
  versus_osars$exact_sign_flip_pvalue[nonreference],
  method = "holm"
)
versus_osars$significant_after_Holm_0_05 <-
  versus_osars$Holm_adjusted_sign_flip_pvalue_vs_OSARS < 0.05
write.csv(
  versus_osars,
  file.path(processed_root, "Top7_Cindex_differences_vs_OSARS.csv"),
  row.names = FALSE
)

summary_table <- data.frame(
  analysis = c(
    "Global comparison across top seven configurations",
    "Pairwise exact sign-flip comparisons across all 21 pairs",
    "Pairwise exact sign-flip comparisons versus OSARS"
  ),
  method = c(
    "Friedman rank-sum test with fold as paired block",
    "Exact paired sign-flip test; Holm correction across 21 pairs",
    "Exact paired sign-flip test; Holm correction across 6 non-reference comparisons"
  ),
  statistic = c(
    unname(as.numeric(friedman_result$statistic)),
    NA_real_,
    NA_real_
  ),
  degrees_of_freedom = c(unname(as.numeric(friedman_result$parameter)), NA_real_, NA_real_),
  pvalue_or_min_adjusted_pvalue = c(
    friedman_result$p.value,
    min(pairwise$Holm_adjusted_sign_flip_pvalue_all_21),
    min(versus_osars$Holm_adjusted_sign_flip_pvalue_vs_OSARS[nonreference])
  ),
  significant_comparisons_0_05 = c(
    as.integer(friedman_result$p.value < 0.05),
    sum(pairwise$significant_after_Holm_0_05),
    sum(versus_osars$significant_after_Holm_0_05[nonreference])
  ),
  note = c(
    sprintf("Kendall's W = %.3f", kendall_w),
    "No adjusted pairwise P value was below 0.05",
    "No adjusted comparison versus OSARS was below 0.05"
  ),
  stringsAsFactors = FALSE
)
write.csv(
  summary_table,
  file.path(processed_root, "Top7_Cindex_statistical_summary.csv"),
  row.names = FALSE
)

cat(sprintf("Friedman chi-squared(%d) = %.4f, P = %.6f\n",
            as.integer(friedman_result$parameter),
            as.numeric(friedman_result$statistic),
            friedman_result$p.value))
cat(sprintf("Holm-significant pairwise comparisons: %d/21\n",
            sum(pairwise$significant_after_Holm_0_05)))
cat(sprintf("Holm-significant comparisons versus OSARS: %d/6\n",
            sum(versus_osars$significant_after_Holm_0_05[nonreference])))
