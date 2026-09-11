#!/usr/bin/env Rscript

# Multi-method bulk projection of the frozen discovery MP4 program.
#
# This analysis uses the pre-existing, frozen OSARS_frozen score and its
# cohort-median High/Low group. The state gene set is fixed before reading bulk
# outcomes: discovery MP4 genes excluding Hallmark ROS genes and genes that
# occur among the 128 frozen OSARS model features. The same state gene set is
# scored by six sample-level approaches. It does not derive a new bulk
# signature from the high/low comparison.

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(GSVA))
suppressPackageStartupMessages(library(jsonlite))

args <- commandArgs(trailingOnly = FALSE)
script_arg <- args[grepl("^--file=", args)]
if (length(script_arg) != 1L) {
  stop("Run this script with Rscript")
}
script_path <- sub("^--file=", "", script_arg)
# This Rscript wrapper encodes spaces in the --file argument as ~+~.
script_path <- gsub("~\\+~", " ", script_path)
script_path <- normalizePath(script_path)
ROOT <- normalizePath(file.path(dirname(script_path), ".."))
GSEA_ROOT <- dirname(ROOT)
REVISION <- dirname(GSEA_ROOT)
MODEL_VALIDATION <- file.path(REVISION, "model_validation")
SINGLE_CELL <- file.path(REVISION, "single_cell")
OUT <- file.path(ROOT, "data", "processed")
TABLES <- file.path(ROOT, "tables")
REFERENCES <- file.path(ROOT, "references")

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLES, recursive = TRUE, showWarnings = FALSE)
dir.create(REFERENCES, recursive = TRUE, showWarnings = FALSE)

SEED <- 20260906L
BOOTSTRAPS <- 2000L

read_expression <- function(path) {
  matrix_data <- read.csv(gzfile(path), row.names = 1, check.names = FALSE)
  matrix_data <- as.matrix(matrix_data)
  storage.mode(matrix_data) <- "double"
  if (anyDuplicated(rownames(matrix_data)) || anyDuplicated(colnames(matrix_data))) {
    stop("Duplicate gene or sample identifiers in ", path)
  }
  if (!all(is.finite(matrix_data))) {
    stop("Non-finite expression values in ", path)
  }
  matrix_data
}

bootstrap_mean_difference <- function(score, group, seed) {
  high <- score[group == "High"]
  low <- score[group == "Low"]
  set.seed(seed)
  draws <- replicate(
    BOOTSTRAPS,
    mean(sample(high, length(high), replace = TRUE)) -
      mean(sample(low, length(low), replace = TRUE))
  )
  unname(quantile(draws, c(0.025, 0.975), names = FALSE, type = 6))
}

bootstrap_standardized_difference <- function(score, group, seed) {
  high <- score[group == "High"]
  low <- score[group == "Low"]
  set.seed(seed)
  draws <- replicate(
    BOOTSTRAPS,
    {
      high_draw <- sample(high, length(high), replace = TRUE)
      low_draw <- sample(low, length(low), replace = TRUE)
      pooled_sd <- sqrt(
        ((length(high_draw) - 1) * var(high_draw) +
           (length(low_draw) - 1) * var(low_draw)) /
          (length(high_draw) + length(low_draw) - 2)
      )
      (mean(high_draw) - mean(low_draw)) / pooled_sd
    }
  )
  draws <- draws[is.finite(draws)]
  unname(quantile(draws, c(0.025, 0.975), names = FALSE, type = 6))
}

bootstrap_spearman <- function(x, y, seed) {
  set.seed(seed)
  n <- length(x)
  draws <- replicate(
    BOOTSTRAPS,
    {
      index <- sample.int(n, n, replace = TRUE)
      suppressWarnings(cor(x[index], y[index], method = "spearman"))
    }
  )
  draws <- draws[is.finite(draws)]
  unname(quantile(draws, c(0.025, 0.975), names = FALSE, type = 6))
}

fractional_rank_module <- function(expression, gene_index) {
  ranks <- apply(expression, 2, rank, ties.method = "average")
  if (is.null(dim(ranks))) {
    ranks <- matrix(ranks, ncol = 1)
  }
  colMeans(ranks[gene_index, , drop = FALSE] / nrow(expression))
}

mean_gene_z_module <- function(expression, gene_index) {
  selected <- expression[gene_index, , drop = FALSE]
  standardized <- t(scale(t(selected)))
  standardized[!is.finite(standardized)] <- 0
  colMeans(standardized)
}

orient_to_rank_module <- function(score, rank_module) {
  association <- suppressWarnings(cor(score, rank_module, method = "spearman"))
  if (!is.finite(association)) {
    stop("Could not orient a sample-level state score")
  }
  if (association < 0) {
    return(list(score = -score, flipped = TRUE, rho = -association))
  }
  list(score = score, flipped = FALSE, rho = association)
}

score_gene_set <- function(expression, gene_set, method) {
  arguments <- list(
    expr = expression,
    gset.idx.list = gene_set,
    method = method,
    min.sz = 5,
    max.sz = 500,
    verbose = FALSE,
    parallel.sz = 1
  )
  if (method %in% c("gsva", "ssgsea")) {
    arguments$kcdf <- "Gaussian"
  }
  if (method == "ssgsea") {
    arguments$ssgsea.norm <- TRUE
  }
  as.numeric(do.call(GSVA::gsva, arguments)[1, ])
}

fit_one_cohort <- function(cohort, expression_path, score_path, clinical_path,
                           state_genes) {
  expression <- read_expression(expression_path)
  scores <- read.csv(score_path, check.names = FALSE)
  clinical <- read.csv(clinical_path, check.names = FALSE)
  rownames(scores) <- as.character(scores$ID)
  rownames(clinical) <- as.character(clinical$ID)
  colnames(expression) <- as.character(colnames(expression))
  if (anyDuplicated(rownames(scores)) || anyDuplicated(rownames(clinical))) {
    stop("Duplicate metadata IDs in ", cohort)
  }
  if (!all(rownames(scores) %in% colnames(expression))) {
    stop("Frozen score IDs missing from expression in ", cohort)
  }
  if (!identical(sort(rownames(scores)), sort(rownames(clinical)))) {
    stop("Frozen score and clinical ID sets differ in ", cohort)
  }
  expression <- expression[, rownames(scores), drop = FALSE]
  clinical <- clinical[rownames(scores), , drop = FALSE]
  group <- as.character(scores$OSARS_frozen_group)
  if (!identical(sort(unique(group)), c("High", "Low"))) {
    stop("Unexpected OSARS group labels in ", cohort)
  }
  risk <- as.numeric(scores$OSARS_frozen)

  upper_expression_genes <- toupper(rownames(expression))
  if (anyDuplicated(upper_expression_genes)) {
    stop("Case-insensitive duplicate gene symbols in ", cohort)
  }
  state_indices <- match(toupper(state_genes), upper_expression_genes)
  state_indices <- state_indices[!is.na(state_indices)]
  if (length(state_indices) < 5L) {
    stop("Too few state genes available in ", cohort)
  }
  state_genes_available <- rownames(expression)[state_indices]
  gene_set <- list(MP4_without_ROS_and_OSARS_features = state_genes_available)
  rank_module <- fractional_rank_module(expression, state_indices)
  raw_scores <- list(
    rank_module = rank_module,
    mean_gene_z = mean_gene_z_module(expression, state_indices),
    GSVA = score_gene_set(expression, gene_set, "gsva"),
    ssGSEA = score_gene_set(expression, gene_set, "ssgsea"),
    PLAGE = score_gene_set(expression, gene_set, "plage"),
    GSVA_zscore = score_gene_set(expression, gene_set, "zscore")
  )
  pca_score <- as.numeric(
    prcomp(t(expression[state_indices, , drop = FALSE]), center = TRUE, scale. = TRUE,
           rank. = 1)$x[, 1]
  )
  raw_scores$PCA_PC1 <- pca_score

  oriented <- lapply(raw_scores, orient_to_rank_module, rank_module = rank_module)
  score_frame <- data.frame(
    cohort = cohort,
    ID = rownames(scores),
    OSARS_frozen = risk,
    OSARS_group = factor(group, levels = c("Low", "High")),
    stage = as.character(clinical$stage),
    check.names = FALSE
  )
  for (method in names(oriented)) {
    score_frame[[method]] <- oriented[[method]]$score
  }
  orientation <- data.frame(
    cohort = cohort,
    method = names(oriented),
    flipped_to_align_with_rank_module = vapply(oriented, function(x) x$flipped, logical(1)),
    spearman_with_rank_module_after_orientation = vapply(
      oriented, function(x) x$rho, numeric(1)
    ),
    stringsAsFactors = FALSE
  )

  method_stats <- lapply(seq_along(names(oriented)), function(index) {
    method <- names(oriented)[index]
    score <- score_frame[[method]]
    high <- score[group == "High"]
    low <- score[group == "Low"]
    wt <- suppressWarnings(wilcox.test(high, low, alternative = "two.sided",
                                       exact = FALSE, conf.int = FALSE))
    pooled_sd <- sqrt(
      ((length(high) - 1) * var(high) + (length(low) - 1) * var(low)) /
        (length(high) + length(low) - 2)
    )
    standardized_difference <- (mean(high) - mean(low)) / pooled_sd
    mean_ci <- bootstrap_mean_difference(score, group, SEED + index)
    standardized_ci <- bootstrap_standardized_difference(score, group, SEED + 50L + index)
    stage_complete <- !is.na(score_frame$stage) & nzchar(score_frame$stage)
    stage_data <- data.frame(
      score = score[stage_complete],
      group = factor(group[stage_complete], levels = c("Low", "High")),
      stage = factor(score_frame$stage[stage_complete])
    )
    stage_fit <- lm(score ~ group + stage, data = stage_data)
    stage_coef <- summary(stage_fit)$coefficients["groupHigh", , drop = FALSE]
    data.frame(
      cohort = cohort,
      method = method,
      n_total = length(score),
      n_low = length(low),
      n_high = length(high),
      mean_low = mean(low),
      mean_high = mean(high),
      mean_difference_high_minus_low = mean(high) - mean(low),
      mean_difference_ci_lower = mean_ci[1],
      mean_difference_ci_upper = mean_ci[2],
      standardized_mean_difference_high_minus_low = standardized_difference,
      standardized_mean_difference_ci_lower = standardized_ci[1],
      standardized_mean_difference_ci_upper = standardized_ci[2],
      wilcoxon_p = wt$p.value,
      n_stage_complete = nrow(stage_data),
      stage_adjusted_coefficient_high_minus_low = stage_coef[1, "Estimate"],
      stage_adjusted_standard_error = stage_coef[1, "Std. Error"],
      stage_adjusted_t = stage_coef[1, "t value"],
      stage_adjusted_p = stage_coef[1, "Pr(>|t|)"],
      stringsAsFactors = FALSE
    )
  })
  method_stats <- do.call(rbind, method_stats)

  correlation_stats <- lapply(seq_along(names(oriented)), function(index) {
    method <- names(oriented)[index]
    score <- score_frame[[method]]
    result <- suppressWarnings(cor.test(risk, score, method = "spearman", exact = FALSE))
    ci <- bootstrap_spearman(risk, score, SEED + 100L + index)
    data.frame(
      cohort = cohort,
      method = method,
      n = length(score),
      spearman_rho_OSARS_vs_state_score = unname(result$estimate),
      spearman_p = result$p.value,
      spearman_ci_lower = ci[1],
      spearman_ci_upper = ci[2],
      stringsAsFactors = FALSE
    )
  })
  correlation_stats <- do.call(rbind, correlation_stats)
  trend_coordinates <- do.call(rbind, lapply(names(oriented), function(method) {
    score <- score_frame[[method]]
    fit <- lm(score ~ risk)
    risk_grid <- seq(min(risk), max(risk), length.out = 100L)
    data.frame(
      cohort = cohort,
      method = method,
      OSARS_frozen = risk_grid,
      fitted_state_score = unname(predict(fit, newdata = data.frame(risk = risk_grid))),
      stringsAsFactors = FALSE
    )
  }))
  coverage <- data.frame(
    cohort = cohort,
    gene = state_genes,
    available = toupper(state_genes) %in% upper_expression_genes,
    stringsAsFactors = FALSE
  )
  list(
    scores = score_frame,
    method_stats = method_stats,
    correlation_stats = correlation_stats,
    trend_coordinates = trend_coordinates,
    coverage = coverage,
    orientation = orientation
  )
}

nmf <- read.csv(file.path(SINGLE_CELL, "data", "processed", "discovery_NMF_genes.csv"))
mp4 <- unique(as.character(nmf$gene[nmf$program == "MP4"]))
ros <- scan(
  file.path(MODEL_VALIDATION, "references", "HALLMARK_ROS_msigdb_v7.0_genes.txt"),
  what = "character", quiet = TRUE
)
osars_features <- read.csv(
  file.path(REVISION, "tables", "Table2_frozen_OSARS_128_genes.csv")
)$gene
strict_state_genes <- sort(setdiff(setdiff(toupper(mp4), toupper(ros)), toupper(osars_features)))

write.csv(
  data.frame(gene = strict_state_genes, stringsAsFactors = FALSE),
  file.path(REFERENCES, "frozen_MP4_without_ROS_and_OSARS_features.csv"),
  row.names = FALSE
)

cohort_specs <- list(
  "TCGA-LIHC" = list(
    expression = file.path(MODEL_VALIDATION, "data", "raw", "TCGA_full_expression_logscale.csv.gz"),
    scores = file.path(MODEL_VALIDATION, "data", "processed", "TCGA-LIHC_all_scores.csv"),
    clinical = file.path(MODEL_VALIDATION, "data", "processed", "TCGA-LIHC_clinical_analysis_rows.csv")
  ),
  "ICGC-LIRI" = list(
    expression = file.path(MODEL_VALIDATION, "data", "raw", "ICGC_full_expression_logscale.csv.gz"),
    scores = file.path(MODEL_VALIDATION, "data", "processed", "ICGC-LIRI_all_scores.csv"),
    clinical = file.path(MODEL_VALIDATION, "data", "processed", "ICGC-LIRI_clinical_analysis_rows.csv")
  )
)

results <- lapply(names(cohort_specs), function(cohort) {
  spec <- cohort_specs[[cohort]]
  fit_one_cohort(cohort, spec$expression, spec$scores, spec$clinical, strict_state_genes)
})
names(results) <- names(cohort_specs)

scores <- do.call(rbind, lapply(results, function(x) x$scores))
method_stats <- do.call(rbind, lapply(results, function(x) x$method_stats))
correlation_stats <- do.call(rbind, lapply(results, function(x) x$correlation_stats))
trend_coordinates <- do.call(rbind, lapply(results, function(x) x$trend_coordinates))
coverage <- do.call(rbind, lapply(results, function(x) x$coverage))
orientation <- do.call(rbind, lapply(results, function(x) x$orientation))

method_stats$wilcoxon_BH_q_within_cohort <- ave(
  method_stats$wilcoxon_p, method_stats$cohort, FUN = p.adjust, method = "BH"
)
method_stats$stage_adjusted_BH_q_within_cohort <- ave(
  method_stats$stage_adjusted_p, method_stats$cohort, FUN = p.adjust, method = "BH"
)
correlation_stats$spearman_BH_q_within_cohort <- ave(
  correlation_stats$spearman_p, correlation_stats$cohort, FUN = p.adjust, method = "BH"
)

method_definitions <- data.frame(
  method = c("rank_module", "mean_gene_z", "GSVA", "ssGSEA", "PLAGE", "GSVA_zscore", "PCA_PC1"),
  definition = c(
    "Mean within-sample fractional expression rank of the fixed state genes",
    "Mean gene-wise z score across samples for the fixed state genes",
    "GSVA score from the GSVA R package",
    "Single-sample GSEA score from the GSVA R package",
    "PLAGE latent-factor score from the GSVA R package; sign oriented to the rank-module score",
    "Gene-set z-score method from the GSVA R package",
    "First principal-component score across the fixed state genes; sign oriented to the rank-module score"
  ),
  stringsAsFactors = FALSE
)

write.csv(scores, file.path(OUT, "per_sample_state_projection_scores.csv"), row.names = FALSE)
write.csv(method_stats, file.path(TABLES, "state_projection_method_statistics.csv"), row.names = FALSE)
write.csv(correlation_stats, file.path(TABLES, "state_projection_continuous_correlations.csv"), row.names = FALSE)
write.csv(trend_coordinates, file.path(OUT, "state_projection_linear_trend_coordinates.csv"), row.names = FALSE)
write.csv(coverage, file.path(OUT, "state_gene_coverage.csv"), row.names = FALSE)
write.csv(orientation, file.path(OUT, "score_orientation_audit.csv"), row.names = FALSE)
write.csv(method_definitions, file.path(REFERENCES, "method_definitions.csv"), row.names = FALSE)

manifest <- list(
  analysis = "Multi-method sample-level projection of frozen discovery MP4 state genes",
  created = "2026-09-06",
  frozen_OSARS_groups = "OSARS_frozen_group from the existing cohort-median risk split; no model refit or cutoff optimization",
  state_gene_set = list(
    source_MP4_genes = 58,
    excluded_Hallmark_ROS_genes = c("ATOX1", "TXN"),
    excluded_OSARS_feature_overlap = 6,
    strict_state_gene_count = length(strict_state_genes)
  ),
  methods = method_definitions,
  inferential_tests = list(
    group_comparison = "Two-sided Wilcoxon rank-sum; BH across seven methods within each cohort",
    stage_adjusted = "Ordinary least squares score ~ OSARS group + categorical stage; BH across seven methods within each cohort",
    continuous_association = "Two-sided Spearman correlation; percentile bootstrap 95% CI with 2,000 sample-row resamples"
  ),
  score_orientation = "PLAGE and PCA signs were aligned to the fixed within-sample rank module score because those factor signs are otherwise arbitrary.",
  no_deconvolution = "This is transcriptional program projection, not cell-fraction deconvolution."
)
write_json(manifest, file.path(OUT, "analysis_manifest.json"), pretty = TRUE, auto_unbox = TRUE)

cat("Completed multi-method state projection.\n")
print(method_stats)
print(correlation_stats)
