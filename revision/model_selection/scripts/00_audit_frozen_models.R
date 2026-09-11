#!/usr/bin/env Rscript

args_all <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args_all[grepl("^--file=", args_all)])
script_path <- normalizePath(script_arg[[1]], mustWork = TRUE)
analysis_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
source(file.path(analysis_root, "scripts", "common.R"))

frozen_result <- Sys.getenv(
  "HCC_OS_FROZEN_RESULT",
  file.path(project_root, "data", "res_no_GSE14520.Rdata")
)
if (!file.exists(frozen_result)) stop("Frozen historical result not found: ", frozen_result)

load(frozen_result)
if (!exists("res") || is.null(res$ml.res) || is.null(res$riskscore)) stop("Malformed frozen result")

historical_names <- names(res$ml.res)
generated <- generate_configurations()
if (length(historical_names) != 117L) stop("Historical object does not contain 117 models")
if (!identical(historical_names, generated$configuration)) {
  missing_generated <- setdiff(historical_names, generated$configuration)
  extra_generated <- setdiff(generated$configuration, historical_names)
  stop(
    "Generated grid does not exactly match frozen model names. Missing: ",
    paste(missing_generated, collapse = ";"), " Extra: ", paste(extra_generated, collapse = ";")
  )
}

extract_historical_scores <- function(model_name, cohort) {
  tab <- as.data.frame(res$riskscore[[model_name]][[cohort]])
  if (!"RS" %in% names(tab)) stop("Historical score table is malformed")
  ids <- if ("ID" %in% names(tab)) as.character(tab$ID) else as.character(seq_len(nrow(tab)))
  stats::setNames(as.numeric(tab$RS), ids)
}

cohorts <- names(res$riskscore[[1]])
pairs <- list()
pair_index <- 0L
for (i in seq_len(length(historical_names) - 1L)) {
  for (j in (i + 1L):length(historical_names)) {
    exact_by_cohort <- logical(length(cohorts))
    rho_by_cohort <- numeric(length(cohorts))
    max_delta_by_cohort <- numeric(length(cohorts))
    for (k in seq_along(cohorts)) {
      a <- extract_historical_scores(historical_names[i], cohorts[k])
      b <- extract_historical_scores(historical_names[j], cohorts[k])
      common <- intersect(names(a), names(b))
      a <- a[common]
      b <- b[common]
      max_delta_by_cohort[k] <- max(abs(a - b), na.rm = TRUE)
      exact_by_cohort[k] <- isTRUE(all.equal(a, b, tolerance = 1e-12, check.attributes = FALSE))
      rho_by_cohort[k] <- suppressWarnings(stats::cor(a, b, method = "spearman", use = "complete.obs"))
    }
    exact_all <- all(exact_by_cohort)
    rank_equivalent_all <- all(is.finite(rho_by_cohort) & abs(rho_by_cohort - 1) < 1e-12)
    if (exact_all || rank_equivalent_all) {
      pair_index <- pair_index + 1L
      pairs[[pair_index]] <- data.frame(
        model_1 = historical_names[i],
        model_2 = historical_names[j],
        exact_predictions_all_historical_cohorts = exact_all,
        rank_equivalent_all_historical_cohorts = rank_equivalent_all,
        max_absolute_difference = max(max_delta_by_cohort),
        minimum_spearman_rho = min(rho_by_cohort),
        stringsAsFactors = FALSE
      )
    }
  }
}

duplicate_table <- if (length(pairs)) do.call(rbind, pairs) else data.frame(
  model_1 = character(), model_2 = character(),
  exact_predictions_all_historical_cohorts = logical(),
  rank_equivalent_all_historical_cohorts = logical(),
  max_absolute_difference = numeric(), minimum_spearman_rho = numeric()
)

historical_listing <- generated
historical_listing$present_in_frozen_object <- historical_listing$configuration %in% historical_names
historical_listing$historical_order_matches <- identical(historical_names, generated$configuration)
write.csv(historical_listing, file.path(processed_root, "all_117_model_configurations_historical_audit.csv"), row.names = FALSE)
write.csv(duplicate_table, file.path(processed_root, "historical_duplicate_predictor_audit.csv"), row.names = FALSE)

historical_c <- as.data.frame(res$Cindex.res)
write.csv(historical_c, file.path(processed_root, "historical_117_model_Cindex_snapshot.csv"), row.names = FALSE)

writeLines(c(
  paste0("Frozen object: ", normalizePath(frozen_result)),
  paste0("Frozen object SHA256: ", sha256_file(frozen_result)),
  paste0("Configurations: ", length(historical_names)),
  paste0("Single algorithms: ", sum(generated$mode == "single")),
  paste0("Combined configurations: ", sum(generated$mode == "double")),
  paste0("Exact/rank-equivalent historical pairs: ", nrow(duplicate_table)),
  "This audit reads the historical object only. Its values are not read by the TCGA CV/model-selection scripts."
), file.path(logs_root, "00_frozen_model_audit.log"))

cat("Historical audit completed: 117 configurations (20 single + 97 combined).\n")
print(duplicate_table)
