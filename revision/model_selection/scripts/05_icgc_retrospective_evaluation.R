#!/usr/bin/env Rscript

args_all <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args_all[grepl("^--file=", args_all)])
script_path <- normalizePath(script_arg[[1]], mustWork = TRUE)
analysis_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
source(file.path(analysis_root, "scripts", "common.R"))

locked_file <- file.path(processed_root, "final_OSARS_model_locked.rds")
lock_file <- file.path(processed_root, "MODEL_LOCKED.sha256")
if (!file.exists(locked_file) || !file.exists(lock_file)) {
  stop("The TCGA-selected model must be fitted and cryptographically locked before external evaluation")
}
lock_parts <- strsplit(readLines(lock_file, n = 1L), "[[:space:]]+")[[1]]
expected_hash <- lock_parts[1]
observed_hash <- sha256_file(locked_file)
if (!identical(expected_hash, observed_hash)) stop("Locked model SHA256 verification failed")

# Only after the lock is verified is the retrospective external cohort opened.
bundle <- readRDS(locked_file)
icgc_file <- file.path(raw_root, "ICGC_original_448_expression.csv")
if (!file.exists(icgc_file)) stop("Frozen retrospective external-evaluation input is missing")
icgc <- read.csv(icgc_file, check.names = FALSE)
required <- c("ID", "OS.time", "OS", bundle$frozen_candidate_genes)
missing <- setdiff(required, names(icgc))
if (length(missing)) stop("ICGC frozen input lacks required fields: ", paste(missing, collapse = ", "))
if (anyDuplicated(icgc$ID)) stop("Duplicate ICGC patient IDs")
if (any(!is.finite(icgc$OS.time)) || any(icgc$OS.time <= 0) || !all(icgc$OS %in% c(0, 1))) {
  stop("Invalid frozen ICGC survival outcome")
}

risk <- predict_locked_bundle(bundle, icgc[, bundle$frozen_candidate_genes, drop = FALSE])
if (any(!is.finite(risk))) stop("Non-finite locked-model predictions in ICGC")

concordance_fit <- survival::concordance(
  survival::Surv(icgc$OS.time, icgc$OS) ~ risk,
  reverse = TRUE,
  timefix = TRUE
)
cindex <- as.numeric(concordance_fit$concordance)
se <- sqrt(as.numeric(concordance_fit$var))
ci_low <- max(0, cindex - stats::qnorm(0.975) * se)
ci_high <- min(1, cindex + stats::qnorm(0.975) * se)

patient_scores <- data.frame(
  patient_id = as.character(icgc$ID),
  OS.time = as.numeric(icgc$OS.time),
  OS = as.integer(icgc$OS),
  locked_OSARS_risk_score = risk,
  evaluation_label = "retrospective external evaluation",
  stringsAsFactors = FALSE
)
write.csv(patient_scores, file.path(processed_root, "ICGC_locked_model_patient_scores.csv"), row.names = FALSE)

evaluation <- data.frame(
  cohort = "ICGC-LIRI-JP",
  evaluation_type = "retrospective external evaluation",
  locked_configuration = bundle$configuration,
  n_patients = nrow(icgc),
  n_events = sum(icgc$OS == 1),
  Cindex = cindex,
  SE = se,
  CI95_low = ci_low,
  CI95_high = ci_high,
  CI_method = "normal approximation using concordance influence-function variance",
  model_sha256 = observed_hash,
  feature_selection_performed_in_ICGC = FALSE,
  algorithm_selection_performed_in_ICGC = FALSE,
  hyperparameter_tuning_performed_in_ICGC = FALSE,
  cutoff_optimization_performed_in_ICGC = FALSE,
  stringsAsFactors = FALSE
)
write.csv(evaluation, file.path(analysis_root, "ICGC_retrospective_evaluation.csv"), row.names = FALSE)

summary_file <- file.path(analysis_root, "final_model_summary.txt")
summary_lines <- if (file.exists(summary_file)) readLines(summary_file, warn = FALSE) else character()
summary_lines <- summary_lines[!grepl("^ICGC retrospective external evaluation", summary_lines)]
summary_lines <- c(
  summary_lines,
  paste0(
    "ICGC retrospective external evaluation: n=", nrow(icgc),
    ", events=", sum(icgc$OS == 1),
    ", C-index=", sprintf("%.6f", cindex),
    " (95% CI ", sprintf("%.6f", ci_low), " to ", sprintf("%.6f", ci_high), ")."
  )
)
writeLines(summary_lines, summary_file)

writeLines(c(
  paste0("Model SHA256 verified before external data access: ", observed_hash),
  paste0("Frozen ICGC input: ", normalizePath(icgc_file)),
  paste0("Frozen ICGC input SHA256: ", sha256_file(icgc_file)),
  paste0("Patients: ", nrow(icgc)),
  paste0("Events: ", sum(icgc$OS == 1)),
  paste0("C-index: ", sprintf("%.6f", cindex)),
  paste0("95% CI: ", sprintf("%.6f to %.6f", ci_low, ci_high)),
  "No ICGC feature selection, algorithm selection, tuning, cutoff optimization, or model modification was performed."
), file.path(logs_root, "05_ICGC_retrospective_evaluation.log"))

print(evaluation)
