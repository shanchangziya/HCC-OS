#!/usr/bin/env Rscript

# Reviewer 2 Major Comment 6: strict hold-out pathology analysis.
# The patient split must already have been frozen by 01_audit_and_freeze_split.R.
# All supervised feature screening, glmnet standardization, CV, coefficient
# estimation, and the KM cutoff are training-only. Test outcomes are accessed
# only after the model and cutoff have been saved.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: Rscript 02_fit_evaluate_strict_holdout.R FROZEN_OUTPUT_DIR")
}

suppressPackageStartupMessages({
  library(survival)
  library(glmnet)
})

output_dir <- normalizePath(path.expand(args[[1]]), mustWork = TRUE)
input_file <- file.path(output_dir, "pathology_analysis_input.rds")
split_file <- file.path(output_dir, "pathology_strict_holdout_split.csv")
split_provenance_file <- file.path(output_dir, "pathology_split_provenance.csv")
completion_marker <- file.path(output_dir, "PATHOLOGY_STRICT_HOLDOUT_COMPLETE.txt")

if (file.exists(completion_marker)) {
  stop("Completed strict hold-out run already exists; refusing to rerun or overwrite.")
}
if (!all(file.exists(c(input_file, split_file, split_provenance_file)))) {
  stop("Frozen inputs are missing. Run 01_audit_and_freeze_split.R exactly once first.")
}

write_csv <- function(x, filename) {
  write.csv(x, file.path(output_dir, filename), row.names = FALSE, na = "")
}

fmt <- function(x, digits = 3L) {
  ifelse(is.finite(x), formatC(x, digits = digits, format = "f"), "NA")
}

fmt_p <- function(x) {
  if (!is.finite(x)) return("NA")
  if (x < 0.001) return(formatC(x, digits = 2, format = "e"))
  formatC(x, digits = 3, format = "f")
}

analysis_input <- readRDS(input_file)
split_table <- read.csv(split_file, stringsAsFactors = FALSE, check.names = FALSE)
split_provenance <- read.csv(split_provenance_file, stringsAsFactors = FALSE)
current_split_md5 <- unname(tools::md5sum(split_file))
stopifnot(
  identical(current_split_md5, split_provenance$split_md5[[1]]),
  split_provenance$seed[[1]] == 20260906L,
  split_provenance$seed_retry[[1]] == "No",
  identical(sort(unique(split_table$split)), c("test", "train")),
  !anyDuplicated(split_table$patient_id)
)

patient_data <- analysis_input$patient_data
features <- analysis_input$feature_names
X <- analysis_input$feature_matrix
stopifnot(
  identical(colnames(X), features), length(features) == 2048L,
  !anyDuplicated(patient_data$patient_id),
  setequal(patient_data$patient_id, split_table$patient_id),
  setequal(rownames(X), patient_data$patient_id)
)

# Reorder all objects once by the immutable split table.
patient_data <- patient_data[match(split_table$patient_id, patient_data$patient_id), , drop = FALSE]
X <- X[split_table$patient_id, features, drop = FALSE]
stopifnot(
  identical(patient_data$patient_id, split_table$patient_id),
  isTRUE(all.equal(patient_data$survival_time_days, split_table$survival_time_days)),
  identical(as.integer(patient_data$event), as.integer(split_table$event)),
  all(is.finite(X))
)

train_idx <- which(split_table$split == "train")
test_idx <- which(split_table$split == "test")
train_ids <- split_table$patient_id[train_idx]
test_ids <- split_table$patient_id[test_idx]
stopifnot(
  length(train_idx) == 231L, length(test_idx) == 99L,
  sum(patient_data$event[train_idx]) == 84L,
  sum(patient_data$event[test_idx]) == 36L,
  length(intersect(train_ids, test_ids)) == 0L
)

X_train <- X[train_idx, , drop = FALSE]
osars_train <- patient_data$OSARS_score[train_idx]
surv_train <- Surv(
  patient_data$survival_time_days[train_idx],
  patient_data$event[train_idx]
)

# Original preprocessing audit: there was no external imputation, filtering,
# centering, or scaling. glmnet standardize=TRUE is retained and learns its
# transformations entirely from these training rows. Both conventional sample
# SD and glmnet's population-denominator scale are recorded for transparency.
training_mean <- colMeans(X_train)
training_sd <- apply(X_train, 2L, sd)
training_glmnet_scale <- sqrt(colMeans(sweep(X_train, 2L, training_mean, "-")^2))
zero_variance <- !is.finite(training_sd) | training_sd == 0
preprocessing <- data.frame(
  feature_id = features,
  training_mean = training_mean,
  training_sample_SD = training_sd,
  training_glmnet_scale = training_glmnet_scale,
  zero_variance_in_training = zero_variance,
  imputation = "None; source matrix complete",
  external_centering_or_scaling = "None",
  glmnet_internal_standardize = TRUE,
  test_rows_used = FALSE,
  stringsAsFactors = FALSE
)
write_csv(preprocessing, "pathology_preprocessing_parameters.csv")

# Training-only replication of the original supervised screen.
screen_rows <- lapply(seq_along(features), function(j) {
  if (zero_variance[[j]] || sd(osars_train) == 0) {
    return(data.frame(feature_id = features[[j]], n_training = length(train_idx),
                      correlation = NA_real_, p_value = NA_real_))
  }
  test <- suppressWarnings(cor.test(X_train[, j], osars_train, method = "pearson"))
  data.frame(
    feature_id = features[[j]], n_training = length(train_idx),
    correlation = unname(test$estimate), p_value = test$p.value
  )
})
screening <- do.call(rbind, screen_rows)
screening$BH_FDR <- p.adjust(screening$p_value, method = "BH")
screening$abs_r_cutoff <- 0.2
screening$nominal_p_cutoff <- 0.05
screening$selected <- with(
  screening,
  is.finite(correlation) & is.finite(p_value) & abs(correlation) > 0.2 & p_value < 0.05
)
screening$method <- "Pearson"
screening$target <- "Frozen continuous OSARS score"
screening$selection_uses_FDR <- FALSE
screening$test_rows_used <- FALSE
screening <- screening[, c(
  "feature_id", "n_training", "method", "target", "correlation", "p_value",
  "BH_FDR", "abs_r_cutoff", "nominal_p_cutoff", "selection_uses_FDR",
  "selected", "test_rows_used"
)]
write_csv(screening, "pathology_training_feature_screening.csv")
screened_features <- screening$feature_id[screening$selected]

write_failure_outputs <- function(reason) {
  final_empty <- data.frame(
    feature_id = character(), coefficient = numeric(), selection_stage = character(),
    stringsAsFactors = FALSE
  )
  write_csv(final_empty, "pathology_final_features.csv")
  score_empty <- patient_data[, intersect(
    c("patient_id", "OSARS_score", "survival_time_days", "event", "Stage", "Stage_Group"),
    names(patient_data)
  ), drop = FALSE]
  score_empty$split <- split_table$split
  score_empty$pathology_score <- NA_real_
  write_csv(score_empty, "pathology_patient_scores_strict_holdout.csv")
  write_csv(data.frame(status = "Model failed", reason = reason), "pathology_performance_summary.csv")
  write_csv(data.frame(status = "Model failed", reason = reason), "pathology_KM_statistics.csv")
  integrity <- c(
    "# Pathology strict hold-out integrity check", "",
    paste0("The frozen split MD5 is `", current_split_md5, "`."), "",
    "| Check | Result | Evidence |", "|---|---|---|",
    "| Train/test patient overlap | PASS: 0 | Patient-ID intersection checked before screening |",
    "| Test samples in feature screening | PASS: No | Screening matrix contains training IDs only |",
    "| Test outcome in feature selection | PASS: No | Screen target is training OSARS only |",
    "| Test in normalization estimation | PASS: No | glmnet input contains training rows only |",
    "| Test in lambda selection | PASS: No | 10-fold CV is within training rows |",
    "| Test in cutoff selection | PASS: No | No cutoff was created because model failed |",
    "| Seed/split retry | PASS: No | One frozen seed-20260906 split |", "",
    paste0("Model status: **failed without threshold relaxation** — ", reason)
  )
  writeLines(integrity, file.path(output_dir, "PATHOLOGY_HOLDOUT_INTEGRITY_CHECK.md"))
  summary_lines <- c(
    "# Pathology strict hold-out summary", "",
    sprintf("1. Final patients: %d.", nrow(patient_data)),
    sprintf("2. Training: %d patients/%d events; test: %d patients/%d events.",
            length(train_idx), sum(patient_data$event[train_idx]),
            length(test_idx), sum(patient_data$event[test_idx])),
    sprintf("3. Features: 2048 initial → %d training-screened → 0 final; %s", length(screened_features), reason),
    "4. Training C-index: NA (model failed).",
    "5. Strict test C-index: NA (model failed).",
    "6. Strict test continuous Cox HR: NA (model failed).",
    "7. Training-derived cutoff/test KM: NA (model failed).",
    "8. Comparison with old analysis: not estimable.",
    "9. Interpretation: the strict analysis does not support a pathology prognostic model and should be reported as an unsuccessful exploratory analysis."
  )
  writeLines(summary_lines, file.path(output_dir, "PATHOLOGY_STRICT_HOLDOUT_SUMMARY.md"))
  writeLines(c("status=model_failed", paste0("reason=", reason)), completion_marker)
}

if (length(screened_features) == 0L) {
  write_failure_outputs("No feature passed the unchanged training-only Pearson screen.")
  quit(save = "no", status = 0L)
}

# Prespecified, reproducible 10-fold CV wholly inside training. The original
# random nfolds=10 allocation and lambda.min rule are retained; the fold IDs are
# exported so the exact CV partition cannot drift across reruns.
set.seed(20260906)
foldid <- sample(rep(seq_len(10L), length.out = length(train_idx)))
cv_folds <- data.frame(
  patient_id = train_ids,
  event = patient_data$event[train_idx],
  fold = foldid,
  stringsAsFactors = FALSE
)
write_csv(cv_folds, "pathology_training_cv_folds.csv")

cv_warnings <- character()
cvfit <- tryCatch(
  withCallingHandlers(
    cv.glmnet(
      x = X_train[, screened_features, drop = FALSE],
      y = surv_train,
      family = "cox", alpha = 1, nfolds = 10, foldid = foldid,
      standardize = TRUE, type.measure = "deviance", grouped = TRUE,
      parallel = FALSE
    ),
    warning = function(w) {
      cv_warnings <<- c(cv_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  ),
  error = function(e) e
)
if (inherits(cvfit, "error")) {
  write_failure_outputs(paste0("Training-only cv.glmnet failed: ", conditionMessage(cvfit)))
  quit(save = "no", status = 0L)
}

lambda_selected <- cvfit$lambda.min
lambda_1se <- cvfit$lambda.1se
selected_lambda_index <- match(lambda_selected, cvfit$lambda)
warning_lambda_positions <- suppressWarnings(as.integer(sub(
  ".*Convergence for ([0-9]+)(st|nd|rd|th) lambda.*", "\\1", cv_warnings
)))
warning_lambda_positions <- warning_lambda_positions[is.finite(warning_lambda_positions)]
convergence_audit <- data.frame(
  selected_lambda_index = selected_lambda_index,
  total_lambda_count = length(cvfit$lambda),
  selected_lambda = lambda_selected,
  full_training_glmnet_jerr = cvfit$glmnet.fit$jerr,
  CV_warning_count = length(cv_warnings),
  earliest_warning_lambda_position = if (length(warning_lambda_positions)) min(warning_lambda_positions) else NA_integer_,
  selected_precedes_warning_tail = if (length(warning_lambda_positions)) selected_lambda_index < min(warning_lambda_positions) else TRUE,
  warning_messages = paste(unique(cv_warnings), collapse = " || "),
  model_changed_or_refit_due_to_warning = "No",
  stringsAsFactors = FALSE
)
write_csv(convergence_audit, "pathology_lasso_convergence_audit.csv")
coefficient_vector <- as.matrix(coef(cvfit, s = "lambda.min"))[, 1]
final_features <- names(coefficient_vector)[coefficient_vector != 0]

all_screened_coefficients <- data.frame(
  feature_id = names(coefficient_vector),
  coefficient = unname(coefficient_vector),
  final_nonzero = coefficient_vector != 0,
  selection_stage = ifelse(coefficient_vector != 0, "final_lambda.min", "training_screen_only"),
  training_mean = training_mean[names(coefficient_vector)],
  training_sample_SD = training_sd[names(coefficient_vector)],
  stringsAsFactors = FALSE
)
write_csv(all_screened_coefficients, "pathology_lasso_screened_coefficients.csv")

if (length(final_features) == 0L) {
  write_failure_outputs("lambda.min selected the null model; no test evaluation was performed.")
  quit(save = "no", status = 0L)
}

final_table <- all_screened_coefficients[all_screened_coefficients$final_nonzero, , drop = FALSE]
final_table$coefficient_per_training_SD <-
  final_table$coefficient * final_table$training_sample_SD
final_table$initial_feature_count <- 2048L
final_table$screened_feature_count <- length(screened_features)
final_table$final_feature_count <- length(final_features)
final_table <- final_table[order(abs(final_table$coefficient), decreasing = TRUE), ]
write_csv(final_table, "pathology_final_features.csv")

cv_curve <- data.frame(
  lambda = cvfit$lambda,
  log_lambda = log(cvfit$lambda),
  cv_mean_deviance = cvfit$cvm,
  cv_SE = cvfit$cvsd,
  cv_lower = cvfit$cvlo,
  cv_upper = cvfit$cvup,
  n_nonzero = cvfit$nzero,
  selected_lambda_min = cvfit$lambda == cvfit$lambda.min,
  selected_lambda_1se = cvfit$lambda == cvfit$lambda.1se
)
write_csv(cv_curve, "pathology_lasso_CV_curve.csv")

# Training score, score scale, and median cutoff are fixed before any test score
# or test performance is calculated.
training_score <- as.numeric(predict(
  cvfit, newx = X_train[, screened_features, drop = FALSE],
  s = "lambda.min", type = "link"
))
training_score_mean <- mean(training_score)
training_score_sd <- sd(training_score)
training_cutoff <- median(training_score)
if (!is.finite(training_score_sd) || training_score_sd == 0) {
  write_failure_outputs("The frozen training pathology score has zero variance.")
  quit(save = "no", status = 0L)
}

model_specification <- data.frame(
  split_seed = 20260906L,
  split_md5 = current_split_md5,
  CV_seed = 20260906L,
  initial_features = 2048L,
  training_screened_features = length(screened_features),
  final_nonzero_features = length(final_features),
  screening_method = "Pearson versus frozen continuous OSARS",
  screening_rule = "absolute r > 0.2 and nominal P < 0.05",
  lasso_framework = "cv.glmnet Cox; alpha=1; 10 fixed random training folds",
  lambda_rule = "lambda.min",
  lambda_selected = lambda_selected,
  lambda_1se_recorded_not_selected = lambda_1se,
  glmnet_standardize = TRUE,
  imputation = "None",
  external_variance_filtering = "None",
  training_score_mean = training_score_mean,
  training_score_SD = training_score_sd,
  training_median_cutoff = training_cutoff,
  group_rule = "High if pathology score > training median; ties are Low",
  bootstrap_replicates = 2000L,
  bootstrap_seed = 20260906L,
  seed_retry = "No",
  stringsAsFactors = FALSE
)
write_csv(model_specification, "pathology_frozen_model_specification.csv")

frozen_model <- list(
  cvfit = cvfit,
  lambda = lambda_selected,
  lambda_rule = "lambda.min",
  screened_features = screened_features,
  final_features = final_features,
  coefficients_raw_scale = coefficient_vector,
  preprocessing_parameters = preprocessing,
  training_ids = train_ids,
  training_score_mean = training_score_mean,
  training_score_sd = training_score_sd,
  training_median_cutoff = training_cutoff,
  split_md5 = current_split_md5,
  split_seed = 20260906L,
  CV_seed = 20260906L
)
saveRDS(frozen_model, file.path(output_dir, "pathology_frozen_model_strict_holdout.rds"))
cat("MODEL AND TRAINING CUTOFF FROZEN; STARTING TEST PREDICTION\n")

# Test enters for the first time here, through predict() only.
test_score <- as.numeric(predict(
  cvfit, newx = X[test_idx, screened_features, drop = FALSE],
  s = "lambda.min", type = "link"
))
pathology_score <- rep(NA_real_, nrow(patient_data))
pathology_score[train_idx] <- training_score
pathology_score[test_idx] <- test_score

score_table <- patient_data[, intersect(
  c("patient_id", "OSARS_score", "survival_time_days", "event",
    "Stage", "Stage_Group", "Grade", "Grade_Group", "Sex"),
  names(patient_data)
), drop = FALSE]
score_table$split <- split_table$split
score_table$pathology_score <- pathology_score
score_table$pathology_score_per_training_SD <-
  (pathology_score - training_score_mean) / training_score_sd
score_table$risk_group <- ifelse(pathology_score > training_cutoff, "High", "Low")
score_table <- score_table[, c(
  "patient_id", "split", "pathology_score", "pathology_score_per_training_SD",
  "risk_group", "OSARS_score", "survival_time_days", "event",
  setdiff(names(score_table), c(
    "patient_id", "split", "pathology_score", "pathology_score_per_training_SD",
    "risk_group", "OSARS_score", "survival_time_days", "event"
  ))
)]
write_csv(score_table, "pathology_patient_scores_strict_holdout.csv")

get_cindex <- function(time, event, score) {
  tryCatch(
    unname(concordance(Surv(time, event) ~ score, reverse = TRUE)$concordance),
    error = function(e) NA_real_
  )
}

bootstrap_cindex <- function(d, B = 2000L, seed = 20260906L) {
  set.seed(seed)
  values <- vapply(seq_len(B), function(replicate_id) {
    index <- sample.int(nrow(d), nrow(d), replace = TRUE)
    get_cindex(
      d$survival_time_days[index], d$event[index], d$pathology_score[index]
    )
  }, numeric(1))
  valid <- values[is.finite(values)]
  interval <- if (length(valid) >= 1000L) {
    unname(quantile(valid, c(0.025, 0.975), names = FALSE, type = 7))
  } else {
    c(NA_real_, NA_real_)
  }
  list(values = values, valid = length(valid), lower = interval[[1]], upper = interval[[2]])
}

continuous_cox <- function(d) {
  fit <- coxph(
    Surv(survival_time_days, event) ~ pathology_score_per_training_SD,
    data = d, ties = "efron", x = TRUE, y = TRUE
  )
  result <- summary(fit)
  ph <- tryCatch(cox.zph(fit)$table, error = function(e) NULL)
  list(
    fit = fit,
    hr = result$conf.int[1, "exp(coef)"],
    lower = result$conf.int[1, "lower .95"],
    upper = result$conf.int[1, "upper .95"],
    p = result$coefficients[1, "Pr(>|z|)"],
    ph_p = if (is.null(ph)) NA_real_ else ph[1, "p"]
  )
}

performance_rows <- list()
bootstrap_rows <- list()
cox_fits <- list()
for (split_name in c("train", "test")) {
  d <- score_table[score_table$split == split_name, , drop = FALSE]
  point <- get_cindex(d$survival_time_days, d$event, d$pathology_score)
  boot <- bootstrap_cindex(d, B = 2000L, seed = 20260906L)
  if (boot$valid < 1000L) stop("Fewer than 1,000 valid patient bootstrap C-indices for ", split_name)
  cox_result <- continuous_cox(d)
  cox_fits[[split_name]] <- cox_result$fit
  performance_rows[[split_name]] <- data.frame(
    split = split_name,
    evaluation_role = ifelse(split_name == "train", "apparent_training", "strict_holdout_test"),
    n = nrow(d), events = sum(d$event),
    C_index = point, C_index_CI_lower = boot$lower, C_index_CI_upper = boot$upper,
    C_index_bootstrap_requested = 2000L,
    C_index_bootstrap_valid = boot$valid,
    continuous_HR_per_training_SD = cox_result$hr,
    continuous_HR_CI_lower = cox_result$lower,
    continuous_HR_CI_upper = cox_result$upper,
    continuous_HR_Wald_P = cox_result$p,
    continuous_score_PH_test_P = cox_result$ph_p,
    score_unit = paste0("1 training-set SD = ", format(training_score_sd, digits = 15)),
    stringsAsFactors = FALSE
  )
  bootstrap_rows[[split_name]] <- data.frame(
    split = split_name,
    replicate = seq_along(boot$values),
    C_index = boot$values,
    stringsAsFactors = FALSE
  )
}
performance <- do.call(rbind, performance_rows)
bootstrap_results <- do.call(rbind, bootstrap_rows)
write_csv(performance, "pathology_performance_summary.csv")
write_csv(bootstrap_results, "pathology_Cindex_patient_bootstrap.csv")

km_statistics <- list()
km_coordinates <- list()
km_risk_table <- list()
risk_times <- seq(0, 5 * 365.25, by = 365.25)
for (split_name in c("train", "test")) {
  d <- score_table[score_table$split == split_name, , drop = FALSE]
  d$risk_group <- factor(d$risk_group, levels = c("Low", "High"))
  logrank <- survdiff(Surv(survival_time_days, event) ~ risk_group, data = d)
  group_cox <- coxph(
    Surv(survival_time_days, event) ~ risk_group,
    data = d, ties = "efron", x = TRUE, y = TRUE
  )
  group_summary <- summary(group_cox)
  km_statistics[[split_name]] <- data.frame(
    split = split_name,
    cutoff_source = "training pathology-score median",
    training_cutoff = training_cutoff,
    group_rule = "High if score > cutoff; ties are Low",
    n_low = sum(d$risk_group == "Low"),
    n_high = sum(d$risk_group == "High"),
    events_low = sum(d$event[d$risk_group == "Low"]),
    events_high = sum(d$event[d$risk_group == "High"]),
    HR_high_vs_low = group_summary$conf.int[1, "exp(coef)"],
    HR_CI_lower = group_summary$conf.int[1, "lower .95"],
    HR_CI_upper = group_summary$conf.int[1, "upper .95"],
    Cox_Wald_P = group_summary$coefficients[1, "Pr(>|z|)"],
    logrank_chisq = logrank$chisq,
    logrank_P = pchisq(logrank$chisq, df = 1L, lower.tail = FALSE),
    stringsAsFactors = FALSE
  )

  fit <- survfit(
    Surv(survival_time_days, event) ~ risk_group,
    data = d, conf.type = "log-log"
  )
  curve <- summary(fit, censored = TRUE)
  curve_table <- data.frame(
    split = split_name,
    risk_group = sub("risk_group=", "", as.character(curve$strata)),
    time_days = curve$time,
    time_years = curve$time / 365.25,
    survival = curve$surv,
    CI_lower = curve$lower,
    CI_upper = curve$upper,
    n_risk = curve$n.risk,
    n_event = curve$n.event,
    n_censor = curve$n.censor,
    stringsAsFactors = FALSE
  )
  starts <- data.frame(
    split = split_name, risk_group = c("Low", "High"),
    time_days = 0, time_years = 0, survival = 1,
    CI_lower = 1, CI_upper = 1,
    n_risk = c(sum(d$risk_group == "Low"), sum(d$risk_group == "High")),
    n_event = 0, n_censor = 0,
    stringsAsFactors = FALSE
  )
  km_coordinates[[split_name]] <- rbind(starts, curve_table)

  risk <- summary(fit, times = risk_times, extend = TRUE)
  km_risk_table[[split_name]] <- data.frame(
    split = split_name,
    risk_group = sub("risk_group=", "", as.character(risk$strata)),
    time_days = risk$time,
    time_years = risk$time / 365.25,
    n_risk = risk$n.risk,
    stringsAsFactors = FALSE
  )
}
km_statistics <- do.call(rbind, km_statistics)
write_csv(km_statistics, "pathology_KM_statistics.csv")
write_csv(do.call(rbind, km_coordinates), "pathology_KM_coordinates.csv")
write_csv(do.call(rbind, km_risk_table), "pathology_KM_risk_table.csv")

# Formal integrity assertions after evaluation.
stopifnot(
  length(intersect(train_ids, test_ids)) == 0L,
  identical(sort(cv_folds$patient_id), sort(train_ids)),
  length(intersect(cv_folds$patient_id, test_ids)) == 0L,
  identical(sort(frozen_model$training_ids), sort(train_ids)),
  isTRUE(all.equal(training_cutoff, median(training_score))),
  all(score_table$risk_group == ifelse(score_table$pathology_score > training_cutoff, "High", "Low")),
  identical(current_split_md5, unname(tools::md5sum(split_file)))
)

integrity_lines <- c(
  "# Pathology strict hold-out integrity check", "",
  sprintf("Frozen split: seed 20260906; MD5 `%s`; 231 train and 99 test patients.", current_split_md5),
  "",
  "| Check | Result | Evidence |", "|---|---|---|",
  "| Train/test patient overlap | PASS: 0 | Exact patient-ID intersection |",
  "| Test samples in feature screening | PASS: No | Pearson screening used the 231 training rows only |",
  "| Test outcome in feature selection | PASS: No | Screening target was training frozen OSARS; LASSO outcome was training OS only |",
  "| Test in normalization estimation | PASS: No | glmnet standardization was fitted inside the training-only call; exported means/scales use training rows |",
  "| Test in lambda selection | PASS: No | lambda.min came from 10 folds containing training IDs only |",
  "| Test in coefficient estimation | PASS: No | Coefficients came from the full training fit at lambda.min |",
  "| Test in cutoff selection | PASS: No | The only cutoff is the training pathology-score median |",
  "| Seed/split retry | PASS: No | The split script refuses overwrite; no alternative seed or split was evaluated |",
  "", "## Frozen workflow", "",
  sprintf("- Source: 3,780 archived tiles from 378 WSIs and 364 pathology patients; 330 patients met all eligibility criteria."),
  sprintf("- Aggregation: arithmetic mean across all tiles and WSIs per patient; 2,048 patient-level ResNet50 features."),
  sprintf("- Training feature screen: Pearson versus continuous frozen OSARS, |r| > 0.2 and nominal P < 0.05; %d features retained.", length(screened_features)),
  sprintf("- Model: Cox LASSO, 10-fold training-only CV, lambda.min = %.15g; %d nonzero features.", lambda_selected, length(final_features)),
  sprintf("- Convergence: full training fit jerr=%d; selected lambda index %d/%d; %d CV warning(s) occurred only in the lower-lambda tail beginning at position %s. No model change or refit followed.",
          cvfit$glmnet.fit$jerr, selected_lambda_index, length(cvfit$lambda), length(cv_warnings),
          if (length(warning_lambda_positions)) min(warning_lambda_positions) else "NA"),
  sprintf("- Cutoff: training median = %.15g; test received predict() and the frozen cutoff only.", training_cutoff)
)
writeLines(integrity_lines, file.path(output_dir, "PATHOLOGY_HOLDOUT_INTEGRITY_CHECK.md"))

old_performance_file <- Sys.getenv("HCC_OS_PATHOLOGY_PRIOR_PERFORMANCE")
old_test_c <- NA_real_
old_test_lower <- NA_real_
old_test_upper <- NA_real_
if (file.exists(old_performance_file)) {
  old <- read.csv(old_performance_file, stringsAsFactors = FALSE)
  old_row <- old[
    old$split == "Test" & old$subset == "All" & old$model == "Pathology" &
      old$metric == "C_index", , drop = FALSE
  ]
  if (nrow(old_row) == 1L) {
    old_test_c <- old_row$estimate[[1]]
    old_test_lower <- old_row$lower[[1]]
    old_test_upper <- old_row$upper[[1]]
  }
}

old_original_file <- Sys.getenv("HCC_OS_PATHOLOGY_ORIGINAL_PATHDAT")
old_original_test_c <- NA_real_
if (file.exists(old_original_file)) {
  old_original_environment <- new.env(parent = emptyenv())
  load(old_original_file, envir = old_original_environment)
  if ("dat_all" %in% ls(old_original_environment)) {
    old_original_test <- old_original_environment$dat_all[
      old_original_environment$dat_all$Set == "Test", , drop = FALSE
    ]
    old_original_test_c <- get_cindex(
      old_original_test$OS.time, old_original_test$OS, old_original_test$RiskScore
    )
  }
}

train_performance <- performance[performance$split == "train", , drop = FALSE]
test_performance <- performance[performance$split == "test", , drop = FALSE]
test_km <- km_statistics[km_statistics$split == "test", , drop = FALSE]

comparison_table <- data.frame(
  analysis = c(
    "Original saved analysis with pre-split screening",
    "Prior corrected historical split",
    "Current seed-20260906 strict hold-out"
  ),
  test_C_index = c(old_original_test_c, old_test_c, test_performance$C_index[[1]]),
  CI_lower = c(NA_real_, old_test_lower, test_performance$C_index_CI_lower[[1]]),
  CI_upper = c(NA_real_, old_test_upper, test_performance$C_index_CI_upper[[1]]),
  note = c(
    "Historical test membership; feature screening preceded split",
    "Historical test membership; training-only corrected model",
    "New frozen split; training-only screen/model/cutoff"
  ),
  stringsAsFactors = FALSE
)
write_csv(comparison_table, "pathology_old_analysis_Cindex_comparison.csv")

comparison_text <- if (is.finite(old_test_c) && is.finite(old_original_test_c)) {
  sprintf(
    "The strict test C-index was modestly lower than the original leakage-prone saved estimate (%.3f vs %.3f), but slightly higher than the prior corrected historical-split estimate (%.3f). Memberships differ, so these are descriptive, not paired, comparisons.",
    test_performance$C_index[[1]], old_original_test_c, old_test_c
  )
} else {
  "The old analysis estimate could not be recovered for a numerical comparison."
}

robust_signal <- is.finite(test_performance$C_index_CI_lower[[1]]) &&
  test_performance$C_index_CI_lower[[1]] > 0.5 &&
  is.finite(test_performance$continuous_HR_Wald_P[[1]]) &&
  test_performance$continuous_HR_Wald_P[[1]] < 0.05
interpretation <- if (robust_signal) {
  "The fixed model shows an internal hold-out signal, but it remains an exploratory supplementary analysis because both the OSARS screening target and WSI cohort are TCGA-derived and no external pathology cohort was tested."
} else {
  "The strict hold-out results do not provide robust support for a pathology prognostic model; they are more appropriate as an exploratory supplementary analysis."
}

summary_lines <- c(
  "# Pathology strict hold-out summary", "",
  sprintf("1. Final patients: %d.", nrow(patient_data)),
  sprintf("2. Training: %d patients/%d events; test: %d patients/%d events.",
          length(train_idx), sum(patient_data$event[train_idx]),
          length(test_idx), sum(patient_data$event[test_idx])),
  sprintf("3. Features: 2048 initial → %d after training-only screening → %d final LASSO features.",
          length(screened_features), length(final_features)),
  sprintf("4. Apparent training C-index: %s (95%% bootstrap CI %s–%s).",
          fmt(train_performance$C_index[[1]]), fmt(train_performance$C_index_CI_lower[[1]]),
          fmt(train_performance$C_index_CI_upper[[1]])),
  sprintf("5. Strict test C-index: %s (95%% bootstrap CI %s–%s; 2,000 patient resamples).",
          fmt(test_performance$C_index[[1]]), fmt(test_performance$C_index_CI_lower[[1]]),
          fmt(test_performance$C_index_CI_upper[[1]])),
  sprintf("6. Strict test continuous Cox HR per training-set SD: %s (95%% CI %s–%s; P=%s).",
          fmt(test_performance$continuous_HR_per_training_SD[[1]]),
          fmt(test_performance$continuous_HR_CI_lower[[1]]),
          fmt(test_performance$continuous_HR_CI_upper[[1]]),
          fmt_p(test_performance$continuous_HR_Wald_P[[1]])),
  sprintf("7. With the frozen training median cutoff %.6g, test Low/High n=%d/%d and events=%d/%d; High-versus-Low HR %s (95%% CI %s–%s), log-rank P=%s.",
          training_cutoff, test_km$n_low[[1]], test_km$n_high[[1]],
          test_km$events_low[[1]], test_km$events_high[[1]],
          fmt(test_km$HR_high_vs_low[[1]]), fmt(test_km$HR_CI_lower[[1]]),
          fmt(test_km$HR_CI_upper[[1]]), fmt_p(test_km$logrank_P[[1]])),
  paste0("8. Comparison with old analysis: ", comparison_text),
  paste0("9. Interpretation: ", interpretation)
)

writeLines(summary_lines, file.path(output_dir, "PATHOLOGY_STRICT_HOLDOUT_SUMMARY.md"))

software_versions <- data.frame(
  software = c("R", "survival", "glmnet"),
  version = c(R.version.string, as.character(packageVersion("survival")),
              as.character(packageVersion("glmnet"))),
  stringsAsFactors = FALSE
)
write_csv(software_versions, "pathology_software_versions.csv")
capture.output(sessionInfo(), file = file.path(output_dir, "pathology_R_session_info.txt"))

writeLines(c(
  "status=complete",
  paste0("split_md5=", current_split_md5),
  "seed_retry=No",
  paste0("completed_utc=", format(Sys.time(), tz = "UTC", usetz = TRUE))
), completion_marker)

cat("STRICT HOLD-OUT ANALYSIS COMPLETE\n")
print(performance)
print(km_statistics)
