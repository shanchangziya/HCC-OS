options(stringsAsFactors = FALSE)

if (!exists("analysis_root", inherits = FALSE)) {
  stop("analysis_root must be defined before sourcing common.R")
}

analysis_root <- normalizePath(analysis_root, mustWork = TRUE)
repo_root <- normalizePath(file.path(analysis_root, "../.."), mustWork = TRUE)
project_root <- repo_root
raw_root <- Sys.getenv(
  "HCC_OS_MODEL_RAW_DIR",
  file.path(repo_root, "revision", "model_validation", "data", "raw")
)
processed_root <- file.path(analysis_root, "data", "processed")
logs_root <- file.path(analysis_root, "logs")
figures_root <- file.path(analysis_root, "figures")
local_lib <- file.path(analysis_root, ".Rlib")
extra_lib <- Sys.getenv("HCC_OS_R_LIBRARY")
candidate_libs <- c(local_lib, extra_lib)
candidate_libs <- candidate_libs[nzchar(candidate_libs) & dir.exists(candidate_libs)]
.libPaths(unique(c(candidate_libs, .libPaths())))

master_seed <- 5201314L
outer_folds <- 5L
univariable_p_cutoff <- 0.01
near_zero_variance_cutoff <- 1e-6

required_model_packages <- c(
  "survival", "randomForestSRC", "glmnet", "plsRcox", "superpc",
  "gbm", "CoxBoost", "survivalsvm", "dplyr", "tibble",
  "miscTools", "compareC", "mixOmics"
)

assert_packages <- function() {
  missing <- required_model_packages[
    !vapply(required_model_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
  ]
  if (length(missing)) {
    stop("Missing R packages: ", paste(missing, collapse = ", "))
  }
  invisible(TRUE)
}

generate_configurations <- function() {
  out <- list()
  add <- function(configuration, mode, single_ml = NA_character_,
                  double_ml1 = NA_character_, double_ml2 = NA_character_,
                  alpha = NA_real_, direction = NA_character_) {
    out[[length(out) + 1L]] <<- data.frame(
      configuration_id = length(out) + 1L,
      configuration = configuration,
      mode = mode,
      single_ml = single_ml,
      double_ml1 = double_ml1,
      double_ml2 = double_ml2,
      alpha = alpha,
      direction = direction,
      stage_count = ifelse(mode == "single", 1L, 2L),
      stringsAsFactors = FALSE
    )
  }

  add("RSF", "single", single_ml = "RSF")
  add("RSF + CoxBoost", "double", double_ml1 = "RSF", double_ml2 = "CoxBoost")
  for (a in seq(0.1, 0.9, 0.1)) {
    add(sprintf("RSF + Enet[α=%.1f]", a), "double", double_ml1 = "RSF", double_ml2 = "Enet", alpha = a)
  }
  add("RSF + GBM", "double", double_ml1 = "RSF", double_ml2 = "GBM")
  add("RSF + Lasso", "double", double_ml1 = "RSF", double_ml2 = "Lasso")
  add("RSF + plsRcox", "double", double_ml1 = "RSF", double_ml2 = "plsRcox")
  add("RSF + Ridge", "double", double_ml1 = "RSF", double_ml2 = "Ridge")
  for (d in c("both", "backward", "forward")) {
    add(sprintf("RSF + StepCox[%s]", d), "double", double_ml1 = "RSF", double_ml2 = "StepCox", direction = d)
  }
  add("RSF + SuperPC", "double", double_ml1 = "RSF", double_ml2 = "superpc")
  add("RSF + survival-SVM", "double", double_ml1 = "RSF", double_ml2 = "survivalsvm")

  for (a in seq(0.1, 0.9, 0.1)) {
    add(sprintf("Enet[α=%.1f]", a), "single", single_ml = "Enet", alpha = a)
  }

  for (d1 in c("both", "backward", "forward")) {
    add(sprintf("StepCox[%s]", d1), "single", single_ml = "StepCox", direction = d1)
  }
  for (d1 in c("both", "backward", "forward")) {
    add(sprintf("StepCox[%s] + CoxBoost", d1), "double", double_ml1 = "StepCox", double_ml2 = "CoxBoost", direction = d1)
    for (a in seq(0.1, 0.9, 0.1)) {
      add(sprintf("StepCox[%s] + Enet[α=%.1f]", d1, a), "double", double_ml1 = "StepCox", double_ml2 = "Enet", alpha = a, direction = d1)
    }
    for (second in c("GBM", "Lasso", "plsRcox", "Ridge", "RSF", "superpc", "survivalsvm")) {
      label <- switch(second, superpc = "SuperPC", survivalsvm = "survival-SVM", second)
      add(sprintf("StepCox[%s] + %s", d1, label), "double", double_ml1 = "StepCox", double_ml2 = second, direction = d1)
    }
  }

  add("CoxBoost", "single", single_ml = "CoxBoost")
  for (a in seq(0.1, 0.9, 0.1)) {
    add(sprintf("CoxBoost + Enet[α=%.1f]", a), "double", double_ml1 = "CoxBoost", double_ml2 = "Enet", alpha = a)
  }
  for (second in c("GBM", "Lasso", "plsRcox", "Ridge")) {
    add(sprintf("CoxBoost + %s", second), "double", double_ml1 = "CoxBoost", double_ml2 = second)
  }
  for (d2 in c("both", "backward", "forward")) {
    add(sprintf("CoxBoost + StepCox[%s]", d2), "double", double_ml1 = "CoxBoost", double_ml2 = "StepCox", direction = d2)
  }
  add("CoxBoost + SuperPC", "double", double_ml1 = "CoxBoost", double_ml2 = "superpc")
  add("CoxBoost + survival-SVM", "double", double_ml1 = "CoxBoost", double_ml2 = "survivalsvm")

  add("plsRcox", "single", single_ml = "plsRcox")
  add("SuperPC", "single", single_ml = "superpc")
  add("GBM", "single", single_ml = "GBM")
  add("survival - SVM", "single", single_ml = "survivalsvm")
  add("Ridge", "single", single_ml = "Ridge")
  add("Lasso", "single", single_ml = "Lasso")

  for (second in c("CoxBoost", "GBM", "plsRcox", "RSF")) {
    add(sprintf("Lasso + %s", second), "double", double_ml1 = "Lasso", double_ml2 = second)
  }
  for (d2 in c("both", "backward", "forward")) {
    add(sprintf("Lasso + StepCox[%s]", d2), "double", double_ml1 = "Lasso", double_ml2 = "StepCox", direction = d2)
  }
  add("Lasso + SuperPC", "double", double_ml1 = "Lasso", double_ml2 = "superpc")
  add("Lasso + survival-SVM", "double", double_ml1 = "Lasso", double_ml2 = "survivalsvm")

  configs <- do.call(rbind, out)
  stopifnot(nrow(configs) == 117L, sum(configs$mode == "single") == 20L, sum(configs$mode == "double") == 97L)
  if (anyDuplicated(configs$configuration)) stop("Configuration names are not unique")
  configs
}

load_mime_core <- function() {
  if (exists("ML.Dev.Prog.Sig", mode = "function", inherits = TRUE)) return(invisible(TRUE))
  assert_packages()
  source_file <- file.path(analysis_root, "references", "Mime_ML.Dev.Prog.Sig_upstream_20250923.R")
  txt <- readLines(source_file, warn = FALSE)

  # These packages are loaded upstream but are not referenced by the modelling function.
  unused <- "^\\s*library\\((BART|ggplot2|ggsci|tidyr|ggbreak|data\\.table)\\)\\s*$"
  txt <- txt[!grepl(unused, txt)]

  # Preserve the post-first-stage predictor set for fold-wise complexity auditing.
  replacement <- paste0(
    "'Sig.genes' = pre_var, 'Input.features' = ",
    "if (exists('rid', inherits = FALSE)) as.character(rid) else as.character(pre_var))"
  )
  txt <- gsub("'Sig.genes' = pre_var)", replacement, txt, fixed = TRUE)
  txt <- gsub("'Sig.genes'=pre_var)", replacement, txt, fixed = TRUE)

  # Avoid an infinite retry loop if SuperPC cannot be fitted for a particular fold/configuration.
  txt <- gsub("repeat {", "for (.mime_retry in seq_len(3L)) {", txt, fixed = TRUE)
  eval(parse(text = txt, keep.source = FALSE), envir = .GlobalEnv)
  if (!exists("ML.Dev.Prog.Sig", mode = "function", inherits = TRUE)) stop("Failed to load Mime modelling core")
  invisible(TRUE)
}

sha256_file <- function(path) {
  out <- system2("shasum", c("-a", "256", shQuote(path)), stdout = TRUE, stderr = TRUE)
  if (!length(out)) return(NA_character_)
  strsplit(out[[1]], "[[:space:]]+")[[1]][1]
}

stratified_folds <- function(status, k = 5L, seed = master_seed) {
  set.seed(seed)
  fold <- integer(length(status))
  for (s in sort(unique(status))) {
    idx <- sample(which(status == s))
    fold[idx] <- rep(seq_len(k), length.out = length(idx))
  }
  # Shuffle fold labels once so no status group is systematically assigned first.
  map <- sample(seq_len(k))
  map[fold]
}

fit_preprocessor <- function(x, variance_cutoff = near_zero_variance_cutoff) {
  x <- as.data.frame(x, check.names = FALSE)
  q01 <- vapply(x, stats::quantile, probs = 0.01, na.rm = TRUE, names = FALSE, FUN.VALUE = numeric(1))
  q99 <- vapply(x, stats::quantile, probs = 0.99, na.rm = TRUE, names = FALSE, FUN.VALUE = numeric(1))
  capped <- x
  for (j in seq_along(capped)) {
    v <- as.numeric(capped[[j]])
    v[!is.finite(v)] <- NA_real_
    v <- pmax(q01[j], pmin(q99[j], v))
    capped[[j]] <- v
  }
  means <- vapply(capped, mean, na.rm = TRUE, FUN.VALUE = numeric(1))
  means[!is.finite(means)] <- 0
  for (j in seq_along(capped)) {
    v <- capped[[j]]
    v[is.na(v)] <- means[j]
    capped[[j]] <- v
  }
  variances <- vapply(capped, stats::var, na.rm = TRUE, FUN.VALUE = numeric(1))
  keep <- is.finite(variances) & variances >= variance_cutoff
  list(q01 = q01, q99 = q99, means = means, variances = variances,
       keep = names(x)[keep], removed = names(x)[!keep])
}

apply_preprocessor <- function(x, prep) {
  x <- as.data.frame(x, check.names = FALSE)
  missing <- setdiff(names(prep$q01), names(x))
  if (length(missing)) stop("Missing expression features: ", paste(missing, collapse = ", "))
  x <- x[, names(prep$q01), drop = FALSE]
  for (j in seq_along(x)) {
    v <- as.numeric(x[[j]])
    v[!is.finite(v)] <- NA_real_
    v <- pmax(prep$q01[j], pmin(prep$q99[j], v))
    v[is.na(v)] <- prep$means[j]
    x[[j]] <- v
  }
  x[, prep$keep, drop = FALSE]
}

univariable_cox_screen <- function(x, time, status, p_cutoff = univariable_p_cutoff) {
  stopifnot(nrow(x) == length(time), length(time) == length(status))
  rows <- lapply(names(x), function(gene) {
    v <- as.numeric(x[[gene]])
    fit <- tryCatch(
      survival::coxph(survival::Surv(time, status) ~ v, ties = "efron"),
      error = function(e) NULL
    )
    if (is.null(fit)) return(data.frame(gene = gene, HR = NA_real_, z = NA_real_, pvalue = NA_real_))
    sm <- summary(fit)
    data.frame(
      gene = gene,
      HR = unname(sm$coefficients[1, "exp(coef)"]),
      z = unname(sm$coefficients[1, "z"]),
      pvalue = unname(sm$coefficients[1, "Pr(>|z|)"])
    )
  })
  table <- do.call(rbind, rows)
  selected <- table$gene[is.finite(table$pvalue) & table$pvalue < p_cutoff]
  list(table = table, selected = selected)
}

harrell_c <- function(time, status, risk) {
  ok <- is.finite(time) & !is.na(status) & is.finite(risk)
  if (sum(ok) < 3L || length(unique(risk[ok])) < 2L) return(NA_real_)
  fit <- tryCatch(
    survival::concordance(
      survival::Surv(time[ok], status[ok]) ~ risk[ok],
      reverse = TRUE,
      timefix = TRUE
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) NA_real_ else as.numeric(fit$concordance)
}

extract_score <- function(result, cohort, expected_ids) {
  if (is.null(result$riskscore) || !length(result$riskscore)) stop("No risk scores returned")
  model_scores <- result$riskscore[[1]]
  if (!cohort %in% names(model_scores)) stop("Missing cohort score: ", cohort)
  tab <- as.data.frame(model_scores[[cohort]])
  if (!all(c("ID", "RS") %in% names(tab))) stop("Malformed risk-score table")
  idx <- match(expected_ids, as.character(tab$ID))
  if (anyNA(idx)) stop("Risk-score IDs do not match input IDs")
  as.numeric(tab$RS[idx])
}

model_parameter_string <- function(configuration_row, result) {
  pieces <- c(
    paste0("mode=", configuration_row$mode),
    "outer_seed=5201314",
    "nodesize=10",
    "inner_folds=10"
  )
  if (is.finite(configuration_row$alpha)) pieces <- c(pieces, paste0("alpha=", configuration_row$alpha))
  if (!is.na(configuration_row$direction)) pieces <- c(pieces, paste0("direction=", configuration_row$direction))
  fit <- if (length(result$ml.res)) result$ml.res[[1]] else NULL
  if (is.list(fit) && !is.null(fit$best)) pieces <- c(pieces, paste0("gbm_best_trees=", fit$best))
  if (inherits(fit, "cv.glmnet")) pieces <- c(pieces, paste0("lambda.min=", signif(fit$lambda.min, 8)))
  if (inherits(fit, "CoxBoost")) {
    if (!is.null(fit$stepno)) pieces <- c(pieces, paste0("stepno=", fit$stepno))
    if (!is.null(fit$penalty)) pieces <- c(pieces, paste0("penalty=", signif(fit$penalty, 8)))
  }
  paste(pieces, collapse = ";")
}

run_mime_configuration <- function(config, train_data, heldout_data, candidate_genes,
                                   log_file, seed = master_seed, gbm_cores = 1L,
                                   keep_model = FALSE) {
  load_mime_core()
  args <- list(
    train_data = train_data,
    list_train_vali_Data = list(TCGA_training = train_data, TCGA_heldout = heldout_data),
    candidate_genes = candidate_genes,
    unicox.filter.for.candi = FALSE,
    unicox_p_cutoff = univariable_p_cutoff,
    mode = config$mode,
    nodesize = 10,
    seed = seed,
    cores_for_parallel = gbm_cores
  )
  if (config$mode == "single") {
    args$single_ml <- config$single_ml
  } else {
    args$double_ml1 <- config$double_ml1
    args$double_ml2 <- config$double_ml2
  }
  if (is.finite(config$alpha)) args$alpha_for_Enet <- config$alpha
  if (!is.na(config$direction)) args$direction_for_stepcox <- config$direction

  warnings <- character()
  con <- file(log_file, open = "at")
  sink(con, type = "output")
  sink(con, type = "message")
  on.exit({
    sink(type = "message")
    sink(type = "output")
    close(con)
  }, add = TRUE)
  cat("\n===== ", config$configuration, " | ", format(Sys.time()), " =====\n", sep = "")
  started <- Sys.time()
  result <- withCallingHandlers(
    do.call(ML.Dev.Prog.Sig, args),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))

  train_score_raw <- extract_score(result, "TCGA_training", train_data$ID)
  heldout_score_raw <- extract_score(result, "TCGA_heldout", heldout_data$ID)
  train_c_raw <- harrell_c(train_data$OS.time, train_data$OS, train_score_raw)
  orientation <- if (is.finite(train_c_raw) && train_c_raw < 0.5) -1 else 1
  train_score <- orientation * train_score_raw
  heldout_score <- orientation * heldout_score_raw

  input_features <- result$Input.features
  if (is.null(input_features) || !length(input_features)) input_features <- candidate_genes
  input_features <- intersect(as.character(input_features), candidate_genes)
  if (!length(input_features)) input_features <- candidate_genes

  output <- list(
    configuration_id = config$configuration_id,
    configuration = config$configuration,
    fold_cindex = harrell_c(heldout_data$OS.time, heldout_data$OS, heldout_score),
    training_cindex = harrell_c(train_data$OS.time, train_data$OS, train_score),
    orientation = orientation,
    heldout = data.frame(
      patient_id = heldout_data$ID,
      OS.time = heldout_data$OS.time,
      OS = heldout_data$OS,
      risk_score = heldout_score,
      raw_risk_score = heldout_score_raw,
      stringsAsFactors = FALSE
    ),
    input_features = input_features,
    input_feature_count = length(input_features),
    univariable_feature_count = length(candidate_genes),
    parameter_string = model_parameter_string(config, result),
    warnings = unique(warnings),
    elapsed_seconds = elapsed,
    upstream_model_name = names(result$ml.res)[1]
  )
  if (isTRUE(keep_model)) output$model_result <- result
  output
}

prepare_locked_model_matrix <- function(bundle, expression_data) {
  x <- apply_preprocessor(expression_data, bundle$preprocessor)
  missing_selected <- setdiff(bundle$univariable_selected_genes, names(x))
  if (length(missing_selected)) stop("Selected genes unavailable after preprocessing")
  x <- x[, bundle$univariable_selected_genes, drop = FALSE]
  names(x) <- gsub("-", ".", names(x), fixed = TRUE)
  missing_model <- setdiff(bundle$model_input_features, names(x))
  if (length(missing_model)) stop("Locked model features missing: ", paste(missing_model, collapse = ", "))
  x[, bundle$model_input_features, drop = FALSE]
}

predict_locked_bundle <- function(bundle, expression_data) {
  x <- prepare_locked_model_matrix(bundle, expression_data)
  config <- bundle$configuration_row
  kind <- if (config$mode == "single") config$single_ml else config$double_ml2
  fit <- bundle$fit
  dummy <- data.frame(OS.time = rep(1, nrow(x)), OS = rep(0, nrow(x)), x, check.names = FALSE)

  raw <- switch(
    kind,
    RSF = as.numeric(predict(fit, newdata = dummy)$predicted),
    Enet = as.numeric(predict(fit, type = "link", newx = as.matrix(x), s = fit$lambda.min)),
    StepCox = as.numeric(predict(fit, type = "risk", newdata = dummy)),
    CoxBoost = as.numeric(predict(
      fit, newdata = as.matrix(x), newtime = dummy$OS.time,
      newstatus = dummy$OS, type = "lp"
    )),
    plsRcox = as.numeric(predict(fit, type = "lp", newdata = x)),
    superpc = {
      model <- fit$fit
      cv.fit <- fit$cv.fit
      train_x <- bundle$training_model_data[, bundle$model_input_features, drop = FALSE]
      train <- list(
        x = t(train_x), y = bundle$training_model_data$OS.time,
        censoring.status = bundle$training_model_data$OS,
        featurenames = bundle$model_input_features
      )
      test <- list(
        x = t(x), y = rep(1, nrow(x)), censoring.status = rep(0, nrow(x)),
        featurenames = bundle$model_input_features
      )
      threshold <- cv.fit$thresholds[which.max(cv.fit[["scor"]][1, ])]
      as.numeric(superpc::superpc.predict(model, train, test, threshold = threshold, n.components = 1)$v.pred)
    },
    GBM = as.numeric(predict(fit$fit, dummy, n.trees = fit$best, type = "link")),
    survivalsvm = as.numeric(predict(fit, dummy)$predicted),
    Ridge = {
      if (is.list(fit) && !inherits(fit, "cv.glmnet") && !is.null(fit$cv.fit)) {
        as.numeric(predict(fit$fit, type = "response", newx = as.matrix(x), s = fit$cv.fit$lambda.min))
      } else {
        as.numeric(predict(fit, type = "response", newx = as.matrix(x), s = fit$lambda.min))
      }
    },
    Lasso = as.numeric(predict(fit, type = "response", newx = as.matrix(x), s = fit$lambda.min)),
    stop("Unsupported locked model kind: ", kind)
  )
  as.numeric(bundle$orientation) * raw
}

extract_linear_coefficients <- function(bundle) {
  config <- bundle$configuration_row
  kind <- if (config$mode == "single") config$single_ml else config$double_ml2
  fit <- bundle$fit
  coef_vector <- NULL
  if (kind %in% c("Enet", "Lasso") && inherits(fit, "cv.glmnet")) {
    coef_vector <- as.matrix(stats::coef(fit, s = fit$lambda.min))[, 1]
  } else if (kind == "Ridge") {
    if (is.list(fit) && !inherits(fit, "cv.glmnet") && !is.null(fit$cv.fit)) {
      coef_vector <- as.matrix(stats::coef(fit$fit, s = fit$cv.fit$lambda.min))[, 1]
    } else {
      coef_vector <- as.matrix(stats::coef(fit, s = fit$lambda.min))[, 1]
    }
  } else if (kind == "StepCox" && inherits(fit, "coxph")) {
    coef_vector <- stats::coef(fit)
  } else if (kind == "CoxBoost") {
    coef_vector <- as.numeric(stats::coef(fit))
    names(coef_vector) <- bundle$model_input_features
  }
  if (is.null(coef_vector)) return(NULL)
  coef_vector[is.finite(coef_vector) & abs(coef_vector) > 0]
}
