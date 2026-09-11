suppressPackageStartupMessages(library(survival))

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]
analysis_dir <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg)), ".."), mustWork = FALSE)
revision_dir <- normalizePath(file.path(analysis_dir, ".."), mustWork = FALSE)
input_dir <- file.path(revision_dir, "model_validation/data/processed")
out_dir <- file.path(analysis_dir, "results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cohorts <- c("TCGA-LIHC", "ICGC-LIRI")
result_rows <- list()
ph_rows <- list()
linearity_rows <- list()
audit_rows <- list()

for (cohort in cohorts) {
  path <- file.path(input_dir, paste0(cohort, "_clinical_analysis_rows.csv"))
  d0 <- read.csv(path, check.names = FALSE)
  required <- c("ID", "time", "event", "OSARS_SD", "stage_advanced", "grade_high")
  stopifnot(all(required %in% colnames(d0)), !anyDuplicated(d0$ID))

  complete <- complete.cases(d0[, required])
  d <- d0[complete, required]
  stopifnot(all(d$time > 0), all(d$event %in% c(0, 1)))
  stopifnot(all(d$stage_advanced %in% c(0, 1)), all(d$grade_high %in% c(0, 1)))

  expected <- if (cohort == "TCGA-LIHC") c(n = 309, events = 104) else c(n = 172, events = 29)
  stopifnot(nrow(d) == expected[["n"]], sum(d$event) == expected[["events"]])

  fit <- coxph(
    Surv(time, event) ~ OSARS_SD + stage_advanced + grade_high,
    data = d,
    ties = "efron",
    x = TRUE,
    y = TRUE,
    model = TRUE
  )
  s <- summary(fit)
  coef_table <- s$coefficients
  ci_table <- s$conf.int
  result_rows[[cohort]] <- data.frame(
    cohort = cohort,
    term = rownames(coef_table),
    beta = coef_table[, "coef"],
    SE = coef_table[, "se(coef)"],
    HR = ci_table[, "exp(coef)"],
    ci_lower = ci_table[, "lower .95"],
    ci_upper = ci_table[, "upper .95"],
    z = coef_table[, "z"],
    p = coef_table[, "Pr(>|z|)"],
    n = nrow(d),
    events = sum(d$event),
    formula = "Surv(time,event) ~ OSARS_SD + stage_advanced + grade_high",
    ties = "efron",
    stringsAsFactors = FALSE
  )

  ph <- cox.zph(fit, transform = "km")$table
  ph_rows[[cohort]] <- data.frame(
    cohort = cohort,
    term = rownames(ph),
    chisq = ph[, "chisq"],
    df = ph[, "df"],
    p = ph[, "p"],
    transform = "km",
    stringsAsFactors = FALSE
  )

  fit_spline <- coxph(
    Surv(time, event) ~ splines::ns(OSARS_SD, df = 3) + stage_advanced + grade_high,
    data = d,
    ties = "efron"
  )
  lr <- 2 * (fit_spline$loglik[2] - fit$loglik[2])
  linearity_rows[[cohort]] <- data.frame(
    cohort = cohort,
    LR_chisq = lr,
    df = 2,
    p = pchisq(lr, df = 2, lower.tail = FALSE),
    comparison = "3-df natural spline versus linear OSARS term, both adjusted for stage and grade",
    stringsAsFactors = FALSE
  )

  audit_rows[[cohort]] <- data.frame(
    cohort = cohort,
    n_total = nrow(d0),
    events_total = sum(d0$event),
    n_complete = nrow(d),
    events_complete = sum(d$event),
    n_excluded = nrow(d0) - nrow(d),
    stage_definition = "stage I-II vs III-IV",
    grade_definition = "grade I-II/G1-G2 vs III-IV/G3-G4; ambiguous interval labels retained as missing",
    OSARS_scale = "per cohort SD, using the frozen project definition",
    stringsAsFactors = FALSE
  )
}

write.csv(do.call(rbind, result_rows), file.path(out_dir, "multivariable_Cox_stage_grade.csv"), row.names = FALSE)
write.csv(do.call(rbind, ph_rows), file.path(out_dir, "multivariable_Cox_stage_grade_PH.csv"), row.names = FALSE)
write.csv(do.call(rbind, linearity_rows), file.path(out_dir, "multivariable_Cox_stage_grade_linearity.csv"), row.names = FALSE)
write.csv(do.call(rbind, audit_rows), file.path(out_dir, "multivariable_Cox_stage_grade_cohort_audit.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "R_sessionInfo.txt"))
