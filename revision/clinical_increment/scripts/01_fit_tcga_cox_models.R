suppressPackageStartupMessages(library(survival))

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]
analysis_dir <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg)), ".."), mustWork = FALSE)
revision_dir <- normalizePath(file.path(analysis_dir, ".."), mustWork = FALSE)
input_file <- file.path(revision_dir, "model_validation/data/processed", "TCGA-LIHC_clinical_analysis_rows.csv")
out_dir <- file.path(analysis_dir, "data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

d <- read.csv(input_file, check.names = FALSE)
required <- c("ID", "time", "event", "stage_num", "ROS_mean_rank", "OSARS_frozen")
stopifnot(all(required %in% colnames(d)), !anyDuplicated(d$ID))
d <- d[complete.cases(d[, required]), required]
stopifnot(nrow(d) == 321L, sum(d$event) == 110L)
stopifnot(all(d$stage_num %in% 1:4), all(d$time > 0), all(d$event %in% c(0, 1)))

forms <- list(
  "ROS only" = Surv(time, event) ~ ROS_mean_rank,
  "OSARS only" = Surv(time, event) ~ OSARS_frozen,
  "Stage" = Surv(time, event) ~ stage_num,
  "Stage + ROS" = Surv(time, event) ~ stage_num + ROS_mean_rank,
  "Stage + OSARS" = Surv(time, event) ~ stage_num + OSARS_frozen,
  "Stage + ROS + OSARS" = Surv(time, event) ~ stage_num + ROS_mean_rank + OSARS_frozen
)

prediction_rows <- list()
coefficient_rows <- list()
fit_rows <- list()

for (model_name in names(forms)) {
  fit <- coxph(forms[[model_name]], data = d, ties = "efron", x = TRUE, y = TRUE)
  lp <- predict(fit, newdata = d, type = "lp", reference = "zero")
  c_stat <- concordance(Surv(d$time, d$event) ~ lp, reverse = TRUE)$concordance

  prediction_rows[[model_name]] <- data.frame(
    ID = d$ID,
    time = d$time,
    event = d$event,
    model = model_name,
    linear_predictor = as.numeric(lp),
    evaluation = "TCGA apparent",
    stringsAsFactors = FALSE
  )
  coefficient_rows[[model_name]] <- data.frame(
    model = model_name,
    term = names(coef(fit)),
    coefficient = unname(coef(fit)),
    stringsAsFactors = FALSE
  )
  fit_rows[[model_name]] <- data.frame(
    model = model_name,
    formula = paste(deparse(forms[[model_name]]), collapse = " "),
    n = nrow(d),
    events = sum(d$event),
    apparent_Harrell_C = as.numeric(c_stat),
    ties = "efron",
    stringsAsFactors = FALSE
  )
}

write.csv(d, file.path(out_dir, "TCGA_complete_case_analysis_cohort.csv"), row.names = FALSE)
write.csv(
  do.call(rbind, prediction_rows),
  file.path(out_dir, "TCGA_stage_ROS_OSARS_model_predictions.csv"),
  row.names = FALSE
)
write.csv(
  do.call(rbind, coefficient_rows),
  file.path(out_dir, "TCGA_stage_ROS_OSARS_model_coefficients.csv"),
  row.names = FALSE
)
write.csv(
  do.call(rbind, fit_rows),
  file.path(out_dir, "TCGA_stage_ROS_OSARS_model_fits.csv"),
  row.names = FALSE
)

writeLines(
  capture.output(sessionInfo()),
  file.path(out_dir, "R_sessionInfo.txt")
)
