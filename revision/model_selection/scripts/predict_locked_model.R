#!/usr/bin/env Rscript

args_all <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args_all[grepl("^--file=", args_all)])
script_path <- normalizePath(script_arg[[1]], mustWork = TRUE)
analysis_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
source(file.path(analysis_root, "scripts", "common.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) {
  stop("Usage: predict_locked_model.R <expression.csv> [output.csv] [locked_model.rds]")
}
expression_file <- normalizePath(args[[1]], mustWork = TRUE)
output_file <- if (length(args) >= 2L) args[[2]] else file.path(getwd(), "OSARS_locked_risk_scores.csv")
model_file <- if (length(args) >= 3L) args[[3]] else file.path(processed_root, "final_OSARS_model_locked.rds")
if (!file.exists(model_file)) stop("Locked model not found: ", model_file)

bundle <- readRDS(model_file)
input <- read.csv(expression_file, check.names = FALSE)
if (!"ID" %in% names(input)) stop("Expression file must contain an ID column")
score <- predict_locked_bundle(bundle, input[, setdiff(names(input), "ID"), drop = FALSE])
write.csv(data.frame(ID = input$ID, OSARS_risk_score = score, stringsAsFactors = FALSE), output_file, row.names = FALSE)
cat("Wrote ", normalizePath(output_file, mustWork = FALSE), "\n", sep = "")
