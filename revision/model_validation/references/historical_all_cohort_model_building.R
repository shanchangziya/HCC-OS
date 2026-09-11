#!/usr/bin/env Rscript

# Historical provenance only. This reconstructs the original short entry point
# in which TCGA-LIHC, ICGC-LIRI, and GSE14520 were all supplied to Mime1 during
# model development/ranking. It is not a leakage-controlled validation workflow
# and must not be used to report prospective or untouched external performance.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop(
    "Usage: Rscript historical_all_cohort_model_building.R ",
    "INPUT_QWEN0208_RDATA OUTPUT_RES_RDATA"
  )
}

input_file <- normalizePath(args[[1]], mustWork = TRUE)
output_file <- normalizePath(args[[2]], mustWork = FALSE)

suppressPackageStartupMessages(library(Mime1))
load(input_file)

if (!exists("list_train_vali_Data") || !exists("g")) {
  stop("Input must contain list_train_vali_Data and candidate-gene vector g")
}
if (length(list_train_vali_Data) != 3L) {
  stop("Historical entry point expects exactly three cohort data frames")
}

names(list_train_vali_Data) <- c("TCGA-LIHC", "ICGC-LIRI", "GSE14520")
res <- ML.Dev.Prog.Sig(
  train_data = list_train_vali_Data$`TCGA-LIHC`,
  list_train_vali_Data = list_train_vali_Data,
  unicox.filter.for.candi = TRUE,
  unicox_p_cutoff = 0.01,
  candidate_genes = g,
  mode = "all",
  nodesize = 10,
  seed = 5201314
)
save(res, file = output_file)
