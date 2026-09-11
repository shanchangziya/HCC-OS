#!/usr/bin/env Rscript

# OSARS-high bulk biology in frozen TCGA-LIHC data
#
# Primary question: what bulk transcriptional state is represented by tumors
# classified as OSARS-high using the pre-existing frozen score/cutoff?
#
# This script does not refit OSARS, alter its cutoff, or use the prior
# single-cell-derived OS-high program. It performs (i) gene-level limma
# differential expression and (ii) competitive pathway testing with CAMERA,
# using all Hallmark and Reactome gene sets available in msigdbr. A stage-
# adjusted analysis is a prespecified sensitivity analysis.

suppressPackageStartupMessages({
  library(limma)
  library(msigdbr)
  library(GSVA)
})

script_arg <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", script_arg[grepl("^--file=", script_arg)][1])
# The local R launcher encodes spaces this way when passing --file.
script_path <- gsub("~\\+~", " ", script_path)
script_path <- normalizePath(script_path)
script_dir <- dirname(script_path)
analysis_dir <- normalizePath(file.path(script_dir, ".."))
revision_dir <- normalizePath(file.path(analysis_dir, ".."))
model_dir <- file.path(revision_dir, "model_validation")

dir.create(file.path(analysis_dir, "data", "processed"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(analysis_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(analysis_dir, "references"), recursive = TRUE, showWarnings = FALSE)

expr_path <- file.path(model_dir, "data", "raw", "TCGA_full_expression_logscale.csv.gz")
score_path <- file.path(model_dir, "data", "processed", "TCGA-LIHC_all_scores.csv")
clinical_path <- file.path(model_dir, "data", "processed", "TCGA-LIHC_clinical_analysis_rows.csv")

stopifnot(file.exists(expr_path), file.exists(score_path), file.exists(clinical_path))

read_expression <- function(path) {
  raw <- read.csv(gzfile(path), check.names = FALSE, stringsAsFactors = FALSE)
  stopifnot("gene" %in% colnames(raw))
  genes <- toupper(trimws(raw$gene))
  if (anyDuplicated(genes)) stop("Duplicated gene symbols after upper-casing.")
  mat <- as.matrix(raw[, setdiff(colnames(raw), "gene"), drop = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- genes
  mat
}

expr <- read_expression(expr_path)
scores <- read.csv(score_path, check.names = FALSE, stringsAsFactors = FALSE)
clinical <- read.csv(clinical_path, check.names = FALSE, stringsAsFactors = FALSE)

required_score_cols <- c("ID", "OSARS_frozen", "OSARS_frozen_group")
if (!all(required_score_cols %in% colnames(scores))) {
  stop("The frozen-score table is missing required columns.")
}
if (!all(c("ID", "stage") %in% colnames(clinical))) {
  stop("The clinical table is missing ID or stage.")
}

sample_ids <- colnames(expr)
score_idx <- match(sample_ids, scores$ID)
if (anyNA(score_idx)) stop("Expression samples could not all be matched to frozen scores.")
score_aligned <- scores[score_idx, required_score_cols, drop = FALSE]
stopifnot(identical(sample_ids, score_aligned$ID))

clinical_idx <- match(sample_ids, clinical$ID)
stage <- clinical$stage[clinical_idx]
stage <- trimws(as.character(stage))
stage[stage %in% c("", "NA", "[Not Available]", "not reported", "Not Reported")] <- NA_character_

group <- factor(score_aligned$OSARS_frozen_group, levels = c("Low", "High"))
if (anyNA(group) || any(table(group) == 0L)) stop("Frozen OSARS groups must contain both Low and High samples.")

manifest <- data.frame(
  sample_id = sample_ids,
  OSARS_frozen = score_aligned$OSARS_frozen,
  OSARS_frozen_group = as.character(group),
  stage = stage,
  stringsAsFactors = FALSE
)
write.csv(manifest, file.path(analysis_dir, "data", "processed", "TCGA_OSARS_frozen_sample_manifest.csv"), row.names = FALSE)

make_gene_sets <- function(tbl, collection_label) {
  keep <- !is.na(tbl$gs_name) & !is.na(tbl$gene_symbol) & nzchar(tbl$gene_symbol)
  tbl <- tbl[keep, c("gs_name", "gene_symbol"), drop = FALSE]
  tbl$gene_symbol <- toupper(trimws(tbl$gene_symbol))
  sets <- split(tbl$gene_symbol, tbl$gs_name)
  sets <- lapply(sets, unique)
  attr(sets, "collection") <- collection_label
  sets
}

filter_gene_sets <- function(sets, available_genes, min_size = 15L, max_size = 500L) {
  used <- lapply(sets, function(x) intersect(x, available_genes))
  coverage <- data.frame(
    pathway = names(sets),
    n_genes_defined = lengths(sets),
    n_genes_used = lengths(used),
    stringsAsFactors = FALSE
  )
  coverage$tested <- coverage$n_genes_used >= min_size & coverage$n_genes_used <= max_size
  coverage <- coverage[order(coverage$pathway), , drop = FALSE]
  keep_names <- coverage$pathway[coverage$tested]
  list(sets = used[keep_names], coverage = coverage)
}

# Gene sets are frozen by package/database version in the output metadata.
hallmark_raw <- msigdbr(species = "Homo sapiens", collection = "H")
reactome_raw <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME")
hallmark_all <- make_gene_sets(hallmark_raw, "Hallmark")
reactome_all <- make_gene_sets(reactome_raw, "Reactome")
hallmark <- filter_gene_sets(hallmark_all, rownames(expr))
reactome <- filter_gene_sets(reactome_all, rownames(expr))

set_coverage <- rbind(
  transform(hallmark$coverage, collection = "Hallmark"),
  transform(reactome$coverage, collection = "Reactome")
)
set_coverage <- set_coverage[, c("collection", "pathway", "n_genes_defined", "n_genes_used", "tested")]
write.csv(set_coverage, file.path(analysis_dir, "tables", "gene_set_coverage.csv"), row.names = FALSE)

write_membership <- function(sets, collection_label) {
  do.call(rbind, lapply(names(sets), function(nm) {
    data.frame(collection = collection_label, pathway = nm, gene = sets[[nm]], stringsAsFactors = FALSE)
  }))
}
membership <- rbind(
  write_membership(hallmark$sets, "Hallmark"),
  write_membership(reactome$sets, "Reactome")
)
write.csv(membership, file.path(analysis_dir, "references", "tested_gene_set_membership.csv"), row.names = FALSE)

run_camera <- function(expression, sets, design, contrast, coverage, collection_label, analysis_label) {
  indices <- ids2indices(sets, rownames(expression), remove.empty = TRUE)
  res <- camera(expression, index = indices, design = design, contrast = contrast, sort = FALSE)
  res$pathway <- rownames(res)
  res$collection <- collection_label
  res$analysis <- analysis_label
  cov <- coverage[match(res$pathway, coverage$pathway), c("n_genes_defined", "n_genes_used"), drop = FALSE]
  res$n_genes_defined <- cov$n_genes_defined
  res$n_genes_used <- cov$n_genes_used
  res$FDR <- p.adjust(res$PValue, method = "BH")
  res$signed_log10_FDR <- ifelse(res$Direction == "Up", 1, -1) * -log10(pmax(res$FDR, .Machine$double.xmin))
  res <- res[, c("analysis", "collection", "pathway", "NGenes", "n_genes_defined", "n_genes_used", "Direction", "PValue", "FDR", "signed_log10_FDR")]
  res[order(res$FDR, res$PValue, res$pathway), , drop = FALSE]
}

run_limma_de <- function(expression, design, contrast, coef_name) {
  fit <- lmFit(expression, design)
  fit <- contrasts.fit(fit, contrast)
  fit <- eBayes(fit, trend = TRUE)
  out <- topTable(fit, coef = coef_name, number = Inf, sort.by = "none")
  out$gene <- rownames(out)
  out <- out[, c("gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]
  out[order(out$adj.P.Val, out$P.Value, decreasing = FALSE), , drop = FALSE]
}

# Primary analysis: frozen OSARS-high vs frozen OSARS-low, without covariates.
design_primary <- model.matrix(~ 0 + group)
colnames(design_primary) <- levels(group)
contrast_primary <- matrix(c(-1, 1), ncol = 1, dimnames = list(colnames(design_primary), "High_vs_Low"))

de_primary <- run_limma_de(expr, design_primary, contrast_primary, "High_vs_Low")
write.csv(de_primary, file.path(analysis_dir, "tables", "TCGA_OSARS_high_vs_low_limma_DE.csv"), row.names = FALSE)

camera_primary_h <- run_camera(expr, hallmark$sets, design_primary, contrast_primary[, 1], hallmark$coverage, "Hallmark", "primary_unadjusted")
camera_primary_r <- run_camera(expr, reactome$sets, design_primary, contrast_primary[, 1], reactome$coverage, "Reactome", "primary_unadjusted")
camera_primary <- rbind(camera_primary_h, camera_primary_r)
write.csv(camera_primary, file.path(analysis_dir, "tables", "TCGA_OSARS_high_vs_low_CAMERA_all.csv"), row.names = FALSE)
write.csv(camera_primary_h, file.path(analysis_dir, "tables", "TCGA_OSARS_high_vs_low_CAMERA_hallmark.csv"), row.names = FALSE)
write.csv(camera_primary_r, file.path(analysis_dir, "tables", "TCGA_OSARS_high_vs_low_CAMERA_reactome.csv"), row.names = FALSE)

# Sensitivity analysis: preserve the frozen group while adjusting for clinical stage.
stage_keep <- !is.na(stage) & !is.na(group)
group_stage <- droplevels(group[stage_keep])
stage_factor <- factor(stage[stage_keep])
design_stage <- model.matrix(~ group_stage + stage_factor)
group_coef <- grep("^group_stageHigh$", colnames(design_stage), value = TRUE)
if (length(group_coef) != 1L) stop("Could not identify the OSARS group coefficient in the stage-adjusted design.")
contrast_stage <- rep(0, ncol(design_stage))
names(contrast_stage) <- colnames(design_stage)
contrast_stage[group_coef] <- 1
contrast_stage_mat <- matrix(contrast_stage, ncol = 1, dimnames = list(colnames(design_stage), "High_vs_Low_stage_adjusted"))

de_stage <- run_limma_de(expr[, stage_keep, drop = FALSE], design_stage, contrast_stage_mat, "High_vs_Low_stage_adjusted")
write.csv(de_stage, file.path(analysis_dir, "tables", "TCGA_OSARS_high_vs_low_stage_adjusted_limma_DE.csv"), row.names = FALSE)

camera_stage_h <- run_camera(expr[, stage_keep, drop = FALSE], hallmark$sets, design_stage, contrast_stage, hallmark$coverage, "Hallmark", "stage_adjusted")
camera_stage_r <- run_camera(expr[, stage_keep, drop = FALSE], reactome$sets, design_stage, contrast_stage, reactome$coverage, "Reactome", "stage_adjusted")
camera_stage <- rbind(camera_stage_h, camera_stage_r)
write.csv(camera_stage, file.path(analysis_dir, "tables", "TCGA_OSARS_high_vs_low_stage_adjusted_CAMERA_all.csv"), row.names = FALSE)
write.csv(camera_stage_h, file.path(analysis_dir, "tables", "TCGA_OSARS_high_vs_low_stage_adjusted_CAMERA_hallmark.csv"), row.names = FALSE)
write.csv(camera_stage_r, file.path(analysis_dir, "tables", "TCGA_OSARS_high_vs_low_stage_adjusted_CAMERA_reactome.csv"), row.names = FALSE)

# Sample-level Hallmark scores are included solely for transparent descriptive
# visualization of pathways selected from the all-pathway CAMERA result. The
# inferential pathway result remains the competitive CAMERA test above.
gsva_hallmark <- suppressWarnings(gsva(expr, hallmark$sets, method = "gsva", kcdf = "Gaussian", verbose = FALSE, parallel.sz = 1))
gsva_samples <- data.frame(
  sample_id = sample_ids,
  OSARS_frozen = score_aligned$OSARS_frozen,
  OSARS_frozen_group = as.character(group),
  stage = stage,
  t(gsva_hallmark),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
write.csv(gsva_samples, file.path(analysis_dir, "data", "processed", "TCGA_Hallmark_GSVA_scores.csv"), row.names = FALSE)

run_gsva_summary <- function(score_matrix, group_factor, stage_vector) {
  rows <- lapply(rownames(score_matrix), function(pathway) {
    y <- as.numeric(score_matrix[pathway, ])
    d_primary <- data.frame(y = y, group = group_factor)
    fit_primary <- lm(y ~ group, data = d_primary)
    p_primary <- coef(summary(fit_primary))["groupHigh", "Pr(>|t|)"]
    beta_primary <- coef(fit_primary)["groupHigh"]
    keep <- !is.na(stage_vector)
    d_stage <- data.frame(y = y[keep], group = droplevels(group_factor[keep]), stage = factor(stage_vector[keep]))
    fit_stage <- lm(y ~ group + stage, data = d_stage)
    p_stage <- coef(summary(fit_stage))["groupHigh", "Pr(>|t|)"]
    beta_stage <- coef(fit_stage)["groupHigh"]
    data.frame(
      pathway = pathway,
      mean_low = mean(y[group_factor == "Low"]),
      mean_high = mean(y[group_factor == "High"]),
      median_low = median(y[group_factor == "Low"]),
      median_high = median(y[group_factor == "High"]),
      beta_high_vs_low = beta_primary,
      PValue = p_primary,
      beta_high_vs_low_stage_adjusted = beta_stage,
      PValue_stage_adjusted = p_stage,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$FDR <- p.adjust(out$PValue, method = "BH")
  out$FDR_stage_adjusted <- p.adjust(out$PValue_stage_adjusted, method = "BH")
  out[order(out$FDR, out$PValue), , drop = FALSE]
}

gsva_summary <- run_gsva_summary(gsva_hallmark, group, stage)
write.csv(gsva_summary, file.path(analysis_dir, "tables", "TCGA_Hallmark_GSVA_group_summary.csv"), row.names = FALSE)

# Representative pathways are selected mechanically only for a descriptive
# panel: the two smallest FDRs in each CAMERA direction, within Hallmark.
select_representative <- function(x, direction, n = 2L) {
  y <- x[x$Direction == direction, , drop = FALSE]
  y <- y[order(y$FDR, y$PValue, y$pathway), , drop = FALSE]
  head(y, n)
}
representative <- rbind(
  select_representative(camera_primary_h, "Up"),
  select_representative(camera_primary_h, "Down")
)
representative$selection_rule <- "Two smallest primary unadjusted Hallmark CAMERA FDRs within direction; descriptive GSVA display only"
write.csv(representative, file.path(analysis_dir, "tables", "representative_hallmark_pathways.csv"), row.names = FALSE)

db_versions <- unique(c(as.character(hallmark_raw$db_version), as.character(reactome_raw$db_version)))
db_versions <- db_versions[!is.na(db_versions) & nzchar(db_versions)]
metadata <- data.frame(
  field = c(
    "cohort", "expression_file", "score_file", "n_expression_samples", "n_OSARS_low", "n_OSARS_high",
    "n_stage_adjusted", "n_hallmark_tested", "n_reactome_tested", "pathway_test", "multiple_testing",
    "gene_set_source", "msigdbr_version", "MSigDB_database_version", "limma_version", "GSVA_version"
  ),
  value = c(
    "TCGA-LIHC", basename(expr_path), basename(score_path), ncol(expr), sum(group == "Low"), sum(group == "High"),
    sum(stage_keep), length(hallmark$sets), length(reactome$sets), "limma CAMERA competitive gene-set test", "Benjamini-Hochberg within each library and analysis",
    "msigdbr", as.character(packageVersion("msigdbr")), paste(db_versions, collapse = ";"), as.character(packageVersion("limma")), as.character(packageVersion("GSVA"))
  ),
  stringsAsFactors = FALSE
)
write.csv(metadata, file.path(analysis_dir, "references", "analysis_metadata.csv"), row.names = FALSE)

capture.output(sessionInfo(), file = file.path(analysis_dir, "references", "R_sessionInfo.txt"))

message("Completed TCGA OSARS-high bulk biology enrichment.")
