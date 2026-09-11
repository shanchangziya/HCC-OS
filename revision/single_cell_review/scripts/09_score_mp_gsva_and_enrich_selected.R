#!/usr/bin/env Rscript
# Panel D analysis for Reviewer 2 Major Comment 1.
#
# Selection and interpretation are intentionally separated:
#   1) score all eight validation GeneNMF metaprograms by GSVA in independent
#      patient-by-OS-quintile pseudobulks and identify the program most
#      increased in OS-high cells;
#   2) only after that selection, annotate the selected program by ORA.
#
# The primary selection removes the 49 frozen Hallmark ROS genes from every
# metaprogram before GSVA.  Full metaprograms are retained as a sensitivity
# analysis.  This avoids selecting a program merely because it directly
# contains genes used to define the OS score.

suppressPackageStartupMessages({
  library(Matrix)
  library(GSVA)
  library(msigdbr)
  library(BiocParallel)
  library(jsonlite)
})

set.seed(123L)

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
if (length(file_arg) != 1L) stop("Run this file with Rscript.")
script_path <- normalizePath(file_arg)
out_root <- normalizePath(file.path(dirname(script_path), ".."))
revision_root <- normalizePath(file.path(out_root, ".."))
input_dir <- file.path(out_root, "data", "genenmf_input")
source_dir <- file.path(revision_root, "single_cell", "data", "processed")
reference_dir <- file.path(revision_root, "model_validation", "references")

matrix_path <- file.path(input_dir, "GSE149614_primary_malignant_2000genes_counts.mtx.gz")
genes_path <- file.path(input_dir, "GSE149614_primary_malignant_2000genes_genes.csv")
cells_path <- file.path(input_dir, "GSE149614_primary_malignant_2000genes_cells.csv")
scores_path <- file.path(source_dir, "validation_primary_scores.csv")
mp_path <- file.path(out_root, "GSE149614_GeneNMF_selected_metaprogram_genes.csv")
ros_path <- file.path(reference_dir, "HALLMARK_ROS_msigdb_v7.0_genes.txt")

required <- c(matrix_path, genes_path, cells_path, scores_path, mp_path, ros_path)
if (any(!file.exists(required))) {
  stop("Missing required input(s): ", paste(required[!file.exists(required)], collapse = ", "))
}

# -------------------------------------------------------------------------
# Frozen input and alignment audit
# -------------------------------------------------------------------------

genes <- read.csv(genes_path, stringsAsFactors = FALSE)$gene
cells <- read.csv(cells_path, stringsAsFactors = FALSE, check.names = FALSE)
scores <- read.csv(scores_path, stringsAsFactors = FALSE, check.names = FALSE)
mp_table <- read.csv(mp_path, stringsAsFactors = FALSE, check.names = FALSE)
ros_genes <- unique(toupper(trimws(readLines(ros_path))))
ros_genes <- ros_genes[nzchar(ros_genes)]

stopifnot(length(genes) == 2000L, length(unique(genes)) == 2000L)
stopifnot(nrow(cells) == 13691L, nrow(scores) == 13691L)
stopifnot(all(c("cell", "patient", "nCount_RNA", "OS_composite") %in% names(scores)))
stopifnot(all(c("metaprogram", "gene", "rank") %in% names(mp_table)))

idx <- match(cells$cell, scores$cell)
if (anyNA(idx) || anyDuplicated(scores$cell) || anyDuplicated(cells$cell)) {
  stop("Cell identifiers are incomplete or non-unique during alignment.")
}
scores <- scores[idx, , drop = FALSE]
if (!identical(as.character(cells$cell), as.character(scores$cell))) {
  stop("Cell-order audit failed.")
}
if (!all(as.character(cells$patient) == as.character(scores$patient))) {
  stop("Patient labels disagree between frozen inputs.")
}

con <- gzfile(matrix_path, open = "rt")
counts <- Matrix::readMM(con)
close(con)
counts <- as(counts, "CsparseMatrix")
if (!identical(dim(counts), c(2000L, 13691L))) stop("Unexpected count-matrix dimensions.")
rownames(counts) <- genes
colnames(counts) <- cells$cell

# -------------------------------------------------------------------------
# Patient-by-OS-quintile pseudobulks
# -------------------------------------------------------------------------

# Equal-frequency OS quintiles are assigned independently within each patient,
# so patient-specific OS-score offsets cannot determine group membership.
os_quintile <- integer(nrow(scores))
for (patient_id in sort(unique(scores$patient))) {
  take <- which(scores$patient == patient_id)
  n <- length(take)
  ord <- order(scores$OS_composite[take], scores$cell[take])
  ranks <- integer(n)
  ranks[ord] <- seq_len(n)
  os_quintile[take] <- floor((ranks - 1L) * 5L / n) + 1L
}
scores$OS_quintile <- os_quintile
scores$group_id <- paste0(scores$patient, "_Q", scores$OS_quintile)

group_levels <- unlist(lapply(sort(unique(scores$patient)), function(x) paste0(x, "_Q", 1:5)))
if (!setequal(unique(scores$group_id), group_levels)) stop("Incomplete patient-by-quintile grid.")
group_index <- match(scores$group_id, group_levels)
membership <- sparseMatrix(
  i = seq_len(nrow(scores)), j = group_index, x = 1,
  dims = c(nrow(scores), length(group_levels))
)

pb_counts <- as.matrix(counts %*% membership)
colnames(pb_counts) <- group_levels
rownames(pb_counts) <- genes

group_sum <- function(x) {
  summed <- rowsum(matrix(x, ncol = 1L), group = group_index, reorder = TRUE)
  if (!identical(rownames(summed), as.character(seq_along(group_levels)))) {
    stop("Pseudobulk summary order does not match the fixed group order.")
  }
  as.numeric(summed[, 1])
}
group_meta <- data.frame(
  group_id = group_levels,
  patient = sub("_Q[1-5]$", "", group_levels),
  OS_quintile = as.integer(sub("^.*_Q", "", group_levels)),
  n_cells = as.integer(tabulate(group_index, nbins = length(group_levels))),
  library_size = group_sum(scores$nCount_RNA),
  mean_OS_composite = group_sum(scores$OS_composite) /
    as.integer(tabulate(group_index, nbins = length(group_levels))),
  stringsAsFactors = FALSE
)
if (any(group_meta$n_cells < 1L) || any(group_meta$library_size <= 0)) {
  stop("Empty or zero-library pseudobulk detected.")
}

# The denominator is the full-transcriptome frozen nCount_RNA sum, not merely
# the counts in the 2,000-gene GeneNMF universe.
log_cpm <- log2(sweep(pb_counts, 2L, group_meta$library_size, "/") * 1e6 + 1)

# -------------------------------------------------------------------------
# Custom GSVA scoring of the eight GeneNMF metaprograms
# -------------------------------------------------------------------------

mp_table$gene <- toupper(trimws(mp_table$gene))
mp_table <- mp_table[order(as.integer(sub("MP", "", mp_table$metaprogram)), mp_table$rank), ]
mp_sets_full <- lapply(split(mp_table$gene, mp_table$metaprogram), unique)
mp_order <- paste0("MP", 1:8)
mp_sets_full <- mp_sets_full[mp_order]
mp_sets_full <- lapply(mp_sets_full, intersect, y = toupper(genes))
mp_sets_no_ros <- lapply(mp_sets_full, setdiff, y = ros_genes)

if (any(vapply(mp_sets_full, length, integer(1)) < 10L) ||
    any(vapply(mp_sets_no_ros, length, integer(1)) < 10L)) {
  stop("A custom metaprogram gene set is too small for prespecified GSVA filtering.")
}

run_gsva <- function(gene_sets) {
  suppressWarnings(GSVA::gsva(
    expr = log_cpm,
    gset.idx.list = gene_sets,
    method = "gsva",
    kcdf = "Gaussian",
    abs.ranking = FALSE,
    min.sz = 10,
    max.sz = 500,
    mx.diff = TRUE,
    parallel.sz = 1,
    verbose = FALSE,
    BPPARAM = BiocParallel::SerialParam(progressbar = FALSE)
  ))
}

gsva_full <- run_gsva(mp_sets_full)
gsva_no_ros <- run_gsva(mp_sets_no_ros)
if (!identical(dim(gsva_full), c(8L, 50L)) || !identical(dim(gsva_no_ros), c(8L, 50L))) {
  stop("Unexpected GSVA result dimensions.")
}

long_scores <- function(mat, analysis) {
  pieces <- lapply(mp_order, function(mp) {
    raw <- as.numeric(mat[mp, group_levels])
    z <- as.numeric(scale(raw))
    data.frame(
      group_meta,
      metaprogram = mp,
      gene_set_variant = analysis,
      GSVA_score = raw,
      GSVA_z = z,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, pieces)
}
gsva_long <- rbind(
  long_scores(gsva_no_ros, "OS-gene-excluded"),
  long_scores(gsva_full, "full")
)
gsva_long$metaprogram <- factor(gsva_long$metaprogram, levels = mp_order)
gsva_long <- gsva_long[order(gsva_long$gene_set_variant, gsva_long$metaprogram,
                             gsva_long$patient, gsva_long$OS_quintile), ]
gsva_long$metaprogram <- as.character(gsva_long$metaprogram)
write.csv(gsva_long,
          file.path(out_root, "GSE149614_MP_GSVA_pseudobulk_scores.csv"),
          row.names = FALSE, quote = TRUE)

# Patient is the inferential unit.  For each MP, compare the within-patient
# highest and lowest OS quintiles after standardizing that MP's GSVA scores
# across the 50 independent pseudobulks.
association_one <- function(dat) {
  q1 <- dat[dat$OS_quintile == 1L, c("patient", "GSVA_z")]
  q5 <- dat[dat$OS_quintile == 5L, c("patient", "GSVA_z")]
  names(q1)[2] <- "q1"
  names(q5)[2] <- "q5"
  paired <- merge(q1, q5, by = "patient", sort = TRUE)
  delta <- paired$q5 - paired$q1
  if (length(delta) != 10L) stop("Expected ten patient-level Q5-Q1 contrasts.")
  tst <- t.test(delta, mu = 0, alternative = "two.sided")

  patient_rho <- vapply(split(dat, dat$patient), function(x) {
    cor(x$mean_OS_composite, x$GSVA_z, method = "spearman")
  }, numeric(1))
  rho_test <- t.test(patient_rho, mu = 0, alternative = "two.sided")

  data.frame(
    n_patients = length(delta),
    mean_delta_Q5_Q1 = mean(delta),
    lower_95CI = unname(tst$conf.int[1]),
    upper_95CI = unname(tst$conf.int[2]),
    t_statistic = unname(tst$statistic),
    P_paired_t = tst$p.value,
    median_patient_delta = median(delta),
    positive_patients = sum(delta > 0),
    mean_patient_spearman = mean(patient_rho),
    spearman_lower_95CI = unname(rho_test$conf.int[1]),
    spearman_upper_95CI = unname(rho_test$conf.int[2]),
    P_mean_spearman = rho_test$p.value,
    median_patient_spearman = median(patient_rho),
    positive_patient_spearman = sum(patient_rho > 0),
    patient_deltas = paste(sprintf("%s:%.6f", paired$patient, delta), collapse = ";"),
    patient_spearman = paste(sprintf("%s:%.6f", names(patient_rho), patient_rho), collapse = ";"),
    stringsAsFactors = FALSE
  )
}

assoc <- do.call(rbind, lapply(unique(gsva_long$gene_set_variant), function(variant) {
  do.call(rbind, lapply(mp_order, function(mp) {
    dat <- gsva_long[gsva_long$gene_set_variant == variant &
                       gsva_long$metaprogram == mp, ]
    cbind(
      data.frame(gene_set_variant = variant, metaprogram = mp, stringsAsFactors = FALSE),
      association_one(dat)
    )
  }))
}))
assoc$BH_FDR_Q5_Q1 <- ave(assoc$P_paired_t, assoc$gene_set_variant,
                          FUN = function(x) p.adjust(x, method = "BH"))
assoc$BH_FDR_mean_spearman <- ave(assoc$P_mean_spearman, assoc$gene_set_variant,
                                  FUN = function(x) p.adjust(x, method = "BH"))
assoc$rank_by_delta <- ave(-assoc$mean_delta_Q5_Q1, assoc$gene_set_variant,
                           FUN = function(x) rank(x, ties.method = "first"))
assoc$rank_by_delta <- as.integer(assoc$rank_by_delta)
assoc$rank_by_spearman <- ave(-assoc$mean_patient_spearman, assoc$gene_set_variant,
                              FUN = function(x) rank(x, ties.method = "first"))
assoc$rank_by_spearman <- as.integer(assoc$rank_by_spearman)

primary <- assoc[assoc$gene_set_variant == "OS-gene-excluded", ]
primary <- primary[order(-primary$mean_patient_spearman, primary$metaprogram), ]
selected_mp <- primary$metaprogram[[1]]
assoc$selected <- assoc$metaprogram == selected_mp
assoc <- assoc[order(assoc$gene_set_variant, assoc$rank_by_spearman), ]
write.csv(assoc,
          file.path(out_root, "GSE149614_MP_OS_association.csv"),
          row.names = FALSE, quote = TRUE)

full_top <- assoc[assoc$gene_set_variant == "full", ]
full_top <- full_top$metaprogram[which.max(full_top$mean_patient_spearman)]

overlap_tab <- data.frame(
  metaprogram = mp_order,
  full_gene_n = vapply(mp_sets_full, length, integer(1)),
  OS_gene_overlap_n = vapply(mp_sets_full, function(x) length(intersect(x, ros_genes)), integer(1)),
  OS_gene_overlap = vapply(mp_sets_full, function(x) paste(sort(intersect(x, ros_genes)), collapse = ";"), character(1)),
  OS_excluded_gene_n = vapply(mp_sets_no_ros, length, integer(1)),
  stringsAsFactors = FALSE
)
write.csv(overlap_tab,
          file.path(out_root, "GSE149614_MP_OS_gene_overlap.csv"),
          row.names = FALSE, quote = TRUE)

selection <- data.frame(
  selected_metaprogram = selected_mp,
  primary_gene_set_variant = "OS-gene-excluded",
  primary_selection_metric = "largest mean within-patient Spearman correlation between OS composite and GSVA score across five OS strata",
  n_patients = 10L,
  primary_rank = 1L,
  primary_mean_patient_spearman = primary$mean_patient_spearman[[1]],
  primary_spearman_lower_95CI = primary$spearman_lower_95CI[[1]],
  primary_spearman_upper_95CI = primary$spearman_upper_95CI[[1]],
  primary_spearman_BH_FDR = primary$BH_FDR_mean_spearman[[1]],
  primary_mean_delta_Q5_Q1 = primary$mean_delta_Q5_Q1[[1]],
  primary_lower_95CI = primary$lower_95CI[[1]],
  primary_upper_95CI = primary$upper_95CI[[1]],
  primary_delta_BH_FDR = primary$BH_FDR_Q5_Q1[[1]],
  full_gene_set_top_metaprogram = full_top,
  rank_concordant_with_full_gene_set = selected_mp == full_top,
  stringsAsFactors = FALSE
)
write.csv(selection,
          file.path(out_root, "GSE149614_OS_related_MP_selection.csv"),
          row.names = FALSE, quote = TRUE)

# -------------------------------------------------------------------------
# Biological annotation only after the OS-related MP is locked
# -------------------------------------------------------------------------

universe <- unique(toupper(genes))
query <- unique(intersect(mp_sets_full[[selected_mp]], universe))
N <- length(universe)
n <- length(query)

hallmark <- msigdbr(species = "Homo sapiens", collection = "H")
reactome <- msigdbr(species = "Homo sapiens", collection = "C2",
                    subcollection = "CP:REACTOME")

make_sets <- function(x, resource) {
  x <- as.data.frame(x[, c("db_version", "gs_name", "gene_symbol")])
  x$gene_symbol <- toupper(trimws(x$gene_symbol))
  x <- x[x$gene_symbol %in% universe, , drop = FALSE]
  split_sets <- lapply(split(x$gene_symbol, x$gs_name), unique)
  versions <- tapply(x$db_version, x$gs_name, function(z) as.character(z[[1]]))
  keep <- vapply(split_sets, length, integer(1)) >= 5L &
    vapply(split_sets, length, integer(1)) <= 500L
  split_sets <- split_sets[keep]
  data.frame(
    resource = resource,
    msigdbr_db_version = unname(versions[names(split_sets)]),
    gene_set = names(split_sets),
    size_in_background = vapply(split_sets, length, integer(1)),
    genes = I(unname(split_sets)),
    stringsAsFactors = FALSE
  )
}

sets <- rbind(make_sets(hallmark, "Hallmark"), make_sets(reactome, "Reactome"))
ora_rows <- lapply(seq_len(nrow(sets)), function(i) {
  pathway_genes <- sets$genes[[i]]
  overlap <- sort(intersect(query, pathway_genes))
  k <- length(overlap)
  m <- sets$size_in_background[[i]]
  expected <- n * m / N
  data.frame(
    selected_metaprogram = selected_mp,
    resource = sets$resource[[i]],
    msigdbr_db_version = sets$msigdbr_db_version[[i]],
    gene_set = sets$gene_set[[i]],
    background_n = N,
    query_n = n,
    size_in_background = m,
    overlap_n = k,
    overlap_genes = paste(overlap, collapse = ";"),
    expected_overlap = expected,
    fold_enrichment = if (expected > 0) k / expected else NA_real_,
    P = if (k == 0L) 1 else phyper(k - 1L, m, N - m, n, lower.tail = FALSE),
    stringsAsFactors = FALSE
  )
})
enrichment <- do.call(rbind, ora_rows)
enrichment$BH_FDR_all_resources <- p.adjust(enrichment$P, method = "BH")
enrichment$BH_FDR_within_resource <- NA_real_
for (resource in unique(enrichment$resource)) {
  take <- enrichment$resource == resource
  enrichment$BH_FDR_within_resource[take] <- p.adjust(enrichment$P[take], method = "BH")
}
enrichment$minus_log10_FDR <- -log10(pmax(enrichment$BH_FDR_all_resources, 1e-300))
enrichment$test_method <- "one-sided hypergeometric ORA"
enrichment$background_definition <- "2,000 genes eligible for the frozen validation GeneNMF analysis"
enrichment <- enrichment[order(enrichment$BH_FDR_all_resources, enrichment$P,
                               -enrichment$fold_enrichment, enrichment$gene_set), ]

# Objective redundancy reduction for plotting.  Hallmark and Reactome are
# represented separately (at most two and four terms, respectively), and terms
# are greedily retained by global FDR only when the MP-overlap genes have
# Jaccard < 0.50 with every already retained term from that resource.  This
# prevents many database labels driven by the same small HSP/ribosomal subset
# from crowding out distinct biological functions.  No keyword filter is used.
candidate_idx <- which(enrichment$BH_FDR_all_resources < 0.05 & enrichment$overlap_n >= 2L)
select_nonredundant <- function(idx, limit) {
  retained_local <- integer(0)
  for (i in idx) {
    candidate_genes <- strsplit(enrichment$overlap_genes[[i]], ";", fixed = TRUE)[[1]]
    redundant <- FALSE
    if (length(retained_local)) {
      redundant <- any(vapply(retained_local, function(j) {
        retained_genes <- strsplit(enrichment$overlap_genes[[j]], ";", fixed = TRUE)[[1]]
        length(intersect(candidate_genes, retained_genes)) /
          length(union(candidate_genes, retained_genes)) >= 0.50
      }, logical(1)))
    }
    if (!redundant) retained_local <- c(retained_local, i)
    if (length(retained_local) == limit) break
  }
  retained_local
}
hallmark_idx <- candidate_idx[enrichment$resource[candidate_idx] == "Hallmark"]
reactome_idx <- candidate_idx[enrichment$resource[candidate_idx] == "Reactome"]
retained <- c(select_nonredundant(hallmark_idx, 2L),
              select_nonredundant(reactome_idx, 4L))
enrichment$selected_for_panel_D <- seq_len(nrow(enrichment)) %in% retained
enrichment$redundancy_rule <- "within resource, greedy by global BH FDR; Jaccard of MP-overlap genes < 0.50; maximum 2 Hallmark plus 4 Reactome terms"
write.csv(enrichment,
          file.path(out_root, "GSE149614_selected_MP_pathway_enrichment.csv"),
          row.names = FALSE, quote = TRUE)

panel_pathways <- enrichment[enrichment$selected_for_panel_D, , drop = FALSE]
if (nrow(panel_pathways) < 3L) {
  stop("Fewer than three non-redundant FDR-significant pathways are available for Panel D.")
}
panel_pathways$display_order <- seq_len(nrow(panel_pathways))
write.csv(panel_pathways,
          file.path(out_root, "GSE149614_Panel_D_pathways_review_v2.csv"),
          row.names = FALSE, quote = TRUE)

manifest <- list(
  analysis = "GSE149614 custom GeneNMF-metaprogram GSVA and selected-program ORA",
  generated = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  input_cells = nrow(scores),
  patients = length(unique(scores$patient)),
  pseudobulks = nrow(group_meta),
  grouping = "five equal-frequency OS-composite strata within each patient",
  normalization = "log2(CPM + 1), CPM denominator from summed frozen full-transcriptome nCount_RNA",
  GSVA = list(
    package_version = as.character(packageVersion("GSVA")),
    method = "gsva",
    kcdf = "Gaussian",
    min_size = 10L,
    max_size = 500L,
    primary_gene_sets = "GeneNMF consensus genes after removal of the 49 frozen Hallmark ROS genes",
    sensitivity_gene_sets = "full GeneNMF consensus genes"
  ),
  selection = as.list(selection[1, ]),
  enrichment = list(
    method = "one-sided hypergeometric ORA",
    resources = c("MSigDB Hallmark", "MSigDB Reactome"),
    msigdbr_version = as.character(packageVersion("msigdbr")),
    database_version = unique(enrichment$msigdbr_db_version),
    background_n = N,
    query_n = n,
    correction = "BH across all tested Hallmark and Reactome terms"
  )
)
write_json(manifest,
           file.path(out_root, "GSE149614_Panel_D_analysis_manifest.json"),
           pretty = TRUE, auto_unbox = TRUE, digits = NA)

cat("PANEL_D_ANALYSIS_COMPLETE\n")
cat("Selected metaprogram:", selected_mp, "\n")
cat("Primary mean patient Spearman rho:", sprintf("%.4f", selection$primary_mean_patient_spearman), "\n")
cat("Primary 95% CI:", sprintf("%.4f to %.4f", selection$primary_spearman_lower_95CI,
                               selection$primary_spearman_upper_95CI), "\n")
cat("Primary BH FDR:", format(selection$primary_spearman_BH_FDR, digits = 4), "\n")
cat("Full-set top metaprogram:", full_top, "\n")
cat("Panel pathways:", nrow(panel_pathways), "\n")
