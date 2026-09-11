#!/usr/bin/env Rscript
# SPDX-License-Identifier: GPL-3.0-only
# This file includes adapted GeneNMF v0.6.2 functions; see
# THIRD_PARTY_NOTICES.md and licenses/GeneNMF-GPL-3.0.txt.
# Refit the GSE149614 malignant-cell programs and select the number of
# metaprograms without using OS scores or pathway annotations.
#
# The analysis reproduces the GeneNMF v0.6.2 core implementation used by the
# frozen 2026-09-05 analysis.  The required functions below are transcribed
# from the official GPL-3 source at tag v0.6.2:
#   https://github.com/carmonalab/GeneNMF/tree/v0.6.2
# `lsa::cosine()` is evaluated by its algebraically identical normalized
# cross-product because `lsa` is unavailable in the frozen local R runtime.

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(cluster)
  library(RcppML)
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
result_dir <- file.path(out_root, "data", "genenmf_review")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------------
# GeneNMF v0.6.2 core functions
# -------------------------------------------------------------------------

normVector <- function(vector) {
  s <- sum(vector)
  if (s > 0) vector <- vector / s
  vector
}

weightCumul <- function(vector, weight.explained = 0.5) {
  x.sorted <- sort(vector, decreasing = TRUE)
  cs <- cumsum(x.sorted)
  norm.cs <- normVector(cs)
  norm.cs <- norm.cs / max(norm.cs)
  x.sorted[norm.cs < weight.explained]
}

wgtLoad <- function(matrix, w) {
  rownorm <- apply(matrix, 1, normVector)
  spec <- apply(rownorm, 2, max)
  spec.w <- spec^w
  matrix <- matrix * spec.w
  apply(matrix, 2, normVector)
}

weightedLoadings <- function(nmf.res, specificity.weight = 5) {
  lapply(nmf.res, function(model) wgtLoad(model$w, w = specificity.weight))
}

geneList2table <- function(gene.vectors) {
  vector.names <- names(gene.vectors)
  gene.vectors <- lapply(vector.names, function(n) {
    g <- gene.vectors[[n]]
    colnames(g) <- paste(n, seq_len(ncol(g)), sep = ".")
    g
  })
  Reduce(f = cbind, x = gene.vectors)
}

cosineSimilarity <- function(gene.table) {
  x <- as.matrix(gene.table)
  den <- sqrt(colSums(x^2))
  sim <- crossprod(x)
  sim <- sim / outer(den, den)
  sim[!is.finite(sim)] <- 0
  diag(sim) <- 1
  sim
}

getNMFgenes <- function(nmf.res, specificity.weight = 5,
                        weight.explained = 0.5, max.genes = 200) {
  if (!is.null(specificity.weight)) {
    nmf.res <- weightedLoadings(nmf.res, specificity.weight = specificity.weight)
  }
  nmf.genes <- lapply(nmf.res, function(model) {
    gene.pass <- apply(model, 2, function(x) {
      weightCumul(x, weight.explained = weight.explained)
    })
    m <- lapply(gene.pass, function(g) head(g, min(length(g), max.genes)))
    isna <- lapply(m, function(x) all(is.na(x)))
    m <- m[!as.numeric(isna)]
    names(m) <- seq_len(length(m))
    m
  })
  unlist(nmf.genes, recursive = FALSE)
}

get_metaprogram_consensus <- function(nmf.wgt, nMP = 10,
                                      min.confidence = 0.5,
                                      weight.explained = 0.5,
                                      max.genes = 200,
                                      cl_members = NULL) {
  nmf.genes.single <- getNMFgenes(
    nmf.res = nmf.wgt,
    specificity.weight = NULL,
    weight.explained = 0.8,
    max.genes = 1000
  )
  markers.consensus <- lapply(seq_len(nMP), function(cluster_id) {
    which.samples <- names(cl_members)[cl_members == cluster_id]
    gene.table <- geneList2table(nmf.wgt)[, which.samples, drop = FALSE]
    genes.avg <- apply(as.matrix(gene.table), 1, function(x) {
      this.mean <- mean(x)
      this.sd <- sd(x)
      x.out <- x[x > this.mean - 2 * this.sd & x < this.mean + 2 * this.sd]
      mean(x.out)
    })
    genes.avg <- sort(genes.avg, decreasing = TRUE)
    genes.pass <- weightCumul(genes.avg, weight.explained = weight.explained)
    this <- nmf.genes.single[which.samples]
    genes.only <- lapply(this, names)
    genes.sum <- sort(table(unlist(genes.only)), decreasing = TRUE)
    genes.confidence <- genes.sum / length(this)
    genes.confidence <- genes.confidence[genes.confidence > min.confidence]
    genes.pass <- genes.pass[names(genes.pass) %in% names(genes.confidence)]
    head(genes.pass, min(length(genes.pass), max.genes))
  })
  names(markers.consensus) <- paste0("MetaProgram", seq_len(nMP))
  markers.consensus
}

get_metaprogram_metrics <- function(J = NULL, Jdist = NULL,
                                    markers.consensus = NULL,
                                    cl_members = NULL) {
  nMP <- length(markers.consensus)
  all.samples <- unique(gsub("\\.k\\d+\\.\\d+", "", colnames(J)))
  sample.coverage <- lapply(seq_len(nMP), function(cluster_id) {
    which.samples <- names(cl_members)[cl_members == cluster_id]
    ss <- gsub("\\.k\\d+\\.\\d+", "", which.samples)
    ss <- factor(ss, levels = all.samples)
    ss.tab <- table(ss)
    sum(ss.tab > 0) / length(ss.tab)
  })
  names(sample.coverage) <- paste0("MetaProgram", seq_len(nMP))

  sil <- cluster::silhouette(cl_members, dist = Jdist)
  sil.widths <- summary(sil)$clus.avg.widths
  names(sil.widths) <- paste0("MetaProgram", seq_len(nMP))

  clusterSim <- rep(NA_real_, nMP)
  for (i in seq_len(nMP)) {
    selectMP <- which(cl_members == i)
    if (length(selectMP) > 1) {
      selectJ <- J[selectMP, selectMP, drop = FALSE]
      value <- round(mean(selectJ[upper.tri(selectJ)]), 3)
    } else {
      value <- 0
    }
    clusterSim[i] <- value
  }

  metaprograms.length <- unlist(lapply(markers.consensus, length))
  metaprograms.size <- as.character(table(cl_members))
  ans <- data.frame(
    sampleCoverage = unlist(sample.coverage),
    silhouette = sil.widths,
    meanSimilarity = clusterSim,
    numberGenes = metaprograms.length,
    numberPrograms = metaprograms.size
  )
  rownames(ans) <- paste0("MetaProgram", seq_len(nMP))
  ans
}

getMetaPrograms <- function(nmf.res, nMP = 10,
                            specificity.weight = 5,
                            weight.explained = 0.5,
                            max.genes = 200,
                            metric = c("cosine", "jaccard"),
                            hclust.method = "ward.D2",
                            min.confidence = 0.5,
                            remove.empty = TRUE) {
  metric <- metric[1]
  nmf.wgt <- weightedLoadings(nmf.res = nmf.res,
                              specificity.weight = specificity.weight)
  if (metric != "cosine") stop("This audit uses the prespecified cosine metric.")
  J <- cosineSimilarity(geneList2table(nmf.wgt))
  Jdist <- as.dist(1 - J)
  tree <- hclust(Jdist, method = hclust.method)
  cl_members <- cutree(tree, k = nMP)
  markers.consensus <- get_metaprogram_consensus(
    nmf.wgt = nmf.wgt,
    nMP = nMP,
    min.confidence = min.confidence,
    weight.explained = weight.explained,
    max.genes = max.genes,
    cl_members = cl_members
  )
  metaprograms.metrics <- get_metaprogram_metrics(
    J = J,
    Jdist = Jdist,
    markers.consensus = markers.consensus,
    cl_members = cl_members
  )
  if (remove.empty) {
    keep <- metaprograms.metrics$numberGenes > 0
    metaprograms.metrics <- metaprograms.metrics[keep, , drop = FALSE]
    markers.consensus <- markers.consensus[keep]
  }
  ord <- order(
    metaprograms.metrics$sampleCoverage,
    metaprograms.metrics$silhouette,
    decreasing = TRUE
  )
  old.names <- names(markers.consensus)[ord]
  new.names <- paste0("MP", seq_along(ord))
  markers.consensus <- markers.consensus[ord]
  names(markers.consensus) <- new.names
  metaprograms.metrics <- metaprograms.metrics[ord, , drop = FALSE]
  rownames(metaprograms.metrics) <- new.names
  map.index <- seq_along(old.names)
  names(map.index) <- as.numeric(gsub("MetaProgram", "", old.names))
  cl_members.new <- map.index[as.character(cl_members)]
  names(cl_members.new) <- names(cl_members)
  markers.consensus <- lapply(markers.consensus, function(m) m / sum(m))
  markers.genes <- lapply(markers.consensus, names)
  list(
    metaprograms.genes = markers.genes,
    metaprograms.genes.weights = markers.consensus,
    metaprograms.metrics = metaprograms.metrics,
    programs.similarity = J,
    programs.tree = tree,
    programs.clusters = cl_members.new
  )
}

getDataMatrix <- function(obj, assay = "RNA", slot = "data", hvg = NULL,
                          center = FALSE, scale = FALSE,
                          non_negative = TRUE) {
  mat <- Seurat::GetAssayData(obj, assay = assay, layer = slot)
  if (!is.null(hvg)) mat <- mat[hvg, ]
  mat <- Matrix::t(scale(Matrix::t(mat), center = center, scale = scale))
  if (scale) mat[is.na(mat)] <- 0
  if (non_negative) mat[mat < 0] <- 0
  mat
}

multiNMF <- function(obj.list, assay = "RNA", slot = "data", k = 5:6,
                     hvg = NULL, nfeatures = 2000, L1 = c(0, 0),
                     min.exp = 0.01, max.exp = 3.0,
                     center = FALSE, scale = FALSE,
                     min.cells.per.sample = 10,
                     hvg.blocklist = NULL, seed = 123) {
  set.seed(seed)
  nc <- unlist(lapply(obj.list, ncol))
  obj.list <- obj.list[nc > min.cells.per.sample]
  if (is.null(hvg)) stop("This frozen audit requires the saved 2,000-gene universe.")
  nmf.res <- lapply(obj.list, function(this) {
    mat <- getDataMatrix(
      obj = this,
      assay = assay,
      slot = slot,
      hvg = hvg,
      center = center,
      scale = scale
    )
    res.k <- lapply(k, function(k.this) {
      model <- RcppML::nmf(mat, k = k.this, L1 = L1,
                          verbose = FALSE, seed = seed)
      rownames(model$h) <- paste0("pattern", seq_len(nrow(model$h)))
      colnames(model$h) <- colnames(mat)
      rownames(model$w) <- rownames(mat)
      colnames(model$w) <- paste0("pattern", seq_len(ncol(model$w)))
      model
    })
    names(res.k) <- paste0("k", k)
    res.k
  })
  unlist(nmf.res, recursive = FALSE)
}

# -------------------------------------------------------------------------
# Frozen input and GeneNMF refit
# -------------------------------------------------------------------------

matrix_path <- file.path(
  input_dir, "GSE149614_primary_malignant_2000genes_counts.mtx.gz"
)
genes_path <- file.path(
  input_dir, "GSE149614_primary_malignant_2000genes_genes.csv"
)
cells_path <- file.path(
  input_dir, "GSE149614_primary_malignant_2000genes_cells.csv"
)

message("Reading frozen 2,000-gene count matrix...")
counts <- Matrix::readMM(gzfile(matrix_path))
genes <- read.csv(genes_path, stringsAsFactors = FALSE)$gene
cell_meta <- read.csv(cells_path, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(
  identical(dim(counts), c(2000L, 13691L)),
  length(genes) == 2000L,
  nrow(cell_meta) == 13691L,
  !anyDuplicated(genes),
  !anyDuplicated(cell_meta$cell),
  length(unique(cell_meta$patient)) == 10L
)
rownames(counts) <- genes
colnames(counts) <- cell_meta$cell
rownames(cell_meta) <- cell_meta$cell

message("Log-normalizing counts and splitting the ten patients...")
v <- CreateSeuratObject(counts = counts, assay = "RNA", meta.data = cell_meta)
v <- NormalizeData(
  v,
  normalization.method = "LogNormalize",
  scale.factor = 1e4,
  verbose = FALSE
)
objects <- SplitObject(v, split.by = "patient")
stopifnot(length(objects) == 10L)

message("Running GeneNMF: ten patients x k=4:9...")
programs <- multiNMF(
  objects,
  assay = "RNA",
  slot = "data",
  k = 4:9,
  hvg = genes,
  min.exp = 0.05,
  seed = 123
)
stopifnot(length(programs) == 60L)

# -------------------------------------------------------------------------
# Determine nMP using clustering quality only (no OS/pathway information)
# -------------------------------------------------------------------------

candidate_nmp <- 3:12
message("Auditing candidate metaprogram counts: ", paste(candidate_nmp, collapse = ", "))
candidate_objects <- lapply(candidate_nmp, function(nmp) {
  getMetaPrograms(
    programs,
    metric = "cosine",
    specificity.weight = 3,
    weight.explained = 0.3,
    nMP = nmp
  )
})
names(candidate_objects) <- as.character(candidate_nmp)

sweep_rows <- lapply(candidate_nmp, function(nmp) {
  mp <- candidate_objects[[as.character(nmp)]]
  raw_clusters <- cutree(mp$programs.tree, k = nmp)
  sil <- cluster::silhouette(raw_clusters, as.dist(1 - mp$programs.similarity))
  metrics <- mp$metaprograms.metrics
  data.frame(
    requested_nMP = nmp,
    retained_nMP = nrow(metrics),
    global_mean_silhouette = mean(sil[, "sil_width"]),
    median_MP_silhouette = median(metrics$silhouette),
    mean_MP_similarity = mean(metrics$meanSimilarity),
    mean_sample_coverage = mean(metrics$sampleCoverage),
    minimum_sample_coverage = min(metrics$sampleCoverage),
    n_MP_coverage_ge_0_5 = sum(metrics$sampleCoverage >= 0.5),
    n_MP_positive_silhouette = sum(metrics$silhouette > 0),
    stringsAsFactors = FALSE
  )
})
sweep <- do.call(rbind, sweep_rows)
sweep$rank_silhouette <- rank(-sweep$global_mean_silhouette, ties.method = "min")
sweep$selected <- FALSE

# Parsimonious plateau rule: first retain solutions reaching at least 95% of
# the best global mean silhouette, then prefer the solution with greatest mean
# patient coverage (and finally the smaller nMP).  This avoids choosing 11 MPs
# for a negligible silhouette gain while never using OS/pathway information.
silhouette_floor <- 0.95 * max(sweep$global_mean_silhouette)
sweep$within_95pct_silhouette_plateau <-
  sweep$global_mean_silhouette >= silhouette_floor
eligible <- sweep[sweep$global_mean_silhouette >= silhouette_floor, , drop = FALSE]
selection_order <- order(
  -eligible$mean_sample_coverage,
  -eligible$global_mean_silhouette,
  eligible$requested_nMP
)
selected_nmp <- eligible$requested_nMP[selection_order[1]]
sweep$selected[sweep$requested_nMP == selected_nmp] <- TRUE
selected <- candidate_objects[[as.character(selected_nmp)]]

write.csv(
  sweep,
  file.path(out_root, "GSE149614_GeneNMF_nMP_sweep.csv"),
  row.names = FALSE
)

selected_metrics <- selected$metaprograms.metrics
selected_metrics$metaprogram <- rownames(selected_metrics)
selected_metrics <- selected_metrics[, c(
  "metaprogram", "sampleCoverage", "silhouette", "meanSimilarity",
  "numberGenes", "numberPrograms"
)]
write.csv(
  selected_metrics,
  file.path(out_root, "GSE149614_GeneNMF_selected_metaprogram_metrics.csv"),
  row.names = FALSE
)

gene_rows <- do.call(rbind, lapply(names(selected$metaprograms.genes.weights), function(mp) {
  w <- selected$metaprograms.genes.weights[[mp]]
  data.frame(
    metaprogram = mp,
    rank = seq_along(w),
    gene = names(w),
    normalized_weight = as.numeric(w),
    stringsAsFactors = FALSE
  )
}))
write.csv(
  gene_rows,
  file.path(out_root, "GSE149614_GeneNMF_selected_metaprogram_genes.csv"),
  row.names = FALSE
)

membership <- data.frame(
  individual_program = names(selected$programs.clusters),
  metaprogram = paste0("MP", as.integer(selected$programs.clusters)),
  stringsAsFactors = FALSE
)
membership$patient <- sub("\\.k[0-9]+\\.[0-9]+$", "", membership$individual_program)
membership$k <- as.integer(sub("^.*\\.k([0-9]+)\\.[0-9]+$", "\\1", membership$individual_program))
membership$factor <- as.integer(sub("^.*\\.([0-9]+)$", "\\1", membership$individual_program))
write.csv(
  membership,
  file.path(out_root, "GSE149614_GeneNMF_individual_program_membership.csv"),
  row.names = FALSE
)

similarity <- selected$programs.similarity
write.csv(
  data.frame(individual_program = rownames(similarity), similarity,
             check.names = FALSE),
  gzfile(file.path(out_root, "GSE149614_GeneNMF_program_similarity.csv.gz")),
  row.names = FALSE
)

tree_table <- data.frame(
  order_index = seq_along(selected$programs.tree$order),
  individual_program = selected$programs.tree$labels[selected$programs.tree$order],
  stringsAsFactors = FALSE
)
write.csv(
  tree_table,
  file.path(out_root, "GSE149614_GeneNMF_heatmap_order.csv"),
  row.names = FALSE
)

# Reproduction audit against the historical nMP=4 export.
historical_genes_path <- file.path(
  revision_root, "single_cell", "data", "processed",
  "validation_de_novo_NMF_genes.csv"
)
historical_metrics_path <- file.path(
  revision_root, "single_cell", "data", "processed",
  "validation_de_novo_NMF_metrics.csv"
)
historical_genes <- read.csv(historical_genes_path, stringsAsFactors = FALSE)
historical_sets <- split(historical_genes$gene, historical_genes$program)
recomputed_four <- candidate_objects[["4"]]
recomputed_sets <- recomputed_four$metaprograms.genes
reproduction <- do.call(rbind, lapply(names(historical_sets), function(old_mp) {
  vals <- vapply(recomputed_sets, function(new_set) {
    length(intersect(historical_sets[[old_mp]], new_set)) /
      length(union(historical_sets[[old_mp]], new_set))
  }, numeric(1))
  best <- names(which.max(vals))
  data.frame(
    historical_metaprogram = old_mp,
    recomputed_metaprogram = best,
    gene_set_jaccard = unname(max(vals)),
    historical_n_genes = length(historical_sets[[old_mp]]),
    recomputed_n_genes = length(recomputed_sets[[best]]),
    stringsAsFactors = FALSE
  )
}))
write.csv(
  reproduction,
  file.path(out_root, "GSE149614_GeneNMF_nMP4_reproduction_audit.csv"),
  row.names = FALSE
)

historical_metrics <- read.csv(historical_metrics_path, stringsAsFactors = FALSE,
                               check.names = FALSE)
write.csv(
  historical_metrics,
  file.path(result_dir, "historical_nMP4_metrics_copy.csv"),
  row.names = FALSE
)

saveRDS(
  list(
    programs = programs,
    selected_metaprograms = selected,
    candidate_metaprograms = candidate_objects,
    selected_nMP = selected_nmp,
    candidate_sweep = sweep,
    implementation = "GeneNMF v0.6.2 official core functions",
    seed = 123L,
    k = 4:9,
    hvg = genes
  ),
  file.path(result_dir, "GSE149614_GeneNMF_refit.rds"),
  compress = FALSE
)

selection_manifest <- list(
  selected_nMP = selected_nmp,
  candidate_nMP = as.list(candidate_nmp),
  selection_rule = paste(
    "Among solutions reaching at least 95% of the maximum global mean",
    "silhouette across the 390 individual GeneNMF programs, select the",
    "solution with the greatest mean patient coverage; remaining ties are",
    "resolved by silhouette and then the smaller nMP. OS scores and pathway",
    "annotations were not used."
  ),
  silhouette_floor = silhouette_floor,
  geneNMF_implementation = "Official GeneNMF v0.6.2 core algorithm",
  geneNMF_source = "https://github.com/carmonalab/GeneNMF/tree/v0.6.2",
  cosine_note = paste(
    "lsa::cosine was replaced by the algebraically identical normalized",
    "cross-product because lsa was absent locally."
  ),
  input_cells = ncol(v),
  input_genes = nrow(v),
  patients = length(objects),
  individual_nmf_fits = length(programs),
  individual_programs = ncol(selected$programs.similarity),
  k = as.list(4:9),
  seed = 123L,
  specificity_weight = 3,
  weight_explained = 0.3,
  metric = "cosine",
  hclust_method = "ward.D2"
)
write_json(
  selection_manifest,
  file.path(out_root, "GSE149614_GeneNMF_selection_manifest.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)

writeLines(
  capture.output(sessionInfo()),
  file.path(result_dir, "GeneNMF_refit_sessionInfo.txt")
)

cat("GENENMF_REFIT_COMPLETE\n")
cat("Selected nMP:", selected_nmp, "\n")
print(sweep)
cat("nMP=4 reproduction audit:\n")
print(reproduction)
