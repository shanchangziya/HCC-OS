#!/usr/bin/env Rscript
# Independent GSE149614 validation-program enrichment.
# This script only reads frozen, processed single-cell outputs and writes new
# source-data tables for Reviewer 2 Major Comment 1.  It does not refit NMF,
# alter cell annotations, or change the OS scoring definition.

suppressPackageStartupMessages(library(msigdbr))

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
if (length(file_arg) != 1L) stop("Run this file with Rscript.")
script_path <- normalizePath(file_arg)
out_root <- normalizePath(file.path(dirname(script_path), ".."))
revision_root <- normalizePath(file.path(out_root, ".."))
source_dir <- file.path(revision_root, "single_cell", "data", "processed")

nmf_genes <- read.csv(file.path(source_dir, "validation_de_novo_NMF_genes.csv"),
                      stringsAsFactors = FALSE, check.names = FALSE)
universe <- read.csv(file.path(source_dir, "NMF_eligible_gene_universe.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)$gene
similarity <- read.csv(file.path(source_dir, "cross_cohort_NMF_overlap.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)

universe <- unique(toupper(trimws(universe)))
nmf_genes$gene <- toupper(trimws(nmf_genes$gene))
stopifnot(length(universe) == 2000L,
          all(nmf_genes$gene %in% universe),
          all(c("discovery", "validation", "jaccard", "FDR") %in% names(similarity)))

# The frozen overlap table determines the validation program to foreground; this
# is not selected after enrichment.
mp4_matches <- similarity[similarity$discovery == "MP4", ]
mp4_matches <- mp4_matches[order(-mp4_matches$jaccard, mp4_matches$FDR), ]
key_validation_program <- as.character(mp4_matches$validation[[1]])

hallmark <- msigdbr(species = "Homo sapiens", collection = "H")
reactome <- msigdbr(species = "Homo sapiens", collection = "C2",
                    subcollection = "CP:REACTOME")

# Construct gene-set lists without relying on an online enrichment server.  The
# test universe is the 2,000 genes actually eligible for the saved validation
# NMF W matrices, rather than the whole genome.
build_sets <- function(x, resource) {
  x <- as.data.frame(x[, c("db_version", "gs_name", "gene_symbol")])
  x$gene_symbol <- toupper(trimws(x$gene_symbol))
  x <- x[x$gene_symbol %in% universe, , drop = FALSE]
  sets <- split(x$gene_symbol, x$gs_name)
  db <- tapply(x$db_version, x$gs_name, function(z) as.character(z[[1]]))
  data.frame(
    resource = resource,
    msigdbr_db_version = unname(db[names(sets)]),
    gene_set = names(sets),
    size_in_background = vapply(sets, function(z) length(unique(z)), integer(1)),
    stringsAsFactors = FALSE
  ) -> tab
  tab$genes <- unname(lapply(sets, unique)[tab$gene_set])
  tab
}

sets <- rbind(build_sets(hallmark, "Hallmark"),
              build_sets(reactome, "Reactome"))
# Exclude poorly powered and overly broad terms before testing.  This fixed rule
# is applied identically to every NMF program and is recorded in the output.
sets <- sets[sets$size_in_background >= 5L & sets$size_in_background <= 500L, , drop = FALSE]
N <- length(universe)

ora_one_program <- function(program, genes) {
  query <- unique(intersect(genes, universe))
  n <- length(query)
  rows <- lapply(seq_len(nrow(sets)), function(i) {
    set_genes <- sets$genes[[i]]
    k <- length(intersect(query, set_genes))
    m <- sets$size_in_background[[i]]
    data.frame(
      validation_program = program,
      resource = sets$resource[[i]],
      msigdbr_db_version = sets$msigdbr_db_version[[i]],
      gene_set = sets$gene_set[[i]],
      background_n = N,
      query_n = n,
      size_in_background = m,
      overlap_n = k,
      overlap_genes = paste(sort(intersect(query, set_genes)), collapse = ";"),
      P = if (k == 0L) 1 else phyper(k - 1L, m, N - m, n, lower.tail = FALSE),
      stringsAsFactors = FALSE
    )
  })
  ans <- do.call(rbind, rows)
  ans$FDR_within_program_all_resources <- p.adjust(ans$P, method = "BH")
  ans$FDR_within_program_resource <- NA_real_
  for (resource in unique(ans$resource)) {
    idx <- which(ans$resource == resource)
    ans$FDR_within_program_resource[idx] <- p.adjust(ans$P[idx], method = "BH")
  }
  ans
}

program_sets <- split(nmf_genes$gene, nmf_genes$program)
enrichment <- do.call(rbind, lapply(names(program_sets), function(program) {
  ora_one_program(program, program_sets[[program]])
}))
enrichment$test_method <- "one-sided hypergeometric ORA"
enrichment$background_definition <- "2,000 genes actually eligible for saved validation NMF W matrices"
enrichment$program_selected_for_panel_D <- enrichment$validation_program == key_validation_program
enrichment <- enrichment[order(enrichment$validation_program,
                               enrichment$FDR_within_program_all_resources,
                               enrichment$P, enrichment$gene_set), ]
write.csv(enrichment,
          file.path(out_root, "GSE149614_NMF_program_enrichment.csv"),
          row.names = FALSE, quote = TRUE)

# Supply a transparent, predeclared candidate table for the plotting script.
# The plotting script can select a non-redundant subset from the significant
# results, while this table preserves every tested pathway and its FDR.
panel_d_candidates <- enrichment[
  enrichment$validation_program == key_validation_program &
    enrichment$FDR_within_program_all_resources < 0.05 &
    enrichment$overlap_n >= 2L,
  ]
write.csv(panel_d_candidates,
          file.path(out_root, "GSE149614_key_MP_pathway_candidates.csv"),
          row.names = FALSE, quote = TRUE)

selection <- data.frame(
  discovery_program = "MP4",
  selected_validation_program = key_validation_program,
  selection_rule = "Highest raw measured-gene Jaccard among the four frozen discovery MP4-to-validation NMF comparisons",
  raw_measured_jaccard = mp4_matches$jaccard[[1]],
  eligible_universe_jaccard = mp4_matches$jaccard_eligible[[1]],
  overlap_genes = mp4_matches$overlap[[1]],
  Fisher_P_eligible_background = mp4_matches$P[[1]],
  BH_FDR_all_16_comparisons = mp4_matches$FDR[[1]],
  stringsAsFactors = FALSE
)
write.csv(selection, file.path(out_root, "GSE149614_key_MP_selection.csv"),
          row.names = FALSE, quote = TRUE)

cat("VALIDATION_NMF_ENRICHMENT_COMPLETE\n")
cat("Key validation program:", key_validation_program, "\n")
