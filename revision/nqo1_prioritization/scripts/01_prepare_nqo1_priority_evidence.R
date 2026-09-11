#!/usr/bin/env Rscript

# Build auditable, plot-ready evidence for NQO1 post-model prioritization.
#
# Scope: This is not a re-fit of OSARS and not a new test of a single-cell
# state. It reconstructs the selection evidence from frozen project files:
# (1) membership of the frozen 128-gene OSARS feature list,
# (2) overlap with the original HALLMARK ROS collection,
# (3) descriptive discovery marker effect size/specificity, and
# (4) independent paired PDC000198 protein abundance.

script_arg <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", script_arg[grepl("^--file=", script_arg)][1])
script_path <- gsub("~\\+~", " ", script_path)
script_path <- normalizePath(script_path)
script_dir <- dirname(script_path)
analysis_dir <- normalizePath(file.path(script_dir, ".."))
revision_dir <- normalizePath(file.path(analysis_dir, ".."))

dir.create(file.path(analysis_dir, "data", "processed"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(analysis_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(analysis_dir, "references"), recursive = TRUE, showWarnings = FALSE)

osars_table_path <- file.path(revision_dir, "tables", "Table2_frozen_OSARS_128_genes.csv")
marker_path <- file.path(revision_dir, "model_validation", "data", "raw", "discovery_deg_os.csv")
importance_path <- file.path(revision_dir, "model_validation", "data", "raw", "frozen_variable_importance.csv")
protein_dir <- file.path(revision_dir, "protein_support", "data", "processed")
pairs_path <- file.path(protein_dir, "nqo1_verified_pairs.csv")
paired_statistics_path <- file.path(protein_dir, "paired_NQO1_statistics.csv")

stopifnot(
  file.exists(osars_table_path), file.exists(marker_path), file.exists(importance_path),
  file.exists(pairs_path), file.exists(paired_statistics_path)
)

osars <- read.csv(osars_table_path, check.names = FALSE, stringsAsFactors = FALSE)
markers <- read.csv(marker_path, check.names = FALSE, stringsAsFactors = FALSE)
importance <- read.csv(importance_path, check.names = FALSE, stringsAsFactors = FALSE)

required_osars <- c("gene", "in_Hallmark_ROS_v7_0", "in_frozen_discovery_MP4")
required_markers <- c("gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")
stopifnot(all(required_osars %in% names(osars)), all(required_markers %in% names(markers)), all(c("var", "rel.inf") %in% names(importance)))

osars$in_Hallmark_ROS_v7_0 <- tolower(as.character(osars$in_Hallmark_ROS_v7_0)) == "true"
osars$in_frozen_discovery_MP4 <- tolower(as.character(osars$in_frozen_discovery_MP4)) == "true"
markers$detect_gap <- markers$pct.1 - markers$pct.2
markers$marker_rank_avg_log2FC <- rank(-markers$avg_log2FC, ties.method = "min")
markers$marker_rank_detect_gap <- rank(-markers$detect_gap, ties.method = "min")
markers$marker_rank_adjusted_p <- rank(markers$p_val_adj, ties.method = "min")
importance$model_importance_rank <- rank(-importance$rel.inf, ties.method = "min")

# These are all the frozen OSARS genes that also belong to the frozen
# HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY source collection.
candidates <- osars[osars$in_Hallmark_ROS_v7_0, , drop = FALSE]
candidates <- merge(candidates, markers, by = "gene", all.x = TRUE, sort = FALSE)
candidates <- merge(candidates, importance, by.x = "gene", by.y = "var", all.x = TRUE, sort = FALSE)
if (!identical(sort(candidates$gene), sort(c("NQO1", "PRDX1", "TXNRD1")))) {
  stop("Unexpected OSARS-HALLMARK ROS overlap; do not continue without auditing the frozen source table.")
}
candidates$selected_for_follow_up <- candidates$gene == "NQO1"
candidates$selection_basis <- ifelse(
  candidates$selected_for_follow_up,
  "Largest descriptive marker effect and detection gap among the OSARS-ROS overlap",
  "Not selected for follow-up"
)
candidates <- candidates[order(-candidates$avg_log2FC, candidates$gene), , drop = FALSE]
write.csv(candidates, file.path(analysis_dir, "data", "processed", "OSARS_ROS_overlap_candidate_evidence.csv"), row.names = FALSE)

# Whole marker population is used only as a visual rank context. It retains no
# new inferential claim: the original marker screen is descriptive/post-cluster.
markers$OSARS_ROS_overlap_gene <- markers$gene %in% candidates$gene
markers$selected_for_follow_up <- markers$gene == "NQO1"
write.csv(markers, file.path(analysis_dir, "data", "processed", "discovery_marker_rank_context.csv"), row.names = FALSE)

pair_data <- read.csv(pairs_path, check.names = FALSE, stringsAsFactors = FALSE)
pair_stats <- read.csv(paired_statistics_path, check.names = FALSE, stringsAsFactors = FALSE)
pair_stats <- pair_stats[pair_stats$subset == "All verified pairs", , drop = FALSE]
stopifnot(nrow(pair_stats) == 1L, nrow(pair_data) == pair_stats$n_pairs[1])
write.csv(pair_data, file.path(analysis_dir, "data", "processed", "PDC000198_NQO1_verified_pairs.csv"), row.names = FALSE)
write.csv(pair_stats, file.path(analysis_dir, "data", "processed", "PDC000198_NQO1_paired_statistics.csv"), row.names = FALSE)

n_markers <- nrow(markers)
n_osars <- nrow(osars)
n_ros_osars <- nrow(candidates)
nq <- candidates[candidates$gene == "NQO1", , drop = FALSE]
summary_tbl <- data.frame(
  metric = c(
    "Frozen OSARS feature count", "Frozen OSARS features overlapping HALLMARK ROS", "Discovery positive marker count",
    "NQO1 average-log2FC rank among discovery markers", "NQO1 detection-gap rank among discovery markers",
    "NQO1 adjusted-P rank among discovery markers", "NQO1 GBM relative-importance rank among frozen OSARS features",
    "NQO1 membership in frozen discovery MP4", "Verified PDC000198 tumor-adjacent protein pairs"
  ),
  value = c(
    n_osars, n_ros_osars, n_markers,
    nq$marker_rank_avg_log2FC, nq$marker_rank_detect_gap, nq$marker_rank_adjusted_p,
    nq$model_importance_rank, nq$in_frozen_discovery_MP4, nrow(pair_data)
  ),
  stringsAsFactors = FALSE
)
write.csv(summary_tbl, file.path(analysis_dir, "tables", "NQO1_priority_evidence_summary.csv"), row.names = FALSE)

sources <- data.frame(
  source_role = c("Frozen OSARS feature table", "Discovery marker screen", "Frozen GBM variable importance", "PDC000198 verified protein pairs", "PDC000198 paired statistics"),
  source_file = c(
    "revision/tables/Table2_frozen_OSARS_128_genes.csv",
    "revision/model_validation/data/raw/discovery_deg_os.csv",
    "revision/model_validation/data/raw/frozen_variable_importance.csv",
    "revision/protein_support/data/processed/nqo1_verified_pairs.csv",
    "revision/protein_support/data/processed/paired_NQO1_statistics.csv"
  ),
  stringsAsFactors = FALSE
)
write.csv(sources, file.path(analysis_dir, "references", "source_files.csv"), row.names = FALSE)

message("Prepared NQO1 candidate-prioritization and protein evidence tables.")
