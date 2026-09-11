#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(SingleCellExperiment)
  library(zellkonverter)
})

args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]
module_dir <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg)), ".."), mustWork = FALSE)
base_dir <- if (length(args) >= 1) args[1] else Sys.getenv("HCC_OS_SPATIAL_WORKDIR", module_dir)
pipeline_dir <- if (length(args) >= 2) args[2] else Sys.getenv("HCC_OS_DISCOVERY_ROOT")
if (!nzchar(pipeline_dir)) stop("Provide the discovery object root as argument 2 or HCC_OS_DISCOVERY_ROOT.")
input_dir <- file.path(base_dir, "inputs")
log_dir <- file.path(base_dir, "logs")
dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(123)
write_h5ad <- identical(Sys.getenv("WRITE_H5AD"), "1")
max_cells_per_type <- as.integer(Sys.getenv("MAX_CELLS_PER_TYPE", "3000"))
skip_sc_export <- identical(Sys.getenv("SKIP_SC_EXPORT"), "1")

log_file <- file.path(log_dir, "01_export_inputs.log")
sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)

cat("=== Export inputs for cell2location + MISTY ===\n")
cat("Start:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

get_counts <- function(obj, assay) {
  tryCatch(
    GetAssayData(obj, assay = assay, slot = "counts"),
    error = function(e) GetAssayData(obj, assay = assay, layer = "counts")
  )
}

write_sparse_bundle <- function(counts, meta, var, out_dir, prefix) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(counts, file.path(out_dir, paste0(prefix, "_counts.rds")))
  write.table(
    meta,
    gzfile(file.path(out_dir, paste0(prefix, "_obs.tsv.gz"))),
    sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA
  )
  write.table(
    var,
    gzfile(file.path(out_dir, paste0(prefix, "_var.tsv.gz"))),
    sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA
  )
}

make_sce <- function(counts, meta, var) {
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts),
    colData = S4Vectors::DataFrame(meta),
    rowData = S4Vectors::DataFrame(var)
  )
  sce
}

sc_out <- file.path(input_dir, "single_cell_reference")
if (skip_sc_export && file.exists(file.path(sc_out, "sc_reference_counts.rds"))) {
  cat("Skipping single-cell export because SKIP_SC_EXPORT=1 and reference bundle exists.\n")
} else {
  cat("Loading single-cell object...\n")
  load(file.path(pipeline_dir, "19_sce.all.integrated.Rdata"))
  sc <- sce.all.integrated
  rm(sce.all.integrated)
  DefaultAssay(sc) <- "RNA"

  if (!"final_celltype" %in% colnames(sc@meta.data)) {
    stop("final_celltype column not found in single-cell metadata")
  }

  sc <- subset(sc, subset = final_celltype != "Normal_Epithelial")
  sc$cell2location_label <- as.character(sc$final_celltype)

  min_cells <- 25
  cell_counts <- table(sc$cell2location_label)
  keep_types <- names(cell_counts[cell_counts >= min_cells])
  sc <- subset(sc, subset = cell2location_label %in% keep_types)

  if (is.finite(max_cells_per_type) && max_cells_per_type > 0) {
    cells_use <- unlist(lapply(split(colnames(sc), sc$cell2location_label), function(cells) {
      if (length(cells) <= max_cells_per_type) cells else sample(cells, max_cells_per_type)
    }), use.names = FALSE)
    sc <- subset(sc, cells = cells_use)
    cat("Applied per-celltype downsampling: max", max_cells_per_type, "cells/type\n")
  }

  cat("Single-cell dimensions after filtering:", paste(dim(sc), collapse = " x "), "\n")
  print(sort(table(sc$cell2location_label), decreasing = TRUE))

  sc_counts <- get_counts(sc, "RNA")
  sc_counts <- as(sc_counts, "dgCMatrix")
  sc_meta <- sc@meta.data
  sc_meta$barcode <- rownames(sc_meta)
  sc_var <- data.frame(gene_symbol = rownames(sc_counts), row.names = rownames(sc_counts))

  write_sparse_bundle(sc_counts, sc_meta, sc_var, sc_out, "sc_reference")

  if (write_h5ad) {
    cat("Writing single-cell h5ad...\n")
    sc_sce <- make_sce(sc_counts, sc_meta, sc_var)
    zellkonverter::writeH5AD(sc_sce, file.path(sc_out, "sc_reference.h5ad"), X_name = "counts")
    rm(sc_sce)
  } else {
    cat("Skipping single-cell h5ad. Set WRITE_H5AD=1 to enable.\n")
  }
  rm(sc, sc_counts)
  gc()
}

cat("\nLoading spatial object...\n")
load(file.path(pipeline_dir, "visium_ST_processed.RData"))
sp <- st_merged
rm(st_merged)

assay_names <- names(sp@assays)
image_names <- names(sp@images)
assay_use <- if ("Spatial" %in% assay_names) "Spatial" else DefaultAssay(sp)
DefaultAssay(sp) <- assay_use
sp_counts <- get_counts(sp, assay_use)
sp_counts <- as(sp_counts, "dgCMatrix")
sp_meta <- sp@meta.data
sp_meta$barcode <- rownames(sp_meta)

cat("Spatial dimensions:", paste(dim(sp), collapse = " x "), "\n")
cat("Spatial assay:", assay_use, "\n")
cat("Spatial images:", paste(image_names, collapse = ", "), "\n")
cat("Spatial samples:\n")
sample_col <- if ("sample" %in% colnames(sp_meta)) "sample" else "orig.ident"
print(table(as.character(sp_meta[[sample_col]])))

sp_var <- data.frame(gene_symbol = rownames(sp_counts), row.names = rownames(sp_counts))

coord_for_image <- function(obj, image_name) {
  coords <- tryCatch(
    GetTissueCoordinates(obj, image = image_name),
    error = function(e) obj@images[[image_name]]@coordinates
  )
  coords$barcode <- rownames(coords)
  coords$image <- image_name
  coords
}

all_coords <- do.call(rbind, lapply(image_names, function(img) coord_for_image(sp, img)))
write.table(
  all_coords,
  gzfile(file.path(input_dir, "spatial_all_coordinates.tsv.gz")),
  sep = "\t", quote = FALSE, row.names = FALSE
)

samples <- sort(unique(as.character(sp_meta[[sample_col]])))
manifest <- data.frame()

for (sample_id in samples) {
  cat("\nExporting spatial sample:", sample_id, "\n")
  cells <- rownames(sp_meta)[as.character(sp_meta[[sample_col]]) == sample_id]
  counts_i <- sp_counts[, cells, drop = FALSE]
  meta_i <- sp_meta[cells, , drop = FALSE]
  coords_i <- all_coords[all_coords$barcode %in% cells, , drop = FALSE]

  out_i <- file.path(input_dir, "spatial", sample_id)
  write_sparse_bundle(counts_i, meta_i, sp_var, out_i, "spatial")
  write.table(
    coords_i,
    gzfile(file.path(out_i, "spatial_coords.tsv.gz")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  if (write_h5ad) {
    sce_i <- make_sce(counts_i, meta_i, sp_var)
    xy_cols <- if (all(c("imagecol", "imagerow") %in% colnames(coords_i))) {
      c("imagecol", "imagerow")
    } else if (all(c("x", "y") %in% colnames(coords_i))) {
      c("x", "y")
    } else {
      setdiff(colnames(coords_i), c("barcode", "image", "cell", "tissue"))[1:2]
    }
    coord_mat <- as.matrix(coords_i[match(colnames(sce_i), coords_i$barcode), xy_cols, drop = FALSE])
    colnames(coord_mat) <- c("x", "y")
    SingleCellExperiment::reducedDim(sce_i, "spatial") <- coord_mat
    zellkonverter::writeH5AD(sce_i, file.path(out_i, "spatial.h5ad"), X_name = "counts")
  }

  manifest <- rbind(
    manifest,
    data.frame(
      sample = sample_id,
      n_spots = length(cells),
      n_genes = nrow(counts_i),
      input_dir = out_i,
      h5ad = if (write_h5ad) file.path(out_i, "spatial.h5ad") else NA_character_,
      stringsAsFactors = FALSE
    )
  )
}

write.csv(manifest, file.path(input_dir, "spatial_manifest.csv"), row.names = FALSE)

cat("\nDone. Manifest:\n")
print(manifest)
cat("End:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
