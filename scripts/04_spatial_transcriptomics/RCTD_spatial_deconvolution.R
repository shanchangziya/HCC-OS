#!/usr/bin/env Rscript

# RCTD Spatial Transcriptomics Deconvolution
# Using single-cell reference with downsampling (30%)
# 3 Macrophage subtypes (HSP_TAM, SLC40A1_TAM, TREM2_TAM) retained
# Tumor OS-high/OS-low retained, other cell types merged to major categories

library(Seurat)
library(spacexr)
library(Matrix)
library(ggplot2)

set.seed(123)

cat("=== RCTD Spatial Deconvolution Analysis ===\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# Create output directory
output_dir <- "RCTD_deconvolution"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# Step 1: Load spatial transcriptomics data
# ============================================================
cat("Step 1: Loading spatial transcriptomics data...\n")
load("visium_ST_processed.RData")

# Identify the spatial object (could be named differently)
spatial_obj_names <- ls()[sapply(ls(), function(x) inherits(get(x), "Seurat"))]
cat("Available Seurat objects:", paste(spatial_obj_names, collapse=", "), "\n")

# Assume the first Seurat object is the spatial data
spatial_obj <- get(spatial_obj_names[1])
cat("Spatial object:", spatial_obj_names[1], "\n")
cat("Dimensions:", dim(spatial_obj), "\n\n")

# ============================================================
# Step 2: Load single-cell reference data
# ============================================================
cat("Step 2: Loading single-cell reference data...\n")
load("注释sceall.Rdata")

# Check object name
sc_obj_names <- ls()[sapply(ls(), function(x) inherits(get(x), "Seurat") & x != spatial_obj_names[1])]
cat("Single-cell object candidates:", paste(sc_obj_names, collapse=", "), "\n")

# Try common names
if(exists("sce.all")) {
  sc_ref <- sce.all
  cat("Using sce.all as reference\n")
} else if(exists("sceall")) {
  sc_ref <- sceall
  cat("Using sceall as reference\n")
} else {
  sc_ref <- get(sc_obj_names[1])
  cat("Using", sc_obj_names[1], "as reference\n")
}

cat("Reference dimensions:", dim(sc_ref), "\n")
cat("Cell types in reference:\n")
print(table(sc_ref$final_celltype))
cat("\n")

# ============================================================
# Step 3: Load macrophage subtype annotations
# ============================================================
cat("Step 3: Loading macrophage subtype annotations...\n")
load("08_Mac_annotated.Rdata")

# Find macrophage object
mac_obj_names <- ls()[sapply(ls(), function(x) inherits(get(x), "Seurat") &
                                    !x %in% c(spatial_obj_names, sc_obj_names))]
cat("Macrophage object candidates:", paste(mac_obj_names, collapse=", "), "\n")

if(exists("sce.mac")) {
  mac_obj <- sce.mac
  cat("Using sce.mac for macrophage annotations\n")
} else if(exists("mac")) {
  mac_obj <- mac
  cat("Using mac for macrophage annotations\n")
} else if(length(mac_obj_names) > 0) {
  mac_obj <- get(mac_obj_names[1])
  cat("Using", mac_obj_names[1], "for macrophage annotations\n")
} else {
  stop("Cannot find macrophage object in 08_Mac_annotated.Rdata")
}

cat("Macrophage object dimensions:", dim(mac_obj), "\n")

# Check which metadata column exists for cell type annotation
mac_celltype_col <- NULL
if("final_celltype" %in% colnames(mac_obj@meta.data)) {
  mac_celltype_col <- "final_celltype"
} else if("celltype" %in% colnames(mac_obj@meta.data)) {
  mac_celltype_col <- "celltype"
} else if("Mac_subtype" %in% colnames(mac_obj@meta.data)) {
  mac_celltype_col <- "Mac_subtype"
} else {
  # Find any column with "type" in the name
  type_cols <- grep("type", colnames(mac_obj@meta.data), value = TRUE, ignore.case = TRUE)
  if(length(type_cols) > 0) {
    mac_celltype_col <- type_cols[1]
  }
}

if(!is.null(mac_celltype_col)) {
  cat("Macrophage subtypes (using column:", mac_celltype_col, "):\n")
  print(table(mac_obj@meta.data[[mac_celltype_col]]))
} else {
  cat("Warning: Could not find cell type annotation column in macrophage object\n")
  cat("Available metadata columns:", paste(colnames(mac_obj@meta.data), collapse=", "), "\n")
}
cat("\n")

# ============================================================
# Step 4: Integrate macrophage annotations into main object
# ============================================================
cat("Step 4: Integrating macrophage annotations...\n")

# Get macrophage cell barcodes and their subtypes
mac_cells <- colnames(mac_obj)

# Use the detected cell type column
if(!is.null(mac_celltype_col)) {
  mac_annotations <- mac_obj@meta.data[[mac_celltype_col]]
  names(mac_annotations) <- mac_cells

  # Update annotations in main object
  sc_ref$celltype_for_rctd <- sc_ref$final_celltype

  # Replace macrophage annotations with subtypes
  for(cell in mac_cells) {
    if(cell %in% colnames(sc_ref)) {
      sc_ref$celltype_for_rctd[cell] <- mac_annotations[cell]
    }
  }
  cat("Integrated", length(mac_cells), "macrophage cells with subtype annotations\n")
} else {
  # If no macrophage subtype column found, just use the main annotation
  cat("Warning: Using main cell type annotations for macrophages (no subtype column found)\n")
  sc_ref$celltype_for_rctd <- sc_ref$final_celltype
}


cat("Updated cell type distribution:\n")
print(table(sc_ref$celltype_for_rctd))
cat("\n")

# ============================================================
# Step 5: Remove normal epithelial cells
# ============================================================
cat("Step 5: Removing Normal_Epithelial cells...\n")
cells_before <- ncol(sc_ref)

# Remove Normal_Epithelial
if("Normal_Epithelial" %in% sc_ref$celltype_for_rctd) {
  sc_ref <- subset(sc_ref, subset = celltype_for_rctd != "Normal_Epithelial")
  cat("Removed Normal_Epithelial cells\n")
}

cells_after <- ncol(sc_ref)
cat("Cells before:", cells_before, "| Cells after:", cells_after, "\n\n")

# ============================================================
# Step 6: Map cell types to final categories
# 3 Macrophage subtypes: HSP_TAM, SLC40A1_TAM, TREM2_TAM (keep)
# Tumor: Tumor_OS-high, Tumor_OS-low (keep)
# Other cells: merge to major categories
# ============================================================
cat("Step 6: Mapping cell types to final categories...\n")

# Define macrophage subtypes to keep (3 subtypes)
mac_subtypes <- c("HSP_TAM", "SLC40A1_TAM", "TREM2_TAM")

# Define tumor subtypes to keep (OS stratification)
tumor_subtypes <- c("Tumor_OS-high", "Tumor_OS-low")

# Define mapping for other cell types to major categories
celltype_mapping <- list(
  # T cells -> T
  "CD8_Trm" = "T",
  "CD8_Tem" = "T",
  "CD8_Tex" = "T",
  "CD8_Naive" = "T",
  "CD4_Cyto" = "T",
  "CD4_Treg" = "T",
  "CD4_Tfh" = "T",
  "CD4_Naive" = "T",
  "CD4_Memory" = "T",
  "T:CD8_Trm" = "T",
  "T:CD8_Tem" = "T",
  "T:CD8_Tex" = "T",
  "T:CD4_Cyto" = "T",
  "T:CD4_Treg" = "T",
  "T:CD4_Tfh" = "T",

  # NK cells -> NK
  "CD56dim_NK" = "NK",
  "CD56bright_NK" = "NK",
  "Activated_NK" = "NK",
  "Exhausted_NK" = "NK",
  "NK:CD56dim_NK" = "NK",
  "NK:CD56bright_NK" = "NK",
  "NK:Activated_NK" = "NK",
  "NK:Exhausted_NK" = "NK",

  # B cells -> B
  "Naive_B" = "B",
  "Memory_B" = "B",
  "Activated_B" = "B",
  "IgG_Plasma" = "B",
  "IgA_Plasma" = "B",
  "Plasma" = "B",

  # DC cells -> DC
  "cDC1" = "DC",
  "cDC2" = "DC",
  "pDC" = "DC",
  "moDC" = "DC",

  # CAF subtypes -> CAF
  "myCAF" = "CAF",
  "iCAF" = "CAF",
  "apCAF" = "CAF",
  "MMP_CAF" = "CAF",
  "VEGFA_CAF" = "CAF",

  # Endothelial subtypes -> Endothelial
  "VEC" = "Endothelial",
  "LSEC" = "Endothelial",
  "LEC" = "Endothelial",
  "Tip_Cell" = "Endothelial",
  "Tumor_Endo" = "Endothelial",

  # Neutrophil subtypes -> Neutrophil
  "Neutrophil" = "Neutrophil",
  "PMN_MDSC" = "Neutrophil",
  "M_MDSC" = "Neutrophil",
  "TAN" = "Neutrophil",
  "Inflam_Neutrophil" = "Neutrophil",

  # Other macrophage subtypes -> Macrophage (merge non-selected subtypes)
  "SPP1_TAM" = "Macrophage",
  "CXCL9_TAM" = "Macrophage",
  "C1QC_TAM" = "Macrophage",
  "FOLR2_TAM" = "Macrophage",
  "MARCO_TAM" = "Macrophage",
  "APOE_TAM" = "Macrophage"
)

# Apply mapping
for(i in 1:ncol(sc_ref)) {
  current_type <- sc_ref$celltype_for_rctd[i]

  # If it's one of the 3 macrophage subtypes, keep it
  if(current_type %in% mac_subtypes) {
    next
  }

  # If it's a tumor subtype (OS-high/OS-low), keep it
  if(current_type %in% tumor_subtypes) {
    next
  }

  # Otherwise, map to major category
  if(current_type %in% names(celltype_mapping)) {
    sc_ref$celltype_for_rctd[i] <- celltype_mapping[[current_type]]
  }
}

cat("Final cell type distribution:\n")
print(table(sc_ref$celltype_for_rctd))
cat("\n")

# ============================================================
# Step 7: Downsample 30% of cells per cell type
# ============================================================
cat("Step 7: Downsampling 30% of cells per cell type...\n")

cell_types <- unique(sc_ref$celltype_for_rctd)
sampled_cells <- c()

for(ct in cell_types) {
  ct_cells <- colnames(sc_ref)[sc_ref$celltype_for_rctd == ct]
  n_sample <- ceiling(length(ct_cells) * 0.3)
  sampled <- sample(ct_cells, n_sample)
  sampled_cells <- c(sampled_cells, sampled)
  cat(ct, ": ", length(ct_cells), " -> ", n_sample, " cells\n", sep="")
}

sc_ref_sampled <- subset(sc_ref, cells = sampled_cells)
cat("\nTotal cells after downsampling:", ncol(sc_ref_sampled), "\n")
cat("Cell type distribution after downsampling:\n")
print(table(sc_ref_sampled$celltype_for_rctd))
cat("\n")

# ============================================================
# Step 8: Prepare reference for RCTD
# ============================================================
cat("Step 8: Preparing reference for RCTD...\n")

# Extract counts matrix - handle both Seurat V3/V4 and V5
# Check if this is Seurat V5 (has @layers) or V3/V4 (has @counts)
if("layers" %in% slotNames(sc_ref_sampled@assays$RNA)) {
  # Seurat V5
  if("counts" %in% names(sc_ref_sampled@assays$RNA@layers)) {
    counts_matrix <- sc_ref_sampled@assays$RNA@layers$counts
    rownames(counts_matrix) <- rownames(sc_ref_sampled)
    colnames(counts_matrix) <- colnames(sc_ref_sampled)
  } else {
    counts_matrix <- GetAssayData(sc_ref_sampled, slot = "counts", assay = "RNA")
  }
} else {
  # Seurat V3/V4
  if(!is.null(sc_ref_sampled@assays$RNA@counts)) {
    counts_matrix <- sc_ref_sampled@assays$RNA@counts
  } else {
    counts_matrix <- GetAssayData(sc_ref_sampled, slot = "counts", assay = "RNA")
  }
}

# Create cell type factor
cell_types_factor <- factor(sc_ref_sampled$celltype_for_rctd)
names(cell_types_factor) <- colnames(sc_ref_sampled)

# Create Reference object
cat("Creating RCTD Reference object...\n")
reference <- Reference(counts_matrix, cell_types_factor)

cat("Reference created successfully\n")
cat("Number of cell types:", length(unique(cell_types_factor)), "\n")
cat("Number of cells:", ncol(counts_matrix), "\n")
cat("Number of genes:", nrow(counts_matrix), "\n\n")

# ============================================================
# Step 9: Prepare spatial data for RCTD
# ============================================================
cat("Step 9: Preparing spatial data for RCTD...\n")

# Extract spatial counts - handle both Seurat V3/V4 and V5
# First try Spatial assay, then RNA assay
spatial_assay <- NULL
if("Spatial" %in% names(spatial_obj@assays)) {
  spatial_assay <- "Spatial"
} else if("RNA" %in% names(spatial_obj@assays)) {
  spatial_assay <- "RNA"
} else {
  spatial_assay <- names(spatial_obj@assays)[1]
}

cat("Using assay:", spatial_assay, "\n")

# Check if this is Seurat V5 (has @layers) or V3/V4 (has @counts)
if("layers" %in% slotNames(spatial_obj@assays[[spatial_assay]])) {
  # Seurat V5
  if("counts" %in% names(spatial_obj@assays[[spatial_assay]]@layers)) {
    spatial_counts <- spatial_obj@assays[[spatial_assay]]@layers$counts
    rownames(spatial_counts) <- rownames(spatial_obj)
    colnames(spatial_counts) <- colnames(spatial_obj)
  } else {
    spatial_counts <- GetAssayData(spatial_obj, slot = "counts", assay = spatial_assay)
  }
} else {
  # Seurat V3/V4
  if(!is.null(spatial_obj@assays[[spatial_assay]]@counts)) {
    spatial_counts <- spatial_obj@assays[[spatial_assay]]@counts
  } else {
    spatial_counts <- GetAssayData(spatial_obj, slot = "counts", assay = spatial_assay)
  }
}

# Get spatial coordinates
if(!is.null(spatial_obj@images)) {
  image_name <- names(spatial_obj@images)[1]
  coords <- spatial_obj@images[[image_name]]@coordinates[, c("imagerow", "imagecol")]
  colnames(coords) <- c("x", "y")
} else {
  # Try to get from metadata
  coords <- spatial_obj@meta.data[, c("row", "col")]
  colnames(coords) <- c("x", "y")
}

cat("Spatial data dimensions:", dim(spatial_counts), "\n")
cat("Coordinates dimensions:", dim(coords), "\n\n")

# Create SpatialRNA object
cat("Creating RCTD SpatialRNA object...\n")
spatial_rna <- SpatialRNA(coords, spatial_counts, colSums(spatial_counts))

cat("SpatialRNA object created successfully\n")
cat("Number of spots:", ncol(spatial_counts), "\n")
cat("Number of genes:", nrow(spatial_counts), "\n\n")

# ============================================================
# Step 10: Run RCTD
# ============================================================
cat("Step 10: Running RCTD deconvolution...\n")
cat("This may take a while...\n\n")

# Create RCTD object
myRCTD <- create.RCTD(spatial_rna, reference, max_cores = 20)

# Run RCTD in full mode
myRCTD <- run.RCTD(myRCTD, doublet_mode = "full")

cat("RCTD deconvolution completed!\n\n")

# ============================================================
# Step 11: Extract and save results
# ============================================================
cat("Step 11: Extracting and saving results...\n")

# Get cell type proportions
results <- myRCTD@results

# Extract weights (cell type proportions)
weights <- normalize_weights(results$weights)

# Save weights matrix
write.csv(weights, file.path(output_dir, "RCTD_cell_type_proportions.csv"))
cat("Saved cell type proportions to:", file.path(output_dir, "RCTD_cell_type_proportions.csv"), "\n")

# Add proportions to spatial object
for(ct in colnames(weights)) {
  spatial_obj@meta.data[[paste0("RCTD_", ct)]] <- weights[, ct]
}

# Save updated spatial object
save(spatial_obj, file = file.path(output_dir, "spatial_with_RCTD.RData"))
cat("Saved spatial object with RCTD results to:", file.path(output_dir, "spatial_with_RCTD.RData"), "\n\n")

# Save RCTD object
save(myRCTD, file = file.path(output_dir, "RCTD_object.RData"))
cat("Saved RCTD object to:", file.path(output_dir, "RCTD_object.RData"), "\n\n")

# ============================================================
# Step 12: Generate visualizations
# ============================================================
cat("Step 12: Generating visualizations...\n")

# Plot cell type proportions
pdf(file.path(output_dir, "RCTD_cell_type_proportions.pdf"), width = 12, height = 10)

# Heatmap of proportions
library(pheatmap)
pheatmap(t(weights),
         main = "RCTD Cell Type Proportions",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_colnames = FALSE,
         color = colorRampPalette(c("white", "red"))(100))

# Boxplot of proportions per cell type
weights_long <- reshape2::melt(weights)
colnames(weights_long) <- c("Spot", "CellType", "Proportion")

p <- ggplot(weights_long, aes(x = CellType, y = Proportion, fill = CellType)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Distribution of Cell Type Proportions Across Spots",
       x = "Cell Type", y = "Proportion") +
  theme(legend.position = "none")
print(p)

dev.off()
cat("Saved visualizations to:", file.path(output_dir, "RCTD_cell_type_proportions.pdf"), "\n\n")

# Spatial plots for each cell type
pdf(file.path(output_dir, "RCTD_spatial_plots.pdf"), width = 15, height = 12)

for(ct in colnames(weights)) {
  tryCatch({
    p <- SpatialFeaturePlot(spatial_obj, features = paste0("RCTD_", ct)) +
      ggtitle(paste("RCTD:", ct))
    print(p)
  }, error = function(e) {
    cat("Warning: Could not create spatial plot for", ct, "\n")
  })
}

dev.off()
cat("Saved spatial plots to:", file.path(output_dir, "RCTD_spatial_plots.pdf"), "\n\n")

# ============================================================
# Step 13: Summary statistics
# ============================================================
cat("Step 13: Generating summary statistics...\n")

# Calculate summary statistics
summary_stats <- data.frame(
  CellType = colnames(weights),
  Mean_Proportion = colMeans(weights),
  Median_Proportion = apply(weights, 2, median),
  SD_Proportion = apply(weights, 2, sd),
  Min_Proportion = apply(weights, 2, min),
  Max_Proportion = apply(weights, 2, max)
)

write.csv(summary_stats, file.path(output_dir, "RCTD_summary_statistics.csv"), row.names = FALSE)
cat("Saved summary statistics to:", file.path(output_dir, "RCTD_summary_statistics.csv"), "\n\n")

print(summary_stats)

# ============================================================
# Completion
# ============================================================
cat("\n=== RCTD Analysis Completed Successfully ===\n")
cat("End time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("\nAll results saved to:", output_dir, "/\n")
cat("\nOutput files:\n")
cat("  - RCTD_cell_type_proportions.csv\n")
cat("  - RCTD_summary_statistics.csv\n")
cat("  - RCTD_cell_type_proportions.pdf\n")
cat("  - RCTD_spatial_plots.pdf\n")
cat("  - spatial_with_RCTD.RData\n")
cat("  - RCTD_object.RData\n")
cat("\nCell types analyzed:\n")
cat(paste("  -", colnames(weights), collapse = "\n"), "\n")
