#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]
module_dir <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg)), ".."), mustWork = FALSE)
base_dir <- if (length(args) >= 1) args[1] else Sys.getenv("HCC_OS_SPATIAL_WORKDIR", module_dir)
pipeline_dir <- if (length(args) >= 2) args[2] else Sys.getenv("HCC_OS_DISCOVERY_ROOT")
result_dir <- file.path(base_dir, "results", "colocalization")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

focal_pairs <- data.frame(
  pair = c(
    "OShigh_TREM2",
    "OSlow_TREM2",
    "OShigh_SLC40A1",
    "OShigh_HSP"
  ),
  tumor = c("Tumor_OS-high", "Tumor_OS-low", "Tumor_OS-high", "Tumor_OS-high"),
  mac = c("Mac_TREM2_TAM", "Mac_TREM2_TAM", "Mac_SLC40A1_TAM", "Mac_HSP_TAM"),
  stringsAsFactors = FALSE
)

sanitize_names <- function(x) {
  x <- gsub("\\.", "-", x)
  x <- gsub("^Tumor_OS-high$", "Tumor_OS-high", x)
  x <- gsub("^Tumor_OS-low$", "Tumor_OS-low", x)
  x
}

read_cell2location_outputs <- function() {
  spatial_dir <- file.path(base_dir, "results", "cell2location_spatial")
  summary_file <- file.path(spatial_dir, "spatial_run_summary.csv")
  if (!file.exists(summary_file)) return(NULL)

  summary <- read.csv(summary_file, stringsAsFactors = FALSE)
  all <- lapply(seq_len(nrow(summary)), function(i) {
    abund <- read.csv(summary$output[i], row.names = 1, check.names = FALSE)
    abund$spot <- rownames(abund)
    abund$sample <- summary$sample[i]
    abund
  })
  bind_rows(all)
}

read_rctd_baseline <- function() {
  f <- file.path(pipeline_dir, "RCTD_deconvolution", "RCTD_cell_type_proportions.csv")
  if (!file.exists(f)) return(NULL)
  x <- read.csv(f, row.names = 1, check.names = FALSE)
  colnames(x) <- sanitize_names(colnames(x))
  x$spot <- rownames(x)
  x$sample <- sub("_.*$", "", x$spot)
  x
}

abund <- read_cell2location_outputs()
source <- "cell2location"
if (is.null(abund)) {
  abund <- read_rctd_baseline()
  source <- "RCTD_baseline"
}
if (is.null(abund)) stop("No cell2location outputs or RCTD baseline found")

colnames(abund) <- sanitize_names(colnames(abund))

needed <- unique(c(focal_pairs$tumor, focal_pairs$mac))
missing <- setdiff(needed, colnames(abund))
if (length(missing) > 0) {
  available <- paste(setdiff(colnames(abund), c("spot", "sample")), collapse = ", ")
  stop(
    "Missing abundance columns: ", paste(missing, collapse = ", "),
    "\nAvailable abundance columns: ", available,
    "\nThis usually means the baseline RCTD result used merged macrophages. ",
    "Run cell2location with macrophage subtypes retained before this core test."
  )
}

spot_cor <- bind_rows(lapply(split(abund, abund$sample), function(df) {
  bind_rows(lapply(seq_len(nrow(focal_pairs)), function(i) {
    p <- focal_pairs[i, ]
    ct1 <- df[[p$tumor]]
    ct2 <- df[[p$mac]]
    ok <- is.finite(ct1) & is.finite(ct2)
    pear <- suppressWarnings(cor.test(ct1[ok], ct2[ok], method = "pearson"))
    spear <- suppressWarnings(cor.test(ct1[ok], ct2[ok], method = "spearman"))
    data.frame(
      sample = unique(df$sample),
      pair = p$pair,
      tumor = p$tumor,
      mac = p$mac,
      n_spots = sum(ok),
      pearson_r = unname(pear$estimate),
      pearson_p = pear$p.value,
      spearman_rho = unname(spear$estimate),
      spearman_p = spear$p.value,
      stringsAsFactors = FALSE
    )
  }))
}))

spot_cor$pearson_fdr <- p.adjust(spot_cor$pearson_p, method = "BH")
spot_cor$spearman_fdr <- p.adjust(spot_cor$spearman_p, method = "BH")
write.csv(spot_cor, file.path(result_dir, paste0(source, "_spot_correlations.csv")), row.names = FALSE)

zscore <- function(x) as.numeric(scale(x))
enrich <- bind_rows(lapply(split(abund, abund$sample), function(df) {
  bind_rows(lapply(seq_len(nrow(focal_pairs)), function(i) {
    p <- focal_pairs[i, ]
    tumor_z <- zscore(df[[p$tumor]])
    mac_z <- zscore(df[[p$mac]])
    high <- tumor_z >= stats::quantile(tumor_z, 0.75, na.rm = TRUE)
    low <- tumor_z <= stats::quantile(tumor_z, 0.25, na.rm = TRUE)
    wt <- suppressWarnings(wilcox.test(mac_z[high], mac_z[low]))
    data.frame(
      sample = unique(df$sample),
      pair = p$pair,
      tumor = p$tumor,
      mac = p$mac,
      mac_mean_tumor_high_q = mean(mac_z[high], na.rm = TRUE),
      mac_mean_tumor_low_q = mean(mac_z[low], na.rm = TRUE),
      delta_z = mean(mac_z[high], na.rm = TRUE) - mean(mac_z[low], na.rm = TRUE),
      wilcox_p = wt$p.value,
      stringsAsFactors = FALSE
    )
  }))
}))
enrich$wilcox_fdr <- p.adjust(enrich$wilcox_p, method = "BH")
write.csv(enrich, file.path(result_dir, paste0(source, "_tumor_high_quartile_enrichment.csv")), row.names = FALSE)

p1 <- ggplot(spot_cor, aes(pair, spearman_rho, fill = pair == "OShigh_TREM2")) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey55") +
  geom_col(width = 0.7) +
  facet_wrap(~sample, scales = "free_y") +
  scale_fill_manual(values = c("TRUE" = "#C43B3B", "FALSE" = "#777777"), guide = "none") +
  labs(x = NULL, y = "Spearman rho", title = paste(source, "spot-level co-localization")) +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(result_dir, paste0(source, "_spot_correlation_barplot.pdf")), p1, width = 8.5, height = 4.8)

p2 <- ggplot(enrich, aes(pair, delta_z, fill = pair == "OShigh_TREM2")) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey55") +
  geom_col(width = 0.7) +
  facet_wrap(~sample, scales = "free_y") +
  scale_fill_manual(values = c("TRUE" = "#C43B3B", "FALSE" = "#777777"), guide = "none") +
  labs(x = NULL, y = "Delta Mac abundance z-score", title = paste(source, "TAM enrichment in high tumor-abundance spots")) +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(result_dir, paste0(source, "_quartile_enrichment_barplot.pdf")), p2, width = 8.5, height = 4.8)

message("Colocalization source: ", source)
message("Wrote: ", result_dir)
