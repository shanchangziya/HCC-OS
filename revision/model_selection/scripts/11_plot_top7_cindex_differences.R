#!/usr/bin/env Rscript

args_all <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args_all[grepl("^--file=", args_all)])
script_path <- normalizePath(script_arg[[1]], mustWork = TRUE)
analysis_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

suppressPackageStartupMessages(library(ggplot2))

input_file <- file.path(analysis_root, "data", "processed", "Top7_Cindex_differences_vs_OSARS.csv")
summary_file <- file.path(analysis_root, "data", "processed", "Top7_Cindex_statistical_summary.csv")
output_root <- file.path(analysis_root, "figures", "top7_comparison")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

comparison <- read.csv(input_file, check.names = FALSE, stringsAsFactors = FALSE)
statistical_summary <- read.csv(summary_file, check.names = FALSE, stringsAsFactors = FALSE)
if (nrow(comparison) != 7L) stop("Expected seven model configurations")
if (sum(comparison$is_OSARS) != 1L) stop("Expected one OSARS reference")

osars_configuration <- "StepCox[forward] + GBM"
comparison <- comparison[order(comparison$CV_rank, -comparison$mean_5fold_Cindex), , drop = FALSE]
comparison$display_label <- comparison$configuration
comparison$display_label[comparison$is_OSARS] <- paste0(osars_configuration, " (OSARS)")
comparison$display_label <- factor(comparison$display_label, levels = rev(comparison$display_label))

global_p <- statistical_summary$pvalue_or_min_adjusted_pvalue[
  grepl("Global comparison", statistical_summary$analysis, fixed = TRUE)
]
pairwise_significant <- statistical_summary$significant_comparisons_0_05[
  grepl("all 21 pairs", statistical_summary$analysis, fixed = TRUE)
]
if (length(global_p) != 1L || length(pairwise_significant) != 1L) {
  stop("Statistical summary is incomplete")
}

accent <- "#B64342"
neutral_ci <- "#8F979F"
neutral_point <- "#343A40"

others <- comparison[!comparison$is_OSARS, , drop = FALSE]
osars <- comparison[comparison$is_OSARS, , drop = FALSE]

base_theme <- theme_classic(base_family = "Arial", base_size = 7.0) +
  theme(
    axis.line = element_line(linewidth = 0.35, colour = "#202124"),
    axis.ticks = element_line(linewidth = 0.3, colour = "#202124"),
    axis.ticks.length = grid::unit(1.6, "pt"),
    axis.text.x = element_text(size = 6.4, colour = "#202124"),
    axis.text.y = element_text(size = 5.7, colour = "#202124", margin = margin(r = 2.5)),
    axis.title.x = element_text(size = 7.0, colour = "#202124", margin = margin(t = 4)),
    plot.title = element_text(size = 8.2, hjust = 0.5, face = "plain", margin = margin(b = 1.5)),
    plot.subtitle = element_text(size = 5.9, hjust = 0.5, colour = "#5F6368", margin = margin(b = 4)),
    plot.tag = element_text(
      family = "Times New Roman", face = "bold", size = 9.2,
      colour = "#111111", hjust = 0, vjust = 1
    ),
    plot.tag.position = c(0.012, 0.992),
    plot.margin = margin(t = 3, r = 4, b = 3, l = 3, unit = "pt")
  )

pB <- ggplot(comparison, aes(y = display_label)) +
  geom_vline(xintercept = 0, colour = "#A9AEB4", linewidth = 0.4, linetype = "22") +
  geom_segment(
    data = others,
    aes(x = CI95_low, xend = CI95_high, yend = display_label),
    linewidth = 0.7, colour = neutral_ci
  ) +
  geom_point(
    data = others,
    aes(x = mean_difference),
    size = 1.65, shape = 21, stroke = 0.35,
    fill = "white", colour = neutral_point
  ) +
  geom_segment(
    data = osars,
    aes(x = CI95_low, xend = CI95_high, yend = display_label),
    linewidth = 1.0, colour = accent
  ) +
  geom_point(
    data = osars,
    aes(x = mean_difference),
    size = 2.25, shape = 21, stroke = 0.4,
    fill = accent, colour = "white"
  ) +
  scale_x_continuous(
    limits = c(-0.082, 0.082),
    breaks = c(-0.08, -0.04, 0, 0.04, 0.08),
    labels = c("-0.08", "-0.04", "0", "0.04", "0.08"),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    tag = "B",
    title = "Top seven model configurations",
    subtitle = sprintf(
      "Friedman P = %.3f; Holm-adjusted P < 0.05: %d/21",
      global_p, pairwise_significant
    ),
    x = "Mean ΔC-index vs OSARS (95% CI)",
    y = NULL
  ) +
  base_theme

stem <- file.path(output_root, "Panel_B_Top7_Cindex_differences_vs_OSARS")

if (Sys.getenv("SKIP_SVG", unset = "0") != "1") {
  if (!requireNamespace("svglite", quietly = TRUE)) {
    stop("svglite is required for editable-text SVG output")
  }
  svg_file <- paste0(stem, ".svg")
  ggsave(
    svg_file, pB,
    width = 8.5, height = 7.0, units = "cm",
    device = svglite::svglite, bg = "white"
  )
  svg_text <- readLines(svg_file, warn = FALSE)
  svg_text <- gsub('font-family: "Nimbus Sans"', 'font-family: "Arial"', svg_text, fixed = TRUE)
  svg_text <- gsub('font-family: "Nimbus Roman"', 'font-family: "Times New Roman"', svg_text, fixed = TRUE)
  writeLines(svg_text, svg_file, useBytes = TRUE)
}

if (Sys.getenv("ONLY_SVG", unset = "0") != "1") {
  pdf_device <- if (requireNamespace("Cairo", quietly = TRUE)) Cairo::CairoPDF else grDevices::cairo_pdf
  png_device <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png"
  ggsave(
    paste0(stem, ".pdf"), pB,
    width = 8.5, height = 7.0, units = "cm",
    device = pdf_device, bg = "white"
  )
  ggsave(
    paste0(stem, ".png"), pB,
    width = 8.5, height = 7.0, units = "cm",
    dpi = 600, device = png_device, bg = "white"
  )
}

cat(sprintf("Top-seven comparison figure written; Friedman P = %.6f.\n", global_p))
