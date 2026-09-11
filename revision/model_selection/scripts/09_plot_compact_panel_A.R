#!/usr/bin/env Rscript

args_all <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args_all[grepl("^--file=", args_all)])
script_path <- normalizePath(script_arg[[1]], mustWork = TRUE)
analysis_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

suppressPackageStartupMessages(library(ggplot2))

input_file <- file.path(analysis_root, "all_117_models_5fold_performance.csv")
output_root <- file.path(analysis_root, "figures", "revised_compact_A")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

performance <- read.csv(input_file, check.names = FALSE, stringsAsFactors = FALSE)
required_columns <- c(
  "configuration", "mean_Cindex", "CI95_low", "CI95_high",
  "selected_final_configuration"
)
if (!all(required_columns %in% names(performance))) {
  stop("Missing required columns: ", paste(setdiff(required_columns, names(performance)), collapse = ", "))
}
if (nrow(performance) != 117L) stop("Expected 117 model configurations")

performance <- performance[
  order(-performance$mean_Cindex, performance$configuration),
  , drop = FALSE
]
performance$model_rank_order <- seq_len(nrow(performance))
performance$selected_final_configuration <- as.logical(performance$selected_final_configuration)
highlight_configuration <- "StepCox[forward] + GBM"
performance$highlighted_in_panel_A <- performance$configuration == highlight_configuration
if (sum(performance$highlighted_in_panel_A) != 1L) stop("Highlighted configuration is not unique")

selected <- performance[performance$highlighted_in_panel_A, , drop = FALSE]
others <- performance[!performance$highlighted_in_panel_A, , drop = FALSE]

accent <- "#B64342"
neutral_ci <- "#C4CBD2"
neutral_line <- "#525A63"
neutral_point <- "#343A40"

base_theme <- theme_classic(base_family = "Arial", base_size = 7.2) +
  theme(
    axis.line = element_line(linewidth = 0.35, colour = "#202124"),
    axis.ticks = element_line(linewidth = 0.3, colour = "#202124"),
    axis.ticks.length = grid::unit(1.6, "pt"),
    axis.text = element_text(size = 6.6, colour = "#202124"),
    axis.title = element_text(size = 7.2, colour = "#202124"),
    plot.title = element_text(size = 8.2, hjust = 0.5, face = "plain", margin = margin(b = 1.5)),
    plot.subtitle = element_text(size = 6.5, hjust = 0.5, colour = "#5F6368", margin = margin(b = 4)),
    plot.tag = element_text(
      family = "Times New Roman", face = "bold", size = 9.2,
      colour = "#111111", hjust = 0, vjust = 1
    ),
    plot.tag.position = c(0.012, 0.992),
    plot.margin = margin(t = 3, r = 3, b = 3, l = 3, unit = "pt")
  )

pA <- ggplot(performance, aes(x = model_rank_order, y = mean_Cindex)) +
  geom_hline(
    yintercept = 0.5, colour = "#A9AEB4", linewidth = 0.3,
    linetype = "22"
  ) +
  geom_linerange(
    data = others,
    aes(ymin = CI95_low, ymax = CI95_high),
    linewidth = 0.25, colour = neutral_ci, alpha = 0.62
  ) +
  geom_line(linewidth = 0.38, colour = neutral_line, alpha = 0.9) +
  geom_point(
    data = others,
    size = 0.62, stroke = 0, colour = neutral_point, alpha = 0.9
  ) +
  geom_linerange(
    data = selected,
    aes(ymin = CI95_low, ymax = CI95_high),
    linewidth = 1.0, colour = accent
  ) +
  geom_point(
    data = selected,
    shape = 21, size = 2.45, stroke = 0.45,
    fill = accent, colour = "white"
  ) +
  annotate(
    "text",
    x = selected$model_rank_order + 3.0,
    y = selected$mean_Cindex + 0.014,
    label = highlight_configuration,
    family = "Arial", fontface = "plain", size = 2.15,
    hjust = 0, vjust = 0.5, colour = accent
  ) +
  scale_x_continuous(
    breaks = c(1, 30, 60, 90, 117),
    limits = c(1, 117),
    expand = expansion(mult = c(0.01, 0.018))
  ) +
  scale_y_continuous(
    breaks = c(0.50, 0.60, 0.70, 0.80),
    labels = sprintf("%.2f", c(0.50, 0.60, 0.70, 0.80)),
    limits = c(0.50, 0.81),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    tag = "A",
    title = "5-fold cross-validation within TCGA",
    subtitle = "Mean C-index (95% CI)",
    x = "Model rank",
    y = "Harrell's C-index"
  ) +
  base_theme

stem <- file.path(output_root, "Panel_A_TCGA_117_models_compact_CI")

if (Sys.getenv("SKIP_SVG", unset = "0") != "1") {
  if (!requireNamespace("svglite", quietly = TRUE)) {
    stop("svglite is required for editable-text SVG output")
  }
  svg_file <- paste0(stem, ".svg")
  ggsave(
    svg_file, pA,
    width = 8.5, height = 8.0, units = "cm",
    device = svglite::svglite, bg = "white"
  )
  # Fontconfig substitutes these metrically compatible families on the Linux
  # compute server. Keep journal-requested family names in the editable SVG.
  svg_text <- readLines(svg_file, warn = FALSE)
  svg_text <- gsub('font-family: "Nimbus Sans"', 'font-family: "Arial"', svg_text, fixed = TRUE)
  svg_text <- gsub('font-family: "Nimbus Roman"', 'font-family: "Times New Roman"', svg_text, fixed = TRUE)
  writeLines(svg_text, svg_file, useBytes = TRUE)
}

if (Sys.getenv("ONLY_SVG", unset = "0") != "1") {
  pdf_device <- if (requireNamespace("Cairo", quietly = TRUE)) Cairo::CairoPDF else grDevices::cairo_pdf
  png_device <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png"

  ggsave(
    paste0(stem, ".pdf"), pA,
    width = 8.5, height = 8.0, units = "cm",
    device = pdf_device, bg = "white"
  )
  ggsave(
    paste0(stem, ".png"), pA,
    width = 8.5, height = 8.0, units = "cm",
    dpi = 600, device = png_device, bg = "white"
  )
}

write.csv(
  performance[, c(
    "model_rank_order", "configuration", "mean_Cindex", "CI95_low",
    "CI95_high", "selected_final_configuration", "highlighted_in_panel_A"
  )],
  file.path(output_root, "Panel_A_TCGA_117_models_compact_CI_source_data.csv"),
  row.names = FALSE
)

cat(
  sprintf(
    "Compact Panel A written at 85 mm x 80 mm; selected model: %s.\n",
    highlight_configuration
  )
)
