#!/usr/bin/env Rscript

args_all <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args_all[grepl("^--file=", args_all)])
script_path <- normalizePath(script_arg[[1]], mustWork = TRUE)
analysis_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
source(file.path(analysis_root, "scripts", "common.R"))

library(ggplot2)
dir.create(figures_root, recursive = TRUE, showWarnings = FALSE)

performance <- read.csv(file.path(analysis_root, "all_117_models_5fold_performance.csv"), check.names = FALSE)
folds <- read.csv(file.path(processed_root, "all_model_fold_level_results.csv"), check.names = FALSE)
selected <- performance$configuration[performance$selected_final_configuration]
if (length(selected) != 1L) stop("Expected one selected model")

ordered_names <- performance$configuration[order(performance$mean_Cindex, na.last = TRUE)]
performance$configuration_factor <- factor(performance$configuration, levels = ordered_names)
folds$configuration_factor <- factor(folds$configuration, levels = ordered_names)
folds$selected_final_configuration <- folds$configuration == selected

range_rows <- do.call(rbind, lapply(split(folds, folds$configuration), function(d) data.frame(
  configuration = d$configuration[1],
  min_Cindex = min(d$Cindex, na.rm = TRUE),
  max_Cindex = max(d$Cindex, na.rm = TRUE),
  stringsAsFactors = FALSE
)))
range_rows$configuration_factor <- factor(range_rows$configuration, levels = ordered_names)

base_theme <- theme_classic(base_family = "Arial", base_size = 8) +
  theme(
    axis.line = element_line(linewidth = 0.35, colour = "black"),
    axis.ticks = element_line(linewidth = 0.3, colour = "black"),
    plot.margin = margin(5, 6, 5, 6),
    plot.title = element_text(size = 9, hjust = 0.5, face = "plain"),
    axis.title = element_text(size = 8),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 7)
  )

pA <- ggplot() +
  annotate(
    "rect", xmin = -Inf, xmax = Inf,
    ymin = which(ordered_names == selected) - 0.45,
    ymax = which(ordered_names == selected) + 0.45,
    fill = "#F6CFCB", alpha = 0.38
  ) +
  geom_segment(
    data = range_rows,
    aes(x = min_Cindex, xend = max_Cindex, y = configuration_factor, yend = configuration_factor),
    linewidth = 0.25, colour = "#B8BEC6"
  ) +
  geom_point(
    data = folds[!folds$selected_final_configuration, , drop = FALSE],
    aes(x = Cindex, y = configuration_factor),
    size = 0.85, colour = "#A9B1BA", alpha = 0.75,
    position = position_jitter(height = 0.11, width = 0, seed = master_seed)
  ) +
  geom_point(
    data = folds[folds$selected_final_configuration, , drop = FALSE],
    aes(x = Cindex, y = configuration_factor),
    size = 1.05, colour = "#B64342", alpha = 0.85,
    position = position_jitter(height = 0.11, width = 0, seed = master_seed)
  ) +
  geom_point(
    data = performance,
    aes(x = mean_Cindex, y = configuration_factor, colour = selected_final_configuration),
    size = 1.45, shape = 18
  ) +
  scale_colour_manual(values = c(`FALSE` = "#252A31", `TRUE` = "#C43C39"), guide = "none") +
  scale_x_continuous(limits = c(
    max(0.45, min(folds$Cindex, na.rm = TRUE) - 0.02),
    min(0.95, max(folds$Cindex, na.rm = TRUE) + 0.03)
  ), expand = expansion(mult = c(0, 0))) +
  labs(x = "Harrell's C-index", y = NULL, title = "5-fold cross-validation within TCGA") +
  base_theme +
  theme(
    axis.text.y = element_text(size = 6.2, colour = "#252A31", margin = margin(r = 2)),
    axis.text.x = element_text(size = 7),
    plot.tag = element_text(family = "Times New Roman", face = "bold", size = 10),
    plot.tag.position = c(0, 1)
  ) +
  labs(tag = "A")

selected_folds <- folds[folds$configuration == selected, , drop = FALSE]
selected_mean <- performance$mean_Cindex[performance$configuration == selected]
pB <- ggplot(selected_folds, aes(x = factor(fold), y = Cindex, group = 1)) +
  geom_hline(yintercept = selected_mean, linewidth = 0.45, linetype = 2, colour = "#C43C39") +
  geom_line(linewidth = 0.45, colour = "#6B737C") +
  geom_point(size = 2.2, shape = 21, stroke = 0.4, fill = "#C43C39", colour = "black") +
  annotate(
    "text", x = 5.3, y = selected_mean, label = sprintf("Mean %.3f", selected_mean),
    hjust = 0, vjust = -0.55, family = "Arial", size = 2.6, colour = "#B64342"
  ) +
  scale_x_discrete(expand = expansion(add = c(0.35, 1.25))) +
  scale_y_continuous(limits = range(c(selected_folds$Cindex, selected_mean)) + c(-0.035, 0.035)) +
  labs(x = "Fold", y = "Harrell's C-index", title = selected) +
  base_theme +
  theme(
    plot.tag = element_text(family = "Times New Roman", face = "bold", size = 10),
    plot.tag.position = c(0, 1),
    legend.position = "none"
  ) +
  labs(tag = "B")

save_panel <- function(plot, stem, width_cm, height_cm) {
  svg_device <- if (requireNamespace("svglite", quietly = TRUE)) svglite::svglite else grDevices::svg
  ggsave(file.path(figures_root, paste0(stem, ".svg")), plot,
         width = width_cm, height = height_cm, units = "cm", device = svg_device)
  ggsave(file.path(figures_root, paste0(stem, ".pdf")), plot,
         width = width_cm, height = height_cm, units = "cm", device = cairo_pdf)
  ggsave(file.path(figures_root, paste0(stem, ".png")), plot,
         width = width_cm, height = height_cm, units = "cm", dpi = 300, bg = "white")
}
save_panel(pA, "Panel_A_117_models_TCGA_5fold_Cindex", 17, 52)
save_panel(pB, "Panel_B_selected_model_fold_Cindex", 17, 8)

draw_combined <- function(device_fun, filename, width_cm = 17, height_cm = 60) {
  device_fun(filename, width = width_cm / 2.54, height = height_cm / 2.54)
  grid::grid.newpage()
  layout <- grid::grid.layout(nrow = 2, ncol = 1, heights = grid::unit(c(52, 8), "cm"))
  vp <- grid::viewport(layout = layout)
  grid::pushViewport(vp)
  print(pA, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(pB, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
  grid::popViewport()
  grDevices::dev.off()
}
draw_combined(grDevices::cairo_pdf, file.path(figures_root, "SuppFig_117_models_TCGA_5fold_Cindex.pdf"))
png_fun <- function(filename, width, height) grDevices::png(filename, width = width, height = height, units = "in", res = 300, type = "cairo")
draw_combined(png_fun, file.path(figures_root, "SuppFig_117_models_TCGA_5fold_Cindex.png"))

svg_fun <- if (requireNamespace("svglite", quietly = TRUE)) {
  function(filename, width, height) svglite::svglite(filename, width = width, height = height)
} else {
  function(filename, width, height) grDevices::svg(filename, width = width, height = height)
}
draw_combined(svg_fun, file.path(figures_root, "SuppFig_117_models_TCGA_5fold_Cindex.svg"))

write.csv(folds[, c("configuration", "fold", "Cindex", "final_input_feature_count")],
          file.path(processed_root, "SuppFig_model_selection_source_data.csv"), row.names = FALSE)

cat("Supplementary panels written at 17-cm width.\n")
