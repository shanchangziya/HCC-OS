#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(pheatmap)
  library(RColorBrewer)
  library(gridExtra)
  library(ggpubr)
  library(ComplexHeatmap)
  library(circlize)
})

# Load the data
cat("Loading res_fixed.Rdata...\n")
load("res_fixed.Rdata")

# Create output directory
dir.create("res_fixed_visualizations", showWarnings = FALSE)

# ============================================================
# 1. C-index Analysis (excluding GSE14520)
# ============================================================
cat("\n=== Processing C-index data ===\n")

# Filter out GSE14520
cindex_filtered <- res$Cindex.res %>%
  filter(ID != "GSE14520")

cat("Total observations after filtering:", nrow(cindex_filtered), "\n")
cat("Datasets included:", unique(cindex_filtered$ID), "\n")
cat("Number of models:", length(unique(cindex_filtered$Model)), "\n")

# Summary statistics
cindex_summary <- cindex_filtered %>%
  group_by(ID) %>%
  summarise(
    Mean_Cindex = mean(Cindex, na.rm = TRUE),
    Median_Cindex = median(Cindex, na.rm = TRUE),
    SD_Cindex = sd(Cindex, na.rm = TRUE),
    Min_Cindex = min(Cindex, na.rm = TRUE),
    Max_Cindex = max(Cindex, na.rm = TRUE),
    N_models = n()
  )

print(cindex_summary)

# Find top 10 models for each dataset
top_models <- cindex_filtered %>%
  group_by(ID) %>%
  arrange(desc(Cindex)) %>%
  slice_head(n = 10) %>%
  ungroup()

cat("\nTop 10 models per dataset:\n")
print(top_models)

# ============================================================
# 2. Visualization 1: C-index comparison boxplot
# ============================================================
cat("\n=== Creating C-index boxplot ===\n")

p1 <- ggplot(cindex_filtered, aes(x = ID, y = Cindex, fill = ID)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 1, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3,
               fill = "red", color = "black") +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "C-index Distribution Across Datasets",
       subtitle = "Red diamond = mean, Box = median ± IQR",
       x = "Dataset", y = "C-index") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 10))

ggsave("res_fixed_visualizations/01_cindex_boxplot.pdf", p1,
       width = 8, height = 6)
cat("Saved: 01_cindex_boxplot.pdf\n")

# ============================================================
# 3. Visualization 2: Top 20 models heatmap
# ============================================================
cat("\n=== Creating top models heatmap ===\n")

# Get top 20 models by average C-index across datasets
top20_models <- cindex_filtered %>%
  group_by(Model) %>%
  summarise(Mean_Cindex = mean(Cindex, na.rm = TRUE)) %>%
  arrange(desc(Mean_Cindex)) %>%
  slice_head(n = 20) %>%
  pull(Model)

# Create matrix for heatmap
heatmap_data <- cindex_filtered %>%
  filter(Model %in% top20_models) %>%
  select(Model, ID, Cindex) %>%
  pivot_wider(names_from = ID, values_from = Cindex) %>%
  column_to_rownames("Model") %>%
  as.matrix()

# Create heatmap
pdf("res_fixed_visualizations/02_top20_models_heatmap.pdf",
    width = 10, height = 12)
pheatmap(heatmap_data,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         breaks = seq(min(heatmap_data, na.rm = TRUE),
                     max(heatmap_data, na.rm = TRUE),
                     length.out = 101),
         display_numbers = TRUE,
         number_format = "%.3f",
         fontsize = 10,
         fontsize_number = 8,
         main = "Top 20 Models C-index Heatmap")
dev.off()
cat("Saved: 02_top20_models_heatmap.pdf\n")

# ============================================================
# 4. Visualization 3: Model performance comparison
# ============================================================
cat("\n=== Creating model performance comparison ===\n")

# Extract base model type (before "+")
cindex_filtered$BaseModel <- sapply(strsplit(cindex_filtered$Model, " \\+ "),
                                    function(x) x[1])

# Calculate average C-index by base model
base_model_summary <- cindex_filtered %>%
  group_by(BaseModel, ID) %>%
  summarise(Mean_Cindex = mean(Cindex, na.rm = TRUE),
            SD_Cindex = sd(Cindex, na.rm = TRUE),
            N = n(),
            .groups = "drop") %>%
  arrange(desc(Mean_Cindex))

p2 <- ggplot(base_model_summary,
             aes(x = reorder(BaseModel, Mean_Cindex),
                 y = Mean_Cindex, fill = ID)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8),
           alpha = 0.8) +
  geom_errorbar(aes(ymin = Mean_Cindex - SD_Cindex,
                   ymax = Mean_Cindex + SD_Cindex),
               position = position_dodge(width = 0.8),
               width = 0.3, alpha = 0.6) +
  coord_flip() +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "Base Model Performance Comparison",
       subtitle = "Mean C-index ± SD across datasets",
       x = "Base Model", y = "Mean C-index", fill = "Dataset") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 9),
        legend.position = "bottom")

ggsave("res_fixed_visualizations/03_base_model_comparison.pdf", p2,
       width = 10, height = 8)
cat("Saved: 03_base_model_comparison.pdf\n")

# ============================================================
# 5. Visualization 4: Scatter plot TCGA vs ICGC
# ============================================================
cat("\n=== Creating TCGA vs ICGC scatter plot ===\n")

# Reshape data for comparison
comparison_data <- cindex_filtered %>%
  select(Model, ID, Cindex) %>%
  pivot_wider(names_from = ID, values_from = Cindex)

# Calculate correlation
cor_value <- cor(comparison_data$`TCGA-LIHC`,
                comparison_data$`ICGC-LIRI`,
                use = "complete.obs")

p3 <- ggplot(comparison_data,
             aes(x = `TCGA-LIHC`, y = `ICGC-LIRI`)) +
  geom_point(alpha = 0.6, size = 2, color = "steelblue") +
  geom_smooth(method = "lm", color = "red", se = TRUE, alpha = 0.2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed",
              color = "gray50") +
  annotate("text", x = min(comparison_data$`TCGA-LIHC`, na.rm = TRUE) + 0.05,
           y = max(comparison_data$`ICGC-LIRI`, na.rm = TRUE) - 0.05,
           label = sprintf("Pearson r = %.3f", cor_value),
           hjust = 0, size = 5, fontface = "bold") +
  labs(title = "Model Performance Correlation",
       subtitle = "C-index comparison between TCGA-LIHC and ICGC-LIRI",
       x = "C-index (TCGA-LIHC)", y = "C-index (ICGC-LIRI)") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 10))

ggsave("res_fixed_visualizations/04_tcga_vs_icgc_scatter.pdf", p3,
       width = 8, height = 7)
cat("Saved: 04_tcga_vs_icgc_scatter.pdf\n")

# ============================================================
# 6. Visualization 5: Model complexity analysis
# ============================================================
cat("\n=== Creating model complexity analysis ===\n")

# Count number of components in each model
cindex_filtered$ModelComplexity <- sapply(
  strsplit(cindex_filtered$Model, " \\+ "),
  length
)

complexity_summary <- cindex_filtered %>%
  group_by(ModelComplexity, ID) %>%
  summarise(Mean_Cindex = mean(Cindex, na.rm = TRUE),
            SD_Cindex = sd(Cindex, na.rm = TRUE),
            N = n(),
            .groups = "drop")

p4 <- ggplot(complexity_summary,
             aes(x = factor(ModelComplexity),
                 y = Mean_Cindex, fill = ID)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8),
           alpha = 0.8) +
  geom_errorbar(aes(ymin = Mean_Cindex - SD_Cindex,
                   ymax = Mean_Cindex + SD_Cindex),
               position = position_dodge(width = 0.8),
               width = 0.3) +
  geom_text(aes(label = N), position = position_dodge(width = 0.8),
           vjust = -0.5, size = 3) +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "Model Complexity vs Performance",
       subtitle = "Numbers indicate sample size per group",
       x = "Number of Model Components", y = "Mean C-index",
       fill = "Dataset") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 10),
        legend.position = "bottom")

ggsave("res_fixed_visualizations/05_model_complexity.pdf", p4,
       width = 10, height = 6)
cat("Saved: 05_model_complexity.pdf\n")

# ============================================================
# 7. Risk Score Analysis
# ============================================================
cat("\n=== Analyzing risk scores ===\n")

# Extract risk scores for TCGA-LIHC and ICGC-LIRI
riskscore_data <- list()

for (model_name in names(res$riskscore)) {
  model_risk <- res$riskscore[[model_name]]

  # Check structure and extract data
  if (is.list(model_risk)) {
    for (dataset in c("TCGA-LIHC", "ICGC-LIRI")) {
      if (dataset %in% names(model_risk)) {
        risk_values <- model_risk[[dataset]]

        # Handle different data structures
        if (!is.null(risk_values)) {
          # Try to convert to numeric vector
          if (is.list(risk_values)) {
            risk_values <- unlist(risk_values)
          }

          if (is.numeric(risk_values) && length(risk_values) > 0) {
            riskscore_data[[paste(model_name, dataset, sep = "___")]] <-
              data.frame(
                Model = model_name,
                Dataset = dataset,
                RiskScore = risk_values,
                stringsAsFactors = FALSE
              )
          }
        }
      }
    }
  }
}

# Combine all risk scores
if (length(riskscore_data) > 0) {
  riskscore_df <- do.call(rbind, riskscore_data)
  rownames(riskscore_df) <- NULL

  cat("Risk score data dimensions:", dim(riskscore_df), "\n")
  cat("Column names:", paste(colnames(riskscore_df), collapse = ", "), "\n")

  # Select top 10 models for visualization
  top10_for_risk <- top_models %>%
    filter(ID != "GSE14520") %>%
    pull(Model) %>%
    unique() %>%
    head(10)

  riskscore_top10 <- riskscore_df %>%
    filter(Model %in% top10_for_risk)

  if (nrow(riskscore_top10) > 0) {
    # Visualization 6: Risk score distribution
    p5 <- ggplot(riskscore_top10,
                 aes(x = RiskScore, fill = Dataset)) +
      geom_density(alpha = 0.5) +
      facet_wrap(~Model, scales = "free", ncol = 2) +
      scale_fill_brewer(palette = "Set1") +
      labs(title = "Risk Score Distribution - Top 10 Models",
           x = "Risk Score", y = "Density") +
      theme_bw(base_size = 10) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            legend.position = "bottom",
            strip.text = element_text(size = 8))

    ggsave("res_fixed_visualizations/06_riskscore_distribution.pdf", p5,
           width = 12, height = 14)
    cat("Saved: 06_riskscore_distribution.pdf\n")
  } else {
    cat("No risk score data available for top 10 models\n")
  }
} else {
  cat("No risk score data extracted\n")
}

# ============================================================
# 8. Significant Genes Analysis
# ============================================================
cat("\n=== Analyzing significant genes ===\n")

sig_genes <- res$Sig.genes
cat("Number of significant genes:", length(sig_genes), "\n")
cat("First 20 genes:", paste(head(sig_genes, 20), collapse = ", "), "\n")

# Save gene list
write.table(sig_genes,
           "res_fixed_visualizations/significant_genes.txt",
           row.names = FALSE, col.names = FALSE, quote = FALSE)
cat("Saved: significant_genes.txt\n")

# ============================================================
# 9. Summary Statistics and Export
# ============================================================
cat("\n=== Creating summary report ===\n")

# Overall best models
best_models_overall <- cindex_filtered %>%
  group_by(Model) %>%
  summarise(
    Mean_Cindex = mean(Cindex, na.rm = TRUE),
    SD_Cindex = sd(Cindex, na.rm = TRUE),
    TCGA_Cindex = Cindex[ID == "TCGA-LIHC"][1],
    ICGC_Cindex = Cindex[ID == "ICGC-LIRI"][1]
  ) %>%
  arrange(desc(Mean_Cindex)) %>%
  head(20)

write.csv(best_models_overall,
         "res_fixed_visualizations/top20_models_summary.csv",
         row.names = FALSE)
cat("Saved: top20_models_summary.csv\n")

# Dataset-specific best models
best_by_dataset <- cindex_filtered %>%
  group_by(ID) %>%
  arrange(desc(Cindex)) %>%
  slice_head(n = 10) %>%
  ungroup()

write.csv(best_by_dataset,
         "res_fixed_visualizations/top10_models_by_dataset.csv",
         row.names = FALSE)
cat("Saved: top10_models_by_dataset.csv\n")

# ============================================================
# 10. Combined Summary Plot
# ============================================================
cat("\n=== Creating combined summary plot ===\n")

# Create a multi-panel summary figure
p_summary1 <- ggplot(cindex_filtered,
                     aes(x = ID, y = Cindex, fill = ID)) +
  geom_violin(alpha = 0.6) +
  geom_boxplot(width = 0.2, alpha = 0.8) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "C-index Distribution", x = "", y = "C-index") +
  theme_bw(base_size = 10) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"))

p_summary2 <- ggplot(comparison_data,
                     aes(x = `TCGA-LIHC`, y = `ICGC-LIRI`)) +
  geom_point(alpha = 0.5, size = 1.5, color = "steelblue") +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(title = "Dataset Correlation",
       x = "TCGA-LIHC", y = "ICGC-LIRI") +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

p_summary3 <- ggplot(complexity_summary,
                     aes(x = factor(ModelComplexity),
                         y = Mean_Cindex, fill = ID)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "Model Complexity",
       x = "Components", y = "Mean C-index") +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "bottom",
        legend.title = element_blank())

# Top 10 models bar plot
top10_avg <- cindex_filtered %>%
  group_by(Model) %>%
  summarise(Mean_Cindex = mean(Cindex, na.rm = TRUE)) %>%
  arrange(desc(Mean_Cindex)) %>%
  head(10)

p_summary4 <- ggplot(top10_avg,
                     aes(x = reorder(Model, Mean_Cindex),
                         y = Mean_Cindex)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
  coord_flip() +
  labs(title = "Top 10 Models", x = "", y = "Mean C-index") +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.y = element_text(size = 7))

combined_plot <- ggarrange(p_summary1, p_summary2, p_summary3, p_summary4,
                          ncol = 2, nrow = 2,
                          labels = c("A", "B", "C", "D"))

ggsave("res_fixed_visualizations/00_summary_combined.pdf", combined_plot,
       width = 14, height = 12)
cat("Saved: 00_summary_combined.pdf\n")

# ============================================================
# 11. Final Summary Report
# ============================================================
cat("\n" , rep("=", 60), "\n", sep = "")
cat("VISUALIZATION COMPLETE\n")
cat(rep("=", 60), "\n\n", sep = "")

cat("Summary Statistics:\n")
cat("-------------------\n")
cat("Total models analyzed:", length(unique(cindex_filtered$Model)), "\n")
cat("Datasets included: TCGA-LIHC, ICGC-LIRI\n")
cat("Dataset excluded: GSE14520\n\n")

cat("C-index Statistics:\n")
print(cindex_summary)

cat("\n\nTop 5 Models (by average C-index):\n")
print(head(best_models_overall, 5))

cat("\n\nSignificant Genes:", length(sig_genes), "genes\n")
cat("First 10:", paste(head(sig_genes, 10), collapse = ", "), "\n")

cat("\n\nOutput Files Generated:\n")
cat("------------------------\n")
cat("1. 00_summary_combined.pdf - Multi-panel overview\n")
cat("2. 01_cindex_boxplot.pdf - C-index distribution\n")
cat("3. 02_top20_models_heatmap.pdf - Top 20 models heatmap\n")
cat("4. 03_base_model_comparison.pdf - Base model performance\n")
cat("5. 04_tcga_vs_icgc_scatter.pdf - Dataset correlation\n")
cat("6. 05_model_complexity.pdf - Complexity analysis\n")
cat("7. 06_riskscore_distribution.pdf - Risk score distributions\n")
cat("8. top20_models_summary.csv - Top 20 models data\n")
cat("9. top10_models_by_dataset.csv - Dataset-specific top models\n")
cat("10. significant_genes.txt - List of significant genes\n")

cat("\nAll files saved in: res_fixed_visualizations/\n")
cat(rep("=", 60), "\n", sep = "")

cat("\nVisualization script completed successfully!\n")
