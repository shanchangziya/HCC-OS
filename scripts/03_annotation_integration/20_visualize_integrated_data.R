#==== Header ====
# 脚本名称: 20_visualize_integrated_data.R
# 目的: 为整合数据生成详细可视化
# 日期: 2026-02-02

rm(list = ls())
gc()
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

cat("=== 加载整合数据并生成可视化 ===\n")

# 加载数据
load("19_sce.all.integrated.Rdata")

cat("总细胞数:", ncol(sce.all.integrated), "\n")
cat("总基因数:", nrow(sce.all.integrated), "\n")

# 读取统计数据
major_summary <- read.csv("19_major_celltype_summary.csv")
detailed_summary <- read.csv("19_detailed_celltype_summary.csv")

# 创建配色方案
major_colors <- c(
  "Tumor" = "#E15759",
  "Normal_Epithelial" = "#76B7B2",
  "T_CD8_Trm" = "#4E79A7",
  "T_CD4_Cytotoxic" = "#F28E2B",
  "T_CD4_Tfh" = "#59A14F",
  "T_CD4_Treg" = "#EDC948",
  "NK" = "#B07AA1",
  "Macrophage" = "#FF9DA7",
  "B" = "#9C755F",
  "DC" = "#BAB0AC",
  "Neutrophil" = "#17BECF",
  "CAF" = "#BCBD22",
  "Endothelial" = "#8C564B"
)

#==== 生成可视化 ====
pdf("20_integrated_data_visualization.pdf", width = 16, height = 12)

# 1. 主要类型柱状图
cat("\n生成主要类型柱状图...\n")
p1 <- ggplot(major_summary, aes(x = reorder(Major_Type, -N_Cells_Major), y = N_Cells_Major)) +
  geom_bar(stat = "identity", fill = "steelblue", width = 0.7) +
  geom_text(aes(label = format(N_Cells_Major, big.mark = ",")),
            vjust = -0.5, size = 3.5) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(x = "Major Cell Type", y = "Number of Cells",
       title = "HCC Integrated Dataset: Major Cell Types",
       subtitle = paste0("Total: ", format(sum(major_summary$N_Cells_Major), big.mark = ","), " cells")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)),
                    labels = scales::comma)

print(p1)

# 2. 主要类型饼图
cat("生成主要类型饼图...\n")
major_summary$Percentage <- major_summary$N_Cells_Major / sum(major_summary$N_Cells_Major) * 100
major_summary$Label <- paste0(major_summary$Major_Type, "\n",
                              round(major_summary$Percentage, 1), "%")

p2 <- ggplot(major_summary, aes(x = "", y = N_Cells_Major, fill = Major_Type)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = major_colors) +
  theme_void(base_size = 12) +
  labs(title = "Cell Type Distribution") +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold"))

print(p2)

# 3. 详细亚型柱状图（Top 30）
cat("生成详细亚型柱状图...\n")
top30 <- head(detailed_summary, 30)

p3 <- ggplot(top30, aes(y = reorder(Final_Celltype, N_Cells), x = N_Cells)) +
  geom_bar(stat = "identity", fill = "coral", width = 0.75) +
  geom_text(aes(label = format(N_Cells, big.mark = ",")),
            hjust = -0.1, size = 3) +
  theme_classic(base_size = 11) +
  labs(y = "Cell Subtype", x = "Number of Cells",
       title = "Top 30 Cell Subtypes") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15)),
                    labels = scales::comma)

print(p3)

# 4. 肿瘤细胞氧化应激分组
cat("生成肿瘤细胞氧化应激分组图...\n")
tumor_df <- detailed_summary[grepl("^Tumor", detailed_summary$Final_Celltype), ]

if(nrow(tumor_df) > 0) {
  p4 <- ggplot(tumor_df, aes(x = Final_Celltype, y = N_Cells, fill = Final_Celltype)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = format(N_Cells, big.mark = ",")), vjust = -0.5) +
    scale_fill_manual(values = c("Tumor_OS-high" = "#D62728",
                                  "Tumor_OS-low" = "#FFA07A")) +
    theme_classic(base_size = 14) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 0)) +
    labs(x = NULL, y = "Number of Cells",
         title = "Tumor Cells: Oxidative Stress Stratification") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)),
                      labels = scales::comma)

  print(p4)
}

# 5. T细胞亚型分布
cat("生成T细胞亚型分布图...\n")
t_df <- detailed_summary[grepl("^CD4_|^CD8_", detailed_summary$Final_Celltype), ]

if(nrow(t_df) > 0) {
  t_df$T_major <- case_when(
    grepl("CD8_Trm", t_df$Final_Celltype) ~ "CD8 Trm",
    grepl("CD4_Cyto", t_df$Final_Celltype) ~ "CD4 Cytotoxic",
    grepl("CD4_Tfh", t_df$Final_Celltype) ~ "CD4 Tfh",
    grepl("CD4_Treg", t_df$Final_Celltype) ~ "CD4 Treg",
    TRUE ~ "Other"
  )

  p5 <- ggplot(t_df, aes(x = T_major, y = N_Cells, fill = Final_Celltype)) +
    geom_bar(stat = "identity") +
    theme_classic(base_size = 13) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "right",
          legend.text = element_text(size = 9)) +
    labs(x = "T Cell Major Type", y = "Number of Cells",
         title = "T Cell Subtype Distribution",
         fill = "Detailed Subtype") +
    scale_y_continuous(labels = scales::comma)

  print(p5)
}

# 6. 免疫细胞总览
cat("生成免疫细胞总览图...\n")
immune_df <- major_summary[major_summary$Major_Type %in% c(
  "T_CD8_Trm", "T_CD4_Cytotoxic", "T_CD4_Tfh", "T_CD4_Treg",
  "NK", "Macrophage", "B", "DC", "Neutrophil"
), ]

if(nrow(immune_df) > 0) {
  p6 <- ggplot(immune_df, aes(x = reorder(Major_Type, -N_Cells_Major),
                              y = N_Cells_Major,
                              fill = Major_Type)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = format(N_Cells_Major, big.mark = ",")),
              vjust = -0.5, size = 3.5) +
    scale_fill_manual(values = major_colors) +
    theme_classic(base_size = 13) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") +
    labs(x = "Immune Cell Type", y = "Number of Cells",
         title = "Immune Cell Composition") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)),
                      labels = scales::comma)

  print(p6)
}

# 7. 样本分布（如果有orig.ident）
if("orig.ident" %in% colnames(sce.all.integrated@meta.data)) {
  cat("生成样本分布图...\n")

  sample_celltype <- table(sce.all.integrated$major_celltype,
                           sce.all.integrated$orig.ident)
  sample_df <- as.data.frame(sample_celltype)
  colnames(sample_df) <- c("Celltype", "Sample", "Count")

  # 堆叠柱状图
  p7 <- ggplot(sample_df, aes(x = Sample, y = Count, fill = Celltype)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = major_colors) +
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "right") +
    labs(y = "Number of Cells",
         title = "Cell Type Distribution Across Samples") +
    scale_y_continuous(labels = scales::comma)

  print(p7)

  # 比例柱状图
  p8 <- ggplot(sample_df, aes(x = Sample, y = Count, fill = Celltype)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_manual(values = major_colors) +
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "right") +
    labs(y = "Proportion",
         title = "Cell Type Proportions Across Samples") +
    scale_y_continuous(labels = scales::percent)

  print(p8)
}

# 8. 数据质量统计
if(all(c("nCount_RNA", "nFeature_RNA") %in% colnames(sce.all.integrated@meta.data))) {
  cat("生成数据质量统计图...\n")

  qc_df <- data.frame(
    Major_Type = sce.all.integrated$major_celltype,
    nCount = sce.all.integrated$nCount_RNA,
    nFeature = sce.all.integrated$nFeature_RNA
  )

  p9 <- ggplot(qc_df, aes(x = Major_Type, y = nCount, fill = Major_Type)) +
    geom_violin(scale = "width") +
    geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", alpha = 0.5) +
    scale_fill_manual(values = major_colors) +
    scale_y_log10(labels = scales::comma) +
    theme_classic(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") +
    labs(x = NULL, y = "UMI Count (log10)",
         title = "UMI Count Distribution by Cell Type")

  print(p9)

  p10 <- ggplot(qc_df, aes(x = Major_Type, y = nFeature, fill = Major_Type)) +
    geom_violin(scale = "width") +
    geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", alpha = 0.5) +
    scale_fill_manual(values = major_colors) +
    scale_y_log10(labels = scales::comma) +
    theme_classic(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") +
    labs(x = NULL, y = "Gene Count (log10)",
         title = "Gene Count Distribution by Cell Type")

  print(p10)
}

dev.off()

cat("\n可视化完成！已保存: 20_integrated_data_visualization.pdf\n")

#==== 生成总结报告 ====
cat("\n=== 生成总结报告 ===\n")

sink("20_integration_summary_report.txt")

cat("="  ,rep("=", 70), "=\n", sep = "")
cat("   HCC Single-Cell RNA-seq Integrated Dataset Summary Report\n")
cat("="  ,rep("=", 70), "=\n", sep = "")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("1. OVERALL STATISTICS\n")
cat(rep("-", 72), "\n", sep = "")
cat(sprintf("Total Cells:            %s\n", format(ncol(sce.all.integrated), big.mark = ",")))
cat(sprintf("Total Genes:            %s\n", format(nrow(sce.all.integrated), big.mark = ",")))
cat(sprintf("Major Cell Types:       %d\n", length(unique(sce.all.integrated$major_celltype))))
cat(sprintf("Detailed Subtypes:      %d\n\n", length(unique(sce.all.integrated$final_celltype))))

if("orig.ident" %in% colnames(sce.all.integrated@meta.data)) {
  cat(sprintf("Number of Samples:      %d\n", length(unique(sce.all.integrated$orig.ident))))
  cat("Sample IDs:            ", paste(sort(unique(sce.all.integrated$orig.ident)), collapse = ", "), "\n\n")
}

cat("\n2. MAJOR CELL TYPE COMPOSITION\n")
cat(rep("-", 72), "\n", sep = "")
major_sorted <- major_summary[order(-major_summary$N_Cells_Major), ]
for(i in 1:nrow(major_sorted)) {
  pct <- major_sorted$N_Cells_Major[i] / sum(major_summary$N_Cells_Major) * 100
  cat(sprintf("%-25s: %8s cells (%5.2f%%)\n",
              major_sorted$Major_Type[i],
              format(major_sorted$N_Cells_Major[i], big.mark = ","),
              pct))
}

cat("\n3. EPITHELIAL CELLS (TUMOR & NORMAL)\n")
cat(rep("-", 72), "\n", sep = "")
epi_cells <- detailed_summary[grepl("^Tumor|^Normal_Epi", detailed_summary$Final_Celltype), ]
for(i in 1:nrow(epi_cells)) {
  pct <- epi_cells$N_Cells[i] / sum(epi_cells$N_Cells) * 100
  cat(sprintf("%-25s: %8s cells (%5.2f%%)\n",
              epi_cells$Final_Celltype[i],
              format(epi_cells$N_Cells[i], big.mark = ","),
              pct))
}

cat("\n4. T CELL SUBTYPES\n")
cat(rep("-", 72), "\n", sep = "")
t_cells <- detailed_summary[grepl("^CD4_|^CD8_", detailed_summary$Final_Celltype), ]
t_cells <- t_cells[order(-t_cells$N_Cells), ]
for(i in 1:nrow(t_cells)) {
  pct <- t_cells$N_Cells[i] / sum(t_cells$N_Cells) * 100
  cat(sprintf("%-30s: %8s cells (%5.2f%%)\n",
              t_cells$Final_Celltype[i],
              format(t_cells$N_Cells[i], big.mark = ","),
              pct))
}

cat("\n5. MYELOID CELLS\n")
cat(rep("-", 72), "\n", sep = "")
myeloid <- detailed_summary[grepl("^Mac_|^DC_|^Neu_", detailed_summary$Final_Celltype), ]
myeloid <- myeloid[order(-myeloid$N_Cells), ]
for(i in 1:nrow(myeloid)) {
  cat(sprintf("%-25s: %8s cells\n",
              myeloid$Final_Celltype[i],
              format(myeloid$N_Cells[i], big.mark = ",")))
}

cat("\n6. OTHER IMMUNE CELLS\n")
cat(rep("-", 72), "\n", sep = "")
other_immune <- detailed_summary[grepl("^B_|^NK_", detailed_summary$Final_Celltype), ]
other_immune <- other_immune[order(-other_immune$N_Cells), ]
for(i in 1:nrow(other_immune)) {
  cat(sprintf("%-25s: %8s cells\n",
              other_immune$Final_Celltype[i],
              format(other_immune$N_Cells[i], big.mark = ",")))
}

cat("\n7. STROMAL CELLS\n")
cat(rep("-", 72), "\n", sep = "")
stromal <- detailed_summary[grepl("^CAF_|^Endo_", detailed_summary$Final_Celltype), ]
stromal <- stromal[order(-stromal$N_Cells), ]
for(i in 1:nrow(stromal)) {
  cat(sprintf("%-25s: %8s cells\n",
              stromal$Final_Celltype[i],
              format(stromal$N_Cells[i], big.mark = ",")))
}

cat("\n")
cat(rep("=", 72), "\n", sep = "")
cat("End of Report\n")
cat(rep("=", 72), "\n", sep = "")

sink()

cat("总结报告已保存: 20_integration_summary_report.txt\n")

cat("\n=== 全部完成 ===\n")
cat("\n输出文件:\n")
cat("  - 20_integrated_data_visualization.pdf (详细可视化)\n")
cat("  - 20_integration_summary_report.txt (文本总结报告)\n")
