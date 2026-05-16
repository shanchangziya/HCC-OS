#==== Header ====
# 脚本名称: 19_integrate_all_cells.R
# 目的: 整合所有细胞群形成最终的单细胞对象
# 包含:
#   - 高/低氧化应激肿瘤细胞
#   - 正常上皮细胞
#   - 所有细分的免疫细胞
#   - 所有细分的基质细胞
# 日期: 2026-02-02

rm(list = ls())
gc()
library(Seurat)
library(dplyr)
library(ggplot2)

cat("=== 开始整合所有细胞群 ===\n\n")

# 用于收集所有对象的列表
cell_list <- list()

#==== 1. 上皮细胞（肿瘤+正常，含氧化应激分组）====
cat("=== 1. 处理上皮细胞 ===\n")
load("05osnmf.Rdata")

# 检查seu对象
cat("上皮细胞总数:", ncol(seu), "\n")
cat("\nOS_group_valley分布:\n")
print(table(seu$OS_group_valley))

# 检查cnv_class
if("cnv_class" %in% colnames(seu@meta.data)) {
  cat("\ncnv_class分布:\n")
  print(table(seu$cnv_class, useNA = "ifany"))

  # 创建最终的细胞类型分组
  # 正常细胞 vs 肿瘤细胞 (高/低OS)
  seu$final_celltype <- NA_character_

  # 正常细胞
  seu$final_celltype[seu$cnv_class == "CNV-low"] <- "Normal_Epithelial"

  # 肿瘤细胞 - 高OS
  seu$final_celltype[seu$cnv_class %in% c("CNV-mid", "CNV-high(malignant)") &
                      seu$OS_group_valley == "OS-high"] <- "Tumor_OS-high"

  # 肿瘤细胞 - 低OS
  seu$final_celltype[seu$cnv_class %in% c("CNV-mid", "CNV-high(malignant)") &
                      seu$OS_group_valley == "OS-low"] <- "Tumor_OS-low"

} else {
  # 如果没有cnv_class，仅使用OS分组
  cat("\n警告: 没有cnv_class，仅按OS分组\n")
  seu$final_celltype <- paste0("Epithelial_", seu$OS_group_valley)
}

cat("\n最终上皮细胞分类:\n")
print(table(seu$final_celltype, useNA = "ifany"))

# 保留上皮细胞（去除NA）
seu <- subset(seu, subset = final_celltype != "NA")
seu$final_celltype[is.na(seu$final_celltype)] <- "Epithelial_Other"

cell_list$Epithelial <- seu
cat("上皮细胞准备完成:", ncol(seu), "cells\n")

rm(seu, geneNMF.metaprograms)
gc()

#==== 2. NK细胞 ====
cat("\n=== 2. NK细胞 ===\n")
load("07_NK_annotated.Rdata")
cat("NK细胞数:", ncol(sce.nk), "\n")
cat("NK亚型:\n")
print(table(sce.nk$NK_subtype))

sce.nk$final_celltype <- paste0("NK_", sce.nk$NK_subtype)
cell_list$NK <- sce.nk

rm(sce.nk)
gc()

#==== 3. T细胞（整合细分版本）====
cat("\n=== 3. T细胞 ===\n")

# 3.1 Trm CD8细分
cat("\n3.1 Trm CD8细分\n")
load("17_TrmCD8_refined_annotated.Rdata")
cat("  Trm CD8细胞数:", ncol(sce.trm), "\n")
print(table(sce.trm$Trm_subtype_refined))

sce.trm$final_celltype <- paste0("CD8_Trm_", sce.trm$Trm_subtype_refined)
cell_list$Trm_CD8 <- sce.trm

rm(sce.trm)
gc()

# 3.2 Cytotoxic CD4细分
cat("\n3.2 Cytotoxic CD4细分\n")
load("17_CytoCD4_refined_annotated.Rdata")
cat("  Cytotoxic CD4细胞数:", ncol(sce.cyto), "\n")
print(table(sce.cyto$CytoCD4_subtype_refined))

sce.cyto$final_celltype <- paste0("CD4_Cyto_", sce.cyto$CytoCD4_subtype_refined)
cell_list$Cyto_CD4 <- sce.cyto

rm(sce.cyto)
gc()

# 3.3 其他T细胞亚型（Tfh和Treg）
cat("\n3.3 其他T细胞亚型\n")
load("07_T_annotated.Rdata")
cat("  总T细胞数:", ncol(sce.tc), "\n")
cat("  T亚型分布:\n")
print(table(sce.tc$T_subtype))

# 提取Tfh和Treg
sce.tfh <- subset(sce.tc, subset = T_subtype == "Tfh")
sce.treg <- subset(sce.tc, subset = T_subtype == "Treg")

cat("  Tfh细胞数:", ncol(sce.tfh), "\n")
cat("  Treg细胞数:", ncol(sce.treg), "\n")

sce.tfh$final_celltype <- "CD4_Tfh"
sce.treg$final_celltype <- "CD4_Treg"

cell_list$Tfh <- sce.tfh
cell_list$Treg <- sce.treg

rm(sce.tc, sce.tfh, sce.treg)
gc()

#==== 4. 其他免疫细胞 ====
cat("\n=== 4. 其他免疫细胞 ===\n")

# 4.1 巨噬细胞
cat("\n4.1 巨噬细胞\n")
load("08_Mac_annotated.Rdata")
cat("  巨噬细胞数:", ncol(sce.mac), "\n")
print(table(sce.mac$Mac_subtype))

sce.mac$final_celltype <- paste0("Mac_", sce.mac$Mac_subtype)
cell_list$Macrophage <- sce.mac

rm(sce.mac)
gc()

# 4.2 B细胞
cat("\n4.2 B细胞\n")
load("09_B_annotated.Rdata")
cat("  B细胞数:", ncol(sce.b), "\n")
print(table(sce.b$B_subtype))

sce.b$final_celltype <- paste0("B_", sce.b$B_subtype)
cell_list$B_cell <- sce.b

rm(sce.b)
gc()

# 4.3 DC细胞
cat("\n4.3 DC细胞\n")
load("10_DC_annotated.Rdata")
cat("  DC细胞数:", ncol(sce.dc), "\n")
print(table(sce.dc$DC_subtype))

sce.dc$final_celltype <- paste0("DC_", sce.dc$DC_subtype)
cell_list$DC <- sce.dc

rm(sce.dc)
gc()

# 4.4 中性粒细胞
cat("\n4.4 中性粒细胞\n")
load("13_Neu_annotated.Rdata")
cat("  中性粒细胞数:", ncol(sce.neu), "\n")
print(table(sce.neu$Neu_subtype))

sce.neu$final_celltype <- paste0("Neu_", sce.neu$Neu_subtype)
cell_list$Neutrophil <- sce.neu

rm(sce.neu)
gc()

#==== 5. 基质细胞 ====
cat("\n=== 5. 基质细胞 ===\n")

# 5.1 CAF
cat("\n5.1 CAF\n")
load("11_CAF_annotated.Rdata")
cat("  CAF细胞数:", ncol(sce.caf), "\n")
print(table(sce.caf$CAF_subtype))

sce.caf$final_celltype <- paste0("CAF_", sce.caf$CAF_subtype)
cell_list$CAF <- sce.caf

rm(sce.caf)
gc()

# 5.2 内皮细胞
cat("\n5.2 内皮细胞\n")
load("12_Endo_annotated.Rdata")
cat("  内皮细胞数:", ncol(sce.endo), "\n")
print(table(sce.endo$Endo_subtype))

sce.endo$final_celltype <- paste0("Endo_", sce.endo$Endo_subtype)
cell_list$Endothelial <- sce.endo

rm(sce.endo)
gc()

#==== 6. 整合所有细胞 ====
cat("\n\n=== 6. 整合所有细胞 ===\n")

cat("\n整合前各细胞群统计:\n")
for(name in names(cell_list)) {
  cat(sprintf("  %s: %d cells\n", name, ncol(cell_list[[name]])))
}

# 合并所有对象
cat("\n开始合并...\n")
sce.all.integrated <- merge(
  x = cell_list[[1]],
  y = cell_list[-1],
  add.cell.ids = names(cell_list),
  project = "HCC_Integrated"
)

cat("\n合并完成!\n")
cat("总细胞数:", ncol(sce.all.integrated), "\n")

#==== 7. 添加主要细胞类型分组 ====
cat("\n=== 7. 添加主要细胞类型分组 ===\n")

# 创建主要类型
sce.all.integrated$major_celltype <- case_when(
  grepl("^Tumor", sce.all.integrated$final_celltype) ~ "Tumor",
  grepl("^Normal_Epithelial", sce.all.integrated$final_celltype) ~ "Normal_Epithelial",
  grepl("^CD8_Trm", sce.all.integrated$final_celltype) ~ "T_CD8_Trm",
  grepl("^CD4_Cyto", sce.all.integrated$final_celltype) ~ "T_CD4_Cytotoxic",
  grepl("^CD4_Tfh", sce.all.integrated$final_celltype) ~ "T_CD4_Tfh",
  grepl("^CD4_Treg", sce.all.integrated$final_celltype) ~ "T_CD4_Treg",
  grepl("^NK", sce.all.integrated$final_celltype) ~ "NK",
  grepl("^Mac", sce.all.integrated$final_celltype) ~ "Macrophage",
  grepl("^B_", sce.all.integrated$final_celltype) ~ "B",
  grepl("^DC", sce.all.integrated$final_celltype) ~ "DC",
  grepl("^Neu", sce.all.integrated$final_celltype) ~ "Neutrophil",
  grepl("^CAF", sce.all.integrated$final_celltype) ~ "CAF",
  grepl("^Endo", sce.all.integrated$final_celltype) ~ "Endothelial",
  TRUE ~ "Other"
)

cat("\n主要细胞类型分布:\n")
print(table(sce.all.integrated$major_celltype))

cat("\n详细细胞亚型分布:\n")
celltype_table <- table(sce.all.integrated$final_celltype)
print(sort(celltype_table, decreasing = TRUE))

#==== 8. 保存结果 ====
cat("\n=== 8. 保存结果 ===\n")

save(sce.all.integrated, file = "19_sce.all.integrated.Rdata")
cat("已保存: 19_sce.all.integrated.Rdata\n")

# 保存细胞类型统计
celltype_summary <- data.frame(
  Major_Type = names(table(sce.all.integrated$major_celltype)),
  N_Cells_Major = as.numeric(table(sce.all.integrated$major_celltype))
)

detailed_summary <- data.frame(
  Final_Celltype = names(celltype_table),
  N_Cells = as.numeric(celltype_table)
)
detailed_summary <- detailed_summary[order(detailed_summary$N_Cells, decreasing = TRUE), ]

write.csv(celltype_summary, "19_major_celltype_summary.csv", row.names = FALSE)
write.csv(detailed_summary, "19_detailed_celltype_summary.csv", row.names = FALSE)

cat("\n细胞类型统计已保存:\n")
cat("  - 19_major_celltype_summary.csv\n")
cat("  - 19_detailed_celltype_summary.csv\n")

#==== 9. 生成可视化 ====
cat("\n=== 9. 生成初步可视化 ===\n")

# 简单的统计图
pdf("19_integration_summary.pdf", width = 14, height = 10)

# 主要类型柱状图
p1 <- ggplot(celltype_summary, aes(x = reorder(Major_Type, -N_Cells_Major), y = N_Cells_Major)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = N_Cells_Major), vjust = -0.5, size = 3) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Major Cell Type", y = "Number of Cells",
       title = "Integrated Dataset: Major Cell Types") +
  ylim(0, max(celltype_summary$N_Cells_Major) * 1.1)

print(p1)

# 详细亚型柱状图（Top 30）
top30 <- head(detailed_summary, 30)
p2 <- ggplot(top30, aes(x = reorder(Final_Celltype, N_Cells), y = N_Cells)) +
  geom_bar(stat = "identity", fill = "coral") +
  geom_text(aes(label = N_Cells), hjust = -0.2, size = 2.5) +
  coord_flip() +
  theme_classic(base_size = 10) +
  labs(x = "Cell Subtype", y = "Number of Cells",
       title = "Integrated Dataset: Top 30 Cell Subtypes") +
  xlim(0, max(top30$N_Cells) * 1.15)

print(p2)

# 样本分布（如果有orig.ident）
if("orig.ident" %in% colnames(sce.all.integrated@meta.data)) {
  sample_celltype <- table(sce.all.integrated$major_celltype,
                           sce.all.integrated$orig.ident)
  sample_df <- as.data.frame(sample_celltype)
  colnames(sample_df) <- c("Celltype", "Sample", "Count")

  p3 <- ggplot(sample_df, aes(x = Sample, y = Count, fill = Celltype)) +
    geom_bar(stat = "identity", position = "fill") +
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "right") +
    labs(y = "Proportion", title = "Cell Type Distribution Across Samples") +
    scale_y_continuous(labels = scales::percent)

  print(p3)
}

dev.off()

cat("可视化已保存: 19_integration_summary.pdf\n")

#==== 10. 总结 ====
cat("\n\n=== 整合完成总结 ===\n")
cat("总细胞数:", ncol(sce.all.integrated), "\n")
cat("主要细胞类型数:", length(unique(sce.all.integrated$major_celltype)), "\n")
cat("详细亚型数:", length(unique(sce.all.integrated$final_celltype)), "\n")

if("orig.ident" %in% colnames(sce.all.integrated@meta.data)) {
  cat("样本数:", length(unique(sce.all.integrated$orig.ident)), "\n")
}

cat("\n输出文件:\n")
cat("  - 19_sce.all.integrated.Rdata (整合后的Seurat对象)\n")
cat("  - 19_major_celltype_summary.csv (主要类型统计)\n")
cat("  - 19_detailed_celltype_summary.csv (详细亚型统计)\n")
cat("  - 19_integration_summary.pdf (可视化汇总)\n")

cat("\n=== 全部完成! ===\n")
