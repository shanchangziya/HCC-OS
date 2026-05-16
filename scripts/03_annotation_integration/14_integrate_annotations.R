#==== Header ====
# 脚本名称: 14_integrate_annotations.R
# 目的: 整合所有细胞类型的亚群注释回sce.all
# 日期: 2026-02-01
# 说明: 加载各细胞类型的分析结果，提取亚型注释，写入sce.all的celltype_fine列

#==== 环境设置 ====
rm(list = ls())
gc()
library(Seurat)
library(dplyr)
library(ggplot2)
library(SCP)

mycols <- c(
  "#4E79A7","#F28E2B","#E15759","#76B7B2","#59A14F",
  "#EDC948","#B07AA1","#FF9DA7","#9C755F","#BAB0AC",
  "#1F77B4","#FF7F0E","#2CA02C","#D62728","#9467BD",
  "#8C564B","#E377C2","#7F7F7F","#BCBD22","#17BECF",
  "#2E4057","#048A81","#F7E1E4","#E4572E","#315C75",
  "#3D405B","#81B29A","#F2CC8F","#A77B93","#5B4A9E"
)

#==== 1. 加载主对象 ====
cat("=== 加载主对象 ===\n")
load("sceall.Rdata")
cat("总细胞数:", ncol(sce.all), "\n")
print(table(sce.all$celltype))

# 初始化celltype_fine
sce.all$celltype_fine <- sce.all$celltype  # 默认保持大类注释

#==== 2. 逐一加载各细胞类型注释并映射 ====
cat("\n=== 加载各细胞类型亚群注释 ===\n")

# --- NK细胞 ---
cat("\n[1/7] 加载NK细胞注释...\n")
if(file.exists("07_NK_annotated.Rdata")) {
  load("07_NK_annotated.Rdata")
  cat("NK细胞数:", ncol(sce.nk), "\n")
  cat("NK亚型分布:\n")
  print(table(sce.nk$NK_subtype))

  # 映射到sce.all
  common_nk <- intersect(colnames(sce.all), colnames(sce.nk))
  sce.all$celltype_fine[common_nk] <- paste0("NK:", sce.nk$NK_subtype[common_nk])
  cat("映射到sce.all:", length(common_nk), "个细胞\n")
  rm(sce.nk); gc()
} else {
  cat("警告: 07_NK_annotated.Rdata 未找到\n")
}

# --- T细胞 ---
cat("\n[2/7] 加载T细胞注释...\n")
if(file.exists("07_T_annotated.Rdata")) {
  load("07_T_annotated.Rdata")
  cat("T细胞数:", ncol(sce.tc), "\n")
  cat("T亚型分布:\n")
  print(table(sce.tc$T_subtype))

  common_tc <- intersect(colnames(sce.all), colnames(sce.tc))
  sce.all$celltype_fine[common_tc] <- paste0("T:", sce.tc$T_subtype[common_tc])
  cat("映射到sce.all:", length(common_tc), "个细胞\n")
  rm(sce.tc); gc()
} else {
  cat("警告: 07_T_annotated.Rdata 未找到\n")
}

# --- 巨噬细胞 ---
cat("\n[3/7] 加载巨噬细胞注释...\n")
if(file.exists("08_Mac_annotated.Rdata")) {
  load("08_Mac_annotated.Rdata")
  cat("Mac细胞数:", ncol(sce.mac), "\n")
  cat("Mac亚型分布:\n")
  print(table(sce.mac$Mac_subtype))

  common_mac <- intersect(colnames(sce.all), colnames(sce.mac))
  sce.all$celltype_fine[common_mac] <- sce.mac$Mac_subtype[common_mac]
  cat("映射到sce.all:", length(common_mac), "个细胞\n")
  rm(sce.mac); gc()
} else {
  cat("警告: 08_Mac_annotated.Rdata 未找到\n")
}

# --- B细胞 ---
cat("\n[4/7] 加载B细胞注释...\n")
if(file.exists("09_B_annotated.Rdata")) {
  load("09_B_annotated.Rdata")
  cat("B细胞数:", ncol(sce.b), "\n")
  cat("B亚型分布:\n")
  print(table(sce.b$B_subtype))

  common_b <- intersect(colnames(sce.all), colnames(sce.b))
  sce.all$celltype_fine[common_b] <- sce.b$B_subtype[common_b]
  cat("映射到sce.all:", length(common_b), "个细胞\n")
  rm(sce.b); gc()
} else {
  cat("警告: 09_B_annotated.Rdata 未找到\n")
}

# --- DC ---
cat("\n[5/7] 加载DC注释...\n")
if(file.exists("10_DC_annotated.Rdata")) {
  load("10_DC_annotated.Rdata")
  cat("DC细胞数:", ncol(sce.dc), "\n")
  cat("DC亚型分布:\n")
  print(table(sce.dc$DC_subtype))

  common_dc <- intersect(colnames(sce.all), colnames(sce.dc))
  sce.all$celltype_fine[common_dc] <- sce.dc$DC_subtype[common_dc]
  cat("映射到sce.all:", length(common_dc), "个细胞\n")
  rm(sce.dc); gc()
} else {
  cat("警告: 10_DC_annotated.Rdata 未找到\n")
}

# --- CAF ---
cat("\n[6/7] 加载CAF注释...\n")
if(file.exists("11_CAF_annotated.Rdata")) {
  load("11_CAF_annotated.Rdata")
  cat("CAF细胞数:", ncol(sce.caf), "\n")
  cat("CAF亚型分布:\n")
  print(table(sce.caf$CAF_subtype))

  common_caf <- intersect(colnames(sce.all), colnames(sce.caf))
  sce.all$celltype_fine[common_caf] <- sce.caf$CAF_subtype[common_caf]
  cat("映射到sce.all:", length(common_caf), "个细胞\n")
  rm(sce.caf); gc()
} else {
  cat("警告: 11_CAF_annotated.Rdata 未找到\n")
}

# --- 内皮细胞 ---
cat("\n[7/7] 加载内皮细胞注释...\n")
if(file.exists("12_Endo_annotated.Rdata")) {
  load("12_Endo_annotated.Rdata")
  cat("Endo细胞数:", ncol(sce.endo), "\n")
  cat("Endo亚型分布:\n")
  print(table(sce.endo$Endo_subtype))

  common_endo <- intersect(colnames(sce.all), colnames(sce.endo))
  sce.all$celltype_fine[common_endo] <- sce.endo$Endo_subtype[common_endo]
  cat("映射到sce.all:", length(common_endo), "个细胞\n")
  rm(sce.endo); gc()
} else {
  cat("警告: 12_Endo_annotated.Rdata 未找到\n")
}

# --- 中性粒细胞 ---
cat("\n[+1] 加载中性粒细胞注释...\n")
if(file.exists("13_Neu_annotated.Rdata")) {
  load("13_Neu_annotated.Rdata")
  cat("Neu细胞数:", ncol(sce.neu), "\n")
  cat("Neu亚型分布:\n")
  print(table(sce.neu$Neu_subtype))

  common_neu <- intersect(colnames(sce.all), colnames(sce.neu))
  sce.all$celltype_fine[common_neu] <- sce.neu$Neu_subtype[common_neu]
  cat("映射到sce.all:", length(common_neu), "个细胞\n")
  rm(sce.neu); gc()
} else {
  cat("警告: 13_Neu_annotated.Rdata 未找到\n")
}

#==== 3. 整合结果统计 ====
cat("\n=== 整合结果 ===\n")
cat("\n细胞类型_fine 全局分布:\n")
print(sort(table(sce.all$celltype_fine), decreasing = TRUE))

cat("\n每种细胞类型的亚群数量:\n")
for(ct in unique(sce.all$celltype)) {
  subtypes <- unique(sce.all$celltype_fine[sce.all$celltype == ct])
  cat(sprintf("  %s: %d个亚群 -> %s\n", ct, length(subtypes), paste(subtypes, collapse = ", ")))
}

#==== 4. 可视化 ====
cat("\n=== 生成可视化 ===\n")
pdf("14_all_subtypes_overview.pdf", width = 20, height = 16)

# 全局UMAP by celltype_fine
p1 <- DimPlot(sce.all, reduction = "umap", group.by = "celltype_fine",
              label = TRUE, repel = TRUE, label.size = 2, cols = mycols) +
  ggtitle("All cell subtypes") +
  theme(legend.position = "right", legend.text = element_text(size = 8))
print(p1)

# 全局UMAP by 大类
p2 <- DimPlot(sce.all, reduction = "umap", group.by = "celltype",
              label = TRUE, repel = TRUE, cols = mycols) +
  ggtitle("Major cell types")
print(p2)

# 样本分布
CellStatPlot(sce.all, stat.by = "orig.ident", group.by = "celltype_fine", plot_type = "bar")

# 各大类内部的亚群比例
CellStatPlot(sce.all, stat.by = "celltype_fine", group.by = "celltype", plot_type = "bar")

dev.off()

#==== 5. 保存最终对象 ====
cat("\n=== 保存最终对象 ===\n")
save(sce.all, file = "sce.all.annotated.Rdata")

# 也保存一个celltype_fine映射表
annotation_table <- sce.all@meta.data[, c("celltype", "celltype_fine", "orig.ident")]
annotation_table <- annotation_table %>%
  group_by(celltype, celltype_fine) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  arrange(celltype, desc(n_cells))
write.csv(annotation_table, "cell_annotation_summary.csv", row.names = FALSE)
cat("注释汇总表已保存: cell_annotation_summary.csv\n")

cat("\n=== 全部分析完成！===\n")
cat("最终输出:\n")
cat("  - sce.all.annotated.Rdata: 带有celltype_fine注释的完整对象\n")
cat("  - cell_annotation_summary.csv: 所有亚群的细胞数统计\n")
cat("  - 14_all_subtypes_overview.pdf: 整合可视化图\n")
