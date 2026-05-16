#==== Header ====
# 脚本名称: 13_Neu_subdivision.R
# 目的: 中性粒细胞/髓系细胞亚群细分（全自动注释）
# 日期: 2026-02-01
# 参考文献:
#   - Nature 2025: Liver tumour immune subtypes and neutrophil heterogeneity
#   - Scientific Reports 2025: NET and immune genes in HCC prognosis
#   - PMC 2023: Myeloid marker genes in HCC

#==== 环境设置 ====
rm(list = ls())
gc()
library(Seurat)
library(harmony)
library(clustree)
library(SCP)
library(dplyr)

mycols <- c(
  "#4E79A7","#F28E2B","#E15759","#76B7B2","#59A14F",
  "#EDC948","#B07AA1","#FF9DA7","#9C755F","#BAB0AC",
  "#1F77B4","#FF7F0E","#2CA02C","#D62728","#9467BD",
  "#8C564B","#E377C2","#7F7F7F","#BCBD22","#17BECF"
)

annotate_clusters <- function(sce, marker_list, resolution_col) {
  exp_mat <- GetAssayData(sce, slot = "data")
  clusters <- as.character(sce[[resolution_col]][[1]])
  unique_clusters <- unique(clusters)

  score_mat <- matrix(0, nrow = length(unique_clusters), ncol = length(marker_list))
  rownames(score_mat) <- unique_clusters
  colnames(score_mat) <- names(marker_list)

  for(i in seq_along(marker_list)) {
    markers <- marker_list[[i]][marker_list[[i]] %in% rownames(exp_mat)]
    if(length(markers) == 0) next
    for(j in seq_along(unique_clusters)) {
      cells <- which(clusters == unique_clusters[j])
      score_mat[j, i] <- mean(colMeans(exp_mat[markers, cells, drop = FALSE]))
    }
  }

  cat("各cluster得分矩阵:\n")
  print(round(score_mat, 3))

  annotations <- apply(score_mat, 1, function(x) names(which.max(x)))
  cat("Cluster → 亚型映射:\n")
  print(annotations)

  return(annotations[clusters])
}

#==== 1. 加载和提取中性粒细胞 ====
cat("=== 加载数据 ===\n")
load("sceall.Rdata")
sce.neu <- subset(sce.all, subset = celltype == "Neu")
rm(sce.all); gc()
cat("Neu细胞数:", ncol(sce.neu), "\n")
print(table(sce.neu$orig.ident))

#==== 2. 标准化和降维 ====
cat("\n=== 标准化和降维 ===\n")
sce.neu <- NormalizeData(sce.neu, normalization.method = "LogNormalize", scale.factor = 1e4)
sce.neu <- FindVariableFeatures(sce.neu)
sce.neu <- ScaleData(sce.neu)
sce.neu <- RunPCA(sce.neu, features = VariableFeatures(object = sce.neu))
sce.neu <- RunHarmony(sce.neu, "orig.ident")

#==== 3. 聚类 ====
cat("\n=== 聚类 ===\n")
sce.neu <- FindNeighbors(sce.neu, reduction = "harmony", dims = 1:20)
sce.neu <- FindClusters(sce.neu, resolution = c(0.05, 0.1, 0.2))
sce.neu <- RunUMAP(sce.neu, reduction = "harmony", dims = 1:20, seed.use = 42)

# 质量过滤
Idents(sce.neu) <- "RNA_snn_res.0.1"
tab <- table(Idents(sce.neu))
cat("\nNeu Resolution 0.1 簇分布:\n")
print(tab)
rm.clusters <- names(tab[tab < 50])
if(length(rm.clusters) > 0) {
  cat("删除小簇:", rm.clusters, "\n")
  sce.neu <- subset(sce.neu, idents = rm.clusters, invert = TRUE)
}

#==== 4. 自动亚型注释 ====
cat("\n=== 中性粒细胞亚型注释 ===\n")

# 髓系/中性粒细胞亚型marker（基于文献）
neu_markers_list <- list(
  "Neutrophil"         = c("FCGR3B", "CSF3R", "CXCR2", "FPR1", "ELANE", "MPO"),
  "PMN_MDSC"          = c("S100A8", "S100A9", "S100A12", "CXCR2", "VCAN", "CEACAM8"),
  "M_MDSC"            = c("CD14", "CCR2", "S100A12", "VCAN", "LYZ", "CSF1R"),
  "TAN"               = c("VEGFA", "MMP9", "ARG1", "CXCL8", "IL8", "TIGIT"),
  "Inflam_Neutrophil"  = c("TNF", "IL1B", "CXCL8", "CCL20", "IL6", "PTGS2")
)

sce.neu$Neu_subtype <- annotate_clusters(sce.neu, neu_markers_list, "RNA_snn_res.0.1")

cat("\n中性粒细胞亚型注释结果:\n")
print(table(sce.neu$Neu_subtype))
print(table(sce.neu$Neu_subtype, sce.neu$orig.ident))

#==== 5. 可视化 ====
pdf("13_Neu_subtyping.pdf", width = 14, height = 10)

clustree(sce.neu@meta.data, prefix = "RNA_snn_res.")
DimPlot(sce.neu, group.by = "Neu_subtype", label = TRUE, repel = TRUE, cols = mycols) +
  ggtitle("Neutrophil / Myeloid subtypes")
DimPlot(sce.neu, group.by = "RNA_snn_res.0.1", label = TRUE)

neu_show <- c("FCGR3B", "CSF3R", "CXCR2", "S100A8", "S100A9", "S100A12",
              "CD14", "CCR2", "LYZ", "VEGFA", "MMP9", "ARG1", "TNF", "IL1B")
neu_show <- neu_show[neu_show %in% rownames(GetAssayData(sce.neu))]
DotPlot(sce.neu, features = neu_show, group.by = "Neu_subtype") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Neutrophil/Myeloid subtype markers")

VlnPlot(sce.neu, features = neu_show[1:min(8, length(neu_show))],
        group.by = "Neu_subtype", pt.size = 0, ncol = 4)
dev.off()

#==== 6. 保存 ====
save(sce.neu, file = "13_Neu_annotated.Rdata")
cat("\n=== 中性粒细胞分析完成！保存: 13_Neu_annotated.Rdata ===\n")
