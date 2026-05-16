#==== Header ====
# 脚本名称: 10_DC_subdivision.R
# 目的: 树突状细胞亚群细分（全自动注释）
# 日期: 2026-02-01
# 参考文献:
#   - MDPI Vaccines 2025: DC function in HCC immunotherapy
#   - Nature Immunology 2025: DC subsets in tumor tissues
#   - Cellular & Molecular Immunology 2025: HCC immune microenvironment

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

#==== 1. 加载和提取DC ====
cat("=== 加载数据 ===\n")
load("sceall.Rdata")
sce.dc <- subset(sce.all, subset = celltype == "DC")
rm(sce.all); gc()
cat("DC细胞数:", ncol(sce.dc), "\n")
print(table(sce.dc$orig.ident))

#==== 2. 标准化和降维 ====
cat("\n=== 标准化和降维 ===\n")
sce.dc <- NormalizeData(sce.dc, normalization.method = "LogNormalize", scale.factor = 1e4)
sce.dc <- FindVariableFeatures(sce.dc)
sce.dc <- ScaleData(sce.dc)
sce.dc <- RunPCA(sce.dc, features = VariableFeatures(object = sce.dc))
sce.dc <- RunHarmony(sce.dc, "orig.ident")

#==== 3. 聚类 ====
cat("\n=== 聚类 ===\n")
sce.dc <- FindNeighbors(sce.dc, reduction = "harmony", dims = 1:20)
sce.dc <- FindClusters(sce.dc, resolution = c(0.05, 0.1, 0.2))
sce.dc <- RunUMAP(sce.dc, reduction = "harmony", dims = 1:20, seed.use = 42)

# 质量过滤
Idents(sce.dc) <- "RNA_snn_res.0.1"
tab <- table(Idents(sce.dc))
cat("\nDC Resolution 0.1 簇分布:\n")
print(tab)
rm.clusters <- names(tab[tab < 50])
if(length(rm.clusters) > 0) {
  cat("删除小簇:", rm.clusters, "\n")
  sce.dc <- subset(sce.dc, idents = rm.clusters, invert = TRUE)
}

#==== 4. 自动亚型注释 ====
cat("\n=== DC亚型注释 ===\n")

# DC亚型marker（基于2025文献）
dc_markers_list <- list(
  "cDC1"  = c("CLEC9A", "XCR1", "BATF3", "IRF8", "PCLCE2", "DNGR1"),
  "cDC2"  = c("CD1C", "FCER1A", "CLEC10A", "SIRPA", "CLEC4A4", "CYP1B1"),
  "pDC"   = c("IL3RA", "CLEC4C", "TCF4", "IRF7", "TLR7", "TLR9"),
  "moDC"  = c("CD14", "FCGR1A", "CD209", "S100A8", "S100A9", "MAFB")
)

sce.dc$DC_subtype <- annotate_clusters(sce.dc, dc_markers_list, "RNA_snn_res.0.1")

cat("\nDC亚型注释结果:\n")
print(table(sce.dc$DC_subtype))
print(table(sce.dc$DC_subtype, sce.dc$orig.ident))

#==== 5. 可视化 ====
pdf("10_DC_subtyping.pdf", width = 14, height = 10)

clustree(sce.dc@meta.data, prefix = "RNA_snn_res.")
DimPlot(sce.dc, group.by = "DC_subtype", label = TRUE, repel = TRUE, cols = mycols) +
  ggtitle("Dendritic cell subtypes")
DimPlot(sce.dc, group.by = "RNA_snn_res.0.1", label = TRUE)

dc_show <- c("CLEC9A", "XCR1", "CD1C", "FCER1A", "CLEC10A",
             "IL3RA", "CLEC4C", "TCF4", "CD14", "CD209", "S100A8")
dc_show <- dc_show[dc_show %in% rownames(GetAssayData(sce.dc))]
DotPlot(sce.dc, features = dc_show, group.by = "DC_subtype") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("DC subtype markers")

VlnPlot(sce.dc, features = dc_show[1:min(8, length(dc_show))],
        group.by = "DC_subtype", pt.size = 0, ncol = 4)
dev.off()

#==== 6. 保存 ====
save(sce.dc, file = "10_DC_annotated.Rdata")
cat("\n=== DC分析完成！保存: 10_DC_annotated.Rdata ===\n")
