#==== Header ====
# 脚本名称: 08_Mac_subdivision.R
# 目的: 巨噬细胞亚群细分（全自动注释）
# 日期: 2026-02-01
# 参考文献:
#   - MDPI Genes 2025: TAM plasticity in HCC (9 macrophage subtypes)
#   - Frontiers Cell Dev Biol 2023: TME cell subpopulations in HCC
#   - Cell Discovery 2020: Global immune characterization of HCC

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

# 自动注释函数
annotate_clusters <- function(sce, marker_list, resolution_col) {
  exp_mat <- GetAssayData(sce, slot = "data")
  all_markers <- unlist(marker_list)
  existing <- all_markers[all_markers %in% rownames(exp_mat)]
  cat("使用marker:", length(existing), "/", length(all_markers), "\n")

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

  cat("\n各cluster得分矩阵:\n")
  print(round(score_mat, 3))

  annotations <- apply(score_mat, 1, function(x) names(which.max(x)))
  cat("\nCluster → 亚型映射:\n")
  print(annotations)

  result <- annotations[clusters]
  return(result)
}

#==== 1. 加载和提取巨噬细胞 ====
cat("=== 加载数据 ===\n")
load("sceall.Rdata")
cat("总细胞数:", ncol(sce.all), "\n")

cat("\n=== 提取巨噬细胞 ===\n")
sce.mac <- subset(sce.all, subset = celltype == "Mac")
rm(sce.all); gc()
cat("Mac细胞数:", ncol(sce.mac), "\n")
print(table(sce.mac$orig.ident))

#==== 2. 标准化和降维 ====
cat("\n=== 标准化和降维 ===\n")
sce.mac <- NormalizeData(sce.mac, normalization.method = "LogNormalize", scale.factor = 1e4)
sce.mac <- FindVariableFeatures(sce.mac)
sce.mac <- ScaleData(sce.mac)
sce.mac <- RunPCA(sce.mac, features = VariableFeatures(object = sce.mac))
sce.mac <- RunHarmony(sce.mac, "orig.ident")

#==== 3. 聚类 ====
cat("\n=== 聚类 ===\n")
sce.mac <- FindNeighbors(sce.mac, reduction = "harmony", dims = 1:20)
sce.mac <- FindClusters(sce.mac, resolution = c(0.1, 0.2, 0.3))
sce.mac <- RunUMAP(sce.mac, reduction = "harmony", dims = 1:20, seed.use = 42)

# 质量过滤
Idents(sce.mac) <- "RNA_snn_res.0.2"
tab <- table(Idents(sce.mac))
cat("\nMac Resolution 0.2 簇分布:\n")
print(tab)
rm.clusters <- names(tab[tab < 100])
if(length(rm.clusters) > 0) {
  cat("删除小簇:", rm.clusters, "\n")
  sce.mac <- subset(sce.mac, idents = rm.clusters, invert = TRUE)
}

#==== 4. 自动亚型注释 ====
cat("\n=== 巨噬细胞亚型注释 ===\n")

# 9种TAM亚型marker（基于2025文献）
mac_markers_list <- list(
  "SPP1_TAM"       = c("SPP1", "CD44", "CSTB", "APOE"),
  "CXCL9_TAM"      = c("CXCL9", "CXCL10", "CXCL11", "STAT1", "GBP1"),
  "SLC40A1_TAM"    = c("SLC40A1", "HMOX1", "FTH1", "FTMT", "TXNIP"),
  "TREM2_TAM"      = c("TREM2", "APOE", "APOC1", "CD9", "PSAP"),
  "HSP_TAM"        = c("HSPA1A", "HSPA1B", "DNAJB1", "HSPA1L", "HSPB1"),
  "Kupffer"        = c("MARCO", "TIMD4", "CD5L", "VCAM1", "MERTK"),
  "MDSC_like"      = c("S100A8", "S100A9", "S100A12", "VCAN", "CD14"),
  "Proliferating"  = c("STMN1", "MKI67", "TOP2A", "BIRC5", "PCNA")
)

sce.mac$Mac_subtype <- annotate_clusters(sce.mac, mac_markers_list, "RNA_snn_res.0.2")

cat("\nMac亚型注释结果:\n")
print(table(sce.mac$Mac_subtype))
print(table(sce.mac$Mac_subtype, sce.mac$orig.ident))

#==== 5. 可视化 ====
pdf("08_Mac_subtyping.pdf", width = 16, height = 14)

clustree(sce.mac@meta.data, prefix = "RNA_snn_res.")
DimPlot(sce.mac, group.by = "Mac_subtype", label = TRUE, repel = TRUE, cols = mycols) +
  ggtitle("Macrophage subtypes")
DimPlot(sce.mac, group.by = "RNA_snn_res.0.2", label = TRUE)

# Marker基因DotPlot
mac_show <- c("SPP1", "CD44", "CXCL9", "CXCL10", "SLC40A1", "HMOX1",
              "TREM2", "APOE", "HSPA1A", "MARCO", "TIMD4",
              "S100A8", "S100A9", "STMN1", "MKI67")
mac_show <- mac_show[mac_show %in% rownames(GetAssayData(sce.mac))]
DotPlot(sce.mac, features = mac_show, group.by = "Mac_subtype") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Macrophage subtype markers")

VlnPlot(sce.mac, features = mac_show[1:8], group.by = "Mac_subtype", pt.size = 0, ncol = 4)

# 样本比例
CellStatPlot(sce.mac, stat.by = "orig.ident", group.by = "Mac_subtype", plot_type = "bar")

dev.off()

#==== 6. 保存 ====
save(sce.mac, file = "08_Mac_annotated.Rdata")
cat("\n=== 巨噬细胞分析完成！保存: 08_Mac_annotated.Rdata ===\n")
