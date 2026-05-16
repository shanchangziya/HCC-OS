#==== Header ====
# 脚本名称: 11_CAF_subdivision.R
# 目的: 成纤维细胞(CAF)亚群细分（全自动注释）
# 日期: 2026-02-01
# 参考文献:
#   - J Translational Medicine 2025: Fibroblast subtypes in HCC
#   - Frontiers Oncology 2025: CAFs in liver cancer
#   - Cell Discovery 2023: CD36+ CAFs in HCC

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

#==== 1. 加载和提取CAF ====
cat("=== 加载数据 ===\n")
load("sceall.Rdata")
sce.caf <- subset(sce.all, subset = celltype == "fibo")
rm(sce.all); gc()
cat("CAF细胞数:", ncol(sce.caf), "\n")
print(table(sce.caf$orig.ident))

#==== 2. 标准化和降维 ====
cat("\n=== 标准化和降维 ===\n")
sce.caf <- NormalizeData(sce.caf, normalization.method = "LogNormalize", scale.factor = 1e4)
sce.caf <- FindVariableFeatures(sce.caf)
sce.caf <- ScaleData(sce.caf)
sce.caf <- RunPCA(sce.caf, features = VariableFeatures(object = sce.caf))
sce.caf <- RunHarmony(sce.caf, "orig.ident")

#==== 3. 聚类 ====
cat("\n=== 聚类 ===\n")
sce.caf <- FindNeighbors(sce.caf, reduction = "harmony", dims = 1:20)
sce.caf <- FindClusters(sce.caf, resolution = c(0.05, 0.1, 0.2, 0.3))
sce.caf <- RunUMAP(sce.caf, reduction = "harmony", dims = 1:20, seed.use = 42)

# 质量过滤
Idents(sce.caf) <- "RNA_snn_res.0.2"
tab <- table(Idents(sce.caf))
cat("\nCAF Resolution 0.2 簇分布:\n")
print(tab)
rm.clusters <- names(tab[tab < 50])
if(length(rm.clusters) > 0) {
  cat("删除小簇:", rm.clusters, "\n")
  sce.caf <- subset(sce.caf, idents = rm.clusters, invert = TRUE)
}

#==== 4. 自动亚型注释 ====
cat("\n=== CAF亚型注释 ===\n")

# CAF亚型marker（基于2025文献）
caf_markers_list <- list(
  "myCAF"        = c("ACTA2", "MYH11", "TAGLN", "RGS5", "ACTG2", "DES"),
  "iCAF"         = c("IL6", "CXCL12", "CCL2", "PDGFRA", "CFD", "EPHA3"),
  "apCAF"        = c("HLA-DRA", "HLA-DRB1", "CD74", "HLA-DRB5", "B2M", "ANTIGEN_PRESENTING"),
  "MMP_CAF"      = c("MMP11", "MMP2", "COL1A1", "COL1A2", "COL3A1", "LOXL2"),
  "VEGFA_CAF"    = c("VEGFA", "PDGFRB", "CXCL12", "ANGPT1", "NTN4", "FGF2")
)
# 修正apCAF中不合法的基因名
caf_markers_list[["apCAF"]] <- c("HLA-DRA", "HLA-DRB1", "CD74", "HLA-DRB5", "B2M", "HLA-DPB1")

sce.caf$CAF_subtype <- annotate_clusters(sce.caf, caf_markers_list, "RNA_snn_res.0.2")

cat("\nCAF亚型注释结果:\n")
print(table(sce.caf$CAF_subtype))
print(table(sce.caf$CAF_subtype, sce.caf$orig.ident))

#==== 5. 可视化 ====
pdf("11_CAF_subtyping.pdf", width = 14, height = 12)

clustree(sce.caf@meta.data, prefix = "RNA_snn_res.")
DimPlot(sce.caf, group.by = "CAF_subtype", label = TRUE, repel = TRUE, cols = mycols) +
  ggtitle("Cancer-associated fibroblast subtypes")
DimPlot(sce.caf, group.by = "RNA_snn_res.0.2", label = TRUE)

caf_show <- c("ACTA2", "MYH11", "TAGLN", "IL6", "CXCL12", "PDGFRA",
              "HLA-DRA", "HLA-DRB1", "CD74", "MMP11", "MMP2", "COL1A1", "VEGFA", "PDGFRB")
caf_show <- caf_show[caf_show %in% rownames(GetAssayData(sce.caf))]
DotPlot(sce.caf, features = caf_show, group.by = "CAF_subtype") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("CAF subtype markers")

VlnPlot(sce.caf, features = caf_show[1:min(8, length(caf_show))],
        group.by = "CAF_subtype", pt.size = 0, ncol = 4)
dev.off()

#==== 6. 保存 ====
save(sce.caf, file = "11_CAF_annotated.Rdata")
cat("\n=== CAF分析完成！保存: 11_CAF_annotated.Rdata ===\n")
