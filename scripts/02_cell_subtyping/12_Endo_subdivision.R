#==== Header ====
# 脚本名称: 12_Endo_subdivision.R
# 目的: 内皮细胞亚群细分（全自动注释）
# 日期: 2026-02-01
# 参考文献:
#   - Discover Oncology 2025: Tumor endothelial cell signatures in HCC
#   - Modern Pathology: Novel endothelial markers in HCC
#   - PLOS One 2025: Endothelial markers in HCC prognosis

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

#==== 1. 加载和提取内皮细胞 ====
cat("=== 加载数据 ===\n")
load("sceall.Rdata")
sce.endo <- subset(sce.all, subset = celltype == "Endo")
rm(sce.all); gc()
cat("Endo细胞数:", ncol(sce.endo), "\n")
print(table(sce.endo$orig.ident))

#==== 2. 标准化和降维 ====
cat("\n=== 标准化和降维 ===\n")
sce.endo <- NormalizeData(sce.endo, normalization.method = "LogNormalize", scale.factor = 1e4)
sce.endo <- FindVariableFeatures(sce.endo)
sce.endo <- ScaleData(sce.endo)
sce.endo <- RunPCA(sce.endo, features = VariableFeatures(object = sce.endo))
sce.endo <- RunHarmony(sce.endo, "orig.ident")

#==== 3. 聚类 ====
cat("\n=== 聚类 ===\n")
sce.endo <- FindNeighbors(sce.endo, reduction = "harmony", dims = 1:20)
sce.endo <- FindClusters(sce.endo, resolution = c(0.05, 0.1, 0.2))
sce.endo <- RunUMAP(sce.endo, reduction = "harmony", dims = 1:20, seed.use = 42)

# 质量过滤
Idents(sce.endo) <- "RNA_snn_res.0.1"
tab <- table(Idents(sce.endo))
cat("\nEndo Resolution 0.1 簇分布:\n")
print(tab)
rm.clusters <- names(tab[tab < 50])
if(length(rm.clusters) > 0) {
  cat("删除小簇:", rm.clusters, "\n")
  sce.endo <- subset(sce.endo, idents = rm.clusters, invert = TRUE)
}

#==== 4. 自动亚型注释 ====
cat("\n=== 内皮细胞亚型注释 ===\n")

# 内皮细胞亚型marker（基于2025文献）
endo_markers_list <- list(
  "VEC"              = c("PECAM1", "VWF", "KDR", "CDH5", "ICAM1", "VCAM1"),
  "LSEC"            = c("CLEC4G", "CLEC4M", "FCN2", "FCN3", "STAB2", "SPARCL1"),
  "LEC"              = c("PROX1", "LYVE1", "FLT4", "PDPN", "VEGFC", "CCBE1"),
  "Tip_Cell"         = c("ANGPT2", "ESM1", "APLN", "CXCL12", "DLL4", "VEGFA"),
  "Tumor_Endo"       = c("PODXL", "RGS5", "VEGFA", "NOTCH4", "COL4A1", "MCAM")
)

sce.endo$Endo_subtype <- annotate_clusters(sce.endo, endo_markers_list, "RNA_snn_res.0.1")

cat("\n内皮细胞亚型注释结果:\n")
print(table(sce.endo$Endo_subtype))
print(table(sce.endo$Endo_subtype, sce.endo$orig.ident))

#==== 5. 可视化 ====
pdf("12_Endo_subtyping.pdf", width = 14, height = 10)

clustree(sce.endo@meta.data, prefix = "RNA_snn_res.")
DimPlot(sce.endo, group.by = "Endo_subtype", label = TRUE, repel = TRUE, cols = mycols) +
  ggtitle("Endothelial cell subtypes")
DimPlot(sce.endo, group.by = "RNA_snn_res.0.1", label = TRUE)

endo_show <- c("PECAM1", "VWF", "KDR", "CLEC4G", "CLEC4M", "FCN3",
               "PROX1", "LYVE1", "FLT4", "ANGPT2", "ESM1", "PODXL", "RGS5")
endo_show <- endo_show[endo_show %in% rownames(GetAssayData(sce.endo))]
DotPlot(sce.endo, features = endo_show, group.by = "Endo_subtype") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Endothelial subtype markers")

VlnPlot(sce.endo, features = endo_show[1:min(8, length(endo_show))],
        group.by = "Endo_subtype", pt.size = 0, ncol = 4)
dev.off()

#==== 6. 保存 ====
save(sce.endo, file = "12_Endo_annotated.Rdata")
cat("\n=== 内皮细胞分析完成！保存: 12_Endo_annotated.Rdata ===\n")
