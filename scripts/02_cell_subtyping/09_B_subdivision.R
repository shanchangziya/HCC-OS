#==== Header ====
# 脚本名称: 09_B_subdivision.R
# 目的: B细胞亚群细分（全自动注释）
# 日期: 2026-02-01
# 参考文献:
#   - Cancer Cell International 2024: B cell landscape in HCC
#   - Frontiers in Oncology 2020: Plasma cells in HCC with cirrhosis

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

  result <- annotations[clusters]
  return(result)
}

#==== 1. 加载和提取B细胞 ====
cat("=== 加载数据 ===\n")
load("sceall.Rdata")
sce.b <- subset(sce.all, subset = celltype == "B")
rm(sce.all); gc()
cat("B细胞数:", ncol(sce.b), "\n")
print(table(sce.b$orig.ident))

#==== 2. 标准化和降维 ====
cat("\n=== 标准化和降维 ===\n")
sce.b <- NormalizeData(sce.b, normalization.method = "LogNormalize", scale.factor = 1e4)
sce.b <- FindVariableFeatures(sce.b)
sce.b <- ScaleData(sce.b)
sce.b <- RunPCA(sce.b, features = VariableFeatures(object = sce.b))
sce.b <- RunHarmony(sce.b, "orig.ident")

#==== 3. 聚类 ====
cat("\n=== 聚类 ===\n")
sce.b <- FindNeighbors(sce.b, reduction = "harmony", dims = 1:20)
sce.b <- FindClusters(sce.b, resolution = c(0.05, 0.1, 0.2))
sce.b <- RunUMAP(sce.b, reduction = "harmony", dims = 1:20, seed.use = 42)

# 质量过滤
Idents(sce.b) <- "RNA_snn_res.0.1"
tab <- table(Idents(sce.b))
cat("\nB Resolution 0.1 簇分布:\n")
print(tab)
rm.clusters <- names(tab[tab < 50])  # B细胞数较少，放宽至50
if(length(rm.clusters) > 0) {
  cat("删除小簇:", rm.clusters, "\n")
  sce.b <- subset(sce.b, idents = rm.clusters, invert = TRUE)
}

#==== 4. 自动亚型注释 ====
cat("\n=== B细胞亚型注释 ===\n")

# B细胞亚型marker（基于文献）
b_markers_list <- list(
  "Naive_B"         = c("TCL1A", "MS4A1", "IGHD", "IGHM", "NAIVE1", "KLF4"),
  "Memory_B"        = c("MS4A1", "IGHG1", "AIM2", "LMNA", "CD19", "CD27"),
  "Activated_B"     = c("EGR1", "JUNB", "NR4A2", "FOS", "JUN", "NR4A1"),
  "IgG_Plasma"      = c("IGHG1", "IGHG2", "IGHG3", "IGHG4", "MZB1", "PRDM1"),
  "IgA_Plasma"      = c("IGHA1", "IGHA2", "MZB1", "PRDM1", "SDC1", "BLIMP1"),
  "Prolif_Plasma"   = c("MKI67", "STMN1", "TOP2A", "IGHG1", "MZB1", "PCNA")
)

sce.b$B_subtype <- annotate_clusters(sce.b, b_markers_list, "RNA_snn_res.0.1")

cat("\nB亚型注释结果:\n")
print(table(sce.b$B_subtype))
print(table(sce.b$B_subtype, sce.b$orig.ident))

#==== 5. 可视化 ====
pdf("09_B_subtyping.pdf", width = 14, height = 12)

clustree(sce.b@meta.data, prefix = "RNA_snn_res.")
DimPlot(sce.b, group.by = "B_subtype", label = TRUE, repel = TRUE, cols = mycols) +
  ggtitle("B cell subtypes")
DimPlot(sce.b, group.by = "RNA_snn_res.0.1", label = TRUE)

b_show <- c("MS4A1", "TCL1A", "IGHD", "IGHM", "IGHG1", "IGHA1",
            "MZB1", "PRDM1", "EGR1", "JUNB", "MKI67", "AIM2", "LMNA")
b_show <- b_show[b_show %in% rownames(GetAssayData(sce.b))]
DotPlot(sce.b, features = b_show, group.by = "B_subtype") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("B cell subtype markers")

VlnPlot(sce.b, features = b_show[1:min(8, length(b_show))],
        group.by = "B_subtype", pt.size = 0, ncol = 4)
dev.off()

#==== 6. 保存 ====
save(sce.b, file = "09_B_annotated.Rdata")
cat("\n=== B细胞分析完成！保存: 09_B_annotated.Rdata ===\n")
