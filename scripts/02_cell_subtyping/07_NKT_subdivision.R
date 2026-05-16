#==== Header ====
# 脚本名称: 07_NKT_subdivision.R
# 目的: NK-T细胞分离与亚群细分（全自动注释）
# 日期: 2026-02-01
# 参考文献:
#   - Scientific Reports 2025: Single-cell landscape of CD8+ T cells in HCC
#   - Nature Communications 2022: NK cells in HCC
#   - Journal of Hepatology 2024: scRNA-seq signatures for HCC

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

# 自动注释函数：对每个cluster打分，分配给得分最高的亚型
annotate_clusters <- function(sce, marker_list, resolution_col, default_name = "Other") {
  # marker_list: 命名列表，每个元素是一种亚型的marker基因向量
  # 返回: 带有注释结果的字符向量

  exp_mat <- GetAssayData(sce, slot = "data")

  # 检查哪些marker在表达矩阵中存在
  all_markers <- unlist(marker_list)
  existing_markers <- all_markers[all_markers %in% rownames(exp_mat)]
  cat("使用的marker基因:", length(existing_markers), "/", length(all_markers), "\n")

  clusters <- as.character(sce[[resolution_col]][[1]])
  unique_clusters <- unique(clusters)

  # 对每个cluster计算各亚型得分
  score_mat <- matrix(0, nrow = length(unique_clusters), ncol = length(marker_list))
  rownames(score_mat) <- unique_clusters
  colnames(score_mat) <- names(marker_list)

  for(i in seq_along(marker_list)) {
    markers <- marker_list[[i]]
    markers <- markers[markers %in% rownames(exp_mat)]
    if(length(markers) == 0) next

    for(j in seq_along(unique_clusters)) {
      clust_cells <- which(clusters == unique_clusters[j])
      if(length(markers) == 1) {
        score_mat[j, i] <- mean(exp_mat[markers, clust_cells])
      } else {
        score_mat[j, i] <- mean(colMeans(exp_mat[markers, clust_cells, drop = FALSE]))
      }
    }
  }

  cat("\n各cluster得分矩阵:\n")
  print(round(score_mat, 3))

  # 分配注释：选择得分最高的亚型
  annotations <- character(length(unique_clusters))
  for(j in seq_along(unique_clusters)) {
    max_idx <- which.max(score_mat[j, ])
    annotations[j] <- names(marker_list)[max_idx]
  }
  names(annotations) <- unique_clusters

  cat("\nCluster → 亚型映射:\n")
  print(annotations)

  # 映射到细胞级别
  result <- annotations[clusters]

  return(result)
}

#==== 1. 加载数据 ====
cat("=== 加载数据 ===\n")
load("sceall.Rdata")
cat("总细胞数:", ncol(sce.all), "\n")
print(table(sce.all$celltype))

#==== 2. 提取NK-T细胞 ====
cat("\n=== 提取NK-T细胞 ===\n")
sce.nkt <- subset(sce.all, subset = celltype == "NK-T")
rm(sce.all); gc()
cat("NK-T细胞数:", ncol(sce.nkt), "\n")
print(table(sce.nkt$orig.ident))

#==== 3. 标准化和降维 ====
cat("\n=== 标准化和降维 ===\n")
sce.nkt <- NormalizeData(sce.nkt, normalization.method = "LogNormalize", scale.factor = 1e4)
sce.nkt <- FindVariableFeatures(sce.nkt)
sce.nkt <- ScaleData(sce.nkt)
sce.nkt <- RunPCA(sce.nkt, features = VariableFeatures(object = sce.nkt))
sce.nkt <- RunHarmony(sce.nkt, "orig.ident")

#==== 4. 初步聚类 (分离NK和T) ====
cat("\n=== 初步聚类识别NK vs T ===\n")
sce.nkt <- FindNeighbors(sce.nkt, reduction = "harmony", dims = 1:20)
sce.nkt <- FindClusters(sce.nkt, resolution = c(0.05, 0.1))
sce.nkt <- RunUMAP(sce.nkt, reduction = "harmony", dims = 1:20, seed.use = 42)

#==== 5. 分离NK和T细胞 ====
cat("\n=== 分离NK和T细胞 ===\n")
exp_mat <- GetAssayData(sce.nkt, slot = "data")
nk_markers_sep <- c("NCAM1", "KLRF1", "KLRD1")
t_markers_sep <- c("CD3D", "CD3E", "CD3G")

nk_markers_sep <- nk_markers_sep[nk_markers_sep %in% rownames(exp_mat)]
t_markers_sep <- t_markers_sep[t_markers_sep %in% rownames(exp_mat)]

sce.nkt$NK_score <- colMeans(exp_mat[nk_markers_sep, , drop = FALSE])
sce.nkt$T_score <- colMeans(exp_mat[t_markers_sep, , drop = FALSE])

# 基于cluster级别的平均分来决定：每个cluster归为NK还是T
Idents(sce.nkt) <- "RNA_snn_res.0.1"
cluster_nk_score <- tapply(sce.nkt$NK_score, Idents(sce.nkt), mean)
cluster_t_score <- tapply(sce.nkt$T_score, Idents(sce.nkt), mean)

cat("\nCluster NK/T 得分:\n")
score_df <- data.frame(
  Cluster = names(cluster_nk_score),
  NK_score = cluster_nk_score,
  T_score = cluster_t_score,
  Lineage = ifelse(cluster_nk_score > cluster_t_score, "NK", "T")
)
print(score_df)

# 映射到细胞
lineage_map <- setNames(score_df$Lineage, score_df$Cluster)
sce.nkt$cell_lineage <- lineage_map[as.character(Idents(sce.nkt))]

cat("\n分离结果:\n")
print(table(sce.nkt$cell_lineage))

# 可视化分离
pdf("07_NKT_overview.pdf", width = 16, height = 10)
p1 <- DimPlot(sce.nkt, group.by = "RNA_snn_res.0.1", label = TRUE) + ggtitle("Resolution 0.1")
p2 <- DimPlot(sce.nkt, group.by = "cell_lineage", cols = c("NK" = "#E15759", "T" = "#4E79A7")) + ggtitle("NK vs T")
p3 <- FeaturePlot(sce.nkt, features = "NK_score")
p4 <- FeaturePlot(sce.nkt, features = "T_score")
print(p1 | p2)
print(p3 | p4)
dev.off()

# 提取
sce.nk <- subset(sce.nkt, subset = cell_lineage == "NK")
sce.tc <- subset(sce.nkt, subset = cell_lineage == "T")
rm(sce.nkt); gc()

cat("NK细胞:", ncol(sce.nk), "\n")
cat("T细胞:", ncol(sce.tc), "\n")

#==============================================================
#==== 8. NK细胞亚群细分 ====
#==============================================================
cat("\n=== NK细胞亚群细分 ===\n")

# 重新分析
sce.nk <- NormalizeData(sce.nk, normalization.method = "LogNormalize", scale.factor = 1e4)
sce.nk <- FindVariableFeatures(sce.nk)
sce.nk <- ScaleData(sce.nk)
sce.nk <- RunPCA(sce.nk, features = VariableFeatures(object = sce.nk))
sce.nk <- RunHarmony(sce.nk, "orig.ident")
sce.nk <- FindNeighbors(sce.nk, reduction = "harmony", dims = 1:20)
sce.nk <- FindClusters(sce.nk, resolution = c(0.1, 0.2, 0.3))
sce.nk <- RunUMAP(sce.nk, reduction = "harmony", dims = 1:20, seed.use = 42)

# 质量过滤：删除<100细胞的簇
Idents(sce.nk) <- "RNA_snn_res.0.2"
tab <- table(Idents(sce.nk))
cat("\nNK Resolution 0.2 簇分布:\n")
print(tab)
rm.clusters <- names(tab[tab < 100])
if(length(rm.clusters) > 0) {
  cat("删除小簇:", rm.clusters, "\n")
  sce.nk <- subset(sce.nk, idents = rm.clusters, invert = TRUE)
}

# NK亚型marker定义（基于文献）
nk_markers_list <- list(
  "CD56dim_NK"    = c("FCGR3A", "GZMB", "PRF1", "FGFBP2", "NKG7"),
  "CD56bright_NK" = c("SELL", "CD27", "XCL1", "XCL2", "IL7R"),
  "Activated_NK"  = c("IFNG", "TNF", "CCL3", "CCL4", "GZMK"),
  "Exhausted_NK"  = c("PDCD1", "LAG3", "HAVCR2", "TIGIT", "ENTPD1")
)

# 自动注释
sce.nk$NK_subtype <- annotate_clusters(sce.nk, nk_markers_list, "RNA_snn_res.0.2")

cat("\nNK亚型注释结果:\n")
print(table(sce.nk$NK_subtype))
print(table(sce.nk$NK_subtype, sce.nk$orig.ident))

# 可视化
pdf("07_NK_subtyping.pdf", width = 14, height = 12)
clustree(sce.nk@meta.data, prefix = "RNA_snn_res.")
DimPlot(sce.nk, group.by = "NK_subtype", label = TRUE, repel = TRUE, cols = mycols) +
  ggtitle("NK cell subtypes")
DimPlot(sce.nk, group.by = "RNA_snn_res.0.2", label = TRUE)

nk_show_markers <- c("FCGR3A", "GZMB", "SELL", "XCL1", "IFNG", "TNF", "PDCD1", "LAG3")
nk_show_markers <- nk_show_markers[nk_show_markers %in% rownames(GetAssayData(sce.nk))]
DotPlot(sce.nk, features = nk_show_markers, group.by = "NK_subtype") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
VlnPlot(sce.nk, features = nk_show_markers, group.by = "NK_subtype", pt.size = 0, ncol = 4)
dev.off()

# 保存
save(sce.nk, file = "07_NK_annotated.Rdata")
cat("NK细胞分析完成！保存: 07_NK_annotated.Rdata\n")
rm(sce.nk); gc()

#==============================================================
#==== 9. T细胞亚群细分 ====
#==============================================================
cat("\n=== T细胞亚群细分 ===\n")

# 重新分析
sce.tc <- NormalizeData(sce.tc, normalization.method = "LogNormalize", scale.factor = 1e4)
sce.tc <- FindVariableFeatures(sce.tc)
sce.tc <- ScaleData(sce.tc)
sce.tc <- RunPCA(sce.tc, features = VariableFeatures(object = sce.tc))
sce.tc <- RunHarmony(sce.tc, "orig.ident")
sce.tc <- FindNeighbors(sce.tc, reduction = "harmony", dims = 1:20)
sce.tc <- FindClusters(sce.tc, resolution = c(0.1, 0.2, 0.3))
sce.tc <- RunUMAP(sce.tc, reduction = "harmony", dims = 1:20, seed.use = 42)

# 质量过滤
Idents(sce.tc) <- "RNA_snn_res.0.3"
tab <- table(Idents(sce.tc))
cat("\nT Resolution 0.3 簇分布:\n")
print(tab)
rm.clusters <- names(tab[tab < 100])
if(length(rm.clusters) > 0) {
  cat("删除小簇:", rm.clusters, "\n")
  sce.tc <- subset(sce.tc, idents = rm.clusters, invert = TRUE)
}

# 区分CD4/CD8 (cluster级别)
exp_mat_tc <- GetAssayData(sce.tc, slot = "data")
Idents(sce.tc) <- "RNA_snn_res.0.3"

cd4_genes <- c("CD4")[c("CD4") %in% rownames(exp_mat_tc)]
cd8_genes <- c("CD8A", "CD8B")[c("CD8A", "CD8B") %in% rownames(exp_mat_tc)]

if(length(cd4_genes) > 0) sce.tc$CD4_exp <- colMeans(exp_mat_tc[cd4_genes, , drop = FALSE]) else sce.tc$CD4_exp <- 0
if(length(cd8_genes) > 0) sce.tc$CD8_exp <- colMeans(exp_mat_tc[cd8_genes, , drop = FALSE]) else sce.tc$CD8_exp <- 0

clust_cd4 <- tapply(sce.tc$CD4_exp, Idents(sce.tc), mean)
clust_cd8 <- tapply(sce.tc$CD8_exp, Idents(sce.tc), mean)

cat("\nCluster CD4/CD8得分:\n")
cd_df <- data.frame(Cluster = names(clust_cd4), CD4 = clust_cd4, CD8 = clust_cd8,
                    Lineage = ifelse(clust_cd4 > clust_cd8, "CD4", "CD8"))
print(cd_df)

cd_map <- setNames(cd_df$Lineage, cd_df$Cluster)
sce.tc$T_lineage <- cd_map[as.character(Idents(sce.tc))]
cat("\nCD4/CD8分类:\n")
print(table(sce.tc$T_lineage))

# T细胞全局亚型marker（CD4和CD8共用框架）
t_markers_list <- list(
  "Naive_CD4"       = c("CCR7", "TCF7", "LEF1", "SELL", "CD4"),
  "Th1"             = c("IFNG", "TBX21", "STAT4", "CD4"),
  "Tfh"             = c("CXCL13", "PDCD1", "BCL6", "ICOS", "CD4"),
  "Treg"            = c("FOXP3", "CTLA4", "IL2RA", "IKZF2", "CD4"),
  "Cytotoxic_CD4"   = c("GZMK", "NKG7", "CD4"),
  "Naive_CD8"       = c("CCR7", "TCF7", "LEF1", "SELL", "CD8A"),
  "Effector_CD8"    = c("GZMK", "GZMH", "PRF1", "IFNG", "CD8A"),
  "Exhausted_CD8"   = c("PDCD1", "TIGIT", "LAG3", "HAVCR2", "TOX", "CD8A"),
  "Trm_CD8"         = c("CD69", "ITGAE", "CXCR6", "CD8A")
)

# 自动注释
sce.tc$T_subtype <- annotate_clusters(sce.tc, t_markers_list, "RNA_snn_res.0.3")

cat("\nT亚型注释结果:\n")
print(table(sce.tc$T_subtype))
print(table(sce.tc$T_subtype, sce.tc$T_lineage))
print(table(sce.tc$T_subtype, sce.tc$orig.ident))

# 可视化
pdf("07_T_subtyping.pdf", width = 16, height = 16)
clustree(sce.tc@meta.data, prefix = "RNA_snn_res.")
DimPlot(sce.tc, group.by = "T_lineage", cols = c("CD4" = "#F28E2B", "CD8" = "#4E79A7")) +
  ggtitle("CD4 vs CD8")
DimPlot(sce.tc, group.by = "T_subtype", label = TRUE, repel = TRUE, cols = mycols) +
  ggtitle("T cell subtypes")
DimPlot(sce.tc, group.by = "RNA_snn_res.0.3", label = TRUE)

t_show_markers <- c("CD4", "CD8A", "CCR7", "TCF7", "FOXP3", "CXCL13",
                     "GZMK", "GZMH", "PRF1", "PDCD1", "TIGIT", "LAG3",
                     "HAVCR2", "TOX", "CD69", "ITGAE", "TBX21", "IFNG")
t_show_markers <- t_show_markers[t_show_markers %in% rownames(GetAssayData(sce.tc))]

DotPlot(sce.tc, features = t_show_markers, group.by = "T_subtype") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("T cell subtype markers")

DotPlot(sce.tc, features = t_show_markers, group.by = "T_subtype",
        split.by = "T_lineage") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

VlnPlot(sce.tc, features = c("CCR7", "GZMK", "PDCD1", "FOXP3", "CXCL13", "TOX", "CD69"),
        group.by = "T_subtype", pt.size = 0, ncol = 4)
dev.off()

# 保存
save(sce.tc, file = "07_T_annotated.Rdata")
cat("T细胞分析完成！保存: 07_T_annotated.Rdata\n")

cat("\n=== 07 NK-T分析全部完成 ===\n")
