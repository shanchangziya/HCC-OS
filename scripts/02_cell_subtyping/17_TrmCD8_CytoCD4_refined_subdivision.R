#==== Header ====
# 脚本名称: 17_TrmCD8_CytoCD4_refined_subdivision.R
# 目的: Trm CD8和Cytotoxic CD4细胞精细化细分（使用差异表达和更高分辨率）
# 日期: 2026-02-02
# 策略:
#   1. 使用更高的分辨率
#   2. 基于差异表达基因进行注释
#   3. 使用相对得分策略而非绝对得分

#==== 环境设置 ====
rm(list = ls())
gc()
library(Seurat)
library(harmony)
library(clustree)
library(SCP)
library(dplyr)
library(ggplot2)
library(patchwork)

mycols <- c(
  "#4E79A7","#F28E2B","#E15759","#76B7B2","#59A14F",
  "#EDC948","#B07AA1","#FF9DA7","#9C755F","#BAB0AC",
  "#1F77B4","#FF7F0E","#2CA02C","#D62728","#9467BD",
  "#8C564B","#E377C2","#7F7F7F","#BCBD22","#17BECF"
)

# 改进的注释函数：使用相对得分和差异表达
annotate_clusters_refined <- function(sce, marker_list, resolution_col, method = "relative") {
  exp_mat <- GetAssayData(sce, slot = "data")

  # 检查marker基因
  all_markers <- unlist(marker_list)
  existing_markers <- all_markers[all_markers %in% rownames(exp_mat)]
  cat("使用的marker基因:", length(existing_markers), "/", length(all_markers), "\n")

  clusters <- as.character(sce[[resolution_col]][[1]])
  unique_clusters <- unique(clusters)

  # 计算得分矩阵
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

  cat("\n各cluster原始得分矩阵:\n")
  print(round(score_mat, 3))

  if(method == "relative") {
    # 标准化得分：每个cluster的得分除以该cluster所有亚型得分的总和
    score_mat_norm <- score_mat / rowSums(score_mat)
    cat("\n标准化后的得分矩阵:\n")
    print(round(score_mat_norm, 3))

    # 计算每个cluster与全局平均的差异
    global_mean <- colMeans(score_mat_norm)
    cat("\n全局平均得分:\n")
    print(round(global_mean, 3))

    # 相对富集得分
    enrichment_mat <- t(apply(score_mat_norm, 1, function(x) x - global_mean))
    cat("\n富集得分矩阵 (相对于全局平均):\n")
    print(round(enrichment_mat, 3))

    # 基于最高富集得分分配亚型
    annotations <- character(length(unique_clusters))
    for(j in seq_along(unique_clusters)) {
      max_idx <- which.max(enrichment_mat[j, ])
      # 如果富集得分太低,标记为混合型
      if(enrichment_mat[j, max_idx] < 0.02) {
        annotations[j] <- paste0("Mixed_", names(marker_list)[max_idx])
      } else {
        annotations[j] <- names(marker_list)[max_idx]
      }
    }
  } else {
    # 使用原始最高得分策略
    annotations <- character(length(unique_clusters))
    for(j in seq_along(unique_clusters)) {
      max_idx <- which.max(score_mat[j, ])
      annotations[j] <- names(marker_list)[max_idx]
    }
  }

  names(annotations) <- unique_clusters

  cat("\nCluster → 亚型映射:\n")
  print(annotations)

  # 映射到细胞级别
  result <- annotations[clusters]

  return(result)
}

#==== 1. 加载T细胞数据 ====
cat("=== 加载T细胞数据 ===\n")
load("07_T_annotated.Rdata")
cat("总T细胞数:", ncol(sce.tc), "\n")

#==============================================================
#==== 2. Trm CD8细胞精细化细分 ====
#==============================================================
cat("\n\n====================================\n")
cat("=== Trm CD8细胞精细化细分 ===\n")
cat("====================================\n")

# 提取Trm CD8细胞
sce.trm <- subset(sce.tc, subset = T_subtype == "Trm_CD8")
cat("\nTrm CD8细胞数:", ncol(sce.trm), "\n")

if(ncol(sce.trm) < 200) {
  cat("警告: Trm CD8细胞数量过少(<200), 跳过细分\n")
} else {
  # 重新标准化和降维
  cat("\n--- 标准化和降维 ---\n")
  sce.trm <- NormalizeData(sce.trm, normalization.method = "LogNormalize", scale.factor = 1e4)
  sce.trm <- FindVariableFeatures(sce.trm, selection.method = "vst", nfeatures = 2000)
  sce.trm <- ScaleData(sce.trm)
  sce.trm <- RunPCA(sce.trm, features = VariableFeatures(object = sce.trm), verbose = FALSE)

  # Harmony整合
  if(length(unique(sce.trm$orig.ident)) > 1) {
    sce.trm <- RunHarmony(sce.trm, "orig.ident", verbose = FALSE)
    reduction_use <- "harmony"
  } else {
    reduction_use <- "pca"
  }

  # 使用更高的分辨率
  cat("\n--- 聚类分析（更高分辨率）---\n")
  sce.trm <- FindNeighbors(sce.trm, reduction = reduction_use, dims = 1:30)
  sce.trm <- FindClusters(sce.trm, resolution = c(0.3, 0.5, 0.7, 0.9, 1.0, 1.2))
  sce.trm <- RunUMAP(sce.trm, reduction = reduction_use, dims = 1:30, seed.use = 42, verbose = FALSE)

  # 尝试多个分辨率
  for(res in c(0.5, 0.7, 0.9)) {
    res_col <- paste0("RNA_snn_res.", res)
    Idents(sce.trm) <- res_col
    tab <- table(Idents(sce.trm))
    cat("\nResolution", res, "簇分布:\n")
    print(tab)
    cat("簇数量:", length(tab), "\n")
  }

  # 选择res=0.7作为基准
  Idents(sce.trm) <- "RNA_snn_res.0.7"
  tab <- table(Idents(sce.trm))

  # 删除小簇(<100细胞)
  rm.clusters <- names(tab[tab < 100])
  if(length(rm.clusters) > 0) {
    cat("\n删除小簇:", paste(rm.clusters, collapse = ", "), "\n")
    sce.trm <- subset(sce.trm, idents = rm.clusters, invert = TRUE)
    Idents(sce.trm) <- "RNA_snn_res.0.7"
  }

  # 优化的Trm CD8亚型marker定义
  trm_markers_list <- list(
    "Precursor_Trm"     = c("TCF7", "LEF1", "SELL", "IL7R"),
    "Mature_Trm"        = c("ITGAE", "CXCR6", "ITGA1", "ZNF683"),
    "Cytotoxic_Trm"     = c("GZMB", "PRF1", "IFNG", "NKG7"),
    "Activated_Trm"     = c("CD69", "TNFRSF9", "CD38", "HLA-DRA"),
    "Exhausted_Trm"     = c("PDCD1", "TIGIT", "LAG3", "HAVCR2", "TOX"),
    "Proliferating_Trm" = c("MKI67", "TOP2A", "STMN1")
  )

  # 使用改进的注释方法
  cat("\n--- 使用相对得分注释Trm亚型 ---\n")
  sce.trm$Trm_subtype_refined <- annotate_clusters_refined(
    sce.trm, trm_markers_list, "RNA_snn_res.0.7", method = "relative"
  )

  cat("\n\nTrm精细化亚型注释结果:\n")
  print(table(sce.trm$Trm_subtype_refined))
  cat("\nTrm亚型样本分布:\n")
  print(table(sce.trm$Trm_subtype_refined, sce.trm$orig.ident))

  # 运行差异表达分析来验证
  cat("\n--- 差异表达分析验证 ---\n")
  Idents(sce.trm) <- "Trm_subtype_refined"
  if(length(unique(Idents(sce.trm))) > 1) {
    trm.markers <- FindAllMarkers(sce.trm,
                                   only.pos = TRUE,
                                   min.pct = 0.25,
                                   logfc.threshold = 0.25,
                                   verbose = FALSE)

    cat("\n各亚型Top5 marker基因:\n")
    top_markers <- trm.markers %>%
      group_by(cluster) %>%
      top_n(n = 5, wt = avg_log2FC)
    print(top_markers)

    # 保存marker基因
    write.csv(trm.markers, "17_TrmCD8_DEG_markers.csv", row.names = FALSE)
  }

  # 可视化
  cat("\n--- 生成可视化 ---\n")
  pdf("17_TrmCD8_refined_subtyping.pdf", width = 20, height = 18)

  # Clustree
  print(clustree(sce.trm@meta.data, prefix = "RNA_snn_res."))

  # UMAP - 多分辨率比较
  p1 <- DimPlot(sce.trm, group.by = "RNA_snn_res.0.5", label = TRUE) +
    ggtitle("Resolution 0.5") + NoLegend()
  p2 <- DimPlot(sce.trm, group.by = "RNA_snn_res.0.7", label = TRUE) +
    ggtitle("Resolution 0.7") + NoLegend()
  p3 <- DimPlot(sce.trm, group.by = "RNA_snn_res.0.9", label = TRUE) +
    ggtitle("Resolution 0.9") + NoLegend()
  print(p1 | p2 | p3)

  # 亚型UMAP
  print(DimPlot(sce.trm, group.by = "Trm_subtype_refined",
                label = TRUE, repel = TRUE, cols = mycols) +
          ggtitle("Trm CD8 refined subtypes (Resolution 0.7)"))

  print(DimPlot(sce.trm, group.by = "orig.ident") + ggtitle("Samples"))

  # Marker基因展示
  trm_show_markers <- c("CD8A", "CD69", "ITGAE", "CXCR6", "TCF7", "LEF1", "IL7R",
                        "GZMB", "PRF1", "IFNG", "NKG7", "MKI67", "TOP2A",
                        "PDCD1", "TIGIT", "LAG3", "HAVCR2", "TOX",
                        "TNFRSF9", "CD38", "ZNF683", "ITGA1")
  trm_show_markers <- trm_show_markers[trm_show_markers %in% rownames(GetAssayData(sce.trm))]

  if(length(unique(Idents(sce.trm))) > 1) {
    print(DotPlot(sce.trm, features = trm_show_markers, group.by = "Trm_subtype_refined") +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            ggtitle("Trm CD8 refined subtype markers"))

    # VlnPlot关键marker
    key_markers <- c("CD69", "ITGAE", "TCF7", "GZMB", "PDCD1", "MKI67")
    key_markers <- key_markers[key_markers %in% rownames(GetAssayData(sce.trm))]
    print(VlnPlot(sce.trm, features = key_markers, group.by = "Trm_subtype_refined",
                  pt.size = 0, ncol = 3))

    # Heatmap - Top markers
    if(exists("top_markers") && nrow(top_markers) > 0) {
      top10 <- trm.markers %>%
        group_by(cluster) %>%
        top_n(n = 10, wt = avg_log2FC)
      print(DoHeatmap(sce.trm, features = unique(top10$gene), size = 3) +
              NoLegend() +
              theme(axis.text.y = element_text(size = 6)))
    }
  }

  # 样本间比较
  if(length(unique(sce.trm$orig.ident)) > 1 && length(unique(sce.trm$Trm_subtype_refined)) > 1) {
    prop_df <- as.data.frame.matrix(prop.table(table(sce.trm$Trm_subtype_refined, sce.trm$orig.ident), margin = 2))
    prop_df$Subtype <- rownames(prop_df)
    prop_df_long <- reshape2::melt(prop_df, id.vars = "Subtype", variable.name = "Sample", value.name = "Proportion")

    p_prop <- ggplot(prop_df_long, aes(x = Sample, y = Proportion, fill = Subtype)) +
      geom_bar(stat = "identity", position = "fill") +
      scale_fill_manual(values = mycols) +
      theme_classic() +
      ggtitle("Trm CD8 refined subtype proportions across samples") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    print(p_prop)
  }

  dev.off()

  # 保存
  save(sce.trm, file = "17_TrmCD8_refined_annotated.Rdata")
  cat("\nTrm CD8精细化细分完成！保存: 17_TrmCD8_refined_annotated.Rdata\n")
}

rm(sce.trm); gc()

#==============================================================
#==== 3. Cytotoxic CD4细胞精细化细分 ====
#==============================================================
cat("\n\n====================================\n")
cat("=== Cytotoxic CD4细胞精细化细分 ===\n")
cat("====================================\n")

# 提取Cytotoxic CD4细胞
sce.cyto <- subset(sce.tc, subset = T_subtype == "Cytotoxic_CD4")
cat("\nCytotoxic CD4细胞数:", ncol(sce.cyto), "\n")

if(ncol(sce.cyto) < 200) {
  cat("警告: Cytotoxic CD4细胞数量过少(<200), 跳过细分\n")
} else {
  # 重新标准化和降维
  cat("\n--- 标准化和降维 ---\n")
  sce.cyto <- NormalizeData(sce.cyto, normalization.method = "LogNormalize", scale.factor = 1e4)
  sce.cyto <- FindVariableFeatures(sce.cyto, selection.method = "vst", nfeatures = 2000)
  sce.cyto <- ScaleData(sce.cyto)
  sce.cyto <- RunPCA(sce.cyto, features = VariableFeatures(object = sce.cyto), verbose = FALSE)

  # Harmony整合
  if(length(unique(sce.cyto$orig.ident)) > 1) {
    sce.cyto <- RunHarmony(sce.cyto, "orig.ident", verbose = FALSE)
    reduction_use <- "harmony"
  } else {
    reduction_use <- "pca"
  }

  # 使用更高的分辨率
  cat("\n--- 聚类分析（更高分辨率）---\n")
  sce.cyto <- FindNeighbors(sce.cyto, reduction = reduction_use, dims = 1:30)
  sce.cyto <- FindClusters(sce.cyto, resolution = c(0.3, 0.5, 0.7, 0.9, 1.0, 1.2))
  sce.cyto <- RunUMAP(sce.cyto, reduction = reduction_use, dims = 1:30, seed.use = 42, verbose = FALSE)

  # 尝试多个分辨率
  for(res in c(0.5, 0.7, 0.9)) {
    res_col <- paste0("RNA_snn_res.", res)
    Idents(sce.cyto) <- res_col
    tab <- table(Idents(sce.cyto))
    cat("\nResolution", res, "簇分布:\n")
    print(tab)
    cat("簇数量:", length(tab), "\n")
  }

  # 选择res=0.7
  Idents(sce.cyto) <- "RNA_snn_res.0.7"
  tab <- table(Idents(sce.cyto))

  # 删除小簇
  rm.clusters <- names(tab[tab < 100])
  if(length(rm.clusters) > 0) {
    cat("\n删除小簇:", paste(rm.clusters, collapse = ", "), "\n")
    sce.cyto <- subset(sce.cyto, idents = rm.clusters, invert = TRUE)
    Idents(sce.cyto) <- "RNA_snn_res.0.7"
  }

  # 优化的Cytotoxic CD4亚型marker定义
  cytocd4_markers_list <- list(
    "Th1_Cyto"          = c("TBX21", "IFNG", "STAT1", "CXCR3"),
    "CTL_like"          = c("GZMB", "PRF1", "GNLY", "NKG7"),
    "Effector_Cyto"     = c("GZMA", "GZMK", "CCL5", "CCL4"),
    "Memory_Cyto"       = c("IL7R", "CCR7", "TCF7", "SELL"),
    "Exhausted_Cyto"    = c("PDCD1", "TIGIT", "LAG3", "HAVCR2", "TOX"),
    "Proliferating_Cyto"= c("MKI67", "TOP2A", "STMN1")
  )

  # 使用改进的注释方法
  cat("\n--- 使用相对得分注释Cytotoxic CD4亚型 ---\n")
  sce.cyto$CytoCD4_subtype_refined <- annotate_clusters_refined(
    sce.cyto, cytocd4_markers_list, "RNA_snn_res.0.7", method = "relative"
  )

  cat("\n\nCytotoxic CD4精细化亚型注释结果:\n")
  print(table(sce.cyto$CytoCD4_subtype_refined))
  cat("\nCytotoxic CD4亚型样本分布:\n")
  print(table(sce.cyto$CytoCD4_subtype_refined, sce.cyto$orig.ident))

  # 运行差异表达分析
  cat("\n--- 差异表达分析验证 ---\n")
  Idents(sce.cyto) <- "CytoCD4_subtype_refined"
  if(length(unique(Idents(sce.cyto))) > 1) {
    cyto.markers <- FindAllMarkers(sce.cyto,
                                    only.pos = TRUE,
                                    min.pct = 0.25,
                                    logfc.threshold = 0.25,
                                    verbose = FALSE)

    cat("\n各亚型Top5 marker基因:\n")
    top_markers <- cyto.markers %>%
      group_by(cluster) %>%
      top_n(n = 5, wt = avg_log2FC)
    print(top_markers)

    # 保存marker基因
    write.csv(cyto.markers, "17_CytoCD4_DEG_markers.csv", row.names = FALSE)
  }

  # 可视化
  cat("\n--- 生成可视化 ---\n")
  pdf("17_CytoCD4_refined_subtyping.pdf", width = 20, height = 18)

  # Clustree
  print(clustree(sce.cyto@meta.data, prefix = "RNA_snn_res."))

  # UMAP - 多分辨率比较
  p1 <- DimPlot(sce.cyto, group.by = "RNA_snn_res.0.5", label = TRUE) +
    ggtitle("Resolution 0.5") + NoLegend()
  p2 <- DimPlot(sce.cyto, group.by = "RNA_snn_res.0.7", label = TRUE) +
    ggtitle("Resolution 0.7") + NoLegend()
  p3 <- DimPlot(sce.cyto, group.by = "RNA_snn_res.0.9", label = TRUE) +
    ggtitle("Resolution 0.9") + NoLegend()
  print(p1 | p2 | p3)

  # 亚型UMAP
  print(DimPlot(sce.cyto, group.by = "CytoCD4_subtype_refined",
                label = TRUE, repel = TRUE, cols = mycols) +
          ggtitle("Cytotoxic CD4 refined subtypes (Resolution 0.7)"))

  print(DimPlot(sce.cyto, group.by = "orig.ident") + ggtitle("Samples"))

  # Marker基因展示
  cyto_show_markers <- c("CD4", "GZMK", "GZMA", "GZMB", "PRF1", "NKG7", "GNLY",
                         "TBX21", "IFNG", "STAT1", "CXCR3",
                         "CCL5", "CCL4", "IL7R", "CCR7", "TCF7",
                         "PDCD1", "TIGIT", "LAG3", "HAVCR2", "TOX",
                         "MKI67", "TOP2A")
  cyto_show_markers <- cyto_show_markers[cyto_show_markers %in% rownames(GetAssayData(sce.cyto))]

  if(length(unique(Idents(sce.cyto))) > 1) {
    print(DotPlot(sce.cyto, features = cyto_show_markers, group.by = "CytoCD4_subtype_refined") +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            ggtitle("Cytotoxic CD4 refined subtype markers"))

    # VlnPlot
    key_markers <- c("CD4", "GZMB", "IFNG", "CCL5", "PDCD1", "MKI67")
    key_markers <- key_markers[key_markers %in% rownames(GetAssayData(sce.cyto))]
    print(VlnPlot(sce.cyto, features = key_markers, group.by = "CytoCD4_subtype_refined",
                  pt.size = 0, ncol = 3))

    # Heatmap
    if(exists("top_markers") && nrow(top_markers) > 0) {
      top10 <- cyto.markers %>%
        group_by(cluster) %>%
        top_n(n = 10, wt = avg_log2FC)
      print(DoHeatmap(sce.cyto, features = unique(top10$gene), size = 3) +
              NoLegend() +
              theme(axis.text.y = element_text(size = 6)))
    }
  }

  # 样本间比较
  if(length(unique(sce.cyto$orig.ident)) > 1 && length(unique(sce.cyto$CytoCD4_subtype_refined)) > 1) {
    prop_df <- as.data.frame.matrix(prop.table(table(sce.cyto$CytoCD4_subtype_refined, sce.cyto$orig.ident), margin = 2))
    prop_df$Subtype <- rownames(prop_df)
    prop_df_long <- reshape2::melt(prop_df, id.vars = "Subtype", variable.name = "Sample", value.name = "Proportion")

    p_prop <- ggplot(prop_df_long, aes(x = Sample, y = Proportion, fill = Subtype)) +
      geom_bar(stat = "identity", position = "fill") +
      scale_fill_manual(values = mycols) +
      theme_classic() +
      ggtitle("Cytotoxic CD4 refined subtype proportions across samples") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    print(p_prop)
  }

  dev.off()

  # 保存
  save(sce.cyto, file = "17_CytoCD4_refined_annotated.Rdata")
  cat("\nCytotoxic CD4精细化细分完成！保存: 17_CytoCD4_refined_annotated.Rdata\n")
}

rm(sce.cyto); gc()

#==== 4. 生成总结报告 ====
cat("\n\n====================================\n")
cat("=== 精细化分析总结 ===\n")
cat("====================================\n")

summary_info <- data.frame(
  Cell_Type = character(),
  Total_Cells = numeric(),
  N_Subtypes = numeric(),
  Subtypes = character(),
  Output_File = character(),
  stringsAsFactors = FALSE
)

if(file.exists("17_TrmCD8_refined_annotated.Rdata")) {
  load("17_TrmCD8_refined_annotated.Rdata")
  subtypes_list <- unique(sce.trm$Trm_subtype_refined)
  summary_info <- rbind(summary_info, data.frame(
    Cell_Type = "Trm_CD8",
    Total_Cells = ncol(sce.trm),
    N_Subtypes = length(subtypes_list),
    Subtypes = paste(subtypes_list, collapse = "; "),
    Output_File = "17_TrmCD8_refined_annotated.Rdata"
  ))
  rm(sce.trm)
}

if(file.exists("17_CytoCD4_refined_annotated.Rdata")) {
  load("17_CytoCD4_refined_annotated.Rdata")
  subtypes_list <- unique(sce.cyto$CytoCD4_subtype_refined)
  summary_info <- rbind(summary_info, data.frame(
    Cell_Type = "Cytotoxic_CD4",
    Total_Cells = ncol(sce.cyto),
    N_Subtypes = length(subtypes_list),
    Subtypes = paste(subtypes_list, collapse = "; "),
    Output_File = "17_CytoCD4_refined_annotated.Rdata"
  ))
  rm(sce.cyto)
}

cat("\n精细化细分结果汇总:\n")
print(summary_info)

write.csv(summary_info, "17_refined_subdivision_summary.csv", row.names = FALSE)

cat("\n=== 17 精细化细分全部完成 ===\n")
cat("\n输出文件:\n")
cat("  - 17_TrmCD8_refined_subtyping.pdf\n")
cat("  - 17_TrmCD8_refined_annotated.Rdata\n")
cat("  - 17_TrmCD8_DEG_markers.csv\n")
cat("  - 17_CytoCD4_refined_subtyping.pdf\n")
cat("  - 17_CytoCD4_refined_annotated.Rdata\n")
cat("  - 17_CytoCD4_DEG_markers.csv\n")
cat("  - 17_refined_subdivision_summary.csv\n")
