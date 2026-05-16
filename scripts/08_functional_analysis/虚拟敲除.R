rm(list = ls())
gc()
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(igraph)
  library(ggrepel)
  library(patchwork)
})
load("tcos.Rdata")
ast_ctrl <- tcos
DefaultAssay(ast_ctrl) <- "RNA"
ast_ctrl <- NormalizeData(ast_ctrl, verbose = FALSE)
ast_ctrl <- FindVariableFeatures(ast_ctrl, nfeatures = 5000, verbose = FALSE)

genes_use <- unique(c(VariableFeatures(ast_ctrl), "NQO1"))
mat <- GetAssayData(ast_ctrl, slot = "counts")[genes_use, , drop = FALSE]

if (!"NQO1" %in% rownames(mat)) stop("NQO1 不在当前对象的基因名中，请检查大小写/物种命名。")

# 转换为标准矩阵并确保是数值型
mat <- as.matrix(mat)
storage.mode(mat) <- "numeric"

# 检查矩阵维度
cat("Matrix dimensions:", dim(mat), "\n")
cat("Matrix class:", class(mat), "\n")

## ===== 2) 安装/加载 scTenifoldKnk =====
if (!requireNamespace("scTenifoldKnk", quietly = TRUE)) {
  install.packages("scTenifoldKnk")
}
library(scTenifoldKnk)

## ===== 3) 跑 virtual KO =====
set.seed(1)
res_knk <- scTenifoldKnk(
  countMatrix   = mat,
  gKO           = "NQO1",
  qc            = FALSE,   # 关闭QC，因为数据已在Seurat中QC过
  nc_nNet       = 20,      # 网络数：更稳（默认10）
  nc_nCells     = min(500, ncol(mat)), # 每次子采样细胞数（默认500）
  nc_nComp      = 10,      # PCA主成分数：增加到10（原来5太少）
  td_K          = 3,
  ma_nDim       = 3,       # manifold维度：增加到3（原来2可能导致维度不足）
  nCores        = max(1, parallel::detectCores() - 1)
)
save(res_knk, file = "NQO1res_knk.Rdata")

## ===== 4) 提取结果 =====
cat("\n=== Virtual Knockout Results ===\n")
cat("Number of differentially regulated genes:", nrow(res_knk), "\n")
if (nrow(res_knk) > 0) {
  cat("Top 10 genes affected by NQO1 knockout:\n")
  print(head(res_knk[order(abs(res_knk$FC), decreasing = TRUE), ], 10))

  # 保存完整结果
  write.csv(res_knk, file = "NQO1_knockout_results.csv", row.names = TRUE)
  cat("\nFull results saved to: NQO1_knockout_results.csv\n")
}
load("NQO1res_knk.Rdata")

colnames(res_knk[["diffRegulation"]])
res_knk[["diffRegulation"]][1:2,]
