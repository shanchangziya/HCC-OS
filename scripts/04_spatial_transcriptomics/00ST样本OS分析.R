library(Seurat)
library(ggplot2)
load("data/02fastcnvST.Rdata")
obj1
suppressPackageStartupMessages({
  library(Seurat)
  library(msigdbr)
  library(dplyr)
  library(GSVA)
})
sp <- obj1

# 1) 取 Hallmark ROS 基因集
os_genes <- msigdbr(species = "Homo sapiens", category = "H") %>%
  filter(gs_name == "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY") %>%
  pull(gene_symbol) %>% unique()

os_use_sp <- intersect(os_genes, rownames(sp))
stopifnot(length(os_use_sp) >= 10)  # 太少说明基因ID不匹配，需要先转SYMBOL

# 2) 准备表达矩阵（genes x spots）
# 建议用 "data"（log-normalized），更稳定
expr <- GetAssayData(sp, assay = DefaultAssay(sp), slot = "data")
# 3) ssGSEA 打分（每个spot一个分数）
# 注意：GSVA 需要矩阵是 genes x samples（samples=spots）
os_list <- list(OxStress = os_use_sp)
os_ssgsea <- gsva(
  expr = as.matrix(expr),
  gset.idx.list = os_list,
  method = "ssgsea",
  kcdf = "Gaussian",        # 适合连续值（log-normalized）
  abs.ranking = TRUE,
  verbose = FALSE
)
# 4) 写回 Seurat meta.data（列名自定）
sp$OS_ssGSEA <- as.numeric(os_ssgsea["OxStress", colnames(sp)])
# （可选）z-score 标准化，方便不同切片/样本比较
sp$OS_ssGSEA_z <- as.numeric(scale(sp$OS_ssGSEA))
# 5) 画空间图
SpatialFeaturePlot(sp, features = "OS_ssGSEA_z", alpha = c(0.2, 1), pt.size.factor = 4000)
rm(list = ls())
gc()
dev.off()
