library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(SCP)
load("19_sce.all.integrated.Rdata")
table(sce.all.integrated$final_celltype)
sce.all.integrated
#load("sceall.Rdata")
sce.all
colnames(sce.all)[1:5]
colnames(sce.all.integrated)[1:5]
# 1) 从 sce.all.integrated 的列名中提取“真实barcode”（取最后一个 "_" 后面的部分）
bar_int <- sub("^.*_", "", colnames(sce.all.integrated))  # e.g. Epithelial_AAAC...-1 -> AAAC...-1

# 2) 建立 integrated -> barcode 的映射（注意：同一barcode如果出现重复，会出问题；下面会检查）
if (any(duplicated(bar_int))) {
  dup <- unique(bar_int[duplicated(bar_int)])
  stop(paste0("integrated 去前缀后出现重复 barcode（示例）：", paste(head(dup, 10), collapse = ", ")))
}

# 3) 找交集
common_cells <- intersect(colnames(sce.all), bar_int)
length(common_cells)

# 4) 只保留交集细胞（在 sce.all 中）
sce.all <- subset(sce.all, cells = common_cells)

# 5) 把 integrated 的 final_celltype 映射到 sce.all
ct_map <- sce.all.integrated$final_celltype
names(ct_map) <- bar_int   # 用“去前缀后的barcode”做名字

sce.all$final_celltype <- unname(ct_map[colnames(sce.all)])

# 6) 检查
sum(is.na(sce.all$final_celltype))
table(sce.all$final_celltype, useNA = "ifany")
DimPlot(sce.all, reduction = "umap", group.by = "final_celltype", label = TRUE)
CellDimPlot(sce.all, group.by = "final_celltype",show_stat = F)
save(sce.all,file = "注释sceall.Rdata")
rm(list = ls())
dev.off()
gc()
