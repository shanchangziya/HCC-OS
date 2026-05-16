library(Seurat)
#remotes::install_github("must-bioinfo/fastCNV")
library(fastCNV)
#载入全部的对象
load("sceall.Rdata") #全部单细胞对象
table(sce.all$celltype)
load("02ep.Rdata") #上皮细胞对象
# 保存
#saveRDS(epi.obj, file = "sce_all_epi.rds")
ep[["RNA5"]] <- as(object = ep[["RNA"]], Class = "Assay5")
DefaultAssay(ep) <- 'RNA5'
table(ep$group)
sce2 <- ep
table(Idents(sce2))
table(sce2$RNA_snn_res.0.5)
DimPlot(sce2,group.by = "RNA_snn_res.0.5")
Idents(sce2) <- sce2$RNA_snn_res.0.5
#转为V5对象
sce.all[["RNA5"]] <- as(object = sce.all[["RNA"]], Class = "Assay5")
DefaultAssay(sce.all) <- 'RNA5'
table(sce.all$celltype)

## 2) 从 sce.all 抽取 T/NK 的 4000 个细胞
set.seed(123)
tnk_cells_all <- colnames(sce.all)[sce.all$celltype == "NK-T"]
stopifnot(length(tnk_cells_all) > 0)
tnk_cells_pick <- sample(tnk_cells_all, size = min(4000, length(tnk_cells_all)), replace = FALSE)
tnk_obj <- subset(sce.all, cells = tnk_cells_pick)
## 3) 确保两边都有 RNA5 assay，并设为默认
if (!"RNA5" %in% Assays(tnk_obj)) {
  tnk_obj[["RNA5"]] <- as(object = tnk_obj[["RNA"]], Class = "Assay5")
}
DefaultAssay(tnk_obj) <- "RNA5"

if (!"RNA5" %in% Assays(sce2)) {
  sce2[["RNA5"]] <- as(object = sce2[["RNA"]], Class = "Assay5")
}
DefaultAssay(sce2) <- "RNA5"

## 4) 合并到 sce2（加前缀避免细胞名重复）
sce2_merge <- merge(
  x = sce2,
  y = tnk_obj,
  add.cell.ids = c("Hep", "TNK"),
  merge.data = TRUE
)

## 5) 检查合并结果
table(sce2_merge$celltype)
dim(sce2_merge)

table(sce2_merge$RNA_snn_res.0.3)


# 1) 找出合并后哪些是原 sce2 的肝细胞（带 Hep_ 前缀）
hep_cells <- grep("^Hep_", colnames(sce2_merge), value = TRUE)
tnk_cells <- grep("^TNK_", colnames(sce2_merge), value = TRUE)  # 如果你当时用的是 "TNK"

# 2) 建一个“合并后细胞名 -> 原 sce2 的 cluster标签”的映射
#    注意：sce2 里原细胞名会变成 Hep_ + 原名
hep_map <- setNames(as.character(sce2$RNA_snn_res.0.5),
                    paste0("Hep_", colnames(sce2)))

# 3) 新建一个细粒度标签列（推荐保留原 celltype 以防后续要用）
sce2_merge$celltype_fine <- NA_character_
sce2_merge$celltype_fine[hep_cells] <- hep_map[hep_cells]
sce2_merge$celltype_fine[tnk_cells] <- "NK-T"

# 4) 检查
table(sce2_merge$celltype_fine, useNA = "ifany")

sampleNames <- "Hepatocyte"
referencelabels <- c("NK-T")
sce2_merge <- fastCNV(sce2_merge, sampleNames, referenceVar = "celltype_fine", 
                      referenceLabel = referencelabels, printPlot = F,getCNVClusters = F,
                      outputType = "pdf",
                      savePath = ".")
library(ggplot2)

common_theme <- theme(
  plot.title = element_text(size = 10),
  legend.text = element_text(size = 8),
  legend.title = element_text(size = 8),
  axis.title = element_text(size = 8),
  axis.text = element_text(size = 6)
)
# 1) 用 T/NK 定阈值
ref <- subset(sce2_merge, subset = celltype_fine == "NK-T")@meta.data$cnv_fraction
T99 <- as.numeric(quantile(ref, 0.99, na.rm = TRUE))
T95 <- as.numeric(quantile(ref, 0.95, na.rm = TRUE))

# 2) 只对肝细胞打标签（避免把免疫细胞误分到肿瘤）
sce2_merge$cnv_class <- NA_character_
is_hep <- sce2_merge$celltype == "Hepatocyte" | sce2_merge$celltype_fine %in% as.character(0:10)

sce2_merge$cnv_class[is_hep & sce2_merge$cnv_fraction <= T95] <- "CNV-low"
sce2_merge$cnv_class[is_hep & sce2_merge$cnv_fraction >  T95 & sce2_merge$cnv_fraction <= T99] <- "CNV-mid"
sce2_merge$cnv_class[is_hep & sce2_merge$cnv_fraction >  T99] <- "CNV-high(malignant)"

# 免疫细胞标注为 reference
sce2_merge$cnv_class[!is_hep] <- "reference"

# 3) 看每个 hepatocyte cluster 的 malignant 占比
prop.table(table(sce2_merge$celltype_fine[is_hep], sce2_merge$cnv_class[is_hep]), 1)
table(sce2$RNA_snn_res.0.5)
DimPlot(sce2,group.by = "RNA_snn_res.0.5",label = T)
sce2_merge <- CNVClassification(sce2_merge)
table(sce2_merge$cnv_class)
head(rownames(sce2_merge@meta.data))

## 1) 找出 sce2_merge 中所有含 cnv 的 metadata 列
cnv_cols <- grep("cnv", colnames(sce2_merge@meta.data), ignore.case = TRUE, value = TRUE)
cnv_cols
table(cnv_cols)
# 例如可能是：cnv_fraction, cnv_class, cnv_leiden, cnv_score...

## 2) 取出这些列，并保留行名（细胞名）
md_cnv <- sce2_merge@meta.data[, cnv_cols, drop = FALSE]
md_cnv$cell_merge <- rownames(md_cnv)
## 3) 只保留 Hep_ 前缀的细胞（也就是原 sce2 的肝细胞部分）
hep_mask <- grepl("^Hep_", md_cnv$cell_merge)
md_cnv_hep <- md_cnv[hep_mask, , drop = FALSE]
## 4) 去掉 Hep_ 前缀，得到原 sce2 的细胞名
md_cnv_hep$cell_sce2 <- sub("^Hep_", "", md_cnv_hep$cell_merge)
## 5) 以 sce2 细胞名为行名，方便直接写回 meta.data
rownames(md_cnv_hep) <- md_cnv_hep$cell_sce2

## 6) 对齐到 sce2 的细胞顺序，并写回 sce2@meta.data
md_cnv_hep <- md_cnv_hep[colnames(sce2), cnv_cols, drop = FALSE]

## 如果你想避免覆盖 sce2 里同名列，可以先加前缀：
# colnames(md_cnv_hep) <- paste0("infercnv_", colnames(md_cnv_hep))
sce2 <- AddMetaData(sce2, metadata = md_cnv_hep)
## 7) 检查：sce2 现在是否有 cnv 列、是否有 NA
colnames(sce2@meta.data)[grep("cnv", colnames(sce2@meta.data), ignore.case=TRUE)]
colSums(is.na(sce2@meta.data[, cnv_cols, drop=FALSE]))
FeaturePlot(sce2,features = "cnv_fraction")
epcnv <- sce2
save(epcnv,file = "03epcnv.Rdata")
#save(sce2_merge,file = "data/scrna/HCC/fastcnv.Rdata")
# 可视化 ---------------------------------------------------------------------
library(ggplot2)
# 1) 取数据
df <- FetchData(sce2_merge, vars = c("celltype_fine", "cnv_fraction"))
df$celltype_fine <- as.character(df$celltype_fine)
# 2) 固定顺序：0-10 + T/NK
df$cluster <- factor(df$celltype_fine, levels = c(as.character(0:10), "NK-T"))
df$xlab <- ifelse(df$cluster == "NK-T", "Ref cell", paste0("Hep_", as.character(df$cluster)))
df$xlab <- factor(df$xlab, levels = c(paste0("Hep_", 0:10), "Ref cell"))

# 3) 设定分组：把 3 也归为“正常/ref”
df$type <- ifelse(df$cluster %in% c("3", "NK-T"), "ref", "hepatocyte")
# 4) 计算 T/NK 的 99th 阈值
thr <- as.numeric(quantile(df$cnv_fraction[df$cluster == "NK-T"], 0.99, na.rm = TRUE))
# 5) 
p <- ggplot(df, aes(x = xlab, y = cnv_fraction)) +
  geom_boxplot(aes(fill = type),
               width = 0.60,
               outlier.size = 0.35,
               outlier.alpha = 0.55,
               color = "black") +
  geom_jitter(aes(color = type),
              width = 0.18,
              size = 0.30,
              alpha = 0.18,
              show.legend = FALSE) +
  geom_hline(yintercept = thr, linetype = "dashed", linewidth = 0.6) +
  annotate("text",
           x = length(levels(df$xlab)) - 1, y = thr,
           label = paste0("T/NK 99th=", sprintf("%.3f", thr)),
           hjust = 0, vjust = -0.6, size = 3.2) +
  scale_fill_manual(values = c(hepatocyte = "#F8766D", ref = "#00BFC4"), guide = "none") +
  scale_color_manual(values = c(hepatocyte = "#F8766D", ref = "#00BFC4"), guide = "none") +
  labs(x = NULL, y = "cnv_fraction") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title.y = element_text(color = "black"),
    legend.position = "none"
  )

p
save(epcnv,file = "03epcnv.Rdata")

##将hep3群剔除，只留下肿瘤细胞
# 只去除 3 群，保留其余群 -> tc
Idents(epcnv) <- "RNA_snn_res.0.5"
tc <- subset(epcnv, idents = "3", invert = TRUE)
# 检查
table(Idents(tc))
save(tc,file = "03tc.Rdata")

rm(list = ls())
gc()
dev.off()
