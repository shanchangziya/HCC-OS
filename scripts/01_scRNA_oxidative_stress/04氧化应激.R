library(Seurat)
library(SCP)
# 给定标签 --------------------------------------------------------------------
load("02ep.Rdata")#所有的上皮细胞
DimPlot(ep,group.by = "RNA_snn_res.0.5")
table(ep$RNA_snn_res.0.5)
load("03epcnv.Rdata")
ep$malignant_status
FeatureDimPlot(epcnv,features = "cnv_fraction")
CellDimPlot(ep,group.by = "malignant_status")
# 确保分群是字符，避免因因子水平导致判断出错
ep$RNA_snn_res.0.5 <- as.character(ep$RNA_snn_res.0.5)

# 生成肿瘤/正常标签：3群=Normal，其余=Tumor
ep$malignant_status <- ifelse(ep$RNA_snn_res.0.5 == "3", "Normal", "Tumor")
ep$malignant_status <- factor(ep$malignant_status, levels = c("Normal", "Tumor"))

# 检查一下分配是否正确
table(ep$RNA_snn_res.0.5, ep$malignant_status)
table(ep$malignant_status,ep$group)

# 肿瘤细胞和正常细胞进行评分 -----------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
})

if (!requireNamespace("msigdbr", quietly = TRUE)) {
  message("msigdbr 未安装，使用手工氧化应激/抗氧化基因集作为 fallback。")
  os_genes <- c(
    "HMOX1","NQO1","SOD1","SOD2","CAT","GPX1","GPX3","PRDX1","PRDX2","PRDX4",
    "TXN","TXNRD1","GCLC","GCLM","GSR","FTH1","FTL","SQSTM1","HSPA1A","HSPA1B",
    "DDIT3","ATF4","JUN","FOS","MAFK"
  )
} else {
  suppressPackageStartupMessages({
    library(msigdbr)
    library(dplyr)
  })
  os_genes <- msigdbr(species = "Homo sapiens", category = "H") %>%
    dplyr::filter(gs_name == "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY") %>%
    dplyr::pull(gene_symbol) %>% unique()
}
sce.all <- ep
## 基因名匹配到你的对象（非常关键）
os_genes_use <- intersect(os_genes, rownames(sce.all))
message("OxStress genes in object: ", length(os_genes_use), " / ", length(os_genes))


## ========= 1) 评分法一：Seurat AddModuleScore（module score） =========
sce.all <- AddModuleScore(
  object = sce.all,
  features = list(OxStress = os_genes_use),
  name = "OS_AddModuleScore_",
  assay = "RNA"
)
## 生成的列名通常为：OS_AddModuleScore_1
VlnPlot(sce.all,group.by = "malignant_status",features ="OS_AddModuleScore_1" )

## ========= 2) 评分法二：UCell（rank-based，稳健，强烈推荐） =========
if (!requireNamespace("UCell", quietly = TRUE)) {
  stop("请先安装 UCell：install.packages('UCell') 或 remotes::install_github('carmonalab/UCell')")
}
sce.all <- UCell::AddModuleScore_UCell(
  sce.all,
  features = list(OxStress = os_genes_use),
  name = "OS_UCell_",
  assay = "RNA"
)
## 常见生成列名：OxStressOS_UCell_
sce.all$OxStressOS_UCell_

## ========= 3) 评分法三：AUCell（AUC enrichment） =========
if (!requireNamespace("AUCell", quietly = TRUE)) {
  stop("请先安装 AUCell：BiocManager::install('AUCell')")
}
expr_mat <- GetAssayData(sce.all, assay = "RNA", slot = "data")  ## genes x cells（log-normalized）
rankings <- AUCell::AUCell_buildRankings(expr_mat, plotStats = FALSE, nCores = 20)
auc_obj  <- AUCell::AUCell_calcAUC(list(OxStress = os_genes_use), rankings, nCores = 20)
auc_vec  <- as.numeric(AUCell::getAUC(auc_obj)["OxStress", ])
names(auc_vec) <- colnames(sce.all)
sce.all <- AddMetaData(sce.all, metadata = auc_vec, col.name = "OS_AUCell_OxStress")


## ========= 4) 看哪个细胞群更高（按 seurat_clusters；也可换成你的 celltype 列） =========
## 小提琴图：三种评分并排看
VlnPlot(
  sce.all,
  features = c("OS_AddModuleScore_1", "OxStressOS_UCell_", "OS_AUCell_OxStress"),
  group.by = "malignant_status",
  pt.size = 0
)

## 计算每个cluster的平均分（快速表格）
score_df <- FetchData(
  sce.all,
  vars = c("malignant_status", "OS_AddModuleScore_1", "OxStressOS_UCell_", "OS_AUCell_OxStress")
)
cluster_mean <- aggregate(
  score_df[, c("OS_AddModuleScore_1","OxStressOS_UCell_","OS_AUCell_OxStress")],
  by = list(cluster = score_df$malignant_status),
  FUN = mean
)
cluster_mean[order(cluster_mean$OxStressOS_UCell_, decreasing = TRUE), ]

FeatureDimPlot(sce.all$malignant_status,
               features =c("OS_AddModuleScore_1","OxStressOS_UCell_","OS_AUCell_OxStress"))
features <- c("OS_AddModuleScore_1", "OxStressOS_UCell_", "OS_AUCell_OxStress")

# 确认列都存在
stopifnot(all(features %in% colnames(sce.all@meta.data)))

md <- sce.all@meta.data[, features, drop = FALSE]

# 1) 直接平均（raw mean）
sce.all$OS_score_mean_raw <- rowMeans(md, na.rm = TRUE)

# 2) 标准化后平均（推荐：避免不同评分尺度不一致）
md_z <- scale(md)  # 每列做 (x-mean)/sd
sce.all$OS_score_mean_z <- rowMeans(md_z, na.rm = TRUE)

# 快速检查
summary(sce.all$OS_score_mean_raw)
summary(sce.all$OS_score_mean_z)

# 可视化（按你的cluster列名自行替换 group.by）
VlnPlot(sce.all, features = c("OS_score_mean_raw", "OS_score_mean_z"), group.by = "malignant_status", pt.size = 0)
CellDimPlot(sce.all,group.by = "celltype",label = T)
FeatureStatPlot(sce.all,group.by = "malignant_status",stat.by = "OS_score_mean_z")

rm(list = ls())
gc()
# 准备单纯给肿瘤细胞评分 ---------------------------------------------------------------
load("03tc.Rdata")#只有肿瘤细胞
table(Idents(tc))
sce <- tc
sce <- NormalizeData(sce, 
                     normalization.method = "LogNormalize",
                     scale.factor = 1e4) 
sce <- FindVariableFeatures(sce)
sce <- ScaleData(sce)
sce <- RunPCA(sce, features = VariableFeatures(object = sce))
library(harmony)
sce <- RunHarmony(sce, "orig.ident")
# 统一：Neighbors / Clusters / UMAP 都用 harmony，同一组 dims
sce <- FindNeighbors(sce, reduction = "harmony", dims = 1:20)

CellDensityPlot(sce,features = "OS_score_mean_z")

# 肿瘤细胞和正常细胞进行评分 -----------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
})

if (!requireNamespace("msigdbr", quietly = TRUE)) {
  message("msigdbr 未安装，使用手工氧化应激/抗氧化基因集作为 fallback。")
  os_genes <- c(
    "HMOX1","NQO1","SOD1","SOD2","CAT","GPX1","GPX3","PRDX1","PRDX2","PRDX4",
    "TXN","TXNRD1","GCLC","GCLM","GSR","FTH1","FTL","SQSTM1","HSPA1A","HSPA1B",
    "DDIT3","ATF4","JUN","FOS","MAFK"
  )
} else {
  suppressPackageStartupMessages({
    library(msigdbr)
    library(dplyr)
  })
  os_genes <- msigdbr(species = "Homo sapiens", category = "H") %>%
    dplyr::filter(gs_name == "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY") %>%
    dplyr::pull(gene_symbol) %>% unique()
}
sce.all <- sce
## 基因名匹配到你的对象（非常关键）
os_genes_use <- intersect(os_genes, rownames(sce.all))
message("OxStress genes in object: ", length(os_genes_use), " / ", length(os_genes))


## ========= 1) 评分法一：Seurat AddModuleScore（module score） =========
sce.all <- AddModuleScore(
  object = sce.all,
  features = list(OxStress = os_genes_use),
  name = "OS_AddModuleScore_",
  assay = "RNA"
)
## 生成的列名通常为：OS_AddModuleScore_1


## ========= 2) 评分法二：UCell（rank-based，稳健，强烈推荐） =========
if (!requireNamespace("UCell", quietly = TRUE)) {
  stop("请先安装 UCell：install.packages('UCell') 或 remotes::install_github('carmonalab/UCell')")
}
sce.all <- UCell::AddModuleScore_UCell(
  sce.all,
  features = list(OxStress = os_genes_use),
  name = "OS_UCell_",
  assay = "RNA"
)
## 常见生成列名：OxStressOS_UCell_
sce.all$OxStressOS_UCell_

## ========= 3) 评分法三：AUCell（AUC enrichment） =========
if (!requireNamespace("AUCell", quietly = TRUE)) {
  stop("请先安装 AUCell：BiocManager::install('AUCell')")
}
expr_mat <- GetAssayData(sce.all, assay = "RNA", slot = "data")  ## genes x cells（log-normalized）
rankings <- AUCell::AUCell_buildRankings(expr_mat, plotStats = FALSE, nCores = 20)
auc_obj  <- AUCell::AUCell_calcAUC(list(OxStress = os_genes_use), rankings, nCores = 20)
auc_vec  <- as.numeric(AUCell::getAUC(auc_obj)["OxStress", ])
names(auc_vec) <- colnames(sce.all)
sce.all <- AddMetaData(sce.all, metadata = auc_vec, col.name = "OS_AUCell_OxStress")
stopifnot(all(features %in% colnames(sce.all@meta.data)))
md <- sce.all@meta.data[, features, drop = FALSE]
# 1) 直接平均（raw mean）
sce.all$OS_score_mean_raw <- rowMeans(md, na.rm = TRUE)
# 2) 标准化后平均（推荐：避免不同评分尺度不一致）
md_z <- scale(md)  # 每列做 (x-mean)/sd
sce.all$OS_score_mean_z <- rowMeans(md_z, na.rm = TRUE)
# 快速检查
summary(sce.all$OS_score_mean_raw)
summary(sce.all$OS_score_mean_z)
FeatureDimPlot(sce.all,features  = "OS_score_mean_z")
CellDensityPlot(sce.all,features  = "OS_score_mean_z")
save(sce.all, file = "tcos.Rdata")
