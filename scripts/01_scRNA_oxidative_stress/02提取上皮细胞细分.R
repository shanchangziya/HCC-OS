rm(list = ls())
gc()
library(Seurat)
load("sceall.Rdata")
table(sce.all$celltype)
# 提取 epi 细胞为新的 Seurat 对象
sce <- subset(sce.all, subset = celltype == "epi")
table(sce@meta.data$celltype)
#####step2重跑单细胞####
sce
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
sce <- FindClusters(sce, resolution = 0.1)
sce <- FindClusters(sce, resolution = 0.2)
sce <- FindClusters(sce, resolution = 0.3)
sce <- FindClusters(sce, resolution = 0.05)
sce <- FindClusters(sce, resolution = 0.5)
sce <- RunUMAP(sce, reduction = "harmony", dims = 1:20, seed.use = 42)
DimPlot(sce, group.by = "RNA_snn_res.0.1", label = TRUE)
DimPlot(sce, group.by = "RNA_snn_res.0.05", label = TRUE)
DimPlot(sce, group.by = "RNA_snn_res.0.3", label = TRUE)
#分辨率画树
library(Seurat)
library(clustree)
library(SCP)
clustree(sce@meta.data, prefix = "RNA_snn_res.")
table(sce$RNA_snn_res.0.1)
###删除细胞数少于100的细胞群,且异常的细胞群####
# 删除细胞数 < 100 的簇（以 RNA_snn_res.0.1 为例）
Idents(sce) <- "RNA_snn_res.0.1"
# 找到要删除的簇
tab <- table(Idents(sce))
rm.clusters <- names(tab[tab < 100])
# 子集保留其余簇
sce.filt <- subset(sce, idents = rm.clusters, invert = TRUE)
#接下来分析，按照分辨率为0.1进行 
DimPlot(sce.filt,group.by = "RNA_snn_res.0.3")
DimPlot(sce.filt,group.by = "RNA_snn_res.0.5")
table(sce$RNA_snn_res.0.5)
ep <- sce.filt
save(ep,file = "02ep.Rdata")

