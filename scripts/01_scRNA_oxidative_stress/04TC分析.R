rm(list = ls())
gc()
library(Seurat)
library(SCP)
load("sceall.Rdata")
mycols <- c(
  "#4E79A7","#F28E2B","#E15759","#76B7B2","#59A14F",
  "#EDC948","#B07AA1","#FF9DA7","#9C755F","#BAB0AC",
  "#1F77B4","#FF7F0E","#2CA02C","#D62728","#9467BD",
  "#8C564B","#E377C2","#7F7F7F","#BCBD22","#17BECF",
  "#3B5B92","#E08D3C","#C23B3A","#2F9C95","#3E8E41",
  "#D6B11E","#8B5FBF","#F06A97","#7A5C44","#9EA3A8",
  "#0F4C81","#FFB000","#6FBE6B","#D64550","#7E57C2",
  "#00A6D6","#F5A6C6","#6B7280","#B8DE29","#FF6F61"
)
CellDimPlot(sce.all,group.by = "celltype",palcolor = mycols,label = T)
load("02ep.Rdata")
ep$RNA_snn_res.0.5
CellDimPlot(ep,group.by = "RNA_snn_res.0.5",palcolor = mycols)
dev.off()
load("03tc.Rdata")
DefaultAssay(epcnv) <- "RNA"
#epcnv$RNA_snn_res.0.5
#Idents(epcnv) <- "RNA_snn_res.0.5"
#tc <- subset(epcnv, idents = "3", invert = TRUE)
# 检查
table(Idents(tc))
seu <- tc


# NMF寻找可靠的功能 ---------------------------------------------------------------------
library(GeneNMF)
library(Seurat)
library(ggplot2)
library(UCell)
library(patchwork)
library(tidyr)
library(dplyr)
library(RColorBrewer)
library(viridis)
DefaultAssay(seu) <- "RNA"
seu$orig.ident
seu.list <- SplitObject(seu, split.by = "orig.ident")
geneNMF.programs <- multiNMF(seu.list, assay="RNA", k=4:9, min.exp = 0.05)
geneNMF.metaprograms <- getMetaPrograms(geneNMF.programs,
                                        metric = "cosine",
                                        specificity.weight = 3,
                                        weight.explained = 0.3,
                                        nMP=10)
plotMetaPrograms(geneNMF.metaprograms)
geneNMF.metaprograms$metaprograms.metrics
#lapply(geneNMF.metaprograms.filtered$metaprograms.genes, head)
geneNMF.metaprograms$metaprograms.genes
mp.genes <- geneNMF.metaprograms$metaprograms.genes
seu <- AddModuleScore_UCell(seu, features = mp.genes, ncores=4, name = "")
seu
# UMAP降维聚类 ----------------------------------------------------------------
sce <- seu
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
DimPlot(sce, group.by = "RNA_snn_res.0.2", label = TRUE)
DimPlot(sce, group.by = "RNA_snn_res.0.3", label = TRUE)
RunSlingshot(sce,group.by = "RNA_snn_res.0.1")
sce$cnv_fraction
FeaturePlot(sce,features = "cnv_fraction")
FeaturePlot(sce,features = "cnv_fraction")
#分辨率画树
library(clustree)
library(SCP)
clustree(sce@meta.data, prefix = "RNA_snn_res.")
table(sce$RNA_snn_res.0.1)
table(sce$RNA_snn_res.0.3)
Idents(sce) <- "RNA_snn_res.0.1"
FeaturePlot(sce,features = c("MP3","cnv_fraction","MP1"))
VlnPlot(sce,features = "MP3",group.by = "RNA_snn_res.0.05")
save(sce,geneNMF.metaprograms,file = "04nmfsce.Rdata")

load("04nmfsce.Rdata")





