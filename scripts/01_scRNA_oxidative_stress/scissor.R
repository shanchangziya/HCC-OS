library(Seurat)
library(SCP)
#devtools::install_github('sunduanchen/Scissor')
library(Scissor)
load("tcos.Rdata")
load("bulk/exp1surv1.rdata")
#检查对应顺序
all(colnames(bulk_dataset) == bulk_survival$TCGA_patient_barcode)
#TRUE
phenotype <- surv1[,2:3]
colnames(phenotype)[1] <- c("status")
head(phenotype)
infos1 <- Scissor(exp1, tcos, phenotype, alpha = 0.05,
                  family = "cox", Save_file = 'Scissor_LIHC_survival.RData')
#Scissor+ Scissor- 定义
Scissor_select <- rep("Background cells", ncol(tcos))
names(Scissor_select) <- colnames(tcos)
Scissor_select[infos1$Scissor_pos] <- "Scissor+ cell"
Scissor_select[infos1$Scissor_neg] <- "Scissor- cell"
#metadata 中添加 Scissor信息
tcos <- AddMetaData(tcos, metadata = Scissor_select, col.name = "scissor")
DimPlot(tcos, reduction = 'umap', group.by = 'scissor', cols = c('grey','royalblue','indianred1'))
library(SCP)
CellDimPlot(tcos,group.by = 'scissor',show_stat = F,label = T)
library(scales)
library(ggplot2)
library(dplyr)

# 分组顺序一定要和 UMAP 一致
lvl <- c("Background cells", "Scissor- cell", "Scissor+ cell")  # 你的三类
cols_umap_like <- c("#a6cee3","#1f78b4","#b2df8a")
cols_umap_like

p_bar2 <- p_bar +
  scale_fill_manual(values = cols_umap_like) +
  guides(fill = guide_legend(title = NULL))

p_bar2

rm(list = ls())
gc()
