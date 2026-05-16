load("sceall.Rdata")
table(sce.all$celltype)
library(Seurat)
library(ggplot2)
library(SCP)
Idents(sce.all) <- sce.all$celltype
DimPlot(sce.all,label=T,group.by = "celltype")
genes_to_check = c('EPCAM', 'ALDH1A1', 'ALB',#EPI
                   'MS4A1', 'CD79A',#B
                   'CD3D', 'CD3E', #T
                   'NCAM1', 'FGFBP2', # NK,
                   'CD33','ITGAM',#MDSC
                   'CD68','CD163',"CD14",#Mono/mac
                   'ITGAX',"CD1E","CD1C",#DC
                   'COL1A2', 'ACTA2',  #fibroblasts  
                   'CD34','PECAM1',#EC
                   "PTPRC","S100A9","FCGR3B","CSF3R","S100A8",
                   "MKI67"
)
DotPlot(sce.all,group.by = "celltype",features = genes_to_check)+coord_flip()

#可视化在后面进行
sce.all
CellDimPlot(sce.all,group.by = "celltype",label = T)

# 评分 ----------------------------------------------------------------------

## ========= 0) 准备：氧化应激基因集（优先用 MSigDB hallmark ROS） =========
## 推荐用：HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY
## 若你没有 msigdbr 包，会自动用一个常用的抗氧化/应激手工基因集作为 fallback

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


## ========= 4) 看哪个细胞群更高（按 seurat_clusters；也可换成你的 celltype 列） =========
## 小提琴图：三种评分并排看
VlnPlot(
  sce.all,
  features = c("OS_AddModuleScore_1", "OxStressOS_UCell_", "OS_AUCell_OxStress"),
  group.by = "celltype",
  pt.size = 0
)

## 计算每个cluster的平均分（快速表格）
score_df <- FetchData(
  sce.all,
  vars = c("celltype", "OS_AddModuleScore_1", "OxStressOS_UCell_", "OS_AUCell_OxStress")
)
cluster_mean <- aggregate(
  score_df[, c("OS_AddModuleScore_1","OxStressOS_UCell_","OS_AUCell_OxStress")],
  by = list(cluster = score_df$celltype),
  FUN = mean
)
cluster_mean[order(cluster_mean$OxStressOS_UCell_, decreasing = TRUE), ]

FeaturePlot(sce.all,
            features =c("OS_AddModuleScore_1","OxStressOS_UCell_","OS_AUCell_OxStress"))
library(SCP)
FeatureDimPlot(sce.all,
               features =c("OS_AddModuleScore_1","OxStressOS_UCell_","OS_AUCell_OxStress"))
FeatureDimPlot(sce.all,features = c("OxStressOS_UCell_"),show_stat = F)
FeatureDimPlot(sce.all,features = c("OS_AddModuleScore_1"),show_stat = F)
FeatureDimPlot(sce.all,features = c("OS_AUCell_OxStress"),show_stat = F)
save(sce.all,file = "01sceall.Rdata")

# 绘图 ----------------------------------------------------------------------
CellDimPlot(sce.all,group.by = "celltype",label = T)
FeatureStatPlot(sce.all,group.by = "celltype",stat.by = "OxStressOS_UCell_")
FeatureStatPlot(sce.all,group.by = "celltype",stat.by = "OS_AddModuleScore_1")
FeatureStatPlot(sce.all,group.by = "celltype",stat.by = "OS_AUCell_OxStress")



rm(list = ls())
gc()
dev.off()
