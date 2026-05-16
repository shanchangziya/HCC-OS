#####0 样本信息 ####
#使用GSE202642数据集做发现集
#过滤前110658个细胞
#初始过滤，标准细胞至少有5个基因表达才会被保留，基因至少在300个细胞表达才会被保留
#后续过滤，线粒体基因表达量小于20%
#过滤后98854个细胞,分辨率寻找0.3
rm(list=ls())
library(Seurat)
library(dplyr)
library(SCP)
load("data/03注释.Rdata")
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
DimPlot(sce.all, reduction = "umap", group.by = "RNA_snn_res.0.3",label = T) 
sel.clust = "RNA_snn_res.0.3"
sce.all <- SetIdent(sce.all, value = sel.clust)
genes_to_check = c('EPCAM', 'ALDH1A1', 'ALB',#EPI
                   'MS4A1', 'CD79A',#B
                   'CD3D', 'CD3E', #T
                   'NCAM1', 'FGFBP2', # NK,
                   'CD33','ITGAM',#MDSC
                   'CD68','CD163',"CD14",#Mono/mac
                   'ITGAX',"CD1E","CD1C",#DC
                   'COL1A2', 'ACTA2',  #fibroblasts  
                   'CD34','PECAM1',#EC
                   "PTPRC","S100A9","FCGR3B","CSF3R","S100A8"
)
DotPlot(sce.all,features = genes_to_check)
# 9：cycling
# 0：NK/T
# 5:Mono
# 4:DC
# 1,11,:Mac
# 2: Endo
# 3,6: epi
# 8,10: B
# 7:CAF
celltype=data.frame(ClusterID=0:11,
                    celltype= 0:11) 
#定义细胞亚群 
celltype[celltype$ClusterID %in% c(9),2]='Cycling'  
celltype[celltype$ClusterID %in% c(4),2]='DC'
celltype[celltype$ClusterID %in% c(0),2]='NK-T' 
celltype[celltype$ClusterID %in% c(2),2]='Endo' 
celltype[celltype$ClusterID %in% c(8,10),2]='B'   
celltype[celltype$ClusterID %in% c(1,11),2]='Mac' 
celltype[celltype$ClusterID %in% c(3,6),2]='epi' 
celltype[celltype$ClusterID %in% c(7),2]='fibo' 
celltype[celltype$ClusterID %in% c(5),2]='Neu'
table(celltype$celltype)
sce.all@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  sce.all@meta.data[which(sce.all@meta.data$RNA_snn_res.0.3 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(sce.all@meta.data$celltype)
CellDimPlot(sce.all,group.by = "celltype",palcolor = mycols,label = T)
#CellStatPlot(sce.all,group.by = "celltype",stat.by = "group",palcolor = mycols,label = T)
CellStatPlot(sce.all, stat.by = "group", group.by = "celltype", plot_type = "trend",palcolor = mycols
             )
CellDimPlot(sce.all,group.by = "celltype",palcolor = mycols)
save(sce.all,file = "sceall.Rdata")
dev.off()
gc()
