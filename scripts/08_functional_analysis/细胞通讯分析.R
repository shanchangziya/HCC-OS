library(Seurat)
load("注释sceall.Rdata")
set.seed(123)

# 查看原始细胞类型分布
cat("\n原始细胞类型统计：\n")
print(table(sce.all$final_celltype))

# 移除正常上皮细胞
cat("\n移除正常上皮细胞...\n")
sce.noNE <- subset(sce.all, subset = final_celltype != "Normal_Epithelial")

# 创建合并后的大类细胞类型
cat("\n创建合并后的大类细胞类型标签...\n")

# 定义细胞类型合并规则
merge_celltype <- function(celltype) {
  # 保持肿瘤细胞OS分组
  if (grepl("^Tumor_OS-high", celltype)) {
    return("Tumor_OS-high")
  }
  if (grepl("^Tumor_OS-low", celltype)) {
    return("Tumor_OS-low")
  }
  
  # T细胞合并（包括前缀T:的情况）
  if (grepl("^T:|CD[48]_|Treg|Tfh|Th17|MAIT|gdT|Cyto.*T|Naive_T|Memory_T|Effector_T|Trm", celltype)) {
    return("T")
  }
  
  # NK细胞合并（包括前缀NK:的情况）
  if (grepl("^NK:|_NK$|CD56", celltype)) {
    return("NK")
  }
  
  # 巨噬细胞合并
  if (grepl("Macrophage|TAM|Mac_|^Mac$", celltype, ignore.case = TRUE)) {
    return("Macrophage")
  }
  
  # B细胞合并
  if (grepl("^B_|_B$|Plasma|^B$|Naive_B|Memory_B|Activated_B", celltype)) {
    return("B")
  }
  
  # DC细胞合并
  if (grepl("^DC|cDC|pDC|moDC|Dendritic", celltype)) {
    return("DC")
  }
  
  # CAF合并
  if (grepl("CAF|Fibroblast", celltype)) {
    return("CAF")
  }
  
  # 内皮细胞合并
  if (grepl("Endothelial|VEC|LSEC|LEC|Endo", celltype)) {
    return("Endothelial")
  }
  
  # 中性粒细胞合并
  if (grepl("Neutrophil|MDSC|TAN|PMN|Inflam_Neu", celltype)) {
    return("Neutrophil")
  }
  
  # 其他细胞保持原名
  return(celltype)
}

# 应用合并规则
sce.noNE$major_celltype <- sapply(sce.noNE$final_celltype, merge_celltype)

# 显示合并后的细胞类型统计
cat("\n合并后的细胞类型统计：\n")
print(table(sce.noNE$major_celltype))


set.seed(123)

# 抽样比例
prop <- 0.30

cells_all <- colnames(sce.noNE)
n_keep <- floor(length(cells_all) * prop)

cells_use <- sample(cells_all, size = n_keep, replace = FALSE)

sce.noNE_30 <- subset(sce.noNE, cells = cells_use)

# 检查抽样后细胞数与各类分布
dim(sce.noNE_30)
print(table(sce.noNE_30$major_celltype))
prop.table(table(sce.noNE_30$major_celltype))


# 细胞通讯分析 ------------------------------------------------------------------

library(CellChat)
seurat_obj <- sce.noNE_30
# 1. 用 RNA 的 counts（如果有 data 也可以改成 "data"）
DefaultAssay(seurat_obj) <- "RNA"

data.input <- GetAssayData(
  seurat_obj,
  assay = "RNA",
  slot  = "counts"   # 如果你确认 data 不为空，也可以改成 "data"
)

dim(data.input)      # 确认这里是 “很多行 × 很多列”，不能是 0 行

# 2. 准备分组信息，用刚才合好的 celltype_refined
meta <- data.frame(
  celltype = seurat_obj$major_celltype,
  row.names = colnames(seurat_obj)
)

# 3. 构建 CellChat 对象
cellchat <- createCellChat(
  object  = data.input,
  meta    = meta,
  group.by = "celltype"
)

# 4. 绑定小鼠数据库
data("CellChatDB.human")
CellChatDB <- CellChatDB.human
cellchat@DB <- CellChatDB

# 只用分泌信号通路（可按需调整）
cellchat@DB$interaction <- subset(CellChatDB$interaction,
                                  annotation == "Secreted Signaling")
cellchat@DB$complex  <- CellChatDB$complex
cellchat@DB$cofactor <- CellChatDB$cofactor

# 5. 子集 + 检查是否还有基因
cellchat <- subsetData(cellchat)
dim(cellchat@data)   # 这里应该是 “>0 行 × 细胞数”

# 6. 后续标准流程
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat@data.raw <- cellchat@data.signaling
cellchat <- computeCommunProb(cellchat, raw.use = T)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
save(cellchat,file = "cellchat.Rdata")
# 三、专门看 CD27+ IgM memory-like B 的通讯 ---------------------------------------
# 所有 outgoing（CD27+ IgM memory-like B 做 sender）
comm_out <- subsetCommunication(cellchat,
                                sources.use = "CD27+ IgM memory-like B")
# 所有 incoming（CD27+ IgM memory-like B 做 receiver）
comm_in  <- subsetCommunication(cellchat,
                                targets.use = "CD27+ IgM memory-like B")
#2. 圈图：展示 CD27+ IgM memory-like B 的出/入信号（图形好看）

#（1）整体细胞通讯圈图（先有一个 big picture）
netVisual_circle(cellchat@net$count,
                 vertex.weight = groupSize,
                 weight.scale  = TRUE,
                 label.edge    = FALSE,
                 title.name    = "Number of interactions")

# 第二张：Interaction weights/strength
netVisual_circle(cellchat@net$weight,
                 vertex.weight = groupSize,
                 weight.scale  = TRUE,
                 label.edge    = FALSE,
                 title.name    = "Interaction weights/strength")


#（2）CD27+ IgM memory-like B 作为 sender 的圈图
netVisual_circle(cellchat@net$count, vertex.weight = groupSize,sources.use = "Tumor_OS-high", 
                 weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$count, vertex.weight = groupSize,targets.use = "Tumor_OS-high", 
                 weight.scale = T, label.edge= F, title.name = "Number of interactions")
#3 细胞点图
cellchat@netP$pathways
levels(cellchat@idents) 
vertex.receiver = c(5,2,3) 
pathways.show <- "CHEMERIN"
netVisual_aggregate(cellchat, signaling = pathways.show,  
                    vertex.receiver = vertex.receiver,layout = "hierarchy")
netVisual_bubble(cellchat, sources.use = c(9),targets.use = c(2,5),remove.isolate = FALSE)
VlnPlot(sce.noNE_30, features = c("CXCR4","CD74",""),group.by = "major_celltype")
table(sce.noNE_30$major_celltype)
#4MIF通路
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
table(meta$celltype)
# (1) show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_bubble(cellchat, sources.use = 9, targets.use = c(5,2,3), remove.isolate = FALSE)
#> Comparing communications on a single object

netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
library(ggplot2)
dotPlot(sce.noNE_30, features = c("CD74","NQO1","MIF"),group.by = "major_celltype")+coord_equal()
table(sce.noNE_30$major_celltype)
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, 
                                  width = 8, height = 2.5, font.size = 10)

save(cellchat,file = "cellchat.Rdata")
rm(list = ls())
dev.off()
gc()

