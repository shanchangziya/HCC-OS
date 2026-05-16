library(SCP)
library(Seurat)
load("08_Mac_annotated.Rdata")
load("注释sceall.Rdata")
print(table(sce.all$final_celltype))
table(sce.mac$Mac_subtype)
table(mac.obj$group_comm)
CellDimPlot(mac.obj,group.by = "group_comm",show_stat = F)
library(Seurat)
library(ggplot2)

Idents(mac.obj) <- mac.obj$group_comm

markers_tam <- list(
  "HSP_TAM"        = c("HSPA1A","HSPA1B","HSP90AA1","DNAJB1","HSPD1"),
  "SLC40A1_TAM"    = c("SLC40A1","FTH1","FTL","HMOX1","SOD2"),
  "TREM2_TAM"      = c("TREM2","APOE","LST1","C1QC","CTSD")
)

# 展开成向量（并保留分组顺序）
features_use <- unique(unlist(markers_tam))

p_dot <- DotPlot(
  mac.obj,
  features = features_use,
  group.by = "group_comm",
  cols = c("lightgrey", "red"),
  dot.scale = 6
) +
  RotatedAxis() +
  theme_classic(base_size = 12) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4)
  )

p_dot
rm(list = ls())
gc()
dev.off()
# 提出三个巨噬细胞亚群开始-----------------------------
# A) 从 sce.mac 取 3 个巨噬细胞亚群
# -----------------------------
Idents(sce.mac) <- sce.mac$Mac_subtype

mac_keep <- c("HSP_TAM", "SLC40A1_TAM", "TREM2_TAM")
mac.obj <- subset(sce.mac, idents = mac_keep)

# 给细胞加前缀，避免与 sce.all 的细胞名冲突
mac.obj <- RenameCells(mac.obj, add.cell.id = "MAC")

# 统一一个用于通讯的标签列（建议叫 celltype / group）
mac.obj$group_comm <- mac.obj$Mac_subtype
mac.obj$source_obj <- "sce.mac"


# -----------------------------
# B) 从 sce.all 取 Tumor_OS-high
# 提出三个巨噬细胞亚群结束-----------------------------
Idents(sce.all) <- sce.all$final_celltype

tum.obj <- subset(sce.all, idents = "Tumor_OS-high")
tum.obj <- RenameCells(tum.obj, add.cell.id = "TUMOR")

tum.obj$group_comm <- "Tumor_OS-high"
tum.obj$source_obj <- "sce.all"

sce.comm <- merge(
  x = mac.obj,
  y = tum.obj,
  project = "Mac3plusTumorHighOS"
)

# 把通讯分组设为 Idents（CellChat 常用）
Idents(sce.comm) <- sce.comm$group_comm

table(Idents(sce.comm))
table(sce.comm$source_obj)

set.seed(123)
prop <- 0.30
sce.noNE <- sce.comm
cells_all <- colnames(sce.noNE)
n_keep <- floor(length(cells_all) * prop)

cells_use <- sample(cells_all, size = n_keep, replace = FALSE)

sce.noNE_30 <- subset(sce.noNE, cells = cells_use)

# 检查抽样后细胞数与各类分布
dim(sce.noNE_30)
print(table(sce.noNE_30$group_comm))

# 开始通讯分析 ------------------------------------------------------------------
library(CellChat)
seurat_obj <- sce.noNE_30
DefaultAssay(seurat_obj) <- "RNA"
data.input <- GetAssayData(
  seurat_obj,
  assay = "RNA",
  slot  = "counts"   # 如果你确认 data 不为空，也可以改成 "data"
)

dim(data.input)      # 确认这里是 “很多行 × 很多列”，不能是 0 行
meta <- data.frame(
  celltype = seurat_obj$group_comm,
  row.names = colnames(seurat_obj)
)
cellchat <- createCellChat(
  object  = data.input,
  meta    = meta,
  group.by = "celltype"
)
data("CellChatDB.human")
CellChatDB <- CellChatDB.human
cellchat@DB <- CellChatDB
cellchat@DB$interaction <- subset(CellChatDB$interaction,
                                  annotation == "Secreted Signaling")
cellchat@DB$complex  <- CellChatDB$complex
cellchat@DB$cofactor <- CellChatDB$cofactor

# 5. 子集 + 检查是否还有基因
cellchat <- subsetData(cellchat)
dim(cellchat@data)   # 这里应该是 “>0 行 × 细胞数”


cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat@data.raw <- cellchat@data.signaling
cellchat <- computeCommunProb(cellchat, raw.use = T)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

netVisual_circle(cellchat@net$count,
                 weight.scale  = TRUE,
                 label.edge    = FALSE,
                 title.name    = "Number of interactions")
netVisual_circle(cellchat@net$weight,
                 weight.scale  = TRUE,
                 label.edge    = FALSE,
                 title.name    = "Interaction weights/strength")
netVisual_circle(cellchat@net$count,sources.use = "Tumor_OS-high", 
                 weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$count, vertex.weight = groupSize,targets.use = "Tumor_OS-high", 
                 weight.scale = T, label.edge= F, title.name = "Number of interactions")

netVisual_bubble(cellchat, sources.use = c(4),targets.use = c(1,2,3),remove.isolate = FALSE)


netAnalysis_signalingRole_scatter(cellchat)
pheatmap::pheatmap(cellchat@net$count, border_color = "black", 
                   cluster_cols = F, fontsize = 10, cluster_rows = F,
                   display_numbers = T,number_color="black",number_format = "%.0f")
maccellchat <- cellchat

save(maccellchat,file = "maccellchat.Rdata")
load("maccellchat.Rdata")
cellchat <- maccellchat
load("gokegg.Rdata")
deg_os$gene[1:10]

rm(list = ls())
gc()
dev.off

