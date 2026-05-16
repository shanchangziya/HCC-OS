# 取出结果表
df <- res_knk[["diffRegulation"]]

# 1) 筛选 p.value < 0.1
df_filt <- df[df$p.value < 0.05, , drop = FALSE]

# 2) 去除所有 HLA 相关基因（HLA开头；大小写不敏感）
df_filt <- df_filt[!grepl("^HLA", df_filt$gene, ignore.case = TRUE), , drop = FALSE]

# 3)（可选）按 p.value / p.adj / |Z| 排序，看起来更顺
df_filt <- df_filt[order(df_filt$p.value, -abs(df_filt$Z)), ]

# 4) 提取基因列表（去重）
genes_p01_noHLA <- unique(df_filt$gene)

# 查看数量与前几个
length(genes_p01_noHLA)
head(genes_p01_noHLA, 15)

# 如果你想导出
write.csv(df_filt, "diffRegulation_p_lt_0.1_noHLA.csv", row.names = FALSE)
# write.table(genes_p01_noHLA, "genes_p_lt_0.1_noHLA.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)


genes <- c("NQO1","CHL1","C1QA","C1QB","AIF1","C1QC","TYROBP","CD74",
           "HMGCS2","ALDOB","PXMP2","TM4SF4","FCER1G","PLPP3","IGSF6","MS4A7",
           "ACADM","APOC3","LAPTM5","CYP4A11","IL32","ACSM2B","CYBB","HSD17B6",
           "LGMN","MS4A6A","FCGR3A","AMACR","NFIB","GADD45A")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(clusterProfiler)
  library(enrichplot)
  library(org.Hs.eg.db)     # 人类；若小鼠换 org.Mm.eg.db
  library(ReactomePA)
})



# Symbol -> Entrez
eg <- bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
entrez <- unique(eg$ENTREZID)

# GO Biological Process
ego <- enrichGO(
  gene = entrez, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
  ont = "BP", pAdjustMethod = "BH",
  pvalueCutoff = 0.2, qvalueCutoff = 0.2, readable = TRUE
)


# 先看一下富集到多少条
c(GO_BP = nrow(as.data.frame(ego)))

head(as.data.frame(ego)[, c("ID","Description","p.adjust","geneID")], 20)



library(dplyr)
library(tidyr)
library(enrichplot)
library(Seurat)
library(ggplot2)

library(dplyr)
library(tidyr)
library(stringr)

ego_df <- as.data.frame(ego)

# 保险：有些情况下 ID 在 rownames
if (!"ID" %in% colnames(ego_df)) {
  ego_df <- tibble::rownames_to_column(ego_df, var = "ID")
}

# 你要的三个主题：用 Description 模糊匹配
term_mac  <- ego_df %>% filter(str_detect(tolower(Description), "macrophage activation")) %>% slice_min(p.adjust, n = 1)
term_fa   <- ego_df %>% filter(str_detect(tolower(Description), "fatty acid metabolic process")) %>% slice_min(p.adjust, n = 1)
term_igg  <- ego_df %>% filter(str_detect(tolower(Description), "immunoglobulin mediated immune response")) %>% slice_min(p.adjust, n = 1)

term_mac[, c("ID","Description","p.adjust","geneID")]
term_fa[,  c("ID","Description","p.adjust","geneID")]
term_igg[, c("ID","Description","p.adjust","geneID")]

# geneID 是 "A/B/C" 形式（readable=TRUE 时是 Symbol）
gs_mac <- unique(unlist(strsplit(term_mac$geneID, "/")))
gs_fa  <- unique(unlist(strsplit(term_fa$geneID,  "/")))
gs_igg <- unique(unlist(strsplit(term_igg$geneID, "/")))

gene_sets <- list(
  Macrophage_activation = gs_mac,
  Fatty_acid_metabolism = gs_fa,
  Ig_mediated_response  = gs_igg
)

# 看看每个基因集大小
sapply(gene_sets, length)


DefaultAssay(tcos) <- "RNA"

# 若没有tsne就运行
if (!"tsne" %in% names(tcos@reductions)) {
  # 确保有 PCA
  if (!"pca" %in% names(tcos@reductions)) {
    tcos <- ScaleData(tcos, verbose = FALSE)
    tcos <- RunPCA(tcos, npcs = 30, verbose = FALSE)
  }
  tcos <- RunTSNE(tcos, dims = 1:30, reduction = "pca", check_duplicates = FALSE)
}

library(UCell)
library(ggplot2)

# 只保留 tcos 里存在的基因
gene_sets2 <- lapply(gene_sets, function(gs) intersect(gs, rownames(tcos)))
sapply(gene_sets2, length)

tcos <- UCell::AddModuleScore_UCell(tcos, features = gene_sets2, name = "KOpath")
tcos$Macrophage_activationKOpath

score_cols <- c("Ig_mediated_responseKOpath",
                "Fatty_acid_metabolismKOpath",
                "Macrophage_activationKOpath")
library(Nebulosa)
plot_density(tcos, c("Ig_mediated_responseKOpath"),size = 0.3)
plot_density(tcos, c("Ig_mediated_responseKOpath", "Fatty_acid_metabolismKOpath",
                     "Macrophage_activationKOpath"), joint = TRUE,combine = FALSE) 

rm(list = ls())
dev.off()
gc()
