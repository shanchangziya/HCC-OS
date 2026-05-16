library(Seurat)
library(SCP)
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})
load("tcos.Rdata")
table(tcos$OS_group_valley)
set.seed(123)

celltype_col <- "OS_group_valley"   # ← 改成你的细胞类型列名
n_per_type <- 500

# 用细胞类型作为身份
Idents(tcos) <- celltype_col

# 每个细胞类型抽样
cells_use <- unlist(lapply(levels(Idents(tcos)), function(ct){
  cells <- WhichCells(tcos, idents = ct)
  if (length(cells) <= n_per_type) cells else sample(cells, n_per_type)
}))
# 子集对象
tcos_500 <- subset(tcos, cells = cells_use)
# 检查：每类数量（有的类可能 <500）
table(tcos_500[[celltype_col]][,1])
tcos_500 <- RunDEtest(tcos_500, group_by = "OS_group_valley", fc.threshold = 1, only.pos = FALSE)
VolcanoPlot(tcos_500, group_by = "OS_group_valley")
DEGs <- tcos_500@tools$DEtest_OS_group_valley$AllMarkers_wilcox
DEGs <- DEGs[with(DEGs, avg_log2FC > 1 & p_val_adj < 0.05), ]
table(tcos_500$OS_group_valley)

# 单细胞方法的差异分析 --------------------------------------------------------------------
# 以分组作为身份
Idents(tcos_500) <- "OS_group_valley"
table(Idents(tcos_500))

# 选择 assay（优先 RNA；如果你主要用 SCT 也可以改）
DefaultAssay(tcos_500) <- if ("RNA" %in% Assays(tcos_500)) "RNA" else DefaultAssay(tcos_500)

# 如果你用 SCT 且想在 SCT 上做 DE，需要先 PrepSCTFindMarkers
if (DefaultAssay(tcos_500) == "SCT") {
  tcos_500 <- PrepSCTFindMarkers(tcos_500)
}

deg_os <- FindMarkers(
  tcos_500,
  ident.1 = "OS-high",
  ident.2 = "OS-low",
  test.use = "wilcox",
  min.pct = 0.10,
  logfc.threshold = 0.25,
  only.pos = FALSE
)

deg_os <- deg_os %>%
  tibble::rownames_to_column("gene") %>%
  arrange(p_val_adj, desc(avg_log2FC))

head(deg_os, 20)
# 画气泡图 ---------------------------------------------------------------------
# 确保分组列存在
Idents(tcos_500) <- "OS_group_valley"

deg_os$p_val_adj
top10 <- deg_os %>%
  filter(!is.na(p_val_adj)) %>%
  filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 10) %>%
  pull(gene)
top10

# 2) DotPlot（颜色=平均表达，点大小=表达细胞比例）
p <- DotPlot(
  object = tcos_500,
  features = top10,
  group.by = "OS_group_valley",
  dot.scale = 7
) +
  scale_color_gradientn(
    colours = c("#2C7BB6", "#ABD9E9", "#FFFFBF", "#FDAE61", "#D7191C") # 蓝-白-红
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major = element_line(size = 0.2, color = "grey90"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title = element_blank()
  ) +
  labs(color = "Avg. expression", size = "Pct. expressed")

p
top10





rm(list = ls())
gc()
dev.off()

# 拟转录组差异分析 ----------------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(edgeR)
  library(dplyr)
  library(ggplot2)
})

# 1) 用 RNA raw counts（不要用 SCT/data）
DefaultAssay(tcos_500) <- if ("RNA" %in% Assays(tcos_500)) "RNA" else DefaultAssay(tcos_500)
Idents(tcos_500) <- "OS_group_valley"

# 2) pseudobulk：按 OS_group_valley + orig.ident 聚合
pb <- AggregateExpression(
  tcos_500,
  assays = "RNA",
  group.by = c("OS_group_valley", "orig.ident"),
  slot = "counts",
  return.seurat = FALSE
)$RNA  # genes x (group_sample)

# 3) 解析分组信息（列名一般像 "OS-high_Tumor1"）
pb_meta <- data.frame(pbsample = colnames(pb)) %>%
  mutate(group = sub("_.*$", "", pbsample)) %>%
  mutate(group = factor(group, levels = c("OS-low","OS-high")))

# 4) edgeR 差异分析（QLF）
y <- DGEList(counts = pb, group = pb_meta$group)
keep <- filterByExpr(y, group = pb_meta$group)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)

design <- model.matrix(~ group, data = pb_meta)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = "groupOS-high")  # OS-high vs OS-low
qlf[['table']][1:4,1:4]
library(dplyr)
library(tibble)
deg_pb <- topTags(qlf, n = Inf)$table %>%
  as.data.frame() %>%
  rownames_to_column("gene")

head(deg_pb)

# 5) 火山图
deg_pb <- deg_pb %>%
  mutate(sig = case_when(
    PValue < 0.05 & logFC >  0.25 ~ "Up in OS-high",
    PValue < 0.05 & logFC < -0.25 ~ "Up in OS-low",
    TRUE ~ "NS"
  ))

p_vol <- ggplot(deg_pb, aes(x = logFC, y = -log10(PValue), color = sig)) +
  geom_point(size = 1.2, alpha = 0.85) +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  scale_color_manual(
    values = c(
      "Up in OS-high" = "#D7191C",  # 红：对应 DotPlot 的高表达端
      "Up in OS-low"  = "#2C7BB6",  # 蓝：对应 DotPlot 的低表达端
      "NS"            = "grey80"
    )
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "right",
    panel.grid.major = element_line(size = 0.2, color = "grey90"),
    panel.grid.minor = element_blank()
  ) +
  labs(x = "log2FC (OS-high vs OS-low)", y = "-log10(P)", color = NULL)
p_vol
ggsave("volcano_pseudobulk_edgeR_OS-high_vs_OS-low.pdf", p_vol, width = 6.5, height = 5)




# GOkegg分析 ----------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(clusterProfiler)
  library(enrichplot)
  library(org.Hs.eg.db)
})

# =========================
# 0) 参数 & 配色（和你前面一致的蓝-白-红体系）
# =========================
fdr_cut <- 0.05
lfc_cut <- 0.25
pal5 <- c("#2C7BB6", "#ABD9E9", "#FFFFBF", "#FDAE61", "#D7191C")  # 蓝->红


FDR
deg_pb$PValue
# =========================
# 1) 准备基因集：Up / Down + universe（建议给 universe，减少偏倚）
# =========================
up_genes   <- deg_pb %>% filter(PValue < fdr_cut, logFC >  lfc_cut) %>% pull(gene) %>% unique()
down_genes <- deg_pb %>% filter(PValue < fdr_cut, logFC < -lfc_cut) %>% pull(gene) %>% unique()
universe_genes <- deg_pb %>% pull(gene) %>% unique()

cat("Up genes:", length(up_genes), "\n")
cat("Down genes:", length(down_genes), "\n")

# SYMBOL -> ENTREZID（Bayes/edgeR 常用 SYMBOL；KEGG/GO 更稳用 ENTREZ）
sym2ent <- function(sym){
  bitr(sym, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) %>%
    distinct(SYMBOL, .keep_all = TRUE)
}
map_up <- sym2ent(up_genes)
map_dn <- sym2ent(down_genes)
map_uni <- sym2ent(universe_genes)

up_ent <- unique(map_up$ENTREZID)
dn_ent <- unique(map_dn$ENTREZID)
uni_ent <- unique(map_uni$ENTREZID)

# =========================
# 2) GO 富集（BP/CC/MF 一起跑，便于一张图展示）
# =========================
ego_up <- enrichGO(
  gene          = up_ent,
  universe      = uni_ent,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

ego_dn <- enrichGO(
  gene          = dn_ent,
  universe      = uni_ent,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

# 保存结果表
write.csv(as.data.frame(ego_up), "GO_up_enrich.csv", row.names = FALSE)
write.csv(as.data.frame(ego_dn), "GO_down_enrich.csv", row.names = FALSE)

# GO dotplot（Up/Down 分开画，且按ONTOLOGY分面）
p_go_up <- dotplot(ego_up, showCategory = 5, split = "ONTOLOGY") +
  facet_grid(ONTOLOGY ~ ., scales = "free_y") +
  scale_color_gradientn(colours = rev(pal5)) +  # 让更显著(更小p.adjust)更偏红
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank()) +
  labs(title = "GO enrichment (Up in OS-high)", color = "p.adjust")

p_go_dn <- dotplot(ego_dn, showCategory = 5, split = "ONTOLOGY") +
  facet_grid(ONTOLOGY ~ ., scales = "free_y") +
  scale_color_gradientn(colours = rev(pal5)) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank()) +
  labs(title = "GO enrichment (Up in OS-low)", color = "p.adjust")

ggsave("GO_dotplot_up.pdf", p_go_up, width = 8, height = 7)
ggsave("GO_dotplot_down.pdf", p_go_dn, width = 8, height = 7)

# （可选）GO 网络图：emapplot（更“高级”，但需要 term 相似度）
ego_up2 <- pairwise_termsim(ego_up)
ego_dn2 <- pairwise_termsim(ego_dn)

p_go_emap_up <- emapplot(ego_up2, showCategory = 30) + ggtitle("GO emapplot (Up in OS-high)")
p_go_emap_dn <- emapplot(ego_dn2, showCategory = 30) + ggtitle("GO emapplot (Up in OS-low)")

ggsave("GO_emap_up.pdf", p_go_emap_up, width = 9, height = 7)
ggsave("GO_emap_down.pdf", p_go_emap_dn, width = 9, height = 7)

# =========================
# 3) KEGG 富集（hsa）
# =========================
ekk_up <- enrichKEGG(
  gene          = up_ent,
  universe      = uni_ent,
  organism      = "hsa",
  keyType       = "kegg",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2
)
ekk_dn <- enrichKEGG(
  gene          = dn_ent,
  universe      = uni_ent,
  organism      = "hsa",
  keyType       = "kegg",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2
)

# KEGG 结果转可读基因名
ekk_up <- setReadable(ekk_up, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
ekk_dn <- setReadable(ekk_dn, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

write.csv(as.data.frame(ekk_up), "KEGG_up_enrich.csv", row.names = FALSE)
write.csv(as.data.frame(ekk_dn), "KEGG_down_enrich.csv", row.names = FALSE)

p_kegg_up <- dotplot(ekk_up, showCategory = 10) +
  scale_color_gradientn(colours = rev(pal5)) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank()) +
  labs(title = "KEGG enrichment (Up in OS-high)", color = "p.adjust")

p_kegg_dn <- dotplot(ekk_dn, showCategory = 15) +
  scale_color_gradientn(colours = rev(pal5)) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank()) +
  labs(title = "KEGG enrichment (Up in OS-low)", color = "p.adjust")

ggsave("KEGG_dotplot_up.pdf", p_kegg_up, width = 8, height = 6)
ggsave("KEGG_dotplot_down.pdf", p_kegg_dn, width = 8, height = 6)

# （可选）KEGG emapplot
ekk_up2 <- pairwise_termsim(ekk_up)
ekk_dn2 <- pairwise_termsim(ekk_dn)

p_kegg_emap_up <- emapplot(ekk_up2, showCategory = 5) + ggtitle("KEGG emapplot (Up in OS-high)")
p_kegg_emap_dn <- emapplot(ekk_dn2, showCategory = 5) + ggtitle("KEGG emapplot (Up in OS-low)")

ggsave("KEGG_emap_up.pdf", p_kegg_emap_up, width = 9, height = 7)
ggsave("KEGG_emap_down.pdf", p_kegg_emap_dn, width = 9, height = 7)

# 打印图（交互式会话里可直接看）
p_go_up
p_kegg_up

rm(list = ls())
gc()
dev.off()
save(ekk_up,ego_up,deg_pb,down_genes,up_genes,deg_os,file = "gokegg.Rdata")
