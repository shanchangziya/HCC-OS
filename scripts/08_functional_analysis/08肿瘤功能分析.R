library(Seurat)
library(scMetabolism)
library(dplyr)
library(ggplot2)
library(stringr)
load("tcos.Rdata")
load("gokegg.Rdata")
library(dplyr)
library(ggplot2)
library(stringr)
library(forcats)
library(scales)

## 你的富集结果（如果还没准备）
goresult <- as.data.frame(ego_up)

plot_go_bar <- function(goresult, ontology = "CC", top_n = 10,
                        p_col = c("p.adjust","pvalue")[1],   # 默认用FDR；想用pvalue就改成 "pvalue"
                        wrap_width = 42,
                        show_right = c("none","Count","GeneRatio")[2], # 右侧显示什么
                        fill_color = "#7f4ea8"){
  
  stopifnot(p_col %in% colnames(goresult))
  
  df <- goresult %>%
    filter(ONTOLOGY == ontology, .data[[p_col]] < 0.05) %>%
    arrange(.data[[p_col]]) %>%
    slice_head(n = top_n) %>%
    mutate(
      score = -log10(.data[[p_col]]),
      Description_wrapped = str_wrap(Description, width = wrap_width),
      Description_wrapped = factor(Description_wrapped, levels = rev(unique(Description_wrapped))),
      GeneRatioNum = {
        m <- str_match(GeneRatio, "^(\\d+)\\/(\\d+)$")
        as.numeric(m[,2]) / as.numeric(m[,3])
      }
    )
  
  xmax <- max(df$score, na.rm = TRUE)
  left_pad  <- xmax * 0.02
  right_pad <- xmax * 0.18
  
  p <- ggplot(df, aes(x = score, y = Description_wrapped)) +
    geom_col(fill = alpha(fill_color, 0.35), width = 0.78) +
    ## 左侧 term 文本（固定在接近0的位置，更稳）
    geom_text(aes(x = left_pad, label = Description_wrapped),
              hjust = 0, vjust = 0.5, size = 4.2, lineheight = 0.95) +
    ## 右侧可选注释
    {
      if(show_right == "Count"){
        geom_text(aes(x = score + xmax*0.02, label = paste0("n=", Count)),
                  hjust = 0, vjust = 0.5, size = 3.6)
      } else if(show_right == "GeneRatio"){
        geom_text(aes(x = score + xmax*0.02, label = GeneRatio),
                  hjust = 0, vjust = 0.5, size = 3.6)
      } else NULL
    } +
    ## 让左侧文本有空间，右侧注释也不会被裁切
    coord_cartesian(clip = "off", xlim = c(0, xmax + right_pad)) +
    scale_x_continuous(expand = c(0, 0)) +
    labs(
      x = ifelse(p_col == "p.adjust", expression(-log[10]("FDR (BH)")), expression(-log[10]("p-value"))),
      y = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid   = element_blank(),
      axis.line.x  = element_line(linewidth = 0.8, color = "black"),
      axis.ticks.x = element_line(linewidth = 0.8, color = "black"),
      plot.margin  = margin(t = 10, r = 18, b = 10, l = 10)
    )
  
  p
}

p_mf <- plot_go_bar(goresult, ontology = "CC", top_n = 10,
                    p_col = "p.adjust",        # 推荐：FDR
                    show_right = "Count")      # 右侧显示 Count（可改 GeneRatio / none）
p_mf


ego_up
ekk_up
save(ego_gsea,file = "egobpgsea.Rdata")
colnames(deg_os)
countexp.Seurat <- sc.metabolism.Seurat(obj = tcos,  #Seuratde单细胞object
                                        method = "AUCell", 
                                        imputation = F, 
                                        ncores = 20, 
                                        metabolism.type = "KEGG")
score <- countexp.Seurat@assays$METABOLISM$score
score[1:4,1:4]
score_change <- score %>% 
  select_all(~str_replace_all(., "\\.", "-"))  #基因ID不规范会报错,下划线替换-
#确定细胞barcode椅子
identical(colnames(score_change) , rownames(countexp.Seurat@meta.data))
countexp.Seurat@meta.data <- cbind(countexp.Seurat@meta.data,t(score_change) )
Idents(countexp.Seurat) <- countexp.Seurat$OS_group_valley
# 画图
DotPlot.metabolism(
  obj = countexp.Seurat,
  pathway = input.pathway,
  phenotype = "OS_group_valley",
  norm = "y"
)
# gsea分析 ------------------------------------------------------------------
library(dplyr)

# 去掉 NA、重复基因（保留 |log2FC| 最大的那条）
deg_os2 <- deg_os %>%
  filter(!is.na(gene), !is.na(avg_log2FC)) %>%
  group_by(gene) %>%
  slice_max(order_by = abs(avg_log2FC), n = 1, with_ties = FALSE) %>%
  ungroup()

geneList <- deg_os2$avg_log2FC
names(geneList) <- deg_os2$gene
geneList <- sort(geneList, decreasing = TRUE)

# quick check
head(geneList)
tail(geneList)


library(clusterProfiler)
library(org.Hs.eg.db)   # 小鼠用 org.Mm.eg.db
library(enrichplot)
library(ggplot2)

set.seed(123)

ego_gsea <- gseGO(
  geneList      = geneList,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",          # BP / MF / CC
  minGSSize     = 10,
  maxGSSize     = 500,
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  verbose       = FALSE
)
library(GseaVis)
#devtools::install_github("junjunlab/GseaVis")
gseaNb(object = ego_gsea,
       geneSetID = 'GO:0016064',
       newGsea = T,addPval = T)
head(ego_gsea@result[, c("ID","Description")])

# 可视化（常用）
dotplot(ego_gsea, showCategory = 6, split = ".sign") +
  facet_grid(. ~ .sign) +
  theme_bw(base_size = 12)
ego_gsea@result

show_terms <- c(
  "oxidative phosphorylation",
  "mitochondrial ATP synthesis coupled electron transport",
  "cellular respiration",
  "antigen processing and presentation"
)

res <- as.data.frame(ego_gsea@result)

pick <- res %>%
  filter(Description %in% show_terms) %>%
  arrange(desc(NES))

pick[, c("Description","setSize","NES","pvalue","p.adjust")]

# 方案A：dotplot（按你挑的这些）
dotplot(ego_gsea, showCategory = show_terms) + theme_bw(base_size = 12)

# 方案B：每条来一张 GSEA running curve（更像文章主图）
for (term in show_terms) {
  pid <- res$ID[match(term, res$Description)]
  if (!is.na(pid)) {
    print(gseaplot2(ego_gsea, geneSetID = pid, title = term))
  }
}


# 可视化 ---------------------------------------------------------------------
# 你要画的通路（按代谢逻辑排序）
pathways_keep <- c(
  "Glycolysis / Gluconeogenesis",
  "Pyruvate metabolism",
  "Citrate cycle (TCA cycle)",
  "Oxidative phosphorylation",
  "Pentose phosphate pathway",
  "Fatty acid biosynthesis",
  "Fatty acid elongation",
  "Fatty acid degradation"
)

# 从对象里取实际可用的通路名
all_pathways <- rownames(countexp.Seurat@assays$METABOLISM$score)

# 只保留在对象里确实存在的通路（防止拼写不一致导致空）
input.pathway <- intersect(pathways_keep, all_pathways)

# 如果发现有缺失，打印出来方便你核对通路命名
missing <- setdiff(pathways_keep, all_pathways)
if (length(missing) > 0) {
  message("These pathways were not found in METABOLISM score rownames:\n  ",
          paste(missing, collapse = "\n  "))
}
rm(list = ls())
dev.off()
gc()
