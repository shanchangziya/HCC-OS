library(SCP)
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})
tcos <- sce.all
save(tcos,file = "tcos.Rdata")
load("tcos.Rdata")
# 寻找高低切割点 -------------------------------------------------------------------
epi$OS_score_mean_z
FeatureDimPlot(epi,features = "OS_score_mean_z")
epi <- tcos
thr_epi <- as.numeric(quantile(epi$OS_score_mean_z, probs = 0.80, na.rm = TRUE))
epi$OS_group <- ifelse(epi$OS_score_mean_z >= thr_epi, "OS-high (top20%)", "OS-low (bottom80%)")
epi$OS_group <- factor(epi$OS_group, levels = c("OS-low (bottom80%)", "OS-high (top20%)"))
## 2) 图1：密度分布 + 阈值线（好看且最直观）
df <- data.frame(OS = epi$OS_score_mean_z, Group = epi$OS_group)
p_density <- ggplot(df, aes(x = OS, fill = Group)) +
  geom_density(alpha = 0.45, linewidth = 0.8) +
  geom_vline(xintercept = thr_epi, linetype = 2, linewidth = 0.9) +
  annotate("text", x = thr_epi, y = Inf, label = paste0("Cutoff = ", round(thr_epi, 3)),
           vjust = 1.3, hjust = -0.05, size = 4) +
  labs(x = "OS composite score (z-mean)", y = "Density", fill = NULL,
       title = "Epithelial cells: OS score stratification") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "top")
p_density
tcos <- sce.all
save(tcos,file = "tcos.Rdata")
sce <- RunSlingshot(tcos,group.by = "OS_group_valley")

# 干性分析 --------------------------------------------------------------------
library(CytoTRACE2) 
expression_data <- GetAssayData(sce, assay = "RNA", slot = "counts")
expression_data <- as.matrix(expression_data)
expression_data[1:4,1:4]
cytotrace2_result <- cytotrace2(expression_data,species = "human")
identical(rownames(sce@meta.data),rownames(cytotrace2_result))
sce <- AddMetaData(sce, metadata = cytotrace2_result)
sce$CytoTRACE2_Relative
FeatureDimPlot(sce,features = "CytoTRACE2_Relative")
# NMF分析 ---------------------------------------------------------------------
library(GeneNMF)
library(RColorBrewer)
library(viridis)
library(UCell)
seu <- tcos
DefaultAssay(seu) <- "RNA"
seu$orig.ident
seu.list <- SplitObject(seu, split.by = "orig.ident")
geneNMF.programs <- multiNMF(seu.list, assay="RNA", k=4:9, min.exp = 0.05)
geneNMF.metaprograms <- getMetaPrograms(geneNMF.programs,
                                        metric = "cosine",
                                        specificity.weight = 3,
                                        weight.explained = 0.3,
                                        nMP=4)
plotMetaPrograms(geneNMF.metaprograms)
geneNMF.metaprograms$metaprograms.metrics
geneNMF.metaprograms$metaprograms.genes
mp.genes <- geneNMF.metaprograms$metaprograms.genes
seu <- AddModuleScore_UCell(seu, features = mp.genes, ncores=4, name = "")
seu$MP1
FeaturePlot(seu,features = "MP1")
FeaturePlot(seu,features = "MP4")
save(seu,geneNMF.metaprograms,file = "05osnmf.Rdata")
#load("05osnmf.Rdata")
anno_colors <- brewer.pal(n=10, name="Paired")
names(anno_colors) <- names(geneNMF.metaprograms$metaprograms.genes)
plotMetaPrograms(geneNMF.metaprograms,
                 palette=viridis(100, option = "H", direction = 1),
                 similarity.cutoff = c(0,1),showtree = T)
seu$OS_group_valley

# 计算 Spearman 相关（用于图上标注）
df <- FetchData(seu, vars = c("MP4", "OS_score_mean_z", "OS_group_valley"))
ct <- cor.test(df$MP4, df$OS_score_mean_z, method = "spearman")

p <- FeatureScatter(
  seu,
  feature1 = "MP4",
  feature2 = "OS_score_mean_z",
  group.by = "OS_group_valley",
  pt.size = 0.25
) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8, color = "grey30") +
  labs(
    x = "MP4 (NMF program score)",
    y = "OS composite score (z-mean)",
    color = NULL
  ) +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = paste0("Spearman ρ = ", round(ct$estimate, 2),
                   "\nP = ", format.pval(ct$p.value, digits = 2, eps = 1e-300)),
    hjust = 1.1, vjust = 1.2, size = 4
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )
p
dev.off()

