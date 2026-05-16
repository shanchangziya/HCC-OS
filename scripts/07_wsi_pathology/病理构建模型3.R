rm(list = ls())
gc()
load("Resnet.Rdata")
load("pathdat.Rdata")
load("bulk/exp1surv1.rdata")
dim(exp1)
exp1[1:4,]
dim(dat_all)
dat_all[1:4,]
library(dplyr)

# 1) 只保留 dat_all 中用于匹配和分组的两列：ID + RiskScore
risk_df <- dat_all %>%
  dplyr::select(ID, RiskScore) %>%
  dplyr::filter(!is.na(ID) & !is.na(RiskScore)) %>%
  dplyr::distinct(ID, .keep_all = TRUE)

# 2) 取 exp1 与 risk_df 的交集样本（以 TCGA-..-01A 为键）
common_ids <- intersect(colnames(exp1), risk_df$ID)

length(common_ids)          # 交集样本数
head(common_ids)

# 3) 对齐顺序：让 risk_df 的顺序与 exp1 的列顺序一致
risk_df2 <- risk_df %>%
  dplyr::filter(ID %in% common_ids) %>%
  dplyr::arrange(match(ID, common_ids))

# 4) 子集表达矩阵（只保留交集样本，并按同顺序排列）
exp1_sub <- exp1[, common_ids, drop = FALSE]

# 5) 按 RiskScore 中位数分组（High / Low）
cutoff <- median(risk_df2$RiskScore, na.rm = TRUE)

risk_group <- ifelse(risk_df2$RiskScore >= cutoff, "High", "Low")
risk_group <- factor(risk_group, levels = c("Low", "High"))
names(risk_group) <- risk_df2$ID

table(risk_group)

# 6) 如果你想把分组写进一个 meta 表（常用于差异分析）
meta_risk <- data.frame(
  Sample = risk_df2$ID,
  RiskScore = risk_df2$RiskScore,
  RiskGroup = risk_group,
  stringsAsFactors = FALSE
)
rownames(meta_risk) <- meta_risk$Sample

head(meta_risk)

# 7)（可选）直接得到 High/Low 两个表达矩阵
exp1_high <- exp1_sub[, meta_risk$RiskGroup == "High", drop = FALSE]
exp1_low  <- exp1_sub[, meta_risk$RiskGroup == "Low",  drop = FALSE]

dim(exp1_high); dim(exp1_low)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

# exp1_sub: genes x samples
# meta_risk: rownames = sample, RiskGroup = Low/High

stopifnot(all(colnames(exp1_sub) %in% rownames(meta_risk)))
meta_risk <- meta_risk[colnames(exp1_sub), , drop = FALSE]
meta_risk$RiskGroup <- factor(meta_risk$RiskGroup, levels = c("Low","High"))

# ---- 判断数据类型：是否近似整数counts ----
is_counts <- all(abs(exp1_sub - round(exp1_sub)) < 1e-6, na.rm = TRUE)

if (is_counts) {
  suppressPackageStartupMessages({
    library(DESeq2)
  })
  dds <- DESeqDataSetFromMatrix(
    countData = round(exp1_sub),
    colData = meta_risk,
    design = ~ RiskGroup
  )
  dds <- dds[rowSums(counts(dds)) >= 10, ]
  dds <- DESeq(dds)
  
  res <- results(dds, contrast = c("RiskGroup", "High", "Low"))
  deg <- as.data.frame(res) %>%
    tibble::rownames_to_column("gene") %>%
    rename(logFC = log2FoldChange, p = pvalue, fdr = padj) %>%
    arrange(fdr)
} else {
  suppressPackageStartupMessages({
    library(limma)
  })
  expr <- as.matrix(exp1_sub)
  design <- model.matrix(~ RiskGroup, data = meta_risk)
  fit <- lmFit(expr, design)
  fit <- eBayes(fit, trend = TRUE)
  tt <- topTable(fit, coef = "RiskGroupHigh", number = Inf, sort.by = "P")
  deg <- tt %>%
    tibble::rownames_to_column("gene") %>%
    rename(logFC = logFC, p = P.Value, fdr = adj.P.Val) %>%
    arrange(fdr)
}

# 保存差异结果
write.csv(deg, "DEG_High_vs_Low_Risk.csv", row.names = FALSE)

# 火山图（论文图）
deg$neglog10p <- -log10(pmax(deg$p, 1e-300))
deg$sig <- with(deg, fdr < 0.05 & abs(logFC) > 0.5)

p_volcano <- ggplot(deg, aes(x = logFC, y = neglog10p)) +
  geom_point(aes(alpha = sig), size = 1.1) +
  theme_classic(base_size = 12) +
  labs(x = "log2FC (High vs Low)", y = "-log10(P)", title = "DEGs: High-risk vs Low-risk") +
  guides(alpha = "none")
ggsave("Fig_volcano_DEG.pdf", p_volcano, width = 6.5, height = 5, useDingbats = FALSE)


library(dplyr)

# 设定阈值
fdr_cut <- 0.05
lfc_cut <- 0.5

deg_labeled <- deg %>%
  mutate(
    Direction = case_when(
      !is.na(fdr) & fdr < fdr_cut & logFC >  lfc_cut ~ "Up",
      !is.na(fdr) & fdr < fdr_cut & logFC < -lfc_cut ~ "Down",
      TRUE ~ "NS"
    )
  )

# 看看数量
table(deg_labeled$Direction)

# 导出（网页常用）
write.csv(deg_labeled, "DEG_High_vs_Low_with_direction.csv", row.names = FALSE)

# 如果你只想导出 up/down 基因列表（更适合网页直接粘贴）
write.table(deg_labeled$gene[deg_labeled$Direction == "Up"],
            "DEG_Up_genes.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(deg_labeled$gene[deg_labeled$Direction == "Down"],
            "DEG_Down_genes.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
rm(list = ls())
gc()
dev.off()
