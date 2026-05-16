#!/usr/bin/env Rscript
# TCGA-LIHC 反卷积分析脚本 - 大类细胞类型版本
# 使用BayesPrism方法对TCGA-LIHC bulk RNA-seq数据进行细胞类型反卷积
# 小亚群合并为大类，但肿瘤细胞OS-high/OS-low分组保持分开

# 加载必要的包 ------------------------------------------------------------------
cat("加载必要的R包...\n")
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(SCP)
library(BayesPrism)

# 加载数据 --------------------------------------------------------------------
cat("加载单细胞参考数据...\n")
load("注释sceall.Rdata")

cat("加载TCGA-LIHC bulk数据...\n")
load("bulk/TCGA-LIHC.Rdata")

cat("TCGA-LIHC数据概览：\n")
cat("  表达矩阵维度:", dim(exp), "\n")
cat("  样本数:", ncol(exp), "\n")
cat("  基因数:", nrow(exp), "\n")
cat("  生存数据样本数:", nrow(surv), "\n")

# 清理内存
gc()

# 单细胞数据预处理 --------------------------------------------------------------
cat("\n=== 单细胞数据预处理 ===\n")
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

# 检查是否有未分类的细胞类型
original_types <- unique(sce.noNE$final_celltype)
merged_types <- unique(sce.noNE$major_celltype)
cat("\n原始细胞类型数:", length(original_types), "\n")
cat("合并后细胞类型数:", length(merged_types), "\n")
cat("\n合并后的细胞类型列表:\n")
print(sort(merged_types))

# 设置细胞类型标识
Idents(sce.noNE) <- "major_celltype"

# 进行细胞抽样（每个细胞类型保留30%，减少计算负担）
cat("\n对单细胞数据进行抽样（保留30%细胞）...\n")
cells.keep <- unlist(lapply(levels(Idents(sce.noNE)), function(ct){
  cells <- WhichCells(sce.noNE, idents = ct)
  n_keep <- max(1, floor(0.30 * length(cells)))  # 每类至少留1个
  sample(cells, n_keep)
}))

sce.noNE.30 <- subset(sce.noNE, cells = cells.keep)
cat("抽样完成，各细胞类型数量：\n")
print(table(sce.noNE.30$major_celltype))

# 准备单细胞参考数据
cat("\n准备单细胞参考矩阵...\n")
# 检查Seurat版本并提取counts
if (.hasSlot(sce.noNE.30@assays$RNA, "layers")) {
  # Seurat V5
  cat("检测到Seurat V5格式，使用layers提取counts...\n")
  ref <- sce.noNE.30@assays$RNA$counts
} else {
  # Seurat V3
  cat("检测到Seurat V3格式，使用@counts提取counts...\n")
  ref <- sce.noNE.30@assays$RNA@counts
}
ref <- as.matrix(ref)
ref <- as.data.frame(ref)
ref <- t(ref)  # 转置为 cell x gene
gc()

sc.dat <- ref
cell.type.labels  <- sce.noNE.30$major_celltype
cell.state.labels <- sce.noNE.30$major_celltype  # 使用相同的细胞类型标签

# 准备bulk数据
cat("准备TCGA-LIHC bulk数据...\n")
bk.dat <- t(exp)  # 转置为 sample x gene

cat("  单细胞参考矩阵维度:", dim(sc.dat), "(cells x genes)\n")
cat("  Bulk数据矩阵维度:", dim(bk.dat), "(samples x genes)\n")

# 质量控制 ---------------------------------------------------------------------
cat("\n=== 质量控制分析 ===\n")

# 细胞类型相关性分析
cat("生成细胞类型相关性图...\n")
plot.cor.phi(input=ref,
             input.labels=cell.type.labels,
             title="cell type correlation",
             pdf.prefix="TCGA-LIHC_major_cell.type.correlation",
             cexRow=0.8, cexCol=0.8)

# 单细胞数据离群值检测
cat("检测单细胞数据离群值...\n")
sc.stat <- plot.scRNA.outlier(
  input=sc.dat,
  cell.type.labels=cell.type.labels,
  species="hs",
  return.raw=TRUE,
  pdf.prefix="TCGA-LIHC_major_sc.stat"
)

# Bulk数据离群值检测
cat("检测bulk数据离群值...\n")
bk.stat <- plot.bulk.outlier(
  bulk.input=bk.dat,
  sc.input=ref,
  cell.type.labels=cell.type.labels,
  species="hs",
  return.raw=TRUE,
  pdf.prefix="TCGA-LIHC_major_bk.stat"
)

# 数据清理与过滤 ----------------------------------------------------------------
cat("\n=== 数据清理与过滤 ===\n")

# 过滤低表达基因和特殊基因
cat("过滤低表达基因和特殊基因（Rb, chrM, chrX, chrY等）...\n")
sc.dat.filtered <- cleanup.genes(
  input=sc.dat,
  input.type="count.matrix",
  species="hs",
  gene.group=c("Rb","other_Rb","chrM","chrX","chrY","Mrp","hb","MALAT1"),
  exp.cells=5
)

cat("过滤后单细胞数据维度:", dim(sc.dat.filtered), "\n")

# Bulk vs Single-cell对比
cat("生成bulk vs single-cell对比图...\n")
plot.bulk.vs.sc(
  sc.input=sc.dat.filtered,
  bulk.input=bk.dat,
  pdf.prefix="TCGA-LIHC_major_bulk.vs.sc"
)

# 标记基因选择 ------------------------------------------------------------------
cat("\n=== 计算差异表达统计 ===\n")
cat("这一步可能需要较长时间，使用20个核心...\n")

diff.exp.stat <- get.exp.stat(
  sc.dat=sc.dat[,colSums(sc.dat>0)>3],  # 过滤基因以减少内存使用
  cell.type.labels=cell.type.labels,
  cell.state.labels=cell.state.labels,
  pseudo.count=0.1,  # 10x数据使用0.1
  cell.count.cutoff=50,
  n.cores=20
)

cat("选择标记基因（pval < 0.01, lfc > 0.1）...\n")
sc.dat.filtered.pc.sig <- select.marker(
  sc.dat=sc.dat.filtered,
  stat=diff.exp.stat,
  pval.max=0.01,
  lfc.min=0.1
)

cat("选择的标记基因数:", length(sc.dat.filtered.pc.sig), "\n")

# BayesPrism反卷积 -------------------------------------------------------------
cat("\n=== 构建BayesPrism模型 ===\n")

myPrism <- new.prism(
  reference=sc.dat.filtered,
  mixture=bk.dat,
  input.type="count.matrix",
  cell.type.labels=cell.type.labels,
  cell.state.labels=cell.state.labels,
  key=NULL,
  outlier.cut=0.01,
  outlier.fraction=0.1
)

cat("\n=== 开始运行BayesPrism反卷积 ===\n")
cat("使用20个核心，这可能需要较长时间...\n")
start_time <- Sys.time()

bp.res <- run.prism(prism=myPrism, n.cores=20)

end_time <- Sys.time()
cat("反卷积完成！耗时:", difftime(end_time, start_time, units="mins"), "分钟\n")

# 结果保存 ---------------------------------------------------------------------
cat("\n=== 保存结果 ===\n")

# 提取细胞类型比例
theta <- bp.res@posterior.initial.cellState@theta
cat("细胞类型比例矩阵维度:", dim(theta), "\n")
cat("样本数:", nrow(theta), "\n")
cat("细胞类型数:", ncol(theta), "\n")

# 保存完整结果
cat("保存完整BayesPrism结果到 TCGA-LIHC_major_bp.res.Rdata...\n")
save(bp.res, file="TCGA-LIHC_major_bp.res.Rdata")

# 保存细胞类型比例矩阵
cat("保存细胞类型比例矩阵到 TCGA-LIHC_major_cell_fractions.csv...\n")
write.csv(theta, file="TCGA-LIHC_major_cell_fractions.csv", row.names=TRUE)

# 合并生存数据
cat("合并细胞类型比例与生存数据...\n")
theta_df <- as.data.frame(theta)
theta_df$sample <- rownames(theta_df)

# 匹配生存数据
surv$sample <- surv$sample
merged_data <- merge(theta_df, surv, by="sample", all.x=TRUE)

cat("保存合并数据到 TCGA-LIHC_major_deconv_with_survival.csv...\n")
write.csv(merged_data, file="TCGA-LIHC_major_deconv_with_survival.csv", row.names=FALSE)

# 结果摘要 ---------------------------------------------------------------------
cat("\n=== 反卷积结果摘要 ===\n")
cat("细胞类型列表：\n")
print(colnames(theta))

cat("\n各细胞类型平均比例：\n")
print(sort(colMeans(theta), decreasing=TRUE))

cat("\n前5个样本的细胞类型比例：\n")
print(head(theta, 5))

# 生成可视化 -------------------------------------------------------------------
cat("\n=== 生成可视化图表 ===\n")

# 细胞类型比例箱线图
pdf("TCGA-LIHC_major_cell_fraction_boxplot.pdf", width=10, height=6)
theta_long <- reshape2::melt(theta)
colnames(theta_long) <- c("Sample", "CellType", "Fraction")
p1 <- ggplot(theta_long, aes(x=reorder(CellType, -Fraction, FUN=median),
                              y=Fraction, fill=CellType)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=10),
        legend.position="none") +
  labs(title="TCGA-LIHC Major Cell Type Fractions",
       x="Cell Type",
       y="Fraction")
print(p1)
dev.off()

# 细胞类型比例热图（所有样本）
pdf("TCGA-LIHC_major_cell_fraction_heatmap.pdf", width=12, height=6)
library(pheatmap)
pheatmap(t(theta),
         cluster_rows=TRUE,
         cluster_cols=TRUE,
         show_colnames=FALSE,
         scale="row",
         main="TCGA-LIHC Major Cell Type Fractions (All samples)",
         fontsize_row=10)
dev.off()

# 肿瘤细胞OS分组比较
if ("Tumor_OS-high" %in% colnames(theta) && "Tumor_OS-low" %in% colnames(theta)) {
  pdf("TCGA-LIHC_tumor_OS_comparison.pdf", width=8, height=6)

  # 计算肿瘤细胞总比例
  tumor_data <- data.frame(
    OS_high = theta[, "Tumor_OS-high"],
    OS_low = theta[, "Tumor_OS-low"]
  )
  tumor_data$Total_tumor <- tumor_data$OS_high + tumor_data$OS_low
  tumor_data$OS_high_ratio <- tumor_data$OS_high / (tumor_data$Total_tumor + 1e-10)

  # OS-high vs OS-low散点图
  p2 <- ggplot(tumor_data, aes(x=OS_low, y=OS_high)) +
    geom_point(alpha=0.5) +
    geom_abline(slope=1, intercept=0, linetype="dashed", color="red") +
    theme_bw() +
    labs(title="Tumor OS-high vs OS-low Fractions",
         x="Tumor OS-low Fraction",
         y="Tumor OS-high Fraction")
  print(p2)

  # OS-high比例分布
  p3 <- ggplot(tumor_data, aes(x=OS_high_ratio)) +
    geom_histogram(bins=30, fill="steelblue", alpha=0.7) +
    theme_bw() +
    labs(title="Distribution of OS-high Ratio in Tumors",
         x="OS-high / Total Tumor",
         y="Count")
  print(p3)

  dev.off()

  cat("\n肿瘤细胞OS分组统计：\n")
  cat("  OS-high平均比例:", mean(theta[, "Tumor_OS-high"]), "\n")
  cat("  OS-low平均比例:", mean(theta[, "Tumor_OS-low"]), "\n")
  cat("  OS-high/(OS-high+OS-low)平均值:", mean(tumor_data$OS_high_ratio), "\n")
}

cat("\n=== 全部完成！===\n")
cat("输出文件：\n")
cat("  1. TCGA-LIHC_major_bp.res.Rdata - 完整BayesPrism结果对象\n")
cat("  2. TCGA-LIHC_major_cell_fractions.csv - 细胞类型比例矩阵\n")
cat("  3. TCGA-LIHC_major_deconv_with_survival.csv - 合并生存数据的结果\n")
cat("  4. TCGA-LIHC_major_cell_fraction_boxplot.pdf - 细胞类型比例箱线图\n")
cat("  5. TCGA-LIHC_major_cell_fraction_heatmap.pdf - 细胞类型比例热图\n")
cat("  6. TCGA-LIHC_tumor_OS_comparison.pdf - 肿瘤OS分组比较图\n")
cat("  7. TCGA-LIHC_major_*.pdf - 各种质量控制图\n")

# 清理
gc()
cat("\n分析完成！\n")
