rm(list = ls())
gc()
library(Mime1)
load("QWEN0208_no_GSE14520.Rdata")
load("res_no_GSE14520.Rdata")
pd <- read.csv("resnet50_features.csv",header = T,row.names = 1)
dim(pd)
rp <- res$riskscore$`StepCox[forward] + GBM`
rp <- res$riskscore$`StepCox[forward] + GBM`$`TCGA-LIHC`
dim(rp)
rp[1:4,1:4]
rownames(pd)[1:20]
colnames(pd)[1:5]
library(dplyr)
library(stringr)

# ---- 0) 确保 pd 是 data.frame（如果是 matrix 也可转） ----
pd <- as.data.frame(pd)

# ---- 1) 重命名特征列：resnet0 ~ resnet(2047)（按列数自动） ----
colnames(pd) <- paste0("resnet", seq_len(ncol(pd)) - 1)

# ---- 2) 从行名路径提取 Sample 和 tile ----
pd2 <- pd %>%
  mutate(
    File = rownames(pd),
    Sample = str_extract(File, "TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}"),
    tile   = str_extract(File, "tile_\\d+")
  )

# 可选：检查是否有提取失败
table(is.na(pd2$Sample)); table(is.na(pd2$tile))

# ---- 3) 按样本取平均（所有 resnet* 数值列求均值） ----
resnet2 <- pd2 %>%
  group_by(Sample) %>%
  summarise(across(starts_with("resnet"), ~ mean(.x, na.rm = TRUE))) %>%
  ungroup()

# ---- 4) 转成“样本为行名”的矩阵样式 ----
resnet2 <- as.data.frame(resnet2)
rownames(resnet2) <- resnet2$Sample
resnet2$Sample <- NULL

# resnet2 就是：n_samples x 2048 的样本级特征矩阵
dim(resnet2)
head(rownames(resnet2))
head(colnames(resnet2)[1:10])


tile_count <- pd2 %>% count(Sample, name = "n_tiles")
tile_count[order(tile_count$n_tiles), ]

rp <- res$riskscore$`StepCox[forward] + GBM`$`TCGA-LIHC`

dim(rp)
rp[1:4,1:4]


# 1) 从 rp 的行名/ID 提取到与 resnet2 一致的短样本号：TCGA-XX-XXXX
rp2 <- rp %>%
  mutate(
    Sample_short = str_extract(ID, "TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}")
  ) %>%
  # 同一个短样本可能出现重复（不同 aliquot），保留第一个或取平均都行
  group_by(Sample_short) %>%
  summarise(RS = mean(RS, na.rm = TRUE), .groups = "drop")

# 2) 按 resnet2 的行名顺序对齐 RS
idx <- match(rownames(resnet2), rp2$Sample_short)
RS_vec <- rp2$RS[idx]

# 3) 加到 resnet2（新列叫 RS）
resnet2_RS <- resnet2
resnet2_RS$RS <- RS_vec

# 4) 检查匹配情况
cat("resnet2 samples:", nrow(resnet2), "\n")
cat("matched RS:", sum(!is.na(resnet2_RS$RS)), "\n")
cat("missing RS:", sum(is.na(resnet2_RS$RS)), "\n")

# 看看没匹配上的样本（如果有）
missing_samples <- rownames(resnet2_RS)[is.na(resnet2_RS$RS)]
head(missing_samples, 20)

# 可选：确保没有顺序错位（抽查）
head(resnet2_RS[, c("RS")])

resnet2_RS

# 相关性分析P小于0.05和R绝对值大于0.2 -------------------------------------------------------------------

get_RS_correlated_features <- function(df,
                                       rs_col = "RS",
                                       method = c("pearson", "spearman"),
                                       p_cutoff = 0.05,
                                       r_cutoff = 0.2,
                                       adjust_method = "BH",
                                       use_fdr = FALSE) {
  method <- match.arg(method)
  
  if (!rs_col %in% colnames(df)) {
    stop("rs_col not found in df: ", rs_col)
  }
  if (!is.numeric(df[[rs_col]])) {
    stop("RS column must be numeric.")
  }
  
  # 只取数值列（去掉RS本身）
  num_cols <- names(df)[sapply(df, is.numeric)]
  feat_cols <- setdiff(num_cols, rs_col)
  
  # 对每个特征做相关检验
  rs <- df[[rs_col]]
  res_list <- lapply(feat_cols, function(coln) {
    x <- df[[coln]]
    ok <- is.finite(rs) & is.finite(x)
    if (sum(ok) < 3) {
      return(data.frame(feature = coln, n = sum(ok), r = NA_real_, p = NA_real_))
    }
    ct <- suppressWarnings(cor.test(x[ok], rs[ok], method = method))
    data.frame(
      feature = coln,
      n = sum(ok),
      r = unname(ct$estimate),
      p = ct$p.value
    )
  })
  
  res <- do.call(rbind, res_list)
  res$fdr <- p.adjust(res$p, method = adjust_method)
  
  # 过滤阈值：p 或 fdr
  p_used <- if (use_fdr) res$fdr else res$p
  hits <- res$feature[!is.na(res$r) & !is.na(p_used) &
                        (p_used < p_cutoff) & (abs(res$r) > r_cutoff)]
  
  # 输出子数据框：命中特征 + RS
  filtered_df <- df[, c(hits, rs_col), drop = FALSE]
  
  # 按显著性和效应量排序
  res <- res[order(ifelse(is.na(p_used), Inf, p_used), -abs(res$r)), ]
  
  list(
    hit_features = hits,
    result_table = res,
    filtered_df = filtered_df
  )
}
out <- get_RS_correlated_features(resnet2_RS,
                                  rs_col = "RS",
                                  method = "pearson",
                                  p_cutoff = 0.05,
                                  r_cutoff = 0.2)

length(out$hit_features)
head(out$hit_features, 20)
head(out$result_table, 20)

# 得到只包含“与RS相关的resnet特征 + RS”的数据
resnet_RS_sig <- out$filtered_df
dim(resnet_RS_sig)


# 可视化 ---------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(stringr)

# 取显著特征表
rt <- out$result_table

# 只保留命中（你已经筛过了，也可以再保险一次）
rt_sig <- rt %>% filter(feature %in% out$hit_features)

# 选网络里展示的数量（建议 15+15，太多会乱）
n_show <- 15
pos <- rt_sig %>% filter(r > 0) %>% arrange(p) %>% slice(1:n_show)
neg <- rt_sig %>% filter(r < 0) %>% arrange(p) %>% slice(1:n_show)

rt_net <- bind_rows(pos, neg) %>%
  mutate(sign = ifelse(r > 0, "Positive", "Negative"),
         weight = abs(r))

# --- 构造 igraph ---
if (!requireNamespace("igraph", quietly = TRUE)) install.packages("igraph")
if (!requireNamespace("ggraph", quietly = TRUE)) install.packages("ggraph")
library(igraph)
library(ggraph)

nodes <- data.frame(
  name = c("RS", rt_net$feature),
  type = c("RS", rep("resnet", nrow(rt_net))),
  sign = c("RS", rt_net$sign),
  size = c(max(rt_net$weight) * 1.2, rt_net$weight),
  stringsAsFactors = FALSE
)

edges <- data.frame(
  from = "RS",
  to = rt_net$feature,
  sign = rt_net$sign,
  weight = rt_net$weight,
  stringsAsFactors = FALSE
)

g <- igraph::graph_from_data_frame(edges, directed = FALSE, vertices = nodes)

# --- 手工布局：左负右正，中间RS ---
V(g)$x <- 0
V(g)$y <- 0
V(g)$x[V(g)$name %in% neg$feature] <- -1
V(g)$x[V(g)$name %in% pos$feature] <-  1

# 给左右节点一个纵向排序
V(g)$y[V(g)$name %in% neg$feature] <- seq(from = 1, to = -1, length.out = nrow(neg))
V(g)$y[V(g)$name %in% pos$feature] <- seq(from = 1, to = -1, length.out = nrow(pos))

layout_df <- data.frame(name = V(g)$name, x = V(g)$x, y = V(g)$y)

p_net <- ggraph(g, layout = "manual", x = layout_df$x, y = layout_df$y) +
  geom_edge_link(aes(width = weight, linetype = sign), alpha = 0.7) +
  geom_node_point(aes(size = size, shape = type)) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  theme_void(base_size = 12) +
  labs(title = "Correlation network: RS-associated ResNet features") +
  theme(legend.position = "none")

p_net


# 选前30个最显著特征（或按 abs(r) 最大）
topK <- 30
top_feat <- rt_sig %>% arrange(p) %>% slice(1:topK) %>% pull(feature)

# 提取矩阵：样本 x 特征
mat <- resnet2_RS[, top_feat, drop = FALSE]

# 按 RS 从低到高排序
ord <- order(resnet2_RS$RS, na.last = NA)
mat <- mat[ord, , drop = FALSE]
rs_sorted <- resnet2_RS$RS[ord]

# z-score 标准化每个特征列（更适合热图比较）
mat_z <- scale(mat)

if (!requireNamespace("pheatmap", quietly = TRUE)) install.packages("pheatmap")
library(pheatmap)

# 注释条：RS（连续变量用色带）
ann <- data.frame(RS = rs_sorted)
rownames(ann) <- rownames(mat_z)

pheatmap(mat_z,
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         show_rownames = FALSE,
         main = "Top RS-associated ResNet features (z-score; samples ordered by RS)",
         annotation_row = ann)



library(dplyr)

# 1) 取显著关联的resnet特征结果
rt_sig <- out$result_table %>%
  filter(feature %in% out$hit_features) %>%
  mutate(
    source = "RS",
    target = feature,
    interaction = "correlation",
    sign = ifelse(r >= 0, "positive", "negative"),
    abs_r = abs(r)
  )

# 2) Edge table（边表）
edges <- rt_sig %>%
  select(source, target, interaction, sign, r, p, fdr, abs_r, n)

# 3) Node table（点表）
nodes_resnet <- rt_sig %>%
  transmute(
    id = target,
    node_type = "resnet",
    sign = sign,
    r = r,
    p = p,
    fdr = fdr,
    abs_r = abs_r,
    n = n
  )

node_RS <- data.frame(
  id = "RS",
  node_type = "RS",
  sign = NA,
  r = NA,
  p = NA,
  fdr = NA,
  abs_r = NA,
  n = NA
)

nodes <- bind_rows(node_RS, nodes_resnet) %>%
  distinct(id, .keep_all = TRUE)

# 4) 导出（推荐 TSV，Cytoscape导入更稳）
write.table(edges, file = "Cytoscape_RS_resnet_edges.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(nodes, file = "Cytoscape_RS_resnet_nodes.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

# 5) 快速检查
dim(edges); dim(nodes)
head(edges); head(nodes)
resnet_RS_sig
save(resnet2,resnet_RS_sig,rp,file = "Resnet.Rdata") #rp中包含临床数据
dim(rp)
rp[1:4,1:4]


