rm(list = ls())
gc()
library(Mime1)
load("QWEN0208_no_GSE14520.Rdata")
load("res_no_GSE14520.Rdata")
pd <- read.table("bulk/TCGA临床信息.txt",sep = "\t",header = T)
colnames(pd)
rownames(pd) <- pd$Patient.ID
rp <- res$riskscore
rp <- rp$`StepCox[forward] + GBM`
rp<- rp[1]
rp <- rp[["TCGA-LIHC"]]
rownames(pd)[1:5]
colnames(pd)[1:5]
rownames(rp)[1:5]
colnames(rp)[1:5]


## 假设对象名就是 pd 和 rp，都是 data.frame

# 1) 从 rp 行名提取病人ID（前12字符）
rp$Patient.ID <- substr(rownames(rp), 1, 12)

# 2) 确保 pd 里也有同名匹配列（用行名生成）
pd$Patient.ID <- rownames(pd)

# 3) 可选：去掉 rp 里列名为 NA 的列（你 colnames(rp) 里有 NA）
rp <- rp[, !is.na(colnames(rp)), drop = FALSE]

# 4) 合并（以 rp 为主，保留全部 rp 样本；pd 没匹配到的会是 NA）
merged <- merge(
  x = rp,
  y = pd,
  by = "Patient.ID",
  all.x = TRUE,
  sort = FALSE
)

# 5) 如果你希望合并后行名仍然是原来的 rp 样本名（带 -01A 那种）
#    merge 后行顺序可能变化，所以用 match 对齐恢复
ord <- match(substr(rownames(rp), 1, 12), merged$Patient.ID)
merged <- merged[ord, , drop = FALSE]
rownames(merged) <- rownames(rp)

colnames(merged)
merged[1:5,]


# 绘制饼图 --------------------------------------------------------------------
library(tidyverse)
library(stringr)
library(scales)

dat <- merged

# -----------------------------
# 1) 风险分组（默认中位数；可替换成你的最佳cutoff）
# -----------------------------
cutoff <- median(dat$RS, na.rm = TRUE)
dat <- dat %>%
  mutate(RiskGroup = ifelse(RS >= cutoff, "High-RS", "Low-RS")) %>%
  mutate(RiskGroup = factor(RiskGroup, levels = c("Low-RS", "High-RS")))

# -----------------------------
# 2) Stage 合并：IIIA/IIIB/IIIC -> III；IVB -> IV；再做 I–II vs III–IV
# -----------------------------
stage_col <- "Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code"

dat <- dat %>%
  mutate(
    Stage_raw = if (stage_col %in% colnames(dat)) as.character(.data[[stage_col]]) else NA_character_,
    Stage_raw = str_to_upper(str_trim(Stage_raw)),
    
    # 先把罗马数字/子分期映射到 1/2/3/4
    Stage_num = case_when(
      is.na(Stage_raw) ~ NA_real_,
      
      # IV 与其子分期（IVB等）优先
      str_detect(Stage_raw, "^STAGE\\s*IV") ~ 4,
      
      # III 与其子分期（IIIA/IIIB/IIIC）
      str_detect(Stage_raw, "^STAGE\\s*III") ~ 3,
      
      # II
      str_detect(Stage_raw, "^STAGE\\s*II\\b") ~ 2,
      
      # I
      str_detect(Stage_raw, "^STAGE\\s*I\\b") ~ 1,
      
      TRUE ~ NA_real_
    ),
    
    Stage_12_34 = case_when(
      is.na(Stage_num) ~ NA_character_,
      Stage_num %in% c(1,2) ~ "Stage I–II",
      Stage_num %in% c(3,4) ~ "Stage III–IV",
      TRUE ~ NA_character_
    )
  )

# -----------------------------
# 3) T stage 合并：T3* -> T3；T4* -> T4
# -----------------------------
t_col <- "American.Joint.Committee.on.Cancer.Tumor.Stage.Code"
dat <- dat %>%
  mutate(
    T_raw = if (t_col %in% colnames(dat)) as.character(.data[[t_col]]) else NA_character_,
    T_raw = str_to_upper(str_trim(T_raw)),
    T_simplified = case_when(
      is.na(T_raw) ~ NA_character_,
      str_detect(T_raw, "^T3") ~ "T3",
      str_detect(T_raw, "^T4") ~ "T4",
      str_detect(T_raw, "^T2") ~ "T2",
      str_detect(T_raw, "^T0") ~ "T1",
      TRUE ~ T_raw
    )
  )
table(dat$T_simplified)
# -----------------------------
# 4) 组装要画的临床特征（NA剔除后重新算比例）
# -----------------------------
df_long <- dat %>%
  transmute(
    RiskGroup,
    Sex = Sex,
    `Histologic Grade` = .data[["Neoplasm.Histologic.Grade"]],
    `AJCC Stage (I–II vs III–IV)` = Stage_12_34,
    `T stage` = T_simplified,
    `N stage` = .data[["Neoplasm.Disease.Lymph.Node.Stage.American.Joint.Committee.on.Cancer.Code"]],
    `M stage` = .data[["American.Joint.Committee.on.Cancer.Metastasis.Stage.Code"]]
  ) %>%
  pivot_longer(cols = -RiskGroup, names_to = "Feature", values_to = "Category") %>%
  filter(!is.na(Category) & Category != "") %>%
  mutate(Category = as.character(Category))

# 同类特征排在一起（列顺序）
feature_order <- c(
  "Sex",
  "Histologic Grade",
  "AJCC Stage (I–II vs III–IV)",
  "T stage",
  "N stage",
  "M stage"
)

df_long <- df_long %>%
  mutate(Feature = factor(Feature, levels = feature_order)) %>%
  filter(!is.na(Feature))

# -----------------------------
# 5) 计数与比例（每个 RiskGroup×Feature 内重新计算）
# -----------------------------
df_sum <- df_long %>%
  group_by(RiskGroup, Feature, Category) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(RiskGroup, Feature) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# -----------------------------
# 6) 每个Feature独立配色（区分度大）
# -----------------------------
df_sum <- df_sum %>% mutate(key = paste(Feature, Category, sep = "||"))

pal_df <- df_sum %>%
  distinct(Feature, Category, key) %>%
  group_by(Feature) %>%
  mutate(col = hcl.colors(n(), palette = "Dark3")) %>%
  ungroup()

pal <- setNames(pal_df$col, pal_df$key)

# -----------------------------
# 7) 标签阈值：只标注占比 >= 8%
# -----------------------------
label_thresh <- 0.08
df_sum <- df_sum %>%
  mutate(label = ifelse(prop >= label_thresh,
                        paste0(Category, "\n", percent(prop, accuracy = 1)),
                        ""))

# -----------------------------
# 8) 画图：完整饼图（去边界线避免“缺口”）+ 横向排布
# -----------------------------
p <- ggplot(df_sum, aes(x = 1, y = n, fill = key)) +
  geom_col(width = 1, color = NA) +   # 关键：完整圆饼，不要白边
  coord_polar(theta = "y") +
  facet_grid(RiskGroup ~ Feature) +
  theme_void(base_size = 12) +
  theme(
    strip.text.x = element_text(face = "bold"),
    strip.text.y = element_text(face = "bold"),
    panel.spacing = unit(1.1, "lines"),
    legend.position = "none"
  ) +
  scale_fill_manual(values = pal) +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 3,
    lineheight = 0.95
  )

print(p)
rm(list = ls())
gc()

