load("Resnet.Rdata")
load("pathdat.Rdata")
dim(resnet_RS_sig)
resnet_RS_sig[1:2,1:2]
rp[1:2,]


library(dplyr)
library(stringr)

# 1) 从 rp 的 ID 提取短样本号：TCGA-XX-XXXX（与 resnet_RS_sig 行名一致）
rp2 <- rp %>%
  mutate(Sample_short = str_extract(ID, "TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}")) %>%
  # 去掉 RS 列
  select(-RS) %>%
  # 同一短样本若出现重复（不同aliquot），取第一条（也可改成 summarise 聚合）
  group_by(Sample_short) %>%
  slice(1) %>%
  ungroup()

# 2) 只保留共有样本，并按 resnet_RS_sig 的顺序对齐
common <- intersect(rownames(resnet_RS_sig), rp2$Sample_short)

resnet_sub <- resnet_RS_sig[common, , drop = FALSE]

rp_sub <- rp2 %>%
  filter(Sample_short %in% common) %>%
  arrange(match(Sample_short, common)) %>%
  as.data.frame()

# 3) 设置行名并删除对齐用列
rownames(rp_sub) <- rp_sub$Sample_short
rp_sub$Sample_short <- NULL

# 4) 合并成新数据框（列拼接）
merged_df <- cbind(rp_sub, resnet_sub)

# 5) 检查行名完全一致
stopifnot(identical(rownames(merged_df), rownames(resnet_sub)))

dim(merged_df)
head(merged_df[, 1:6])

# 机器学习 --------------------------------------------------------------------
library(glmnet)
library(survival)
library(dplyr)

set.seed(2026)

# ---- 1) 构建生存对象 ----
dat <- merged_df
dat$OS.time <- as.numeric(dat$OS.time)
dat$OS <- as.numeric(dat$OS)

# ---- 2) 特征矩阵：只用 resnet* 列 ----
x_cols <- grep("^resnet\\d+$", colnames(dat), value = TRUE)
x <- as.matrix(dat[, x_cols])
y <- Surv(dat$OS.time, dat$OS)

# ---- 3) 7:3 划分 ----
n <- nrow(dat)
idx_train <- sample(seq_len(n), size = floor(0.7 * n))
idx_test  <- setdiff(seq_len(n), idx_train)

x_train <- x[idx_train, , drop = FALSE]
x_test  <- x[idx_test,  , drop = FALSE]
y_train <- y[idx_train]
y_test  <- y[idx_test]

dat_train <- dat[idx_train, , drop = FALSE]
dat_test  <- dat[idx_test,  , drop = FALSE]


cvfit <- cv.glmnet(
  x_train, y_train,
  family = "cox",
  alpha = 1,          # LASSO
  nfolds = 10,
  standardize = TRUE
)

# 系数路径图：每条线是一条特征的系数随 log(lambda) 的变化
pdf("LASSO_coefficient_path.pdf", width = 7, height = 5)
plot(cvfit$glmnet.fit, xvar = "lambda", label = FALSE)
abline(v = log(cvfit$lambda.min), lty = 2)  # lambda.min
abline(v = log(cvfit$lambda.1se), lty = 2)  # lambda.1se
dev.off()

# CV曲线：横轴 log(lambda)，纵轴 CV partial likelihood deviance
pdf("LASSO_CV_curve.pdf", width = 7, height = 5)
plot(cvfit)
abline(v = log(cvfit$lambda.min), lty = 2)  # lambda.min
abline(v = log(cvfit$lambda.1se), lty = 2)  # lambda.1se
dev.off()


use_lambda <- "lambda.min"  # 可改成 "lambda.min"
coef_mat <- coef(cvfit, s = use_lambda)
sel <- which(as.numeric(coef_mat) != 0)
sel_genes <- rownames(coef_mat)[sel]
sel_coef  <- as.numeric(coef_mat)[sel]

sel_table <- data.frame(feature = sel_genes, coef = sel_coef) %>%
  arrange(desc(abs(coef)))

sel_table
length(sel_genes)
# 训练集风险评分
risk_train <- as.numeric(predict(cvfit, newx = x_train, s = use_lambda, type = "link"))
# 验证集风险评分
risk_test  <- as.numeric(predict(cvfit, newx = x_test,  s = use_lambda, type = "link"))

dat_train$RiskScore <- risk_train
dat_test$RiskScore  <- risk_test

# CoxPH 用于计算 C-index（也可用 survConcordance）
fit_train <- coxph(Surv(OS.time, OS) ~ RiskScore, data = dat_train)
fit_test  <- coxph(Surv(OS.time, OS) ~ RiskScore, data = dat_test)

cindex_train <- summary(fit_train)$concordance[1]
cindex_test  <- summary(fit_test)$concordance[1]

cindex_train
cindex_test

df_cv <- data.frame(
  lambda = cvfit$lambda,
  loglambda = log(cvfit$lambda),
  cvm = cvfit$cvm,
  cvup = cvfit$cvup,
  cvlo = cvfit$cvlo
)

p_cv <- ggplot(df_cv, aes(x = loglambda, y = cvm)) +
  geom_line() +
  geom_ribbon(aes(ymin=cvlo, ymax=cvup), alpha = 0.2) +
  theme_classic(base_size = 12) +
  labs(x = "log(Lambda)", y = "CV partial likelihood deviance",
       title = "10-fold cross-validation for LASSO Cox") +
  geom_vline(xintercept = log(cvfit$lambda.min), linetype = 2) +
  geom_vline(xintercept = log(cvfit$lambda.1se), linetype = 2)

p_cv

# surv_cutpoint：maximally selected rank statistics
cut_test <- surv_cutpoint(
  dat_test,
  time = "OS.time",
  event = "OS",
  variables = "RiskScore",
  minprop = 0.25  # 防止切出特别小的组，可改 0.2/0.3
)

ggsurvplot(
survfit(Surv(OS.time, OS) ~ RiskGroup, data = dat_train),
data = dat_train,
pval = TRUE, risk.table = F,
title = "Training set: KM by LASSO Cox risk score")
p_km_train
cut_test$cutpoint
best_cut <- cut_test$cutpoint$cutpoint

dat_test$RiskGroup_best <- ifelse(dat_test$RiskScore > best_cut, "High", "Low")
dat_test$RiskGroup_best <- factor(dat_test$RiskGroup_best, levels = c("Low","High"))
table(dat_test$RiskGroup_best)


fit_km_test <- survfit(Surv(OS.time, OS) ~ RiskGroup_best, data = dat_test)
dat_test$RiskGroup_best <- factor(dat_test$RiskGroup_best, levels = c("Low","High"))
# 固定颜色：Low = 蓝色，High = 红色（你也可以换成你论文既定配色）
cols_risk <- c("Low" = "#1F77B4", "High" = "#D62728")
fit_km_test <- survfit(Surv(OS.time, OS) ~ RiskGroup_best, data = dat_test)

p_km_test <- ggsurvplot(
  fit_km_test,
  data = dat_test,
  pval = TRUE,
  conf.int = FALSE,
  risk.table = F,
  risk.table.height = 0.25,
  legend.title = NULL,
  legend.labs = c("Low", "High"),
  palette = cols_risk,   # ✅ 强制颜色映射
  xlab = "Time",
  ylab = "Overall survival probability",
  title = paste0("Validation set OS (best cutpoint = ", round(best_cut, 3), ")")
)

p_km_test
dat_all <- dat
dat_all$RiskScore <- NA_real_
dat_all$RiskScore[idx_train] <- risk_train
dat_all$RiskScore[idx_test]  <- risk_test
dat_all$Set <- "Test"
dat_all$Set[idx_train] <- "Train"

head(dat_all[, c("ID","OS.time","OS","Set","RiskScore")])
save(dat_all,file = "pathdat.Rdata")
rm(list = ls())
gc()
dev.off()
