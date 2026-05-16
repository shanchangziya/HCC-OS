load("TCGA-LIHC_major_bp.res.Rdata")
bp.res
os <- read.csv("TCGA-LIHC_major_deconv_with_survival.csv")
colnames(os)
meta <- read.table("TCGA临床信息.txt",sep="\t",header = T)
meta[1:2,]
dat[1:2,]
library(dplyr)

# 1) dat: sample 截取前12位作为 Patient.ID
dat2 <- dat %>%
  mutate(Patient.ID = substr(sample, 1, 12))

# 2) meta: 第一列重命名为 Patient.ID（若已是 Patient.ID 则无影响）
meta2 <- meta
colnames(meta2)[1] <- "Patient.ID"

# 3) 以 meta 为主表，把 dat 拼上来（保留所有 meta 行）
meta_dat <- meta2 %>%
  left_join(dat2, by = "Patient.ID")

str(meta_dat)
# 以 OS 评分（这里用 Tumor_OS.high 这一列）按中位数分组：High / Low
meta_dat <- meta_dat %>%
  mutate(
    OS_score_median_group = if_else(
      Tumor_OS.high >= median(Tumor_OS.high, na.rm = TRUE),
      "High", "Low"
    ) %>% factor(levels = c("Low","High"))
  )

# 可选：看一下分组人数
table(meta_dat$OS_score_median_group, useNA = "ifany")
table(meta_dat$OS_score_median_group, meta_dat$Sex)
table(meta_dat$OS_score_median_group, meta_dat$Neoplasm.Histologic.Grade)
table(meta_dat$OS_score_median_group, meta_dat$Sex)
write.csv(meta_dat,file = "metaos.csv")


suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
  library(survminer)
})

# -----------------------------
# 1) 清理数据：去掉关键列为 NA 的样本
# -----------------------------
dat <- os %>%
  dplyr::select(sample, Tumor_OS.high, OS.time, OS) %>%
  dplyr::filter(
    !is.na(sample),
    !is.na(Tumor_OS.high),
    !is.na(OS.time),
    !is.na(OS)
  )

# 确保类型正确
dat <- dat %>%
  mutate(
    Tumor_OS.high = as.numeric(Tumor_OS.high),
    OS.time = as.numeric(OS.time),
    OS = as.numeric(OS)   # 0/1
  )

cat("After removing NA, n =", nrow(dat), "\n")


dat <- os %>%
  dplyr::select(sample, Tumor_OS.high, OS.time, OS) %>%
  filter(!is.na(sample), !is.na(Tumor_OS.high), !is.na(OS.time), !is.na(OS)) %>%
  mutate(
    Tumor_OS.high = as.numeric(Tumor_OS.high),
    OS.time = as.numeric(OS.time),
    OS = as.numeric(OS)
  ) %>%
  group_by(sample) %>%
  slice_max(order_by = OS.time, n = 1, with_ties = FALSE) %>%  # 去重
  ungroup()

cat("After NA filter + dedup, n =", nrow(dat), "\n")


# -----------------------------
# 2) 最佳截断值（限制最小组比例，避免极端不平衡）
# -----------------------------
cut <- surv_cutpoint(
  dat,
  time = "OS.time",
  event = "OS",
  variables = "Tumor_OS.high",
  minprop = 0.20   # 可改 0.10/0.25；建议别太小
)

cp <- cut$cutpoint$cutpoint
cat("Best cutpoint (Tumor_OS.high) =", cp, "\n")

dat$TumorHigh_group <- ifelse(dat$Tumor_OS.high > cp, "High", "Low")
dat$TumorHigh_group <- factor(dat$TumorHigh_group, levels = c("Low", "High"))
print(table(dat$TumorHigh_group))

# -----------------------------
# 3) KM 生存曲线 + log-rank
# -----------------------------
fit_km <- survfit(Surv(OS.time, OS) ~ TumorHigh_group, data = dat)

p_km <- ggsurvplot(
  fit_km, data = dat,
  pval = TRUE, conf.int = TRUE,
  risk.table = TRUE, risk.table.height = 0.22,
  legend.title = NULL,
  legend.labs = c("Tumor_OS.high Low", "Tumor_OS.high High")
)
print(p_km)

ggsave("KM_Tumor_OS_high_cutpoint.pdf", p_km$plot, width = 6.8, height = 5.2)

# 中位生存期
print(surv_median(fit_km))

# -----------------------------
# 4) Cox（High vs Low），输出 HR(95%CI)
# -----------------------------
fit_cox <- coxph(Surv(OS.time, OS) ~ TumorHigh_group, data = dat)
print(summary(fit_cox))
print(exp(cbind(HR = coef(fit_cox), confint(fit_cox))))

rm(list = ls())
gc()
dev.off()
