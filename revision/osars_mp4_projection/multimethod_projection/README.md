# OSARS 高低组的多方法状态投射

本补充分析把“OS 评分”明确为冻结的 OSARS_frozen 风险分数。TCGA-LIHC 的 343 个样本已由该分数的原有队列中位数切为 OSARS-low（n = 172）和 OSARS-high（n = 171）；没有重新训练模型或重新选择阈值。ICGC-LIRI（n = 243）作为外部重复队列。

## 为什么这不同于前一轮 GSEA

前一轮 GSEA 把所有基因按 OSARS-high 减 OSARS-low 的差异表达排序，并检验固定状态基因集在排序中的位置。本目录则对每一个肿瘤样本直接计算状态分数，再检验：

1. 高低 OSARS 组的状态分数是否不同；
2. 连续 OSARS 分数与状态分数是否相关；
3. 在模型 score ~ OSARS group + pathological stage 中，关联是否仍存在。

这更直接回答“高风险肿瘤是否保留 OS-high 细胞状态转录特征”。

## 冻结的状态特征

主特征为发现队列的 58 基因 MP4 程序。为避免用风险模型的输入基因反过来验证风险模型，主分析先去除 2 个 Hallmark ROS 重叠基因（ATOX1、TXN），以及 6 个 OSARS 特征重叠基因，得到 50 个固定基因。TCGA 可测 49 个，ICGC 可测 48 个。没有从 TCGA 或 ICGC 的高低组比较中再挑选基因。

这检验的是冻结的 MP4 转录程序，而不是细胞比例；它不能代替真正的细胞组分反卷积。

## 样本级投射方法

- Within-sample rank module：每个样本内的固定状态基因平均表达分位数。
- GSVA。
- ssGSEA。
- PLAGE：基于潜变量的基因集评分。
- 额外完整表格还包括 mean gene z-score、GSVA z-score 和 PCA PC1。单个基因集时 PLAGE 与 PCA PC1、mean z-score 与 GSVA z-score 分别是高度冗余的实现，因此图中只展示每对的一种代表。

PLAGE 与 PCA 的符号本身可任意翻转，已按与 rank module 的正相关方向统一，不改变任何样本间排序或统计检验。

## 主要结果

四种非冗余的样本级投射方法在 TCGA 中都显示 OSARS-high 的固定 MP4 状态分数更低，而不是更高：

| 方法 | 标准化均值差（high minus low） | BH q | 分期校正 BH q |
| --- | ---: | ---: | ---: |
| Rank module | -0.51 [-0.66, -0.35] | <0.001 | 0.004 |
| GSVA | -0.72 [-0.96, -0.50] | <0.001 | <0.001 |
| ssGSEA | -0.52 [-0.68, -0.37] | <0.001 | 0.004 |
| PLAGE | -0.87 [-1.09, -0.67] | <0.001 | <0.001 |

连续关联也均为负向；TCGA 的 Spearman rho 分别为 -0.39、-0.39、-0.40 和 -0.55（所有 BH q < 0.001）。ICGC 的方向一致且效应更强。

因此，多种对同一冻结状态基因集的聚合方式都不支持“OSARS-high 肿瘤保留 OS-high/MP4 状态”的结论。它们是算法稳健性检验，不是相互独立的生物学队列验证。

## 图件

- figures/R1Q9_TCGA_multimethod_state_scores_by_OSARS_group.pdf：TCGA 高低 OSARS 组的逐样本状态投射分数。
- figures/R1Q9_TCGA_continuous_OSARS_state_projection.pdf：TCGA 连续风险分数与状态投射分数的关系。
- figures/R1Q9_multimethod_state_projection_effects_TCGA_ICGC.pdf：TCGA 与 ICGC 的标准化组间效应量及 bootstrap 95% CI；实心点表示分期校正后 BH q < 0.05。

所有图件同时提供 600 dpi PNG。处理后的逐样本分数、统计表、基因覆盖、方法说明和可重运行脚本均保存在本目录；未改动任何原图。
