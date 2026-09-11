# OSARS 与冻结单细胞 MP4 状态的 GSEA

本目录是对审稿意见 9 的独立、可复现补充分析；未改动任何既有原始图。

## 结论

冻结的发现队列 MP4 程序在 OSARS-high 肿瘤中没有正向保留。相反，它在两个批量队列中均显著富集于 OSARS-low 一侧；在每个基因的模型中纳入分类病理分期后，方向和显著性仍然一致。因此，这组冻结数据不能支持“OSARS-high 肿瘤保留 OS-high 恶性上皮状态转录特征”的表述。

主分析的 MP4 程序结果为：

| 分析 | TCGA-LIHC | ICGC-LIRI |
| --- | --- | --- |
| 全队列 | NES = -3.49，BH q < 0.001 | NES = -3.43，BH q < 0.001 |
| 剔除 ROS 与 OSARS 特征基因 | NES = -3.33，BH q < 0.001 | NES = -3.25，BH q < 0.001 |
| 分期校正后、剔除 ROS 与 OSARS 特征基因 | NES = -2.86，BH q < 0.001 | NES = -3.38，BH q < 0.001 |

负 NES 表示该程序富集于 OSARS-low，而非 OSARS-high。这里的“MP4”是冻结的发现队列程序标签；本分析不把它重新定义为已被独立确认的恶性细胞状态。

## 分析设计

- 高、低风险分组：冻结的队列内 OSARS 中位数分组，未重新选择阈值或训练模型。
- GSEA 排序：每个基因的 Welch t 统计量（OSARS-high 减 OSARS-low）。
- 状态程序：发现队列的 58 个 MP4 基因，固定读取自 single_cell/data/processed/discovery_NMF_genes.csv。
- 去循环敏感性：从 MP4 程序中剔除 Hallmark ROS 基因和冻结 OSARS 的 128 个模型特征中的重叠基因（50 个定义基因；TCGA 可用 49 个、ICGC 可用 48 个）。
- 分期敏感性：逐基因拟合 expression ~ OSARS group + categorical stage，并以 OSARS-high 系数的 t 统计量排序后进行 GSEA。TCGA 使用 321 个分期完整病例；ICGC 使用 243 个病例。
- 推断：每项 GSEA 采用 10,000 次与观察到的基因集大小相同的随机基因集置换；每种分析内 6 项检验使用 BH 校正。

## 图件

- figures/R1Q9_primary_GSEA_NES_summary.pdf：两个队列的主分析汇总。
- figures/R1Q9_stage_adjusted_GSEA_NES_summary.pdf：分期校正敏感性汇总。
- figures/R1Q9_TCGA_MP4minusROS_GSEA_running.pdf：TCGA-LIHC 去 ROS/OSARS 重叠后的运行富集曲线。
- figures/R1Q9_ICGC_MP4minusROS_GSEA_running.pdf：ICGC-LIRI 对应运行富集曲线。

所有图件同时提供 600 dpi PNG 预览。完整结果、基因覆盖和每个基因的排序统计量分别保存于 tables/ 和 data/processed/。

## 对反卷积的决定

不建议把反卷积作为本轮回复中用来“补强”该状态延续性的证据。直接的、冻结状态程序投射已经在两个队列中给出相反且分期稳健的方向；反卷积不能将这一结果重新解释为“高风险肿瘤保留 MP4 程序”。若后续另做反卷积，只应作为明确标为探索性的肿瘤组成分析，而不能用于维持原有的状态连接性结论。

## 扩展的样本级投射

随后完成的多方法直接投射位于 multimethod_projection/。它在同一冻结的去重叠 MP4 基因集中采用 rank module、GSVA、ssGSEA、PLAGE、z-score 与 PCA 型评分，并给出 TCGA 高低组、连续风险与分期校正的结果。它同样支持状态程序在 OSARS-high 肿瘤中较低，而不是较高。
