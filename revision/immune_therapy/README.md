# Fig. 5 and treatment supplements — revision candidate, 2026-09-05

**重要限制：归档 CIBERSORT 只有 36/343 例满足 fit P < 0.05（Low 9 / High 27）。全体样本中 Tregs、M0 的组间差异，在该小型、不平衡的过滤子集中未得到显著确认。主图仅作为 inferred immune features 的返修候选，不得称为已经稳健验证的免疫抑制或细胞富集。必须一并提供 `Supplementary_CIBERSORT_fit_sensitivity.pdf` 和完整源表。**

本目录根据真实个体数据重构 Fig. 5，未修改 Article2、原始脚本、原始 Rdata 或旧图。统计与作图分开；Python 作图只读取 `data/processed/`。新增结果不包含模拟患者或编造检验值。

## Reviewable figures

| 文件（`figures/composites/`） | 尺寸 | 用途 |
|---|---|---|
| `Fig5_immune_TIDE_revised.pdf` / `.png` | 170 × 125 mm | A 免疫特征差异；B 六项 TIDE；C 预测 ICB response |
| `Supplementary_therapy_contexts.pdf` / `.png` | 170 × 68 mm | A TACE；B BIOSTORM 辅助 sorafenib RFS-benefit signature；C 黑色素瘤 TIL-ACT |
| `Supplementary_CIBERSORT_fit_sensitivity.pdf` / `.png` | 170 × 63 mm | 全体343例与fit合格36例的 Tregs / M0 敏感性 |
| `Supplementary_DrugReflector.pdf` / `.png` | 85 × 110 mm | 原始排名前20（包含未命名化合物），探索性 |
| `Supplementary_checkpoints.pdf` / `.png` | 170 × 95 mm | 可选八基因；源表达单位尚待核实，不作为当前必须提交项 |

独立面板在 `figures/panels/`，无字母角标。组合步骤才加入 Arial Bold 10 pt 角标；图内正文 Arial 7 pt，坐标标题 8 pt。Low `#0072B2`，High `#D55E00`。PDF 全部为矢量并嵌入 TrueType 字体；PNG 为 300 dpi。固定物理画布避免合图时缩小字体。旧 PNG 没有参与生成。

英文正式图注见 `LEGENDS.md`。补图编号由全文最终排版确定。

## 冻结分组、模型与统计

- TCGA-LIHC 沿用旧 Fig5 的 171 Low / 172 High，High = score ≥ cohort median；精确 cutoff 在 `data/processed/analysis_manifest.json`。将边界改为严格 `>` 的敏感性表另存，不静默改变旧图分组。仅1例边界样本改变分组，主要免疫显著性不变，预测 ICB Fisher P 从 0.001418 变为 0.002044。
- TACE / sorafenib **只读取 2026-09-04 最终 helper 表**，不使用早期 aggregate 分支。患者计数、Fisher P 和冻结汇总完全一致。其 High = score > eligible-cohort median。分别有147和67例。
- GSE100797 使用冻结 StepCox[forward]+GBM 重新直接预测：128输入、1,696树、唯一缺失基因GC用TCGA训练均值补齐，21例全部可用；High = score > median，Low11/High10。直接预测与旧个体分数最大差异 < 1e-10。R同样重新预测343例TCGA，与归档risk最大差异 < 1e-10，记录在`data/source/TCGA_reprediction_audit.csv`。
- 免疫11指标、TIDE6指标、checkpoint8基因为彼此分开的BH家族；CIBERSORT敏感性每个分析集使用22个fraction的BH家族。检验为双侧 Mann–Whitney U，渐近、连续性和ties校正。
- 连续效应为High−Low均值差。主A按全343例每个特征的SD标准化；CIB补图用原始fraction差的百分比点。区间为组内有放回抽样5,000次的非校正95%百分位bootstrap CI，seed=20260905。均值效应CI和秩检验FDR对应不同统计目标，不应要求二者数值等价。
- 所有response比例使用患者表计数、双侧Fisher和Wilson95%CI。三种治疗比较额外提供BH，图内显示原始P，图注同时给出校正值。
- 归档 CIBERSORT 不良fit、TCGA免疫计算结果未经独立临床ICB治疗队列验证，以及跨平台GBM输入均为实际限制。

## 核实后的结果和必须修改的叙述

| 内容 | 核实结果 |
|---|---|
| TIDE predicted ICB | Low84/171(49.1%)，High55/172(32.0%)，P=.001418；仅预测标签 |
| TIDE方向 | High组TIDE、Exclusion、MDSC高，**Dysfunction低**；CAF/TAM M2不显著 |
| 11个免疫指标 | FDR显著仅Tregs、M0、NK、Stromal；其余不能写成显著改变；Tregs/M0受fit限制 |
| TACE | 51/74 vs30/73，P=.0008926，三治疗BH=.002678 |
| BIOSTORM | 15/34 vs6/33，P=.03444，三治疗BH=.05166；**RFS-benefit基因签名标签关联**，不是肿瘤缩小或OSARS治疗获益预测验证 |
| Melanoma TIL-ACT | 6/11 vs2/10，P=.1827；非显著、跨癌种，不是HCC CAR-T |
| DrugReflector | 仅探索排名；模型概率不等于患者药物有效概率；输入为26,280基因，不是只用128基因 |

可用于替换 Results 的短段落：

> In TCGA-LIHC, OSARS was associated with differences in inferred immune features and TIDE-derived scores. Higher OSARS was associated with higher TIDE, exclusion, and MDSC scores but lower dysfunction scores. TIDE classified fewer OSARS-high tumors as predicted ICB responders (55/172, 32.0%) than OSARS-low tumors (84/171, 49.1%; Fisher P = 0.001418). Although the all-sample CIBERSORT analysis indicated differences in Treg and M0 fractions, these associations were not statistically confirmed after restricting to the 36 samples with fit P < 0.05. The findings therefore describe computational immune associations and do not establish observed immunotherapy benefit.

治疗段落只可称 endpoint associations；sorafenib经三比较校正后不满足FDR<.05，ACT也不显著。

## Reproduce

在仓库根目录运行。受限原始资料不入库：用 `HCC_OS_IMMUNE_INPUT_DIR`
指向 `ssgsea_result.Rdata`、`GSE100797.Rdata`、CIBERSORT/ESTIMATE/TIDE 等源表所在目录，
用 `HCC_OS_LEGACY_PROJECT_ROOT` 指向含冻结 Fig. 5 中间表的原分析根目录：

```sh
Rscript revision/immune_therapy/scripts/export_r_objects.R
python revision/immune_therapy/scripts/prepare_data.py
python revision/immune_therapy/scripts/verify_source_labels.py
python revision/immune_therapy/scripts/plot_panels.py
python revision/immune_therapy/scripts/assemble_and_qc.py
```

Python环境包含 numpy、pandas、scipy、matplotlib、pymupdf。原运行环境为 R 4.3.3 及 gbm 2.3.1；若需附加 R 库路径，可设置 `HCC_OS_R_LIBRARY`。原 R 对象仅被读取。R 导出的全部中间输入存于 `data/source/`。`prepare_data.py`进行统计；绘图不导入 scipy。

源文件位置和SHA256在`data/processed/analysis_manifest.json`，公开标签原件hash在`qc/public_source_manifest.json`，字体/矢量/尺寸检查在`qc/vector_export_qc.json`。原始资料网络记录保存在`data/source/`；GEO 67个Sor样本标签逐例完全匹配，记录见`qc/GEO_label_validation.txt`。

主要源表为`immune_patient_values.csv`、`immune_statistics.csv`、`tide_patient_values.csv`、`tide_statistics.csv`、`tide_response_statistics.csv`、`therapy_patients.csv`、`therapy_response_statistics.csv`、`cibersort_22_fraction_sensitivity.csv`。每种箱线图的五数摘要独立保存，绘图只渲染已计算摘要。
