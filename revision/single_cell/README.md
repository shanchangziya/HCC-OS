# 单细胞返修图及独立验证

2026-09-05；所有原稿、原模型和原始对象保留。此目录含实际重算的数据、单幅面板、组合 PDF/PNG 和可重运行脚本。

## 先看哪些文件

- `figures/composites/Fig1_discovery_definition_reconstructed.pdf`：发现对象、评分公式、历史分组与返修 top-20% 的区别、组织样本组成。
- `figures/composites/Fig2_independent_single_cell_validation.pdf`：独立图谱、冻结 MP4 评分、患者级主检验及敏感性、技术调整和独立 NMF overlap。
- `figures/composites/Supplementary_discovery_definition_audit.pdf`：严格 top-20% 的发现集 UMAP 及定义一致性，供主图择优取舍；与 Figure 1 部分重叠，不建议重复放入最终论文。
- `LEGENDS.md`：逐图英文图注和方法限制。
- `audit/independent_scientific_review.md`：独立方法审核；后续已修复事项见下。

## 结论边界

**预设主分析没有证实跨患者稳定的 OS-high/MP4 富集。** GSE149614 中，全队列 composite top-20% 的 MP4 患者级差异 P=0.46094，BH q=0.76823。去掉 ATOX1/TXN 后结果相同。供体内 UCell 分组提供敏感性支持（q=0.00488），但改变了评分器和分组问题，不能替代主分析。控制测序深度与线粒体比例后，平均 MP4−ROS partial rho=0.14859，signed-rank P=0.08398。

NMF 的冻结 discovery MP4 与 validation MP2 有 34 个重叠基因（描述性 Jaccard=0.26772）。这支持部分基因结构重现，不能替代患者级主检验或证明机制。

发现集保存对象实际上是 **10,413 个上皮细胞、11 个组织样本**，包括来自肿瘤组织的 9,412 个细胞和邻近组织的 1,001 个细胞；包含不同 CNV 分类。不能写成 11 名已核实独立患者或 10,413 个已确认恶性细胞。历史 OS-high=2,437；新 top-20%=2,083；其中历史 high 的 85.31% 来自 p5。原模型仍追溯到历史分组，没有用新标签重新训练。

## 输入、分组与软件

发现集：服务器原 `analysis/current_pipeline/05osnmf.Rdata` 中的保存对象、NMF 基因列表和评分组件；重新计算 composite 与原值最大差为 0。原 MP4 UCell 重算与保存值一致。原标签 `OS_group_valley` 完整保留。

独立验证：[GSE149614](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149614)，[Lu et al., Nature Communications 2022](https://www.nature.com/articles/s41467-022-32283-3)。使用作者 fresh counts 与 updated metadata；输入 SHA256 见 `audit/atlas_input_manifest.json`、`audit/metadata_qc.json`。全图谱 71,915 cells / 21 samples / 10 patients；按照作者恶性簇 3,4,12,15,17,19,22,24,27,42,43,47 且 site=Tumor 纳入主原发肿瘤集合，13,691 cells / 10 patients。不是根据 OS score 选择恶性身份。

可视化：NormalizeTotal(1e4)、log1p、2,500 batch-aware HVGs、30 PCs、Harmony(sample)、20-neighbor graph、UMAP(min_dist=.35)，seed=20260905。Harmony 只用于展示，不用于基因程序评分。

评分：Seurat 4.3.0.1、UCell 2.4.0、AUCell 1.22.0；未校正 RNA log-normalized 数据。composite = 三个评分各自按队列 z 标准化后的均值。验证集重新应用相同的相对定义，80th percentile=0.79170559；没有移植发现集数值 cutoff=1.01733225。验证高组=2,739 cells。冻结 MSigDB v7.0 Hallmark ROS 49 genes 中可测 47；MP4 58 中可测 57；MP4−ROS 可测 55/56。缺失基因列于 `audit/gene_coverage.json`。

统计：在每患者内求 high−low 的 UCell 均值差，每组至少 5 个细胞。主分析满足条件 8 人，两个患者因高组不足排除；供体内分位数敏感性为 10 人。效应为患者差值的等权均值，CI 为 5,000 次患者 percentile bootstrap；P 来自双侧 Wilcoxon signed-rank；在每个队列×分组定义的 5 个程序内 BH。均值 CI 和 signed-rank P 的估计目标不相同，不能由其中一个替换另一个判断显著性。相关先在患者内计算，控制 library-size ranks 和 mitochondrial-percentage ranks 后给 partial rank correlation。

QC 敏感性已完成：保留原评分和标签，限定 200≤nFeature≤6000、percent.mt<20 的 12,694 cells；主分析仍阴性（去 ROS MP4 mean difference=0.00994，P=0.46094，q=0.76823）。该阈值没有用于主图谱或主 NMF，不得将其称为严格过滤后的主结果。

NMF：GeneNMF `multiNMF(k=4:9,min.exp=.05)`，10 个原发肿瘤样本共 60 个 NMF fits；`getMetaPrograms(nMP=4, specificity=3, weight.explained=.3)`。实际内部 seed=123，已显式写入重现脚本。评分/UMAP seed=20260905 不应混写成 NMF seed。原模型未因审核而重跑。

NMF overlap 修复：以保存 W 的实际 2,000 个候选基因为条件背景，两侧 gene sets 均裁剪后，16 项比较做 BH；MP4→validation MP2 的 eligible n=44/104，overlap=34，P=1.157e−37，q=1.851e−36。原 25,712 全基因背景结果仅归档，禁止用于图中星号。主图显示未改变的 measured-gene Jaccard=.26772；eligible Jaccard=.29825 是不同分母。没有重建 discovery W 的候选全集，不能把本检验描述为对称共同候选空间检验。

## 重现顺序

`00_probe_objects.R` → `01_build_validation_atlas.py`（本地）→ `02_validate_frozen_programs.R`（服务器）→ `03_summarize_validation.py`（本地）→ `05_validation_de_novo_NMF.R`（服务器，内部调用 06）→ `04_plot_scRNA_figures.py` + `07_plot_discovery_definition.py`（本地）。06 可以单独读取保存 NMF 重算 overlap，无须重新拟合。

使用 `python` 和 `Rscript` 运行；可通过 `HCC_OS_RSCRIPT`、`HCC_OS_DISCOVERY_ROOT`、`HCC_OS_GSE149614_RDS` 等变量配置环境与受限输入（见仓库 `config/paths.example.env`）。大型 `validation_scored.rds` 和 NMF RDS 不随代码仓库分发；可公开的小型评分、基因定义、统计表与审计记录按数据清单管理。

## 审核闭环

已修复：NMF 候选背景；实际内部 seed 记录；统计协议遗漏 UCell 敏感性；新增保留标签的 technical-QC 敏感性；字体统一最低 7 pt；Figure 2 D 标签与 C 色标冲突；Figure 1 不再把保存上皮对象泛称恶性细胞。

仍然存在的研究限制：发现组织样本与 donor 的映射未核实，邻近组织高 CNV 比例需要回查，历史 high 被 p5 主导；未重新构建一个去除这些发现集问题的新模型。当前图可以用于真实报告返修分析，不能据此宣称原“稳定恶性 OS-high 状态→OSARS”主线得到完整验证。
