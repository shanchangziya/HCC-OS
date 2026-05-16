# TCGA-LIHC反卷积分析 - 大类细胞类型版本

## 概述

这是TCGA-LIHC反卷积分析的改进版本，将单细胞数据中的小亚群合并为大类进行反卷积，同时保持肿瘤细胞的OS-high/OS-low分组。

## 与原版本的区别

### 原版本 (TCGA-LIHC_deconvolution.R)
- 使用详细的细胞亚型（~30种细胞类型）
- 包含所有细分类型：CD8_Trm, CD4_Cyto, SPP1_TAM, CXCL9_TAM等

### 新版本 (TCGA-LIHC_deconv_major_celltype.R)
- 将细胞亚型合并为大类（~10种细胞类型）
- **保留肿瘤细胞OS分组**：Tumor_OS-high 和 Tumor_OS-low
- 去除Normal_Epithelial

## 细胞类型合并策略

| 原始细胞亚型 | 合并后的大类 |
|-------------|-------------|
| Tumor_OS-high | **Tumor_OS-high** (保持) |
| Tumor_OS-low | **Tumor_OS-low** (保持) |
| Normal_Epithelial | **移除** |
| CD8_Trm, CD4_Cyto, CD4_Treg, CD4_Tfh, T:* 等所有T细胞亚型 | T |
| CD56dim_NK, CD56bright_NK, NK:* 等所有NK细胞亚型 | NK |
| SPP1_TAM, CXCL9_TAM, SLC40A1_TAM, TREM2_TAM 等所有巨噬细胞亚型 | Macrophage |
| Naive_B, Memory_B, Activated_B, IgG_Plasma, IgA_Plasma 等所有B细胞亚型 | B |
| cDC1, cDC2, pDC, moDC | DC |
| myCAF, iCAF, apCAF, MMP_CAF, VEGFA_CAF 等所有CAF亚型 | CAF |
| VEC, LSEC, LEC, Tip_Cell, Tumor_Endo 等所有内皮细胞亚型 | Endothelial |
| Neutrophil, PMN_MDSC, M_MDSC, TAN, Inflam_Neutrophil 等 | Neutrophil |

## 预期输出细胞类型

反卷积后预计包含以下大类细胞类型（约10种）：

1. **Tumor_OS-high** - 高氧化应激肿瘤细胞
2. **Tumor_OS-low** - 低氧化应激肿瘤细胞
3. **T** - T细胞（所有亚型合并）
4. **NK** - NK细胞（所有亚型合并）
5. **Macrophage** - 巨噬细胞（包括所有TAM亚型）
6. **B** - B细胞和浆细胞
7. **DC** - 树突状细胞
8. **CAF** - 癌相关成纤维细胞
9. **Endothelial** - 内皮细胞
10. **Neutrophil** - 中性粒细胞和髓系细胞

## 使用方法

### 1. 运行反卷积分析

```bash
Rscript TCGA-LIHC_deconv_major_celltype.R
```

如需后台运行，可在服务器环境中自行用 `nohup` 或任务调度系统包装上述命令。

### 2. 监控分析进度

```bash
tail -f <your_log_file>
```

### 3. 检查输出文件

```bash
ls -lh TCGA-LIHC_major_*
```

## 输出文件

| 文件名 | 说明 |
|--------|------|
| `TCGA-LIHC_major_bp.res.Rdata` | 完整的BayesPrism结果对象 |
| `TCGA-LIHC_major_cell_fractions.csv` | 细胞类型比例矩阵（样本×细胞类型） |
| `TCGA-LIHC_major_deconv_with_survival.csv` | 合并生存数据的细胞类型比例 |
| `TCGA-LIHC_major_cell_fraction_boxplot.pdf` | 细胞类型比例箱线图 |
| `TCGA-LIHC_major_cell_fraction_heatmap.pdf` | 细胞类型比例热图 |
| `TCGA-LIHC_tumor_OS_comparison.pdf` | 肿瘤OS-high vs OS-low比较图 |
| `TCGA-LIHC_major_*.pdf` | 质量控制图（相关性、离群值等） |

## 肿瘤OS分组分析

新版本特别关注肿瘤细胞的氧化应激分组：

- **Tumor_OS-high**: 高氧化应激肿瘤细胞（原80th百分位以上）
- **Tumor_OS-low**: 低氧化应激肿瘤细胞（原80th百分位以下）

输出文件 `TCGA-LIHC_tumor_OS_comparison.pdf` 包含：
1. OS-high vs OS-low散点图
2. OS-high比例在所有肿瘤中的分布直方图

可用于分析：
- TCGA样本中肿瘤细胞OS状态的异质性
- OS-high/OS-low比例与患者生存的关系
- 不同肿瘤样本的氧化应激特征

## 参数说明

与原版本相同的参数：
- 细胞抽样比例：30%
- 核心数：20
- Marker基因选择：pval < 0.01, lfc > 0.1
- 过滤基因组：Rb, chrM, chrX, chrY, Mrp, hb, MALAT1

## 分析流程时间估计

基于424个TCGA样本和约95K单细胞：

| 步骤 | 预计时间 |
|------|---------|
| 数据加载与预处理 | 5-10分钟 |
| 质量控制分析 | 10-15分钟 |
| 差异表达统计（get.exp.stat） | 30-60分钟 |
| BayesPrism反卷积 | 60-120分钟 |
| 结果保存与可视化 | 5-10分钟 |
| **总计** | **2-4小时** |

## 故障排查

### 内存不足
如果遇到内存问题，可以调整抽样比例：
- 修改第45行：`floor(0.30 * length(cells))` → `floor(0.20 * length(cells))`

### Seurat版本兼容性
脚本已包含V3/V5兼容性检查，自动处理不同版本的assay结构。

### 细胞类型合并问题
检查合并后的细胞类型统计输出，确保：
- Tumor_OS-high 和 Tumor_OS-low 都存在
- Normal_Epithelial 已被移除
- 所有其他细胞类型已正确合并

## 下游分析建议

使用反卷积结果可以进行：

1. **生存分析**
   - 肿瘤OS-high/OS-low比例与生存的关系
   - 各大类细胞比例与预后的关联

2. **相关性分析**
   - Tumor_OS-high比例与免疫细胞浸润的关系
   - T细胞、NK细胞与肿瘤比例的相关性

3. **分组比较**
   - 根据OS-high比例高低分组比较临床特征
   - 不同免疫浸润模式的生存差异

4. **Cox回归**
   - 多因素分析各细胞类型的独立预后价值
   - OS-high比例作为连续变量或分类变量

## 引用

如果使用此脚本进行分析，请引用BayesPrism方法：
> Chu, T., Wang, Z., Pe'er, D. et al. Cell type and gene expression deconvolution with BayesPrism enables Bayesian integrative analysis across bulk and single-cell RNA sequencing in oncology. Nat Cancer 3, 505–517 (2022).

## 相关脚本

- 查看原始反卷积脚本对比：`TCGA-LIHC_deconvolution.R`
