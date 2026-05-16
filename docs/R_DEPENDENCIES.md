# R Dependencies

Packages detected in the included scripts:

```r
cran_packages <- c(
  "aplot", "circlize", "clustree", "dplyr", "forcats", "gbm",
  "ggplot2", "ggpubr", "ggraph", "ggrepel", "glmnet", "gridExtra",
  "igraph", "Matrix", "msigdbr", "patchwork", "pheatmap",
  "RColorBrewer", "scales", "stringr", "survival", "survminer",
  "tibble", "tidyr", "tidyverse", "viridis"
)

bioc_packages <- c(
  "AUCell", "ComplexHeatmap", "DESeq2", "edgeR", "enrichplot",
  "GSVA", "limma", "org.Hs.eg.db", "ReactomePA"
)

special_packages <- c(
  "BayesPrism", "CellChat", "CytoTRACE2", "fastCNV", "GeneNMF",
  "GseaVis", "harmony", "Mime1", "Nebulosa", "SCP", "Scissor",
  "scMetabolism", "scTenifoldKnk", "Seurat", "spacexr", "UCell"
)
```

Install source varies by package version. For reproducible manuscript work, record the exact package versions used in the analysis environment with:

```r
sessionInfo()
```

Some packages may require GitHub installation or package-specific setup:

- `BayesPrism`
- `CellChat`
- `CytoTRACE2`
- `fastCNV`
- `GeneNMF`
- `GseaVis`
- `Mime1`
- `Scissor`
- `SCP`
- `spacexr`

The scripts were kept close to the original analysis code; dependency installation is not performed automatically inside most scripts.
