# R dependencies

The complete repository spans single-cell, survival, enrichment, spatial, and
historical analyses, so one monolithic R environment is not recommended. Use a
separate environment for the 117-model workflow and record `sessionInfo()` for
every full reproduction.

## Core revision packages

```r
cran <- c(
  "data.table", "dplyr", "gbm", "ggplot2", "glmnet", "jsonlite",
  "Matrix", "msigdbr", "patchwork", "survival", "survminer", "tidyr"
)

bioconductor <- c(
  "AUCell", "GSVA", "limma", "org.Hs.eg.db", "SingleCellExperiment",
  "UCell", "zellkonverter"
)

single_cell_special <- c("GeneNMF", "RcppML", "Seurat")
model_special <- c(
  "compareC", "CoxBoost", "miscTools", "mixOmics", "plsRcox",
  "randomForestSRC", "superpc", "survivalsvm"
)
```

The historical archive additionally references BayesPrism, CellChat,
CytoTRACE2, fastCNV, GseaVis, Mime1, Nebulosa, SCP, Scissor,
scMetabolism, scTenifoldKnk, and spacexr. Those packages are needed only when
reproducing the corresponding historical step; several historical functional
panels are not part of the revised evidence chain.

## Recorded reference versions

- Model-selection server: survival 3.5.7, randomForestSRC 3.3.1, glmnet 4.1.8,
  plsRcox 1.7.7, superpc 1.12, gbm 2.2.2, CoxBoost 1.5,
  survivalsvm 0.0.5, dplyr 1.1.3, tibble 3.2.1, miscTools 0.6.28,
  compareC 1.3.2, and mixOmics 6.24.0.
- Model-validation/immune reference run: R 4.3.3, survival 3.8-3,
  Matrix 1.6-5, and gbm 2.3.1 where applicable.
- Independent GeneNMF review: GeneNMF 0.6.2 core behavior, GSVA 1.50.0,
  and msigdbr 26.1.0 / MSigDB 2026.1.Hs.
- Corrected pathology reference run: R 4.3.1, survival 3.5-7,
  glmnet 4.1-8, jsonlite 1.8.7, and timeROC 0.4.

Exact session records are outputs of the module scripts. A complete portable
`renv.lock` has not been reconstructed and remains a pre-release limitation.

## Third-party source

The model-selection wrapper vendors one unmodified Mime1 MIT source file. The
reviewer GeneNMF script contains GPL-3.0-only adapted functions. See
`THIRD_PARTY_NOTICES.md` before choosing a repository-level license.
