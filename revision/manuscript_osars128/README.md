# Manuscript OSARS: the 128-gene GBM

Figure 3 and Supplementary Data 4 report the transcriptomic OSARS derived
with StepCox[forward] + GBM. The frozen Cox-loss GBM was fitted in 343
TCGA-LIHC patients with 124 overall-survival events. It scores the expression
of 128 genes using 1,696 trees. Higher link-scale scores indicate higher
predicted risk.

| Model property | Value |
| --- | --- |
| Ordered input | [128 gene symbols](../model_validation/data/raw/frozen_128_genes.txt) |
| Training objective | Cox partial-likelihood loss |
| Interaction depth | 3 |
| Learning rate | 0.001 |
| Minimum node size | 10 |
| Bag fraction | 0.5 |
| Initial score | 0 |

The score is the sum of the terminal-node contributions across all 1,696
trees. The stored contributions already include the learning rate. The gene
influences in Supplementary Data 4 describe the fitted model. Exact scoring
uses its tree splits and terminal-node values.

## Score the frozen model

Download the [scoring bundle](https://github.com/shanchangziya/HCC-OS/raw/refs/heads/main/revision/manuscript_osars128/OSARS128_reviewer_scoring_bundle.zip)
or use the copy in this directory. It contains the frozen tree splits and leaf
values, the ordered genes, a Python scorer, and two synthetic profiles with
expected scores. The ZIP SHA-256 is
`6dfb1364d9dba971847401d8c72ab3a3d1277aa4c41a4cebfaff7d7f12ce400d`.
With Python 3, run:

```sh
python3 -m zipfile -e OSARS128_reviewer_scoring_bundle.zip osars128-review
cd osars128-review
python3 score_osars.py frozen_osars128 example_input.csv example_scores.csv
```

Compare the two synthetic scores with [the expected output](example_expected.csv)
using an absolute tolerance of 1e−12. The profile with all 128 gene values
equal to 6 gives −0.4081560724653323. The corresponding R prediction on the
original fitted object is:

```r
predict(gbm_model$gbm_fit, newdata = x, n.trees = 1696, type = "link")
```

For researchers who already have the original trusted R fitting archive,
[the exporter](export_frozen_model.R) regenerates the three model parameter
files from `StepCox_GBM_model_object.Rdata` (archive SHA-256
`a48c7de11e6e39b03bf86820f5816c98ee086e98fe272fafe494dcd32689c121`).
The scorer checks each parameter file's hash before prediction.

For your own profiles, replace `example_input.csv` with a UTF-8 CSV
containing one row per sample, a unique `sample_id`, and 128 named gene
columns. Use finite, preprocessed expression values on the analysis scale.
The scorer matches genes by name and writes one continuous OSARS per sample.

For the reported TCGA and ICGC analyses, each gene was winsorized at its
1st and 99th percentiles within the corresponding cohort. Infinite values
were treated as missing before winsorization. Missing values were then
filled with the within-cohort gene mean after winsorization. ICGC inputs
were log2(FPKM + 1). Pass these values directly to the GBM without further
gene-wise z scoring. Risk strata are assigned separately for each analysis.
The external treatment analyses used their documented missing-data procedure.

## Related analyses

- The [TCGA-only five-fold comparison](../model_selection/README.md)
  evaluated 117 configurations and selected a **49-feature StepCox[backward] +
  plsRcox** for that comparison. Its PLS coefficients and C-indices belong to
  the comparison.
- Supplementary Data 5 reports additional internal validation of the
  original StepCox[forward] + GBM **configuration** and of the 117-configuration
  selection procedure. Both refit models within the training partitions.
- The CPTAC proteomic surrogate and the image-based pathology model are
  separate analyses. This scoring program takes transcriptomic inputs.
