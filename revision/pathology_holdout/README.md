# Reproduction

Run the scripts from the project root with R 4.3, `survival`, and `glmnet` installed. Use a new empty output directory when independently reproducing the analysis; the frozen result directory intentionally refuses overwrite or split regeneration.

```bash
Rscript revision/pathology_holdout/scripts/01_audit_and_freeze_split.R SOURCE_DIR CLINICAL_CSV FRESH_OUTPUT_DIR
Rscript revision/pathology_holdout/scripts/02_fit_evaluate_strict_holdout.R FRESH_OUTPUT_DIR
python revision/pathology_holdout/scripts/03_plot_diagnostics.py FRESH_OUTPUT_DIR
```

`SOURCE_DIR` contains `Resnet.Rdata`, `pathdat.Rdata`, and `resnet50_features.csv`; `CLINICAL_CSV` is the matched TCGA clinical table. Optional fourth and fifth arguments to the first script may identify the historical feature/model scripts for checksum provenance. The first script audits the archived tile and patient objects and then creates exactly one seed-20260906 event-stratified patient split. The second script performs the unchanged training-only Pearson screen and `cv.glmnet` Cox LASSO with `lambda.min`, freezes the model and training median, and only then evaluates fixed predictions in test patients. The third script draws the two standalone diagnostics from frozen CSV files and does not refit any model.

Optional descriptive comparisons with prior pathology results use
`HCC_OS_PATHOLOGY_PRIOR_PERFORMANCE` and `HCC_OS_PATHOLOGY_ORIGINAL_PATHDAT`.
They do not alter the strict model, split, or cutoff when absent.
