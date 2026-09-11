# Locked exploratory pathology reanalysis

This reanalysis repairs the previously documented all-patient feature selection and outcome-informed test cutoff. No original model, feature file, or submitted figure is overwritten.

1. Inspect the original tile-level ResNet50 feature matrix and saved model objects. Aggregate the 2048 image features to TCGA patient means exactly as in the original pipeline. Validate unique patient identifiers and survival units.
2. Reuse the saved 230 training / 100 test patient identities if recoverable; otherwise document why that split cannot be reproduced before locking any replacement. Never choose a split using its performance.
3. Correlate image features with the existing molecular OSARS score using training patients only. Preserve the original absolute Pearson r > 0.2 and nominal P < 0.05 eligibility rule for comparability; report Benjamini–Hochberg FDR for all 2048 tests and disclose that nominal selection is exploratory.
4. Fit LASSO Cox regression in training patients. Standardization and cross-validation use training data only. Fix CV folds and seed; choose the prespecified minimum CV deviance lambda. All predictors and the training median risk cutoff are frozen before accessing test outcomes for evaluation.
5. Apply the training-derived coefficients, transformations, and median cutoff to test patients. Use High for risk > training median and Low otherwise. Do not optimize a test cutoff.
6. Report per-split N/events, selected features, coefficients, lambda, cutoff, Kaplan–Meier estimates with confidence limits/censor coordinates/risk tables, log-rank P, group HR with 95% CI, and continuous-score C-index with 95% CI. Estimate 1/3/5-year IPCW AUC with uncertainty when follow-up supports estimation; flag horizons with inadequate support. Clinical covariate adjustment is conditional on source availability and predeclared covariate coding, never test-based variable selection.

Training performance is apparent performance. The historical test cohort has already been used in the original analysis, so this is a corrected internal reassessment, not a newly untouched holdout. Moreover, molecular OSARS was developed using TCGA and used to screen pathology features; its target construction still creates a non-independent biological target even after pathology screening is restricted to training patients. This cannot support a claim of external validation.

Outputs are numerical plotting-ready tables only. Rendering specification: Arial 7 pt; Low #0072B2; High #D55E00; all panels identified as exploratory internal TCGA reassessment.
