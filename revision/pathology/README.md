# Corrected exploratory pathology analysis — 2026-09-05

The corrected analysis **does not establish significant test-set survival separation or incremental value over the clinical model**. It should replace the outcome-optimized pathology result and be presented as exploratory supplementary evidence.

## Frozen result

| Quantity | Training | Historical internal test |
|---|---:|---:|
| Patients / OS events | 230 / 86 | 100 / 34 |
| Low / High | 115 / 115 | 42 / 58 |
| Common training cutoff | 1.001726611104115 | 1.001726611104115 |
| High-versus-Low HR (95% CI) | 4.339 (2.695–6.984) | 1.687 (0.825–3.448) |
| Log-rank P | 8.03 × 10⁻¹¹ | 0.1472 |
| C-index (95% bootstrap CI) | 0.744 (0.686–0.807) | 0.590 (0.478–0.704) |
| 1-year AUC (95% CI) | 0.817 (0.733–0.893) | 0.598 (0.413–0.792) |
| 3-year AUC (95% CI) | 0.791 (0.711–0.869) | 0.689 (0.508–0.822) |
| 5-year AUC (95% CI) | 0.775 (0.680–0.866) | 0.611 (0.420–0.782) |

The original 230/100 patient membership was recovered from `pathdat.Rdata:dat_all$Set`; `set.seed(2026)` reproduced the same membership. All 3,780 tiles and 364 patient-level means were recovered. The saved full 2,048-feature matrix matched fresh per-patient tile means to a maximum absolute difference of 3.55 × 10⁻¹⁵. Thirty-four patients lacked membership in the original matched cohort and were recorded, leaving 330.

Training-only Pearson correlation screening retained 113 features at |r| > 0.2 and nominal P < 0.05. BH-adjusted P values are supplied for all 2,048 tests. The final LASSO-Cox model retained 13 nonzero features. The selected lambda was 0.06711232781216821, corresponding to 0.3486365227678085 of the training lambda maximum. The alternative 1-SE lambda is recorded for transparency and was **not** fitted or evaluated as a competing test-selected model.

Ten training folds were stratified by event status with seed 20260905. Correlation screening and feature standardization were repeated within each fold's training portion. A fixed 60-point lambda-fraction grid from 1 to 0.001 was compared using held-out Breslow partial-likelihood deviance per event. The minimum mean validation deviance chose the fraction. Full training data then determined the final screen, coefficients, and median cutoff. No test-set cutpoint search was performed.

The retained features are resnet110, resnet193, resnet258, resnet261, resnet349, resnet350, resnet450, resnet546, resnet1293, resnet1357, resnet1548, resnet1901, and resnet2002. Full-precision raw-scale coefficients, training means, and SDs are in `data/processed/lasso_coefficients.csv`. The risk score is the sum of raw feature values multiplied by these coefficients. A score strictly greater than the training median is High; ties are Low.

## Clinical comparison

The pre-existing TCGA clinical table provided stage, grade, and sex for 301 patients: 207 training patients with 73 events and 94 test patients with 29 events. All matched survival outcomes and molecular OSARS values agreed exactly with the pathology source. The baseline formula was `Surv(time_days,event) ~ stage_advanced + grade_high + sex_male`, where advanced stage is III–IV versus I–II, high grade is G3–G4 versus G1–G2, and male is versus female. The combined model added the frozen continuous pathology score. Both models were fitted in complete-case training patients and applied unchanged to complete-case test patients. Age was not available in that source table.

In the same 94 test patients, clinical and combined C-indices were 0.642 (95% CI 0.541–0.757) and 0.646 (0.529–0.759). The paired difference was **0.0036 (95% CI −0.1312 to 0.1258)**. Test AUC differences at 1, 3, and 5 years were 0.0105, −0.0049, and −0.0383; every difference CI included zero.

## Statistical definitions and QA

- Kaplan–Meier coordinates use `survival::survfit` with log-log 95% confidence intervals. Times are in days; figure years use 365.25 days. Censor coordinates and annual risk counts are exported separately. Group HRs use univariate `coxph`, with Low as reference; P values use two-sided log-rank tests.
- Harrell C-index uses `survival::concordance(..., reverse=TRUE)` because larger frozen scores represent higher risk. [R survival documentation](https://stat.ethz.ch/R-manual/R-devel/library/survival/html/concordance.html).
- Cumulative/dynamic AUC uses marginal Kaplan–Meier censoring weights: cases are deaths by the horizon, controls survive past it; cases receive 1/G(T−), ties receive half credit. The independent implementation agreed with `timeROC` to within 2.77 × 10⁻⁵. [timeROC documentation](https://search.r-project.org/CRAN/refmans/timeROC/html/timeROC.html).
- C-index/AUC CIs are percentile intervals from 500 patient bootstraps with seed 20260906, with the fitted models and predictions fixed. The same resampled patients were used across clinical and combined models. All requested replicates were estimable. These CIs describe evaluation-sample uncertainty conditional on the fitted model, not full model-development uncertainty.
- The final glmnet path converged (`jerr=0`), including the selected lambda. An initial lambda-max probe produced a numerical warning only at its smallest lambda; no model was selected there. This warning remains in the original analysis log.
- Schoenfeld tests did not flag pathology-score PH violations in training (P=0.475) or test (P=0.365). Small test-event counts limit diagnostic power. Five-year test AUC has only 14 controls remaining and a wide interval.
- `independent_output_qa.json` records that coefficients reconstruct all scores within 6.00 × 10⁻¹⁵; all groups and risk tables agree with patient records; curves are monotone; train/test IDs are disjoint; and source SHA-256 hashes were unchanged.

## Figure-ready files and draft legend

All files below are under `data/processed/`.

| File | Purpose |
|---|---|
| `km_coordinates.csv` | Post-step KM survival/CI coordinates by split and group |
| `km_censor_coordinates.csv` | Censor marks at the corresponding KM height |
| `km_risk_table.csv` | Annual patient counts at risk, 0–5 years |
| `survival_group_statistics.csv` | N, events, common cutoff, log-rank P, HR and CI |
| `performance_metrics.csv` | C-index/AUC and bootstrap CI; keep `All` and `Clinical_complete` distinct |
| `paired_incremental_performance.csv` | Paired CI for combined-minus-clinical performance |
| `cv_curve.csv`, `coefficient_path.csv` | Training CV and coefficient paths |
| `frozen_patient_predictions.csv` | Patient IDs, split, outcomes, scores, groups, and clinical covariates |
| `source_manifest.json`, `input_qa.json` | Provenance, inclusion counts, hashes |

**Draft legend.** Corrected exploratory evaluation of a ResNet50-derived pathology score in TCGA-LIHC. **(A,B)** Overall-survival Kaplan–Meier curves in 230 training patients and 100 historical internal test patients, respectively. Screening, standardization, cross-validation, and risk-threshold determination were confined to training patients. The same training median threshold (1.0017266) classified both sets. Curves show log-log 95% confidence intervals, censor marks, and numbers at risk; HRs compare High with Low, and P values are two-sided log-rank tests. **(C)** Training and test Harrell C-indices for the frozen pathology score. **(D)** C-indices of the pathology, clinical, and combined models in the same 94 complete-case test patients. Clinical models were fitted only in training patients and included stage, grade, and sex. The displayed incremental value is combined-minus-clinical C-index. Error bars are 95% percentile intervals from 500 patient bootstraps of frozen predictions, paired across models. Test survival separation and incremental performance were not statistically established. This is a corrected internal reassessment of a previously examined TCGA cohort. Molecular OSARS itself was developed in TCGA and served as the pathology screening target; residual target non-independence remains. It is not an external pathology validation.

The visualizer generated `figures/Supplementary_pathology_corrected_internal.pdf` and `.png`, with individual A–D panels under `figures/panels/`. The statistical backend did not generate plots.

Rendering: Arial 7 pt; Low `#0072B2`; High `#D55E00`. Figures are produced by the designated visualizer; these statistical scripts do not draw plots.

## Reproduction

Run `scripts/run_pathology_analysis.py` with Python containing NumPy/Pandas and an Rscript containing survival, glmnet, jsonlite, and optionally timeROC. Defaults read `data/raw/`; provide the existing clinical CSV through `--clinical-source`. The reference run used R 4.3.1, survival 3.5-7, glmnet 4.1-8, jsonlite 1.8.7, and timeROC 0.4. `scripts/audit_processed.py` checks the exported results without fitting any model. A frozen model RDS and all intermediate training feature tests are retained.
