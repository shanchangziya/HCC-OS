# TCGA five-fold comparison of 117 configurations

The manuscript OSARS is the [128-gene StepCox[forward] + GBM model](../manuscript_osars128/README.md) with 1,696 boosting trees. This directory records a separate TCGA five-fold comparison of 117 configurations. Its ranking selected the 49-feature PLS-Cox predictor described below. The additional nested analysis reported in Supplementary Data 5 is separate from this comparison.

## Scope and headline result

The frozen modelling grid contains **117 configurations (20 single algorithms + 97 two-stage configurations)**, not 101. Model development and algorithm selection were repeated using stratified 5-fold cross-validation entirely within TCGA-LIHC. The selected configuration was **StepCox[backward] + plsRcox**, with a mean held-out-fold Harrell C-index of **0.7016** (SD 0.0699).

ICGC-LIRI-JP was opened only after the TCGA-selected model had been refitted on complete TCGA and saved with a SHA256 lock. Because the historical analysis had already inspected ICGC while choosing the original model, its present result is labelled **retrospective external evaluation**, not untouched or fully independent validation.

## Frozen inputs and cohort definition

- TCGA expression: `TCGA_unclipped_448_expression.csv` (343 patients; 448 frozen candidate genes).
- TCGA outcome: `TCGA-LIHC_frozen_risks.csv` (343 patients; 124 deaths/events). Only `ID`, `OS.time`, and `OS` were used; the historical risk score was not used.
- ICGC retrospective evaluation: `ICGC_original_448_expression.csv` (243 patients).
- Patient inclusion, survival-time/status definitions, and expression matrices were not redefined.
- Input file SHA256 hashes are recorded in `data/processed/TCGA_input_manifest.csv` and the external-evaluation log.
- Server pre-lock raw-input inventory: TCGA-LIHC_frozen_risks.csv	14916 bytes; TCGA_unclipped_448_expression.csv	2605255 bytes.
- The TCGA development run was locked at 2026-09-07T14:42:10Z; the ICGC file was transferred into the isolated server workspace only afterwards, at 2026-09-07T14:42:48Z.

## Five-fold comparison workflow

1. A fixed master seed (`5201314`) generated event-stratified 5-fold assignments. Held-out folds contain 67–69 patients and event fractions of 35.8%–36.2%.
2. For each outer fold, 1st/99th percentile winsorization bounds, mean-imputation values, and near-zero-variance filtering were estimated using only the 80% training partition and then applied unchanged to its held-out partition.
3. Univariable Cox screening (`P < 0.01`) was repeated independently inside each training partition. The five folds retained
   78, 152, 97, 106, 68 genes, respectively.
4. All model fitting and all algorithm-specific tuning (10-fold internal tuning where used by the frozen Mime implementation) used only that outer training partition.
5. Risk-score orientation was determined from training-partition predictions only and then fixed before scoring the held-out patients. Held-out outcomes were used only to calculate Harrell's C-index.
6. The 117 configurations were ranked by mean outer-fold C-index. ICGC was not used as a tie-breaker.
7. The pre-specified tie rule was: highest mean C-index; exact numerical ties then favour fewer input features, then a single-stage configuration, then frozen configuration order.
8. The chosen configuration was refitted once on complete TCGA, serialized, and SHA256-locked. Only then was the external-evaluation script allowed to open ICGC.

## Leakage audit

| Potential leakage | Five-fold comparison |
|---|---|
| Univariable Cox on all TCGA before CV | Prevented; repeated in each outer training fold |
| Global TCGA winsorization/imputation before CV | Prevented; parameters estimated in each outer training fold |
| Held-out fold used for tuning | Prevented; tuning calls receive training-fold data only |
| Risk direction inferred from held-out outcomes | Prevented; orientation fixed from training predictions |
| Model changed after viewing held-out performance | Prevented; one frozen grid and one pre-specified ranking rule |
| ICGC used for feature/algorithm/tuning/cutoff decisions | Prevented in revised scripts; model hash verified before ICGC is opened |

No external cohort path or object is referenced by the fold-preparation, fold-fitting, aggregation/selection, or complete-TCGA locking scripts.

## Manuscript development versus five-fold comparison

| Item | Original OSARS development | Separate five-fold comparison |
|---|---|---|
| Configuration count | Earlier draft stated 101; frozen object contains 117 | Audited and reports all 117 |
| Univariable screen | Complete TCGA before the reported multi-model analysis | Repeated inside each outer training fold |
| Expression capping/imputation | Estimated cohort-wise before model evaluation | Estimated only in each TCGA training fold; full-TCGA values locked for later use |
| Algorithm selection | ICGC results were available and used when prioritising the historical StepCox[forward] + GBM/GBM predictor | Mean 5-fold CV C-index within TCGA only |
| External cohort label | Previously treated as independent validation | Retrospective external evaluation |
| Risk direction | Historical evaluator fitted `Surv ~ RS` separately in each evaluated cohort | Fixed using training data before held-out scoring |

## Duplicate-predictor audit

The frozen historical object contains 20 exact or rank-equivalent predictor pairs. In particular, `StepCox[forward] + GBM` and `GBM` were exactly identical in the historical TCGA and ICGC score vectors. This is expected from the implementation because StepCox `forward` starts from the full Cox model and therefore may retain the complete input set.
In the revised out-of-fold predictions, 55 exact or rank-equivalent pairs were detected; they are retained as separate configurations and marked in `all_117_models_5fold_performance.csv`.

## Models within 0.01 of the best

| Rank | Configuration | Mean C-index | SD | Mean input features | Duplicate group |
|---:|---|---:|---:|---:|---|
| 1 | StepCox[backward] + plsRcox | 0.7016 | 0.0699 | 55.6 |  |
| 2 | StepCox[forward] + GBM | 0.6986 | 0.0660 | 100.2 | D15 |
| 2 | GBM | 0.6986 | 0.0660 | 100.2 | D15 |
| 4 | StepCox[both] + plsRcox | 0.6972 | 0.0756 | 54.8 |  |
| 5 | StepCox[backward] + GBM | 0.6954 | 0.0647 | 55.6 |  |
| 6 | StepCox[both] + GBM | 0.6950 | 0.0646 | 54.8 |  |
| 7 | StepCox[backward] | 0.6922 | 0.0744 | 100.2 |  |

## Predictor selected in the comparison

- Configuration: `StepCox[backward] + plsRcox`.
- Mean 5-fold C-index within TCGA: 0.7016.
- Fold values: 0.6358, 0.6808, 0.6443, 0.7495, 0.7977.
- Fold SD: 0.0699; range 0.6358–0.7977 (存在较明显的fold间波动).
- Complete-TCGA univariable genes: 128.
- Comparison model input features: 49.
- Scoring for the selected comparison predictor is implemented in `scripts/predict_locked_model.R`; preprocessing vectors, feature order, model object, parameters, and risk orientation are all stored in `data/processed/final_OSARS_model_locked.rds`.
- The explicit nine-component PLS-Cox calculation and its parameter tables are exported in `comparison_model_risk_formula.txt` and `data/processed/final_model_*`. These reproduce serialized predictions to machine precision without changing the model hash.

The archived filename `final_OSARS_model_locked.rds` and its SHA256 belong to the comparison predictor. With the archived object and cohort inputs staged, run `Rscript scripts/predict_locked_model.R expression.csv comparison_scores.csv` for its scores.

## ICGC retrospective external evaluation

The selected comparison predictor yielded a Harrell C-index of **0.6297** (95% CI 0.5269–0.7324; n=243, events=44).
No ICGC feature selection, algorithm selection, hyperparameter tuning, risk-direction selection, cutoff optimisation, or model modification was performed.

## Interpretation of performance change

The manuscript GBM had an apparent TCGA C-index of 0.8705, whereas the selected comparison predictor had a mean five-fold C-index of 0.7016. Their difference (-0.1689) combines changes in the predictor and in its evaluation.
The ICGC C-index was 0.7664 for the manuscript GBM and 0.6297 for the comparison predictor (difference -0.1368). This is a retrospective comparison because ICGC had been examined historically.
For reference only, the refitted locked model's complete-TCGA apparent C-index is 0.8534.

## Statistical limitation

The requested five outer folds provide honest held-out predictions for every individual configuration and prevent fold-level feature/tuning leakage. However, the same set of five fold estimates is used to rank 117 configurations and to quote the selected configuration's mean; therefore some winner's-curse/model-selection optimism can remain. These values are explicitly called **internal cross-validation/model-selection performance**, not independent test performance. The additional nested assessment reported in Supplementary Data 5 evaluates the configuration-selection procedure separately.

## Software fidelity and recorded deviations

- The modelling code vendors the upstream `ML.Dev.Prog.Sig` implementation used to reconstruct the frozen grid and keeps the original algorithm parameters (RSF 1,000 trees/nodesize 10; internal 10-fold tuning; GBM up to 10,000 trees; CoxBoost penalty/step tuning; and the published StepCox directions).
- Unused plotting-only package imports were removed from the vendored function; this does not alter model calculations.
- SuperPC's upstream infinite retry loop was bounded to three attempts so a failed configuration is reported rather than hanging indefinitely.
- Server-side package versions were recorded in the locked object: survival 3.5.7; randomForestSRC 3.3.1; glmnet 4.1.8; plsRcox 1.7.7; superpc 1.12; gbm 2.2.2; CoxBoost 1.5; survivalsvm 0.0.5; dplyr 1.1.3; tibble 3.2.1; miscTools 0.6.28; compareC 1.3.2; mixOmics 6.24.0.
- GBM was run with one worker to keep parallel execution deterministic; its statistical parameters are unchanged.

## Output map

- `all_117_models_5fold_performance.csv`: ranked fold-level and summary C-indices, CI, completion, feature counts, and duplicate groups.
- `cv_patient_predictions.csv`: long-format out-of-fold risk scores for every patient/configuration.
- `feature_selection_by_fold.csv`: train-only preprocessing/univariable screening counts.
- `comparison_model_summary.txt`: locked algorithm, CV performance, parameters, and risk-score method.
- `comparison_model_risk_formula.txt` and `data/processed/final_model_*`: explicit preprocessing, PLS projection, Cox coefficients, and formula validation.
- `ICGC_retrospective_evaluation.csv`: locked-model retrospective external result.
- `data/processed/duplicate_predictor_audit.csv`: revised exact/rank-equivalent pairs.
- `figures/`: 17-cm-wide combined and individual PDF/SVG/PNG panels.

Comparison results generated: 2026-09-07 22:51:49 CST
