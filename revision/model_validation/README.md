# Model revision statistics — 5 September 2026

This directory contains auditable statistics and plotting tables for the revised model and clinical figures. Statistical scripts do not generate plots. Figure files and `plot*.py` are owned by the visualization agent. The primary model remains the frozen legacy OSARS; the new bounded training strategy is a separate analysis.

## Findings that determine the manuscript wording

| Analysis | Result | Interpretation |
|---|---|---|
| Frozen OSARS, TCGA | C-index 0.8705, 95% CI 0.8452–0.8944 | Training apparent performance |
| New bounded TCGA strategy | Mean five outer-fold C-index 0.6789, 0.6283–0.7267 | Nested feature and model selection; CI conditional on the fitted outer models |
| Frozen OSARS, ICGC | C-index 0.7664, 0.6886–0.8343 | Retrospective evaluation in a historically consulted cohort |
| TCGA-trained new strategy, ICGC | C-index 0.7329, 0.6499–0.8063 | New model trained without ICGC outcomes; ICGC itself is not pristine |
| ICGC common benchmark subset | OSARS 0.7772; Hong 0.7437; Ma 0.7261 | All three use the same 141 samples / 26 events; harmonized formula implementations |
| ICGC stage increment | C-index 0.6921 to 0.7907; paired difference 0.0985, 0.0204–0.1847 | Coefficients trained in TCGA; retrospective discrimination, not demonstrated clinical utility |
| Verified PDC000198 proteome, OS | C-index 0.6351, 0.5604–0.7089; adjusted HR/SD 1.5134, 1.1521–1.9879 | Frozen correlation-weighted protein surrogate; 158 verified tumors |
| Verified PDC000198 proteome, RFS | C-index 0.5424, 0.4781–0.6112; adjusted HR/SD 1.0353, 0.8264–1.2969 | No clear association |
| TCGA MP4–OSARS, ROS/stage/purity adjusted | Partial rank correlation −0.4395, −0.5253 to −0.3481 | Negative association persists after these measured covariates |
| ICGC MP4–OSARS, ROS/stage adjusted | Partial rank correlation −0.5404, −0.6313 to −0.4328 | Does not support preservation of the MP4 state by OSARS |

On the common ICGC benchmark subset, paired C-index differences were 0.0335 (−0.0512 to 0.1249) versus Hong and 0.0511 (−0.0451 to 0.1658) versus Ma. Both intervals cross zero; superiority over these signatures is not demonstrated.

The apparent-minus-nested difference (0.1916, CI 0.1529–0.2355) compares different fitting strategies. It is **not** an estimate of the pure optimism of the frozen GBM. The 8-configuration nested experiment does not repeat the original 117-model search.

## Audit conclusions

1. The original result object contains 117 models, 20 single algorithms and 97 two-stage configurations, and 234 cohort/model C-index rows. `StepCox[forward] + GBM` and the single GBM use the same 128 features and identical scores. The saved final GBM has 1,696 trees, interaction depth 3, shrinkage 0.001, minimum node size 10 and bag fraction 0.5. Ten-fold CV tunes trees after the full-TCGA univariate filter; it is not nested validation of the complete pipeline.
2. The original source sorts model performance using ICGC. The selected model is then hard-coded. No prospectively fixed selection rule was located. ICGC must not be described as an untouched independent validation set.
3. The candidate origin was reconstructed exactly: 657 positive OS-high scRNA markers, intersected with genes available in TCGA, GSE14520 and ICGC, reproduce the **ordered** archived 448-gene vector. No survival filter or Scissor output was identified in this candidate step. The nested experiment is conditional on this archived, independently derived candidate set and does not re-estimate scRNA discovery.
4. Original preprocessing clips each cohort independently at its 1st/99th percentiles and uses cohort-specific means for imputation. Frozen-score auditing preserves this historical behavior. The new nested strategy estimates clipping, imputation and scaling from the current training fold only.
5. TCGA full expression and original candidate expression agree to approximately 5e-14. `log2(FPKM+1)` applied to the correctly identified ICGC matrix agrees with all 448 original candidate values to 1.8e-15; IDs and outcomes match exactly. Filenames were not taken as proof of cohort identity.

## Cohort and assumption issues

- ICGC has 243 sample IDs, but 11 of these represent five groups of identical full expression profiles. Donor IDs are unavailable. A sensitivity retaining one row per identical profile gives 237 samples / 43 events and nearly unchanged frozen C-index 0.7670. This is an expression-profile sensitivity, not confirmed patient deduplication. One identical-expression cluster has conflicting recorded stage. See `ICGC_duplicate_expression_clinical.csv`.
- The published Hong formula requires MT3, which has unobserved raw ICGC values in 99 samples; SLC7A11 has four unobserved values. These were not replaced by zero, imputed or refitted. The common benchmark subset therefore has 141 samples. The full-cohort Ma and OSARS results remain separately available.
- ICGC has only two subjects observed beyond five years; the complete benchmark subset has none. Five-year outputs are retained as flagged audit data, not primary figure evidence. Display one- and three-year AUCs. PDC follow-up also does not support five-year AUC.
- The archived `CPTAC_*` files are retained for traceability, but the cohort is the HBV-HCC proteome dataset **PDC000198**. Official PDC biospecimen mapping identifies T724 as normal tissue and P723 as tumor, in conflict with the local prefixes. Primary survival excludes T724 (158 tumors); it does not silently swap samples or invent P723 survival. The original 159-sample files are archival sensitivity results.
- Per-SD Cox HRs summarize association. TCGA score-specific Schoenfeld test P=0.0381 and log-time interaction P=0.00884 indicate a time-varying association. ICGC stage violates PH (P=0.0122); stratifying the baseline hazard by early/advanced stage retains an OSARS association, HR/SD 1.5980 (1.2586–2.0291). Adjusted spline-versus-linear tests are significant in TCGA (P=0.000960) and ICGC (P=0.00642). These diagnostics prevent treating one constant HR as a complete model of risk.
- Clinical adjustment uses age, sex and stage for TCGA/ICGC, and age/sex for PDC because stage is absent in the archived protein clinical table. Treatment, etiology, liver function, performance status and all other potential confounders are not uniformly available. ESTIMATE purity is computationally inferred and does not establish cell-intrinsic expression.

## Published comparators and provenance

- [Hong and Cai, 2022, Disease Markers](https://doi.org/10.1155/2022/6201987): eight coefficients are explicitly reported, but the exact expression unit/transformation/standardization is insufficiently specified. The primary implementation applies them to archived log-scale expression; within-cohort z-scoring is a sensitivity. Official publisher, PubMed and Crossref were checked during this revision; no retraction or update relation was returned. The Crossref response is saved in `references/Hong2022_crossref_status.json`. This records the check, not a guarantee against future status changes.
- [Ma et al., 2024, Journal of Hepatocellular Carcinoma](https://doi.org/10.2147/JHC.S465592): three coefficients are explicit. “Standardized expression” is operationalized here as within-cohort gene z-scores; unstandardized log expression is a sensitivity because the original precise standardization formula is not reported. The original model concerns TACE-related oxidative stress; ICGC general survival extends that setting.
- [Wang and Liu, 2023, Frontiers in Genetics](https://doi.org/10.3389/fgene.2022.975211): the seven genes are identifiable, but exact coefficients were not available in the inspected main text. It is documented as not evaluable and was not reconstructed by fitting new coefficients.
- Exact formulas and preprocessing limitations are in `references/signature_definitions.json`; retrieved full text and metadata are preserved. TCGA/ICGC overlap with the development or evaluation cohorts of published signatures, so this is retrospective benchmarking, not unbiased superiority testing.
- The 49-gene Hallmark ROS set is frozen from the [Broad ssGSEA2.0 MSigDB v7.0 release](https://raw.githubusercontent.com/broadinstitute/ssGSEA2.0/master/db/msigdb/h.all.v7.0.symbols.gmt). MP4 comes from the discovery NMF export. Bulk scores are mean within-sample fractional expression ranks over the same 11,222-gene universe: 45/49 ROS, 55/58 MP4 and 53/56 MP4-minus-ROS genes. ATOX1 and TXN are the two ROS–MP4 overlaps. These scores are not UCell or ssGSEA.

## Plotting handoff

| Figure content | Statistical input |
|---|---|
| New nested training strategy / old apparent contrast | `nested_Cindex_summary_CI.csv`, `nested_outer_fold_Cindex_CI.csv`, `nested_outer_fold_performance.csv` |
| Published signature C-index / one- and three-year AUC | `discrimination_metrics_with_CI.csv`; select TCGA-LIHC and ICGC-LIRI_benchmark_completecases; OSARS_frozen, Hong2022_8gene_log, Ma2024_TR_OSRG_3gene_zscore |
| ROC curves / paired comparisons | `ROC_coordinates.csv`, `paired_discrimination_differences.csv` |
| TCGA / ICGC survival curves and risk tables | `KM_coordinates.csv`, `KM_risk_tables.csv`, `KM_logrank_tests.csv`; frozen OSARS and explicit cutoff rule |
| PDC158 OS/RFS | `PDC_verified158_KM_coordinates.csv`, `PDC_verified158_KM_risk_tables.csv`, `PDC_verified158_KM_logrank_tests.csv`, `PDC_verified158_discrimination_CI.csv` |
| Adjusted HR forest | `Cox_forest_long.csv`, model multivariable_primary, term OSARS_SD; PDC158 from `Cox_sensitivity_forest_long.csv` |
| Clinical stage increments | `clinical_incremental_discrimination_CI.csv`, `clinical_incremental_paired_differences.csv` |
| ROS / MP4 linkage | `OSARS_program_correlations.csv`, `MP4_partial_rank_correlations.csv`, `MP4_partial_rank_residual_coordinates.csv` |
| Cox assumptions | `Cox_PH_diagnostics.csv`, `Cox_adjusted_linearity_diagnostics.csv`, `Cox_adjusted_spline_coordinates.csv`, `TCGA_OSARS_timevarying_HR.csv` |
| Follow-up | `cohort_survival_summary.csv`, `PDC_verified158_cohort_survival_summary.csv` |

All times are days. Divide by 365 for year axes. KM risk tables count `time >= landmark`; AUC controls require `time > horizon`. The primary median rule here is **score > median → High**; equality is Low. The old immune figure used `>=`, so its frozen groups differ by one boundary patient. Do not silently mix these conventions. Pooled OOF KM/ROC is descriptive and uses fold-training percentiles / training median thresholds; do not attach an ordinary pooled log-rank P as an unbiased model-validation test.

## Reproduce

Run from the repository root. Python dependencies are recorded in `data/processed/python_versions.json`; R information is written to `logs/clinical_R_sessionInfo.txt`. Set `HCC_OS_RSCRIPT` if `Rscript` is not on `PATH`, and configure source locations using `config/paths.example.env`.

```sh
python revision/model_validation/scripts/run_statistics.py
```

This reruns statistics from `data/raw/` plus explicitly referenced discovery MP4, protein mapping, ESTIMATE and frozen protein-score exports. The sequence is 04→05→06→07→08→09→10. Extraction scripts 01–03b document upstream Rdata/full-matrix transformations; scripts 01/02 accept source and output directories as their first two arguments. Restricted inputs are not committed and must be supplied through the documented environment variables. Output model pickles require the recorded scikit-survival/scikit-learn versions.

`source_and_code_manifest.csv` records source and code SHA-256 values. `verification_report.json` records 26 passed checks, including numerical agreement with six original timeROC AUCs, agreement of vectorized Efron Cox with R on actual data, every nested split's isolation, exactly one held-out prediction per TCGA sample, exact optimized-versus-direct bootstrap C-index agreement and sample-eligibility checks. Bootstrap intervals use 1,000 paired sample-row resamples, seed 20260905, with fitted scores held fixed. They do not cover model-training or model-selection uncertainty.

## Reviewer follow-up: stage plus grade, and competing events

**R2-M3 is addressed by the existing TCGA grade-adjusted sensitivity**, originally exported by script 06 as `multivariable_plus_grade` in `Cox_forest_long.csv`. No model was newly selected or refitted for this follow-up. The model includes frozen OSARS per cohort SD, age per ten years, sex, stage III/IV versus I/II, and grade G3/G4 versus G1/G2. It uses 309 complete cases with 104 deaths from the original 343 TCGA rows. Grade is missing in 32 rows and stage in 22, with overlap; 34 rows are excluded by the combined completeness requirement. The grade field matches the source `FIG4_clinical_utility/Fig4_TCGA_OSARS_clinical_all_matched.csv` by unique sample ID.

The adjusted OSARS HR/SD is **4.2427 (95% CI 3.4124–5.2750), P=1.13×10⁻³⁸**. The OSARS-specific Schoenfeld test is P=0.0283, although the global test is P=0.2650. This remains a training-apparent association summary with a score-specific proportional-hazards limitation; it is not independent validation or proof of a constant hazard ratio. The existing main model and figures are unchanged. `reviewer_R2M3_TCGA_grade_sensitivity.csv` includes all coefficients and PH tests; `reviewer_R2M3_TCGA_grade_completecase_inclusion.csv` records eligibility and missing fields for every TCGA sample.

**PDC000198 competing-event readiness:** the archived clinical source does contain `ostime`, `osevent`, `rsftime` and `rsfevent`; therefore, it would be inaccurate to state that death fields are absent. It lacks a verified recurrence/death/censor event-type variable, separate dated events, and an endpoint data dictionary confirming whether the recorded RFS event is recurrence only and how death without recurrence was handled. In the verified 158-row primary cohort, 17 rows have a death event but no recorded RFS event, 39 have both recorded events, 41 have an RFS event without recorded death, and 61 have neither. The 17 deaths could be candidate competing events only after confirming the endpoint definitions. T615 also has RFS time 38.93 months exceeding OS censoring time 37.37 months, leaving the joint follow-up basis unresolved. No Fine–Gray model or cumulative-incidence analysis was performed from these unverified assumptions; the current RFS results are not competing-risk analyses.

Field availability and row-level ambiguities are saved in `PDC_competing_event_field_audit.csv`, `PDC_competing_event_crossclassification.csv` and `PDC_competing_event_ambiguities.csv`. The examined `肝癌和正常组织/cliincal.txt` in the local protein tree actually contains ICGC sample IDs and fields, and was not used as PDC clinical evidence. Script `12_reviewer_grade_and_event_audit.py` reproduces this bounded audit without refitting models. `reviewer_followup_source_manifest.csv` records the specific source/result hashes used.
