# Figure-to-code crosswalk

This crosswalk follows the current six-main-figure revision. The former Figure
7 NQO1 analysis is now supplementary. Final supplementary numbering still
requires author review.

| Figure | Current content | Primary code | Important boundary |
|---|---|---|---|
| Figure 1 | Discovery object, oxidative-stress definition, historical versus top-20% labels | `revision/single_cell/scripts/00_probe_objects.R`, `02_validate_frozen_programs.R`, `03_summarize_validation.py`, `07_plot_discovery_definition.py` | The saved object contains epithelial cells from 11 tissue samples, not 11 verified independent patients or a purely confirmed malignant set |
| Figure 2 | Independent GSE149614 atlas, patient-level tests, technical sensitivity, and NMF correspondence | `revision/single_cell/scripts/01_build_validation_atlas.py` through `06_recompute_NMF_overlap.R`; reviewer refinements in `revision/single_cell_review/scripts/` | The prespecified global-composite primary test was not significant; within-patient UCell is a sensitivity |
| Figure 3 | OSARS/program direction, 117-model audit and TCGA-only selection, published-signature comparison | `revision/model_validation/scripts/04_nested_tcga_cv.py`, `05_prepare_scores_and_metrics.py`, `07_nested_and_incremental_summaries.py`; `revision/model_selection/scripts/`; published-comparator submodule; historical three-cohort provenance in `revision/model_validation/references/historical_all_cohort_model_building.R` | TCGA historical C-index is apparent; ICGC was historically consulted; OSARS–MP4 is negative in both bulk cohorts |
| Figure 4 | TCGA/ICGC/PDC survival, adjusted HRs, assumptions, and clinical discrimination increments | `revision/model_validation/scripts/05_prepare_scores_and_metrics.py` through `12_reviewer_grade_and_event_audit.py`; stage/grade sensitivity modules | PDC uses a 123-protein correlation-weighted surrogate, not the mRNA GBM; HRs require PH/nonlinearity caveats |
| Figure 5 | Inferred immune features, TIDE scores, and predicted ICB response | `revision/immune_therapy/scripts/export_r_objects.R`, `prepare_data.py`, `verify_source_labels.py`, `plot_panels.py`, `assemble_and_qc.py` | TIDE is predictive only; CIBERSORT fit sensitivity is weak/unbalanced; treatment cohorts are exploratory contexts |
| Figure 6 | Four-section spatial abundance associations | `revision/spatial_association/upstream_scripts/`, `scripts/01_build_source_manifest.py`, `02_export_descriptive_spatial.py`, and plotting scripts | Section-level heterogeneity is descriptive; four sections are not automatically four verified patients; no communication or causal claim |

## Proposed supplementary crosswalk

| Provisional item | Code source |
|---|---|
| NQO1 prioritization and orthogonal PDC000198 protein evidence (former Figure 7) | `revision/nqo1_prioritization/` and `revision/protein_support/` |
| Published-signature 1- and 3-year AUC | `revision/model_validation/` and its `published_signatures/` submodule |
| Cutoff transfer, protein RFS, and Cox sensitivities | `revision/model_validation/scripts/06_clinical_survival.R`, `08_sensitivity_statistics.py`, `09_cox_assumption_sensitivities.R` |
| Cox nonlinearity/time-varying effects | `revision/model_validation/scripts/09_cox_assumption_sensitivities.R` |
| Treatment contexts | `revision/immune_therapy/` |
| CIBERSORT fit sensitivity | `revision/immune_therapy/scripts/prepare_data.py` and plotting code |
| DrugReflector exploratory ranking | `revision/immune_therapy/`; optional, not efficacy evidence |
| Corrected internal pathology assessment | `revision/pathology/` |
| Strict hold-out pathology audit | `revision/pathology_holdout/`; additional sensitivity |
| Bulk OSARS-high biology | `revision/bulk_biology/`; optional supporting figure |
| MP4 projection/GSEA | `revision/osars_mp4_projection/`; required if the historical MP4-preservation claim is discussed |

## Excluded historical evidence

The current main-figure map does not use the legacy raw-count CellChat networks,
the mislabeled NQO1 correlation inset, or the selected virtual-knockout pathway
panels in `scripts/08_functional_analysis/`. Their provenance and validity
limitations are documented under `revision/spatial_association/audit/`.
