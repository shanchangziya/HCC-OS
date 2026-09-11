# Code manifest

This manifest separates the current revision workflow from the historical
snapshot. File-level hashes are generated in `docs/CODE_INVENTORY.csv`.

## Current revision code

| Module | Purpose | Status |
|---|---|---|
| `revision/single_cell/` | Discovery-definition audit, GSE149614 atlas, frozen program scoring, patient-level summaries, de novo NMF, and Figures 1–2 | Core |
| `revision/single_cell_review/` | Independent GeneNMF refit, program correspondence, GSVA/ORA, and reviewer-focused panels | Core/sensitivity |
| `revision/model_selection/` | Historical 117-model audit, leakage-controlled five-fold TCGA selection, model lock, and later ICGC evaluation | Core |
| `revision/model_validation/` | Frozen-score audit, nested comparator strategy, signature benchmarking, survival/clinical analyses, assumptions, and Figures 3–4 | Core |
| `revision/osars_mp4_projection/` | Frozen MP4 projection, GSEA, and multi-method sensitivity | Supporting; negative-direction result |
| `revision/bulk_biology/` | TCGA bulk pathway phenotype associated with frozen OSARS | Supporting |
| `revision/osars_ros_correlation/` | Direct ROS–OSARS correlation figure | Supporting |
| `revision/osars_ros_overlap/` | Frozen OSARS/Hallmark ROS overlap diagram | Supporting |
| `revision/osars_discrimination/` | Compact frozen OSARS discrimination figure | Supporting |
| `revision/clinical_increment/` | TCGA stage/ROS/OSARS apparent discrimination comparisons | Sensitivity |
| `revision/clinical_cox/` | Stage- and grade-adjusted Cox models | Sensitivity |
| `revision/immune_therapy/` | Inferred immune features, TIDE predictions, treatment contexts, CIBERSORT-fit sensitivity, and Figure 5 | Core, exploratory clinical interpretation |
| `revision/spatial_association/` | cell2location input/export code, four-section descriptive correlations, and Figure 6 | Core, descriptive only |
| `revision/protein_support/` | PDC000198 pairing, NQO1 abundance, continuous survival, and surrogate audit | Supplementary support |
| `revision/nqo1_prioritization/` | Transparent post-model NQO1 prioritization and paired-protein source tables | Supplementary figure |
| `revision/pathology/` | Corrected historical TCGA pathology reassessment and supplementary figure | Core supplementary |
| `revision/pathology_holdout/` | One-time strict patient split and training-only sensitivity | Audit/sensitivity |
| `revision/scripts/build_summary_tables.py` | Assemble aggregate manuscript Tables 1–3 from audited module outputs | Reporting |
| `revision/model_validation/references/historical_all_cohort_model_building.R` | Portable preservation of the original three-cohort Mime1 call | Historical audit only |

## Historical code

| Module | Historical role | Current handling |
|---|---|---|
| `scripts/01_scRNA_oxidative_stress/` | Initial single-cell oxidative-stress, CNV, NMF, Scissor, and enrichment workflow | Retained for provenance |
| `scripts/02_cell_subtyping/` | Immune/stromal subtype analyses | Retained for provenance |
| `scripts/03_annotation_integration/` | Integrated annotation and QC | Retained for provenance |
| `scripts/04_spatial_transcriptomics/` | Earlier spatial scoring/RCTD workflow | Superseded for current Figure 6 |
| `scripts/05_bulk_deconvolution/` | Earlier BayesPrism analyses | Not used to reverse the negative MP4 projection result |
| `scripts/06_prognostic_modeling/` | Historical Mime1 development and extraction | Audited/superseded by 117-model revision code |
| `scripts/07_wsi_pathology/` | Historical WSI feature/prognosis workflow | Superseded by corrected pathology modules |
| `scripts/08_functional_analysis/` | Historical tumor, CellChat, and virtual-knockout analyses | Legacy only; not current mechanistic evidence |

## Inclusion rule

The candidate includes executable source code and small interpretation-critical
definitions. It excludes raw/controlled data, large derived matrices, serialized
objects, virtual environments, package libraries, caches, generated figures,
and logs. No placeholder CT-radiomics, docking, or unverified analysis code was
invented where a real source script was unavailable.

See `SOURCE_CODE_COVERAGE.md` for the explicit comparison with the local
manuscript-revision bundle and the reasons for the few non-analysis exclusions.
