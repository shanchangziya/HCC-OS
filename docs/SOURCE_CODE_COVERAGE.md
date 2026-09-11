# Source-code coverage review

The manuscript-revision bundle was compared with this candidate before
handoff. Local package libraries, virtual environments, caches, generated
figures, and data files were excluded from the comparison.

The filtered revision source contained 103 project-authored R/Python/shell
files. One hundred scientific or provenance scripts are retained or explicitly
represented under the cleaned semantic directory names described below. Three
internal review-only utilities are deliberately excluded.

## Retained or deliberately adapted

- All analysis scripts in the current single-cell, model-selection,
  model-validation, bulk, clinical-sensitivity, immune, protein, pathology,
  and NQO1 modules are retained.
- The four spatial upstream scripts were moved out of `data/raw/scripts/` to
  `revision/spatial_association/upstream_scripts/` and made configurable.
- The workstation/server download helper was replaced by
  `revision/spatial_association/scripts/01_build_source_manifest.py`; the
  public workflow validates a user-staged relative input layout and contains
  no personal SSH destination.
- `revision/scripts/build_summary_tables.py` retains the manuscript-table
  assembly step.
- The original all-cohort Mime1 call is represented by
  `revision/model_validation/references/historical_all_cohort_model_building.R`
  and explicitly marked as historical, non-leakage-controlled provenance.
- The longer historical TCGA/ICGC model-building script is byte-identical to
  `scripts/06_prognostic_modeling/构建预后模型_no_GSE14520.R` and is therefore
  not duplicated under `revision/model_validation/references/`.
- The separate strict pathology hold-out scripts are included in
  `revision/pathology_holdout/`.

## Deliberately not released as analysis code

- `package_figure_review.py` only creates local review PDFs/ZIPs from generated
  figures; figures and review bundles are excluded from the code release.
- `图文说明/build_guide.py` and `图文说明/font_probe.py` create internal visual
  review aids and probe workstation fonts; they do not generate scientific
  results.
- Installed `.Rlib` package source, `.venv` content, caches, and package demos
  are dependencies rather than project-authored code and are excluded.

This is a scope audit, not proof of numerical reproduction. Clean-environment
execution remains required before public release.
