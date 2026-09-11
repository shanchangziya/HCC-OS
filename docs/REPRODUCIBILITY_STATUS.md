# Reproducibility and interpretation status

## What is ready

- Current revision scripts are collected under stable module paths.
- Personal workstation/server paths have been replaced with relative paths,
  command-line arguments, or documented environment variables.
- Small frozen gene definitions, the 117-configuration grid, explicit locked
  PLS-Cox parameters, and formula audit tables are included.
- Raw data, large objects, generated outputs, environments, and caches are
  excluded by default.
- Third-party Mime1 and GeneNMF source components are identified with license
  copies and exact upstream commits.
- A file-level inventory/checksum generator and a stage-level entry point are included.
- Static checks passed on 11 September 2026: 85 R, 58 Python, and 2 shell
  files; 145 code files total. The hygiene audit found no known personal path,
  obvious secret, banned cache/environment, or public candidate file over 5 MB.

## What still requires a clean-room run

The full scientific workflow has not been rerun from public downloads in a
fresh environment because several inputs are large, derived, controlled, or
stored outside this repository. Syntax checks and static repository checks do
not prove numerical end-to-end reproduction. Before release, at least one clean
environment should execute every retained primary stage using hashes from the
approved data manifest.

Known environment gaps are an exact portable R lock across all modules and an
exact lock for the historical GPU cell2location run. Serialized scikit-learn/
scikit-survival objects also require the frozen model-validation environment.

## Scientific reporting boundaries

This package preserves results that challenge the original narrative:

- The independent GSE149614 global-composite primary comparison did not show
  the prespecified MP4 enrichment; within-patient UCell is a sensitivity.
- OSARS and MP4 were negatively associated in TCGA and ICGC, including limited
  adjusted analyses. A claim that OSARS-high tumors preserve MP4-high biology
  is not supported.
- The historical TCGA C-index is apparent. The frozen object contains 117, not
  101, model configurations, and ICGC had been consulted historically.
- The revised TCGA-only five-fold selection chose StepCox[backward] + plsRcox;
  its post-selection performance remains internal cross-validation, not an
  untouched external test.
- The PDC score is a correlation-weighted 123-protein surrogate and must not be
  described as direct deployment of the mRNA GBM.
- TIDE outputs are predictions; the treatment-context analyses do not establish
  clinical treatment benefit.
- Spatial associations vary strongly across four sections and are descriptive.
  The legacy CellChat and virtual-knockout panels do not establish communication
  or causality.
- NQO1 paired protein abundance supports a tissue difference, not independent
  prognosis or causal macrophage regulation.
- The corrected pathology hold-out result is weak and is an internal
  reassessment of a previously examined TCGA cohort, not external validation.

These boundaries should travel with any code release, revised manuscript,
figure legend, and response letter.

## Reproducibility levels

| Level | Meaning | Current state |
|---|---|---|
| Static | Code parses/compiles; repository contains no known personal paths or obvious secrets | Passed for this local review candidate; repeat after any edit |
| Structural | Entry points, code/figure map, inputs, and dependencies are documented | Ready for author review |
| Formula/data-definition | Interpretation-critical gene sets, formulas, model grid, and selected audit tables are present | Ready, subject to data-rights review |
| Numerical | Every retained stage reproduces frozen outputs from staged inputs | Previously run in project environments; not yet clean-room verified from this package |
| Independent | A third party can reproduce from public/approved sources alone | Not yet established |
