# Published HCC prognostic-signature comparison

This folder is an isolated, reproducible response-analysis output. It consumes
the project's frozen TCGA-LIHC and ICGC-LIRI survival/expression inputs and
does not modify any frozen source file or retrain any model.

## What is compared

The primary same-patient Harrell C-index comparison includes OSARS plus the
five requested articles that publish numerical Cox formulas:

- Liu et al. 2021, FAIS signature 1, the paper's selected best-performing
  signature.
- Fu and Song 2021.
- Lin et al. 2021.
- Hong and Cai 2022.
- Ma et al. 2024, using STC2. The printed equation's SCT2 is treated as an
  apparent typographical error because the selected gene and surrounding text
  state STC2.

Tian et al. 2021 identifies five genes but does not report the fitted
multivariable Cox coefficients in the main article or accessible supplement.
It is therefore deliberately unscored, rather than refitted on the present
cohorts.

Wang et al. 2022 does not report the fitted numerical score coefficients or
its fitted TPM/SVA preprocessing parameters. Its rounded multivariable hazard
ratios in Figure 3F are used only for a separate sensitivity proxy with
coefficient equal to ln(HR); that table is not an exact reproduction of the
published model.

## Inputs and consistency

The primary table uses shared complete cases for every included predictor:

| Cohort | Primary shared n (events) | Wang-proxy sensitivity n (events) |
|---|---:|---:|
| TCGA-LIHC | 343 (124) | 343 (124) |
| ICGC-LIRI | 124 (24) | 74 (16) |

The frozen filtered ICGC expression matrix lacks some required gene rows
because its original export dropped genes with any missing sample. Required
genes were recovered without imputation from the same archived ICGC FPKM
source and transformed as log2(FPKM + 1). For gene rows already present in the
frozen matrix, the recovery agrees to less than 1e-10; see the file named
ICGC_source_recovery_consistency.csv in data/processed. GSDME is recovered
from its accepted historical symbol DFNA5 and is explicitly recorded in the
mapping file.

## Outputs

- tables/head_to_head_cindex_summary.md: compact reviewer-facing results.
- tables/head_to_head_cindex.csv: point estimate, 95% bootstrap interval,
  common n, events, and pair count.
- tables/paired_cindex_difference_vs_OSARS.csv: paired same-patient
  OSARS-minus-comparator C-index contrasts.
- references/formula_audit.csv: formula, source location, preprocessing
  limitations, and reproducibility status for all seven requested articles.
- references/published_signature_coefficients.csv: coefficient-level audit,
  including the Wang sensitivity-proxy derivation.
- data/processed/OSARS_cindex_reproduction_check.csv: confirms identical
  full-cohort OSARS C-index to the frozen pipeline.
- data/processed/source_manifest.json: checksums of all consumed inputs.

## Re-run

From the repository root, set `HCC_OS_ICGC_FPKM` to the unfiltered ICGC
expression source when it is not staged under `revision/model_validation/data/raw/`, then run:

    python revision/model_validation/published_signatures/scripts/run_published_signature_head_to_head.py

The script uses 1,000 paired sample-row bootstrap replicates (seed 20260906).
It applies fixed published scores, never retrains a published signature.

## Scope limitation

TCGA-LIHC overlaps development data for several signature papers, and
ICGC-LIRI overlaps historical validation/setting data for some. This is a
descriptive head-to-head recalculation on the requested frozen cohorts, not an
independent external validation or a claim of post-publication generalization.
