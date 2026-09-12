# HCC-OS

Code supporting the study **“A machine learning-driven framework integrating
an oxidative stress-associated risk signature for prognostic stratification
and therapy-response prediction in hepatocellular carcinoma.”**

> **Pre-release code snapshot (12 September 2026).** This repository contains
> the author-approved code upload for the current revision. The six-main-figure
> map remains provisional until the revised manuscript and supplementary
> numbering are finalized.

## Repository map

| Location | Role | Release status |
|---|---|---|
| `revision/` | Corrected and expanded analyses used for the current six-main-figure revision candidate | Primary code for review |
| `scripts/` | Historical analysis snapshot already present in the public repository | Retained for provenance; not all scripts support the revised claims |
| `workflow/run_pipeline.py` | Single public entry point for listing, checking, and launching revision stages | Ready for dry-run/preflight use |
| `config/paths.example.env` | Names every external-data path expected by the portable scripts | Copy locally; never commit real paths |
| `docs/` | Code inventory, figure crosswalk, data manifest, dependencies, limitations, and release checklist | Review before release |

Large omics matrices, clinical tables, serialized R/Python objects, model
checkpoints, whole-slide images, and generated figures are intentionally not
tracked. Small frozen gene definitions, formula tables, and audit resources
needed to interpret the code are included.

## Quick start

```bash
python workflow/run_pipeline.py --list
python workflow/run_pipeline.py --check
python tools/audit_repository.py --check-only
```

To inspect a stage without executing it:

```bash
python workflow/run_pipeline.py --stage single-cell
```

After staging the documented inputs, add `--execute` to run that stage. Some
stages are intentionally separate—for example, TCGA-only model selection and
the later ICGC retrospective evaluation—to preserve the model-lock boundary.

## Current revision modules

- Independent GSE149614 single-cell assessment and GeneNMF sensitivity.
- Audit and TCGA-only re-evaluation of all 117 prognostic-model configurations.
- Frozen OSARS discrimination, survival, published-signature comparisons,
  clinical adjustment, proportional-hazards, and nonlinearity analyses.
- Bulk ROS/MP4 projection and pathway-level sensitivity analyses.
- Inferred immune/TIDE and treatment-context analyses.
- Four-section descriptive spatial association analysis.
- Supplementary NQO1 prioritization, PDC000198 pairing, and protein analyses.
- Corrected internal pathology evaluation and strict hold-out sensitivity.

The revision results do not support every historical claim. In particular,
the independent single-cell primary test was not significant, bulk OSARS–MP4
associations were negative, the protein score is a distinct surrogate rather
than the mRNA GBM, spatial findings are heterogeneous/descriptive, and the
pathology hold-out result is weak. See `docs/REPRODUCIBILITY_STATUS.md` before
reusing manuscript language.

## Documentation

- `docs/FIGURE_CODE_MAP.md` — current Figure 1–6 and supplementary crosswalk.
- `docs/CODE_MANIFEST.md` — module ownership and historical/current status.
- `docs/SOURCE_CODE_COVERAGE.md` — source-bundle comparison and explicit exclusions.
- `docs/CODE_INVENTORY.csv` — file-level generated inventory and SHA-256 values.
- `docs/DATA_MANIFEST.md` — inputs, access status, and environment variables.
- `docs/R_DEPENDENCIES.md` and `docs/PYTHON_DEPENDENCIES.md` — software setup.
- `docs/PRE_RELEASE_CHECKLIST.md` — decisions still requiring author approval.
- `THIRD_PARTY_NOTICES.md` — Mime1 MIT and GeneNMF GPL-3 components.

## Citation and license

Draft citation metadata are in `CITATION.cff` and require author confirmation.
The project-level license is deliberately not finalized. Third-party code keeps
its own license; in particular, one adapted GeneNMF script is GPL-3.0-only.
Choose a compatible release strategy before replacing the pre-release
`LICENSE` notice.
