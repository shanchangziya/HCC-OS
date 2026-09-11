# Data manifest

The repository is code-first. Raw data, controlled clinical data, large derived
matrices, and serialized model objects are not distributed. This document tells
reviewers exactly what must be re-staged and which small definitions are tracked.

## External and restricted inputs

| Analysis | Expected source | Configuration | Repository policy |
|---|---|---|---|
| Discovery scRNA-seq | Archived discovery root containing `05osnmf.Rdata` and related frozen objects (study accession GSE202642) | `HCC_OS_DISCOVERY_ROOT` | Not bundled; regenerate from accession/source workflow where permitted |
| Independent scRNA-seq | GSE149614 fresh count matrix and updated metadata; intermediate Seurat RDS | `HCC_OS_GSE149614_COUNTS`, `HCC_OS_GSE149614_METADATA_CSV`, `HCC_OS_GSE149614_METADATA`, `HCC_OS_GSE149614_RDS` | Public-source downloads and large derived objects are not bundled |
| Historical OSARS objects | `QWEN0208_no_GSE14520.Rdata`, `res_no_GSE14520.Rdata`, `StepCox_GBM_model_object.Rdata` | `HCC_OS_DISCOVERY_ROOT`, `HCC_OS_FROZEN_RESULT` | Not bundled; required for audit/extraction |
| TCGA-only 117-model selection | Frozen TCGA expression and survival tables in an isolated input directory | `HCC_OS_MODEL_RAW_DIR` | Not bundled; ICGC files must be absent during TCGA model selection |
| ICGC-LIRI | Archived expression/clinical objects and unfiltered FPKM matrix | `HCC_OS_EXTERNAL_DATA_ROOT`, `HCC_OS_ICGC_FPKM` | Not bundled; use according to ICGC access terms |
| Model/clinical revision | Original project-level frozen score, clinical, ESTIMATE, and protein-surrogate tables | `HCC_OS_LEGACY_PROJECT_ROOT`, `HCC_OS_TCGA_CLINICAL_GRADE` | Not bundled; filenames are documented by module README files |
| Immune/TIDE/treatment | Frozen Fig. 5 input directory and project-level treatment tables | `HCC_OS_IMMUNE_INPUT_DIR`, `HCC_OS_LEGACY_PROJECT_ROOT` | Large matrices and derived patient tables are not bundled |
| Spatial transcriptomics | Four section coordinate/abundance exports and cell2location run summary | Stage under `revision/spatial_association/data/raw/`; optional `HCC_OS_SPATIAL_WORKDIR` | Not bundled; relative file layout is validated by `01_build_source_manifest.py` |
| PDC000198 proteomics | Tumor/adjacent RData, unshared-ratio table, preprocessing script, and frozen surrogate weights/scores | `HCC_OS_PROTEOMICS_ROOT`, `HCC_OS_PDC_SURROGATE_WEIGHTS`, `HCC_OS_PDC_SURROGATE_SCORES` | Not bundled; PDC biospecimen mapping is queried by accession |
| TCGA pathology | `Resnet.Rdata`, `pathdat.Rdata`, optional tile-level `resnet50_features.csv`, and matched clinical CSV | `HCC_OS_PATHOLOGY_SOURCE_DIR`, `HCC_OS_PATHOLOGY_CLINICAL` | WSI/features and patient tables are not bundled |

Copy `config/paths.example.env` to `config/paths.env`, replace placeholders,
and source it locally. `config/paths.env` is ignored by git.

## Tracked small resources

The following are intentionally included because they define the analysis or
permit formula-level review without exposing patient-level data:

- Frozen 128-gene, 448-candidate-gene, pseudobulk, Hallmark ROS, and NMF gene lists.
- NQO1 discovery marker and frozen variable-importance audit tables.
- The complete 117-configuration grid and explicit locked PLS-Cox formula parameters.
- Published-signature coefficient/formula audit tables.
- Small reviewer-facing GeneNMF program/pathway summaries.
- Frozen pathology model specification, PDC GraphQL query text, and public-source hashes.

These exceptions are enumerated explicitly in `.gitignore`; all other
`revision/**/data/**` files are ignored by default.

## Generated outputs

Each module writes only beneath its own directory. Common output locations are
`data/processed/`, `tables/`, `figures/`, `logs/`, and `results/`. Generated
figures, logs, serialized models, caches, and large tables remain untracked.
Before deciding to release any newly generated table, confirm that it contains
no sample-level restricted fields and that its source terms permit redistribution.

## Integrity

Module scripts record SHA-256 or MD5 checksums of consumed files using source
identifiers or repository-relative paths. Personal filesystem roots are not
written into new public manifests. Run `python tools/audit_repository.py` to
refresh the code inventory and `python tools/audit_repository.py --check-only`
for a non-writing hygiene check.
