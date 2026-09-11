# Spatial association: frozen descriptive re-export

2026-09-05. Source models were read only. No cell2location, CellChat, MISTY or virtual-knockout model was refitted. The original reference labels and original posterior means were preserved. This folder separates copied source files (`data/raw`), export code (`scripts`), plotting-ready tables (`data/processed`), and source/QA records (`audit`). Figure rendering is assigned separately.

## Figure-ready data

| File in `data/processed` | Content |
|---|---|
| `spatial_spot_coordinates_and_abundance.csv` | 16,535 spots: `section`, `barcode`, `imagecol`, `imagerow`, `image`, `Tumor_OS-high`, three `Mac_*_TAM` columns; also OS-low and total of 28 reference states for transparency |
| `spatial_abundance_long.csv` | The same unchanged four requested abundance vectors in long form |
| `section_TAM_descriptive_spearman.csv` | All 12 section × TAM comparisons; `spearman_rho`, `n_spots`, state identities, units and descriptive-only definition |
| `section_TAM_descriptive_spearman_matrix.csv` | Four-section × three-TAM effect-size matrix |
| `section_OSlow_TREM2_descriptive_control.csv` | Four original OS-low–TREM2 descriptive correlations; optional contextual table |
| `section_metadata_and_source.csv` | Spot counts, original file paths, model gene counts/epochs, coordinate directions, observed spacing and units |
| `pooled_display_quantiles.csv` | Unmodified abundance range and pooled display percentiles per state, computed across all four sections |

Abundance is the **cell2location posterior mean cell abundance**, `means_cell_abundance_w_sf`, from the saved `cell_abundance_mean.csv`. It is an inferred cell abundance, not a measured cell count, fraction, TPM or z-score. Values are not row-normalized or transformed in these exports. Model run summaries record 300 epochs on GPU for each section; this audit does not establish convergence from epoch count alone.

The historical single-cell reference used `final_celltype` after exclusion of `Normal_Epithelial`, a minimum of 25 cells per state, and a maximum of 3,000 cells per state. Its exported metadata contains 56,753 reference cells. The 2,363 `Tumor_OS-high` reference cells all have historical `OS_group_valley=OS-high`; the 3,000 retained `Tumor_OS-low` cells all have historical `OS_group_valley=OS-low`. These are **not the new revision top-20% labels**, and reference-cell counts are not patient counts. No state was reclassified in this re-export.

## Numerical findings and interpretive scope

| Section | Spots | OS-high–TREM2 | OS-high–SLC40A1 | OS-high–HSP |
|---|---:|---:|---:|---:|
| HCC-1T | 3184 | 0.1313 | 0.2670 | 0.3361 |
| HCC-2T | 4733 | 0.7447 | 0.7202 | 0.6863 |
| HCC-3T | 4456 | 0.1664 | 0.5493 | 0.5515 |
| HCC-4T | 4162 | 0.0480 | 0.1407 | 0.2925 |

The TREM2 association is positive but heterogeneous and is strongest in HCC-2T. It exceeds both alternative TAM correlations only in that section. These results do not establish preferential TREM2 association across HCC, communication, recruitment or causal interaction. The exported files establish four section identifiers, not independently verified patient identities.

Spearman correlations are within-section summaries of all finite paired posterior-mean abundances. Original coefficients were copied and independently reproduced using average ranks; maximum absolute difference was 4.72 × 10⁻¹⁶. The original spot-level P values/FDR values are retained only in the copied source table for audit. **No P values or confidence intervals appear in the processed plotting tables.** Spatial autocorrelation, uncertain deconvolution, total cellularity and section-level heterogeneity are not addressed by treating thousands of spots as independent observations. No cross-section significance test was performed.

## Map conventions

Use x=`imagecol`, y=`imagerow`, equal aspect ratio and y increasing downward. These are coordinates exported by Seurat `GetTissueCoordinates`; their physical micrometre calibration is not established here. Show only the observed tissue spots, without interpolation or a newly invented tissue mask. `map_point_radius_suggestion` is 0.42 times the median nearest-neighbour spot spacing in these coordinate units; it is only a display suggestion, not the measured physical spot radius.

For visual comparability, use a fixed scale for each state across all four sections. If the pooled 1st–99th percentile limits are used, explicitly state that colour values are clipped for display and that correlations used full original values. A separate colour range for each state does not permit comparison of absolute abundance between states. Do not make a tissue map look like a measured NQO1 expression map: these are inferred cell-state abundances.

Recommended main figure: A, the four-section × three-TAM correlation matrix with numeric rho values and no significance stars; B–E, paired OS-high/TREM2 maps for HCC-1T through HCC-4T. All sections must be included. Additional TAM maps or the historical OS-low control may be placed in a supplement without selecting sections by strength.

## Proposed legend

**Figure 6. Section-dependent spatial associations between historical oxidative-stress states and macrophage states.** (A) Within-section Spearman correlations between posterior-mean abundances of the historical OS-high tumor reference state and TREM2, SLC40A1 or HSP macrophage states in four HCC tissue sections. All 12 comparisons are shown as descriptive effect sizes. (B–E) Spatial distributions of the historical OS-high tumor reference state and TREM2 macrophage state in HCC-1T, HCC-2T, HCC-3T and HCC-4T, respectively, comprising 3,184, 4,733, 4,456 and 4,162 analyzed spots. Abundances are frozen cell2location posterior means; the single-cell reference retains the historical valley-based OS-high labels, not the newly evaluated top-20% definition. All eight maps use the same 0–0.5 abundance colour scale, with no colour clipping; all observed values lie within this range. Spots are spatially dependent, and no spot-level P values or confidence intervals are shown. Association does not establish direct cell communication or causality.

## Reproduction and QA

Stage the restricted files in `data/raw/` as documented in the repository data manifest, then run `scripts/01_build_source_manifest.py` followed by `scripts/02_export_descriptive_spatial.py` using Python with pandas, NumPy and SciPy. The first script records source hashes; the second validates section membership, nonnegative finite abundances, exact spot-coordinate alignment, model spot counts and all 16 original descriptive correlations before exporting. `audit/export_qa.json` records PASS. Source files remain unchanged. The independent Figure 6/7 audit documents the original CellChat input and virtual-KO limitations; those old quantitative panels are not validated by this re-export.
