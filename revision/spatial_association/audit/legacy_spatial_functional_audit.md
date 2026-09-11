# Legacy spatial and functional source audit

2026-09-05. Scope: Article2 figure captions/methods, local `../FIG/FIG6`, local `../虚拟敲除/NQO1ST`, submitted functional-analysis scripts, and read-only inspection of original server objects and saved spatial outputs. No new CellChat/KO/MISTY inference or enrichment was performed. Original panel letters are caption-based mappings; exact final assembly should be checked against the source figure before reuse.

## Decisions for original panels

| Original caption panel | Traceable source | Recommendation |
|---|---|---|
| Fig6 A macrophage UMAP | `../FIG/FIG6/maccelltype.pdf`; source `08_Mac_annotated.Rdata` through `mac细胞通讯分析.R` | May retain as a descriptive embedding of annotated cells, with original subtype labels and cell counts. It does not supply independent-patient evidence. |
| Fig6 B marker expression | `../FIG/FIG6/MACdot.pdf`; SLC40A1, TREM2 and heat-shock marker sets | May retain as annotation support. Marker selection and display are descriptive; do not add cell-level significance as patient evidence. |
| Fig6 C interaction scatter and correlation inset | `通讯点图.pdf` / outgoing-incoming strength; `NQO1_TREM2_correlation_KWftw.pdf` | Hold interaction strengths pending normalized-input CellChat reanalysis. **Do not retain the inset as an NQO1 correlation:** its actual x-axis is `log2(10 Signatures TPM)`, despite the filename. |
| Fig6 D–E communication networks | `通讯.pdf`, `interaction.pdf`, `weight.pdf`, `细胞通讯点图.pdf`; `maccellchat.Rdata` | Hold quantitative weights, ligand-receptor P values and TREM2-preference claims. Saved object confirms raw-count input, pooled cells and no donor metadata. |
| Fig6 F cell2location spatial deconvolution | `Fig1A_HCC2T_STdeconvolve_style_spatial_pies.pdf`; `11_make_hcc2t_revised_fig1ab.py` | May retain only as an illustrative HCC-2T posterior-abundance composition map. Filename means STdeconvolve **visual style**; algorithm is cell2location. Pies sum 28 states into 13 display groups and normalize within each spot. They are not measured cell proportions. |
| Fig6 G spatial co-localization | HCC-2T maps and/or `Fig1C_HCC2T_selected_MISTY_{intra,juxta.120,para.360}_heatmap.pdf` | Prefer replacement with all four sections. If MISTY remains in supplement, label unsigned model importance and view/radius units; do not call it positive interaction, directionality or a replicated effect. |
| Fig6 H correlation/enrichment summary | `cell2location_spot_correlations.csv`, `cell2location_tumor_high_quartile_enrichment.csv` | Reuse descriptive rho for **all four** sections and all three TAM controls. Remove spot-level P/FDR and naive spot-group Wilcoxon claims. Original strongest-section selection cannot support a general preference. |
| Fig6 I abundance maps | `Fig1B_HCC2T_OShigh_and_three_macrophage_spatial_maps.pdf` | Data may be retained, but main figure should include four sections with consistent within-state scales. Original script independently clips each map at its own 1st/99th percentile; apparent contrasts are not comparable across maps unless disclosed. |
| Legacy NQO1 maps | `../虚拟敲除/NQO1ST/HCC{2,3,4}/HE-*.png`, `annotation-celltype-p*.png`, `1728.png` | Potential descriptive tissue illustrations only after source identification. Current gene-map filenames are numeric; available exports alone do not establish gene identity, expression unit/normalization, scale or exact section correspondence. Do not relabel as verified NQO1 maps from filename alone. |
| Legacy KO-perturbed network | `res_knk.Rdata`, `虚拟敲除.R`, `虚拟敲除可视化.R` | Existing selected-gene network should not support a mechanistic main claim. Only NQO1 and CHL1 pass stored adjusted P<0.05. A fully disclosed algorithmic sensitivity result could be supplemental; target self-perturbation is expected by construction. |
| Legacy literature/GO enrichment | Hard-coded 30-gene list, GO at P/q cutoffs 0.2, selected named pathways | Do not retain confirmatory enrichment claims. Nominal prefiltering, manual HLA removal, expanded hand-picked genes and non-explicit tested universe break the claimed unbiased downstream chain. |
| Legacy pathway density / STRING networks | UCell/Nebulosa on original unchanged `tcos`; local STRING interaction TSV/SVG files | Do not label as post-KO pathway activation or observed transcriptional change. These are baseline pathway-score maps plus prior protein-interaction annotations. At most relocate as explicitly exploratory annotations after selection is fully disclosed. |

## Spatial source and limitations

Original analysis root was external to this repository. Re-stage the restricted inputs under `revision/spatial_association/data/raw/` using the relative layout in the data manifest.

- Input reference: `19_sce.all.integrated.Rdata`; exported reference metadata identifies 56,753 cells after state-level filtering/downsampling. All 2,363 `Tumor_OS-high` reference cells have historical `OS_group_valley=OS-high`; all 3,000 retained OS-low reference cells have historical low labels. This does not implement the revised top-20% definition.
- Spatial object: `visium_ST_processed.RData` (`st_merged`); section identities HCC-1T to HCC-4T. Independent patient identities are not established by the inspected exported files.
- Saved outputs: `results/cell2location_spatial/<section>/cell_abundance_mean.csv`; posterior means from `means_cell_abundance_w_sf`. These are neither the q05 posterior abundance nor spot proportions. Four section run summaries record GPU and 300 epochs; convergence requires diagnostics beyond this summary.
- Coordinates: `inputs/spatial/<section>/spatial_coords.tsv.gz` from Seurat `GetTissueCoordinates`, x=imagecol, y=imagerow, increasing downward. Physical micrometre calibration is not verified.
- Original `04_colocalization_analysis.R` runs correlations and quartile-group Wilcoxon tests on spots. Neither test accounts for spatial autocorrelation. Tiny spot-level P values must not be treated as independent-patient replication. Posterior uncertainty and common tissue cellularity are also unpropagated.
- Original `07_make_positive_sample_figures.py` fixes HCC-1T, HCC-2T and HCC-3T, excluding HCC-4T. The HCC-2T-only main examples represent the strongest TREM2 correlation, not a prespecified typical section.

| Section | Spots | OS-high/TREM2 rho | OS-high/SLC40A1 rho | OS-high/HSP rho | OS-low/TREM2 rho |
|---|---:|---:|---:|---:|---:|
| HCC-1T | 3184 | 0.13131 | 0.26701 | 0.33606 | −0.12439 |
| HCC-2T | 4733 | 0.74470 | 0.72015 | 0.68635 | 0.60002 |
| HCC-3T | 4456 | 0.16643 | 0.54931 | 0.55146 | 0.04915 |
| HCC-4T | 4162 | 0.04795 | 0.14072 | 0.29247 | −0.26782 |

TREM2 is weaker than both comparison TAM correlations in three sections. Thus the evidence supports heterogeneous, section-dependent abundance associations. It does **not** support a general preferential TREM2 association, recruitment, cross-cell signaling, or an NQO1→TREM2 causal route. Deconvolution from related reference signatures also requires caution about common latent factors; no reference-label uncertainty is represented by these descriptive coefficients.

MISTY uses standardized abundance views, intra/juxta/para representations and radii 120/360 in exported coordinates. Forest importance is unsigned and prediction-based; high importance alone gives no positive/negative relationship or molecular causal direction. The inspected workflow does not introduce patient-level validation or spatially blocked validation. Retain MISTY only as exploratory model structure, if needed; it adds little to a conservative main figure once all sections are visible.

## CellChat: confirmed implementation issue

`mac细胞通讯分析.R` merges three TAM subtypes with `Tumor_OS-high`, samples 30% of pooled cells with seed 123 and supplies `GetAssayData(...,slot="counts")` to `createCellChat`. Metadata contains only `celltype`. It restricts the database to secreted signaling, uses `computeCommunProb(raw.use=TRUE)` and `filterCommunication(min.cells=10)`.

The saved `maccellchat.Rdata` confirms 27,778 genes × 4,989 cells: HSP_TAM 276, SLC40A1_TAM 3,392, TREM2_TAM 609, Tumor_OS-high 712. All stored nonzero data values are integers (range 1–9,997), so this is an actual raw-count run, not merely an unused code line. Object options are single RNA mode, triMean, nboot=100, seed.use=1 and population.size=FALSE. Equal or unequal cell counts should not be interpreted as population-adjusted effects.

The official [createCellChat documentation](https://github.com/sqjin/CellChat/blob/master/man/createCellChat.Rd) specifies a normalized expression matrix rather than counts; the [current vignette](https://github.com/jinworks/CellChat/blob/main/tutorial/CellChat-vignette.Rmd) uses normalized expression. `raw.use=TRUE` selects unprojected expression, not permission to pass unnormalized counts. Hence a normalized-input rerun is required before presenting these weights as valid quantitative communication evidence. Even a corrected pooled run remains inferred signaling potential rather than donor-replicated or experimental communication.

The apparent NQO1 correlation inset is a separate provenance problem. The three local PDFs named for NQO1 have x-axis `log2(10 Signatures TPM)`: TREM2 r=0.49, HSPA1A r=0.33, SLC40A1 r=0.025. Neither the filename nor the old caption can establish what “10 Signatures” means. Remove NQO1 attribution until exact input and source are recovered.

## Virtual knockout: saved object and selection chain

Source `tcos.Rdata` contains **10,413 epithelial cells**, with 9,412 from `group=HCCtumor` and 1,001 from `group=Normal`. All are `celltype=epi`, and the object has 9,883 CNV-high, 42 CNV-low, 134 CNV-mid and 354 reference-class cells. Samples p5 (4,832 cells) and p7 (4,087) dominate. These cells are not a verified pure malignant-cell set from 11 independent patients. Immune marker responses may reflect composition, labeling, contamination or model structure; this audit does not identify which explanation applies.

The original KO script uses normalized data only for selecting 5,000 variable genes, then supplies raw counts for those genes (plus target if absent) to scTenifoldKnk. Counts are appropriate for this tool. However, reported “default parameters” are inaccurate: seed 1, qc=FALSE, nc_nNet=20, nc_nCells=min(500,nCells), nc_nComp=10, td_K=3 and ma_nDim=3 are explicitly set. The completed log records a 5,000 × 10,413 input. The saved object is valid even though a later `nrow(res_knk)` check errors because the object is a list.

The stored `diffRegulation` has 5,000 rows with gene, distance, Z, FC, nominal P and adjusted P. Nineteen rows have nominal P<0.05; only **NQO1 and CHL1** have adjusted P<0.05. C1QA, C1QB, AIF1, C1QC, TYROBP, CD74 and other highlighted immune genes have stored adjusted P=1. The zero stored P values for the two top hits are numerical underflow and should not be typeset as exact P=0.

The visualization script first filters nominal P<0.05 and removes genes beginning `HLA`; it then introduces a hard-coded list of 30 genes including others outside the nominal threshold. GO uses P/q cutoffs 0.2 and no explicit tested-gene universe, then manually selects macrophage activation, fatty-acid metabolism and immunoglobulin-mediated response terms. This does not justify confirmatory pathway enrichment or a claim that NQO1 coordinates those functions.

The [scTenifoldKnk primary publication](https://pmc.ncbi.nlm.nih.gov/articles/PMC9058914/) describes manifold-based differential-regulation distance after an in-silico network perturbation. This is not experimental gene knockout, and the algorithmic FC/distance should not be read as a signed expression fold change. The later pathway scores are calculated from the unchanged original `tcos`, so they show baseline score localization, not a predicted or observed post-knockout expression state. STRING maps are prior interaction annotations rather than experimentally established links induced by this intervention.

## Safest revised Figure 6

Use a coherent five-panel figure: A, all four sections × all three TAM states, showing the 12 descriptive Spearman coefficients without stars; B–E, side-by-side historical OS-high and TREM2 abundance maps for each of the four sections. Show section spot counts, keep within-state colour limits consistent across sections, and state any display clipping. Preserve the historic reference label; do not relabel the original cell2location output as revised top-20% validation.

Macrophage UMAP/marker annotation can remain as supplementary support. Hold the old count-input CellChat networks, mismatched NQO1 inset and selected virtual-KO pathway panels from the mechanistic main evidence chain. New independent paired NQO1 protein expression may support orthogonal expression evidence in its own clearly identified panel/figure, but its negative survival findings and weak self-removed surrogate correlation should remain visible. Protein abundance plus spatial association still does not establish NQO1-driven TAM recruitment.

Following explicit authorization, `../data/processed` now contains the full frozen spatial coordinate/abundance export and descriptive correlation tables for the five-panel revision. This does not repair or endorse the original CellChat/KO analyses. No original source, selected gene list, spatial model or output from another agent was modified.
