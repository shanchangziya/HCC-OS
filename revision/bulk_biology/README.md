# TCGA bulk biology associated with frozen OSARS status

This self-contained analysis addresses a deliberately narrower question than direct single-cell state projection: **what bulk transcriptional phenotype is associated with frozen OSARS-high status in TCGA-LIHC?** It does not refit OSARS, alter its cutoff, or use the previously tested scRNA-derived OS-high signature.

## Design

- Cohort: the frozen TCGA-LIHC log-scale expression matrix (`n=343`).
- Exposure: the pre-existing `OSARS_frozen_group` classification (`Low n=172`, `High n=171`).
- Primary pathway analysis: limma competitive gene-set testing (CAMERA) for every available Hallmark (`n=50`) and testable Reactome (`n=1,132`) set. `Up` means higher in OSARS-high; `Down` means lower in OSARS-high.
- Sensitivity analysis: identical pathway testing after adjustment for pathologic stage (`n=321`).
- Multiple testing: Benjamini–Hochberg FDR within each gene-set library and analysis.
- Presentation-only sample scores: GSVA for Hallmark pathways. The four displayed pathways were selected mechanically as the two smallest primary CAMERA FDRs in each direction; their boxplots are descriptive and do not replace the all-pathway CAMERA inference.

The exact package/database versions, gene-set membership, sample manifest, and R session information are preserved in `references/`.

## Principal result

In TCGA-LIHC, OSARS-high tumors show a strong **proliferative** bulk program (E2F targets, G2M checkpoint, MYC targets, DNA repair and mitotic pathways) together with lower **hepatocyte-associated metabolic** programs (bile acid metabolism, xenobiotic metabolism, fatty-acid metabolism and coagulation). The principal Hallmark signals remain after pathologic-stage adjustment. This supports a bulk-level, proliferative/metabolically attenuated phenotype; it does **not** establish that a specific scRNA-seq epithelial state is preserved in bulk tumors or that the effect is tumor-cell intrinsic.

## Key outputs

- `figures/Supplementary_Figure_OSARS_high_bulk_biology_TCGA.pdf` and `.png`: composite manuscript-ready supplementary figure.
- `figures/Panel_A_Hallmark_CAMERA_TCGA.*`, `Panel_B_Reactome_CAMERA_TCGA.*`, and `Panel_C_Representative_Hallmark_GSVA_TCGA.*`: separately reusable panels.
- `tables/TCGA_OSARS_high_vs_low_CAMERA_hallmark.csv` and `..._reactome.csv`: primary inference tables.
- `tables/TCGA_OSARS_high_vs_low_stage_adjusted_CAMERA_*.csv`: stage-adjusted sensitivity tables.
- `tables/TCGA_OSARS_high_vs_low_limma_DE.csv`: full gene-level differential-expression table.
- `scripts/01_bulk_biology_enrichment.R` and `scripts/02_plot_bulk_biology_figure.py`: complete reproducible workflow.

No pre-existing manuscript figure was modified.
