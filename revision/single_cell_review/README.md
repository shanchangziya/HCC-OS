# GSE149614 independent single-cell validation — concise record

## Scope

This folder addresses only Reviewer 2 Major Comment 1(a): independent single-cell assessment of the OS-related malignant-cell programme in GSE149614. It does not alter the manuscript main figures, bulk analyses, OSARS performance analyses, or other reviewer-response modules.

## Cohort use

| Dataset | Patient/donor number | Samples | Total cells | Cells used for OS/NMF |
|---|---:|---:|---:|---:|
| GSE202642 (frozen discovery object) | Not verifiably mapped | 11 tissue samples | 10,413 epithelial cells | 10,413 historical epithelial cells |
| GSE149614 (independent validation) | 10 patients | 21 samples (10 primary-tumour analysis samples) | 71,915 | 13,691 published malignant hepatocytes |

GSE202642 used 9,412 tumour-source and 1,001 adjacent-source epithelial cells. Its frozen analysis object contains 9,883 CNV-high cells but was not a pure confirmed-malignant set; the 11 tissue samples are therefore not reported as 11 independent patients.

## Reused analysis and score definition

The existing GSE149614 results passed the current revision audit and were reused without refitting NMF. Author-annotated primary-tumour malignant-hepatocyte clusters were selected before scoring (13,691 cells from 10 patients). The OS composite was the frozen within-cohort mean of separately z-standardised AddModuleScore, UCell, and AUCell scores on log-normalised uncorrected RNA; 47/49 frozen HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY genes were measurable.

The composite remained continuous across the independent malignant-cell compartment (median 0.05; range -4.30 to 2.70). Its within-cohort 80th-percentile label comprised 2,739 cells (20.0%). Panel B is deliberately continuous rather than a dichotomised survival-style display.

## Independent NMF correspondence

GSE149614 malignant cells underwent the pre-existing de novo GeneNMF workflow independently of GSE202642 (`multiNMF`, k=4–9, min.exp=0.05, seed=123; four metaprograms retained). Validation MP2 was the closest match to frozen discovery MP4: 34 shared genes, raw measured-gene Jaccard=0.268, and eligible-background Jaccard=0.298. The Fisher overlap result was calculated using the 2,000 genes actually eligible for the saved validation NMF fits (BH FDR across all 16 discovery–validation comparisons=1.85e-36); it is reported in source data rather than used as a visual significance badge.

## Enrichment of the key validation programme

Panel D shows pre-specified representative, non-redundant significant pathways for MP2, selected after MP2 had been fixed by the MP4 similarity comparison: Coagulation (FDR=4.24e-07); Drug ADME (FDR=2.16e-06); Xenobiotic metabolism (FDR=1.08e-04); Complement cascade (FDR=9.33e-04); Heme scavenging (FDR=9.78e-04); Biological oxidations (FDR=0.002). Over-representation analysis used MSigDB 2026.1.Hs Hallmark and Reactome gene sets with the same 2,000-gene NMF-eligible background, one-sided hypergeometric tests, and BH correction within each programme across both resources. The direct HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY term was not enriched for MP2 (0 overlapping genes; FDR=1); the displayed detoxification/oxidation terms are therefore presented as oxidative-stress-adjacent biology, not as a substituted direct ROS result.

## Figure legend

**Supplementary Fig. X | Independent single-cell assessment in GSE149614.** **A**, all-cell UMAP coloured by major cell type. **B**, continuous OS composite in published primary-tumour malignant hepatocytes. **C**, independently derived validation NMF metaprogram scores; MP2 (orange title) is the validation programme with the greatest measured-gene overlap with frozen discovery MP4. **D**, representative FDR-significant pathways for MP2; colours denote Hallmark or Reactome resource. UMAP panels were generated within GSE149614 only.

## Interpretation boundary and suggested response wording

GSE149614 was analysed independently with harmonised cell-selection and score criteria; residual cross-cohort heterogeneity cannot be completely excluded. The figure supports heterogeneous OS-composite activity and partially conserved NMF programme features, but not a claim of complete cross-patient reproduction of a fixed OS-high state. In the pre-existing patient-level primary global-top-20% comparison, frozen MP4 enrichment was not demonstrated (n=8 comparable patients; P=0.46, BH FDR=0.77); this result is retained in the audit record rather than obscured or redefined.

Suggested restrained sentence: *“An oxidative-stress-adjacent malignant-cell programme could be observed in the independent GSE149614 cohort, with partially conserved NMF metaprogram features.”*

## Output map

- `SuppFig_GSE149614_validation.pdf`, `.svg`, `.png`: 170-mm-wide complete supplementary figure (222 mm high).
- `figures/panels/*Panel_[A–D].pdf`: independently editable panel PDFs, with corresponding SVG and PNG files.
- `GSE202642_usage_summary.csv` and `GSE149614_usage_summary.csv`: cohort counts separated by patient/donor, sample, and cell scope.
- `GSE149614_OS_score_metadata.csv`: cell-level source data for Panels B/C.
- `GSE149614_MP_similarity_to_discovery.csv`: full 4×4 programme correspondence table.
- `GSE149614_NMF_program_enrichment.csv`: all tested Hallmark/Reactome ORA terms.
