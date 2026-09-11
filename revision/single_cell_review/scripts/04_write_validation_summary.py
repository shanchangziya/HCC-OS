"""Write the concise, evidence-bound Markdown companion for Supplementary Fig. X."""
from pathlib import Path
import pandas as pd


HERE = Path(__file__).resolve()
OUT = HERE.parents[1]
REVISION = OUT.parent
SOURCE = REVISION / "single_cell" / "data" / "processed"


gse202642 = pd.read_csv(OUT / "GSE202642_usage_summary.csv").iloc[0]
gse149614 = pd.read_csv(OUT / "GSE149614_usage_summary.csv").iloc[0]
os_summary = pd.read_csv(OUT / "GSE149614_OS_score_summary.csv").set_index("metric")
scores = pd.read_csv(OUT / "GSE149614_OS_score_metadata.csv")
sim = pd.read_csv(OUT / "GSE149614_MP_similarity_to_discovery.csv")
panel_d = pd.read_csv(OUT / "GSE149614_panel_D_pathways.csv")
enrichment = pd.read_csv(OUT / "GSE149614_NMF_program_enrichment.csv")
patient_tests = pd.read_csv(SOURCE / "patient_level_program_tests.csv")

mp4_to_mp2 = sim[(sim["discovery_program"] == "MP4") & (sim["validation_program"] == "MP2")].iloc[0]
primary_global = patient_tests[
    (patient_tests["dataset"] == "GSE149614")
    & (patient_tests["definition"] == "OS_high_top20")
    & (patient_tests["program"] == "MP4")
].iloc[0]
ros_mp2 = enrichment[
    (enrichment["validation_program"] == "MP2")
    & (enrichment["gene_set"] == "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY")
].iloc[0]

def q(x: float) -> str:
    return f"{x:.2g}" if x >= 0.001 else f"{x:.2e}"


terms = "; ".join(
    f"{row.display_term} (FDR={q(row.FDR_within_program_all_resources)})"
    for row in panel_d.sort_values("display_order").itertuples()
)

text = f"""# GSE149614 independent single-cell validation — concise record

## Scope

This folder addresses only Reviewer 2 Major Comment 1(a): independent single-cell assessment of the OS-related malignant-cell programme in GSE149614. It does not alter the manuscript main figures, bulk analyses, OSARS performance analyses, or other reviewer-response modules.

## Cohort use

| Dataset | Patient/donor number | Samples | Total cells | Cells used for OS/NMF |
|---|---:|---:|---:|---:|
| GSE202642 (frozen discovery object) | Not verifiably mapped | {int(gse202642.tissue_or_sample_n)} tissue samples | {int(gse202642.total_cells_in_frozen_analysis_object):,} epithelial cells | {int(gse202642.cells_actually_used_for_OS_and_NMF):,} historical epithelial cells |
| GSE149614 (independent validation) | {int(gse149614.patient_n)} patients | {int(gse149614.sample_n_all_cells)} samples ({int(gse149614.primary_tumor_malignant_hepatocyte_sample_n)} primary-tumour analysis samples) | {int(gse149614.total_cells_all_samples):,} | {int(gse149614.primary_tumor_malignant_hepatocyte_cell_n):,} published malignant hepatocytes |

GSE202642 used {int(gse202642.tumor_source_epithelial_cells):,} tumour-source and {int(gse202642.adjacent_source_epithelial_cells):,} adjacent-source epithelial cells. Its frozen analysis object contains {int(gse202642.CNV_high_malignant_epithelial_cells_in_object):,} CNV-high cells but was not a pure confirmed-malignant set; the 11 tissue samples are therefore not reported as 11 independent patients.

## Reused analysis and score definition

The existing GSE149614 results passed the current revision audit and were reused without refitting NMF. Author-annotated primary-tumour malignant-hepatocyte clusters were selected before scoring (13,691 cells from 10 patients). The OS composite was the frozen within-cohort mean of separately z-standardised AddModuleScore, UCell, and AUCell scores on log-normalised uncorrected RNA; 47/49 frozen HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY genes were measurable.

The composite remained continuous across the independent malignant-cell compartment (median {os_summary.loc['OS_composite', '50%']:.2f}; range {os_summary.loc['OS_composite', 'min']:.2f} to {os_summary.loc['OS_composite', 'max']:.2f}). Its within-cohort 80th-percentile label comprised {int(scores.OS_high_top20.sum()):,} cells ({100 * scores.OS_high_top20.mean():.1f}%). Panel B is deliberately continuous rather than a dichotomised survival-style display.

## Independent NMF correspondence

GSE149614 malignant cells underwent the pre-existing de novo GeneNMF workflow independently of GSE202642 (`multiNMF`, k=4–9, min.exp=0.05, seed=123; four metaprograms retained). Validation MP2 was the closest match to frozen discovery MP4: 34 shared genes, raw measured-gene Jaccard={mp4_to_mp2.raw_measured_gene_jaccard:.3f}, and eligible-background Jaccard={mp4_to_mp2.eligible_background_jaccard:.3f}. The Fisher overlap result was calculated using the 2,000 genes actually eligible for the saved validation NMF fits (BH FDR across all 16 discovery–validation comparisons={q(mp4_to_mp2.BH_FDR_all_16_comparisons)}); it is reported in source data rather than used as a visual significance badge.

## Enrichment of the key validation programme

Panel D shows pre-specified representative, non-redundant significant pathways for MP2, selected after MP2 had been fixed by the MP4 similarity comparison: {terms}. Over-representation analysis used MSigDB 2026.1.Hs Hallmark and Reactome gene sets with the same 2,000-gene NMF-eligible background, one-sided hypergeometric tests, and BH correction within each programme across both resources. The direct HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY term was not enriched for MP2 (0 overlapping genes; FDR={q(ros_mp2.FDR_within_program_all_resources)}); the displayed detoxification/oxidation terms are therefore presented as oxidative-stress-adjacent biology, not as a substituted direct ROS result.

## Figure legend

**Supplementary Fig. X | Independent single-cell assessment in GSE149614.** **A**, all-cell UMAP coloured by major cell type. **B**, continuous OS composite in published primary-tumour malignant hepatocytes. **C**, independently derived validation NMF metaprogram scores; MP2 (orange title) is the validation programme with the greatest measured-gene overlap with frozen discovery MP4. **D**, representative FDR-significant pathways for MP2; colours denote Hallmark or Reactome resource. UMAP panels were generated within GSE149614 only.

## Interpretation boundary and suggested response wording

GSE149614 was analysed independently with harmonised cell-selection and score criteria; residual cross-cohort heterogeneity cannot be completely excluded. The figure supports heterogeneous OS-composite activity and partially conserved NMF programme features, but not a claim of complete cross-patient reproduction of a fixed OS-high state. In the pre-existing patient-level primary global-top-20% comparison, frozen MP4 enrichment was not demonstrated (n={int(primary_global.n_paired)} comparable patients; P={q(primary_global.P)}, BH FDR={q(primary_global.FDR_within_definition)}); this result is retained in the audit record rather than obscured or redefined.

Suggested restrained sentence: *“An oxidative-stress-adjacent malignant-cell programme could be observed in the independent GSE149614 cohort, with partially conserved NMF metaprogram features.”*

## Output map

- `SuppFig_GSE149614_validation.pdf`, `.svg`, `.png`: 170-mm-wide complete supplementary figure (222 mm high).
- `figures/panels/*Panel_[A–D].pdf`: independently editable panel PDFs, with corresponding SVG and PNG files.
- `GSE202642_usage_summary.csv` and `GSE149614_usage_summary.csv`: cohort counts separated by patient/donor, sample, and cell scope.
- `GSE149614_OS_score_metadata.csv`: cell-level source data for Panels B/C.
- `GSE149614_MP_similarity_to_discovery.csv`: full 4×4 programme correspondence table.
- `GSE149614_NMF_program_enrichment.csv`: all tested Hallmark/Reactome ORA terms.
"""

(OUT / "GSE149614_validation_summary.md").write_text(text, encoding="utf-8")
print("VALIDATION_SUMMARY_WRITTEN")
