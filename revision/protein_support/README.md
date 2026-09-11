# NQO1 orthogonal protein support

**NQO1 protein abundance was higher in verified tumor–adjacent tissue pairs, but continuous NQO1 protein did not associate with OS or RFS.** Its association with the protein OSARS surrogate after removing NQO1's own contribution was weak. These findings support protein-level differential expression, not an independent prognostic claim for NQO1.

## Source and pairing audit

The original RData matrices contain 10,249 proteins, 165 tumor columns, and 165 adjacent-liver columns; a tumor-only object has survival metadata for 159 samples. They derive from the **HBV-Related Hepatocellular Carcinoma — Proteome** study, PDC000198. Pairing was established from the PDC `case_id` shared by tumor and normal aliquots, not column adjacency or numeric ID arithmetic. The official biospecimen response, full query, and mapping table are preserved. [PDC study](https://pdc.cancer.gov/pdc/study/PDC000198), [original Gao et al. study](https://pubmed.ncbi.nlm.nih.gov/31585088/).

The PDC mapping contains 160 cases with both tissue types and 10 specimens whose counterpart is not assigned to the same PDC case. One additional case has a tissue-label conflict: local **T724** maps to PDC **Solid Tissue Normal**, while local **P723** maps to **Primary Tumor**. Original labels remain intact; neither specimen was relabeled. Primary paired analyses exclude this ambiguous case and the 10 specimens without a verifiable shared-case counterpart, leaving **159 verified pairs**. All 318 NQO1 values in these pairs were present in the original unshared log-ratio matrix before imputation. The same ambiguous T724 specimen was excluded from primary survival/correlation analyses, leaving **158 tumor samples**, with 56 OS and 80 RFS events. The unchanged original 159-sample survival analysis is supplied as an explicit sensitivity analysis.

The local file `cliincal.txt` was rejected as pairing evidence because it contains ICGC SA identifiers belonging to another dataset. All exclusions are recorded in `nqo1_unpaired_excluded_from_paired_test.csv` and `pdc_source_tissue_discrepancies.csv`.

The supplied preprocessing script first used `impute.knn`, then `limma::normalizeBetweenArrays`. Values are normalized unshared log2 protein ratios, not absolute concentrations. Unshared ratios use uniquely assigned peptides. [PDC quantitation documentation](https://pdc-docs.cancer.gov/pdc-docs/data-analysis-guides).

## Results

| Analysis | Estimate | 95% CI | P |
|---|---:|---:|---:|
| Paired tumor-minus-adjacent mean difference, n=159 | 0.952 | 0.775 to 1.139 | — |
| Paired median difference | 0.662 | 0.509 to 0.793 | — |
| Hodges–Lehmann paired shift | 0.815 | 0.648 to 0.998 | 1.91 × 10⁻²⁰, Wilcoxon |
| Paired standardized effect, dᶻ | 0.795 | 0.655 to 0.963 | — |
| OS HR per 1 SD NQO1, n=158, 56 events | 0.960 | 0.742 to 1.242 | 0.756 |
| OS HR, age/sex adjusted | 0.988 | 0.765 to 1.275 | 0.926 |
| RFS HR per 1 SD NQO1, n=158, 80 events | 1.016 | 0.823 to 1.255 | 0.879 |
| RFS HR, age/sex adjusted | 1.016 | 0.828 to 1.247 | 0.880 |
| Spearman NQO1 vs surrogate excluding NQO1 | 0.185 | 0.031 to 0.330 | 0.0198 |
| Pearson correlation, sensitivity | 0.144 | 0.001 to 0.292¹ | 0.0715 |

¹ The Pearson correlation CI is a percentile bootstrap interval while its P value is the usual analytical correlation test; their slight disagreement near zero reflects different inferential procedures. Spearman correlation is the rank-based primary association. The magnitude is weak and should not be overstated.

Tumor/adjacent medians were −1.261/−1.983 in the normalized source scale. There were 138 positive, 19 negative, and 2 zero within-pair differences. Nonparametric effect CIs used 2,000 paired patient bootstraps; the Wilcoxon test was two-sided with tie/continuity handling. Its Hodges–Lehmann interval is the test-derived interval. A paired t-test is supplied as a sensitivity result. Since all paired NQO1 values were originally measured, the measured-only sensitivity uses the same pairs.

The two unadjusted OS/RFS primary P values were BH-adjusted (both q=0.879). Age and sex were prespecified available covariates; stage was not present in the original tumor metadata. No biomarker dichotomization or survival-optimized cutoff was used. NQO1 was standardized within the relevant survival analysis sample. Schoenfeld diagnostics and a 3-df natural-spline versus linear likelihood-ratio check did not flag NQO1 PH violation or nonlinearity; those checks have limited power. The original survival time fields are retained without conversion; their unit is not encoded explicitly in the RData and is not required for interpreting Cox HRs under a common time rescaling.

The existing 123-protein surrogate was reconstructed exactly (maximum difference 7.11 × 10⁻¹⁴). The NQO1 weight was 0.290638874813801. We subtracted `NQO1_z × weight` from that fixed score before correlation, avoiding a direct self-inclusion correlation. Remaining association is cross-sectional and does not demonstrate mechanism or causality.

## Figure-ready files

All files are under `data/processed/`.

- `nqo1_verified_pairs.csv`: one row per verified case, tumor/adjacent abundance and paired difference. Use this for paired dots/lines. It already excludes the conflicting case.
- `paired_NQO1_statistics.csv`: effect estimates, CIs, P values and pair counts. Use the `All verified pairs` row for the primary annotation.
- `NQO1_continuous_cox.csv`: OS/RFS HR forest data; primary rows have `cohort = Tissue-label verified primary` and `term = NQO1_z`.
- `NQO1_surrogate_scatter_coordinates.csv`: NQO1 versus protein surrogate with its own contribution removed, n=158.
- `NQO1_surrogate_without_self_correlation.csv`: Spearman primary and Pearson sensitivity coefficients/intervals.
- `nqo1_all_samples_with_verified_mapping.csv`, `pdc_biospecimen_mapping.csv`, and exclusion tables: auditable tissue/case provenance. Do not plot all 330 samples as paired.

**Draft legend.** Orthogonal assessment of NQO1 protein abundance in the PDC000198 HBV-related HCC proteome. Paired comparisons included 159 patient pairs with a shared PDC case identifier and concordant source/PDC tissue labels. One case with conflicting tissue annotations and specimens lacking a verifiable counterpart were excluded. Each line connects a tumor and adjacent-liver specimen from the same verified case; values are normalized unshared log2 protein ratios. P values are from a two-sided paired Wilcoxon signed-rank test. Forest plots show NQO1 HRs per 1 SD, with Wald 95% confidence intervals, in 158 tissue-label-verified tumor samples; adjusted models include age and sex. The scatterplot relates NQO1 abundance to the existing protein OSARS surrogate after removing NQO1's own weighted contribution; the annotation is Spearman correlation. NQO1 abundance differed between tissues but did not show a prognostic association in this cohort.

Run `scripts/run_protein_support.py` using Python with NumPy/Pandas and an Rscript with survival/jsonlite. Supply the restricted proteomics source root with `--proteomics-root` (or `HCC_OS_PROTEOMICS_ROOT`) and, when not staged in `data/raw/`, the frozen surrogate tables with `HCC_OS_PDC_SURROGATE_WEIGHTS` and `HCC_OS_PDC_SURROGATE_SCORES`. No plots are generated by the statistical scripts. Source files are read-only and SHA-256 hashes are recorded without personal filesystem paths in `manifest.json`.
