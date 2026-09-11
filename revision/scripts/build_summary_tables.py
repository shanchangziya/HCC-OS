#!/usr/bin/env python3
"""Assemble manuscript-support tables from audited outputs; do no fitting."""

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "tables"
OUTPUT.mkdir(exist_ok=True)

genes = (ROOT / "model_validation/data/raw/frozen_128_genes.txt").read_text().split()
ros = set((ROOT / "model_validation/references/HALLMARK_ROS_msigdb_v7.0_genes.txt").read_text().split())
mp = pd.read_csv(ROOT / "single_cell/data/processed/discovery_NMF_genes.csv")
mp4 = set(mp.loc[mp.program == "MP4", "gene"])
assert len(genes) == 128 and len(set(genes)) == 128

pd.DataFrame(
    {
        "input_order": range(1, 129),
        "gene": genes,
        "in_Hallmark_ROS_v7_0": [gene in ros for gene in genes],
        "in_frozen_discovery_MP4": [gene in mp4 for gene in genes],
    }
).to_csv(OUTPUT / "Table2_frozen_OSARS_128_genes.csv", index=False)

cohort_rows = [
    ("GSE202642", "scRNA-seq", "10,413 cells; 11 tissue samples", "Historical epithelial discovery collection", "Includes 9,412 tumor-source and 1,001 adjacent-source cells; donor mapping unresolved; mixed CNV labels", "single_cell/data/processed/discovery_historical_and_top20_counts.csv"),
    ("GSE149614", "scRNA-seq", "71,915 atlas cells; 21 samples; 10 patients", "Independent frozen-program assessment", "13,691 primary-tumor malignant-cluster cells; primary paired comparison 8 patients; QC sensitivity 12,694 cells", "single_cell/data/processed/validation_primary_summary.csv"),
    ("TCGA-LIHC", "Bulk RNA-seq", "343 samples; 124 OS events", "Legacy model development; new nested strategy", "Legacy performance apparent; age/sex/stage complete cases 321; plus-grade sensitivity 309", "model_validation/data/processed/cohort_survival_summary.csv"),
    ("ICGC-LIRI", "Bulk RNA-seq", "243 sample IDs; 44 OS events", "Retrospective evaluation; historically consulted during ranking", "Common published-signature subset 141; distinct-expression sensitivity 237, not confirmed donor deduplication", "model_validation/data/processed/cohort_survival_summary.csv"),
    ("PDC000198", "Proteomics", "159 verified pairs; survival 158 tumors", "Orthogonal protein expression and frozen protein surrogate", "Paired and survival subsets differ; 56 OS / 80 RFS events; T724 excluded from survival; 123-protein surrogate is not mRNA GBM", "protein_support/README.md"),
    ("HCCDB spatial HCC-1T to HCC-4T", "Spatial transcriptomics", "4 sections; 16,535 spots", "Descriptive historical-reference abundance associations", "Section IDs are not verified independent patient IDs; no spot-level significance used", "spatial_association/data/raw/inputs/spatial_manifest.csv"),
    ("GSE104580", "Expression; GEO platform GPL570", "147 samples", "TACE response-context association", "Low74/High73; retrospective endpoint association", "immune_therapy/data/source/GSE104580_GEO.txt"),
    ("GSE109211", "Expression; GEO platform GPL13938", "67 sorafenib-treated samples", "BIOSTORM adjuvant RFS-benefit signature label", "Not tumor-shrinkage response; source series contains additional placebo cases not included in this analysis", "immune_therapy/data/source/GSE109211_GEO.txt"),
    ("GSE100797", "Expression; GEO platform GPL11154", "21 samples", "Exploratory melanoma TIL adoptive-cell therapy", "Cross-cancer; not HCC CAR-T; null comparison", "immune_therapy/data/source/GSE100797_GEO.txt"),
    ("TCGA-LIHC WSI", "ResNet50 pathology features", "330 matched patients; 230 training / 100 test", "Corrected historical internal reassessment", "Overlaps TCGA molecular cohort; feature selection and cutoff within training; not external WSI validation", "pathology/data/processed/patient_manifest.csv"),
]
pd.DataFrame(
    cohort_rows,
    columns=["dataset", "modality_or_source_platform", "analysis_size", "role", "limitation", "source"],
).to_csv(OUTPUT / "Table1_analysis_cohorts.csv", index=False)

processed = ROOT / "model_validation/data/processed"
discrimination = pd.read_csv(processed / "discrimination_metrics_with_CI.csv")
discrimination = discrimination[
    discrimination.cohort.isin(["TCGA-LIHC", "ICGC-LIRI", "ICGC-LIRI_benchmark_completecases"])
].copy()
discrimination = discrimination[
    (discrimination.metric == "C-index") | discrimination.time_days.isin([365, 1095])
]
pdc = pd.read_csv(processed / "PDC_verified158_discrimination_CI.csv")
pdc = pdc[(pdc.metric == "C-index") | pdc.time_days.isin([365, 1095])]
performance = pd.concat([discrimination, pdc], ignore_index=True)
performance["interpretation"] = performance.cohort.map(
    lambda cohort: "training apparent / retrospective benchmarking"
    if cohort == "TCGA-LIHC"
    else "retrospective historically consulted cohort"
    if cohort.startswith("ICGC")
    else "verified protein-surrogate cohort; not mRNA GBM"
)
performance.to_csv(OUTPUT / "Table3_retrospective_discrimination.csv", index=False)

for filename in (
    "nested_Cindex_summary_CI.csv",
    "original_117_model_performance.csv",
    "reviewer_R2M3_TCGA_grade_sensitivity.csv",
):
    source = processed / filename
    if source.exists():
        pd.read_csv(source).to_csv(OUTPUT / filename, index=False)

print(
    "TABLES_COMPLETE",
    len(performance),
    "performance rows; ROS overlap",
    [gene for gene in genes if gene in ros],
)
