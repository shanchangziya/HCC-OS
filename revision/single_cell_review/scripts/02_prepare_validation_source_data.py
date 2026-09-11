"""Prepare only source-data tables for the GSE149614 validation supplement.

Inputs are the pre-existing, audited single-cell outputs.  This script does
not compute scores, select malignant cells, fit NMF, or calculate statistics.
It makes the tables consumed by the plotting script and the requested cohort
summaries traceable in one self-contained Reviewer 2 output folder.
"""
from pathlib import Path
import numpy as np
import pandas as pd


HERE = Path(__file__).resolve()
OUT = HERE.parents[1]
REVISION = OUT.parent
SOURCE = REVISION / "single_cell" / "data" / "processed"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


all_cells = pd.read_csv(SOURCE / "GSE149614_all_cells_embedding.csv")
malignant = pd.read_csv(SOURCE / "GSE149614_primary_malignant_embedding.csv")
scores = pd.read_csv(SOURCE / "validation_primary_scores.csv")
nmf_scores = pd.read_csv(SOURCE / "validation_de_novo_NMF_scores.csv")
nmf_metrics = pd.read_csv(SOURCE / "validation_de_novo_NMF_metrics.csv")
similarity = pd.read_csv(SOURCE / "cross_cohort_NMF_overlap.csv")
discovery = pd.read_csv(SOURCE / "discovery_validation_ready.csv")
enrichment = pd.read_csv(OUT / "GSE149614_NMF_program_enrichment.csv")
selection = pd.read_csv(OUT / "GSE149614_key_MP_selection.csv")

require(len(all_cells) == 71915, "Unexpected GSE149614 all-cell count.")
require(len(malignant) == 13691, "Unexpected GSE149614 primary malignant-cell count.")
require(malignant["patient"].nunique() == 10, "Unexpected GSE149614 patient count.")
require(malignant["sample"].nunique() == 10, "Unexpected GSE149614 analysis-sample count.")
require(set(malignant["Cell"]) == set(scores["cell"]), "OS-score cells do not match malignant embedding.")
require(set(malignant["Cell"]) == set(nmf_scores["Cell"]), "NMF-score cells do not match malignant embedding.")
require((malignant["site"] == "Tumor").all(), "Analysis cells must be primary-tumor cells.")
require((malignant["published_malignancy"] == "Malignant hepatocyte").all(),
        "Analysis cells must retain the published malignant-hepatocyte annotation.")
require(len(discovery) == 10413, "Unexpected frozen discovery analysis-object cell count.")
require(selection.loc[0, "selected_validation_program"] == "MP2",
        "The audited MP4-to-validation mapping should select validation MP2.")

# Table 1: explicitly distinguish the historical analysis set from confirmed
# malignancy.  The frozen GSE202642 object is epithelial-only and contains 530
# non-CNV-high cells; it must not be relabelled as a pure malignant cohort.
discovery_tumor = int(discovery["group"].astype(str).str.contains("tumor", case=False, na=False).sum())
discovery_adjacent = int(len(discovery) - discovery_tumor)
cnv_high = int((discovery["cnv_class"] == "CNV-high(malignant)").sum())
gse202642_usage = pd.DataFrame([{
    "dataset": "GSE202642",
    "analysis_role": "frozen discovery analysis object",
    "patient_or_donor_n": np.nan,
    "patient_or_donor_reporting": "Not verifiably available in the frozen object; 11 tissue samples are not treated as 11 independent patients.",
    "tissue_or_sample_n": int(discovery["sample"].nunique()),
    "total_cells_in_frozen_analysis_object": int(len(discovery)),
    "epithelial_cells_in_frozen_analysis_object": int((discovery["celltype"] == "epi").sum()),
    "tumor_source_epithelial_cells": discovery_tumor,
    "adjacent_source_epithelial_cells": discovery_adjacent,
    "CNV_high_malignant_epithelial_cells_in_object": cnv_high,
    "cells_actually_used_for_OS_and_NMF": int(len(discovery)),
    "usage_definition": "Historical epithelial-cell set used by the frozen OS/NMF analysis; it includes CNV-high, CNV-mid, CNV-low, and reference cells.",
}])
gse202642_usage.to_csv(OUT / "GSE202642_usage_summary.csv", index=False)

# Table 2: all published GSE149614 cells versus the independently defined,
# primary-tumor malignant-hepatocyte analysis subset.
gse149614_usage = pd.DataFrame([{
    "dataset": "GSE149614",
    "analysis_role": "independent validation cohort",
    "patient_n": int(all_cells["patient"].nunique()),
    "sample_n_all_cells": int(all_cells["sample"].nunique()),
    "total_cells_all_samples": int(len(all_cells)),
    "primary_tumor_malignant_hepatocyte_sample_n": int(malignant["sample"].nunique()),
    "primary_tumor_malignant_hepatocyte_cell_n": int(len(malignant)),
    "cells_actually_used_for_OS_and_de_novo_NMF": int(len(malignant)),
    "selection_definition": "Published malignant-hepatocyte clusters 3, 4, 12, 15, 17, 19, 22, 24, 27, 42, 43, and 47 restricted to site=Tumor; selection was independent of OS score.",
}])
gse149614_usage.to_csv(OUT / "GSE149614_usage_summary.csv", index=False)

# All-cell UMAP source data for panel A.
all_plot = all_cells[["Cell", "sample", "patient", "site", "celltype", "published_malignancy", "UMAP1", "UMAP2"]].copy()
all_plot = all_plot.rename(columns={"Cell": "cell_id", "celltype": "major_cell_type"})
all_plot.to_csv(OUT / "GSE149614_all_cells_UMAP_metadata.csv", index=False)

# All score components (and de novo NMF scores used in panel C) are retained at
# the cell level.  The values are read directly from the pre-existing scoring
# tables; standardized scores are not recalculated here.
score_cols = [
    "cell", "patient", "sample", "site", "stage", "virus", "res.3",
    "published_malignancy", "technical_qc_pass", "UMAP1", "UMAP2",
]
score_meta = malignant.set_index("Cell").join(
    scores.set_index("cell")[["OS_AddModuleScore", "OS_UCell", "OS_AUCell", "OS_composite", "OS_high_top20"]],
    how="inner", validate="one_to_one",
).join(
    nmf_scores.set_index("Cell")[["MP1_deNovo", "MP2_deNovo", "MP3_deNovo", "MP4_deNovo"]],
    how="inner", validate="one_to_one",
).reset_index().rename(columns={"Cell": "cell_id", "res.3": "published_cluster"})
score_meta = score_meta.rename(columns={"cell": "source_cell_id"})
score_meta = score_meta[[c for c in [
    "cell_id", "patient", "sample", "site", "stage", "virus", "published_cluster",
    "published_malignancy", "technical_qc_pass", "UMAP1", "UMAP2",
    "OS_AddModuleScore", "OS_UCell", "OS_AUCell", "OS_composite", "OS_high_top20",
    "MP1_deNovo", "MP2_deNovo", "MP3_deNovo", "MP4_deNovo",
] if c in score_meta.columns]]
require(len(score_meta) == 13691 and score_meta["OS_composite"].notna().all(),
        "Invalid source table for OS-score panel.")
score_meta.to_csv(OUT / "GSE149614_OS_score_metadata.csv", index=False)

# The complete 4 × 4 cross-cohort NMF table is retained.  Raw Jaccard is a
# descriptive measured-gene similarity; the Fisher P/FDR use the 2,000-gene
# validation-NMF eligible background and are left as audit data rather than
# decorations on the figure.
sim = similarity.copy().rename(columns={
    "discovery": "discovery_program",
    "validation": "validation_program",
    "jaccard": "raw_measured_gene_jaccard",
    "jaccard_eligible": "eligible_background_jaccard",
    "P": "Fisher_P_eligible_background",
    "FDR": "BH_FDR_all_16_comparisons",
})
sim["similarity_rank_within_discovery"] = (
    sim.sort_values(["discovery_program", "raw_measured_gene_jaccard", "BH_FDR_all_16_comparisons"],
                    ascending=[True, False, True])
       .groupby("discovery_program").cumcount().add(1)
       .reindex(sim.index)
)
sim["closest_validation_program_for_discovery_program"] = sim["similarity_rank_within_discovery"].eq(1)
sim["selected_for_panel_D"] = (
    (sim["discovery_program"] == "MP4") & (sim["validation_program"] == "MP2")
)
sim.to_csv(OUT / "GSE149614_MP_similarity_to_discovery.csv", index=False)

# Representative (non-redundant) significant terms for MP2 are selected after
# the program itself was fixed by the MP4 similarity table.  This panel set is
# intentionally not a filtered claim of direct ROS-pathway enrichment.
panel_d_specs = pd.DataFrame([
    ("HALLMARK_COAGULATION", "Coagulation", 1),
    ("REACTOME_DRUG_ADME", "Drug ADME", 2),
    ("HALLMARK_XENOBIOTIC_METABOLISM", "Xenobiotic metabolism", 3),
    ("REACTOME_COMPLEMENT_CASCADE", "Complement cascade", 4),
    ("REACTOME_SCAVENGING_OF_HEME_FROM_PLASMA", "Heme scavenging", 5),
    ("REACTOME_BIOLOGICAL_OXIDATIONS", "Biological oxidations", 6),
], columns=["gene_set", "display_term", "display_order"])
panel_d = panel_d_specs.merge(
    enrichment[enrichment["validation_program"] == "MP2"], on="gene_set", how="left", validate="one_to_one"
)
require(len(panel_d) == 6 and panel_d["FDR_within_program_all_resources"].lt(0.05).all(),
        "Each selected representative MP2 pathway must be FDR significant.")
panel_d["minus_log10_FDR"] = -np.log10(panel_d["FDR_within_program_all_resources"])
panel_d = panel_d.sort_values("display_order")
panel_d.to_csv(OUT / "GSE149614_panel_D_pathways.csv", index=False)

# Simple, compact score summary used in the prose report and table audit.
score_summary = score_meta[["OS_AddModuleScore", "OS_UCell", "OS_AUCell", "OS_composite"]].describe().T
score_summary.index.name = "metric"
score_summary.reset_index().to_csv(OUT / "GSE149614_OS_score_summary.csv", index=False)
nmf_metrics.to_csv(OUT / "GSE149614_validation_NMF_metrics.csv", index=False)

manifest = pd.DataFrame([
    {"output_file": "GSE202642_usage_summary.csv", "source": "single_cell/data/processed/discovery_validation_ready.csv", "purpose": "Requested discovery-cohort usage summary"},
    {"output_file": "GSE149614_usage_summary.csv", "source": "single_cell/data/processed/GSE149614_all_cells_embedding.csv; GSE149614_primary_malignant_embedding.csv", "purpose": "Requested validation-cohort usage summary"},
    {"output_file": "GSE149614_all_cells_UMAP_metadata.csv", "source": "single_cell/data/processed/GSE149614_all_cells_embedding.csv", "purpose": "Panel A source data"},
    {"output_file": "GSE149614_OS_score_metadata.csv", "source": "single_cell/data/processed/GSE149614_primary_malignant_embedding.csv; validation_primary_scores.csv; validation_de_novo_NMF_scores.csv", "purpose": "Panels B/C source data and requested OS-score metadata"},
    {"output_file": "GSE149614_MP_similarity_to_discovery.csv", "source": "single_cell/data/processed/cross_cohort_NMF_overlap.csv", "purpose": "Requested complete MP correspondence table"},
    {"output_file": "GSE149614_NMF_program_enrichment.csv", "source": "validation_de_novo_NMF_genes.csv; NMF_eligible_gene_universe.csv; msigdbr 2026.1.Hs", "purpose": "Requested complete ORA source data"},
    {"output_file": "GSE149614_panel_D_pathways.csv", "source": "GSE149614_NMF_program_enrichment.csv", "purpose": "Panel D display subset"},
])
manifest.to_csv(OUT / "source_data_manifest.csv", index=False)

print("VALIDATION_SOURCE_DATA_PREPARED")
print(f"Output directory: {OUT}")
