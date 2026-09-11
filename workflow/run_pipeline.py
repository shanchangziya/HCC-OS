#!/usr/bin/env python3
"""List, preflight, or explicitly execute one HCC-OS revision stage."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]


@dataclass(frozen=True)
class Stage:
    description: str
    commands: tuple[tuple[str, ...], ...]
    settings: tuple[str, ...] = ()
    note: str = ""


PY = "{python}"
RS = "{rscript}"

STAGES: dict[str, Stage] = {
    "single-cell": Stage(
        "Discovery audit and independent GSE149614 assessment (Figures 1-2)",
        (
            (RS, "revision/single_cell/scripts/00_probe_objects.R"),
            (PY, "revision/single_cell/scripts/01_build_validation_atlas.py"),
            (RS, "revision/single_cell/scripts/02_validate_frozen_programs.R"),
            (PY, "revision/single_cell/scripts/03_summarize_validation.py"),
            (RS, "revision/single_cell/scripts/05_validation_de_novo_NMF.R"),
            (PY, "revision/single_cell/scripts/04_plot_scRNA_figures.py"),
            (PY, "revision/single_cell/scripts/07_plot_discovery_definition.py"),
        ),
        ("HCC_OS_DISCOVERY_ROOT", "HCC_OS_GSE149614_COUNTS", "HCC_OS_GSE149614_METADATA_CSV", "HCC_OS_GSE149614_METADATA", "HCC_OS_GSE149614_RDS"),
    ),
    "single-cell-review": Stage(
        "Reviewer-focused independent GeneNMF/GSVA/ORA sensitivity",
        tuple(
            ((RS if name.endswith(".R") else PY), f"revision/single_cell_review/scripts/{name}")
            for name in (
                "01_run_validation_nmf_enrichment.R", "02_prepare_validation_source_data.py",
                "03_plot_validation_suppfig.py", "04_write_validation_summary.py",
                "05_plot_panel_B_review_v2.py", "06_export_genenmf_input.py",
                "07_refit_genenmf_and_select_nmp.R", "08_plot_panel_C_genenmf_heatmap.py",
                "09_score_mp_gsva_and_enrich_selected.R", "10_plot_panel_D_gsva_enrichment.py",
                "11_compute_plot_panel_E_cross_cohort_similarity.py", "12_plot_panel_A_compact.py",
            )
        ),
        note="Requires the completed single-cell stage; script 07 contains GPL-3.0-only adapted GeneNMF code.",
    ),
    "model-input-export": Stage(
        "Export frozen model/cohort inputs from archived R objects",
        (
            (RS, "revision/model_validation/scripts/01_inspect_sources.R"),
            (RS, "revision/model_validation/scripts/02_export_analysis_inputs.R"),
            (RS, "revision/model_validation/scripts/03_export_external_sources.R"),
            (PY, "revision/model_validation/scripts/03b_extract_unfiltered_signature_expression.py"),
        ),
        ("HCC_OS_DISCOVERY_ROOT", "HCC_OS_EXTERNAL_DATA_ROOT", "HCC_OS_ICGC_FPKM"),
    ),
    "model-selection-tcga": Stage(
        "Run and lock all 117 configurations using TCGA only",
        (("bash", "revision/model_selection/scripts/run_server_tcga_development.sh"),),
        ("HCC_OS_MODEL_RAW_DIR",),
        "The input directory must not contain ICGC files; this stage is intentionally separate from ICGC evaluation.",
    ),
    "model-selection-icgc": Stage(
        "Evaluate the already locked TCGA-selected model retrospectively in ICGC",
        (
            (RS, "revision/model_selection/scripts/05_icgc_retrospective_evaluation.R"),
            (RS, "revision/model_selection/scripts/08_export_locked_plsrcox_formula.R"),
            (RS, "revision/model_selection/scripts/10_compare_top7_models.R"),
            (RS, "revision/model_selection/scripts/06_plot_supplementary_figure.R"),
            (RS, "revision/model_selection/scripts/09_plot_compact_panel_A.R"),
            (RS, "revision/model_selection/scripts/11_plot_top7_cindex_differences.R"),
        ),
        ("HCC_OS_MODEL_RAW_DIR",),
        "Run only after the TCGA lock exists and the ICGC input has been staged.",
    ),
    "model-validation": Stage(
        "Frozen score, survival, comparator, and clinical analyses (Figures 3-4)",
        (
            (RS, "revision/model_validation/scripts/11_frozen_model_audit.R"),
            (PY, "revision/model_validation/scripts/run_statistics.py"),
            (PY, "revision/model_validation/scripts/12_reviewer_grade_and_event_audit.py"),
            (PY, "revision/model_validation/published_signatures/scripts/run_published_signature_head_to_head.py"),
            (PY, "revision/model_validation/published_signatures/scripts/plot_head_to_head_cindex.py"),
            (PY, "revision/model_validation/scripts/plot_revision_panels.py"),
            (PY, "revision/model_validation/scripts/plot_assemble_revision.py"),
        ),
        ("HCC_OS_FROZEN_RESULT", "HCC_OS_ICGC_FPKM", "HCC_OS_LEGACY_PROJECT_ROOT", "HCC_OS_TCGA_CLINICAL_GRADE"),
    ),
    "bulk-support": Stage(
        "MP4/ROS projection, pathway biology, and compact sensitivity panels",
        (
            (PY, "revision/osars_mp4_projection/scripts/run_osars_os_high_gsea.py"),
            (PY, "revision/osars_mp4_projection/scripts/plot_osars_os_high_gsea.py"),
            (RS, "revision/osars_mp4_projection/multimethod_projection/scripts/run_multimethod_state_projection.R"),
            (PY, "revision/osars_mp4_projection/multimethod_projection/scripts/plot_multimethod_state_projection.py"),
            (RS, "revision/bulk_biology/scripts/01_bulk_biology_enrichment.R"),
            (PY, "revision/bulk_biology/scripts/02_plot_bulk_biology_figure.py"),
            (PY, "revision/osars_ros_correlation/scripts/plot_osars_direct_ros_correlation.py"),
            (PY, "revision/osars_ros_overlap/scripts/plot_venn.py"),
        ),
        note="Requires model-validation processed tables; negative MP4 direction must be retained in reporting.",
    ),
    "clinical-sensitivities": Stage(
        "Stage/ROS discrimination and stage/grade Cox sensitivities",
        (
            (RS, "revision/clinical_increment/scripts/01_fit_tcga_cox_models.R"),
            (PY, "revision/clinical_increment/scripts/02_summarize_and_plot.py"),
            (RS, "revision/clinical_cox/scripts/run_stage_grade_cox.R"),
            (PY, "revision/osars_discrimination/scripts/plot_osars_cindex_bars.py"),
        ),
        note="Requires completed model-validation processed tables.",
    ),
    "immune-therapy": Stage(
        "Inferred immune/TIDE and exploratory treatment contexts (Figure 5)",
        (
            (RS, "revision/immune_therapy/scripts/export_r_objects.R"),
            (PY, "revision/immune_therapy/scripts/prepare_data.py"),
            (PY, "revision/immune_therapy/scripts/verify_source_labels.py"),
            (PY, "revision/immune_therapy/scripts/plot_panels.py"),
            (PY, "revision/immune_therapy/scripts/assemble_and_qc.py"),
        ),
        ("HCC_OS_IMMUNE_INPUT_DIR", "HCC_OS_LEGACY_PROJECT_ROOT"),
    ),
    "spatial-association": Stage(
        "Validate staged cell2location exports and render descriptive Figure 6",
        (
            (PY, "revision/spatial_association/scripts/01_build_source_manifest.py"),
            (PY, "revision/spatial_association/scripts/02_export_descriptive_spatial.py"),
            (PY, "revision/spatial_association/scripts/plot_fig6_panels.py"),
            (PY, "revision/spatial_association/scripts/plot_assemble_fig6.py"),
        ),
        note="Restricted inputs must follow the relative data/raw layout in docs/DATA_MANIFEST.md.",
    ),
    "protein-support": Stage(
        "PDC000198 mapping, paired NQO1, survival, and supplementary support",
        (
            (PY, "revision/protein_support/scripts/run_protein_support.py"),
            (PY, "revision/protein_support/scripts/audit_protein_outputs.py"),
            (PY, "revision/protein_support/scripts/plot_nqo1_panels.py"),
            (PY, "revision/protein_support/scripts/plot_assemble_nqo1.py"),
        ),
        ("HCC_OS_PROTEOMICS_ROOT", "HCC_OS_PDC_SURROGATE_WEIGHTS", "HCC_OS_PDC_SURROGATE_SCORES"),
    ),
    "pathology": Stage(
        "Corrected internal TCGA pathology reassessment",
        (
            (PY, "revision/pathology/scripts/run_pathology_analysis.py"),
            (PY, "revision/pathology/scripts/audit_processed.py"),
            (PY, "revision/pathology/scripts/plot_revision_panels.py"),
            (PY, "revision/pathology/scripts/plot_assemble_revision.py"),
        ),
        ("HCC_OS_PATHOLOGY_SOURCE_DIR", "HCC_OS_PATHOLOGY_CLINICAL"),
    ),
    "nqo1-prioritization": Stage(
        "Post-model NQO1 prioritization and paired-protein supplementary figure",
        (
            (RS, "revision/nqo1_prioritization/scripts/01_prepare_nqo1_priority_evidence.R"),
            (PY, "revision/nqo1_prioritization/scripts/02_plot_nqo1_priority_protein.py"),
            (PY, "revision/nqo1_prioritization/scripts/03_plot_nature_supplementary_figure.py"),
        ),
        note="Requires completed protein-support outputs.",
    ),
    "manuscript-tables": Stage(
        "Assemble aggregate manuscript-support tables from audited outputs",
        ((PY, "revision/scripts/build_summary_tables.py"),),
        note="Run after the single-cell, model-validation, protein-support, and pathology outputs are complete.",
    ),
}


def resolve(command: tuple[str, ...]) -> list[str]:
    return [
        token.replace(PY, sys.executable).replace(RS, os.environ.get("HCC_OS_RSCRIPT", "Rscript"))
        for token in command
    ]


def structural_check() -> int:
    problems: list[str] = []
    for name, stage in STAGES.items():
        for command in stage.commands:
            resolved = resolve(command)
            executable = resolved[0]
            is_path = Path(executable).is_absolute() or "/" in executable or "\\" in executable
            if is_path:
                if not Path(executable).expanduser().is_file():
                    problems.append(f"{name}: executable path not found: {executable}")
            elif shutil.which(executable) is None:
                problems.append(f"{name}: executable not found: {executable}")
            for token in resolved[1:]:
                if token.startswith("revision/") and not (ROOT / token).is_file():
                    problems.append(f"{name}: entry file missing: {token}")
    required_docs = ["docs/DATA_MANIFEST.md", "docs/FIGURE_CODE_MAP.md", "config/paths.example.env"]
    problems.extend(f"required file missing: {path}" for path in required_docs if not (ROOT / path).is_file())
    if problems:
        print("Preflight: FAIL")
        print("\n".join(f"- {problem}" for problem in problems))
        return 1
    print(f"Preflight: PASS ({len(STAGES)} stage entry points checked)")
    return 0


def show_stage(name: str, execute: bool) -> int:
    stage = STAGES[name]
    print(f"[{name}] {stage.description}")
    if stage.settings:
        for setting in stage.settings:
            value = os.environ.get(setting)
            print(f"  {setting}: {'set' if value else 'not set / use documented local default'}")
    if stage.note:
        print(f"  Note: {stage.note}")
    for command in stage.commands:
        resolved = resolve(command)
        print("  $", shlex.join(resolved))
        if execute:
            subprocess.run(resolved, cwd=ROOT, check=True)
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--list", action="store_true", help="List stages")
    parser.add_argument("--check", action="store_true", help="Check entry files and executables")
    parser.add_argument("--stage", choices=sorted(STAGES), help="Show or run one stage")
    parser.add_argument("--execute", action="store_true", help="Actually execute the selected stage")
    args = parser.parse_args()
    if args.execute and not args.stage:
        parser.error("--execute requires --stage")
    if args.list:
        for name, stage in STAGES.items():
            print(f"{name:24} {stage.description}")
    status = structural_check() if args.check else 0
    if args.stage:
        status = max(status, show_stage(args.stage, args.execute))
    if not (args.list or args.check or args.stage):
        parser.print_help()
    return status


if __name__ == "__main__":
    raise SystemExit(main())
