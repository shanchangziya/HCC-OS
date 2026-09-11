#!/usr/bin/env python3
"""Preranked GSEA of the frozen discovery MP4 program in OSARS risk groups.

The question is deliberately narrow: are the transcriptional features of the
frozen single-cell OS-high state enriched in bulk OSARS-high tumors?

The primary test compares the frozen median-split OSARS high and low groups in
each frozen cohort. A predeclared stage-adjusted sensitivity residualizes every
gene on categorical pathological stage before forming the same high-minus-low
rank statistic. No model, gene set, risk cutoff, or expression value is refit
or changed by this script.

The GSEA null is a preranked competitive gene-set permutation null: 10,000
random gene sets with the same observed gene-set size are sampled from each
cohort's tested expression universe. Nominal sign-specific empirical P values
are BH-adjusted across the three fixed gene sets in both cohorts, separately
for the primary and stage-adjusted analyses. The state-specific sensitivity
also removes genes used as OSARS model features, avoiding a direct
gene-component overlap between the tested state program and the classifier.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
REVISION = ROOT.parent
MODEL_VALIDATION = REVISION / "model_validation"
SINGLE_CELL = REVISION / "single_cell"
RAW = MODEL_VALIDATION / "data" / "raw"
MODEL_PROCESSED = MODEL_VALIDATION / "data" / "processed"
OUT = ROOT / "data" / "processed"
REFERENCES = ROOT / "references"
TABLES = ROOT / "tables"
FROZEN_OSARS_FEATURES = REVISION / "tables" / "Table2_frozen_OSARS_128_genes.csv"

SEED = 20260906
PERMUTATIONS = 10_000
WEIGHT_POWER = 1.0


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def bh_adjust(p_values: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjustment for a finite vector of P values."""
    p_values = np.asarray(p_values, dtype=float)
    order = np.argsort(p_values)
    ranked = p_values[order]
    m = len(ranked)
    adjusted = ranked * m / np.arange(1, m + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result = np.empty(m, dtype=float)
    result[order] = np.minimum(adjusted, 1.0)
    return result


def welch_t(values: np.ndarray, high: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Gene-wise Welch t ranks for OSARS-high minus OSARS-low expression."""
    high_values = values[:, high]
    low_values = values[:, ~high]
    high_mean = high_values.mean(axis=1)
    low_mean = low_values.mean(axis=1)
    standard_error = np.sqrt(
        high_values.var(axis=1, ddof=1) / high_values.shape[1]
        + low_values.var(axis=1, ddof=1) / low_values.shape[1]
    )
    statistic = np.divide(
        high_mean - low_mean,
        standard_error,
        out=np.zeros_like(high_mean),
        where=standard_error > 0,
    )
    statistic[~np.isfinite(statistic)] = 0.0
    return statistic, high_mean, low_mean


def stage_adjusted_t(
    values: np.ndarray, high: np.ndarray, stage: pd.Series
) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[str]]:
    """Per-gene OLS t statistic for high risk after adjustment for stage."""
    labels = stage.astype(str).to_numpy()
    levels = sorted(pd.unique(labels).tolist())
    stage_design = np.column_stack(
        [np.ones(len(labels), dtype=float)]
        + [(labels == level).astype(float) for level in levels[1:]]
    )
    design = np.column_stack([stage_design, high.astype(float)])
    degrees_of_freedom = values.shape[1] - design.shape[1]
    if degrees_of_freedom < 2:
        raise ValueError("Insufficient degrees of freedom for stage-adjusted test")
    coefficients = np.linalg.lstsq(design, values.T, rcond=None)[0]
    residual = values - (design @ coefficients).T
    residual_variance = (residual**2).sum(axis=1) / degrees_of_freedom
    high_variance = float(np.linalg.pinv(design.T @ design)[-1, -1])
    standard_error = np.sqrt(residual_variance * high_variance)
    statistic = np.divide(
        coefficients[-1, :],
        standard_error,
        out=np.zeros(values.shape[0], dtype=float),
        where=standard_error > 0,
    )
    statistic[~np.isfinite(statistic)] = 0.0
    high_mean = values[:, high].mean(axis=1)
    low_mean = values[:, ~high].mean(axis=1)
    return statistic, high_mean, low_mean, levels


def enrichment_from_positions(
    rank_scores: np.ndarray, positions: np.ndarray
) -> tuple[float, np.ndarray, int]:
    """Compute weighted running GSEA score from sorted hit positions."""
    total_genes = len(rank_scores)
    hit_count = len(positions)
    if hit_count < 2 or total_genes - hit_count < 2:
        raise ValueError("Gene set is too small or too large for GSEA")
    weights = np.abs(rank_scores[positions]) ** WEIGHT_POWER
    if weights.sum() == 0:
        weights = np.ones(hit_count, dtype=float)
    hit_increment = weights / weights.sum()
    running = np.full(total_genes, -1.0 / (total_genes - hit_count), dtype=float)
    running[positions] = hit_increment
    running = np.cumsum(running)
    positive_peak = int(np.argmax(running))
    negative_peak = int(np.argmin(running))
    if running[positive_peak] >= abs(running[negative_peak]):
        return float(running[positive_peak]), running, positive_peak
    return float(running[negative_peak]), running, negative_peak


def permutation_null(
    rank_scores: np.ndarray, hit_count: int, seed: int
) -> np.ndarray:
    """Same-size gene-set permutation null for a preranked GSEA."""
    rng = np.random.default_rng(seed)
    total = len(rank_scores)
    null = np.empty(PERMUTATIONS, dtype=float)
    for index in range(PERMUTATIONS):
        positions = np.sort(rng.choice(total, size=hit_count, replace=False))
        null[index] = enrichment_from_positions(rank_scores, positions)[0]
    return null


def normalized_score_and_p(es: float, null: np.ndarray) -> tuple[float, float, float]:
    """Return NES, sign-specific empirical P, and null scaling factor."""
    if es >= 0:
        same_sign = null[null >= 0]
        scale = float(np.mean(same_sign))
        p_value = (1 + int(np.sum(null >= es))) / (1 + len(null))
    else:
        same_sign = -null[null < 0]
        scale = float(np.mean(same_sign))
        p_value = (1 + int(np.sum(null <= es))) / (1 + len(null))
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError("Invalid same-sign GSEA null scaling factor")
    return es / scale, p_value, scale


def gene_sets() -> dict[str, set[str]]:
    """Read frozen state signatures, without adapting them to bulk results."""
    nmf = pd.read_csv(SINGLE_CELL / "data" / "processed" / "discovery_NMF_genes.csv")
    mp4 = set(nmf.loc[nmf["program"] == "MP4", "gene"].astype(str))
    osars_features = set(pd.read_csv(FROZEN_OSARS_FEATURES)["gene"].astype(str))
    ros = set(
        (
            MODEL_VALIDATION
            / "references"
            / "HALLMARK_ROS_msigdb_v7.0_genes.txt"
        ).read_text().split()
    )
    return {
        "Frozen discovery MP4 state": mp4,
        "Frozen discovery MP4 state excluding ROS and OSARS features": (
            mp4.difference(ros).difference(osars_features)
        ),
        "Hallmark ROS positive control": ros,
    }


def load_expression_and_metadata(cohort: str) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Load one frozen expression matrix, risk groups, and stage table."""
    expression_file = (
        RAW / "TCGA_full_expression_logscale.csv.gz"
        if cohort == "TCGA-LIHC"
        else RAW / "ICGC_full_expression_logscale.csv.gz"
    )
    score_file = MODEL_PROCESSED / f"{cohort}_all_scores.csv"
    clinical_file = MODEL_PROCESSED / f"{cohort}_clinical_analysis_rows.csv"
    expression = pd.read_csv(expression_file, index_col=0)
    scores = pd.read_csv(score_file).set_index("ID")
    clinical = pd.read_csv(clinical_file).set_index("ID")
    scores.index = scores.index.astype(str)
    clinical.index = clinical.index.astype(str)
    expression.columns = expression.columns.astype(str)
    if expression.index.has_duplicates or expression.columns.has_duplicates:
        raise ValueError(f"Duplicate gene or sample identifier in {expression_file.name}")
    if scores.index.has_duplicates or clinical.index.has_duplicates:
        raise ValueError(f"Duplicate clinical/sample identifier in {cohort}")
    absent = scores.index.difference(expression.columns)
    if len(absent):
        raise ValueError(f"{cohort} expression misses frozen score IDs: {absent[:5].tolist()}")
    if set(scores.index) != set(clinical.index):
        raise ValueError(f"{cohort} score and clinical ID sets differ")
    clinical = clinical.loc[scores.index]
    if set(scores["OSARS_frozen_group"].dropna().unique()) != {"High", "Low"}:
        raise ValueError(f"Unexpected frozen OSARS group labels in {cohort}")
    expression = expression.loc[:, scores.index]
    if not np.isfinite(expression.to_numpy(float)).all():
        raise ValueError(f"Non-finite frozen expression in {cohort}")
    return expression, scores, clinical


def run_one(
    cohort: str,
    analysis: str,
    expression: pd.DataFrame,
    scores: pd.DataFrame,
    clinical: pd.DataFrame,
    sets: dict[str, set[str]],
) -> tuple[list[dict], list[pd.DataFrame], pd.DataFrame, pd.DataFrame]:
    """Run one cohort/analysis GSEA and return all audit-ready records."""
    if analysis == "primary":
        selected_ids = scores.index
        values = expression.to_numpy(float)
        stage_levels = ""
    elif analysis == "stage_adjusted":
        stage = clinical["stage"]
        selected_ids = stage.index[stage.notna()]
        if not len(selected_ids):
            raise ValueError(f"No stage-complete samples in {cohort}")
        values = expression.loc[:, selected_ids].to_numpy(float)
        stage_levels = ""
    else:
        raise ValueError(f"Unknown analysis: {analysis}")

    high = scores.loc[selected_ids, "OSARS_frozen_group"].eq("High").to_numpy()
    if high.sum() < 5 or (~high).sum() < 5:
        raise ValueError(f"Insufficient risk-group sample size in {cohort}/{analysis}")
    if analysis == "primary":
        statistic, high_mean, low_mean = welch_t(values, high)
        rank_statistic = "Welch t statistic for OSARS-high minus OSARS-low"
    else:
        statistic, high_mean, low_mean, stage_levels_list = stage_adjusted_t(
            values, high, clinical.loc[selected_ids, "stage"]
        )
        stage_levels = ";".join(stage_levels_list)
        rank_statistic = (
            "OLS t statistic for OSARS-high coefficient, adjusted for categorical stage"
        )
    order = np.argsort(-statistic, kind="stable")
    ranked_genes = expression.index.to_numpy()[order]
    rank_scores = statistic[order]
    ranked = pd.DataFrame(
        {
            "cohort": cohort,
            "analysis": analysis,
            "rank": np.arange(1, len(ranked_genes) + 1),
            "gene": ranked_genes,
            "rank_statistic_high_minus_low": rank_scores,
            "mean_expression_high": high_mean[order],
            "mean_expression_low": low_mean[order],
        }
    )
    result_rows: list[dict] = []
    coordinate_frames: list[pd.DataFrame] = []
    coverage_rows: list[dict] = []
    for set_position, (set_name, genes) in enumerate(sets.items(), start=1):
        hit = np.isin(ranked_genes, sorted(genes))
        positions = np.flatnonzero(hit)
        es, running, peak = enrichment_from_positions(rank_scores, positions)
        null = permutation_null(
            rank_scores,
            len(positions),
            seed=SEED + (1_000 if cohort == "ICGC-LIRI" else 0)
            + (10_000 if analysis == "stage_adjusted" else 0)
            + set_position,
        )
        nes, p_value, null_scale = normalized_score_and_p(es, null)
        if es >= 0:
            leading = ranked_genes[(np.arange(len(ranked_genes)) <= peak) & hit]
            direction = "OSARS-high"
        else:
            leading = ranked_genes[(np.arange(len(ranked_genes)) >= peak) & hit]
            direction = "OSARS-low"
        result_rows.append(
            {
                "cohort": cohort,
                "analysis": analysis,
                "gene_set": set_name,
                "rank_contrast": "OSARS-high minus OSARS-low",
                "rank_statistic": rank_statistic,
                "n_total": len(selected_ids),
                "n_high": int(high.sum()),
                "n_low": int((~high).sum()),
                "stage_levels_adjusted": stage_levels,
                "ranking_genes": len(ranked_genes),
                "genes_defined": len(genes),
                "genes_used": int(hit.sum()),
                "ES": es,
                "NES": nes,
                "nominal_empirical_p": p_value,
                "null_same_sign_mean_abs_ES": null_scale,
                "permutations": PERMUTATIONS,
                "permutation_scheme": "same-size random gene sets from tested ranked universe",
                "enriched_at": direction,
                "peak_rank": peak + 1,
                "peak_rank_fraction": (peak + 1) / len(ranked_genes),
                "leading_edge_count": len(leading),
                "leading_edge_genes": ";".join(leading.tolist()),
            }
        )
        coordinate_frames.append(
            pd.DataFrame(
                {
                    "cohort": cohort,
                    "analysis": analysis,
                    "gene_set": set_name,
                    "rank": np.arange(1, len(ranked_genes) + 1),
                    "rank_fraction": np.arange(1, len(ranked_genes) + 1)
                    / len(ranked_genes),
                    "gene": ranked_genes,
                    "rank_statistic_high_minus_low": rank_scores,
                    "is_gene_set_hit": hit,
                    "running_ES": running,
                    "is_peak": np.arange(len(ranked_genes)) == peak,
                }
            )
        )
        coverage_rows.extend(
            {
                "cohort": cohort,
                "analysis": analysis,
                "gene_set": set_name,
                "gene": gene,
                "available_in_ranked_universe": gene in set(ranked_genes),
            }
            for gene in sorted(genes)
        )
    return result_rows, coordinate_frames, ranked, pd.DataFrame(coverage_rows)


def main() -> None:
    for directory in (OUT, REFERENCES, TABLES):
        directory.mkdir(parents=True, exist_ok=True)
    manifest_inputs: dict[str, str] = {}
    sets = gene_sets()
    for file in (
        SINGLE_CELL / "data" / "processed" / "discovery_NMF_genes.csv",
        MODEL_VALIDATION / "references" / "HALLMARK_ROS_msigdb_v7.0_genes.txt",
        FROZEN_OSARS_FEATURES,
    ):
        manifest_inputs[file.name] = sha256(file)
    pd.DataFrame(
        [
            {"gene_set": name, "gene": gene}
            for name, genes in sets.items()
            for gene in sorted(genes)
        ]
    ).to_csv(REFERENCES / "frozen_gene_set_definitions.csv", index=False)

    all_results: list[dict] = []
    all_coordinates: list[pd.DataFrame] = []
    all_coverage: list[pd.DataFrame] = []
    stage_rows: list[dict] = []
    for cohort in ("TCGA-LIHC", "ICGC-LIRI"):
        expression, scores, clinical = load_expression_and_metadata(cohort)
        expression_file = (
            RAW / "TCGA_full_expression_logscale.csv.gz"
            if cohort == "TCGA-LIHC"
            else RAW / "ICGC_full_expression_logscale.csv.gz"
        )
        manifest_inputs[expression_file.name] = sha256(expression_file)
        for file in (
            MODEL_PROCESSED / f"{cohort}_all_scores.csv",
            MODEL_PROCESSED / f"{cohort}_clinical_analysis_rows.csv",
        ):
            manifest_inputs[file.name] = sha256(file)
        stage_table = (
            clinical.assign(OSARS_group=scores["OSARS_frozen_group"])
            .groupby(["OSARS_group", "stage"], dropna=False)
            .size()
            .rename("n")
            .reset_index()
        )
        stage_table.insert(0, "cohort", cohort)
        stage_rows.append(stage_table)
        for analysis in ("primary", "stage_adjusted"):
            results, coordinates, ranked, coverage = run_one(
                cohort, analysis, expression, scores, clinical, sets
            )
            all_results.extend(results)
            all_coordinates.extend(coordinates)
            all_coverage.append(coverage)
            ranked.to_csv(
                OUT / f"{cohort}_{analysis}_ranked_gene_statistics.csv.gz",
                index=False,
                compression="gzip",
            )

    results = pd.DataFrame(all_results)
    results["BH_q_within_analysis_6tests"] = np.nan
    for analysis, index in results.groupby("analysis").groups.items():
        indices = list(index)
        results.loc[indices, "BH_q_within_analysis_6tests"] = bh_adjust(
            results.loc[indices, "nominal_empirical_p"].to_numpy(float)
        )
    results.to_csv(TABLES / "GSEA_results.csv", index=False)
    pd.concat(all_coordinates, ignore_index=True).to_csv(
        OUT / "GSEA_running_score_coordinates.csv.gz",
        index=False,
        compression="gzip",
    )
    pd.concat(all_coverage, ignore_index=True).to_csv(
        OUT / "gene_set_coverage.csv", index=False
    )
    pd.concat(stage_rows, ignore_index=True).to_csv(
        OUT / "stage_by_OSARS_group.csv", index=False
    )

    manifest = {
        "analysis": "Preranked GSEA of frozen discovery MP4 state in OSARS risk groups",
        "created": "2026-09-06",
        "risk_groups": "Frozen cohort-median OSARS_frozen_group labels; no cutoff reoptimization",
        "rank_statistic": "Gene-wise Welch t statistic, OSARS-high minus OSARS-low",
        "primary": "All frozen cohort samples",
        "stage_adjusted_sensitivity": "Fit each gene as expression ~ frozen OSARS group + categorical stage, then rank by the adjusted high-group t statistic",
        "gene_sets": {name: len(genes) for name, genes in sets.items()},
        "permutations": PERMUTATIONS,
        "null": "Same-size random gene-set permutations from tested ranked gene universe",
        "multiple_testing": "Benjamini-Hochberg within the 6 tests per analysis",
        "raw_inputs_sha256": manifest_inputs,
        "state_gene_sources": {
            "MP4": str(SINGLE_CELL / "data" / "processed" / "discovery_NMF_genes.csv"),
            "Hallmark_ROS": str(
                MODEL_VALIDATION
                / "references"
                / "HALLMARK_ROS_msigdb_v7.0_genes.txt"
            ),
            "OSARS_features": str(FROZEN_OSARS_FEATURES),
        },
    }
    (OUT / "analysis_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(results.to_string(index=False))


if __name__ == "__main__":
    main()
