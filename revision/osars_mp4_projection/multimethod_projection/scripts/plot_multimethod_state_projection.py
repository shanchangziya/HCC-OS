#!/usr/bin/env python3
"""Render figures from the processed multi-method state-projection outputs.

This script reads only the previously calculated score and statistics tables.
It does not calculate tests, risk groups, or gene-set scores.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "processed"
TABLES = ROOT / "tables"
FIGURES = ROOT / "figures"

SCORES_FILE = DATA / "per_sample_state_projection_scores.csv"
TRENDS_FILE = DATA / "state_projection_linear_trend_coordinates.csv"
STATS_FILE = TABLES / "state_projection_method_statistics.csv"
CORRELATIONS_FILE = TABLES / "state_projection_continuous_correlations.csv"

METHODS_MAIN = ["rank_module", "GSVA", "ssGSEA", "PLAGE"]
METHODS_FOREST = ["rank_module", "mean_gene_z", "GSVA", "ssGSEA", "PLAGE"]
METHOD_LABELS = {
    "rank_module": "Within-sample rank module",
    "mean_gene_z": "Mean gene z-score",
    "GSVA": "GSVA",
    "ssGSEA": "ssGSEA",
    "PLAGE": "PLAGE",
    "GSVA_zscore": "GSVA z-score",
    "PCA_PC1": "PCA PC1",
}
LOW = "#0072B2"
HIGH = "#D55E00"
AXIS = "#333333"
GRAY = "#6E6E6E"
LIGHT_GRAY = "#C9C9C9"

mpl.rcParams.update(
    {
        "font.family": "Arial",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 7,
        "axes.titlesize": 8,
        "axes.labelsize": 8,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "axes.linewidth": 0.75,
        "xtick.major.width": 0.75,
        "ytick.major.width": 0.75,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "savefig.dpi": 600,
    }
)


def q_label(value: float) -> str:
    if value < 0.001:
        return "BH q < 0.001"
    return f"BH q = {value:.3f}"


def base_axis(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(color=AXIS, labelcolor=AXIS)
    ax.spines["left"].set_color(AXIS)
    ax.spines["bottom"].set_color(AXIS)


def save_figure(fig: plt.Figure, basename: str) -> dict[str, str]:
    FIGURES.mkdir(parents=True, exist_ok=True)
    pdf = FIGURES / f"{basename}.pdf"
    png = FIGURES / f"{basename}.png"
    fig.savefig(pdf, bbox_inches="tight", pad_inches=0.02)
    fig.savefig(png, dpi=600, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    return {"pdf": pdf.name, "png": png.name}


def draw_group_comparisons(scores: pd.DataFrame, stats: pd.DataFrame) -> dict[str, str]:
    """TCGA distribution panels using individual observations and precomputed q values."""
    cohort = "TCGA-LIHC"
    subset = scores.loc[scores["cohort"].eq(cohort)].copy()
    fig, axes = plt.subplots(2, 2, figsize=(7.0, 4.65))
    rng = np.random.default_rng(20260906)
    for ax, method in zip(axes.flat, METHODS_MAIN):
        values = [
            subset.loc[subset["OSARS_group"].eq(group), method].to_numpy(float)
            for group in ("Low", "High")
        ]
        violin = ax.violinplot(
            values,
            positions=[0, 1],
            widths=0.7,
            showmeans=False,
            showmedians=False,
            showextrema=False,
        )
        for body, color in zip(violin["bodies"], (LOW, HIGH)):
            body.set_facecolor(color)
            body.set_edgecolor(color)
            body.set_alpha(0.22)
        box = ax.boxplot(
            values,
            positions=[0, 1],
            widths=0.35,
            patch_artist=True,
            showfliers=False,
            medianprops={"color": AXIS, "linewidth": 1.0},
            whiskerprops={"color": AXIS, "linewidth": 0.75},
            capprops={"color": AXIS, "linewidth": 0.75},
        )
        for patch, color in zip(box["boxes"], (LOW, HIGH)):
            patch.set_facecolor("white")
            patch.set_edgecolor(color)
            patch.set_linewidth(0.9)
        for position, (score, color) in enumerate(zip(values, (LOW, HIGH))):
            jitter = rng.uniform(-0.16, 0.16, len(score))
            ax.scatter(
                np.full(len(score), position) + jitter,
                score,
                s=5,
                color=color,
                alpha=0.32,
                linewidths=0,
                rasterized=True,
            )
        stat = stats.loc[
            (stats["cohort"].eq(cohort)) & (stats["method"].eq(method))
        ].iloc[0]
        lower = min(np.min(values[0]), np.min(values[1]))
        upper = max(np.max(values[0]), np.max(values[1]))
        span = upper - lower if upper > lower else 1.0
        ax.set_ylim(lower - 0.05 * span, upper + 0.23 * span)
        annotation = (
            f"{q_label(float(stat['wilcoxon_BH_q_within_cohort']))}\n"
            f"d = {float(stat['standardized_mean_difference_high_minus_low']):.2f} "
            f"[{float(stat['standardized_mean_difference_ci_lower']):.2f}, "
            f"{float(stat['standardized_mean_difference_ci_upper']):.2f}]"
        )
        ax.text(
            0.5,
            upper + 0.18 * span,
            annotation,
            ha="center",
            va="top",
            fontsize=6.5,
            color=AXIS,
        )
        ax.set_title(METHOD_LABELS[method], pad=4)
        ax.set_xticks([0, 1], ["OSARS-low", "OSARS-high"])
        ax.set_ylabel("Strict MP4 state score")
        base_axis(ax)
    fig.suptitle(
        "TCGA-LIHC: frozen MP4 state projection by frozen OSARS group",
        x=0.5,
        y=0.995,
        fontsize=9,
        fontweight="bold",
    )
    fig.text(
        0.5,
        0.005,
        "Fixed 50-gene MP4 program excluding ROS and OSARS-feature overlap; points are individual tumors.",
        ha="center",
        va="bottom",
        fontsize=6.5,
        color=GRAY,
    )
    fig.tight_layout(rect=(0, 0.045, 1, 0.94), h_pad=2.0, w_pad=1.9)
    return save_figure(fig, "R1Q9_TCGA_multimethod_state_scores_by_OSARS_group")


def draw_continuous_associations(
    scores: pd.DataFrame, trends: pd.DataFrame, correlations: pd.DataFrame
) -> dict[str, str]:
    """TCGA continuous risk versus projection-score scatter panels."""
    cohort = "TCGA-LIHC"
    subset = scores.loc[scores["cohort"].eq(cohort)].copy()
    fig, axes = plt.subplots(2, 2, figsize=(7.0, 4.65))
    for ax, method in zip(axes.flat, METHODS_MAIN):
        low = subset.loc[subset["OSARS_group"].eq("Low")]
        high = subset.loc[subset["OSARS_group"].eq("High")]
        ax.scatter(
            low["OSARS_frozen"], low[method], s=9, color=LOW, alpha=0.48,
            linewidths=0, rasterized=True, label="OSARS-low"
        )
        ax.scatter(
            high["OSARS_frozen"], high[method], s=9, color=HIGH, alpha=0.48,
            linewidths=0, rasterized=True, label="OSARS-high"
        )
        line = trends.loc[
            (trends["cohort"].eq(cohort)) & (trends["method"].eq(method))
        ].sort_values("OSARS_frozen")
        ax.plot(
            line["OSARS_frozen"],
            line["fitted_state_score"],
            color=AXIS,
            linewidth=1.0,
            zorder=3,
        )
        correlation = correlations.loc[
            (correlations["cohort"].eq(cohort))
            & (correlations["method"].eq(method))
        ].iloc[0]
        annotation = (
            f"Spearman rho = {float(correlation['spearman_rho_OSARS_vs_state_score']):.2f}\n"
            f"95% CI [{float(correlation['spearman_ci_lower']):.2f}, "
            f"{float(correlation['spearman_ci_upper']):.2f}]\n"
            f"{q_label(float(correlation['spearman_BH_q_within_cohort']))}"
        )
        ax.text(
            0.03,
            0.97,
            annotation,
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=6.2,
            color=AXIS,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.78, "pad": 1.2},
        )
        ax.set_title(METHOD_LABELS[method], pad=4)
        ax.set_xlabel("Frozen OSARS score")
        ax.set_ylabel("Strict MP4 state score")
        base_axis(ax)
    axes.flat[0].legend(
        frameon=False,
        loc="lower left",
        fontsize=6.2,
        handletextpad=0.3,
        borderpad=0.1,
    )
    fig.suptitle(
        "TCGA-LIHC: continuous association of OSARS with frozen MP4 state scores",
        x=0.5,
        y=0.995,
        fontsize=9,
        fontweight="bold",
    )
    fig.text(
        0.5,
        0.005,
        "Lines are descriptive least-squares fits; inference is the labeled Spearman correlation.",
        ha="center",
        va="bottom",
        fontsize=6.5,
        color=GRAY,
    )
    fig.tight_layout(rect=(0, 0.045, 1, 0.94), h_pad=2.0, w_pad=1.9)
    return save_figure(fig, "R1Q9_TCGA_continuous_OSARS_state_projection")


def draw_effect_forest(stats: pd.DataFrame) -> dict[str, str]:
    """Replication forest plot for standardized High-minus-Low score differences."""
    fig, axes = plt.subplots(1, 2, figsize=(6.8, 3.15), sharey=True)
    y_positions = np.arange(len(METHODS_FOREST))[::-1]
    for ax, cohort in zip(axes, ("TCGA-LIHC", "ICGC-LIRI")):
        subset = (
            stats.loc[
                (stats["cohort"].eq(cohort))
                & (stats["method"].isin(METHODS_FOREST))
            ]
            .set_index("method")
            .loc[METHODS_FOREST]
            .reset_index()
        )
        for y, row in zip(y_positions, subset.itertuples(index=False)):
            effect = float(row.standardized_mean_difference_high_minus_low)
            lower = float(row.standardized_mean_difference_ci_lower)
            upper = float(row.standardized_mean_difference_ci_upper)
            q = float(row.stage_adjusted_BH_q_within_cohort)
            color = LOW if effect < 0 else HIGH
            ax.errorbar(
                effect,
                y,
                xerr=[[effect - lower], [upper - effect]],
                fmt="o",
                color=color,
                markerfacecolor=color if q < 0.05 else "white",
                markeredgecolor=color,
                markersize=5.4,
                capsize=2,
                capthick=0.8,
                linewidth=0.9,
                zorder=3,
            )
            ax.text(
                0.12,
                y,
                q_label(q),
                ha="left",
                va="center",
                fontsize=5.7,
                color=color if q < 0.05 else GRAY,
            )
        ax.axvline(0, color=LIGHT_GRAY, linestyle=(0, (2, 2)), linewidth=0.8, zorder=0)
        ax.set_xlim(-1.55, 0.55)
        ax.set_title(cohort, fontweight="bold", pad=5)
        ax.set_xlabel("Standardized mean difference\n(OSARS-high minus OSARS-low)")
        ax.set_yticks(y_positions, [METHOD_LABELS[method] for method in METHODS_FOREST])
        base_axis(ax)
    fig.text(
        0.5,
        0.005,
        "Error bars: 2,000-sample bootstrap 95% CIs. Filled points: stage-adjusted BH q < 0.05.",
        ha="center",
        va="bottom",
        fontsize=6.2,
        color=GRAY,
    )
    fig.suptitle(
        "Frozen MP4 state score is lower in OSARS-high tumors across scoring frameworks",
        y=0.99,
        fontsize=9,
        fontweight="bold",
    )
    fig.tight_layout(rect=(0, 0.08, 1, 0.91), w_pad=1.7)
    return save_figure(fig, "R1Q9_multimethod_state_projection_effects_TCGA_ICGC")


def main() -> None:
    for file in (SCORES_FILE, TRENDS_FILE, STATS_FILE, CORRELATIONS_FILE):
        if not file.exists():
            raise FileNotFoundError(f"Required processed input is missing: {file}")
    scores = pd.read_csv(SCORES_FILE)
    trends = pd.read_csv(TRENDS_FILE)
    stats = pd.read_csv(STATS_FILE)
    correlations = pd.read_csv(CORRELATIONS_FILE)
    expected = {"TCGA-LIHC", "ICGC-LIRI"}
    if set(scores["cohort"]) != expected or set(stats["cohort"]) != expected:
        raise ValueError("Unexpected cohort set in processed multi-method outputs")
    if not set(METHODS_FOREST).issubset(set(stats["method"])):
        raise ValueError("Expected state-projection methods are absent from statistics")
    manifest = {
        "source_scores": str(SCORES_FILE.relative_to(ROOT)),
        "source_statistics": str(STATS_FILE.relative_to(ROOT)),
        "source_correlations": str(CORRELATIONS_FILE.relative_to(ROOT)),
        "source_trends": str(TRENDS_FILE.relative_to(ROOT)),
        "strict_state_gene_set": "Frozen MP4 genes excluding Hallmark ROS and OSARS-feature overlap",
        "plotted_methods": METHODS_FOREST,
        "note": (
            "PLAGE and PCA PC1 are equivalent single-set latent-factor formulations; "
            "mean gene z-score and GSVA z-score are closely related. The figures show "
            "one representative of each redundant pair, while all seven methods remain "
            "available in the source tables."
        ),
        "figures": {},
    }
    for name, function in (
        ("R1Q9_TCGA_multimethod_state_scores_by_OSARS_group", lambda: draw_group_comparisons(scores, stats)),
        (
            "R1Q9_TCGA_continuous_OSARS_state_projection",
            lambda: draw_continuous_associations(scores, trends, correlations),
        ),
        ("R1Q9_multimethod_state_projection_effects_TCGA_ICGC", lambda: draw_effect_forest(stats)),
    ):
        manifest["figures"][name] = function()
    (FIGURES / "figure_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print("Rendered multi-method state-projection figures")


if __name__ == "__main__":
    main()
