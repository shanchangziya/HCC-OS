#!/usr/bin/env python3
"""Plot transparent NQO1 prioritization plus orthogonal protein evidence.

Statistics are read from the processed source tables created by
01_prepare_nqo1_priority_evidence.R; no quantities are recalculated here.
"""

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


HERE = Path(__file__).resolve()
ANALYSIS_DIR = HERE.parents[1]
DATA = ANALYSIS_DIR / "data" / "processed"
FIGURES = ANALYSIS_DIR / "figures"
FIGURES.mkdir(exist_ok=True)

NQO1 = "#D55E00"
LOW = "#0072B2"
MID = "#6E8FB3"
GREY = "#B8B8B8"
DARK = "#292929"
LIGHT = "#E7E7E7"

mpl.rcParams.update(
    {
        "font.family": "Arial",
        "font.size": 7.2,
        "axes.titlesize": 8.3,
        "axes.labelsize": 7.6,
        "xtick.labelsize": 6.8,
        "ytick.labelsize": 7.1,
        "axes.linewidth": 0.7,
        "xtick.major.width": 0.7,
        "ytick.major.width": 0.7,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.04,
    }
)


def clean_axis(ax):
    ax.spines[["top", "right"]].set_visible(False)


def plot_candidate_evidence(ax_effect, ax_detect, candidates):
    """Panel A: use a data-driven comparison only among 3 eligible genes."""
    cand = candidates.sort_values("avg_log2FC", ascending=False).reset_index(drop=True)
    y = np.arange(len(cand))
    colors = [NQO1 if g == "NQO1" else MID for g in cand["gene"]]

    # Effect-size lollipop
    ax_effect.hlines(y, 0, cand["avg_log2FC"], color=colors, lw=1.3, alpha=0.72)
    ax_effect.scatter(cand["avg_log2FC"], y, s=38, c=colors, edgecolors="white", linewidths=0.6, zorder=3)
    for yi, value in zip(y, cand["avg_log2FC"]):
        ax_effect.text(value + 0.035, yi, f"{value:.2f}", va="center", fontsize=6.5)
    ax_effect.set_yticks(y, cand["gene"])
    ax_effect.invert_yaxis()
    ax_effect.set_xlim(0, max(cand["avg_log2FC"]) * 1.28)
    ax_effect.set_xlabel("Average log2FC")
    ax_effect.set_title("Marker effect size", loc="left", pad=5, fontweight="bold")
    ax_effect.grid(axis="x", color=LIGHT, lw=0.65, zorder=0)
    clean_axis(ax_effect)

    # Detection-rate dumbbells make specificity visible without relying solely
    # on a post-clustering p value.
    for yi, row in cand.iterrows():
        color = NQO1 if row["gene"] == "NQO1" else MID
        ax_detect.plot([row["pct.2"] * 100, row["pct.1"] * 100], [yi, yi], color=color, lw=1.4, alpha=0.8, zorder=1)
        ax_detect.scatter(row["pct.2"] * 100, yi, s=27, color=LOW, edgecolor="white", linewidth=0.45, zorder=2)
        ax_detect.scatter(row["pct.1"] * 100, yi, s=31, color=color, edgecolor="white", linewidth=0.45, zorder=2)
        ax_detect.text(102.5, yi, f"Δ{row['detect_gap'] * 100:.1f}%", va="center", fontsize=6.25, color=DARK)
    ax_detect.set_yticks(y, [""] * len(cand))
    ax_detect.invert_yaxis()
    ax_detect.set_xlim(0, 125)
    ax_detect.set_xlabel("Cells expressing gene (%)")
    ax_detect.set_title("Detection specificity", loc="left", pad=5, fontweight="bold")
    ax_detect.grid(axis="x", color=LIGHT, lw=0.65, zorder=0)
    clean_axis(ax_detect)
    return ax_effect, ax_detect


def plot_marker_context(ax, markers):
    """Panel B: whole discovery marker-space context (descriptive only)."""
    non_candidates = markers.loc[~markers["OSARS_ROS_overlap_gene"]]
    other_candidates = markers.loc[markers["OSARS_ROS_overlap_gene"] & ~markers["selected_for_follow_up"]]
    selected = markers.loc[markers["selected_for_follow_up"]].iloc[0]

    ax.scatter(
        non_candidates["avg_log2FC"],
        non_candidates["detect_gap"] * 100,
        s=4,
        color=GREY,
        alpha=0.42,
        linewidths=0,
        zorder=1,
    )
    ax.scatter(
        other_candidates["avg_log2FC"],
        other_candidates["detect_gap"] * 100,
        s=30,
        color=MID,
        edgecolor="white",
        linewidth=0.45,
        zorder=2,
        label="Other OSARS–ROS overlap genes",
    )
    ax.scatter(
        [selected["avg_log2FC"]],
        [selected["detect_gap"] * 100],
        s=52,
        color=NQO1,
        edgecolor="white",
        linewidth=0.7,
        zorder=3,
        label="NQO1",
    )
    ax.annotate(
        "NQO1\nrank 1/2,185",
        xy=(selected["avg_log2FC"], selected["detect_gap"] * 100),
        xytext=(selected["avg_log2FC"] - 0.43, selected["detect_gap"] * 100 - 10),
        fontsize=6.7,
        fontweight="bold",
        color=NQO1,
        ha="left",
        arrowprops={"arrowstyle": "-", "color": NQO1, "lw": 0.8},
    )
    ax.set_xlabel("Average log2FC, historical OS-high vs OS-low")
    ax.set_ylabel("Detection-rate difference (percentage points)")
    ax.set_title("B  Discovery marker-screen context", loc="left", fontweight="bold", pad=5)
    ax.grid(color=LIGHT, lw=0.65, zorder=0)
    clean_axis(ax)
    ax.legend(frameon=False, fontsize=5.9, loc="lower right", handletextpad=0.3, borderpad=0.1)


def plot_pdc_pairs(ax, pairs, stats):
    """Panel C: orthogonal paired protein evidence requested by Reviewer 2."""
    adjacent = pairs["adjacent"].to_numpy()
    tumor = pairs["tumor"].to_numpy()
    for a, t in zip(adjacent, tumor):
        ax.plot([0, 1], [a, t], color="#BDBDBD", lw=0.45, alpha=0.36, zorder=1)
    rng = np.random.default_rng(20260906)
    for pos, values, color in [(0, adjacent, LOW), (1, tumor, NQO1)]:
        jitter = rng.normal(0, 0.018, len(values))
        ax.scatter(np.full(len(values), pos) + jitter, values, s=9, color=color, alpha=0.78, linewidths=0, zorder=2)
        median = np.median(values)
        ax.hlines(median, pos - 0.17, pos + 0.17, color=DARK, lw=1.35, zorder=3)
    ax.set_xlim(-0.42, 1.42)
    ax.set_xticks([0, 1], ["Adjacent", "Tumor"])
    ax.set_ylabel("NQO1 abundance\n(normalized unshared log2 protein ratio)")
    ax.set_title("C  Independent paired protein evidence", loc="left", fontweight="bold", pad=5)
    ax.text(
        0.02,
        0.97,
        f"PDC000198, n={int(stats['n_pairs'])} verified pairs; paired Wilcoxon P={stats['wilcoxon_p']:.2e}",
        transform=ax.transAxes,
        fontsize=6.1,
        color="#4D4D4D",
        va="top",
        bbox={"boxstyle": "round,pad=0.15", "facecolor": "white", "edgecolor": "none", "alpha": 0.78},
    )
    ax.text(
        0.5,
        0.05,
        f"Median paired shift = {stats['hodges_lehmann_shift']:.2f}\n95% CI {stats['hl_lower']:.2f} to {stats['hl_upper']:.2f}",
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=6.5,
        color=DARK,
        bbox={"boxstyle": "round,pad=0.22", "facecolor": "white", "edgecolor": "#D0D0D0", "linewidth": 0.55, "alpha": 0.95},
    )
    clean_axis(ax)


def make_figure(candidates, markers, pairs, stats):
    # 170 mm wide main-figure candidate; source plots remain vector in PDF.
    fig = plt.figure(figsize=(170 / 25.4, 150 / 25.4))
    grid = fig.add_gridspec(
        2,
        2,
        height_ratios=[0.78, 1.0],
        width_ratios=[1.0, 1.0],
        left=0.095,
        right=0.985,
        top=0.83,
        bottom=0.16,
        hspace=0.72,
        wspace=0.38,
    )
    subgrid_a = grid[0, :].subgridspec(1, 2, width_ratios=[0.85, 1.15], wspace=0.45)
    ax_a1 = fig.add_subplot(subgrid_a[0, 0])
    ax_a2 = fig.add_subplot(subgrid_a[0, 1])
    plot_candidate_evidence(ax_a1, ax_a2, candidates)
    ax_b = fig.add_subplot(grid[1, 0])
    ax_c = fig.add_subplot(grid[1, 1])
    plot_marker_context(ax_b, markers)
    plot_pdc_pairs(ax_c, pairs, stats)

    fig.suptitle("NQO1 prioritization and orthogonal protein support", x=0.095, ha="left", fontsize=11.6, fontweight="bold")
    fig.text(
        0.095,
        0.875,
        "A  Post-model prioritization: frozen OSARS ∩ HALLMARK ROS (n=3)",
        fontsize=9.2,
        fontweight="bold",
        va="bottom",
    )
    fig.text(
        0.095,
        0.055,
        "Discovery marker statistics are descriptive because the groups were defined within the discovery dataset; they do not constitute independent single-cell-state validation.",
        ha="left",
        va="bottom",
        fontsize=5.8,
        color="#454545",
    )
    fig.text(
        0.095,
        0.032,
        "NQO1 was not selected by GBM feature importance (rank 56/128) and is not a frozen discovery MP4 gene. PDC protein evidence supports abundance differences, not prognostic or causal function.",
        ha="left",
        va="bottom",
        fontsize=5.8,
        color="#454545",
    )
    fig.savefig(FIGURES / "Supplementary_NQO1_prioritization_and_protein.pdf")
    fig.savefig(FIGURES / "Supplementary_NQO1_prioritization_and_protein.png", dpi=600)
    plt.close(fig)


if __name__ == "__main__":
    candidates = pd.read_csv(DATA / "OSARS_ROS_overlap_candidate_evidence.csv")
    markers = pd.read_csv(DATA / "discovery_marker_rank_context.csv")
    pairs = pd.read_csv(DATA / "PDC000198_NQO1_verified_pairs.csv")
    stats = pd.read_csv(DATA / "PDC000198_NQO1_paired_statistics.csv").iloc[0]
    make_figure(candidates, markers, pairs, stats)
    print(f"Wrote figure to {FIGURES}")
