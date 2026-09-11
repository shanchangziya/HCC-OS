#!/usr/bin/env python3
"""Nature-style, compact supplementary figure for NQO1 prioritization.

Inputs are the processed tables made by 01_prepare_nqo1_priority_evidence.R.
This plotting-only script performs no statistical testing or score calculation.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


HERE = Path(__file__).resolve()
ANALYSIS_DIR = HERE.parents[1]
DATA = ANALYSIS_DIR / "data" / "processed"
OUT = ANALYSIS_DIR / "figures" / "nature"
OUT.mkdir(parents=True, exist_ok=True)

# Reuse the manuscript-wide high/low semantic mapping, with neutral support.
LOW = "#0072B2"
HIGH = "#D55E00"
NEUTRAL_LIGHT = "#D8D8D8"
NEUTRAL_MID = "#8F8F8F"
NEUTRAL_DARK = "#4D4D4D"
BLACK = "#272727"

# Nature-style export rules: editable SVG, Arial-style sans serif, minimal axes.
mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans"],
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "font.size": 7.5,
        "axes.labelsize": 7.5,
        "axes.titlesize": 8.2,
        "xtick.labelsize": 6.8,
        "ytick.labelsize": 7.0,
        "axes.linewidth": 0.8,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "legend.frameon": False,
        "savefig.bbox": None,
        "savefig.pad_inches": 0.0,
    }
)


def style_axis(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(length=3.2, width=0.8, pad=2.3)


def add_panel_label(ax, label):
    ax.text(
        -0.075,
        1.08,
        label,
        transform=ax.transAxes,
        fontsize=9.4,
        fontweight="bold",
        ha="left",
        va="bottom",
        color=BLACK,
    )


def plot_candidate_panel(ax, candidates):
    """A: Compact dumbbell display of all eligible OSARS-ROS candidates."""
    d = candidates.sort_values("avg_log2FC", ascending=False).reset_index(drop=True)
    y = np.arange(len(d))[::-1]

    for yi, (_, row) in zip(y, d.iterrows()):
        low = row["pct.2"] * 100
        high = row["pct.1"] * 100
        is_nqo1 = row["gene"] == "NQO1"
        line_color = HIGH if is_nqo1 else NEUTRAL_MID
        ax.plot([low, high], [yi, yi], color=line_color, lw=1.5 if is_nqo1 else 1.0, solid_capstyle="round", zorder=1)
        ax.scatter(low, yi, s=25, color=LOW, edgecolor="white", linewidth=0.45, zorder=2)
        ax.scatter(high, yi, s=27 if is_nqo1 else 23, color=HIGH, edgecolor="white", linewidth=0.45, zorder=2)
        ax.text(
            116,
            yi,
            f"{row['avg_log2FC']:.2f}",
            va="center",
            ha="left",
            fontsize=7.0,
            color=HIGH if is_nqo1 else NEUTRAL_DARK,
            fontweight="bold" if is_nqo1 else "normal",
        )

    ax.set_yticks(y, d["gene"])
    # Highlight the selected candidate without a competing legend.
    for tick in ax.get_yticklabels():
        if tick.get_text() == "NQO1":
            tick.set_color(HIGH)
            tick.set_fontweight("bold")
    ax.set_xlim(0, 133)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_ylim(-0.45, len(d) - 0.45)
    ax.set_xlabel("Cells expressing gene (%)")
    ax.set_title("Frozen OSARS ∩ HALLMARK ROS (n = 3)", loc="left", pad=5, fontweight="bold")
    ax.text(0.01, 0.92, "● Low", color=LOW, transform=ax.transAxes, fontsize=6.6)
    ax.text(0.09, 0.92, "● High", color=HIGH, transform=ax.transAxes, fontsize=6.6)
    ax.text(0.87, 0.92, "Avg. log2FC", color=NEUTRAL_DARK, transform=ax.transAxes, fontsize=6.6, ha="center")
    style_axis(ax)
    add_panel_label(ax, "a")


def plot_marker_context(ax, markers):
    """B: Whole marker landscape establishes rank context without duplicate tables."""
    background = markers.loc[~markers["OSARS_ROS_overlap_gene"]]
    other = markers.loc[markers["OSARS_ROS_overlap_gene"] & ~markers["selected_for_follow_up"]]
    nq = markers.loc[markers["selected_for_follow_up"]].iloc[0]

    ax.scatter(
        background["avg_log2FC"],
        background["detect_gap"] * 100,
        s=5.0,
        color=NEUTRAL_LIGHT,
        alpha=0.72,
        linewidths=0,
        zorder=1,
    )
    ax.scatter(
        other["avg_log2FC"],
        other["detect_gap"] * 100,
        s=28,
        color=NEUTRAL_MID,
        edgecolor="white",
        linewidth=0.45,
        zorder=2,
    )
    ax.scatter(
        nq["avg_log2FC"],
        nq["detect_gap"] * 100,
        s=42,
        color=HIGH,
        edgecolor="white",
        linewidth=0.6,
        zorder=3,
    )
    ax.annotate(
        "NQO1\nrank 1/2,185",
        xy=(nq["avg_log2FC"], nq["detect_gap"] * 100),
        xytext=(nq["avg_log2FC"] - 0.40, nq["detect_gap"] * 100 - 9.5),
        fontsize=6.8,
        color=HIGH,
        fontweight="bold",
        ha="left",
        arrowprops={"arrowstyle": "-", "lw": 0.8, "color": HIGH},
    )
    ax.set_xlabel("Average log2FC")
    ax.set_ylabel("Δ detection (%)")
    ax.set_title("Discovery marker screen", loc="left", pad=5, fontweight="bold")
    style_axis(ax)
    add_panel_label(ax, "b")


def plot_protein_panel(ax, pairs, stats):
    """C: Independent pair-level protein abundance support."""
    adjacent = pairs["adjacent"].to_numpy()
    tumor = pairs["tumor"].to_numpy()
    for a, t in zip(adjacent, tumor):
        ax.plot([0, 1], [a, t], color=NEUTRAL_LIGHT, lw=0.42, alpha=0.70, zorder=1)
    rng = np.random.default_rng(20260906)
    for x, values, color, median in [
        (0, adjacent, LOW, stats["adjacent_median"]),
        (1, tumor, HIGH, stats["tumor_median"]),
    ]:
        ax.scatter(
            x + rng.normal(0, 0.017, size=len(values)),
            values,
            s=9.5,
            color=color,
            alpha=0.80,
            linewidths=0,
            zorder=2,
        )
        ax.hlines(median, x - 0.17, x + 0.17, color=BLACK, lw=1.25, zorder=3)

    ax.set_xlim(-0.38, 1.38)
    ax.set_xticks([0, 1], ["Adjacent", "Tumor"])
    ax.set_ylabel("NQO1 protein abundance\n(normalized log2 ratio)")
    ax.set_title("Paired PDC000198 HCC proteome", loc="left", pad=5, fontweight="bold")
    ax.text(
        0.02,
        0.94,
        f"n = {int(stats['n_pairs'])} pairs\nP = {stats['wilcoxon_p']:.2e}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=6.7,
        color=BLACK,
        bbox={"boxstyle": "round,pad=0.16", "facecolor": "white", "edgecolor": "none", "alpha": 0.88},
    )
    style_axis(ax)
    add_panel_label(ax, "c")


def save(fig, stem):
    # SVG is the editable primary export; PDF and high-resolution PNG are secondary.
    # Preserve the requested physical canvas dimensions exactly (170 × 150 mm
    # for the composite), rather than shrinking them with a tight crop.
    fig.savefig(OUT / f"{stem}.svg", bbox_inches=None, pad_inches=0)
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches=None, pad_inches=0)
    fig.savefig(OUT / f"{stem}.png", dpi=600, bbox_inches=None, pad_inches=0)


def make_composite(candidates, markers, pairs, stats):
    # Fixed at 170 mm wide × 150 mm high, safely below the requested 210 mm maximum.
    fig = plt.figure(figsize=(170 / 25.4, 150 / 25.4))
    gs = fig.add_gridspec(
        2,
        2,
        height_ratios=[0.67, 1.0],
        width_ratios=[0.88, 1.12],
        left=0.09,
        right=0.98,
        top=0.94,
        bottom=0.10,
        hspace=0.62,
        wspace=0.45,
    )
    ax_a = fig.add_subplot(gs[0, :])
    ax_b = fig.add_subplot(gs[1, 0])
    ax_c = fig.add_subplot(gs[1, 1])
    plot_candidate_panel(ax_a, candidates)
    plot_marker_context(ax_b, markers)
    plot_protein_panel(ax_c, pairs, stats)
    save(fig, "Supplementary_Figure_NQO1_prioritization_and_protein_Nature")
    plt.close(fig)


def make_standalone_panels(candidates, markers, pairs, stats):
    # Standalone source panels facilitate later manuscript assembly.
    fig, ax = plt.subplots(figsize=(170 / 25.4, 54 / 25.4))
    plot_candidate_panel(ax, candidates)
    save(fig, "Panel_a_NQO1_candidate_prioritization_Nature")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(82 / 25.4, 80 / 25.4))
    plot_marker_context(ax, markers)
    save(fig, "Panel_b_NQO1_marker_context_Nature")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(88 / 25.4, 80 / 25.4))
    plot_protein_panel(ax, pairs, stats)
    save(fig, "Panel_c_NQO1_paired_protein_Nature")
    plt.close(fig)


if __name__ == "__main__":
    candidate_data = pd.read_csv(DATA / "OSARS_ROS_overlap_candidate_evidence.csv")
    marker_data = pd.read_csv(DATA / "discovery_marker_rank_context.csv")
    pair_data = pd.read_csv(DATA / "PDC000198_NQO1_verified_pairs.csv")
    pair_stats = pd.read_csv(DATA / "PDC000198_NQO1_paired_statistics.csv").iloc[0]
    make_composite(candidate_data, marker_data, pair_data, pair_stats)
    make_standalone_panels(candidate_data, marker_data, pair_data, pair_stats)
    print(f"Nature-style figure exports written to {OUT}")
