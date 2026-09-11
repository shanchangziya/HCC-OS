#!/usr/bin/env python3
"""Create two standalone, minimal diagnostic plots from frozen CSV outputs."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


BLUE = "#0072B2"
ORANGE = "#D55E00"
DARK = "#252525"
MID_GREY = "#737373"
LIGHT_GREY = "#D9D9D9"
MM_TO_INCH = 1 / 25.4


def configure_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 7,
            "axes.labelsize": 7,
            "axes.titlesize": 8,
            "xtick.labelsize": 6.5,
            "ytick.labelsize": 6.5,
            "legend.fontsize": 6.5,
            "axes.linewidth": 0.6,
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
            "xtick.major.size": 2.5,
            "ytick.major.size": 2.5,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
        }
    )


def extend_step(data: pd.DataFrame, end: float = 5.0) -> pd.DataFrame:
    data = data.sort_values("time_years").copy()
    data = data[data["time_years"] <= end]
    if data.empty:
        raise ValueError("No KM coordinates within plotting range")
    if data["time_years"].iloc[-1] < end:
        last = data.iloc[-1].copy()
        last["time_years"] = end
        last["time_days"] = end * 365.25
        data = pd.concat([data, last.to_frame().T], ignore_index=True)
    return data


def save_both(fig: mpl.figure.Figure, stem: Path) -> None:
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", pad_inches=0.03)
    fig.savefig(stem.with_suffix(".png"), dpi=600, bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)


def plot_test_km(output_dir: Path) -> None:
    curves = pd.read_csv(output_dir / "pathology_KM_coordinates.csv")
    risk = pd.read_csv(output_dir / "pathology_KM_risk_table.csv")
    stats = pd.read_csv(output_dir / "pathology_KM_statistics.csv")
    curves = curves[curves["split"] == "test"].copy()
    risk = risk[risk["split"] == "test"].copy()
    stat = stats.loc[stats["split"] == "test"].iloc[0]

    fig = plt.figure(figsize=(88 * MM_TO_INCH, 82 * MM_TO_INCH))
    grid = fig.add_gridspec(2, 1, height_ratios=[4.2, 1.05], hspace=0.18)
    ax = fig.add_subplot(grid[0])
    ax_risk = fig.add_subplot(grid[1], sharex=ax)

    group_info = {
        "Low": (BLUE, int(stat["n_low"])),
        "High": (ORANGE, int(stat["n_high"])),
    }
    for group, (color, n_group) in group_info.items():
        d = extend_step(curves[curves["risk_group"] == group])
        x = d["time_years"].astype(float).to_numpy()
        y = d["survival"].astype(float).to_numpy()
        lo = d["CI_lower"].astype(float).to_numpy()
        hi = d["CI_upper"].astype(float).to_numpy()
        ax.fill_between(x, lo, hi, step="post", color=color, alpha=0.10, linewidth=0)
        ax.step(x, y, where="post", color=color, linewidth=1.15, label=f"{group} (n={n_group})")
        censored = d[(d["n_censor"].astype(float) > 0) & (d["time_years"].astype(float) <= 5)]
        ax.scatter(
            censored["time_years"], censored["survival"], marker="|", s=12,
            linewidths=0.65, color=color, zorder=3,
        )

    ax.set_xlim(0, 5)
    ax.set_ylim(0, 1.02)
    ax.set_ylabel("Overall survival")
    ax.set_title("Strict hold-out test", loc="left", fontweight="normal", pad=3)
    ax.legend(frameon=False, loc="lower left", handlelength=1.8, borderaxespad=0.3)
    ax.text(
        0.98, 0.96,
        f"HR {stat['HR_high_vs_low']:.2f} ({stat['HR_CI_lower']:.2f}–{stat['HR_CI_upper']:.2f})\n"
        f"Log-rank P = {stat['logrank_P']:.3f}",
        transform=ax.transAxes, ha="right", va="top", color=DARK, linespacing=1.25,
    )
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(axis="x", labelbottom=False)
    ax.grid(False)

    times = np.arange(0, 6)
    ax_risk.set_xlim(0, 5)
    ax_risk.set_ylim(-0.7, 1.3)
    ax_risk.set_xticks(times)
    ax_risk.set_xlabel("Time (years)")
    ax_risk.set_yticks([1, 0], labels=["Low", "High"])
    ax_risk.tick_params(axis="y", length=0, pad=9)
    for label, color in zip(ax_risk.get_yticklabels(), [BLUE, ORANGE]):
        label.set_color(color)
    for row_y, group in zip([1, 0], ["Low", "High"]):
        d = risk[risk["risk_group"] == group].copy()
        for year in times:
            nearest = d.iloc[(d["time_years"] - year).abs().argmin()]
            if year == 0:
                ax_risk.text(0.03, row_y, f"{int(nearest['n_risk'])}", ha="left", va="center")
            elif year == 5:
                ax_risk.text(4.97, row_y, f"{int(nearest['n_risk'])}", ha="right", va="center")
            else:
                ax_risk.text(year, row_y, f"{int(nearest['n_risk'])}", ha="center", va="center")
    ax_risk.spines[["left", "right", "top"]].set_visible(False)
    ax_risk.spines["bottom"].set_linewidth(0.6)
    ax_risk.grid(False)

    fig.subplots_adjust(left=0.17, right=0.98, top=0.94, bottom=0.14)
    save_both(fig, output_dir / "pathology_strict_test_KM")


def plot_cindex_forest(output_dir: Path) -> None:
    performance = pd.read_csv(output_dir / "pathology_performance_summary.csv")
    performance = performance.set_index("split").loc[["train", "test"]].reset_index()
    labels = ["Training", "Strict test"]
    colors = [MID_GREY, ORANGE]
    markers = ["o", "s"]
    y = np.array([1, 0])

    fig, ax = plt.subplots(figsize=(82 * MM_TO_INCH, 45 * MM_TO_INCH))
    ax.axvline(0.5, color=LIGHT_GREY, linewidth=0.8, linestyle=(0, (2, 2)), zorder=0)
    for idx, row in performance.iterrows():
        estimate = float(row["C_index"])
        lower = float(row["C_index_CI_lower"])
        upper = float(row["C_index_CI_upper"])
        ax.errorbar(
            estimate, y[idx], xerr=[[estimate - lower], [upper - estimate]],
            fmt=markers[idx], color=colors[idx], ecolor=colors[idx],
            markersize=4.2, markeredgewidth=0, elinewidth=1.0, capsize=2.0,
            zorder=2,
        )
        ax.text(estimate, y[idx] + 0.18, f"{estimate:.2f} ({lower:.2f}–{upper:.2f})",
                ha="center", va="bottom", color=DARK)

    ax.set_yticks(y, labels=labels)
    ax.set_xlim(0.44, 0.86)
    ax.set_ylim(-0.62, 1.55)
    ax.set_xticks([0.5, 0.6, 0.7, 0.8])
    ax.set_xlabel("Harrell C-index")
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.grid(False)
    fig.subplots_adjust(left=0.25, right=0.99, top=0.92, bottom=0.28)
    save_both(fig, output_dir / "pathology_Cindex_train_test_forest")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    output_dir = args.output_dir.resolve()
    required = [
        "pathology_KM_coordinates.csv", "pathology_KM_risk_table.csv",
        "pathology_KM_statistics.csv", "pathology_performance_summary.csv",
    ]
    missing = [name for name in required if not (output_dir / name).is_file()]
    if missing:
        raise FileNotFoundError(f"Missing frozen result files: {missing}")
    configure_style()
    plot_test_km(output_dir)
    plot_cindex_forest(output_dir)


if __name__ == "__main__":
    main()
