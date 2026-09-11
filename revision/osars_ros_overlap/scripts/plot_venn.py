#!/usr/bin/env python3
"""Draw the frozen OSARS versus HALLMARK ROS overlap as a compact vector Venn diagram."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Circle


MODULE = Path(__file__).resolve().parents[1]
REVISION = MODULE.parent
OSARS_TABLE = REVISION / "tables" / "Table2_frozen_OSARS_128_genes.csv"
ROS_FILE = (
    REVISION
    / "model_validation"
    / "references"
    / "HALLMARK_ROS_msigdb_v7.0_genes.txt"
)
FIGURE_DIR = MODULE / "figures"
DATA_DIR = MODULE / "data"


def read_sets() -> tuple[set[str], set[str]]:
    with OSARS_TABLE.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    osars = {row["gene"].strip() for row in rows if row["gene"].strip()}
    ros = {
        line.strip()
        for line in ROS_FILE.read_text(encoding="utf-8").splitlines()
        if line.strip()
    }
    overlap = osars & ros
    assert len(osars) == 128, f"Expected 128 frozen OSARS genes, found {len(osars)}"
    assert len(ros) == 49, f"Expected 49 HALLMARK ROS genes, found {len(ros)}"
    assert overlap == {"NQO1", "PRDX1", "TXNRD1"}, overlap
    return osars, ros


def write_source_data(osars: set[str], ros: set[str]) -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    with (DATA_DIR / "OSARS_HALLMARK_ROS_overlap_source.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.writer(handle)
        writer.writerow(["gene", "in_OSARS", "in_HALLMARK_ROS", "venn_region"])
        for gene in sorted(osars | ros):
            in_osars = gene in osars
            in_ros = gene in ros
            region = "overlap" if in_osars and in_ros else "OSARS_only" if in_osars else "HALLMARK_ROS_only"
            writer.writerow([gene, in_osars, in_ros, region])


def draw_venn(osars: set[str], ros: set[str]) -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "font.size": 8,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.spines.right": False,
            "axes.spines.top": False,
            "axes.linewidth": 1.0,
            "legend.frameon": False,
        }
    )

    overlap = sorted(osars & ros)
    left_only = len(osars - ros)
    right_only = len(ros - osars)

    # Compact single-column panel: 85 mm wide, with direct labels and no title.
    fig, ax = plt.subplots(figsize=(85 / 25.4, 62 / 25.4), facecolor="white")
    ax.set_position([0.01, 0.01, 0.98, 0.98])
    ax.set_xlim(0.25, 9.75)
    ax.set_ylim(0.35, 6.55)
    ax.set_aspect("equal")
    ax.axis("off")

    blue = "#3775BA"
    red = "#B64342"
    text_dark = "#272727"

    ax.add_patch(
        Circle(
            (3.85, 3.18),
            2.28,
            facecolor=blue,
            edgecolor=blue,
            linewidth=1.25,
            alpha=0.25,
        )
    )
    ax.add_patch(
        Circle(
            (6.15, 3.18),
            2.28,
            facecolor=red,
            edgecolor=red,
            linewidth=1.25,
            alpha=0.23,
        )
    )

    # Set names are direct labels; no detached legend or figure title is needed.
    ax.text(2.70, 6.10, "OSARS", ha="center", va="center", fontsize=9, fontweight="bold", color=blue)
    ax.text(2.70, 5.72, f"n = {len(osars)}", ha="center", va="center", fontsize=7.3, color=text_dark)
    ax.text(7.30, 6.10, "HALLMARK ROS", ha="center", va="center", fontsize=9, fontweight="bold", color=red)
    ax.text(7.30, 5.72, f"n = {len(ros)}", ha="center", va="center", fontsize=7.3, color=text_dark)

    ax.text(2.76, 3.18, str(left_only), ha="center", va="center", fontsize=15, fontweight="bold", color=text_dark)
    ax.text(7.24, 3.18, str(right_only), ha="center", va="center", fontsize=15, fontweight="bold", color=text_dark)

    ax.text(5.00, 3.85, str(len(overlap)), ha="center", va="center", fontsize=13, fontweight="bold", color=text_dark)
    ax.text(
        5.00,
        2.75,
        "\n".join(overlap),
        ha="center",
        va="center",
        fontsize=7.2,
        fontstyle="italic",
        linespacing=1.28,
        color=text_dark,
    )

    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    stem = FIGURE_DIR / "OSARS_HALLMARK_ROS_Venn"
    save_args = {"bbox_inches": "tight", "pad_inches": 0.02, "facecolor": "white"}
    fig.savefig(stem.with_suffix(".svg"), **save_args)
    fig.savefig(stem.with_suffix(".pdf"), metadata={"Title": ""}, **save_args)
    fig.savefig(stem.with_suffix(".png"), dpi=600, **save_args)
    plt.close(fig)


def main() -> None:
    osars, ros = read_sets()
    write_source_data(osars, ros)
    draw_venn(osars, ros)


if __name__ == "__main__":
    main()
