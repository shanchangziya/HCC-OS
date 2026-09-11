#!/usr/bin/env python3
"""Render the accepted all-cell UMAP as a compact standalone Panel A."""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
OUT_ROOT = SCRIPT_DIR.parent
FIG_DIR = OUT_ROOT / "figures" / "panel_review"
FIG_DIR.mkdir(parents=True, exist_ok=True)

SOURCE = OUT_ROOT / "GSE149614_all_cells_UMAP_metadata.csv"
MM_PER_INCH = 25.4
WIDTH_MM = 98
HEIGHT_MM = 70

PANEL_FONT = "Times New Roman"
TEXT_FONT = "Arial"

CELL_COLORS = {
    "Hepatocyte": "#D55E00",
    "T/NK": "#0072B2",
    "Myeloid": "#009E73",
    "B": "#CC79A7",
    "Endothelial": "#56B4E9",
    "Fibroblast": "#7A7A7A",
}
DRAW_ORDER = ["Fibroblast", "Endothelial", "B", "T/NK", "Myeloid", "Hepatocyte"]
LEGEND_ORDER = ["Hepatocyte", "T/NK", "Myeloid", "B", "Endothelial", "Fibroblast"]


def main() -> None:
    data = pd.read_csv(SOURCE)
    if len(data) != 71915:
        raise ValueError(f"Expected 71,915 cells, found {len(data):,}.")
    if set(data["major_cell_type"]) != set(CELL_COLORS):
        raise ValueError("Unexpected major-cell-type labels.")
    if data[["UMAP1", "UMAP2"]].isna().any().any():
        raise ValueError("Missing UMAP coordinates.")

    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": [TEXT_FONT, "Helvetica", "DejaVu Sans"],
            "font.size": 6.2,
            "legend.fontsize": 6.0,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )

    fig = plt.figure(
        figsize=(WIDTH_MM / MM_PER_INCH, HEIGHT_MM / MM_PER_INCH), facecolor="white"
    )
    # The UMAP remains nearly square; the remaining narrow strip holds the
    # legend.  No unused half-page canvas is retained.
    ax = fig.add_axes([0.105, 0.075, 0.615, 0.845])
    for cell_type in DRAW_ORDER:
        block = data.loc[data["major_cell_type"].eq(cell_type)]
        ax.scatter(
            block["UMAP1"],
            block["UMAP2"],
            s=0.30,
            color=CELL_COLORS[cell_type],
            linewidths=0,
            rasterized=True,
        )

    ax.set_aspect("equal")
    ax.set_anchor("C")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.margins(0.025)
    for spine in ax.spines.values():
        spine.set_visible(False)

    ax.text(
        0,
        1.015,
        "Major cell types",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=7.2,
        color="#222222",
    )
    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            label=cell_type,
            markerfacecolor=CELL_COLORS[cell_type],
            markeredgewidth=0,
            markersize=3.2,
        )
        for cell_type in LEGEND_ORDER
    ]
    ax.legend(
        handles=handles,
        loc="center left",
        bbox_to_anchor=(1.005, 0.50),
        borderaxespad=0,
        frameon=False,
        handletextpad=0.35,
        labelspacing=0.46,
    )

    fig.text(
        0.015,
        0.978,
        "A",
        ha="left",
        va="top",
        fontsize=10.5,
        fontfamily=PANEL_FONT,
        fontweight="bold",
        color="#111111",
    )

    stem = FIG_DIR / "Panel_A_all_cell_UMAP_review_v2"
    fig.savefig(stem.with_suffix(".pdf"), dpi=600, bbox_inches=None)
    fig.savefig(stem.with_suffix(".svg"), dpi=600, bbox_inches=None)
    fig.savefig(stem.with_suffix(".png"), dpi=600, bbox_inches=None)
    plt.close(fig)
    print(f"Saved {stem}.pdf/.svg/.png")
    print(f"Canvas: {WIDTH_MM} x {HEIGHT_MM} mm; cells: {len(data):,}")


if __name__ == "__main__":
    main()
