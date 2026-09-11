#!/usr/bin/env python3
"""Draw the review candidate for Panel C: GeneNMF similarity heatmap."""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
OUT_ROOT = SCRIPT_DIR.parent
FIG_DIR = OUT_ROOT / "figures" / "panel_review"
FIG_DIR.mkdir(parents=True, exist_ok=True)

MM_PER_INCH = 25.4
WIDTH_MM = 85
HEIGHT_MM = 82

PANEL_FONT = "Times New Roman"
TEXT_FONT = "Arial"

MP_COLORS = {
    "MP1": "#3F7F93",
    "MP2": "#D6973C",
    "MP3": "#6E63A8",
    "MP4": "#C95F6F",
    "MP5": "#5D9368",
    "MP6": "#9A704B",
    "MP7": "#5986C4",
    "MP8": "#9A6FAE",
}


def contiguous_blocks(labels: np.ndarray) -> list[tuple[int, int, str]]:
    change = np.flatnonzero(labels[1:] != labels[:-1]) + 1
    starts = np.r_[0, change]
    ends = np.r_[change, len(labels)]
    return [(int(a), int(b), str(labels[a])) for a, b in zip(starts, ends)]


def main() -> None:
    sweep = pd.read_csv(OUT_ROOT / "GSE149614_GeneNMF_nMP_sweep.csv")
    selected_rows = sweep.loc[sweep["selected"].astype(bool)]
    if len(selected_rows) != 1 or int(selected_rows.iloc[0]["requested_nMP"]) != 8:
        raise ValueError("Panel C expects the audited eight-metaprogram solution.")

    similarity = pd.read_csv(
        OUT_ROOT / "GSE149614_GeneNMF_program_similarity.csv.gz", index_col=0
    )
    membership = pd.read_csv(
        OUT_ROOT / "GSE149614_GeneNMF_individual_program_membership.csv"
    ).set_index("individual_program")
    order = (
        pd.read_csv(OUT_ROOT / "GSE149614_GeneNMF_heatmap_order.csv")
        .sort_values("order_index")["individual_program"]
        .astype(str)
        .tolist()
    )

    if similarity.shape != (390, 390) or set(order) != set(similarity.index):
        raise ValueError("Unexpected GeneNMF similarity matrix or tree order.")
    matrix = similarity.loc[order, order].to_numpy(dtype=float)
    labels = membership.loc[order, "metaprogram"].astype(str).to_numpy()
    blocks = contiguous_blocks(labels)
    if len(blocks) != 8 or set(labels) != set(MP_COLORS):
        raise ValueError("Metaprogram blocks are incomplete or non-contiguous.")

    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": [TEXT_FONT, "DejaVu Sans"],
            "font.size": 6.0,
            "axes.linewidth": 0.45,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )

    cmap = LinearSegmentedColormap.from_list(
        "genenmf_similarity",
        ["#FAFAF8", "#D9E5E8", "#86B3C1", "#2F7389", "#143E57"],
        N=256,
    )
    norm = mpl.colors.Normalize(vmin=0.05, vmax=0.80, clip=True)

    fig = plt.figure(
        figsize=(WIDTH_MM / MM_PER_INCH, HEIGHT_MM / MM_PER_INCH),
        facecolor="white",
    )

    # Square heatmap with two narrow MP annotation strips.
    ax = fig.add_axes([0.145, 0.175, 0.705, 0.730])
    ax_top = fig.add_axes([0.145, 0.911, 0.705, 0.018])
    ax_left = fig.add_axes([0.121, 0.175, 0.018, 0.730])
    cax = fig.add_axes([0.548, 0.080, 0.285, 0.022])

    image = ax.imshow(
        matrix,
        cmap=cmap,
        norm=norm,
        origin="upper",
        interpolation="nearest",
        aspect="auto",
        rasterized=True,
    )
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    mp_to_index = {mp: idx for idx, mp in enumerate(MP_COLORS)}
    strip_values = np.array([mp_to_index[x] for x in labels], dtype=int)
    strip_cmap = ListedColormap([MP_COLORS[mp] for mp in MP_COLORS])
    ax_top.imshow(strip_values[np.newaxis, :], cmap=strip_cmap, vmin=0, vmax=7,
                  interpolation="nearest", aspect="auto")
    ax_left.imshow(strip_values[:, np.newaxis], cmap=strip_cmap, vmin=0, vmax=7,
                   interpolation="nearest", aspect="auto")
    for strip_ax in (ax_top, ax_left):
        strip_ax.set_xticks([])
        strip_ax.set_yticks([])
        for spine in strip_ax.spines.values():
            spine.set_visible(False)

    # White separators delineate the consensus blocks while keeping the plot quiet.
    for start, end, mp in blocks[:-1]:
        boundary = end - 0.5
        ax.axvline(boundary, color="white", lw=0.38, alpha=0.95)
        ax.axhline(boundary, color="white", lw=0.38, alpha=0.95)
        ax_top.axvline(boundary, color="white", lw=0.45)
        ax_left.axhline(boundary, color="white", lw=0.45)

    # Direct MP labels follow the dendrogram order; no separate legend is needed.
    for start, end, mp in blocks:
        center = (start + end - 1) / 2
        x_fig = 0.145 + 0.705 * ((center + 0.5) / len(labels))
        fig.text(
            x_fig,
            0.945,
            mp,
            ha="center",
            va="bottom",
            rotation=90,
            fontsize=5.4,
            color=MP_COLORS[mp],
            fontfamily=TEXT_FONT,
            fontweight="bold",
        )

    cb = fig.colorbar(image, cax=cax, orientation="horizontal", ticks=[0.1, 0.4, 0.8])
    cb.outline.set_linewidth(0.35)
    cb.ax.tick_params(length=1.8, width=0.4, pad=1.4, labelsize=5.2)
    cb.ax.set_title("Cosine similarity", fontsize=5.5, pad=2.0, loc="left")

    fig.text(
        0.020,
        0.975,
        "C",
        ha="left",
        va="top",
        fontsize=10.5,
        fontfamily=PANEL_FONT,
        fontweight="bold",
        color="#111111",
    )

    stem = FIG_DIR / "Panel_C_GeneNMF_similarity_heatmap_review_v1"
    fig.savefig(stem.with_suffix(".pdf"), dpi=600, bbox_inches=None)
    fig.savefig(stem.with_suffix(".svg"), dpi=600, bbox_inches=None)
    fig.savefig(stem.with_suffix(".png"), dpi=600, bbox_inches=None)
    plt.close(fig)

    print(f"Saved {stem}.pdf/.svg/.png")
    print("Block order:", ", ".join(f"{mp} ({end-start})" for start, end, mp in blocks))


if __name__ == "__main__":
    main()
