#!/usr/bin/env python3
"""Create revised HCC-2T Fig1A/B panels requested by the user."""

from __future__ import annotations

import gzip
import os
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge, Patch
from matplotlib.collections import PatchCollection
import numpy as np
import pandas as pd


BASE_DIR = Path(
    os.environ.get("HCC_OS_SPATIAL_WORKDIR", Path(__file__).resolve().parents[1])
).expanduser().resolve()
SAMPLE = "HCC-2T"
OUT_DIR = BASE_DIR / "results" / "hcc2t_main_figures"

MAJOR_GROUPS = {
    "Tumor_OS-high": "Tumor OS-high",
    "Tumor_OS-low": "Tumor OS-low",
    "Mac_TREM2_TAM": "Mac TREM2 TAM",
    "Mac_SLC40A1_TAM": "Mac SLC40A1 TAM",
    "Mac_HSP_TAM": "Mac HSP TAM",
    "CAF_myCAF": "CAF",
    "CAF_apCAF": "CAF",
    "CAF_MMP_CAF": "CAF",
    "Endo_LSEC": "Endothelial",
    "Endo_VEC": "Endothelial",
    "DC_cDC1": "DC",
    "DC_cDC2": "DC",
    "DC_moDC": "DC",
    "Neu_PMN_MDSC": "Neutrophil/MDSC",
    "B_Activated_B": "B/Plasma",
    "B_IgG_Plasma": "B/Plasma",
    "B_Naive_B": "B/Plasma",
    "NK_CD56dim_NK": "NK",
    "CD4_Treg": "CD4 T",
    "CD4_Tfh": "CD4 T",
    "CD4_Cyto_Th1_Cyto": "CD4 T",
    "CD4_Cyto_Memory_Cyto": "CD4 T",
    "CD4_Cyto_Exhausted_Cyto": "CD4 T",
    "CD4_Cyto_Effector_Cyto": "CD4 T",
    "CD4_Cyto_CTL_like": "CD4 T",
    "CD8_Trm_Activated_Trm": "CD8 T",
    "CD8_Trm_Cytotoxic_Trm": "CD8 T",
    "CD8_Trm_Precursor_Trm": "CD8 T",
}

GROUP_ORDER = [
    "Tumor OS-high",
    "Tumor OS-low",
    "Mac TREM2 TAM",
    "Mac SLC40A1 TAM",
    "Mac HSP TAM",
    "CAF",
    "Endothelial",
    "CD4 T",
    "CD8 T",
    "B/Plasma",
    "DC",
    "NK",
    "Neutrophil/MDSC",
]

GROUP_COLORS = {
    "Tumor OS-high": "#D55E00",
    "Tumor OS-low": "#E69F00",
    "Mac TREM2 TAM": "#0072B2",
    "Mac SLC40A1 TAM": "#56B4E9",
    "Mac HSP TAM": "#009E73",
    "CAF": "#CC79A7",
    "Endothelial": "#7A7A7A",
    "CD4 T": "#6A3D9A",
    "CD8 T": "#CAB2D6",
    "B/Plasma": "#A6CEE3",
    "DC": "#B15928",
    "NK": "#33A02C",
    "Neutrophil/MDSC": "#8C6D31",
}

mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans"],
        "font.size": 8,
        "axes.labelsize": 8,
        "axes.titlesize": 9,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 6.2,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.linewidth": 0.6,
    }
)


def read_abundance() -> pd.DataFrame:
    path = BASE_DIR / "results" / "cell2location_spatial" / SAMPLE / "cell_abundance_mean.csv"
    return pd.read_csv(path, index_col=0)


def read_coords(index: pd.Index) -> pd.DataFrame:
    path = BASE_DIR / "inputs" / "spatial" / SAMPLE / "spatial_coords.tsv.gz"
    with gzip.open(path, "rt") as handle:
        coords = pd.read_csv(handle, sep="\t").set_index("barcode")
    return coords.loc[index]


def save(fig: mpl.figure.Figure, name: str) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_DIR / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(OUT_DIR / f"{name}.png", dpi=300, bbox_inches="tight")


def panel_label(ax: mpl.axes.Axes, label: str) -> None:
    ax.text(-0.05, 1.04, label, transform=ax.transAxes, fontsize=12, fontweight="bold", va="bottom", ha="right")


def group_abundance(abund: pd.DataFrame) -> pd.DataFrame:
    missing = sorted(set(MAJOR_GROUPS) - set(abund.columns))
    if missing:
        raise RuntimeError(f"Missing abundance columns for grouping: {missing}")
    grouped = pd.DataFrame(index=abund.index)
    for group in GROUP_ORDER:
        cols = [state for state, grp in MAJOR_GROUPS.items() if grp == group]
        grouped[group] = abund[cols].sum(axis=1)
    grouped.to_csv(OUT_DIR / "HCC-2T_spot_major_group_abundance_for_pie.csv")
    return grouped


def make_spatial_pie(abund: pd.DataFrame, coords: pd.DataFrame) -> None:
    grouped = group_abundance(abund)
    fractions = grouped.div(grouped.sum(axis=1), axis=0).fillna(0)

    fig, ax = plt.subplots(figsize=(6.4, 5.6), constrained_layout=True)
    radius = 3.8
    for group in GROUP_ORDER:
        patches = []
        colors = []
        current_start = fractions[GROUP_ORDER[: GROUP_ORDER.index(group)]].sum(axis=1).to_numpy() * 360.0
        current_end = current_start + fractions[group].to_numpy() * 360.0
        keep = current_end > current_start + 0.2
        for x, y, start, end in zip(
            coords.loc[keep, "imagecol"],
            coords.loc[keep, "imagerow"],
            current_start[keep],
            current_end[keep],
        ):
            patches.append(Wedge((x, y), radius, start, end))
            colors.append(GROUP_COLORS[group])
        collection = PatchCollection(patches, facecolor=colors, edgecolor="none", linewidth=0, rasterized=True)
        ax.add_collection(collection)

    ax.set_aspect("equal")
    ax.set_xlim(coords["imagecol"].min() - 15, coords["imagecol"].max() + 15)
    ax.set_ylim(coords["imagerow"].max() + 15, coords["imagerow"].min() - 15)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_title("HCC-2T cell2location deconvolution pies")
    ax.text(
        0.02,
        0.02,
        "Each pie represents spot-level proportions\nfrom all 28 cell states grouped for display.",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=7,
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.75, pad=2),
    )
    legend_handles = [Patch(facecolor=GROUP_COLORS[g], edgecolor="none", label=g) for g in GROUP_ORDER]
    ax.legend(
        handles=legend_handles,
        loc="center left",
        bbox_to_anchor=(1.01, 0.5),
        frameon=False,
        ncol=1,
        handlelength=1.0,
        handletextpad=0.45,
    )
    panel_label(ax, "A")
    save(fig, "Fig1A_HCC2T_STdeconvolve_style_spatial_pies")
    plt.close(fig)


def clipped(values: pd.Series) -> tuple[np.ndarray, float, float]:
    arr = values.to_numpy(dtype=float)
    vmin = float(np.nanpercentile(arr, 1))
    vmax = float(np.nanpercentile(arr, 99))
    return np.clip(arr, vmin, vmax), vmin, vmax


def make_four_spatial_maps(abund: pd.DataFrame, coords: pd.DataFrame) -> None:
    features = [
        ("Tumor_OS-high", "Tumor OS-high", "magma"),
        ("Mac_TREM2_TAM", "Mac TREM2 TAM", "viridis"),
        ("Mac_SLC40A1_TAM", "Mac SLC40A1 TAM", "cividis"),
        ("Mac_HSP_TAM", "Mac HSP TAM", "YlGnBu"),
    ]
    fig, axes = plt.subplots(1, 4, figsize=(8.2, 2.55), constrained_layout=True)
    for ax, (feature, title, cmap) in zip(axes, features):
        vals, vmin, vmax = clipped(abund[feature])
        sc = ax.scatter(
            coords["imagecol"],
            coords["imagerow"],
            c=vals,
            s=5,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            linewidths=0,
            rasterized=True,
        )
        ax.set_title(title)
        ax.set_aspect("equal")
        ax.invert_yaxis()
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        cbar = fig.colorbar(sc, ax=ax, fraction=0.045, pad=0.01)
        cbar.ax.tick_params(length=2, width=0.5)
    panel_label(axes[0], "B")
    save(fig, "Fig1B_HCC2T_OShigh_and_three_macrophage_spatial_maps")
    plt.close(fig)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    abund = read_abundance()
    coords = read_coords(abund.index)
    make_spatial_pie(abund, coords)
    make_four_spatial_maps(abund, coords)
    print(f"Wrote revised Fig1A/B to {OUT_DIR}")


if __name__ == "__main__":
    main()
