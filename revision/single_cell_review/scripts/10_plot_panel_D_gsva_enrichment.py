#!/usr/bin/env python3
"""Draw the review candidate for Panel D.

The left subpanel documents the prespecified GSVA-based selection among all
eight GeneNMF metaprograms.  The right subpanel shows only the biology of the
selected program.  All statistics are read from source-data tables produced by
09_score_mp_gsva_and_enrich_selected.R; this script performs no inference.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
OUT_ROOT = SCRIPT_DIR.parent
FIG_DIR = OUT_ROOT / "figures" / "panel_review"
FIG_DIR.mkdir(parents=True, exist_ok=True)

MM_PER_INCH = 25.4
WIDTH_MM = 170
HEIGHT_MM = 70

PANEL_FONT = "Times New Roman"
TEXT_FONT = "Arial"

SELECTED = "#B64E3C"
OTHER = "#B8BEC1"
ZERO = "#D9DDDF"
TEXT = "#222222"
RESOURCE_COLORS = {"Reactome": "#397B87", "Hallmark": "#C56A43"}

TERM_LABELS = {
    "HALLMARK_COAGULATION": "Coagulation",
    "REACTOME_COMPLEMENT_CASCADE": "Complement cascade",
    "REACTOME_DRUG_ADME": "Drug ADME",
    "REACTOME_RESPONSE_TO_ELEVATED_PLATELET_CYTOSOLIC_CA2": "Platelet Ca2+ response",
    "REACTOME_METABOLISM_OF_FAT_SOLUBLE_VITAMINS": "Fat-soluble vitamin metabolism",
    "HALLMARK_XENOBIOTIC_METABOLISM": "Xenobiotic metabolism",
}


def short_title(ax: mpl.axes.Axes, title: str, color: str = TEXT) -> None:
    ax.text(
        0,
        1.035,
        title,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=7.2,
        color=color,
        fontweight="normal",
    )


def main() -> None:
    association = pd.read_csv(OUT_ROOT / "GSE149614_MP_OS_association.csv")
    pathways = pd.read_csv(OUT_ROOT / "GSE149614_Panel_D_pathways_review_v2.csv")
    selection = pd.read_csv(OUT_ROOT / "GSE149614_OS_related_MP_selection.csv")

    selected_mp = str(selection.loc[0, "selected_metaprogram"])
    if selected_mp not in {f"MP{i}" for i in range(1, 9)}:
        raise ValueError(f"Unexpected audited metaprogram selection: {selected_mp}.")

    assoc = association.loc[
        association["gene_set_variant"].eq("OS-gene-excluded"),
        [
            "metaprogram",
            "mean_patient_spearman",
            "spearman_lower_95CI",
            "spearman_upper_95CI",
            "BH_FDR_mean_spearman",
            "selected",
        ],
    ].copy()
    if len(assoc) != 8 or assoc["metaprogram"].nunique() != 8:
        raise ValueError("The primary GSVA association table must contain eight metaprograms.")
    assoc = assoc.sort_values(
        ["mean_patient_spearman", "metaprogram"], ascending=[False, True]
    ).reset_index(drop=True)
    if assoc.iloc[0]["metaprogram"] != selected_mp:
        raise ValueError("The plotted selection is not the strongest primary GSVA association.")

    if len(pathways) != 6 or not pathways["selected_metaprogram"].eq(selected_mp).all():
        raise ValueError("Panel D must contain six enrichment terms for the selected metaprogram.")
    if not pathways["BH_FDR_all_resources"].lt(0.05).all():
        raise ValueError("A displayed pathway does not pass global BH FDR < 0.05.")
    missing_labels = set(pathways["gene_set"]) - set(TERM_LABELS)
    if missing_labels:
        raise ValueError(f"Missing concise labels for: {sorted(missing_labels)}")

    # Save the exact plotting subsets so every visual mark is traceable.
    assoc.to_csv(OUT_ROOT / "GSE149614_Panel_D_GSVA_association_review_v1.csv", index=False)
    plot_pathways = pathways.copy()
    plot_pathways["display_term"] = plot_pathways["gene_set"].map(TERM_LABELS)
    plot_pathways.to_csv(
        OUT_ROOT / "GSE149614_Panel_D_enrichment_review_v1.csv", index=False
    )

    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": [TEXT_FONT, "Helvetica", "DejaVu Sans"],
            "font.size": 6.2,
            "axes.labelsize": 6.2,
            "xtick.labelsize": 5.8,
            "ytick.labelsize": 6.0,
            "legend.fontsize": 5.7,
            "axes.linewidth": 0.45,
            "xtick.major.width": 0.45,
            "ytick.major.width": 0.45,
            "xtick.major.size": 2.2,
            "ytick.major.size": 0,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )

    fig = plt.figure(
        figsize=(WIDTH_MM / MM_PER_INCH, HEIGHT_MM / MM_PER_INCH),
        facecolor="white",
    )

    # Compact evidence of MP selection.
    ax_left = fig.add_axes([0.065, 0.175, 0.275, 0.705])
    y = np.arange(len(assoc), dtype=float)
    values = assoc["mean_patient_spearman"].to_numpy(float)
    lowers = assoc["spearman_lower_95CI"].to_numpy(float)
    uppers = assoc["spearman_upper_95CI"].to_numpy(float)
    colors = np.where(assoc["metaprogram"].eq(selected_mp), SELECTED, OTHER)

    ax_left.axvline(0, color=ZERO, lw=0.65, zorder=0)
    for yi, lo, hi, val, color in zip(y, lowers, uppers, values, colors):
        ax_left.hlines(yi, lo, hi, color=color, lw=1.05, zorder=2)
        ax_left.scatter(
            val,
            yi,
            s=18 if color == SELECTED else 11,
            color=color,
            edgecolor="white",
            linewidth=0.35,
            zorder=3,
        )

    ax_left.set_yticks(y, assoc["metaprogram"])
    ax_left.invert_yaxis()
    ax_left.set_xlim(-0.95, 0.95)
    ax_left.set_xticks([-0.8, -0.4, 0, 0.4, 0.8])
    ax_left.set_xlabel("Mean patient ρ", labelpad=2.2)
    ax_left.tick_params(axis="x", pad=2)
    ax_left.tick_params(axis="y", pad=3)
    ax_left.spines["top"].set_visible(False)
    ax_left.spines["right"].set_visible(False)
    ax_left.spines["left"].set_visible(False)
    for label in ax_left.get_yticklabels():
        if label.get_text() == selected_mp:
            label.set_color(SELECTED)
            label.set_fontweight("bold")
        else:
            label.set_color("#666B6E")
    short_title(ax_left, "OS association")

    # Biological function of only the locked program.
    ax_right = fig.add_axes([0.585, 0.175, 0.385, 0.705])
    plot_pathways = plot_pathways.sort_values(
        ["BH_FDR_all_resources", "gene_set"], ascending=[True, True]
    ).reset_index(drop=True)
    yp = np.arange(len(plot_pathways), dtype=float)
    fdr_strength = -np.log10(plot_pathways["BH_FDR_all_resources"].to_numpy(float))
    bar_colors = [RESOURCE_COLORS[x] for x in plot_pathways["resource"]]
    ax_right.barh(yp, fdr_strength, height=0.58, color=bar_colors, linewidth=0)
    ax_right.set_yticks(yp, plot_pathways["display_term"])
    ax_right.invert_yaxis()
    ax_right.set_xlabel("−log10(FDR)", labelpad=2.2)
    ax_right.set_xlim(0, max(7.5, float(fdr_strength.max()) + 0.35))
    ax_right.set_xticks([0, 2, 4, 6])
    ax_right.tick_params(axis="x", pad=2)
    ax_right.tick_params(axis="y", pad=3)
    ax_right.spines["top"].set_visible(False)
    ax_right.spines["right"].set_visible(False)
    ax_right.spines["left"].set_visible(False)
    short_title(ax_right, f"{selected_mp} enrichment", color=SELECTED)
    ax_right.legend(
        handles=[
            Patch(facecolor=RESOURCE_COLORS["Reactome"], edgecolor="none", label="Reactome"),
            Patch(facecolor=RESOURCE_COLORS["Hallmark"], edgecolor="none", label="Hallmark"),
        ],
        loc="upper right",
        bbox_to_anchor=(1.0, 1.035),
        ncol=2,
        frameon=False,
        borderaxespad=0,
        handlelength=0.9,
        handleheight=0.7,
        handletextpad=0.35,
        columnspacing=0.9,
    )

    fig.text(
        0.014,
        0.975,
        "D",
        ha="left",
        va="top",
        fontsize=10.5,
        fontfamily=PANEL_FONT,
        fontweight="bold",
        color="#111111",
    )

    stem = FIG_DIR / f"Panel_D_{selected_mp}_GSVA_enrichment_review_v1"
    fig.savefig(stem.with_suffix(".pdf"), dpi=600, bbox_inches=None)
    fig.savefig(stem.with_suffix(".svg"), dpi=600, bbox_inches=None)
    fig.savefig(stem.with_suffix(".png"), dpi=600, bbox_inches=None)
    plt.close(fig)

    print(f"Saved {stem}.pdf/.svg/.png")
    print(
        "Selected:",
        selected_mp,
        "rho=",
        f"{assoc.iloc[0]['mean_patient_spearman']:.2f}",
        "FDR=",
        f"{assoc.iloc[0]['BH_FDR_mean_spearman']:.4f}",
    )


if __name__ == "__main__":
    main()
