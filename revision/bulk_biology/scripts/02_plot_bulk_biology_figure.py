#!/usr/bin/env python3
"""Create a standalone supplementary figure for TCGA OSARS-high bulk biology.

This script reads the precomputed statistics produced by 01_bulk_biology_enrichment.R.
It does not perform statistical tests.  All panels retain the frozen OSARS groups
and use the same pathway tables as their source of inference.
"""

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


HERE = Path(__file__).resolve()
ANALYSIS_DIR = HERE.parents[1]
TABLES = ANALYSIS_DIR / "tables"
DATA = ANALYSIS_DIR / "data" / "processed"
FIGURES = ANALYSIS_DIR / "figures"
FIGURES.mkdir(exist_ok=True)


# Colour-blind-safe Okabe-Ito family; warm = higher in OSARS-high, cool = lower.
HIGH = "#D55E00"
LOW = "#0072B2"
GREY = "#8A8A8A"
BLACK = "#252525"
PALE_GREY = "#E8E8E8"

mpl.rcParams.update(
    {
        "font.family": "Arial",
        "font.size": 8.5,
        "axes.titlesize": 11,
        "axes.labelsize": 9,
        "axes.linewidth": 0.7,
        "xtick.labelsize": 8,
        "ytick.labelsize": 7.4,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.12,
    }
)


def label_pathway(pathway: str) -> str:
    """Make MSigDB names compact while retaining the source ontology."""
    label = pathway.replace("HALLMARK_", "").replace("REACTOME_", "")
    label = label.replace("_", " ").lower().title()
    exact = {
        "Tgf Beta Signaling": "TGF-β signaling",
        "Il6 Jak Stat3 Signaling": "IL6–JAK–STAT3 signaling",
        "Il2 Stat5 Signaling": "IL2–STAT5 signaling",
        "Tnfa Signaling Via Nfkb": "TNFα signaling via NF-κB",
        "Wnt Beta Catenin Signaling": "WNT/β-catenin signaling",
        "P53 Pathway": "p53 pathway",
        "Uv Response Up": "UV response up",
        "Uv Response Dn": "UV response down",
        "Kras Signaling Up": "KRAS signaling up",
        "Kras Signaling Dn": "KRAS signaling down",
        "Interferon Alpha Response": "Interferon-α response",
        "Interferon Gamma Response": "Interferon-γ response",
    }
    if label in exact:
        return exact[label]
    replacements = {
        "E2F": "E2F",
        "G2M": "G2M",
        "Myc": "MYC",
        "Mtorc1": "mTORC1",
        "Pi3K Akt Mtor": "PI3K–AKT–mTOR",
        "Dna": "DNA",
        "Mrna": "mRNA",
        "Rho": "RHO",
        "Apc C": "APC/C",
        "Cdc20": "CDC20",
        "P450": "P450",
        "G1": "G1",
        "S Phase": "S phase",
        "M Phase": "M phase",
        "G2M Checkpoint": "G2M checkpoint",
    }
    for old, new in replacements.items():
        label = label.replace(old, new)
    return label


def fdr_text(value: float) -> str:
    if value < 1e-3:
        return f"q={value:.1e}"
    return f"q={value:.3f}"


def load_results():
    hallmark = pd.read_csv(TABLES / "TCGA_OSARS_high_vs_low_CAMERA_hallmark.csv")
    hallmark_stage = pd.read_csv(TABLES / "TCGA_OSARS_high_vs_low_stage_adjusted_CAMERA_hallmark.csv")
    reactome = pd.read_csv(TABLES / "TCGA_OSARS_high_vs_low_CAMERA_reactome.csv")
    gsva = pd.read_csv(DATA / "TCGA_Hallmark_GSVA_scores.csv")
    representative = pd.read_csv(TABLES / "representative_hallmark_pathways.csv")

    stage_cols = hallmark_stage[["pathway", "Direction", "FDR", "signed_log10_FDR"]].rename(
        columns={
            "Direction": "Direction_stage",
            "FDR": "FDR_stage",
            "signed_log10_FDR": "signed_log10_FDR_stage",
        }
    )
    hallmark = hallmark.merge(stage_cols, on="pathway", how="left", validate="one_to_one")
    hallmark["stage_significant"] = hallmark["FDR_stage"] < 0.05
    hallmark["label"] = hallmark["pathway"].map(label_pathway)
    reactome["label"] = reactome["pathway"].map(label_pathway)
    return hallmark, reactome, gsva, representative


def pathway_lollipop(ax, frame, title, subtitle, *, all_paths=False):
    """Signed pathway enrichment lollipop plot."""
    plot = frame.copy().sort_values("signed_log10_FDR", ascending=True).reset_index(drop=True)
    y = np.arange(len(plot))
    colors = np.where(plot["Direction"].eq("Up"), HIGH, LOW)
    ax.axvline(0, color=BLACK, lw=0.8, zorder=0)
    ax.hlines(y, 0, plot["signed_log10_FDR"], color=colors, lw=1.05, alpha=0.65, zorder=1)

    if "stage_significant" in plot.columns:
        # An outer dark ring makes stage-adjustment robustness visible without
        # reusing a separate redundant panel.
        outer = np.where(plot["stage_significant"], BLACK, "#B8B8B8")
        ax.scatter(
            plot["signed_log10_FDR"],
            y,
            s=30 + 0.55 * plot["NGenes"],
            c=colors,
            edgecolors=outer,
            linewidths=np.where(plot["stage_significant"], 0.75, 0.30),
            zorder=2,
        )
    else:
        ax.scatter(
            plot["signed_log10_FDR"],
            y,
            s=30 + 0.55 * plot["NGenes"],
            c=colors,
            edgecolors="white",
            linewidths=0.4,
            zorder=2,
        )

    ax.set_yticks(y)
    ax.set_yticklabels(plot["label"])
    ax.grid(axis="x", color=PALE_GREY, lw=0.7, zorder=0)
    ax.tick_params(axis="y", length=0, pad=3)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.set_xlabel("Signed -log10(CAMERA FDR)")
    ax.set_title(title, loc="left", pad=9, fontweight="bold")

    xmax = max(abs(plot["signed_log10_FDR"]).max() * 1.12, 1)
    ax.set_xlim(-xmax, xmax)
    top_y = len(plot) - 0.25
    ax.text(-xmax * 0.98, top_y, "lower in OSARS-high", ha="left", va="bottom", color=LOW, fontsize=7.4, fontweight="bold")
    ax.text(xmax * 0.98, top_y, "higher in OSARS-high", ha="right", va="bottom", color=HIGH, fontsize=7.4, fontweight="bold")

    if all_paths:
        legend = [
            Line2D([0], [0], marker="o", color="none", label="higher in OSARS-high", markerfacecolor=HIGH, markeredgecolor="white", markersize=6),
            Line2D([0], [0], marker="o", color="none", label="lower in OSARS-high", markerfacecolor=LOW, markeredgecolor="white", markersize=6),
            Line2D([0], [0], marker="o", color="none", label="retained after stage adjustment (q<0.05)", markerfacecolor="white", markeredgecolor=BLACK, markersize=6),
        ]
        ax.legend(handles=legend, frameon=False, fontsize=7.2, loc="lower left", bbox_to_anchor=(0, -0.075), ncol=1)
    return plot


def plot_gsva_panels(fig, spec, gsva, representative):
    """Show raw sample-level GSVA distributions for mechanically selected terms."""
    ordered = representative.copy()
    ordered["rank"] = np.where(ordered["Direction"].eq("Up"), 0, 1)
    ordered = ordered.sort_values(["rank", "FDR", "pathway"]).reset_index(drop=True)
    grid = spec.subgridspec(2, 2, hspace=0.48, wspace=0.34)
    rng = np.random.default_rng(20260906)
    count_low = int((gsva["OSARS_frozen_group"] == "Low").sum())
    count_high = int((gsva["OSARS_frozen_group"] == "High").sum())

    for i, term in ordered.iterrows():
        ax = fig.add_subplot(grid[i // 2, i % 2])
        pathway = term["pathway"]
        low_values = gsva.loc[gsva["OSARS_frozen_group"] == "Low", pathway].to_numpy()
        high_values = gsva.loc[gsva["OSARS_frozen_group"] == "High", pathway].to_numpy()
        values = [low_values, high_values]
        box = ax.boxplot(
            values,
            positions=[0, 1],
            widths=0.52,
            patch_artist=True,
            showfliers=False,
            medianprops={"color": BLACK, "linewidth": 1.05},
            whiskerprops={"color": "#686868", "linewidth": 0.75},
            capprops={"color": "#686868", "linewidth": 0.75},
        )
        for patch, color in zip(box["boxes"], [LOW, HIGH]):
            patch.set(facecolor=color, alpha=0.42, edgecolor=color, linewidth=1.0)
        for pos, sample_values, color in zip([0, 1], values, [LOW, HIGH]):
            jitter = rng.normal(0, 0.055, size=len(sample_values))
            ax.scatter(np.full(len(sample_values), pos) + jitter, sample_values, s=6, color=color, alpha=0.17, linewidths=0, rasterized=True)
        ax.set_xticks([0, 1])
        ax.set_xticklabels([f"Low\n(n={count_low})", f"High\n(n={count_high})"])
        ax.set_ylabel("GSVA score")
        ax.grid(axis="y", color=PALE_GREY, lw=0.7, zorder=0)
        ax.spines[["top", "right"]].set_visible(False)
        direction = "higher in OSARS-high" if term["Direction"] == "Up" else "lower in OSARS-high"
        ax.set_title(f"{label_pathway(pathway)}\n{direction}", loc="left", fontsize=8.8, fontweight="bold", pad=6)
    return ordered


def make_composite(hallmark, reactome, gsva, representative):
    fig = plt.figure(figsize=(17.7, 14.2))
    outer = fig.add_gridspec(
        2,
        2,
        width_ratios=[1.12, 1.0],
        height_ratios=[1.0, 1.0],
        wspace=0.35,
        hspace=0.35,
        left=0.08,
        right=0.98,
        top=0.94,
        bottom=0.10,
    )
    ax_a = fig.add_subplot(outer[:, 0])
    pathway_lollipop(
        ax_a,
        hallmark,
        "A  Hallmark pathways",
        "All 50 Hallmark gene sets; frozen OSARS-high versus OSARS-low TCGA-LIHC tumors",
        all_paths=True,
    )

    top_up = reactome.loc[reactome["Direction"].eq("Up")].nsmallest(10, "FDR")
    top_down = reactome.loc[reactome["Direction"].eq("Down")].nsmallest(10, "FDR")
    reactome_display = pd.concat([top_down, top_up], ignore_index=True)
    reactome_display = reactome_display.sort_values("signed_log10_FDR", ascending=True)
    reactome_display.to_csv(TABLES / "figure_reactome_top10_per_direction.csv", index=False)
    ax_b = fig.add_subplot(outer[0, 1])
    pathway_lollipop(
        ax_b,
        reactome_display,
        "B  Reactome pathways",
        "Top 10 terms per direction, ranked by primary CAMERA FDR",
    )
    ax_b.tick_params(axis="y", labelsize=6.65)

    plot_gsva_panels(fig, outer[1, 1], gsva, representative)
    fig.text(0.585, 0.513, "C  Representative sample-level pathway scores", ha="left", va="bottom", fontsize=11, fontweight="bold")
    fig.text(
        0.585,
        0.498,
        "Two lowest-FDR Hallmark terms in each direction; distributions are descriptive, while pathway inference is from CAMERA.",
        ha="left",
        va="bottom",
        fontsize=7.2,
        color="#555555",
    )

    fig.suptitle("Bulk transcriptional state associated with frozen OSARS-high status in TCGA-LIHC", x=0.08, ha="left", fontsize=15, fontweight="bold")
    fig.text(
        0.08,
        0.018,
        "TCGA-LIHC: low n=172, high n=171. Competitive pathway analysis: limma CAMERA; Benjamini–Hochberg FDR within each library. "
        "Black outline in A denotes pathways retained after adjustment for pathologic stage (n=321).",
        ha="left",
        va="bottom",
        fontsize=7.5,
        color="#4C4C4C",
    )
    save_figure(fig, "Supplementary_Figure_OSARS_high_bulk_biology_TCGA")
    plt.close(fig)


def save_figure(fig, stem):
    fig.savefig(FIGURES / f"{stem}.pdf")
    fig.savefig(FIGURES / f"{stem}.png", dpi=600)


def make_individual_panels(hallmark, reactome, gsva, representative):
    fig, ax = plt.subplots(figsize=(9.5, 14.2))
    pathway_lollipop(
        ax,
        hallmark,
        "Hallmark pathways associated with frozen OSARS-high status",
        "All 50 Hallmark gene sets; TCGA-LIHC OSARS-high versus OSARS-low",
        all_paths=True,
    )
    fig.text(0.125, 0.015, "Black outline: FDR<0.05 after pathologic-stage adjustment (n=321).", fontsize=7.5, color="#4C4C4C")
    save_figure(fig, "Panel_A_Hallmark_CAMERA_TCGA")
    plt.close(fig)

    top_up = reactome.loc[reactome["Direction"].eq("Up")].nsmallest(10, "FDR")
    top_down = reactome.loc[reactome["Direction"].eq("Down")].nsmallest(10, "FDR")
    reactome_display = pd.concat([top_down, top_up], ignore_index=True).sort_values("signed_log10_FDR")
    fig, ax = plt.subplots(figsize=(10.5, 8.2))
    pathway_lollipop(
        ax,
        reactome_display,
        "Reactome pathways associated with frozen OSARS-high status",
        "Top 10 terms per direction, ranked by primary CAMERA FDR; TCGA-LIHC",
    )
    ax.tick_params(axis="y", labelsize=7.1)
    save_figure(fig, "Panel_B_Reactome_CAMERA_TCGA")
    plt.close(fig)

    fig = plt.figure(figsize=(11.8, 7.6))
    spec = fig.add_gridspec(1, 1, left=0.08, right=0.98, top=0.88, bottom=0.08)
    plot_gsva_panels(fig, spec[0], gsva, representative)
    fig.suptitle("Representative Hallmark scores in frozen OSARS groups", x=0.08, ha="left", fontsize=14, fontweight="bold")
    fig.text(0.08, 0.025, "Terms selected mechanically as the two lowest-FDR Hallmark pathways in each direction; group distributions are descriptive.", fontsize=7.5, color="#4C4C4C")
    save_figure(fig, "Panel_C_Representative_Hallmark_GSVA_TCGA")
    plt.close(fig)


if __name__ == "__main__":
    hallmark_results, reactome_results, gsva_scores, representative_terms = load_results()
    make_composite(hallmark_results, reactome_results, gsva_scores, representative_terms)
    make_individual_panels(hallmark_results, reactome_results, gsva_scores, representative_terms)
    print(f"Figures written to: {FIGURES}")
