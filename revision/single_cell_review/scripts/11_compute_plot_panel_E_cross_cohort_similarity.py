#!/usr/bin/env python3
"""Formal cross-cohort GeneNMF metaprogram similarity for Panel E.

Discovery programs are the frozen GSE202642 consensus gene lists.  Validation
programs are the independently refitted eight-program GSE149614 solution.  All
tests are conditioned on the exact 2,000-gene validation GeneNMF universe.

The primary figure removes every frozen Hallmark ROS scoring gene from both the
gene sets and the background.  A full-gene-set analysis is exported as a
sensitivity analysis.  Pairwise enrichment uses a one-sided Fisher/hypergeometric
test followed by BH correction across all 4 x 8 comparisons within each variant.
"""

from __future__ import annotations

from datetime import datetime
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, hypergeom, spearmanr


SCRIPT_DIR = Path(__file__).resolve().parent
OUT_ROOT = SCRIPT_DIR.parent
REVISION_ROOT = OUT_ROOT.parent
SOURCE_DIR = REVISION_ROOT / "single_cell" / "data" / "processed"
REFERENCE_DIR = REVISION_ROOT / "model_validation" / "references"
FIG_DIR = OUT_ROOT / "figures" / "panel_review"
FIG_DIR.mkdir(parents=True, exist_ok=True)

DISCOVERY_PATH = SOURCE_DIR / "discovery_NMF_genes.csv"
VALIDATION_PATH = OUT_ROOT / "GSE149614_GeneNMF_selected_metaprogram_genes.csv"
UNIVERSE_PATH = SOURCE_DIR / "NMF_eligible_gene_universe.csv"
ROS_PATH = REFERENCE_DIR / "HALLMARK_ROS_msigdb_v7.0_genes.txt"

MM_PER_INCH = 25.4
WIDTH_MM = 170
HEIGHT_MM = 70
PANEL_FONT = "Times New Roman"
TEXT_FONT = "Arial"

SELECTED = "#B64E3C"
OTHER = "#B8BEC1"
TEXT = "#222222"
GRID = "#FFFFFF"


def bh_adjust(pvalues: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted P values, preserving input order."""
    pvalues = np.asarray(pvalues, dtype=float)
    m = len(pvalues)
    order = np.argsort(pvalues)
    ranked = pvalues[order] * m / np.arange(1, m + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adjusted = np.empty(m, dtype=float)
    adjusted[order] = np.minimum(ranked, 1.0)
    return adjusted


def ordered_sets(table: pd.DataFrame, group_col: str, gene_col: str) -> dict[str, list[str]]:
    result: dict[str, list[str]] = {}
    for group, block in table.groupby(group_col, sort=False):
        genes = block[gene_col].astype(str).str.upper().str.strip().tolist()
        if len(genes) != len(set(genes)):
            raise ValueError(f"Duplicated genes in {group_col}={group}.")
        result[str(group)] = genes
    return result


def pairwise_similarity(
    discovery: dict[str, list[str]],
    validation: dict[str, list[str]],
    background: set[str],
    variant: str,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for discovery_mp in [f"MP{i}" for i in range(1, 5)]:
        dset = set(discovery[discovery_mp]) & background
        for validation_mp in [f"MP{i}" for i in range(1, 9)]:
            vset = set(validation[validation_mp]) & background
            shared = sorted(dset & vset)
            union = dset | vset
            overlap_n = len(shared)
            table = np.array(
                [
                    [overlap_n, len(dset) - overlap_n],
                    [len(vset) - overlap_n, len(background) - len(union)],
                ],
                dtype=int,
            )
            odds_ratio, fisher_p = fisher_exact(table, alternative="greater")
            hypergeom_p = hypergeom.sf(
                overlap_n - 1, len(background), len(vset), len(dset)
            )
            if not np.isclose(fisher_p, hypergeom_p, rtol=1e-9, atol=1e-300):
                raise RuntimeError("Fisher and hypergeometric P values disagree.")
            expected = len(dset) * len(vset) / len(background)
            rows.append(
                {
                    "gene_set_variant": variant,
                    "discovery_program": discovery_mp,
                    "validation_program": validation_mp,
                    "background_n": len(background),
                    "n_discovery": len(dset),
                    "n_validation": len(vset),
                    "overlap_n": overlap_n,
                    "overlap_genes": ";".join(shared),
                    "jaccard": overlap_n / len(union) if union else np.nan,
                    "discovery_recall": overlap_n / len(dset) if dset else np.nan,
                    "validation_fraction": overlap_n / len(vset) if vset else np.nan,
                    "expected_overlap": expected,
                    "fold_enrichment": overlap_n / expected if expected else np.nan,
                    "odds_ratio": odds_ratio,
                    "P_one_sided_fisher": fisher_p,
                    "test_method": "one-sided Fisher exact / hypergeometric enrichment",
                    "multiple_testing_family": "all 32 discovery-by-validation comparisons within gene_set_variant",
                }
            )
    result = pd.DataFrame(rows)
    result["BH_FDR_32"] = bh_adjust(result["P_one_sided_fisher"].to_numpy())
    result["significant_FDR_0.05"] = result["BH_FDR_32"] < 0.05
    return result


def draw_panel(primary: pd.DataFrame, summary: pd.Series) -> Path:
    d_order = [f"MP{i}" for i in range(1, 5)]
    v_order = [f"MP{i}" for i in range(1, 9)]
    matrix = (
        primary.pivot(
            index="discovery_program", columns="validation_program", values="jaccard"
        )
        .loc[d_order, v_order]
        .to_numpy(float)
    )
    fdr = (
        primary.pivot(
            index="discovery_program", columns="validation_program", values="BH_FDR_32"
        )
        .loc[d_order, v_order]
        .to_numpy(float)
    )

    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": [TEXT_FONT, "Helvetica", "DejaVu Sans"],
            "font.size": 6.2,
            "axes.labelsize": 6.2,
            "xtick.labelsize": 5.8,
            "ytick.labelsize": 5.9,
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

    cmap = LinearSegmentedColormap.from_list(
        "program_overlap", ["#F7F8F7", "#D7E5E7", "#8CB8BE", "#397B87"], N=256
    )
    norm = Normalize(vmin=0, vmax=0.30, clip=True)
    fig = plt.figure(
        figsize=(WIDTH_MM / MM_PER_INCH, HEIGHT_MM / MM_PER_INCH), facecolor="white"
    )

    # Global 4 x 8 mapping establishes that MP4-to-MP1 is not chosen from a
    # single favorable pair in isolation.
    ax_hm = fig.add_axes([0.085, 0.245, 0.455, 0.585])
    image = ax_hm.imshow(matrix, cmap=cmap, norm=norm, aspect="auto", interpolation="nearest")
    ax_hm.set_xticks(np.arange(8), v_order)
    ax_hm.set_yticks(np.arange(4), d_order)
    ax_hm.set_xlabel("Validation", labelpad=4.0)
    ax_hm.set_ylabel("Discovery", labelpad=4.0)
    ax_hm.tick_params(axis="x", length=0, pad=3)
    ax_hm.tick_params(axis="y", length=0, pad=3)
    for spine in ax_hm.spines.values():
        spine.set_visible(False)

    # Quiet white gutters make each comparison readable without heavy grids.
    ax_hm.set_xticks(np.arange(-0.5, 8, 1), minor=True)
    ax_hm.set_yticks(np.arange(-0.5, 4, 1), minor=True)
    ax_hm.grid(which="minor", color=GRID, linewidth=1.2)
    ax_hm.tick_params(which="minor", bottom=False, left=False)

    for row in range(4):
        for col in range(8):
            value = matrix[row, col]
            label = "–" if value < 0.005 else f"{value:.2f}"
            color = "white" if value >= 0.17 else "#323638"
            weight = "bold" if (row == 3 and col == 0) else "normal"
            ax_hm.text(
                col,
                row,
                label,
                ha="center",
                va="center",
                fontsize=5.8,
                color=color,
                fontweight=weight,
            )
            if fdr[row, col] < 0.05:
                ax_hm.scatter(
                    col + 0.34,
                    row - 0.30,
                    s=3.8,
                    color="white" if value >= 0.17 else "#3B3E40",
                    linewidth=0,
                    zorder=4,
                )

    # Warm outline reserves the accent color for the requested correspondence.
    ax_hm.add_patch(
        Rectangle(
            (-0.5, 2.5),
            1,
            1,
            fill=False,
            edgecolor=SELECTED,
            linewidth=1.35,
            clip_on=False,
        )
    )
    ax_hm.text(
        0,
        1.13,
        "Cross-cohort gene overlap",
        transform=ax_hm.transAxes,
        ha="left",
        va="bottom",
        fontsize=7.2,
        color=TEXT,
    )

    cax = fig.add_axes([0.355, 0.865, 0.185, 0.028])
    cbar = fig.colorbar(image, cax=cax, orientation="horizontal", ticks=[0, 0.15, 0.30])
    cbar.outline.set_linewidth(0.35)
    cbar.ax.tick_params(length=1.8, width=0.4, pad=1.5, labelsize=5.2)
    cbar.ax.set_title("Jaccard", fontsize=5.5, loc="left", pad=1.8)
    fig.text(0.087, 0.183, "•  FDR < 0.05", ha="left", va="center", fontsize=5.4, color="#3B3E40")

    # Formal comparison restricted to the discovery OS-associated MP4 row.
    mp4 = primary.loc[primary["discovery_program"].eq("MP4")].copy()
    mp4["validation_number"] = mp4["validation_program"].str.removeprefix("MP").astype(int)
    mp4 = mp4.sort_values("validation_number")
    strength = -np.log10(np.clip(mp4["BH_FDR_32"].to_numpy(float), 1e-300, 1.0))
    y = np.arange(len(mp4), dtype=float)
    colors = np.where(mp4["validation_program"].eq("MP1"), SELECTED, OTHER)

    ax_bar = fig.add_axes([0.655, 0.175, 0.305, 0.705])
    ax_bar.barh(y, strength, height=0.54, color=colors, linewidth=0)
    ax_bar.axvline(-np.log10(0.05), color="#90979A", lw=0.6, ls=(0, (2, 2)), zorder=0)
    ax_bar.set_yticks(y, mp4["validation_program"])
    ax_bar.invert_yaxis()
    ax_bar.set_xlim(0, 42)
    ax_bar.set_xticks([0, 10, 20, 30, 40])
    ax_bar.set_xlabel("−log10(FDR)", labelpad=2.2)
    ax_bar.tick_params(axis="x", pad=2)
    ax_bar.tick_params(axis="y", pad=3)
    ax_bar.spines["top"].set_visible(False)
    ax_bar.spines["right"].set_visible(False)
    ax_bar.spines["left"].set_visible(False)
    for label in ax_bar.get_yticklabels():
        if label.get_text() == "MP1":
            label.set_color(SELECTED)
            label.set_fontweight("bold")
        else:
            label.set_color("#666B6E")
    ax_bar.text(
        0,
        1.035,
        "Discovery MP4",
        transform=ax_bar.transAxes,
        ha="left",
        va="bottom",
        fontsize=7.2,
        color=SELECTED,
    )
    ax_bar.text(
        strength[0] + 0.7,
        y[0],
        f"{int(summary['overlap_n'])}/{int(summary['n_discovery'])}",
        ha="left",
        va="center",
        fontsize=5.8,
        color=SELECTED,
        fontweight="bold",
    )

    fig.text(
        0.014,
        0.975,
        "E",
        ha="left",
        va="top",
        fontsize=10.5,
        fontfamily=PANEL_FONT,
        fontweight="bold",
        color="#111111",
    )

    stem = FIG_DIR / "Panel_E_discovery_MP4_validation_MP1_similarity_review_v1"
    fig.savefig(stem.with_suffix(".pdf"), dpi=600, bbox_inches=None)
    fig.savefig(stem.with_suffix(".svg"), dpi=600, bbox_inches=None)
    fig.savefig(stem.with_suffix(".png"), dpi=600, bbox_inches=None)
    plt.close(fig)
    return stem


def main() -> None:
    for path in [DISCOVERY_PATH, VALIDATION_PATH, UNIVERSE_PATH, ROS_PATH]:
        if not path.exists():
            raise FileNotFoundError(path)

    discovery_table = pd.read_csv(DISCOVERY_PATH)
    validation_table = pd.read_csv(VALIDATION_PATH)
    universe = set(pd.read_csv(UNIVERSE_PATH)["gene"].astype(str).str.upper().str.strip())
    ros = {
        line.strip().upper()
        for line in ROS_PATH.read_text(encoding="utf-8").splitlines()
        if line.strip()
    }
    if len(universe) != 2000:
        raise ValueError("The frozen validation GeneNMF universe must contain 2,000 genes.")

    discovery = ordered_sets(discovery_table, "program", "gene")
    validation = ordered_sets(validation_table, "metaprogram", "gene")
    if set(discovery) != {f"MP{i}" for i in range(1, 5)}:
        raise ValueError("Expected four frozen discovery metaprograms.")
    if set(validation) != {f"MP{i}" for i in range(1, 9)}:
        raise ValueError("Expected eight independently refitted validation metaprograms.")

    full = pairwise_similarity(discovery, validation, universe, "full")
    no_ros_background = universe - ros
    no_ros = pairwise_similarity(
        discovery, validation, no_ros_background, "OS-score-genes-excluded"
    )
    all_results = pd.concat([no_ros, full], ignore_index=True)
    all_results = all_results.sort_values(
        ["gene_set_variant", "discovery_program", "validation_program"]
    )

    primary = no_ros.copy()
    mp4 = primary.loc[primary["discovery_program"].eq("MP4")].sort_values(
        ["jaccard", "BH_FDR_32"], ascending=[False, True]
    )
    selected = mp4.iloc[0]
    if selected["validation_program"] != "MP1":
        raise RuntimeError("Validation MP1 is not the strongest discovery MP4 match.")
    if int(selected["overlap_n"]) != 38 or float(selected["BH_FDR_32"]) >= 0.05:
        raise RuntimeError("The expected MP4-to-MP1 overlap was not reproduced.")
    global_best = primary.sort_values(
        ["jaccard", "BH_FDR_32"], ascending=[False, True]
    ).iloc[0]
    if not (
        global_best["discovery_program"] == "MP4"
        and global_best["validation_program"] == "MP1"
    ):
        raise RuntimeError("MP4-to-MP1 is not the global maximum of the 4 x 8 matrix.")

    # Rank concordance among the 38 overlapping genes is a descriptive second
    # axis of evidence; the Fisher-BH result remains the formal overlap test.
    discovery_rank = pd.DataFrame(
        {
            "gene": discovery["MP4"],
            "discovery_rank_raw": np.arange(1, len(discovery["MP4"]) + 1),
        }
    )
    discovery_rank = discovery_rank.loc[
        discovery_rank["gene"].isin(no_ros_background)
    ].copy()
    discovery_rank["discovery_rank_eligible"] = np.arange(1, len(discovery_rank) + 1)

    validation_rank = validation_table.loc[
        validation_table["metaprogram"].eq("MP1"),
        ["gene", "rank", "normalized_weight"],
    ].copy()
    validation_rank["gene"] = validation_rank["gene"].astype(str).str.upper().str.strip()
    validation_rank = validation_rank.loc[validation_rank["gene"].isin(no_ros_background)]
    validation_rank = validation_rank.rename(columns={"rank": "validation_rank"})

    shared = discovery_rank.merge(validation_rank, on="gene", how="inner", validate="one_to_one")
    shared = shared.sort_values("discovery_rank_eligible")
    if len(shared) != 38 or shared["gene"].isin(ros).any():
        raise RuntimeError("Shared-gene audit failed.")
    rank_test = spearmanr(
        shared["discovery_rank_eligible"], shared["validation_rank"]
    )
    shared["is_frozen_OS_score_gene"] = shared["gene"].isin(ros)

    second = mp4.iloc[1]
    summary = pd.Series(
        {
            "discovery_program": "MP4",
            "validation_program": "MP1",
            "primary_gene_set_variant": "OS-score-genes-excluded",
            "background_n": int(selected["background_n"]),
            "OS_score_genes_removed_from_background": len(universe & ros),
            "n_discovery": int(selected["n_discovery"]),
            "n_validation": int(selected["n_validation"]),
            "overlap_n": int(selected["overlap_n"]),
            "jaccard": float(selected["jaccard"]),
            "discovery_recall": float(selected["discovery_recall"]),
            "validation_fraction": float(selected["validation_fraction"]),
            "fold_enrichment": float(selected["fold_enrichment"]),
            "odds_ratio": float(selected["odds_ratio"]),
            "P_one_sided_fisher": float(selected["P_one_sided_fisher"]),
            "BH_FDR_32": float(selected["BH_FDR_32"]),
            "next_best_validation_program": str(second["validation_program"]),
            "next_best_jaccard": float(second["jaccard"]),
            "jaccard_margin_over_next_best": float(selected["jaccard"] - second["jaccard"]),
            "shared_gene_rank_spearman_rho": float(rank_test.statistic),
            "shared_gene_rank_spearman_P": float(rank_test.pvalue),
            "shared_frozen_OS_score_gene_n": 0,
        }
    )

    all_results.to_csv(
        OUT_ROOT / "GSE202642_GSE149614_GeneNMF_cross_cohort_similarity.csv", index=False
    )
    shared.to_csv(
        OUT_ROOT / "GSE202642_MP4_GSE149614_MP1_shared_genes.csv", index=False
    )
    summary.to_frame().T.to_csv(
        OUT_ROOT / "GSE202642_MP4_GSE149614_MP1_similarity_summary.csv", index=False
    )

    stem = draw_panel(primary, summary)
    manifest = {
        "analysis": "Formal cross-cohort GeneNMF program correspondence",
        "generated": datetime.now().astimezone().isoformat(timespec="seconds"),
        "discovery_dataset": "GSE202642",
        "validation_dataset": "GSE149614",
        "discovery_programs": 4,
        "validation_programs": 8,
        "primary_background_n": len(no_ros_background),
        "primary_background": "Frozen 2,000-gene validation GeneNMF universe after removing the 17 measurable frozen OS-score genes",
        "primary_gene_set_variant": "OS-score genes excluded from both programs and background",
        "sensitivity_gene_set_variant": "Full program gene sets in the complete 2,000-gene validation universe",
        "test": "One-sided Fisher exact / hypergeometric enrichment",
        "multiple_testing": "Benjamini-Hochberg across all 32 discovery-by-validation comparisons within each gene-set variant",
        "figure": str(stem.with_suffix(".pdf")),
        "selected_pair": summary.to_dict(),
        "limitation": "The test conditions on the validation GeneNMF candidate universe and treats the frozen discovery signatures as fixed; it is not a symmetric reconstruction of the discovery W-matrix candidate universe.",
    }
    (OUT_ROOT / "GSE202642_GSE149614_cross_cohort_similarity_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    print("PANEL_E_CROSS_COHORT_SIMILARITY_COMPLETE")
    print("Selected pair: discovery MP4 -> validation MP1")
    print(
        f"Overlap: {int(summary['overlap_n'])}/{int(summary['n_discovery'])}; "
        f"Jaccard={summary['jaccard']:.4f}; BH FDR={summary['BH_FDR_32']:.3e}"
    )
    print(
        f"Shared-gene rank rho={summary['shared_gene_rank_spearman_rho']:.3f}; "
        f"P={summary['shared_gene_rank_spearman_P']:.3e}"
    )
    print(f"Saved {stem}.pdf/.svg/.png")


if __name__ == "__main__":
    main()
