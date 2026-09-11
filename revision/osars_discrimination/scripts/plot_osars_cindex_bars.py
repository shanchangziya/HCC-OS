#!/usr/bin/env python3
"""Nature-style OSARS C-index bar chart for TCGA-LIHC and ICGC-LIRI."""

from pathlib import Path

import matplotlib as mpl
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


OUT = Path(__file__).resolve().parents[1]
PROCESSED = OUT.parent / "model_validation/data/processed"
FIG_DIR = OUT / "figures"
DATA_DIR = OUT / "data"

ARIAL = "/System/Library/Fonts/Supplemental/Arial.ttf"
TIMES_BOLD = "/System/Library/Fonts/Supplemental/Times New Roman Bold.ttf"


def mm(value: float) -> float:
    return value / 25.4


def setup_style() -> None:
    if Path(ARIAL).is_file():
        fm.fontManager.addfont(ARIAL)
    if Path(TIMES_BOLD).is_file():
        fm.fontManager.addfont(TIMES_BOLD)
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.sans-serif": ["Arial", "DejaVu Sans"],
            "font.size": 7.5,
            "axes.labelsize": 8.2,
            "axes.titlesize": 8.8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 0.75,
            "xtick.labelsize": 7.2,
            "ytick.labelsize": 7.0,
            "xtick.major.width": 0.75,
            "ytick.major.width": 0.75,
            "xtick.major.size": 2.8,
            "ytick.major.size": 2.8,
            "legend.frameon": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )


def load_frozen_results():
    metrics = pd.read_csv(PROCESSED / "discrimination_metrics_with_CI.csv")
    selected = metrics.loc[
        metrics["predictor"].eq("OSARS_frozen")
        & metrics["metric"].eq("C-index")
        & metrics["cohort"].isin(["TCGA-LIHC", "ICGC-LIRI"])
    ].copy()
    selected = selected.set_index("cohort").loc[["TCGA-LIHC", "ICGC-LIRI"]].reset_index()
    if len(selected) != 2:
        raise AssertionError("Expected one frozen OSARS C-index row per plotted cohort")
    expected = {
        "TCGA-LIHC": (343, 124, 0.8704939919893191),
        "ICGC-LIRI": (243, 44, 0.7664233576642335),
    }
    for _, row in selected.iterrows():
        n, events, estimate = expected[row["cohort"]]
        if int(row["n"]) != n or int(row["events"]) != events:
            raise AssertionError(f"Frozen sample manifest changed for {row['cohort']}")
        if not np.isclose(row["estimate"], estimate, rtol=0, atol=1e-12):
            raise AssertionError(f"Frozen C-index changed for {row['cohort']}")
    selected["role"] = [
        "model-development cohort; apparent discrimination",
        "retrospective external evaluation",
    ]

    pdc = pd.read_csv(PROCESSED / "PDC_verified158_discrimination_CI.csv")
    pdc = pdc.loc[
        pdc["cohort"].eq("PDC000198_OS_verified158")
        & pdc["predictor"].eq("OSARS_protein_surrogate")
        & pdc["metric"].eq("C-index")
    ].copy()
    if len(pdc) != 1:
        raise AssertionError("Expected one verified CPTAC/PDC OS C-index row")
    pdc_row = pdc.iloc[0]
    if int(pdc_row["n"]) != 158 or int(pdc_row["events"]) != 56:
        raise AssertionError("Verified CPTAC/PDC tumor cohort changed")
    return selected, pdc


def plot_bars(data: pd.DataFrame) -> plt.Figure:
    setup_style()
    fig, ax = plt.subplots(figsize=(mm(85), mm(82)))
    fig.subplots_adjust(left=0.19, right=0.965, bottom=0.19, top=0.82)

    x = np.arange(len(data))
    estimates = data["estimate"].to_numpy(float)
    yerr = np.vstack(
        [
            estimates - data["ci_lower"].to_numpy(float),
            data["ci_upper"].to_numpy(float) - estimates,
        ]
    )
    colors = ["#3DB7CC", "#E64B35"]
    bars = ax.bar(
        x,
        estimates,
        width=0.58,
        color=colors,
        edgecolor="#222222",
        linewidth=0.8,
        zorder=2,
    )
    ax.errorbar(
        x,
        estimates,
        yerr=yerr,
        fmt="none",
        ecolor="#222222",
        elinewidth=1.0,
        capsize=3.5,
        capthick=1.0,
        zorder=3,
    )
    for bar, estimate, upper in zip(bars, estimates, data["ci_upper"]):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            upper + 0.027,
            f"{estimate:.3f}",
            ha="center",
            va="bottom",
            fontsize=7.7,
            fontweight="bold",
            color="#222222",
        )

    ax.axhline(0.5, color="#B8B8B8", linewidth=0.7, linestyle=(0, (3, 3)), zorder=1)
    ax.set_ylim(0, 1.04)
    ax.set_xlim(-0.55, 1.55)
    ax.set_ylabel("Harrell C-index", labelpad=4)
    ax.set_xticks(x)
    ax.set_xticklabels(["TCGA-LIHC", "ICGC-LIRI"])
    ax.set_yticks(np.arange(0, 1.01, 0.2))
    ax.tick_params(axis="both", direction="out", pad=2.5)
    ax.spines["left"].set_color("#222222")
    ax.spines["bottom"].set_color("#222222")
    ax.grid(False)
    ax.set_title("C-index of OSARS", fontweight="bold", pad=8)

    panel_font = (fm.FontProperties(fname=TIMES_BOLD, size=11) if Path(TIMES_BOLD).is_file()
                  else fm.FontProperties(family="DejaVu Serif", weight="bold", size=11))
    ax.text(
        -0.22,
        1.13,
        "A",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontproperties=panel_font,
        color="black",
        clip_on=False,
    )
    return fig


def save(fig: plt.Figure, stem: Path) -> None:
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches=None)
    fig.savefig(stem.with_suffix(".svg"), bbox_inches=None)
    fig.savefig(stem.with_suffix(".png"), dpi=600, bbox_inches=None)


def main() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    plotted, pdc = load_frozen_results()
    plotted.to_csv(DATA_DIR / "OSARS_Cindex_TCGA_ICGC_source_data.csv", index=False)
    pdc.to_csv(DATA_DIR / "CPTAC_PDC_verified158_OS_Cindex.csv", index=False)

    pdc_row = pdc.iloc[0]
    report = (
        "CPTAC/PDC000198 verified tumor cohort (overall survival; n=158; 56 events): "
        f"C-index = {pdc_row['estimate']:.3f} "
        f"(95% CI, {pdc_row['ci_lower']:.3f}–{pdc_row['ci_upper']:.3f}).\n"
        "The score is the frozen correlation-weighted protein OSARS surrogate, not the mRNA GBM predictor.\n"
    )
    (DATA_DIR / "CPTAC_Cindex_text_report.txt").write_text(report, encoding="utf-8")

    fig = plot_bars(plotted)
    save(fig, FIG_DIR / "Panel_A_OSARS_Cindex_TCGA_ICGC")
    plt.close(fig)


if __name__ == "__main__":
    main()
