#!/usr/bin/env python3
"""Plot frozen OSARS versus direct Hallmark ROS scores in TCGA and ICGC."""

from pathlib import Path

import matplotlib as mpl
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


OUT_DIR = Path(__file__).resolve().parents[1]
INPUT_DIR = OUT_DIR.parent / "model_validation/data/processed"
FIG_DIR = OUT_DIR / "figures"
DATA_DIR = OUT_DIR / "data"

TCGA_FILE = INPUT_DIR / "TCGA-LIHC_all_scores.csv"
ICGC_FILE = INPUT_DIR / "ICGC-LIRI_all_scores.csv"
STATS_FILE = INPUT_DIR / "OSARS_program_correlations.csv"

ARIAL = "/System/Library/Fonts/Supplemental/Arial.ttf"
TIMES_BOLD = "/System/Library/Fonts/Supplemental/Times New Roman Bold.ttf"

POINT_COLOR = "#6F8FAF"
LINE_COLOR = "#B24745"
BAND_COLOR = "#E7C4C1"
AXIS_COLOR = "#222222"


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
            "font.size": 7.5,
            "axes.labelsize": 8.0,
            "axes.titlesize": 8.5,
            "xtick.labelsize": 7.0,
            "ytick.labelsize": 7.0,
            "axes.linewidth": 0.65,
            "xtick.major.width": 0.65,
            "ytick.major.width": 0.65,
            "xtick.major.size": 2.8,
            "ytick.major.size": 2.8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "mathtext.fontset": "custom",
            "mathtext.rm": "Arial",
            "mathtext.it": "Arial:italic",
            "mathtext.bf": "Arial:bold",
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )


def format_p(value: float) -> str:
    if value < 0.001:
        exponent = int(np.floor(np.log10(value)))
        coefficient = value / (10 ** exponent)
        return rf"${coefficient:.2f}\times10^{{{exponent}}}$"
    if value < 0.01:
        return f"{value:.5f}".rstrip("0")
    return f"{value:.3f}"


def regression_with_ci(x: np.ndarray, y: np.ndarray, n_grid: int = 240):
    """Least-squares line and 95% CI for the fitted mean response."""
    slope, intercept = np.polyfit(x, y, 1)
    grid = np.linspace(x.min(), x.max(), n_grid)
    fit = intercept + slope * grid
    fitted = intercept + slope * x
    residual = y - fitted
    dof = len(x) - 2
    residual_se = np.sqrt(np.sum(residual**2) / dof)
    x_centered_ss = np.sum((x - x.mean()) ** 2)
    mean_se = residual_se * np.sqrt(
        1.0 / len(x) + ((grid - x.mean()) ** 2) / x_centered_ss
    )
    critical = stats.t.ppf(0.975, dof)
    return grid, fit, fit - critical * mean_se, fit + critical * mean_se


def load_and_validate():
    cohorts = {
        "TCGA": pd.read_csv(TCGA_FILE),
        "ICGC": pd.read_csv(ICGC_FILE),
    }
    expected = {"TCGA": "TCGA-LIHC", "ICGC": "ICGC-LIRI"}
    summary = pd.read_csv(STATS_FILE)
    summary = summary.loc[summary["program"].eq("ROS_mean_rank")].copy()

    rows = []
    source = []
    for short_name, frame in cohorts.items():
        needed = ["ID", "OSARS_frozen", "ROS_mean_rank"]
        if frame[needed].isna().any().any():
            raise ValueError(f"Missing plotting values in {short_name}")
        if frame["ID"].duplicated().any():
            raise ValueError(f"Duplicated sample IDs in {short_name}")

        x = frame["ROS_mean_rank"].to_numpy(float)
        y = frame["OSARS_frozen"].to_numpy(float)
        rho, p_value = stats.spearmanr(x, y)
        frozen = summary.loc[summary["cohort"].eq(expected[short_name])].iloc[0]
        if len(frame) != int(frozen["n"]):
            raise AssertionError(f"Sample-size mismatch for {short_name}")
        if not np.isclose(rho, frozen["spearman_rho"], rtol=0, atol=1e-12):
            raise AssertionError(f"Spearman rho mismatch for {short_name}")
        if not np.isclose(p_value, frozen["p"], rtol=1e-10, atol=1e-15):
            raise AssertionError(f"Spearman P mismatch for {short_name}")

        rows.append(
            {
                "cohort": short_name,
                "n": len(frame),
                "spearman_rho": rho,
                "p_value": p_value,
                "bootstrap_95CI_lower": frozen["ci_lower"],
                "bootstrap_95CI_upper": frozen["ci_upper"],
                "ROS_genes_used": int(frozen["genes_used"]),
            }
        )
        table = frame[needed].copy()
        table.insert(0, "cohort", short_name)
        source.append(table)

    return cohorts, pd.DataFrame(rows), pd.concat(source, ignore_index=True)


def style_axis(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color(AXIS_COLOR)
    ax.spines["bottom"].set_color(AXIS_COLOR)
    ax.tick_params(axis="both", colors=AXIS_COLOR, direction="out", pad=2.2)
    ax.grid(False)


def draw_panel(
    ax: plt.Axes,
    data: pd.DataFrame,
    stat_row: pd.Series,
    title: str,
    panel_label: str,
) -> None:
    x = data["ROS_mean_rank"].to_numpy(float)
    y = data["OSARS_frozen"].to_numpy(float)
    grid, fit, lower, upper = regression_with_ci(x, y)

    ax.fill_between(grid, lower, upper, color=BAND_COLOR, alpha=0.62, linewidth=0, zorder=1)
    ax.scatter(
        x,
        y,
        s=12,
        c=POINT_COLOR,
        alpha=0.56,
        edgecolors="white",
        linewidths=0.22,
        zorder=2,
    )
    ax.plot(grid, fit, color=LINE_COLOR, linewidth=1.25, zorder=3)

    x_pad = 0.045 * (x.max() - x.min())
    y_pad = 0.065 * (y.max() - y.min())
    ax.set_xlim(x.min() - x_pad, x.max() + x_pad)
    ax.set_ylim(y.min() - y_pad, y.max() + y_pad)
    ax.set_xlabel("Direct ROS score", labelpad=4)
    ax.set_ylabel("OSARS score", labelpad=4)
    ax.set_title(title, fontweight="bold", pad=7)

    rho = stat_row["spearman_rho"]
    ci_low = stat_row["bootstrap_95CI_lower"]
    ci_high = stat_row["bootstrap_95CI_upper"]
    p_value = stat_row["p_value"]
    n = int(stat_row["n"])
    annotation = (
        f"ρ = {rho:.3f} (95% CI {ci_low:.3f}–{ci_high:.3f})\n"
        f"P = {format_p(p_value)}; n = {n}"
    )
    ax.text(
        0.035,
        0.965,
        annotation,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=7.2,
        color=AXIS_COLOR,
        linespacing=1.35,
    )

    panel_font = (fm.FontProperties(fname=TIMES_BOLD, size=11) if Path(TIMES_BOLD).is_file()
                  else fm.FontProperties(family="DejaVu Serif", weight="bold", size=11))
    ax.text(
        -0.16,
        1.12,
        panel_label,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontproperties=panel_font,
        color="black",
        clip_on=False,
    )
    style_axis(ax)


def save_figure(fig: plt.Figure, stem: Path) -> None:
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches=None)
    fig.savefig(stem.with_suffix(".svg"), bbox_inches=None)
    fig.savefig(stem.with_suffix(".png"), dpi=600, bbox_inches=None)


def main() -> None:
    setup_style()
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    cohorts, stats_table, source_table = load_and_validate()
    stats_table.to_csv(DATA_DIR / "OSARS_direct_ROS_correlation_statistics.csv", index=False)
    source_table.to_csv(DATA_DIR / "OSARS_direct_ROS_correlation_source_data.csv", index=False)

    # Combined figure: exact 170-mm journal width with a compact two-panel layout.
    fig, axes = plt.subplots(1, 2, figsize=(mm(170), mm(74)))
    fig.subplots_adjust(left=0.085, right=0.985, bottom=0.205, top=0.84, wspace=0.29)
    for ax, cohort, label in zip(axes, ["TCGA", "ICGC"], ["A", "B"]):
        stat_row = stats_table.loc[stats_table["cohort"].eq(cohort)].iloc[0]
        draw_panel(ax, cohorts[cohort], stat_row, cohort, label)
    save_figure(fig, FIG_DIR / "OSARS_direct_ROS_correlation_TCGA_ICGC")
    plt.close(fig)

    # Standalone panels for later assembly.
    for cohort, label in [("TCGA", "A"), ("ICGC", "B")]:
        fig, ax = plt.subplots(1, 1, figsize=(mm(85), mm(74)))
        fig.subplots_adjust(left=0.18, right=0.97, bottom=0.205, top=0.84)
        stat_row = stats_table.loc[stats_table["cohort"].eq(cohort)].iloc[0]
        draw_panel(ax, cohorts[cohort], stat_row, cohort, label)
        save_figure(fig, FIG_DIR / f"Panel_{label}_{cohort}_OSARS_direct_ROS")
        plt.close(fig)


if __name__ == "__main__":
    main()
