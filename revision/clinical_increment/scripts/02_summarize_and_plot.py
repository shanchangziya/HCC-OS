#!/usr/bin/env python3
"""Summarize and plot TCGA Stage/ROS/OSARS discrimination models."""

from pathlib import Path

import matplotlib as mpl
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ANALYSIS_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = ANALYSIS_DIR / "data"
FIG_DIR = ANALYSIS_DIR / "figures"
B = 1000
SEED = 20260905
ARIAL = "/System/Library/Fonts/Supplemental/Arial.ttf"
TIMES_BOLD = "/System/Library/Fonts/Supplemental/Times New Roman Bold.ttf"

MODEL_ORDER_ALL = [
    "ROS only",
    "OSARS only",
    "Stage",
    "Stage + ROS",
    "Stage + OSARS",
    "Stage + ROS + OSARS",
]
MODEL_ORDER_PLOT = [
    "ROS only",
    "OSARS only",
    "Stage",
    "Stage + ROS",
    "Stage + OSARS",
    "Stage + ROS + OSARS",
]
COLORS = {
    "ROS only": "#4C78A8",
    "OSARS only": "#D55E00",
    "Stage": "#747474",
    "Stage + ROS": "#4C78A8",
    "Stage + OSARS": "#D55E00",
    "Stage + ROS + OSARS": "#9C3B30",
}


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
            "axes.spines.right": False,
            "axes.spines.top": False,
            "axes.linewidth": 0.7,
            "xtick.labelsize": 7.3,
            "ytick.labelsize": 7.3,
            "xtick.major.size": 3.0,
            "ytick.major.size": 3.0,
            "xtick.major.width": 0.7,
            "ytick.major.width": 0.7,
            "legend.frameon": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )


def load_predictions() -> pd.DataFrame:
    path = DATA_DIR / "TCGA_stage_ROS_OSARS_model_predictions.csv"
    long = pd.read_csv(path)
    required = {"ID", "time", "event", "model", "linear_predictor"}
    if not required.issubset(long.columns):
        raise ValueError("Model prediction file is missing required columns")
    wide = long.pivot(
        index=["ID", "time", "event"], columns="model", values="linear_predictor"
    ).reset_index()
    if len(wide) != 321 or int(wide["event"].sum()) != 110:
        raise AssertionError("Frozen TCGA complete-case cohort changed")
    if wide[MODEL_ORDER_ALL].isna().any().any():
        raise AssertionError("Incomplete model predictions")
    return wide


def paired_cindex_bootstrap(
    time: np.ndarray,
    event: np.ndarray,
    score_matrix: np.ndarray,
    model_names: list[str],
):
    """Harrell C-index and paired row-bootstrap draws with fitted scores fixed."""
    time = np.asarray(time, float)
    event = np.asarray(event, bool)
    score_matrix = np.asarray(score_matrix, float)
    comparable = event[:, None] & (
        (time[:, None] < time[None, :])
        | ((time[:, None] == time[None, :]) & ~event[None, :])
    )
    earlier, later = np.where(comparable)
    difference = score_matrix[earlier] - score_matrix[later]
    agreement = (difference > 1e-8) + 0.5 * (np.abs(difference) <= 1e-8)
    point = agreement.mean(axis=0)

    rng = np.random.default_rng(SEED)
    n = len(time)
    draws = np.full((B, len(model_names)), np.nan)
    for b in range(B):
        index = rng.integers(n, size=n)
        counts = np.bincount(index, minlength=n)
        pair_weight = counts[earlier] * counts[later]
        if pair_weight.sum() > 0:
            draws[b] = np.sum(pair_weight[:, None] * agreement, axis=0) / pair_weight.sum()
    return point, draws


def summarize(wide: pd.DataFrame):
    matrix = wide[MODEL_ORDER_ALL].to_numpy(float)
    point, cube = paired_cindex_bootstrap(
        wide["time"].to_numpy(float),
        wide["event"].to_numpy(int),
        matrix,
        MODEL_ORDER_ALL,
    )
    fit_audit = pd.read_csv(DATA_DIR / "TCGA_stage_ROS_OSARS_model_fits.csv").set_index("model")
    metrics_rows = []
    draw_rows = []
    for j, name in enumerate(MODEL_ORDER_ALL):
        valid = cube[:, j][np.isfinite(cube[:, j])]
        low, high = np.quantile(valid, [0.025, 0.975])
        if not np.isclose(point[j], fit_audit.loc[name, "apparent_Harrell_C"], atol=1e-12):
            raise AssertionError(f"C-index implementation mismatch for {name}")
        metrics_rows.append(
            {
                "cohort": "TCGA-LIHC",
                "predictor": name,
                "metric": "C-index",
                "estimate": point[j],
                "ci_lower": low,
                "ci_upper": high,
                "n": len(wide),
                "events": int(wide["event"].sum()),
                "bootstrap_B": B,
                "bootstrap_valid": len(valid),
                "ci_method": "paired sample-row percentile bootstrap; fitted scores fixed",
                "evaluation": "apparent discrimination",
                "plotted": name in MODEL_ORDER_PLOT,
            }
        )
        draw_rows.extend(
            {
                "cohort": "TCGA-LIHC",
                "replicate": b + 1,
                "predictor": name,
                "metric": "C-index",
                "value": value,
            }
            for b, value in enumerate(cube[:, j])
        )
    metrics = pd.DataFrame(metrics_rows)
    draws = pd.DataFrame(draw_rows)

    point = metrics.set_index("predictor")["estimate"].to_dict()
    positions = {name: i for i, name in enumerate(MODEL_ORDER_ALL)}
    contrasts = [
        ("Stage + ROS", "Stage"),
        ("Stage + OSARS", "Stage"),
        ("Stage + ROS + OSARS", "Stage"),
        ("Stage + OSARS", "Stage + ROS"),
        ("Stage + ROS + OSARS", "Stage + OSARS"),
        ("OSARS only", "ROS only"),
    ]
    rows = []
    for model, reference in contrasts:
        values = cube[:, positions[model]] - cube[:, positions[reference]]
        values = values[np.isfinite(values)]
        low, high = np.quantile(values, [0.025, 0.975])
        rows.append(
            {
                "cohort": "TCGA-LIHC",
                "model": model,
                "reference": reference,
                "delta_Cindex": point[model] - point[reference],
                "bootstrap_95CI_lower": low,
                "bootstrap_95CI_upper": high,
                "bootstrap_B": B,
                "bootstrap_valid": len(values),
                "CI_excludes_zero": bool(low > 0 or high < 0),
            }
        )
    contrasts_df = pd.DataFrame(rows)
    return metrics, contrasts_df, draws


def draw_delta_bracket(ax, x1, x2, y, label, color="#555555") -> None:
    height = 0.007
    ax.plot([x1, x1, x2, x2], [y - height, y, y, y - height], color=color, lw=0.75, clip_on=False)
    ax.text((x1 + x2) / 2, y + 0.004, label, ha="center", va="bottom", fontsize=6.7, color=color)


def plot(metrics: pd.DataFrame, contrasts: pd.DataFrame) -> plt.Figure:
    setup_style()
    plot_data = (
        metrics.loc[metrics["predictor"].isin(MODEL_ORDER_PLOT)]
        .set_index("predictor")
        .loc[MODEL_ORDER_PLOT]
        .reset_index()
    )

    fig, ax = plt.subplots(figsize=(mm(170), mm(82)))
    fig.subplots_adjust(left=0.095, right=0.985, bottom=0.27, top=0.82)
    x = np.array([0.0, 1.0, 2.35, 3.35, 4.35, 5.35])

    ax.axhline(0.5, color="#B9B9B9", lw=0.7, ls=(0, (3, 3)), zorder=0)
    for idx, row in plot_data.iterrows():
        model = row["predictor"]
        color = COLORS[model]
        lower = row["estimate"] - row["ci_lower"]
        upper = row["ci_upper"] - row["estimate"]
        ax.errorbar(
            x[idx],
            row["estimate"],
            yerr=np.array([[lower], [upper]]),
            fmt="o",
            ms=6.0,
            mfc=color,
            mec="white",
            mew=0.65,
            ecolor=color,
            elinewidth=1.35,
            capsize=3.4,
            capthick=1.15,
            zorder=3,
        )
        ax.text(
            x[idx],
            row["ci_upper"] + 0.012,
            f"{row['estimate']:.3f}",
            ha="center",
            va="bottom",
            fontsize=7.2,
            fontweight="bold" if "OSARS" in model else "normal",
            color=color,
        )

    # Two clinically relevant nested increments; CIs are available in the source table.
    delta_ros = contrasts.loc[
        contrasts["model"].eq("Stage + ROS") & contrasts["reference"].eq("Stage")
    ].iloc[0]
    delta_both = contrasts.loc[
        contrasts["model"].eq("Stage + ROS + OSARS")
        & contrasts["reference"].eq("Stage + OSARS")
    ].iloc[0]
    draw_delta_bracket(ax, x[2], x[3], 0.746, f"ΔC = {delta_ros['delta_Cindex']:+.3f}")
    draw_delta_bracket(ax, x[4], x[5], 0.925, f"ΔC = {delta_both['delta_Cindex']:+.3f}")

    ax.axvline(1.68, color="#E2E2E2", lw=0.7, zorder=0)
    ax.set_xlim(-0.45, 5.78)
    ax.set_ylim(0.48, 0.955)
    ax.set_ylabel("Apparent Harrell C-index", labelpad=5)
    ax.set_xticks(x)
    ax.set_xticklabels(
        ["ROS only", "OSARS only", "Stage", "Stage + ROS", "Stage + OSARS", "Stage + ROS\n+ OSARS"]
    )
    ax.set_yticks([0.5, 0.6, 0.7, 0.8, 0.9])
    ax.tick_params(axis="both", direction="out", pad=3)
    ax.spines["bottom"].set_color("#222222")
    ax.spines["left"].set_color("#222222")
    ax.grid(False)

    ax.set_title("TCGA-LIHC", fontweight="bold", pad=8)
    ax.text(
        1.0,
        1.035,
        "n = 321; events = 110",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=7.0,
        color="#555555",
    )
    panel_font = (fm.FontProperties(fname=TIMES_BOLD, size=11) if Path(TIMES_BOLD).is_file()
                  else fm.FontProperties(family="DejaVu Serif", weight="bold", size=11))
    ax.text(
        -0.087,
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
    wide = load_predictions()
    metrics, contrasts, draws = summarize(wide)
    metrics.to_csv(DATA_DIR / "TCGA_stage_ROS_OSARS_Cindex.csv", index=False)
    contrasts.to_csv(DATA_DIR / "TCGA_stage_ROS_OSARS_paired_differences.csv", index=False)
    draws.to_csv(DATA_DIR / "TCGA_stage_ROS_OSARS_Cindex_bootstrap_draws.csv", index=False)
    metrics.loc[metrics["plotted"]].to_csv(
        DATA_DIR / "Figure_source_data_TCGA_stage_ROS_OSARS_Cindex.csv", index=False
    )

    fig = plot(metrics, contrasts)
    save(fig, FIG_DIR / "TCGA_stage_ROS_OSARS_Cindex")
    save(fig, FIG_DIR / "Panel_A_TCGA_stage_ROS_OSARS_Cindex")
    plt.close(fig)


if __name__ == "__main__":
    main()
