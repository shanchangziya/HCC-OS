"""Render the standalone GSE149614 validation supplementary figure.

The plotting layer reads only source-data tables prepared in this Reviewer 2
folder.  It performs no score, NMF, enrichment, or inferential calculation.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


HERE = Path(__file__).resolve()
ROOT = HERE.parents[1]
FIGDIR = ROOT / "figures"
PANELDIR = FIGDIR / "panels"
FIGDIR.mkdir(exist_ok=True)
PANELDIR.mkdir(parents=True, exist_ok=True)

# Nature-style sparse layout with editable vector text.  Panel letters alone
# use the explicitly requested Times New Roman typeface.
mpl.rcParams.update({
    "font.family": "Arial",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 6.8,
    "axes.titlesize": 7.2,
    "axes.labelsize": 6.8,
    "xtick.labelsize": 6.2,
    "ytick.labelsize": 6.2,
    "legend.fontsize": 6.1,
    "axes.linewidth": 0.65,
    "xtick.major.width": 0.65,
    "ytick.major.width": 0.65,
    "xtick.major.size": 2.4,
    "ytick.major.size": 2.4,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "figure.dpi": 600,
    "savefig.dpi": 600,
    "savefig.facecolor": "white",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "legend.frameon": False,
})

MM = 1 / 25.4
SCORE_CMAP = mpl.colormaps["magma"]
CELL_COLORS = {
    "Hepatocyte": "#D55E00",
    "T/NK": "#0072B2",
    "Myeloid": "#009E73",
    "B": "#CC79A7",
    "Endothelial": "#56B4E9",
    "Fibroblast": "#7A7A7A",
}
CELL_ORDER = ["Fibroblast", "Endothelial", "B", "T/NK", "Myeloid", "Hepatocyte"]
RESOURCE_COLORS = {"Hallmark": "#C5643F", "Reactome": "#3E7C88"}

ALL = pd.read_csv(ROOT / "GSE149614_all_cells_UMAP_metadata.csv")
OS = pd.read_csv(ROOT / "GSE149614_OS_score_metadata.csv")
PANEL_D = pd.read_csv(ROOT / "GSE149614_panel_D_pathways.csv")
SIM = pd.read_csv(ROOT / "GSE149614_MP_similarity_to_discovery.csv")

if len(ALL) != 71915 or len(OS) != 13691:
    raise RuntimeError("Unexpected source-data dimensions for GSE149614 figure.")
if not {"MP1_deNovo", "MP2_deNovo", "MP3_deNovo", "MP4_deNovo"}.issubset(OS.columns):
    raise RuntimeError("Missing validation NMF score columns.")
if not ((SIM["discovery_program"] == "MP4") & (SIM["validation_program"] == "MP2") &
        SIM["selected_for_panel_D"].astype(bool)).any():
    raise RuntimeError("Panel D must use the frozen MP4-to-validation MP2 mapping.")
if not PANEL_D["FDR_within_program_all_resources"].lt(0.05).all():
    raise RuntimeError("Panel D contains a non-significant pathway.")


def clean_umap(ax, anchor="C"):
    ax.set_aspect("equal")
    ax.set_anchor(anchor)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.margins(0.025)


def panel_letter(fig, x, y, letter):
    fig.text(x, y, letter, ha="left", va="top", fontsize=12.5, fontweight="bold",
             fontfamily="Times New Roman", color="#121212")


def add_short_title(ax, title, color="#222222"):
    ax.text(0, 1.015, title, transform=ax.transAxes, ha="left", va="bottom",
            fontsize=7.3, color=color)


def draw_all_cells(ax):
    for cell_type in CELL_ORDER:
        dat = ALL.loc[ALL["major_cell_type"] == cell_type]
        ax.scatter(dat["UMAP1"], dat["UMAP2"], s=0.32, color=CELL_COLORS[cell_type],
                   linewidths=0, rasterized=True)
    clean_umap(ax, anchor="W")
    add_short_title(ax, "Major cell types")
    handles = [
        Line2D([0], [0], marker="o", color="none", label=ct,
               markerfacecolor=CELL_COLORS[ct], markeredgewidth=0, markersize=3.4)
        for ct in ["Hepatocyte", "T/NK", "Myeloid", "B", "Endothelial", "Fibroblast"]
    ]
    ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.015, 0.50),
              borderaxespad=0, handletextpad=0.35, labelspacing=0.48)


def draw_os_score(ax):
    dat = OS.sort_values("OS_composite")
    points = ax.scatter(dat["UMAP1"], dat["UMAP2"], c=dat["OS_composite"],
                        cmap=SCORE_CMAP, norm=Normalize(vmin=-3.0, vmax=2.5),
                        s=0.62, linewidths=0, rasterized=True)
    clean_umap(ax)
    add_short_title(ax, "Malignant hepatocytes")
    return points


def draw_nmf_maps(fig, subgrid):
    axes = []
    for i, program in enumerate(["MP1", "MP2", "MP3", "MP4"]):
        ax = fig.add_subplot(subgrid[i // 2, i % 2])
        dat = OS.sort_values(f"{program}_deNovo")
        ax.scatter(dat["UMAP1"], dat["UMAP2"], c=dat[f"{program}_deNovo"],
                   cmap=SCORE_CMAP, norm=Normalize(vmin=0, vmax=1),
                   s=0.62, linewidths=0, rasterized=True)
        clean_umap(ax)
        add_short_title(ax, program, color="#C5643F" if program == "MP2" else "#222222")
        axes.append(ax)
    colorbar = fig.colorbar(ScalarMappable(norm=Normalize(0, 1), cmap=SCORE_CMAP),
                            ax=axes, orientation="horizontal", fraction=0.045,
                            pad=0.080, aspect=32)
    colorbar.set_label("UCell score", labelpad=1.5)
    colorbar.ax.tick_params(length=2, pad=1)
    return axes


def draw_enrichment(ax):
    dat = PANEL_D.sort_values("display_order", ascending=False).copy()
    y = np.arange(len(dat))
    colors = [RESOURCE_COLORS[x] for x in dat["resource"]]
    ax.barh(y, dat["minus_log10_FDR"], color=colors, height=0.62, linewidth=0)
    ax.set_yticks(y, dat["display_term"])
    ax.set_xlabel("−log10(FDR)")
    ax.set_xlim(0, max(6.8, dat["minus_log10_FDR"].max() + 0.35))
    ax.tick_params(axis="y", length=0, pad=3)
    ax.tick_params(axis="x", pad=2)
    ax.spines["left"].set_visible(False)
    add_short_title(ax, "MP2 enrichment")
    ax.legend(handles=[Patch(facecolor=RESOURCE_COLORS["Hallmark"], label="Hallmark"),
                       Patch(facecolor=RESOURCE_COLORS["Reactome"], label="Reactome")],
              loc="upper right", ncol=2, handlelength=0.9, columnspacing=0.8,
              borderaxespad=0.05)


def save_figure(fig, stem, directory=FIGDIR):
    for suffix, dpi in (("pdf", None), ("svg", None), ("png", 600)):
        kwargs = {"facecolor": "white", "bbox_inches": None, "pad_inches": 0}
        if dpi is not None:
            kwargs["dpi"] = dpi
        fig.savefig(directory / f"{stem}.{suffix}", **kwargs)
    plt.close(fig)


def make_panel_a():
    fig = plt.figure(figsize=(170 * MM, 66 * MM))
    ax = fig.add_axes([0.075, 0.10, 0.72, 0.84])
    draw_all_cells(ax)
    panel_letter(fig, 0.014, 0.978, "A")
    save_figure(fig, "SuppFig_GSE149614_validation_Panel_A", PANELDIR)


def make_panel_b():
    fig = plt.figure(figsize=(82 * MM, 83 * MM))
    ax = fig.add_axes([0.12, 0.085, 0.70, 0.84])
    points = draw_os_score(ax)
    colorbar = fig.colorbar(points, ax=ax, fraction=0.050, pad=0.025, shrink=0.72)
    colorbar.set_label("OS composite", labelpad=1.5)
    colorbar.ax.tick_params(length=2, pad=1)
    panel_letter(fig, 0.016, 0.978, "B")
    save_figure(fig, "SuppFig_GSE149614_validation_Panel_B", PANELDIR)


def make_panel_c():
    fig = plt.figure(figsize=(82 * MM, 83 * MM))
    grid = fig.add_gridspec(2, 2, left=0.10, right=0.94, bottom=0.14, top=0.91,
                            hspace=0.20, wspace=0.15)
    draw_nmf_maps(fig, grid)
    panel_letter(fig, 0.015, 0.978, "C")
    save_figure(fig, "SuppFig_GSE149614_validation_Panel_C", PANELDIR)


def make_panel_d():
    fig = plt.figure(figsize=(170 * MM, 57 * MM))
    ax = fig.add_axes([0.245, 0.17, 0.73, 0.69])
    draw_enrichment(ax)
    panel_letter(fig, 0.014, 0.978, "D")
    save_figure(fig, "SuppFig_GSE149614_validation_Panel_D", PANELDIR)


def make_composite():
    # Fixed manuscript width: 170 mm.  Height is intentionally unrestricted
    # per the user-specified large-figure convention (here 222 mm).
    fig = plt.figure(figsize=(170 * MM, 222 * MM))
    grid = fig.add_gridspec(
        3, 2, left=0.075, right=0.965, bottom=0.045, top=0.982,
        height_ratios=[0.98, 1.23, 0.75], hspace=0.32, wspace=0.16,
    )
    ax_a = fig.add_subplot(grid[0, :])
    draw_all_cells(ax_a)
    pos_a = ax_a.get_position()
    panel_letter(fig, 0.012, pos_a.y1 + 0.012, "A")

    ax_b = fig.add_subplot(grid[1, 0])
    points = draw_os_score(ax_b)
    colorbar_b = fig.colorbar(points, ax=ax_b, fraction=0.046, pad=0.026, shrink=0.73)
    colorbar_b.set_label("OS composite", labelpad=1.5)
    colorbar_b.ax.tick_params(length=2, pad=1)
    pos_b = ax_b.get_position()
    panel_letter(fig, pos_b.x0 - 0.045, pos_b.y1 + 0.014, "B")

    grid_c = grid[1, 1].subgridspec(2, 2, hspace=0.20, wspace=0.15)
    axes_c = draw_nmf_maps(fig, grid_c)
    pos_c = axes_c[0].get_position()
    panel_letter(fig, pos_c.x0 - 0.040, pos_c.y1 + 0.014, "C")

    d_slot = grid[2, :].get_position(fig)
    # Panel D needs a larger left gutter for full pathway names.  It is kept
    # within the same row while preserving the 170-mm composite canvas.
    ax_d = fig.add_axes([0.225, d_slot.y0, 0.740, d_slot.height])
    draw_enrichment(ax_d)
    pos_d = ax_d.get_position()
    panel_letter(fig, 0.012, pos_d.y1 + 0.014, "D")

    save_figure(fig, "SuppFig_GSE149614_validation", ROOT)


if __name__ == "__main__":
    make_panel_a()
    make_panel_b()
    make_panel_c()
    make_panel_d()
    make_composite()
    print("SUPPLEMENTARY_VALIDATION_FIGURE_RENDERED")
