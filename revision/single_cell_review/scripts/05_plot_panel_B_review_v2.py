"""Render only the revised Panel B candidate for author review.

All malignant hepatocytes form a pale spatial reference.  The frozen
within-cohort top-20% OS group is overlaid with its continuous OS composite.
No grouping, score, or statistic is recalculated in this plotting script.
"""
from pathlib import Path
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize


HERE = Path(__file__).resolve()
ROOT = HERE.parents[1]
OUT = ROOT / "figures" / "panel_review"
OUT.mkdir(parents=True, exist_ok=True)

mpl.rcParams.update({
    "font.family": "Arial",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7,
    "axes.titlesize": 7.5,
    "axes.labelsize": 7,
    "xtick.labelsize": 6.5,
    "ytick.labelsize": 6.5,
    "axes.linewidth": 0.7,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "figure.dpi": 600,
    "savefig.dpi": 600,
    "savefig.facecolor": "white",
})

MM = 1 / 25.4
BACKGROUND = "#E8E8E8"
OS_CMAP = LinearSegmentedColormap.from_list(
    "os_high_warm",
    ["#F6D7B0", "#EEA45F", "#DD6B3D", "#B83B32", "#741F2B"],
)
OS_HIGH_CUTOFF = 0.79170559244785
OS_HIGH_MAX = 2.70203002818365

data = pd.read_csv(ROOT / "GSE149614_OS_score_metadata.csv")
required = {"UMAP1", "UMAP2", "OS_composite", "OS_high_top20", "published_malignancy"}
if not required.issubset(data.columns):
    raise RuntimeError("Panel B source table is incomplete.")
if len(data) != 13691 or not (data["published_malignancy"] == "Malignant hepatocyte").all():
    raise RuntimeError("Panel B must contain the frozen malignant-hepatocyte subset only.")

high = data.loc[data["OS_high_top20"].astype(bool)].sort_values("OS_composite")
if len(high) != 2739 or abs(high["OS_composite"].min() - OS_HIGH_CUTOFF) > 1e-10:
    raise RuntimeError("Frozen OS-high membership or cutoff has changed.")

fig = plt.figure(figsize=(85 * MM, 80 * MM))
ax = fig.add_axes([0.09, 0.075, 0.73, 0.865])

# Pale context first; only the prespecified OS-high cells carry signal colour.
ax.scatter(data["UMAP1"], data["UMAP2"], s=0.48, c=BACKGROUND,
           alpha=0.70, linewidths=0, rasterized=True, zorder=1)
points = ax.scatter(
    high["UMAP1"], high["UMAP2"], c=high["OS_composite"],
    cmap=OS_CMAP, norm=Normalize(OS_HIGH_CUTOFF, OS_HIGH_MAX),
    s=1.15, alpha=0.94, linewidths=0, rasterized=True, zorder=2,
)

ax.set_aspect("equal")
ax.set_anchor("C")
ax.set_xticks([])
ax.set_yticks([])
for spine in ax.spines.values():
    spine.set_visible(False)
ax.margins(0.025)
ax.text(0, 1.015, "OS-high malignant cells", transform=ax.transAxes,
        ha="left", va="bottom", fontsize=7.5, color="#222222")

cbar = fig.colorbar(points, ax=ax, fraction=0.047, pad=0.025, shrink=0.72)
cbar.set_label("OS composite", labelpad=2)
cbar.set_ticks([0.8, 1.4, 2.0, 2.7])
cbar.ax.tick_params(length=2.2, width=0.65, pad=1.5)
cbar.outline.set_linewidth(0.65)

fig.text(0.015, 0.982, "B", ha="left", va="top", fontsize=12.5,
         fontweight="bold", fontfamily="Times New Roman", color="#111111")

stem = OUT / "Panel_B_OS_high_review_v2"
for suffix, dpi in (("pdf", None), ("svg", None), ("png", 600)):
    kwargs = {"facecolor": "white", "bbox_inches": None, "pad_inches": 0}
    if dpi is not None:
        kwargs["dpi"] = dpi
    fig.savefig(stem.with_suffix(f".{suffix}"), **kwargs)
plt.close(fig)

print("PANEL_B_REVIEW_V2_RENDERED")
