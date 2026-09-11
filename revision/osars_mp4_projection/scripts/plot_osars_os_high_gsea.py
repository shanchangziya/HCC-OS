#!/usr/bin/env python3
"""Render standalone, publication-ready GSEA panels from processed results.

This script is intentionally plot-only: it reads the processed GSEA result and
running-score tables created by run_osars_os_high_gsea.py and does not compute
or modify any statistics.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont
from reportlab.lib.colors import HexColor
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfgen import canvas


ROOT = Path(__file__).resolve().parents[1]
TABLE = ROOT / "tables" / "GSEA_results.csv"
COORDINATES = ROOT / "data" / "processed" / "GSEA_running_score_coordinates.csv.gz"
OUTPUT = ROOT / "figures"

ARIAL = Path("/System/Library/Fonts/Supplemental/Arial.ttf")
ARIAL_BOLD = Path("/System/Library/Fonts/Supplemental/Arial Bold.ttf")

SETS = [
    "Frozen discovery MP4 state",
    "Frozen discovery MP4 state excluding ROS and OSARS features",
    "Hallmark ROS positive control",
]
SET_LABELS = {
    "Frozen discovery MP4 state": "Frozen MP4 state",
    "Frozen discovery MP4 state excluding ROS and OSARS features": (
        "MP4 state minus ROS/OSARS genes"
    ),
    "Hallmark ROS positive control": "Hallmark ROS",
}
COLORS = {
    "high": "#D55E00",
    "low": "#0072B2",
    "axis": "#333333",
    "gray": "#707070",
    "light_gray": "#C7C7C7",
    "very_light": "#ECECEC",
    "white": "#FFFFFF",
}


def q_label(value: float) -> str:
    if value < 0.001:
        return "BH q < 0.001"
    return f"BH q = {value:.3f}"


def pdf_fonts() -> tuple[str, str]:
    regular, bold = "Helvetica", "Helvetica-Bold"
    if ARIAL.exists() and ARIAL_BOLD.exists():
        regular, bold = "ArialEmbedded", "ArialEmbedded-Bold"
        if regular not in pdfmetrics.getRegisteredFontNames():
            pdfmetrics.registerFont(TTFont(regular, str(ARIAL)))
        if bold not in pdfmetrics.getRegisteredFontNames():
            pdfmetrics.registerFont(TTFont(bold, str(ARIAL_BOLD)))
    return regular, bold


def pil_font(path: Path, size: float) -> ImageFont.FreeTypeFont:
    if path.exists():
        return ImageFont.truetype(str(path), max(1, round(size)))
    return ImageFont.load_default()


def draw_summary_pdf(results: pd.DataFrame, analysis: str, path: Path) -> None:
    """Draw a two-cohort NES summary with fixed scale and BH labels."""
    width, height = 510.0, 228.0
    label_right = 123.0
    panels = {"TCGA-LIHC": (142.0, 308.0), "ICGC-LIRI": (332.0, 498.0)}
    x_min, x_max = -4.0, 2.0
    row_y = [77.0, 118.0, 159.0]
    axis_y = 184.0
    document = canvas.Canvas(str(path), pagesize=(width, height), pageCompression=1)
    document.setTitle(f"OSARS MP4 GSEA {analysis}")
    regular, bold = pdf_fonts()

    def py(y: float) -> float:
        return height - y

    def px(value: float, left: float, right: float) -> float:
        return left + (value - x_min) / (x_max - x_min) * (right - left)

    document.setFillColor(HexColor(COLORS["axis"]))
    document.setFont(bold, 9)
    document.drawString(142.0, py(17.0), "GSEA: OSARS-high versus OSARS-low")
    document.setFont(regular, 7)
    document.setFillColor(HexColor(COLORS["gray"]))
    descriptor = "Primary analysis" if analysis == "primary" else "Stage-adjusted sensitivity"
    document.drawString(142.0, py(30.0), descriptor)

    for position, gene_set in enumerate(SETS):
        document.setFillColor(HexColor(COLORS["axis"]))
        document.setFont(regular, 7.3)
        document.drawRightString(label_right, py(row_y[position] - 2.5), SET_LABELS[gene_set])

    for cohort, (left, right) in panels.items():
        subset = results.loc[
            (results["cohort"] == cohort) & (results["analysis"] == analysis)
        ].set_index("gene_set").loc[SETS].reset_index()
        n_text = f"n = {int(subset['n_total'].iloc[0])}"
        document.setFillColor(HexColor(COLORS["axis"]))
        document.setFont(bold, 8)
        document.drawCentredString((left + right) / 2, py(48.0), cohort)
        document.setFillColor(HexColor(COLORS["gray"]))
        document.setFont(regular, 6.5)
        document.drawCentredString((left + right) / 2, py(59.0), n_text)
        zero = px(0.0, left, right)
        document.setStrokeColor(HexColor(COLORS["light_gray"]))
        document.setLineWidth(0.7)
        document.setDash(2.0, 2.0)
        document.line(zero, py(64.0), zero, py(axis_y))
        document.setDash()
        for position, (_, row) in enumerate(subset.iterrows()):
            value = float(row["NES"])
            color = COLORS["high"] if value > 0 else COLORS["low"]
            x = px(value, left, right)
            y = py(row_y[position])
            significant = float(row["BH_q_within_analysis_6tests"]) < 0.05
            document.setFillColor(HexColor(color if significant else COLORS["white"]))
            document.setStrokeColor(HexColor(color))
            document.setLineWidth(1.25)
            document.circle(x, y, 3.6, fill=1, stroke=1)
            document.setFillColor(HexColor(color))
            document.setFont(regular, 6.1)
            text_x = x + 6.2 if value < 1.2 else x - 6.2
            if value < 1.2:
                document.drawString(text_x, y - 2.2, q_label(float(row["BH_q_within_analysis_6tests"])))
            else:
                document.drawRightString(text_x, y - 2.2, q_label(float(row["BH_q_within_analysis_6tests"])))

        document.setStrokeColor(HexColor(COLORS["axis"]))
        document.setLineWidth(0.75)
        document.line(left, py(axis_y), right, py(axis_y))
        document.setFont(regular, 6.4)
        document.setFillColor(HexColor(COLORS["axis"]))
        for tick in (-4, -2, 0, 2):
            x = px(tick, left, right)
            document.line(x, py(axis_y), x, py(axis_y) - 3.0)
            document.drawCentredString(x, py(axis_y + 13), str(tick))
        document.setFont(regular, 7)
        document.drawCentredString((left + right) / 2, py(216.0), "Normalized enrichment score")
    document.save()


def draw_summary_png(results: pd.DataFrame, analysis: str, path: Path) -> None:
    """Render the same NES summary as a 600 dpi PNG preview."""
    width, height, scale = 510.0, 228.0, 4
    label_right = 123.0
    panels = {"TCGA-LIHC": (142.0, 308.0), "ICGC-LIRI": (332.0, 498.0)}
    x_min, x_max = -4.0, 2.0
    row_y = [77.0, 118.0, 159.0]
    axis_y = 184.0
    image = Image.new("RGB", (round(width * scale), round(height * scale)), COLORS["white"])
    draw = ImageDraw.Draw(image)
    regular = pil_font(ARIAL, 7.3 * scale)
    bold = pil_font(ARIAL_BOLD, 8 * scale)
    title = pil_font(ARIAL_BOLD, 9 * scale)
    small = pil_font(ARIAL, 6.4 * scale)
    tiny = pil_font(ARIAL, 6.1 * scale)

    def pt(x: float, y: float) -> tuple[int, int]:
        return round(x * scale), round(y * scale)

    def line(x1: float, y1: float, x2: float, y2: float, color: str, line_width: float) -> None:
        draw.line([pt(x1, y1), pt(x2, y2)], fill=color, width=max(1, round(line_width * scale)))

    def right_text(x: float, y: float, text: str, font: ImageFont.FreeTypeFont, color: str) -> None:
        box = draw.textbbox((0, 0), text, font=font)
        draw.text((round(x * scale) - (box[2] - box[0]), round(y * scale)), text, font=font, fill=color)

    def centered(x: float, y: float, text: str, font: ImageFont.FreeTypeFont, color: str) -> None:
        box = draw.textbbox((0, 0), text, font=font)
        draw.text((round(x * scale) - (box[2] - box[0]) / 2, round(y * scale)), text, font=font, fill=color)

    def px(value: float, left: float, right: float) -> float:
        return left + (value - x_min) / (x_max - x_min) * (right - left)

    draw.text(pt(142.0, 7.0), "GSEA: OSARS-high versus OSARS-low", font=title, fill=COLORS["axis"])
    descriptor = "Primary analysis" if analysis == "primary" else "Stage-adjusted sensitivity"
    draw.text(pt(142.0, 23.0), descriptor, font=small, fill=COLORS["gray"])
    for position, gene_set in enumerate(SETS):
        right_text(label_right, row_y[position] - 6.5, SET_LABELS[gene_set], regular, COLORS["axis"])

    for cohort, (left, right) in panels.items():
        subset = results.loc[
            (results["cohort"] == cohort) & (results["analysis"] == analysis)
        ].set_index("gene_set").loc[SETS].reset_index()
        centered((left + right) / 2, 40.0, cohort, bold, COLORS["axis"])
        centered((left + right) / 2, 53.0, f"n = {int(subset['n_total'].iloc[0])}", small, COLORS["gray"])
        zero = px(0.0, left, right)
        for y in np.arange(64, axis_y, 4.0):
            line(zero, float(y), zero, float(min(y + 2.0, axis_y)), COLORS["light_gray"], 0.7)
        for position, (_, row) in enumerate(subset.iterrows()):
            value = float(row["NES"])
            color = COLORS["high"] if value > 0 else COLORS["low"]
            x = px(value, left, right)
            y = row_y[position]
            significant = float(row["BH_q_within_analysis_6tests"]) < 0.05
            radius = 3.6 * scale
            px_x, px_y = pt(x, y)
            draw.ellipse(
                (px_x - radius, px_y - radius, px_x + radius, px_y + radius),
                fill=color if significant else COLORS["white"],
                outline=color,
                width=max(1, round(1.25 * scale)),
            )
            label = q_label(float(row["BH_q_within_analysis_6tests"]))
            box = draw.textbbox((0, 0), label, font=tiny)
            if value < 1.2:
                draw.text((round((x + 6.2) * scale), round((y - 4.5) * scale)), label, font=tiny, fill=color)
            else:
                draw.text((round((x - 6.2) * scale) - (box[2] - box[0]), round((y - 4.5) * scale)), label, font=tiny, fill=color)
        line(left, axis_y, right, axis_y, COLORS["axis"], 0.75)
        for tick in (-4, -2, 0, 2):
            x = px(tick, left, right)
            line(x, axis_y, x, axis_y + 3, COLORS["axis"], 0.75)
            centered(x, axis_y + 6, str(tick), small, COLORS["axis"])
        centered((left + right) / 2, 207.0, "Normalized enrichment score", regular, COLORS["axis"])
    image.save(path, dpi=(600, 600))


def draw_running_pdf(
    coordinates: pd.DataFrame, result: pd.Series, cohort: str, path: Path
) -> None:
    """Draw one GSEA running-score plot for the non-overlapping MP4 sensitivity."""
    width, height = 260.0, 246.0
    left, right, top, bottom = 48.0, 244.0, 62.0, 182.0
    subtitle = (
        f"MP4 state minus ROS/OSARS genes; NES = {float(result['NES']):.2f}; "
        f"{q_label(float(result['BH_q_within_analysis_6tests']))}"
    )
    max_abs = max(abs(float(coordinates["running_ES"].min())), abs(float(coordinates["running_ES"].max())), 0.1)
    y_min, y_max = -max_abs * 1.12, max_abs * 1.12
    regular, bold = pdf_fonts()
    document = canvas.Canvas(str(path), pagesize=(width, height), pageCompression=1)
    document.setTitle(f"{cohort} MP4-minus-ROS-and-OSARS GSEA")

    def py(y: float) -> float:
        return height - y

    def px(rank_fraction: float) -> float:
        return left + rank_fraction * (right - left)

    def yv(value: float) -> float:
        return bottom - (value - y_min) / (y_max - y_min) * (bottom - top)

    document.setFillColor(HexColor(COLORS["axis"]))
    document.setFont(bold, 9)
    document.drawString(left, py(17.0), cohort)
    document.setFont(regular, 6.5)
    document.setFillColor(HexColor(COLORS["gray"]))
    document.drawString(left, py(29.0), subtitle)
    document.drawString(left, py(40.0), f"n = {int(result['n_total'])}; high = {int(result['n_high'])}; low = {int(result['n_low'])}")

    zero_y = yv(0.0)
    document.setStrokeColor(HexColor(COLORS["light_gray"]))
    document.setLineWidth(0.7)
    document.line(left, py(zero_y), right, py(zero_y))
    document.setStrokeColor(HexColor(COLORS["low"]))
    document.setLineWidth(1.35)
    coordinates = coordinates.sort_values("rank")
    stride = max(1, len(coordinates) // 2000)
    sampled = coordinates.iloc[::stride]
    path_object = document.beginPath()
    first = sampled.iloc[0]
    path_object.moveTo(px(float(first["rank_fraction"])), py(yv(float(first["running_ES"]))))
    for _, row in sampled.iloc[1:].iterrows():
        path_object.lineTo(px(float(row["rank_fraction"])), py(yv(float(row["running_ES"]))))
    last = coordinates.iloc[-1]
    path_object.lineTo(px(float(last["rank_fraction"])), py(yv(float(last["running_ES"]))))
    document.drawPath(path_object, stroke=1, fill=0)

    document.setStrokeColor(HexColor(COLORS["axis"]))
    document.setLineWidth(0.55)
    rug_top, rug_bottom = 190.0, 201.0
    for fraction in coordinates.loc[coordinates["is_gene_set_hit"], "rank_fraction"]:
        x = px(float(fraction))
        document.line(x, py(rug_top), x, py(rug_bottom))
    document.setStrokeColor(HexColor(COLORS["axis"]))
    document.setLineWidth(0.75)
    document.line(left, py(bottom), right, py(bottom))
    document.setFont(regular, 6.5)
    document.setFillColor(HexColor(COLORS["axis"]))
    document.drawCentredString(left, py(218.0), "OSARS-high")
    document.drawCentredString(right, py(218.0), "OSARS-low")
    document.drawCentredString((left + right) / 2, py(235.0), "Genes ranked by high minus low expression")
    document.save()


def draw_running_png(
    coordinates: pd.DataFrame, result: pd.Series, cohort: str, path: Path
) -> None:
    width, height, scale = 260.0, 246.0, 4
    left, right, top, bottom = 48.0, 244.0, 62.0, 182.0
    subtitle = (
        f"MP4 state minus ROS/OSARS genes; NES = {float(result['NES']):.2f}; "
        f"{q_label(float(result['BH_q_within_analysis_6tests']))}"
    )
    max_abs = max(abs(float(coordinates["running_ES"].min())), abs(float(coordinates["running_ES"].max())), 0.1)
    y_min, y_max = -max_abs * 1.12, max_abs * 1.12
    image = Image.new("RGB", (round(width * scale), round(height * scale)), COLORS["white"])
    draw = ImageDraw.Draw(image)
    regular = pil_font(ARIAL, 6.5 * scale)
    bold = pil_font(ARIAL_BOLD, 9 * scale)

    def pt(x: float, y: float) -> tuple[int, int]:
        return round(x * scale), round(y * scale)

    def line(x1: float, y1: float, x2: float, y2: float, color: str, line_width: float) -> None:
        draw.line([pt(x1, y1), pt(x2, y2)], fill=color, width=max(1, round(line_width * scale)))

    def centered(x: float, y: float, text: str, font: ImageFont.FreeTypeFont, color: str) -> None:
        box = draw.textbbox((0, 0), text, font=font)
        draw.text((round(x * scale) - (box[2] - box[0]) / 2, round(y * scale)), text, font=font, fill=color)

    def px(rank_fraction: float) -> float:
        return left + rank_fraction * (right - left)

    def yv(value: float) -> float:
        return bottom - (value - y_min) / (y_max - y_min) * (bottom - top)

    draw.text(pt(left, 8.0), cohort, font=bold, fill=COLORS["axis"])
    draw.text(pt(left, 24.0), subtitle, font=regular, fill=COLORS["gray"])
    draw.text(pt(left, 36.0), f"n = {int(result['n_total'])}; high = {int(result['n_high'])}; low = {int(result['n_low'])}", font=regular, fill=COLORS["gray"])
    line(left, yv(0), right, yv(0), COLORS["light_gray"], 0.7)

    coordinates = coordinates.sort_values("rank")
    stride = max(1, len(coordinates) // 2000)
    sampled = coordinates.iloc[::stride]
    points = [pt(px(float(row.rank_fraction)), yv(float(row.running_ES))) for row in sampled.itertuples()]
    points.append(pt(px(float(coordinates.iloc[-1]["rank_fraction"])), yv(float(coordinates.iloc[-1]["running_ES"]))))
    draw.line(points, fill=COLORS["low"], width=max(1, round(1.35 * scale)), joint="curve")
    for fraction in coordinates.loc[coordinates["is_gene_set_hit"], "rank_fraction"]:
        x = px(float(fraction))
        line(x, 190.0, x, 201.0, COLORS["axis"], 0.55)
    line(left, bottom, right, bottom, COLORS["axis"], 0.75)
    centered(left, 211.0, "OSARS-high", regular, COLORS["axis"])
    centered(right, 211.0, "OSARS-low", regular, COLORS["axis"])
    centered((left + right) / 2, 228.0, "Genes ranked by high minus low expression", regular, COLORS["axis"])
    image.save(path, dpi=(600, 600))


def main() -> None:
    if not TABLE.exists() or not COORDINATES.exists():
        raise FileNotFoundError("Run run_osars_os_high_gsea.py before plotting")
    OUTPUT.mkdir(parents=True, exist_ok=True)
    results = pd.read_csv(TABLE)
    coordinates = pd.read_csv(COORDINATES)
    expected = {"TCGA-LIHC", "ICGC-LIRI"}
    if set(results["cohort"]) != expected or set(results["analysis"]) != {"primary", "stage_adjusted"}:
        raise ValueError("Unexpected GSEA result table structure")
    manifest: dict[str, object] = {
        "input_result_table": str(TABLE.relative_to(ROOT)),
        "input_running_coordinates": str(COORDINATES.relative_to(ROOT)),
        "figures": {},
        "colors": {
            "OSARS_high_enrichment": COLORS["high"],
            "OSARS_low_enrichment": COLORS["low"],
        },
        "plot_definition": "Filled point: BH q < 0.05. Open point: BH q >= 0.05.",
    }
    for analysis, basename in [
        ("primary", "R1Q9_primary_GSEA_NES_summary"),
        ("stage_adjusted", "R1Q9_stage_adjusted_GSEA_NES_summary"),
    ]:
        pdf = OUTPUT / f"{basename}.pdf"
        png = OUTPUT / f"{basename}.png"
        draw_summary_pdf(results, analysis, pdf)
        draw_summary_png(results, analysis, png)
        manifest["figures"][basename] = {
            "pdf": pdf.name,
            "png": png.name,
            "analysis": analysis,
            "data": "tables/GSEA_results.csv",
        }
    state_set = "Frozen discovery MP4 state excluding ROS and OSARS features"
    for cohort, basename in [
        ("TCGA-LIHC", "R1Q9_TCGA_MP4minusROS_GSEA_running"),
        ("ICGC-LIRI", "R1Q9_ICGC_MP4minusROS_GSEA_running"),
    ]:
        result = results.loc[
            (results["cohort"] == cohort)
            & (results["analysis"] == "primary")
            & (results["gene_set"] == state_set)
        ]
        if len(result) != 1:
            raise ValueError(f"Missing primary running-score result for {cohort}")
        result_row = result.iloc[0]
        curve = coordinates.loc[
            (coordinates["cohort"] == cohort)
            & (coordinates["analysis"] == "primary")
            & (coordinates["gene_set"] == state_set)
        ].copy()
        pdf = OUTPUT / f"{basename}.pdf"
        png = OUTPUT / f"{basename}.png"
        draw_running_pdf(curve, result_row, cohort, pdf)
        draw_running_png(curve, result_row, cohort, png)
        manifest["figures"][basename] = {
            "pdf": pdf.name,
            "png": png.name,
            "analysis": "primary",
            "gene_set": state_set,
            "data": "data/processed/GSEA_running_score_coordinates.csv.gz",
        }
    (OUTPUT / "figure_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print("Rendered GSEA figures")
    for name, value in manifest["figures"].items():
        print(name, value)


if __name__ == "__main__":
    main()
