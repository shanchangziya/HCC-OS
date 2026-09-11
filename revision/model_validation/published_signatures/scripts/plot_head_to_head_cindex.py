#!/usr/bin/env python3
"""Render separate TCGA and ICGC C-index forest plots from fixed result data.

No statistical quantity is recalculated here. The script reads the primary
complete-formula rows from tables/head_to_head_cindex.csv and renders exactly
two standalone panels. OSARS is the only highlighted marker; published
comparators are identified by PMID rather than article names.
"""

from __future__ import annotations

import hashlib
import html
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
INPUT = ROOT / "tables" / "head_to_head_cindex.csv"
OUTPUT = ROOT / "figures"
PLOT_DATA = OUTPUT / "head_to_head_cindex_plot_data.csv"

ARIAL = Path("/System/Library/Fonts/Supplemental/Arial.ttf")
ARIAL_BOLD = Path("/System/Library/Fonts/Supplemental/Arial Bold.ttf")

ORDER = [
    "OSARS_frozen",
    "Liu2021_FAIS_signature1",
    "FuSong2021_pyroptosis_3gene",
    "Lin2021_inflammatory_8gene",
    "HongCai2022_OSRG_8gene",
    "Ma2024_OS_ER_stress_5gene",
]
LABELS = {
    "OSARS_frozen": "OSARS",
    "Liu2021_FAIS_signature1": "PMID 34970594",
    "FuSong2021_pyroptosis_3gene": "PMID 34820376",
    "Lin2021_inflammatory_8gene": "PMID 33828988",
    "HongCai2022_OSRG_8gene": "PMID 36133439",
    "Ma2024_OS_ER_stress_5gene": "PMID 37957902",
}

COLORS = {
    "osars": "#D55E00",
    "comparator": "#5F5F5F",
    "ci": "#8C8C8C",
    "axis": "#333333",
    "baseline": "#BDBDBD",
    "subtle": "#666666",
    "white": "#FFFFFF",
}

# Dimensions are in PDF/SVG points (1/72 inch), 4.4 x 3.47 inches. This
# allows readable 7-9 pt labels without crowding the six forest rows.
WIDTH = 316.8
HEIGHT = 249.8
LEFT = 108.0
RIGHT = 301.0
TITLE_Y = 17.0
SUBTITLE_Y = 30.0
FIRST_ROW_Y = 58.0
ROW_STEP = 25.0
AXIS_Y = 207.0
TICK_LABEL_Y = 220.0
XLABEL_Y = 241.0
XMIN = 0.40
XMAX = 0.95
XTICKS = (0.4, 0.5, 0.6, 0.7, 0.8, 0.9)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def x_coordinate(value: float) -> float:
    return LEFT + (float(value) - XMIN) / (XMAX - XMIN) * (RIGHT - LEFT)


def layout_rows(frame: pd.DataFrame) -> list[dict]:
    """Produce ordered drawing records from already-calculated C-index rows."""
    records = []
    for position, predictor in enumerate(ORDER):
        candidate = frame.loc[frame["predictor"] == predictor]
        if len(candidate) != 1:
            raise ValueError(f"Expected one {predictor} row, found {len(candidate)}")
        row = candidate.iloc[0]
        records.append(
            {
                "predictor": predictor,
                "label": LABELS[predictor],
                "cindex": float(row["cindex"]),
                "ci_lower": float(row["ci_lower"]),
                "ci_upper": float(row["ci_upper"]),
                "n": int(row["n"]),
                "events": int(row["events"]),
                "row_y": FIRST_ROW_Y + position * ROW_STEP,
                "is_osars": predictor == "OSARS_frozen",
            }
        )
    return records


def register_pdf_fonts() -> tuple[str, str]:
    """Embed the system Arial faces where the PDF runtime permits it."""
    regular = "Helvetica"
    bold = "Helvetica-Bold"
    if ARIAL.exists() and ARIAL_BOLD.exists():
        regular = "ArialEmbedded"
        bold = "ArialEmbedded-Bold"
        if regular not in pdfmetrics.getRegisteredFontNames():
            pdfmetrics.registerFont(TTFont(regular, str(ARIAL)))
        if bold not in pdfmetrics.getRegisteredFontNames():
            pdfmetrics.registerFont(TTFont(bold, str(ARIAL_BOLD)))
    return regular, bold


def draw_pdf(cohort: str, rows: list[dict], path: Path) -> None:
    regular, bold = register_pdf_fonts()
    document = canvas.Canvas(str(path), pagesize=(WIDTH, HEIGHT), pageCompression=1)
    document.setTitle(f"{cohort} C-index forest plot")
    document.setAuthor("OSARS revision analysis")

    def py(y: float) -> float:
        return HEIGHT - y

    document.setFillColor(HexColor(COLORS["axis"]))
    document.setFont(bold, 9)
    document.drawString(LEFT, py(TITLE_Y), cohort)
    document.setFillColor(HexColor(COLORS["subtle"]))
    document.setFont(regular, 6.8)
    document.drawString(LEFT, py(SUBTITLE_Y), f"n = {rows[0]['n']}; OS events = {rows[0]['events']}")

    baseline_x = x_coordinate(0.5)
    document.setStrokeColor(HexColor(COLORS["baseline"]))
    document.setLineWidth(0.65)
    document.setDash(2.2, 2.2)
    document.line(baseline_x, py(42), baseline_x, py(AXIS_Y))
    document.setDash()

    for row in rows:
        y = py(row["row_y"])
        accent = COLORS["osars"] if row["is_osars"] else COLORS["comparator"]
        ci_color = COLORS["osars"] if row["is_osars"] else COLORS["ci"]
        x_low, x_mid, x_high = (
            x_coordinate(row["ci_lower"]),
            x_coordinate(row["cindex"]),
            x_coordinate(row["ci_upper"]),
        )
        document.setStrokeColor(HexColor(ci_color))
        document.setLineWidth(1.35 if row["is_osars"] else 1.05)
        document.line(x_low, y, x_high, y)
        document.setLineWidth(0.95)
        document.line(x_low, y - 3.6, x_low, y + 3.6)
        document.line(x_high, y - 3.6, x_high, y + 3.6)
        document.setFillColor(HexColor(accent))
        document.setStrokeColor(HexColor(COLORS["white"]))
        document.setLineWidth(0.75)
        radius = 4.1 if row["is_osars"] else 3.1
        document.circle(x_mid, y, radius, fill=1, stroke=1)
        document.setFillColor(HexColor(accent))
        document.setFont(bold if row["is_osars"] else regular, 7.6)
        document.drawRightString(LEFT - 9, y - 2.6, row["label"])

    document.setStrokeColor(HexColor(COLORS["axis"]))
    document.setLineWidth(0.75)
    document.line(LEFT, py(AXIS_Y), RIGHT, py(AXIS_Y))
    document.setFont(regular, 6.8)
    document.setFillColor(HexColor(COLORS["axis"]))
    for tick in XTICKS:
        x = x_coordinate(tick)
        document.line(x, py(AXIS_Y), x, py(AXIS_Y) - 3.0)
        document.drawCentredString(x, py(TICK_LABEL_Y), f"{tick:.1f}")
    document.setFont(regular, 8)
    document.drawCentredString((LEFT + RIGHT) / 2, py(XLABEL_Y), "Harrell C-index (95% CI)")
    document.save()


def draw_svg(cohort: str, rows: list[dict], path: Path) -> None:
    """Write an editable, font-preserving SVG version of the same panel."""
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH / 72:.3f}in" height="{HEIGHT / 72:.3f}in" viewBox="0 0 {WIDTH} {HEIGHT}">',
        f'<rect width="{WIDTH}" height="{HEIGHT}" fill="{COLORS["white"]}"/>',
        '<style>text { font-family: Arial, Helvetica, sans-serif; } .axis { fill: #333333; }</style>',
        f'<text x="{LEFT}" y="{TITLE_Y}" font-size="9" font-weight="700" fill="{COLORS["axis"]}">{html.escape(cohort)}</text>',
        f'<text x="{LEFT}" y="{SUBTITLE_Y}" font-size="6.8" fill="{COLORS["subtle"]}">n = {rows[0]["n"]}; OS events = {rows[0]["events"]}</text>',
    ]
    baseline_x = x_coordinate(0.5)
    parts.append(
        f'<line x1="{baseline_x:.3f}" y1="42" x2="{baseline_x:.3f}" y2="{AXIS_Y}" '
        f'stroke="{COLORS["baseline"]}" stroke-width="0.65" stroke-dasharray="2.2 2.2"/>'
    )
    for row in rows:
        accent = COLORS["osars"] if row["is_osars"] else COLORS["comparator"]
        ci_color = COLORS["osars"] if row["is_osars"] else COLORS["ci"]
        x_low, x_mid, x_high = x_coordinate(row["ci_lower"]), x_coordinate(row["cindex"]), x_coordinate(row["ci_upper"])
        y = row["row_y"]
        line_width = 1.35 if row["is_osars"] else 1.05
        radius = 4.1 if row["is_osars"] else 3.1
        weight = "700" if row["is_osars"] else "400"
        parts += [
            f'<line x1="{x_low:.3f}" y1="{y:.3f}" x2="{x_high:.3f}" y2="{y:.3f}" stroke="{ci_color}" stroke-width="{line_width}"/>',
            f'<line x1="{x_low:.3f}" y1="{y - 3.6:.3f}" x2="{x_low:.3f}" y2="{y + 3.6:.3f}" stroke="{ci_color}" stroke-width="0.95"/>',
            f'<line x1="{x_high:.3f}" y1="{y - 3.6:.3f}" x2="{x_high:.3f}" y2="{y + 3.6:.3f}" stroke="{ci_color}" stroke-width="0.95"/>',
            f'<circle cx="{x_mid:.3f}" cy="{y:.3f}" r="{radius}" fill="{accent}" stroke="{COLORS["white"]}" stroke-width="0.75"/>',
            f'<text x="{LEFT - 9}" y="{y + 2.6:.3f}" text-anchor="end" font-size="7.6" font-weight="{weight}" fill="{accent}">{html.escape(row["label"])}</text>',
        ]
    parts.append(
        f'<line x1="{LEFT}" y1="{AXIS_Y}" x2="{RIGHT}" y2="{AXIS_Y}" stroke="{COLORS["axis"]}" stroke-width="0.75"/>'
    )
    for tick in XTICKS:
        x = x_coordinate(tick)
        parts += [
            f'<line x1="{x:.3f}" y1="{AXIS_Y}" x2="{x:.3f}" y2="{AXIS_Y + 3}" stroke="{COLORS["axis"]}" stroke-width="0.75"/>',
            f'<text x="{x:.3f}" y="{TICK_LABEL_Y}" text-anchor="middle" font-size="6.8" fill="{COLORS["axis"]}">{tick:.1f}</text>',
        ]
    parts += [
        f'<text x="{(LEFT + RIGHT) / 2:.3f}" y="{XLABEL_Y}" text-anchor="middle" font-size="8" fill="{COLORS["axis"]}">Harrell C-index (95% CI)</text>',
        "</svg>",
    ]
    path.write_text("\n".join(parts) + "\n")


def _font(path: Path, size: float, fallback: bool = False) -> ImageFont.FreeTypeFont:
    chosen = ARIAL_BOLD if fallback else path
    if not chosen.exists():
        return ImageFont.load_default()
    return ImageFont.truetype(str(chosen), int(round(size)))


def draw_png(cohort: str, rows: list[dict], path: Path) -> None:
    """Render a high-resolution preview in the same physical layout."""
    scale = 4
    image = Image.new("RGB", (round(WIDTH * scale), round(HEIGHT * scale)), COLORS["white"])
    draw = ImageDraw.Draw(image)
    regular = _font(ARIAL, 7.6 * scale)
    bold = _font(ARIAL_BOLD, 7.6 * scale)
    title_font = _font(ARIAL_BOLD, 9 * scale)
    subtitle_font = _font(ARIAL, 6.8 * scale)
    tick_font = _font(ARIAL, 6.8 * scale)
    xlabel_font = _font(ARIAL, 8 * scale)

    def point(x: float, y: float) -> tuple[int, int]:
        return round(x * scale), round(y * scale)

    def line(x1: float, y1: float, x2: float, y2: float, fill: str, width: float) -> None:
        draw.line([point(x1, y1), point(x2, y2)], fill=fill, width=max(1, round(width * scale)))

    def right_text(x: float, y: float, value: str, font: ImageFont.FreeTypeFont, fill: str) -> None:
        bbox = draw.textbbox((0, 0), value, font=font)
        draw.text((round(x * scale) - (bbox[2] - bbox[0]), round(y * scale)), value, font=font, fill=fill)

    draw.text(point(LEFT, TITLE_Y - 9), cohort, font=title_font, fill=COLORS["axis"])
    draw.text(point(LEFT, SUBTITLE_Y - 7), f"n = {rows[0]['n']}; OS events = {rows[0]['events']}", font=subtitle_font, fill=COLORS["subtle"])

    baseline_x = x_coordinate(0.5)
    for y in np.arange(42, AXIS_Y, 4.4):
        line(baseline_x, float(y), baseline_x, float(min(y + 2.2, AXIS_Y)), COLORS["baseline"], 0.65)

    for row in rows:
        y = row["row_y"]
        accent = COLORS["osars"] if row["is_osars"] else COLORS["comparator"]
        ci_color = COLORS["osars"] if row["is_osars"] else COLORS["ci"]
        x_low, x_mid, x_high = x_coordinate(row["ci_lower"]), x_coordinate(row["cindex"]), x_coordinate(row["ci_upper"])
        line(x_low, y, x_high, y, ci_color, 1.35 if row["is_osars"] else 1.05)
        line(x_low, y - 3.6, x_low, y + 3.6, ci_color, 0.95)
        line(x_high, y - 3.6, x_high, y + 3.6, ci_color, 0.95)
        radius = (4.1 if row["is_osars"] else 3.1) * scale
        px, py = point(x_mid, y)
        draw.ellipse((px - radius, py - radius, px + radius, py + radius), fill=accent, outline=COLORS["white"], width=max(1, round(0.75 * scale)))
        label_font = bold if row["is_osars"] else regular
        right_text(LEFT - 9, y - 5.1, row["label"], label_font, accent)

    line(LEFT, AXIS_Y, RIGHT, AXIS_Y, COLORS["axis"], 0.75)
    for tick in XTICKS:
        x = x_coordinate(tick)
        line(x, AXIS_Y, x, AXIS_Y + 3, COLORS["axis"], 0.75)
        label = f"{tick:.1f}"
        bbox = draw.textbbox((0, 0), label, font=tick_font)
        draw.text((round(x * scale) - (bbox[2] - bbox[0]) / 2, round((TICK_LABEL_Y - 6) * scale)), label, font=tick_font, fill=COLORS["axis"])
    xlabel = "Harrell C-index (95% CI)"
    bbox = draw.textbbox((0, 0), xlabel, font=xlabel_font)
    draw.text((round((LEFT + RIGHT) / 2 * scale) - (bbox[2] - bbox[0]) / 2, round((XLABEL_Y - 7) * scale)), xlabel, font=xlabel_font, fill=COLORS["axis"])
    image.save(path, dpi=(600, 600))


def main() -> None:
    if not INPUT.exists():
        raise FileNotFoundError(f"Missing frozen C-index result table: {INPUT}")
    OUTPUT.mkdir(parents=True, exist_ok=True)
    result = pd.read_csv(INPUT)
    primary = result.loc[result["analysis_set"] == "primary_complete_formula"].copy()
    if set(primary["predictor"]) != set(ORDER):
        raise ValueError("Primary C-index table does not contain the expected six predictors")

    plot_rows = []
    manifest = {
        "input": str(INPUT.relative_to(ROOT)),
        "input_sha256": sha256(INPUT),
        "analysis_set": "primary_complete_formula",
        "scope": "OSARS plus published signatures with complete reported numerical formulas",
        "omitted_from_plot": {
            "Tian2021_five_gene": "numerical multivariable coefficients not reported",
            "Wang2022_cuproptosis_5gene": "not an exact rerun; rounded-HR proxy is sensitivity-only",
        },
        "pmid_labels": {LABELS[predictor]: predictor for predictor in ORDER if predictor != "OSARS_frozen"},
        "outputs": {},
    }
    for cohort, shortname in [("TCGA-LIHC", "TCGA_LIHC"), ("ICGC-LIRI", "ICGC_LIRI")]:
        records = layout_rows(primary.loc[primary["cohort"] == cohort])
        for record in records:
            plot_rows.append({"cohort": cohort, **{key: value for key, value in record.items() if key != "row_y"}})
        pdf_path = OUTPUT / f"{shortname}_head_to_head_Cindex.pdf"
        svg_path = OUTPUT / f"{shortname}_head_to_head_Cindex.svg"
        png_path = OUTPUT / f"{shortname}_head_to_head_Cindex.png"
        draw_pdf(cohort, records, pdf_path)
        draw_svg(cohort, records, svg_path)
        draw_png(cohort, records, png_path)
        manifest["outputs"][cohort] = {
            "pdf": pdf_path.name,
            "svg": svg_path.name,
            "png": png_path.name,
            "n": records[0]["n"],
            "events": records[0]["events"],
        }
    pd.DataFrame(plot_rows).to_csv(PLOT_DATA, index=False)
    (OUTPUT / "head_to_head_Cindex_figure_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print("Rendered two C-index forest plots")
    for cohort, details in manifest["outputs"].items():
        print(cohort, details)


if __name__ == "__main__":
    main()
