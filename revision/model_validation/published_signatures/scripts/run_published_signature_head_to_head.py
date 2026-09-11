#!/usr/bin/env python3
"""Recalculate published HCC signature scores on the frozen OSARS cohorts.

This isolated script consumes frozen inputs only and never refits a published
model. The primary comparison includes only articles that report numerical Cox
formulas. Tian 2021 does not report its fitted coefficients. Wang 2022 reports
rounded multivariable hazard ratios, but not fitted coefficients; a clearly
labelled log(rounded-HR) sensitivity proxy is emitted separately.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
from pathlib import Path

import numpy as np
import pandas as pd


ANALYSIS_ROOT = Path(__file__).resolve().parents[1]
MODEL_VALIDATION_ROOT = ANALYSIS_ROOT.parent
RAW = MODEL_VALIDATION_ROOT / "data" / "raw"
FROZEN_PROCESSED = MODEL_VALIDATION_ROOT / "data" / "processed"
OUT = ANALYSIS_ROOT / "data" / "processed"
TABLES = ANALYSIS_ROOT / "tables"
REFERENCES = ANALYSIS_ROOT / "references"
ICGC_UNFILTERED_SOURCE_VALUE = os.environ.get("HCC_OS_ICGC_FPKM")
ICGC_UNFILTERED_SOURCE = (
    Path(ICGC_UNFILTERED_SOURCE_VALUE).expanduser().resolve()
    if ICGC_UNFILTERED_SOURCE_VALUE
    else RAW / "Merge_RNAseq_FKPM_Symbol.txt"
)
B = 1000
SEED = 20260906
TIE_TOLERANCE = 1e-8


MODELS: dict[str, dict] = {
    "Liu2021_FAIS_signature1": {
        "display_name": "Liu et al. 2021 — FAIS signature 1",
        "short_citation": "Liu et al., Front Mol Biosci, 2021",
        "doi": "10.3389/fmolb.2021.766609",
        "url": "https://www.frontiersin.org/journals/molecular-biosciences/articles/10.3389/fmolb.2021.766609/full",
        "reported_formula": "−0.13 × FCN3 + 0.11 × SPP1 + 0.12 × IGFBP3 + 0.16 × MCM7",
        "coefficients": {"FCN3": -0.13, "SPP1": 0.11, "IGFBP3": 0.12, "MCM7": 0.16},
        "coefficient_source": "Results: FAIS signature-1 stepwise Cox formula.",
        "original_preprocessing": "Exact abundance unit and transform are not sufficiently reported.",
        "frozen_operationalization": "Reported coefficients applied to frozen log-scale expression; no refitting.",
        "formula_status": "complete_reported_formula",
        "comparison_tier": "primary",
    },
    "Tian2021_five_gene": {
        "display_name": "Tian et al. 2021 — five-gene signature",
        "short_citation": "Tian et al., Front Med, 2021",
        "doi": "10.3389/fmed.2021.681388",
        "url": "https://www.frontiersin.org/journals/medicine/articles/10.3389/fmed.2021.681388/full",
        "reported_formula": "Multivariable Cox risk value; numerical coefficients not reported.",
        "coefficients": {},
        "reported_genes": ["CDC20", "TOP2A", "RRM2", "UBE2C", "AOX1"],
        "coefficient_source": "Five genes are reported, but the fitted multivariable Cox coefficients are absent.",
        "original_preprocessing": "Not evaluable because the numerical score is unavailable.",
        "frozen_operationalization": "Not scored; coefficients were neither inferred nor refitted.",
        "formula_status": "not_reproducible_coefficients_not_reported",
        "comparison_tier": "not_scored",
    },
    "FuSong2021_pyroptosis_3gene": {
        "display_name": "Fu & Song 2021 — pyroptosis three-gene signature",
        "short_citation": "Fu & Song, Front Cell Dev Biol, 2021",
        "doi": "10.3389/fcell.2021.748039",
        "url": "https://www.frontiersin.org/journals/cell-and-developmental-biology/articles/10.3389/fcell.2021.748039/full",
        "reported_formula": "0.0182 × GSDME + 0.0005 × GPX4 + 0.0188 × SCAF11",
        "coefficients": {"GSDME": 0.0182, "GPX4": 0.0005, "SCAF11": 0.0188},
        "coefficient_source": "Results: explicit three-gene risk-score formula.",
        "original_preprocessing": "The article specifies R scale() for expression, but not archived means/SDs.",
        "frozen_operationalization": "Within-cohort gene-wise z score with sample SD followed by published coefficients; no refitting.",
        "formula_status": "complete_reported_formula",
        "comparison_tier": "primary",
        "score_transform": "within_cohort_zscore",
    },
    "Lin2021_inflammatory_8gene": {
        "display_name": "Lin et al. 2021 — inflammatory-response eight-gene signature",
        "short_citation": "Lin et al., Front Oncol, 2021",
        "doi": "10.3389/fonc.2021.644416",
        "url": "https://www.frontiersin.org/journals/oncology/articles/10.3389/fonc.2021.644416/full",
        "reported_formula": "0.118 × SLC7A1 + 0.114 × RIPK2 + 0.113 × NOD2 + 0.022 × ADORA2B + 0.058 × MEP1A + 0.051 × ITGA5 + 0.016 × P2RX4 + 0.018 × SERPINE1",
        "coefficients": {"SLC7A1": 0.118, "RIPK2": 0.114, "NOD2": 0.113, "ADORA2B": 0.022, "MEP1A": 0.058, "ITGA5": 0.051, "P2RX4": 0.016, "SERPINE1": 0.018},
        "coefficient_source": "Results: explicit risk-score formula.",
        "original_preprocessing": "A normalized expression matrix is cited, but exact unit and transform are incomplete.",
        "frozen_operationalization": "Reported linear predictor applied to frozen log-scale expression; no refitting. Its exp transform is rank-monotone and leaves C-index unchanged.",
        "formula_status": "complete_reported_formula",
        "comparison_tier": "primary",
    },
    "Wang2022_cuproptosis_5gene": {
        "display_name": "Wang et al. 2022 — cuproptosis five-gene signature",
        "short_citation": "Wang et al., Front Mol Biosci, 2022",
        "doi": "10.3389/fmolb.2022.1001788",
        "url": "https://www.frontiersin.org/journals/molecular-biosciences/articles/10.3389/fmolb.2022.1001788/full",
        "reported_formula": "Cu-PS = Σ(coef_i × mRNA_i); numerical coefficients not reported.",
        "coefficients": {},
        "reported_genes": ["C7", "MAGEA6", "HK2", "CYP26B1", "EPO"],
        "reported_hazard_ratios": {"C7": 0.917, "MAGEA6": 1.098, "HK2": 1.147, "CYP26B1": 1.145, "EPO": 1.077},
        "coefficient_source": "Figure 3F reports rounded multivariable HRs, not fitted score coefficients.",
        "original_preprocessing": "TPM conversion plus cross-cohort SVA are described, but parameters and transformed values are not published.",
        "frozen_operationalization": "Not exact-scored. Sensitivity only: ln(rounded Figure 3F HR) on frozen log-scale expression.",
        "formula_status": "not_reproducible_coefficients_not_reported",
        "comparison_tier": "rounded_HR_proxy_sensitivity",
    },
    "HongCai2022_OSRG_8gene": {
        "display_name": "Hong & Cai 2022 — oxidative-stress eight-gene signature",
        "short_citation": "Hong & Cai, Disease Markers, 2022",
        "doi": "10.1155/2022/6201987",
        "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC9484914/",
        "reported_formula": "0.069 × G6PD + 0.177 × MT3 + 0.206 × CBX2 + 0.063 × CDKN2B + 0.078 × CCNA2 + 0.164 × MAPT + 0.248 × EZH2 + 0.213 × SLC7A11",
        "coefficients": {"G6PD": 0.069, "MT3": 0.177, "CBX2": 0.206, "CDKN2B": 0.063, "CCNA2": 0.078, "MAPT": 0.164, "EZH2": 0.248, "SLC7A11": 0.213},
        "coefficient_source": "Results 3.1: explicit risk-score formula.",
        "original_preprocessing": "Exact abundance unit, log transform and standardization details are not sufficiently reported.",
        "frozen_operationalization": "Reported coefficients applied to frozen log-scale expression; no refitting.",
        "formula_status": "complete_reported_formula",
        "comparison_tier": "primary",
    },
    "Ma2024_OS_ER_stress_5gene": {
        "display_name": "Ma et al. 2024 — oxidative/ER-stress five-gene signature",
        "short_citation": "Ma et al., Comb Chem High Throughput Screen, 2024",
        "doi": "10.2174/0113862073257308231026073951",
        "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC11497145/",
        "reported_formula": "−0.2310 × IL18RAP + 0.2303 × ECT2 − 0.1228 × PPARGC1A + 0.1809 × STC2 + 0.0567 × NQO1",
        "coefficients": {"IL18RAP": -0.2310, "ECT2": 0.2303, "PPARGC1A": -0.1228, "STC2": 0.1809, "NQO1": 0.0567},
        "coefficient_source": "Explicit five-gene formula; the equation says SCT2 but the selected gene/text identify STC2.",
        "original_preprocessing": "Exact abundance unit, transform and standardization parameters are incompletely reported.",
        "frozen_operationalization": "Reported coefficients applied to frozen log-scale expression; STC2 resolves the apparent SCT2 typo; no refitting.",
        "formula_status": "complete_reported_formula",
        "comparison_tier": "primary",
    },
}

# The raw ICGC source uses historical symbol DFNA5 for the current GSDME gene.
ICGC_SOURCE_GENE = {"GSDME": "DFNA5"}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def model_genes(model: dict) -> list[str]:
    return list(model.get("coefficients", {})) or list(model.get("reported_genes", []))


def load_frozen_cohort(cohort: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    risk = pd.read_csv(RAW / f"{cohort}_frozen_risks.csv")
    expected = {"ID", "OS.time", "OS", "RS"}
    if set(risk.columns) != expected:
        raise ValueError(f"Unexpected columns in {cohort}: {list(risk.columns)}")
    risk = risk.rename(columns={"OS.time": "time_days", "OS": "event", "RS": "OSARS_frozen"}).set_index("ID")
    risk.index = risk.index.astype(str)
    if risk.index.has_duplicates:
        raise ValueError(f"Duplicate frozen IDs in {cohort}")
    filename = "TCGA_full_expression_logscale.csv.gz" if cohort == "TCGA-LIHC" else "ICGC_full_expression_logscale.csv.gz"
    expression = pd.read_csv(RAW / filename, index_col=0).T
    expression.index = expression.index.astype(str)
    if expression.index.has_duplicates or expression.columns.has_duplicates:
        raise ValueError(f"Duplicate expression IDs/genes in {filename}")
    missing = risk.index.difference(expression.index)
    if len(missing):
        raise ValueError(f"{cohort} expression misses frozen risk IDs: {missing[:5].tolist()}")
    return risk, expression.loc[risk.index]


def parse_fpkm(value: str) -> float:
    if value is None or value.strip() in {"", "NA", "NaN", "nan"}:
        return np.nan
    answer = float(value)
    if answer < 0:
        raise ValueError(f"Negative FPKM: {answer}")
    return answer


def recover_icgc_signature_expression(ids: pd.Index, genes: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    if not ICGC_UNFILTERED_SOURCE.exists():
        raise FileNotFoundError(f"Unfiltered ICGC source absent: {ICGC_UNFILTERED_SOURCE}")
    source_to_canonical: dict[str, list[str]] = {}
    for gene in genes:
        source_to_canonical.setdefault(ICGC_SOURCE_GENE.get(gene, gene), []).append(gene)
    source_values: dict[str, np.ndarray] = {}
    with ICGC_UNFILTERED_SOURCE.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        lookup = {identifier: position for position, identifier in enumerate(header)}
        absent = [identifier for identifier in ids if identifier not in lookup]
        if absent:
            raise ValueError(f"Unfiltered ICGC source misses IDs: {absent[:5]}")
        positions = [lookup[identifier] for identifier in ids]
        for row in reader:
            if not row:
                continue
            source_gene = row[0].strip()
            if source_gene not in source_to_canonical:
                continue
            if source_gene in source_values:
                raise ValueError(f"Duplicate source gene: {source_gene}")
            raw = [parse_fpkm(row[position]) if position < len(row) else np.nan for position in positions]
            source_values[source_gene] = np.log2(np.asarray(raw, dtype=float) + 1.0)
    missing_genes = sorted(set(source_to_canonical).difference(source_values))
    if missing_genes:
        raise ValueError(f"Missing required ICGC source rows: {missing_genes}")
    recovered, mapping = {}, []
    for source_gene, canonical_genes in source_to_canonical.items():
        for gene in canonical_genes:
            recovered[gene] = source_values[source_gene]
            mapping.append({
                "cohort": "ICGC-LIRI",
                "canonical_gene": gene,
                "source_gene": source_gene,
                "mapping_type": "official_historical_alias" if source_gene != gene else "exact_symbol",
                "raw_transform": "log2(FPKM + 1)",
            })
    return pd.DataFrame(recovered, index=ids), pd.DataFrame(mapping)


def cindex(time: np.ndarray, event: np.ndarray, score: np.ndarray) -> float:
    """Harrell C-index for a score increasing with risk."""
    time, event, score = np.asarray(time, float), np.asarray(event, bool), np.asarray(score, float)
    if not (np.isfinite(time).all() and np.isfinite(score).all()):
        raise ValueError("Non-finite time or score")
    comparable = event[:, None] & ((time[:, None] < time[None, :]) | ((time[:, None] == time[None, :]) & ~event[None, :]))
    i, j = np.where(comparable)
    if not len(i):
        return np.nan
    delta = score[i] - score[j]
    agreement = (delta > TIE_TOLERANCE).astype(float) + 0.5 * (np.abs(delta) <= TIE_TOLERANCE)
    return float(agreement.mean())


def bootstrap_cindices(survival: pd.DataFrame, columns: list[str], cohort: str, analysis_set: str) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Fixed-score paired sample-row percentile bootstrap."""
    time = survival["time_days"].to_numpy(float)
    event = survival["event"].to_numpy(bool)
    matrix = survival.loc[:, columns].to_numpy(float)
    if not np.isfinite(matrix).all():
        raise ValueError(f"Non-finite score in {cohort}/{analysis_set}")
    comparable = event[:, None] & ((time[:, None] < time[None, :]) | ((time[:, None] == time[None, :]) & ~event[None, :]))
    i, j = np.where(comparable)
    if not len(i):
        raise ValueError(f"No comparable pairs in {cohort}/{analysis_set}")
    delta = matrix[i] - matrix[j]
    agreement = (delta > TIE_TOLERANCE).astype(float) + 0.5 * (np.abs(delta) <= TIE_TOLERANCE)
    rng = np.random.default_rng(SEED)
    draws = np.full((B, len(columns)), np.nan)
    for replicate in range(B):
        sampled = rng.integers(len(survival), size=len(survival))
        count = np.bincount(sampled, minlength=len(survival))
        weight = count[i] * count[j]
        if weight.sum():
            draws[replicate] = np.sum(weight[:, None] * agreement, axis=0) / weight.sum()
    point = np.asarray([cindex(time, event, matrix[:, k]) for k in range(matrix.shape[1])])
    metric_rows, draw_rows = [], []
    for k, predictor in enumerate(columns):
        finite = draws[:, k][np.isfinite(draws[:, k])]
        low, high = np.quantile(finite, [0.025, 0.975])
        metric_rows.append({
            "analysis_set": analysis_set, "cohort": cohort, "predictor": predictor,
            "cindex": point[k], "ci_lower": low, "ci_upper": high, "n": len(survival),
            "events": int(event.sum()), "comparable_pairs": int(len(i)),
            "bootstrap_B": B, "bootstrap_valid": int(len(finite)),
            "ci_method": "paired sample-row percentile bootstrap; fixed published scores",
            "tie_rule": "higher risk score concordant; ties=0.5; tolerance=1e-8",
        })
        draw_rows.extend({
            "analysis_set": analysis_set, "cohort": cohort, "replicate": replicate + 1,
            "predictor": predictor, "cindex": draws[replicate, k],
        } for replicate in range(B))
    ref = columns.index("OSARS_frozen")
    difference_rows = []
    for k, predictor in enumerate(columns):
        if predictor == "OSARS_frozen":
            continue
        difference = draws[:, ref] - draws[:, k]
        finite = difference[np.isfinite(difference)]
        difference_rows.append({
            "analysis_set": analysis_set, "cohort": cohort, "reference": "OSARS_frozen",
            "comparator": predictor, "difference_OSARS_minus_comparator": point[ref] - point[k],
            "ci_lower": np.quantile(finite, 0.025), "ci_upper": np.quantile(finite, 0.975),
            "n": len(survival), "events": int(event.sum()), "bootstrap_B": B,
            "bootstrap_valid": int(len(finite)),
            "interpretation": "Descriptive paired contrast; development/validation cohort overlap is not removed.",
        })
    return pd.DataFrame(metric_rows), pd.DataFrame(draw_rows), pd.DataFrame(difference_rows)


def score_complete_formula(model: dict, expression: pd.DataFrame) -> pd.Series:
    coefficients = pd.Series(model["coefficients"], dtype=float)
    missing = coefficients.index.difference(expression.columns)
    if len(missing):
        raise KeyError(f"Missing formula genes: {missing.tolist()}")
    values = expression.loc[:, coefficients.index]
    if model.get("score_transform") == "within_cohort_zscore":
        values = (values - values.mean(axis=0)) / values.std(axis=0, ddof=1)
    return values @ coefficients


def score_wang_proxy(model: dict, expression: pd.DataFrame) -> pd.Series:
    coefficients = pd.Series({gene: math.log(hr) for gene, hr in model["reported_hazard_ratios"].items()}, dtype=float)
    return expression.loc[:, coefficients.index] @ coefficients


def formula_audit() -> pd.DataFrame:
    rows = []
    for signature_id, model in MODELS.items():
        rows.append({
            "signature_id": signature_id, "display_name": model["display_name"],
            "citation": model["short_citation"], "doi": model["doi"], "article_url": model["url"],
            "reported_genes": ";".join(model_genes(model)), "reported_formula": model["reported_formula"],
            "coefficient_source": model["coefficient_source"], "original_preprocessing": model["original_preprocessing"],
            "frozen_operationalization": model["frozen_operationalization"],
            "formula_status": model["formula_status"], "comparison_tier": model["comparison_tier"],
        })
    return pd.DataFrame(rows)


def coefficient_audit() -> pd.DataFrame:
    rows = []
    for signature_id, model in MODELS.items():
        coefficients = model.get("coefficients", {})
        for gene in model_genes(model):
            row = {
                "signature_id": signature_id, "gene": gene,
                "published_coefficient": coefficients.get(gene, np.nan),
                "coefficient_type": "reported_cox_coefficient" if gene in coefficients else "not_reported",
                "formula_status": model["formula_status"],
            }
            if signature_id == "Wang2022_cuproptosis_5gene":
                hr = model["reported_hazard_ratios"][gene]
                row.update({
                    "reported_multivariable_HR_figure3F": hr,
                    "sensitivity_proxy_coefficient_ln_rounded_HR": math.log(hr),
                    "coefficient_type": "not_reported; proxy uses ln(rounded HR)",
                })
            rows.append(row)
    return pd.DataFrame(rows)


def gene_availability(cohort: str, expression: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for signature_id, model in MODELS.items():
        for gene in model_genes(model):
            source_gene = ICGC_SOURCE_GENE.get(gene, gene) if cohort == "ICGC-LIRI" else gene
            value = expression[gene] if gene in expression.columns else pd.Series(np.nan, index=expression.index)
            rows.append({
                "cohort": cohort, "signature_id": signature_id, "gene": gene, "source_gene": source_gene,
                "source_mapping": "official_historical_alias" if source_gene != gene else "exact_symbol",
                "gene_row_present": gene in expression.columns, "n": len(expression),
                "available_n": int(value.notna().sum()), "missing_n": int(value.isna().sum()),
                "complete_for_all_patients": bool(value.notna().all()),
            })
    return pd.DataFrame(rows)


def write_summary(metrics: pd.DataFrame) -> None:
    primary = metrics[metrics.analysis_set == "primary_complete_formula"]
    proxy = metrics[metrics.analysis_set == "wang_rounded_hr_proxy_sensitivity"]
    lines = [
        "# Published-signature head-to-head C-index", "",
        "Harrell C-index for scores increasing with risk. Intervals are 1,000 paired sample-row percentile bootstrap replicates.", "",
        "## Primary comparison: numerical published formulas only", "",
        "| Cohort | n (events) | Predictor | C-index (95% CI) |",
        "|---|---:|---|---:|",
    ]
    for _, row in primary.iterrows():
        lines.append(f"| {row.cohort} | {int(row.n)} ({int(row.events)}) | {row.predictor} | {row.cindex:.3f} ({row.ci_lower:.3f}–{row.ci_upper:.3f}) |")
    lines += [
        "", "Tian et al. 2021 is unscored because the five multivariable Cox coefficients are not reported. Wang et al. 2022 is excluded from the primary table because numerical coefficients and fitted preprocessing parameters are not reported.", "",
        "## Sensitivity only: Wang log(rounded-HR) proxy", "",
        "This is not an exact rerun of Wang et al. 2022. It applies ln of rounded Figure 3F multivariable HRs to frozen log-scale expression.", "",
        "| Cohort | n (events) | Predictor | C-index (95% CI) |",
        "|---|---:|---|---:|",
    ]
    for _, row in proxy.iterrows():
        lines.append(f"| {row.cohort} | {int(row.n)} ({int(row.events)}) | {row.predictor} | {row.cindex:.3f} ({row.ci_lower:.3f}–{row.ci_upper:.3f}) |")
    lines += [
        "", "## Interpretation guardrails", "",
        "- No published signature was refitted, recalibrated, or optimized.",
        "- Primary rows use identical complete cases within each cohort.",
        "- TCGA-LIHC is an apparent comparison for signatures developed there; ICGC-LIRI overlaps historical validation cohorts for some signatures. Results are descriptive, not a new independent validation.",
        "- Paired same-patient contrasts are in paired_cindex_difference_vs_OSARS.csv.",
    ]
    (TABLES / "head_to_head_cindex_summary.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    for directory in (OUT, TABLES, REFERENCES):
        directory.mkdir(parents=True, exist_ok=True)
    formula_audit().to_csv(REFERENCES / "formula_audit.csv", index=False)
    coefficient_audit().to_csv(REFERENCES / "published_signature_coefficients.csv", index=False)

    tcga_survival, tcga_expression = load_frozen_cohort("TCGA-LIHC")
    icgc_survival, icgc_frozen_expression = load_frozen_cohort("ICGC-LIRI")
    all_genes = sorted({gene for model in MODELS.values() for gene in model_genes(model)})
    icgc_recovered, icgc_mapping = recover_icgc_signature_expression(icgc_survival.index, all_genes)
    icgc_expression = icgc_frozen_expression.copy()
    for gene in icgc_recovered:
        icgc_expression[gene] = icgc_recovered[gene]
    icgc_mapping.to_csv(OUT / "ICGC_signature_gene_source_mapping.csv", index=False)

    consistency_rows = []
    for gene in all_genes:
        if gene not in icgc_frozen_expression:
            continue
        original = icgc_frozen_expression[gene].to_numpy(float)
        recovered = icgc_recovered[gene].to_numpy(float)
        keep = np.isfinite(original) & np.isfinite(recovered)
        maximum = float(np.max(np.abs(original[keep] - recovered[keep]))) if keep.any() else np.nan
        consistency_rows.append({
            "gene": gene, "source_gene": ICGC_SOURCE_GENE.get(gene, gene),
            "n_compared": int(keep.sum()), "max_abs_difference": maximum,
            "matching_within_1e-10": bool(keep.any() and maximum < 1e-10),
        })
    consistency = pd.DataFrame(consistency_rows)
    consistency.to_csv(OUT / "ICGC_source_recovery_consistency.csv", index=False)
    if len(consistency) and not consistency["matching_within_1e-10"].all():
        raise AssertionError("Recovered ICGC expression does not match retained frozen values")

    pd.concat([
        gene_availability("TCGA-LIHC", tcga_expression),
        gene_availability("ICGC-LIRI", icgc_expression),
    ], ignore_index=True).to_csv(OUT / "input_gene_availability.csv", index=False)

    frames = {}
    for cohort, survival, expression in [
        ("TCGA-LIHC", tcga_survival, tcga_expression),
        ("ICGC-LIRI", icgc_survival, icgc_expression),
    ]:
        frame = survival.copy()
        for signature_id, model in MODELS.items():
            if model["formula_status"] == "complete_reported_formula":
                frame[signature_id] = score_complete_formula(model, expression)
            elif signature_id == "Wang2022_cuproptosis_5gene":
                frame["Wang2022_cuproptosis_5gene_ln_roundedHR_proxy"] = score_wang_proxy(model, expression)
        frame["analysis_role"] = "development-overlap apparent comparison" if cohort == "TCGA-LIHC" else "historical-validation/setting-overlap descriptive comparison"
        frame.to_csv(OUT / f"{cohort}_all_scores.csv", index_label="ID")
        frames[cohort] = frame

    primary = ["OSARS_frozen", "Liu2021_FAIS_signature1", "FuSong2021_pyroptosis_3gene", "Lin2021_inflammatory_8gene", "HongCai2022_OSRG_8gene", "Ma2024_OS_ER_stress_5gene"]
    proxy = "Wang2022_cuproptosis_5gene_ln_roundedHR_proxy"
    sensitivity = primary + [proxy]
    memberships, metric_list, draw_list, difference_list = [], [], [], []
    for cohort, frame in frames.items():
        primary_keep = frame[primary].notna().all(axis=1)
        sensitivity_keep = frame[sensitivity].notna().all(axis=1)
        memberships.extend({
            "cohort": cohort, "ID": identifier,
            "primary_complete_formula": bool(primary_keep.loc[identifier]),
            "wang_rounded_hr_proxy_sensitivity": bool(sensitivity_keep.loc[identifier]),
        } for identifier in frame.index)
        for label, columns, keep in [
            ("primary_complete_formula", primary, primary_keep),
            ("wang_rounded_hr_proxy_sensitivity", sensitivity, sensitivity_keep),
        ]:
            selected = frame.loc[keep, ["time_days", "event"] + columns]
            metrics, draws, differences = bootstrap_cindices(selected, columns, cohort, label)
            metric_list.append(metrics)
            draw_list.append(draws)
            difference_list.append(differences)
    pd.DataFrame(memberships).to_csv(OUT / "analysis_set_membership.csv", index=False)
    metrics = pd.concat(metric_list, ignore_index=True)
    draws = pd.concat(draw_list, ignore_index=True)
    differences = pd.concat(difference_list, ignore_index=True)
    metrics.to_csv(TABLES / "head_to_head_cindex.csv", index=False)
    draws.to_csv(OUT / "bootstrap_cindex_draws.csv.gz", index=False, compression="gzip")
    differences.to_csv(TABLES / "paired_cindex_difference_vs_OSARS.csv", index=False)

    prior_path = FROZEN_PROCESSED / "discrimination_metrics_with_CI.csv"
    prior = pd.read_csv(prior_path) if prior_path.exists() else pd.DataFrame()
    reproduction_rows = []
    for cohort, frame in frames.items():
        observed = cindex(frame.time_days, frame.event, frame.OSARS_frozen)
        old = prior[(prior.cohort == cohort) & (prior.predictor == "OSARS_frozen") & (prior.metric == "C-index")]
        expected = float(old.estimate.iloc[0]) if len(old) else np.nan
        reproduction_rows.append({
            "cohort": cohort, "full_frozen_n": len(frame), "full_frozen_events": int(frame.event.sum()),
            "recalculated_OSARS_cindex": observed, "existing_frozen_pipeline_cindex": expected,
            "absolute_difference": abs(observed - expected) if np.isfinite(expected) else np.nan,
        })
    reproduction = pd.DataFrame(reproduction_rows)
    reproduction.to_csv(OUT / "OSARS_cindex_reproduction_check.csv", index=False)
    delta = reproduction.absolute_difference.dropna()
    if len(delta) and delta.max() > 1e-12:
        raise AssertionError("OSARS C-index fails frozen-pipeline reproduction check")

    manifest = {
        "analysis": "Published HCC prognostic signature head-to-head C-index",
        "created": "2026-09-06", "seed": SEED, "bootstrap_B": B,
        "primary_analysis_set": "OSARS plus five signatures with complete reported formulas",
        "sensitivity_analysis_set": "Primary predictors plus Wang ln(rounded HR) proxy; not exact rerun",
        "unscored_article": "Tian2021_five_gene: numerical multivariable coefficients not reported",
        "raw_inputs": {path.name: sha256(path) for path in [
            RAW / "TCGA-LIHC_frozen_risks.csv", RAW / "ICGC-LIRI_frozen_risks.csv",
            RAW / "TCGA_full_expression_logscale.csv.gz", RAW / "ICGC_full_expression_logscale.csv.gz",
            ICGC_UNFILTERED_SOURCE,
        ]},
    }
    (OUT / "source_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    write_summary(metrics)
    print("Completed published-signature C-index analysis")
    print(metrics.to_string(index=False))


if __name__ == "__main__":
    main()
