#!/usr/bin/env python3
"""Validate staged spatial inputs and record repository-relative SHA-256 values.

Restricted source files are intentionally not downloaded by this repository.
Place them under ``data/raw`` with the layout below, then run this script.
"""

from pathlib import Path
import hashlib
import json


BASE = Path(__file__).resolve().parents[1]
RAW = BASE / "data" / "raw"
FILES = [
    "results/colocalization/cell2location_spot_correlations.csv",
    "results/colocalization/cell2location_tumor_high_quartile_enrichment.csv",
    "results/cell2location_spatial/spatial_run_summary.csv",
    "inputs/spatial_manifest.csv",
]
for section in ["HCC-1T", "HCC-2T", "HCC-3T", "HCC-4T"]:
    FILES.extend(
        [
            f"results/cell2location_spatial/{section}/cell_abundance_mean.csv",
            f"inputs/spatial/{section}/spatial_coords.tsv.gz",
        ]
    )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


missing = [relative for relative in FILES if not (RAW / relative).is_file()]
if missing:
    formatted = "\n  - ".join(missing)
    raise SystemExit(f"Missing required spatial inputs under {RAW}:\n  - {formatted}")

manifest = [
    {
        "relative_path": relative,
        "bytes": (RAW / relative).stat().st_size,
        "sha256": sha256(RAW / relative),
    }
    for relative in FILES
]
(BASE / "audit").mkdir(parents=True, exist_ok=True)
(BASE / "audit" / "source_manifest.json").write_text(
    json.dumps(manifest, indent=2), encoding="utf-8"
)
print(json.dumps({"validated_files": len(manifest)}, indent=2))
