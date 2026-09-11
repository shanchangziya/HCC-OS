#!/usr/bin/env python3
"""Read original pathology sources; export plotting-ready numerical data, no plots.

Python validates identifiers, source hashes, source joins, full feature recovery,
and the saved split. R is used only as a numerical backend for RData, glmnet,
survival estimators, and timeROC. All writes stay beneath --output-dir.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
from datetime import datetime, timezone
import numpy as np
import pandas as pd


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def ids(values):
    result = values.astype(str).str.extract(r"(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4})", expand=False)
    if result.isna().any():
        raise ValueError("Unparseable TCGA patient identifier")
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-dir", type=Path, default=Path(os.environ.get("HCC_OS_PATHOLOGY_SOURCE_DIR", "data/raw")))
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent.parent)
    parser.add_argument("--clinical-source", type=Path, default=Path(os.environ["HCC_OS_PATHOLOGY_CLINICAL"]) if os.environ.get("HCC_OS_PATHOLOGY_CLINICAL") else None)
    parser.add_argument("--rscript", default=os.environ.get("HCC_OS_RSCRIPT", "Rscript"))
    parser.add_argument("--r-library", type=Path)
    parser.add_argument("--bootstrap", type=int, default=500)
    parser.add_argument("--skip-tile-check", action="store_true",
                        help="Use complete saved patient feature matrix; tile audit remains unavailable")
    parser.add_argument("--prepare-only", action="store_true")
    args = parser.parse_args()
    out = args.output_dir.resolve()
    scripts = Path(__file__).resolve().parent
    raw, proc, logs = [out / p for p in ("data/raw", "data/processed", "logs")]
    for folder in (raw, proc, logs):
        folder.mkdir(parents=True, exist_ok=True)
    source = args.source_dir.resolve()
    source_specs = [
        ("resnet_object", source / "Resnet.Rdata"),
        ("pathology_split_object", source / "pathdat.Rdata"),
    ]
    if not args.skip_tile_check:
        source_specs.append(("tile_features", source / "resnet50_features.csv"))
    if args.clinical_source:
        source_specs.append(("clinical_covariates", args.clinical_source.resolve()))
    source_paths = [path for _, path in source_specs]
    for path in source_paths:
        if not path.is_file():
            raise FileNotFoundError(path)
    manifest = {
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "sources": [
            {"role": role, "source_id": path.name, "bytes": path.stat().st_size, "sha256": sha256(path)}
            for role, path in source_specs
        ],
        "bootstrap_replicates": args.bootstrap,
        "tile_audit_requested": not args.skip_tile_check,
        "scope": "exploratory corrected internal TCGA reassessment",
    }
    (proc / "source_manifest.json").write_text(json.dumps(manifest, indent=2))
    with (logs / "export.log").open("w") as log:
        subprocess.run([args.rscript, str(scripts / "export_original_objects.R"),
                        str(source), str(raw)], check=True, stdout=log, stderr=subprocess.STDOUT)
    features = pd.read_csv(raw / "original_patient_features.csv")
    assert len(features.columns) == 2049 and features.patient_id.is_unique
    feature_columns = [f"resnet{i}" for i in range(2048)]
    assert features.columns.tolist()[1:] == feature_columns
    assert np.isfinite(features[feature_columns].to_numpy()).all()
    qa = {"saved_full_feature_patients": len(features), "n_raw_features": 2048,
          "source_features_finite": True, "tile_check_performed": not args.skip_tile_check}
    if not args.skip_tile_check:
        tiles = pd.read_csv(source / "resnet50_features.csv", index_col=0)
        assert tiles.shape[1] == 2048
        assert np.isfinite(tiles.to_numpy()).all()
        tiles.columns = feature_columns
        tiles["patient_id"] = ids(pd.Series(tiles.index, index=tiles.index)).values
        counts = tiles.groupby("patient_id").size().rename("n_tiles")
        aggregate = tiles.groupby("patient_id")[feature_columns].mean()
        saved = features.set_index("patient_id")[feature_columns]
        assert set(aggregate.index) == set(saved.index)
        maxdiff = float(np.max(np.abs(aggregate.loc[saved.index].to_numpy() - saved.to_numpy())))
        assert maxdiff < 1e-6, f"Saved patient mean differs from tile aggregation: {maxdiff}"
        qa.update(n_tiles=len(tiles), tile_aggregate_max_abs_difference=maxdiff)
        features = features.merge(counts, on="patient_id", validate="one_to_one")
        counts.to_csv(proc / "patient_tile_counts.csv")
    else:
        features["n_tiles"] = np.nan
    split = pd.read_csv(raw / "original_pathology_split.csv")
    split["patient_id"] = ids(split.ID)
    assert split.patient_id.is_unique
    assert set(split.Set) == {"Train", "Test"}
    assert split.Set.value_counts().to_dict() == {"Train": 230, "Test": 100}
    assert set(split.patient_id).issubset(set(features.patient_id))
    rp = pd.read_csv(raw / "original_molecular_clinical.csv")
    rp["patient_id"] = ids(rp.ID)
    dup = rp[rp.duplicated("patient_id", keep=False)].copy()
    dup.to_csv(proc / "molecular_duplicate_ids.csv", index=False)
    for col in ("OS.time", "OS"):
        assert rp.groupby("patient_id")[col].nunique(dropna=False).max() == 1
    molecular = rp.groupby("patient_id").agg(
        molecular_OSARS=("RS", "mean"), time_days=("OS.time", "first"), event=("OS", "first")
    ).reset_index()
    cohort = split[["patient_id", "ID", "Set", "OS.time", "OS", "RS", "RiskScore"]].merge(
        molecular, on="patient_id", validate="one_to_one", how="left"
    )
    assert len(cohort) == 330 and cohort.time_days.notna().all()
    assert np.allclose(cohort["OS.time"], cohort.time_days)
    assert np.array_equal(cohort.OS, cohort.event)
    assert np.allclose(cohort.RS, cohort.molecular_OSARS)
    assert set(cohort.event.unique()).issubset({0, 1})
    assert (cohort.time_days > 0).all()
    cohort = cohort.rename(columns={"Set": "split", "RiskScore": "legacy_pathology_score"})
    cohort = cohort[["patient_id", "ID", "split", "time_days", "event", "molecular_OSARS",
                     "legacy_pathology_score"]].merge(features, on="patient_id", validate="one_to_one")
    excluded = features.loc[~features.patient_id.isin(cohort.patient_id), ["patient_id", "n_tiles"]].copy()
    excluded["reason"] = "absent from original matched pathology survival cohort"
    excluded.to_csv(proc / "excluded_patients.csv", index=False)
    # Secondary covariates are fixed before any new test performance is computed.
    if args.clinical_source:
        clinical = pd.read_csv(args.clinical_source)
        clinical["patient_id"] = ids(clinical.ID)
        assert clinical.patient_id.is_unique
        clinical = clinical[["patient_id", "time", "status", "OSARS", "Stage_Group", "Grade_Group", "Sex"]]
        cohort = cohort.merge(clinical, on="patient_id", how="left", validate="one_to_one")
        matched = cohort.time.notna()
        assert np.allclose(cohort.loc[matched, "time"], cohort.loc[matched, "time_days"])
        assert np.array_equal(cohort.loc[matched, "status"], cohort.loc[matched, "event"])
        assert np.allclose(cohort.loc[matched, "OSARS"], cohort.loc[matched, "molecular_OSARS"])
        for column, labels in {"Stage_Group": {"Stage I-II": 0, "Stage III-IV": 1},
                               "Grade_Group": {"Grade G1-G2": 0, "Grade G3-G4": 1},
                               "Sex": {"Female": 0, "Male": 1}}.items():
            target = {"Stage_Group": "stage_advanced", "Grade_Group": "grade_high", "Sex": "sex_male"}[column]
            assert set(cohort[column].dropna()).issubset(labels)
            cohort[target] = cohort[column].map(labels)
        cohort = cohort.drop(columns=["time", "status", "OSARS"])
        qa["clinical_matched"] = int(matched.sum())
        qa["clinical_complete"] = int(cohort[["stage_advanced", "grade_high", "sex_male"]].notna().all(axis=1).sum())
    qa.update(n_matched=len(cohort), n_excluded=len(excluded),
              n_train=int((cohort.split == "Train").sum()), n_test=int((cohort.split == "Test").sum()),
              train_test_patient_overlap=0, original_split_recovered_from="pathdat.Rdata dat_all$Set",
              outcome_unit="days, cross-checked against existing TCGA clinical table",
              legacy_seed=2026, cv_seed=20260905, cutoff_rule="High iff score > training median")
    cohort.to_csv(proc / "analysis_input.csv", index=False)
    cohort[["patient_id", "split", "event", "time_days", "n_tiles"]].to_csv(proc / "patient_manifest.csv", index=False)
    (proc / "input_qa.json").write_text(json.dumps(qa, indent=2))
    print(json.dumps(qa, indent=2), flush=True)
    if args.prepare_only:
        return
    command = [args.rscript, str(scripts / "fit_and_export.R"), str(proc), str(args.bootstrap)]
    if args.r_library:
        command.append(str(args.r_library.resolve()))
    with (logs / "analysis.log").open("w") as log:
        subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)
    for path, original in zip(source_paths, manifest["sources"]):
        assert sha256(path) == original["sha256"], "Original source changed during run"
    manifest["finished_utc"] = datetime.now(timezone.utc).isoformat()
    manifest["sources_unchanged_after_analysis"] = True
    manifest["outputs"] = [{"path": str(p.relative_to(out)), "sha256": sha256(p)}
                           for p in sorted(proc.glob("*")) if p.is_file() and p.name != "source_manifest.json"]
    (proc / "source_manifest.json").write_text(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
