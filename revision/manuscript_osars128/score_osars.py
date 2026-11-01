"""Score already-preprocessed, finite OSARS expression profiles using a frozen tree export."""
import argparse
import csv
import hashlib
import io
import math
import os
import tempfile
from pathlib import Path


MODEL_SHA256 = {
    "model.tsv": "6e4d091e27421729b5cc3a3427ea16dfcdc7162c03325233b5d7a7b75bac8bdc",
    "genes.txt": "e82c4193cebdd9a53ad373943861174befce2fa172517ce5047743b7bf324759",
    "trees.tsv": "e96002949afb4e197c7e783cb8933b170bf573845c016f3a37863fff72461df8",
}


def verified_text(model_dir, name):
    data = (model_dir / name).read_bytes().replace(b"\r\n", b"\n")
    if hashlib.sha256(data).hexdigest() != MODEL_SHA256[name]:
        raise ValueError(f"{name} differs from the manuscript model export.")
    return data.decode("ascii")


def load_model(model_dir):
    with io.StringIO(verified_text(model_dir, "model.tsv")) as stream:
        meta = dict(csv.reader(stream, delimiter="\t"))
    if meta != {"format": "osars-numeric-gbm-v1", "n_features": "128", "n_trees": "1696",
                "initF": "0", "prediction_scale": "link"}:
        raise ValueError("Unexpected model metadata.")
    genes = verified_text(model_dir, "genes.txt").splitlines()
    if len(genes) != 128 or len(set(genes)) != 128 or any(not g for g in genes):
        raise ValueError("Expected 128 unique, ordered gene names.")

    trees = [[] for _ in range(1696)]
    with io.StringIO(verified_text(model_dir, "trees.tsv")) as stream:
        rows = csv.DictReader(stream, delimiter="\t")
        if rows.fieldnames != ["tree", "node", "feature", "split_or_leaf", "left", "right", "missing"]:
            raise ValueError("Unexpected tree columns.")
        for row in rows:
            tree_id = int(row["tree"])
            if not 1 <= tree_id <= len(trees):
                raise ValueError(f"Invalid tree number: {tree_id}")
            tree = trees[tree_id - 1]
            if int(row["node"]) != len(tree):
                raise ValueError(f"Unexpected node order in tree {tree_id}.")
            tree.append((int(row["feature"]), float(row["split_or_leaf"]),
                         int(row["left"]), int(row["right"]), int(row["missing"])))
    for tree in trees:
        if not tree:
            raise ValueError("Missing tree.")
        for node, (feature, value, left, right, missing) in enumerate(tree):
            if not math.isfinite(value):
                raise ValueError("Nonfinite tree value.")
            if feature == -1 and (left, right, missing) == (-1, -1, -1):
                continue
            if not (0 <= feature < len(genes) and
                    all(node < child < len(tree) for child in (left, right, missing))):
                raise ValueError("Invalid split or child node.")
    return genes, trees


def predict(profile, trees):
    score = 0.0  # initF for this Cox GBM
    for tree in trees:
        node = 0
        while True:
            feature, value, left, right, _ = tree[node]
            if feature == -1:
                score += value  # shrinkage is already included in the leaf
                break
            node = left if profile[feature] < value else right
    return score


def score_csv(model_dir, input_csv, output_csv):
    genes, trees = load_model(model_dir)
    scores = {}
    with input_csv.open(newline="", encoding="utf-8-sig") as stream:
        rows = csv.DictReader(stream)
        header = rows.fieldnames or []
        if len(header) != len(genes) + 1 or set(header) != set(genes) | {"sample_id"}:
            raise ValueError("CSV must have sample_id and each of the 128 gene columns exactly once.")
        for line, row in enumerate(rows, start=2):
            if None in row or any(v is None for v in row.values()):
                raise ValueError(f"Wrong number of values at row {line}.")
            sample_id = row["sample_id"]
            if not sample_id.strip() or sample_id in scores:
                raise ValueError(f"Empty or duplicate sample_id at row {line}.")
            profile = [float(row[g]) for g in genes]
            if not all(math.isfinite(x) for x in profile):
                raise ValueError(f"Preprocess or impute nonfinite expression values at row {line}.")
            scores[sample_id] = predict(profile, trees)
    if not scores:
        raise ValueError("The input has no samples.")
    with tempfile.TemporaryDirectory(dir=output_csv.parent) as stage:
        staged = Path(stage) / "scores.csv"
        with staged.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream, lineterminator="\n")
            writer.writerow(["sample_id", "OSARS"])
            writer.writerows(scores.items())
        os.link(staged, output_csv)  # atomically publish without overwriting an existing file


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("model_dir", type=Path, help="Directory from export_frozen_model.R")
    cli.add_argument("input_csv", type=Path, help="CSV with sample_id and 128 named expression columns")
    cli.add_argument("output_csv", type=Path, help="New CSV for link-scale scores")
    args = cli.parse_args()
    try:
        score_csv(args.model_dir, args.input_csv, args.output_csv)
    except (OSError, ValueError) as exc:
        cli.exit(1, f"Cannot score profiles: {exc}\n")
