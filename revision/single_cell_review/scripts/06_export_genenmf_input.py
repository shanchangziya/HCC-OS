#!/usr/bin/env python3
"""Export the frozen GSE149614 malignant-cell GeneNMF input.

Only the 13,691 author-annotated primary-tumour malignant cells and the
2,000 genes that were present in every saved historical GeneNMF W matrix are
exported.  The matrix contains raw integer counts; normalization remains an
explicit step in the R analysis.
"""

from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmwrite


SCRIPT_DIR = Path(__file__).resolve().parent
OUT_ROOT = SCRIPT_DIR.parent
REVISION_ROOT = OUT_ROOT.parent
SOURCE_DIR = REVISION_ROOT / "single_cell" / "data" / "processed"
INPUT_DIR = OUT_ROOT / "data" / "genenmf_input"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    INPUT_DIR.mkdir(parents=True, exist_ok=True)

    h5ad_path = SOURCE_DIR / "GSE149614_counts.h5ad"
    universe_path = SOURCE_DIR / "NMF_eligible_gene_universe.csv"
    matrix_path = INPUT_DIR / "GSE149614_primary_malignant_2000genes_counts.mtx.gz"
    genes_path = INPUT_DIR / "GSE149614_primary_malignant_2000genes_genes.csv"
    cells_path = INPUT_DIR / "GSE149614_primary_malignant_2000genes_cells.csv"
    manifest_path = INPUT_DIR / "GSE149614_primary_malignant_2000genes_manifest.json"

    universe = pd.read_csv(universe_path)["gene"].astype(str).tolist()
    if len(universe) != 2000 or len(set(universe)) != 2000:
        raise ValueError("The frozen GeneNMF universe must contain 2,000 unique genes.")

    adata = ad.read_h5ad(h5ad_path, backed="r")
    if "primary_analysis" not in adata.obs or "patient" not in adata.obs:
        raise KeyError("Required frozen annotations are absent from the H5AD object.")

    selected = adata.obs["primary_analysis"].astype(bool).to_numpy()
    cell_names = adata.obs_names[selected].astype(str)
    if selected.sum() != 13691:
        raise ValueError(f"Expected 13,691 malignant cells, found {selected.sum():,}.")

    var_lookup = pd.Series(np.arange(adata.n_vars), index=adata.var_names.astype(str))
    missing = [gene for gene in universe if gene not in var_lookup.index]
    if missing:
        raise ValueError(f"Frozen GeneNMF genes missing from H5AD: {missing[:10]}")
    gene_indices = var_lookup.loc[universe].to_numpy(dtype=int)

    # AnnData is cells x genes; GeneNMF/Seurat expects genes x cells.
    counts = adata[selected, gene_indices].to_memory().X
    if not sparse.issparse(counts):
        counts = sparse.csr_matrix(counts)
    counts = counts.T.tocsc()
    if counts.shape != (2000, 13691):
        raise ValueError(f"Unexpected matrix shape: {counts.shape}")
    if counts.data.size and (counts.data.min() < 0 or not np.all(counts.data == np.floor(counts.data))):
        raise ValueError("The GeneNMF input is not a non-negative integer count matrix.")

    with gzip.open(matrix_path, "wb", compresslevel=6) as handle:
        mmwrite(handle, counts, field="integer", symmetry="general")

    pd.DataFrame({"gene": universe}).to_csv(genes_path, index=False)
    cell_meta = adata.obs.loc[cell_names, ["patient", "sample", "site"]].copy()
    cell_meta.insert(0, "cell", cell_names)
    cell_meta.to_csv(cells_path, index=False)

    patient_counts = cell_meta["patient"].value_counts().sort_index().astype(int).to_dict()
    manifest = {
        "source_h5ad": str(h5ad_path),
        "frozen_gene_universe": str(universe_path),
        "selection": "primary_analysis == TRUE (author-annotated primary-tumour malignant hepatocytes)",
        "matrix_orientation": "genes_x_cells",
        "matrix_values": "raw_integer_counts",
        "n_genes": int(counts.shape[0]),
        "n_cells": int(counts.shape[1]),
        "n_patients": int(cell_meta["patient"].nunique()),
        "n_samples": int(cell_meta["sample"].nunique()),
        "nnz": int(counts.nnz),
        "count_sum": int(counts.sum()),
        "cells_per_patient": patient_counts,
        "matrix_sha256": sha256(matrix_path),
        "genes_sha256": sha256(genes_path),
        "cells_sha256": sha256(cells_path),
    }
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
