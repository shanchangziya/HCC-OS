#!/usr/bin/env python3
"""Run cell2location deconvolution for each HCC Visium sample."""

from __future__ import annotations

import argparse
import gzip
import os
from pathlib import Path


def require_deps():
    missing = []
    for mod in ("anndata", "scanpy", "cell2location", "torch", "pandas"):
        try:
            __import__(mod)
        except Exception as exc:  # pragma: no cover
            missing.append(f"{mod}: {type(exc).__name__}: {exc}")
    if missing:
        raise SystemExit(
            "Missing Python dependencies for cell2location:\n"
            + "\n".join(missing)
            + "\nCreate/activate the environment from env/cell2location.yml first."
        )


def main() -> None:
    require_deps()

    import anndata as ad
    import numpy as np
    import pandas as pd
    import scanpy as sc
    import torch
    from scipy import io
    from cell2location.models import Cell2location

    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", default="inputs/spatial_manifest.csv")
    parser.add_argument("--signatures", default="results/cell2location_reference/cell_type_signatures.csv")
    parser.add_argument("--output-dir", default="results/cell2location_spatial")
    parser.add_argument("--max-epochs", type=int, default=None)
    parser.add_argument("--sample", default=None, help="Optional single sample to run")
    parser.add_argument("--N-cells-per-location", type=float, default=8.0)
    parser.add_argument("--detection-alpha", type=float, default=20.0)
    parser.add_argument("--seed", type=int, default=123)
    args = parser.parse_args()

    base = Path(__file__).resolve().parent.parent
    manifest_path = (base / args.manifest).resolve() if not os.path.isabs(args.manifest) else Path(args.manifest)
    signatures_path = (base / args.signatures).resolve() if not os.path.isabs(args.signatures) else Path(args.signatures)
    out_dir = (base / args.output_dir).resolve() if not os.path.isabs(args.output_dir) else Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    np.random.seed(args.seed)
    torch.manual_seed(args.seed)
    use_gpu = bool(torch.cuda.is_available())
    max_epochs = args.max_epochs if args.max_epochs is not None else (3000 if use_gpu else 800)

    signatures = pd.read_csv(signatures_path, index_col=0)
    manifest = pd.read_csv(manifest_path)
    if args.sample is not None:
        manifest = manifest.loc[manifest["sample"].astype(str) == args.sample].copy()
        if manifest.empty:
            raise SystemExit(f"Sample not found in manifest: {args.sample}")
    summary_rows = []

    for row in manifest.itertuples(index=False):
        sample = str(row.sample)
        sample_dir = out_dir / sample
        sample_dir.mkdir(parents=True, exist_ok=True)

        h5ad = getattr(row, "h5ad", None)
        if isinstance(h5ad, str) and h5ad != "nan" and Path(h5ad).exists():
            adata_vis = ad.read_h5ad(h5ad)
        else:
            input_dir = Path(row.input_dir)
            mtx = input_dir / "spatial_counts.mtx.gz"
            if not mtx.exists():
                mtx = input_dir / "spatial_counts.mtx"
            print(f"Reading spatial matrix for {sample}: {mtx}", flush=True)
            if str(mtx).endswith(".gz"):
                with gzip.open(mtx, "rb") as fh:
                    counts = io.mmread(fh).tocsr().T
            else:
                counts = io.mmread(mtx).tocsr().T
            obs = pd.read_csv(input_dir / "spatial_obs.tsv.gz", sep="\t", index_col=0)
            var = pd.read_csv(input_dir / "spatial_var.tsv.gz", sep="\t", index_col=0)
            adata_vis = ad.AnnData(X=counts, obs=obs, var=var)
        sc.pp.filter_genes(adata_vis, min_cells=1)
        adata_vis.var["SYMBOL"] = adata_vis.var_names

        common = adata_vis.var_names.intersection(signatures.index)
        adata_vis = adata_vis[:, common].copy()
        sig = signatures.loc[common, :].copy()
        if hasattr(adata_vis.X, "tocsr"):
            adata_vis.X = adata_vis.X.tocsr()
        print(f"{sample}: {adata_vis.n_obs} spots x {adata_vis.n_vars} genes, {sig.shape[1]} signatures", flush=True)

        Cell2location.setup_anndata(adata=adata_vis)
        mod = Cell2location(
            adata_vis,
            cell_state_df=sig,
            N_cells_per_location=args.N_cells_per_location,
            detection_alpha=args.detection_alpha,
        )
        train_kwargs = {"max_epochs": max_epochs, "batch_size": None, "train_size": 1}
        if use_gpu:
            train_kwargs.update({"accelerator": "gpu", "devices": 1})
        print(f"{sample}: training cell2location for {max_epochs} epochs; use_gpu={use_gpu}", flush=True)
        try:
            mod.train(**train_kwargs)
        except TypeError:
            train_kwargs.pop("accelerator", None)
            train_kwargs.pop("devices", None)
            mod.train(**train_kwargs)
        adata_vis = mod.export_posterior(
            adata_vis,
            sample_kwargs={"num_samples": 1000, "batch_size": None},
        )

        mod.save(str(sample_dir / "model"), overwrite=True)
        adata_vis.write_h5ad(sample_dir / "spatial_cell2location.h5ad")

        abund = adata_vis.obsm["q05_cell_abundance_w_sf"]
        abund = pd.DataFrame(abund, index=adata_vis.obs_names)
        abund.columns = [c.replace("q05cell_abundance_w_sf_", "") for c in abund.columns]
        abund.to_csv(sample_dir / "cell_abundance_q05.csv")

        mean_abund = adata_vis.obsm["means_cell_abundance_w_sf"]
        mean_abund = pd.DataFrame(mean_abund, index=adata_vis.obs_names)
        mean_abund.columns = [c.replace("meanscell_abundance_w_sf_", "") for c in mean_abund.columns]
        mean_abund.to_csv(sample_dir / "cell_abundance_mean.csv")

        obs = adata_vis.obs.copy()
        obs.to_csv(sample_dir / "spot_metadata.csv")
        summary_rows.append(
            {
                "sample": sample,
                "n_spots": adata_vis.n_obs,
                "n_genes": adata_vis.n_vars,
                "use_gpu": use_gpu,
                "max_epochs": max_epochs,
                "output": str(sample_dir / "cell_abundance_mean.csv"),
            }
        )

    summary_path = out_dir / "spatial_run_summary.csv"
    new_summary = pd.DataFrame(summary_rows)
    if summary_path.exists():
        old_summary = pd.read_csv(summary_path)
        merged_summary = pd.concat([old_summary, new_summary], ignore_index=True)
        merged_summary = merged_summary.drop_duplicates(subset=["sample"], keep="last")
    else:
        merged_summary = new_summary
    merged_summary.sort_values("sample").to_csv(summary_path, index=False)


if __name__ == "__main__":
    main()
