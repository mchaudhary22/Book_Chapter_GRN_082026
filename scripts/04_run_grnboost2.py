#!/usr/bin/env python
"""Run GRNBoost2 for a single cell type. Must be executed on a compute node
(submitted via sbatch, scripts/04_submit.sh) — never on a login node."""
import os

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"

import argparse
import time

import dask

dask.config.set({"dataframe.query-planning": False})

import numpy as np
import pandas as pd
import anndata as ad
from distributed import Client, LocalCluster
from arboreto.algo import grnboost2

PROCESSED = "data/processed"
RESULTS = "results"
SEED = 42


def slug(celltype):
    return celltype.replace(" ", "_").replace("/", "-")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--celltype", required=True,
                         help="Cell type to run GRNBoost2 on (must match adata.obs['celltype'] exactly)")
    parser.add_argument("--n-workers", type=int, default=64)
    parser.add_argument("--threads-per-worker", type=int, default=1)
    args = parser.parse_args()

    os.makedirs(RESULTS, exist_ok=True)

    print(f"Loading h5ad...", flush=True)
    adata = ad.read_h5ad(f"{PROCESSED}/root_singlets.h5ad")
    print(f"h5ad: {adata.shape[0]} cells x {adata.shape[1]} genes", flush=True)

    if args.celltype not in adata.obs["celltype"].values:
        raise SystemExit(f"FATAL: celltype '{args.celltype}' not found in adata.obs['celltype']")

    sub = adata[adata.obs["celltype"] == args.celltype].copy()
    n_cells = sub.shape[0]
    print(f"Subset to celltype '{args.celltype}': {n_cells} cells x {sub.shape[1]} genes", flush=True)

    with open(f"{PROCESSED}/regulators.txt") as f:
        regulators = [line.strip() for line in f if line.strip()]
    print(f"Loaded {len(regulators)} regulators", flush=True)

    regulators_in_data = [r for r in regulators if r in sub.var_names]
    print(f"Regulators present in data: {len(regulators_in_data)}", flush=True)
    if len(regulators_in_data) != len(regulators):
        missing = len(regulators) - len(regulators_in_data)
        print(f"WARNING: {missing} regulators missing from expression matrix var_names", flush=True)

    expr = pd.DataFrame(
        sub.X.toarray() if hasattr(sub.X, "toarray") else np.asarray(sub.X),
        index=sub.obs_names,
        columns=sub.var_names,
    )

    print(f"Starting dask LocalCluster: n_workers={args.n_workers}, "
          f"threads_per_worker={args.threads_per_worker}", flush=True)
    cluster = LocalCluster(
        n_workers=args.n_workers,
        threads_per_worker=args.threads_per_worker,
        processes=True,
    )
    client = Client(cluster)
    print(client, flush=True)

    t0 = time.time()
    network = grnboost2(
        expression_data=expr,
        tf_names=regulators_in_data,
        client_or_address=client,
        seed=SEED,
    )
    elapsed = time.time() - t0

    client.close()
    cluster.close()

    out_path = f"{RESULTS}/grnboost2_edges_{slug(args.celltype)}.csv"
    network.to_csv(out_path, index=False)

    print("\n=== SUMMARY ===", flush=True)
    print(f"Cell type: {args.celltype}", flush=True)
    print(f"Cells: {n_cells}", flush=True)
    print(f"Regulators used: {len(regulators_in_data)}", flush=True)
    print(f"Edges written: {len(network)}", flush=True)
    print(f"Runtime: {elapsed:.1f} s ({elapsed / 60:.2f} min)", flush=True)
    print(f"Output: {out_path}", flush=True)


if __name__ == "__main__":
    main()
