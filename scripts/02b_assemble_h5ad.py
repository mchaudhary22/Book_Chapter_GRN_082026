#!/usr/bin/env python
"""Assemble data/processed/root_singlets.h5ad from the MatrixMarket export
produced by scripts/02_export_seurat.R (SeuratDisk route unavailable, see
Step 00 / prompt.md)."""
import scipy.io as sio
import pandas as pd
import anndata as ad

PROCESSED = "data/processed"


def main():
    print("Reading MatrixMarket counts...")
    counts = sio.mmread(f"{PROCESSED}/counts.mtx").tocsr()  # features x cells, as written by R

    genes = pd.read_csv(f"{PROCESSED}/genes.csv")["gene"].tolist()
    cells = pd.read_csv(f"{PROCESSED}/cells.csv")["cell"].astype(str).tolist()
    meta = pd.read_csv(f"{PROCESSED}/metadata.csv", index_col=0)
    meta.index = meta.index.astype(str)

    assert counts.shape == (len(genes), len(cells)), (
        f"counts.mtx shape {counts.shape} does not match "
        f"genes ({len(genes)}) x cells ({len(cells)})"
    )

    X = counts.T.tocsr()  # -> cells x genes

    meta = meta.loc[cells]

    var = pd.DataFrame(index=genes)
    adata = ad.AnnData(X=X, obs=meta, var=var)
    adata.obs_names = cells
    adata.var_names = genes

    umap_path = f"{PROCESSED}/umap_embedding.csv"
    try:
        umap = pd.read_csv(umap_path, index_col=0)
        umap.index = umap.index.astype(str)
        umap = umap.loc[cells]
        adata.obsm["X_umap"] = umap.values
        print(f"Attached X_umap, shape {umap.values.shape}")
    except FileNotFoundError:
        print(f"No UMAP embedding file at {umap_path}; skipping.")

    out_path = f"{PROCESSED}/root_singlets.h5ad"
    adata.write_h5ad(out_path)
    print(f"Wrote {out_path}: {adata.shape[0]} cells x {adata.shape[1]} genes")


if __name__ == "__main__":
    main()
