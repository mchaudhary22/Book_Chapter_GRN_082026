#!/usr/bin/env python
"""Drill into one specific GRNBoost2 tie block to test whether identical
importance values arise from identical (or proportional) target expression
vectors, vs. being a more general GBM behavior."""
import numpy as np
import pandas as pd
import anndata as ad

PROCESSED = "data/processed"
RESULTS = "results"

TF = "SbiRTX430.09G017700"
CELLTYPE = "mature-endodermis"
IMPORTANCE_VALUE = 3.156694826102661


def main():
    edges = pd.read_csv(f"{RESULTS}/grnboost2_edges_mature-endodermis.csv")
    tie_group = edges[(edges["TF"] == TF) & (edges["importance"] == IMPORTANCE_VALUE)]
    targets = tie_group["target"].tolist()
    print(f"Tie group: TF={TF}, importance={IMPORTANCE_VALUE}, n_targets={len(targets)}")
    print(targets)

    adata = ad.read_h5ad(f"{PROCESSED}/root_singlets.h5ad")
    sub = adata[adata.obs["celltype"] == CELLTYPE]
    n_nuclei = sub.shape[0]
    print(f"\n{CELLTYPE}: {n_nuclei} nuclei")

    X = sub[:, targets].X
    X = np.asarray(X.todense()) if hasattr(X, "todense") else np.asarray(X)

    print(f"\nPer-target detection / count summary ({len(targets)} tied targets):")
    detected_nuclei_sets = []
    for i, t in enumerate(targets):
        col = X[:, i]
        nz_idx = np.nonzero(col)[0]
        detected_nuclei_sets.append(frozenset(nz_idx.tolist()))
        vals = col[nz_idx]
        print(f"  {t}: detected in {len(nz_idx)}/{n_nuclei} nuclei, "
              f"nonzero values={vals.tolist()}, nuclei_idx={nz_idx.tolist()}")

    unique_patterns = set(detected_nuclei_sets)
    print(f"\nNumber of DISTINCT detection patterns (which nuclei are nonzero) among "
          f"the {len(targets)} tied targets: {len(unique_patterns)}")
    if len(unique_patterns) == 1:
        print("-> ALL 19 tied targets are detected in exactly the SAME nucleus/nuclei.")

    # exact vector identity check
    identical_to_first = np.all(X == X[:, [0]], axis=0)
    print(f"\nColumns exactly identical to target[0] ({targets[0]})'s vector: "
          f"{identical_to_first.sum()}/{len(targets)}")

    # proportionality check (in case counts differ but pattern is scaled)
    ref = X[:, 0].astype(float)
    ref_nonzero = ref[ref != 0]
    proportional = []
    for i in range(1, len(targets)):
        col = X[:, i].astype(float)
        if np.array_equal(np.nonzero(col)[0], np.nonzero(ref)[0]) and len(ref_nonzero) > 0:
            ratios = col[col != 0] / ref[col != 0]
            proportional.append(np.allclose(ratios, ratios[0]))
        else:
            proportional.append(False)
    print(f"Columns proportional to target[0] (same nonzero support, constant ratio): "
          f"{sum(proportional)}/{len(targets)-1}")

    # full matrix summary stats
    print(f"\nFull tie-group submatrix: shape={X.shape}, total nonzero entries={np.count_nonzero(X)}, "
          f"sum={X.sum()}, max value={X.max()}")


if __name__ == "__main__":
    main()
