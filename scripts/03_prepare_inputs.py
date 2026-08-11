#!/usr/bin/env python
"""Prepare the regulator (TF) list for GRNBoost2 from the exported h5ad
and the raw TF annotation file."""
import numpy as np
import pandas as pd
import anndata as ad

PROCESSED = "data/processed"
RAW = "data/raw"

EXPECTED_TF_ROWS = 1827
EXPECTED_PRESENT = 1370
EXPECTED_THRESHOLDED = 1022
DETECTION_MIN_NUCLEI = 50


def main():
    print("Reading h5ad...")
    adata = ad.read_h5ad(f"{PROCESSED}/root_singlets.h5ad")
    print(f"h5ad: {adata.shape[0]} cells x {adata.shape[1]} genes")

    print("\nReading TF list...")
    tf_df = pd.read_csv(
        f"{RAW}/Sbi_RTx430_TF_list.txt",
        sep="\t",
        header=None,
        names=["protein_id", "family"],
    )
    print(f"TF list rows: {len(tf_df)} (expected {EXPECTED_TF_ROWS})")
    if len(tf_df) != EXPECTED_TF_ROWS:
        raise SystemExit(f"FATAL: TF list row count {len(tf_df)} != expected {EXPECTED_TF_ROWS}")

    tf_df["gene_id"] = tf_df["protein_id"].str.replace(r"\.1\.p$", "", regex=True)
    n_dup = tf_df["gene_id"].duplicated().sum()
    print(f"Rows: {len(tf_df)}, unique gene_ids after stripping '.1.p': {tf_df['gene_id'].nunique()} "
          f"(duplicates collapsed: {n_dup})")

    # Filter to TFs present in the expression matrix
    present_mask = tf_df["gene_id"].isin(adata.var_names)
    tf_present = tf_df[present_mask].drop_duplicates(subset="gene_id").copy()
    n_present = len(tf_present)
    print(f"\nTFs present in matrix: {n_present} (expected {EXPECTED_PRESENT})")
    if n_present != EXPECTED_PRESENT:
        raise SystemExit(f"FATAL: present TF count {n_present} != expected {EXPECTED_PRESENT}")

    # Detection threshold: nonzero raw count in >= 50 nuclei
    sub = adata[:, tf_present["gene_id"].values]
    n_detected = np.asarray(sub.X.getnnz(axis=0)).flatten()
    tf_present["n_nuclei_detected"] = n_detected

    tf_thresholded = tf_present[tf_present["n_nuclei_detected"] >= DETECTION_MIN_NUCLEI].copy()
    n_thresholded = len(tf_thresholded)
    print(f"TFs detected in >={DETECTION_MIN_NUCLEI} nuclei: {n_thresholded} (expected {EXPECTED_THRESHOLDED})")
    if n_thresholded != EXPECTED_THRESHOLDED:
        raise SystemExit(f"FATAL: thresholded TF count {n_thresholded} != expected {EXPECTED_THRESHOLDED}")

    with open(f"{PROCESSED}/regulators.txt", "w") as f:
        for gid in tf_thresholded["gene_id"]:
            f.write(gid + "\n")
    print(f"\nWrote {PROCESSED}/regulators.txt ({n_thresholded} genes)")

    tf_thresholded[["gene_id", "family", "n_nuclei_detected"]].sort_values("gene_id").to_csv(
        f"{PROCESSED}/tf_families.csv", index=False
    )
    print(f"Wrote {PROCESSED}/tf_families.csv")

    print("\nFamily breakdown (post-threshold, sorted by count):")
    fam_counts = tf_thresholded["family"].value_counts()
    print(fam_counts.to_string())

    print("\n=== SUMMARY ===")
    print(f"TF list rows:                 {len(tf_df)}")
    print(f"Present in matrix:            {n_present}")
    print(f"After >={DETECTION_MIN_NUCLEI} nuclei threshold: {n_thresholded}")
    print(f"Distinct TF families (final): {tf_thresholded['family'].nunique()}")


if __name__ == "__main__":
    main()
