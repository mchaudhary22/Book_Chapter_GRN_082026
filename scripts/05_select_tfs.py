#!/usr/bin/env python
"""Select top TFs per target cell type from GRNBoost2 edges.

For each cell type: identify top-50 marker genes by Wilcoxon DE (celltype
vs rest, ranked by adjusted p-value), sum GRNBoost2 importance per TF across
edges into those markers, rank, take top 10.

Also characterizes GRNBoost2's exact-tie behavior (GBM feature_importances_
ties when a regulator dominates splits identically across several targets),
since large tie blocks affect how confidently the summed-importance ranking
and the Step 06 network figures can be read.
"""
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

PROCESSED = "data/processed"
RESULTS = "results"

CELLTYPES = ["mature-endodermis", "dividing cells", "exodermis"]
TOP_N_MARKERS = 50
TOP_N_TFS = 10
TIE_MIN_SIZE = 3
BULK_THRESHOLD = 0.5


def slug(celltype):
    return celltype.replace(" ", "_").replace("/", "-")


def main():
    print("Loading h5ad...")
    adata = ad.read_h5ad(f"{PROCESSED}/root_singlets.h5ad")
    print(f"h5ad: {adata.shape[0]} cells x {adata.shape[1]} genes")

    tf_fam = pd.read_csv(f"{PROCESSED}/tf_families.csv").set_index("gene_id")["family"].to_dict()

    # Normalize + log1p for DE testing only; raw counts on disk (h5ad, GRNBoost2
    # inputs) are untouched — this copy exists only for rank_genes_groups.
    adata_norm = adata.copy()
    sc.pp.normalize_total(adata_norm, target_sum=1e4)
    sc.pp.log1p(adata_norm)

    print("Running Wilcoxon rank_genes_groups (celltype vs rest, all 27 groups)...")
    sc.tl.rank_genes_groups(adata_norm, groupby="celltype", method="wilcoxon", pts=False)

    tie_summary_rows = []

    for ct in CELLTYPES:
        print(f"\n{'=' * 70}\n{ct}\n{'=' * 70}")

        de = sc.get.rank_genes_groups_df(adata_norm, group=ct)
        de_sorted = de.sort_values("pvals_adj", ascending=True, kind="stable")
        top_markers = de_sorted.head(TOP_N_MARKERS)["names"].tolist()
        print(f"Identified top {len(top_markers)} markers by adjusted p-value "
              f"(p_adj range: {de_sorted['pvals_adj'].iloc[0]:.3g} - "
              f"{de_sorted['pvals_adj'].iloc[TOP_N_MARKERS - 1]:.3g})")

        edge_path = f"{RESULTS}/grnboost2_edges_{slug(ct)}.csv"
        edges = pd.read_csv(edge_path)
        n_total_edges = len(edges)

        # Exact-tie blocks: >=3 edges from the same TF sharing an identical
        # importance value, computed over the FULL edge set for this cell type.
        grp_sizes = edges.groupby(["TF", "importance"])["target"].transform("size")
        edges["in_tie_block"] = grp_sizes >= TIE_MIN_SIZE
        n_tie_edges = int(edges["in_tie_block"].sum())
        frac_tie = n_tie_edges / n_total_edges
        print(f"\nTie-block edges (>= {TIE_MIN_SIZE} edges/TF at identical importance): "
              f"{n_tie_edges} / {n_total_edges} = {frac_tie:.4%}")

        marker_edges = edges[edges["target"].isin(top_markers)].copy()
        print(f"Edges from regulators into the top-{TOP_N_MARKERS} markers: {len(marker_edges)}")

        tf_sum = marker_edges.groupby("TF")["importance"].sum().sort_values(ascending=False)
        tf_tie_sum = marker_edges[marker_edges["in_tie_block"]].groupby("TF")["importance"].sum()

        top10 = tf_sum.head(TOP_N_TFS).reset_index()
        top10.columns = ["TF", "summed_importance"]
        top10["rank"] = np.arange(1, len(top10) + 1)
        top10["family"] = top10["TF"].map(tf_fam).fillna("NA")
        top10["tie_block_importance"] = top10["TF"].map(tf_tie_sum).fillna(0.0)
        top10["frac_from_tie_blocks"] = top10["tie_block_importance"] / top10["summed_importance"]

        out_path = f"{RESULTS}/selected_tfs_{slug(ct)}.csv"
        top10[["TF", "summed_importance", "rank", "family"]].to_csv(out_path, index=False)
        print(f"Wrote {out_path}")

        print(f"\nTop {TOP_N_TFS} TFs for '{ct}':")
        print(top10[["rank", "TF", "family", "summed_importance", "frac_from_tie_blocks"]]
              .to_string(index=False, float_format=lambda x: f"{x:.4f}"))

        bulk_flag = top10[top10["frac_from_tie_blocks"] >= BULK_THRESHOLD]
        if len(bulk_flag) > 0:
            print(f"\nFLAG: {len(bulk_flag)} of top {TOP_N_TFS} TFs derive "
                  f">={BULK_THRESHOLD:.0%} of summed importance from tie blocks: "
                  f"{bulk_flag['TF'].tolist()}")
        else:
            print(f"\nNo top-{TOP_N_TFS} TF derives >={BULK_THRESHOLD:.0%} of its "
                  f"summed importance from tie blocks.")

        tie_summary_rows.append({
            "celltype": ct,
            "n_total_edges": n_total_edges,
            "n_tie_edges": n_tie_edges,
            "frac_tie_edges": frac_tie,
            "n_top10_bulk_from_ties": len(bulk_flag),
            "top10_bulk_TFs": ";".join(bulk_flag["TF"].tolist()),
        })

    tie_summary_df = pd.DataFrame(tie_summary_rows)
    tie_summary_df.to_csv(f"{RESULTS}/tie_block_summary.csv", index=False)
    print(f"\n\n{'=' * 70}\nTIE-BLOCK SUMMARY (all cell types)\n{'=' * 70}")
    print(tie_summary_df.to_string(index=False))
    print(f"\nWrote {RESULTS}/tie_block_summary.csv")


if __name__ == "__main__":
    main()
