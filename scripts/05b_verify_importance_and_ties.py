#!/usr/bin/env python
"""Verification/audit script for Step 05 findings:
(1) confirm frac_from_tie_blocks == 0 for top-10 TFs is real, not a join bug,
    by auditing one TF's individual edges and a ranks-11-30 control;
(2) empirically check arboreto's per-target importance denormalization;
(3) locus-adjacency / co-expression check for the two mature-endodermis
    rank-1/rank-2 MYB TFs.

Not a new pipeline stage — diagnostic companion to scripts/05_select_tfs.py.
"""
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy.stats import pearsonr, spearmanr

PROCESSED = "data/processed"
RESULTS = "results"

CELLTYPES = ["mature-endodermis", "dividing cells", "exodermis"]
TOP_N_MARKERS = 50
TIE_MIN_SIZE = 3


def slug(ct):
    return ct.replace(" ", "_").replace("/", "-")


def main():
    print("Loading h5ad...")
    adata = ad.read_h5ad(f"{PROCESSED}/root_singlets.h5ad")

    adata_norm = adata.copy()
    sc.pp.normalize_total(adata_norm, target_sum=1e4)
    sc.pp.log1p(adata_norm)
    print("Recomputing Wilcoxon rank_genes_groups (must match Step 05 exactly)...")
    sc.tl.rank_genes_groups(adata_norm, groupby="celltype", method="wilcoxon", pts=False)

    all_ranked = {}
    for ct in CELLTYPES:
        de = sc.get.rank_genes_groups_df(adata_norm, group=ct)
        de_sorted = de.sort_values("pvals_adj", ascending=True, kind="stable")
        top_markers = de_sorted.head(TOP_N_MARKERS)["names"].tolist()

        edges = pd.read_csv(f"{RESULTS}/grnboost2_edges_{slug(ct)}.csv")
        grp_sizes = edges.groupby(["TF", "importance"])["target"].transform("size")
        edges["in_tie_block"] = grp_sizes >= TIE_MIN_SIZE

        marker_edges = edges[edges["target"].isin(top_markers)].copy()
        tf_sum = marker_edges.groupby("TF")["importance"].sum().sort_values(ascending=False)
        tf_tie_sum = marker_edges[marker_edges["in_tie_block"]].groupby("TF")["importance"].sum()

        ranked = tf_sum.reset_index()
        ranked.columns = ["TF", "summed_importance"]
        ranked["rank"] = np.arange(1, len(ranked) + 1)
        ranked["tie_block_importance"] = ranked["TF"].map(tf_tie_sum).fillna(0.0)
        ranked["frac_from_tie_blocks"] = ranked["tie_block_importance"] / ranked["summed_importance"]

        all_ranked[ct] = (ranked, marker_edges, edges, top_markers)

        print(f"\n{'=' * 70}\n{ct}: CONTROL — frac_from_tie_blocks for ranks 11-30\n{'=' * 70}")
        control = ranked.iloc[10:30][["rank", "TF", "summed_importance", "frac_from_tie_blocks"]]
        print(control.to_string(index=False))
        n_nonzero = (control["frac_from_tie_blocks"] > 0).sum()
        print(f"-> {n_nonzero}/20 of ranks 11-30 have nonzero frac_from_tie_blocks "
              f"(nonzero values expected somewhere if the metric is working)")

    # ---- (1) Detailed audit: rank-1 TF in mature-endodermis ----
    ct = "mature-endodermis"
    ranked, marker_edges, edges, top_markers = all_ranked[ct]
    top_tf = ranked.iloc[0]["TF"]
    print(f"\n{'=' * 70}\nDETAILED AUDIT: rank-1 TF {top_tf} in {ct}\n{'=' * 70}")

    tf_marker_edges = marker_edges[marker_edges["TF"] == top_tf].sort_values("importance", ascending=False)
    print(f"\nAll {len(tf_marker_edges)} edges from {top_tf} into the top-{TOP_N_MARKERS} markers:")
    print(tf_marker_edges[["target", "importance", "in_tie_block"]].to_string(index=False))

    dup_check = tf_marker_edges["importance"].value_counts()
    dup_check = dup_check[dup_check >= 2]
    print(f"\nImportance values repeated 2+ times among THESE {len(tf_marker_edges)} marker-edges: "
          f"{len(dup_check)} such values" + (f" -> {dup_check.to_dict()}" if len(dup_check) else " (none)"))

    tf_all_edges = edges[edges["TF"] == top_tf]
    n_tie_total = int(tf_all_edges["in_tie_block"].sum())
    print(f"\nAcross {top_tf}'s FULL network ({len(tf_all_edges)} edges to all targets, not just markers): "
          f"{n_tie_total} edges ({n_tie_total / len(tf_all_edges):.2%}) are in tie blocks elsewhere.")
    tie_edges_this_tf = tf_all_edges[tf_all_edges["in_tie_block"]]
    if len(tie_edges_this_tf) > 0:
        example_val = tie_edges_this_tf["importance"].iloc[0]
        example_group = tf_all_edges[tf_all_edges["importance"] == example_val]
        print(f"Example tie block for {top_tf}: importance={example_val!r}, "
              f"{len(example_group)} tied targets (first 10): {example_group['target'].tolist()[:10]}")
        # Are these tied targets low-expression / sparse? (candidate mechanism)
        sub = adata[adata.obs["celltype"] == ct]
        tied_targets = example_group["target"].tolist()
        present = [t for t in tied_targets if t in sub.var_names]
        if present:
            X = sub[:, present].X
            detect_rate = np.asarray((X > 0).sum(axis=0)).flatten() / sub.shape[0]
            mean_expr = np.asarray(X.mean(axis=0)).flatten()
            print("Detection rate / mean raw count for those tied targets (checking if low-expression "
                  "targets drive ties):")
            for t, d, m in zip(present, detect_rate, mean_expr):
                print(f"  {t}: detected in {d:.1%} of nuclei, mean count {m:.4f}")
        # for context, same stats for the top (non-tied) marker edge
        top_target = tf_marker_edges.iloc[0]["target"]
        if top_target in sub.var_names:
            Xtop = sub[:, [top_target]].X
            d = float(np.asarray((Xtop > 0).sum(axis=0)).flatten()[0]) / sub.shape[0]
            m = float(np.asarray(Xtop.mean(axis=0)).flatten()[0])
            print(f"For comparison, {top_tf}'s TOP (non-tied) marker edge target {top_target}: "
                  f"detected in {d:.1%} of nuclei, mean count {m:.4f}")

    # ---- (2) arboreto denormalization: empirical check ----
    print(f"\n{'=' * 70}\nEMPIRICAL CHECK: per-target total importance (sum over ALL TFs)\n"
          f"arboreto (core.py) uses regressor_type='GBM', SGBM_KWARGS "
          f"(subsample=0.9 < 1.0) -> is_oob_heuristic_supported=True ->\n"
          f"denormalized_importances = trained_regressor.feature_importances_ * n_estimators_used.\n"
          f"sklearn's feature_importances_ sums to 1.0 across predictors for that target's model,\n"
          f"so the per-target TOTAL importance (summed over all TFs) should approx. equal the\n"
          f"number of boosting rounds actually used for that target before early stopping\n"
          f"(an integer-ish value, bounded by n_estimators=5000, and varying target to target).\n{'=' * 70}")
    sample_targets = top_markers[:15]
    per_target_sum = edges[edges["target"].isin(sample_targets)].groupby("target")["importance"].sum()
    per_target_sum = per_target_sum.reindex(sample_targets)
    print("\nTotal importance summed over ALL ~1022 TFs, for the first 15 top-50 markers:")
    print(per_target_sum.to_string())
    print(f"\nRange across these 15 targets: {per_target_sum.min():.2f} - {per_target_sum.max():.2f} "
          f"(mean {per_target_sum.mean():.2f})")
    all_target_sums = edges.groupby("target")["importance"].sum()
    print(f"\nAcross ALL {len(all_target_sums)} targets in this network: "
          f"min={all_target_sums.min():.2f}, max={all_target_sums.max():.2f}, "
          f"mean={all_target_sums.mean():.2f}, median={all_target_sums.median():.2f}")

    # ---- (3) SbiRTX430.09G017700 vs SbiRTX430.09G017600 ----
    print(f"\n{'=' * 70}\nSbiRTX430.09G017700 vs SbiRTX430.09G017600 (mature-endodermis)\n{'=' * 70}")
    g1, g2 = "SbiRTX430.09G017700", "SbiRTX430.09G017600"
    print(f"Locus IDs differ by exactly 100 on chromosome 09 ({g1} vs {g2}) — consistent with "
          f"immediately adjacent genes under this annotation's locus-numbering convention "
          f"(neighboring genes typically increment by 100). No GFF/annotation file is available in "
          f"data/raw to directly confirm physical bp distance or gene orientation — this is an "
          f"inference from ID adjacency, not a confirmed coordinate-based tandem-duplicate call.")

    sub = adata[adata.obs["celltype"] == "mature-endodermis"]
    x1 = np.asarray(sub[:, g1].X.todense()).flatten()
    x2 = np.asarray(sub[:, g2].X.todense()).flatten()
    pr, pp = pearsonr(x1, x2)
    sr, sp = spearmanr(x1, x2)
    print(f"\nRaw counts, {sub.shape[0]} mature-endodermis nuclei:")
    print(f"  Pearson r = {pr:.4f} (p = {pp:.3g})")
    print(f"  Spearman rho = {sr:.4f} (p = {sp:.3g})")
    print(f"  {g1}: mean={x1.mean():.3f}, detected in {int((x1 > 0).sum())}/{len(x1)} nuclei "
          f"({(x1 > 0).mean():.1%})")
    print(f"  {g2}: mean={x2.mean():.3f}, detected in {int((x2 > 0).sum())}/{len(x2)} nuclei "
          f"({(x2 > 0).mean():.1%})")

    subn = adata_norm[adata_norm.obs["celltype"] == "mature-endodermis"]
    y1 = np.asarray(subn[:, g1].X.todense()).flatten()
    y2 = np.asarray(subn[:, g2].X.todense()).flatten()
    pr2, pp2 = pearsonr(y1, y2)
    print(f"\nLog-normalized (normalize_total 1e4 + log1p):")
    print(f"  Pearson r = {pr2:.4f} (p = {pp2:.3g})")


if __name__ == "__main__":
    main()
