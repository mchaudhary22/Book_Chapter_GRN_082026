#!/usr/bin/env python
"""Step 06 figures: UMAP by celltype, TF feature plots, marker dotplot,
and per-celltype TF->target bipartite networks.

Implementation choices not fully pinned down by prompt.md, documented here
and in logs/run_log.md:
  - fig2: one combined figure, 3 rows (celltype) x 10 columns (TF rank 1-10),
    rather than 30 separate files, per "grouped by cell type" in the spec.
  - fig4: top 5 targets per TF (by importance into that celltype's top-50
    markers) shown per network, for legibility with 10 TFs x N targets.
  - TF-family colors are assigned from one global palette shared across all
    figures so the same family reads as the same color everywhere.
  - Any fig4 edge whose target is ALSO one of that celltype's own top-10 TFs
    (a TF-TF edge) is drawn in a distinct style and logged explicitly —
    this generically catches the SbiRTX430.09G017700 -> SbiRTX430.09G017600
    case (rank-1 TF -> rank-2 TF, r=0.56 raw-count co-expression in
    mature-endodermis) plus any other instance, without hardcoding gene IDs.
"""
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib.patches import FancyArrowPatch
from scipy.stats import pearsonr

PROCESSED = "data/processed"
RESULTS = "results"
RAW = "data/raw"
FIGURES = "figures"

CELLTYPES = ["mature-endodermis", "dividing cells", "exodermis"]
TOP_N_MARKERS = 50
TOP_TARGETS_PER_TF = 5
DPI = 300

plt.rcParams.update({
    "font.size": 9,
    "savefig.dpi": DPI,
    "figure.dpi": 100,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


def slug(ct):
    return ct.replace(" ", "_").replace("/", "-")


def save_fig(fig, name):
    fig.savefig(f"{FIGURES}/{name}.png", dpi=DPI, bbox_inches="tight")
    fig.savefig(f"{FIGURES}/{name}.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {FIGURES}/{name}.png and .pdf")


def build_celltype_palette(categories):
    n = len(categories)
    base = plt.get_cmap("tab20").colors + plt.get_cmap("tab20b").colors
    colors = base[:n] if n <= len(base) else [base[i % len(base)] for i in range(n)]
    return dict(zip(categories, colors))


def build_family_palette(families):
    n = len(families)
    cmap = plt.get_cmap("tab20")
    colors = [cmap(i / max(n - 1, 1)) for i in range(n)]
    return dict(zip(sorted(families), colors))


def fig1_umap_celltypes(adata):
    print("\n=== fig1: UMAP by celltype ===")
    cats = sorted(adata.obs["celltype"].cat.categories if hasattr(adata.obs["celltype"], "cat")
                  else adata.obs["celltype"].unique())
    palette = build_celltype_palette(cats)

    coords = adata.obsm["X_umap"]
    fig, ax = plt.subplots(figsize=(9, 7))
    for ct in cats:
        mask = (adata.obs["celltype"] == ct).values
        ax.scatter(coords[mask, 0], coords[mask, 1], s=3, alpha=0.7,
                   color=palette[ct], label=ct, linewidths=0)
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_title(f"Sorghum root atlas — {adata.shape[0]} singlet nuclei, {len(cats)} cell types")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    leg = ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), ncol=1,
                     fontsize=6.5, markerscale=3, frameon=False, title="celltype")
    save_fig(fig, "fig1_umap_celltypes")


def fig2_tf_featureplots(adata, top10_by_ct):
    print("\n=== fig2: TF feature plots, grouped by celltype ===")
    coords = adata.obsm["X_umap"]

    adata_norm = adata.copy()
    sc.pp.normalize_total(adata_norm, target_sum=1e4)
    sc.pp.log1p(adata_norm)

    n_rows = len(CELLTYPES)
    n_cols = 10
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(2.0 * n_cols, 2.0 * n_rows))

    for row, ct in enumerate(CELLTYPES):
        top10 = top10_by_ct[ct]
        for col in range(n_cols):
            ax = axes[row, col]
            if col >= len(top10):
                ax.axis("off")
                continue
            tf_row = top10.iloc[col]
            tf, family, rank = tf_row["TF"], tf_row["family"], int(tf_row["rank"])
            expr = np.asarray(adata_norm[:, tf].X.todense()).flatten()
            order = np.argsort(expr)  # plot high-expressers on top
            sca = ax.scatter(coords[order, 0], coords[order, 1], c=expr[order], s=2,
                              cmap="viridis", linewidths=0)
            ax.set_title(f"#{rank} {tf}\n({family})", fontsize=6.5)
            ax.set_xticks([])
            ax.set_yticks([])
            for spine in ax.spines.values():
                spine.set_visible(False)
            if col == 0:
                ax.set_ylabel(ct, fontsize=8, fontweight="bold")
            cb = fig.colorbar(sca, ax=ax, fraction=0.046, pad=0.02)
            cb.ax.tick_params(labelsize=4)

    fig.suptitle("Top-10 selected TFs per cell type — log-normalized expression on shared UMAP",
                 fontsize=11, y=1.01)
    fig.tight_layout()
    save_fig(fig, "fig2_tf_featureplots")


def fig3_marker_dotplot(adata):
    print("\n=== fig3: marker dotplot (Table S10) ===")
    markers = pd.read_excel(f"{RAW}/41477_2025_2047_MOESM3_ESM.xlsx",
                            sheet_name="TableS10.sorghum_lit_markers")
    n_raw = len(markers)
    markers = markers.drop_duplicates(subset="SbicolorRTx430")
    n_dedup = len(markers)
    present_mask = markers["SbicolorRTx430"].isin(adata.var_names)
    n_present = int(present_mask.sum())
    missing = markers.loc[~present_mask, "SbicolorRTx430"].tolist()
    print(f"Table S10 markers: {n_raw} rows -> {n_dedup} unique gene IDs -> "
          f"{n_present} present in the expression matrix.")
    if missing:
        print(f"NOT present in matrix (excluded from dotplot): {missing}")

    markers_present = markers[present_mask].copy()

    # Disambiguate aliases that are reused across different gene IDs (e.g. "SWEET13"
    # labels two distinct SbiRTX430 loci in this table) so gene_symbols lookups are unique.
    dup_alias_mask = markers_present["alias"].duplicated(keep=False)
    if dup_alias_mask.any():
        dup_aliases = sorted(markers_present.loc[dup_alias_mask, "alias"].unique())
        print(f"Alias values reused across multiple gene IDs (disambiguated with gene ID suffix): "
              f"{dup_aliases}")
    markers_present["display_label"] = markers_present["alias"]
    markers_present.loc[dup_alias_mask, "display_label"] = (
        markers_present.loc[dup_alias_mask, "alias"] + " ("
        + markers_present.loc[dup_alias_mask, "SbicolorRTx430"].str.split(".").str[-1] + ")"
    )
    assert markers_present["display_label"].is_unique, "marker display labels must be unique"

    adata_norm = adata.copy()
    sc.pp.normalize_total(adata_norm, target_sum=1e4)
    sc.pp.log1p(adata_norm)

    # gene_symbols expects var_names given in the SYMBOL namespace (the values of the
    # gene_symbols column), which it then maps back to the real adata.var_names for lookup.
    id_to_label = dict(zip(markers_present["SbicolorRTx430"], markers_present["display_label"]))
    adata_norm.var["marker_label"] = adata_norm.var_names.map(lambda g: id_to_label.get(g, g))
    display_labels = markers_present["display_label"].tolist()

    dp = sc.pl.dotplot(
        adata_norm,
        var_names=display_labels,
        groupby="celltype",
        gene_symbols="marker_label",
        standard_scale="var",
        dendrogram=False,
        show=False,
        figsize=(max(10, n_present * 0.28), 8),
        return_fig=True,
    )
    dp.make_figure()  # return_fig=True hands back the DotPlot object, not a rendered
                       # figure — must force rendering before .fig is populated.
    fig = dp.fig
    fig.suptitle(f"Table S10 literature root markers ({n_present}/{n_raw} present) "
                 f"across all 27 cell types", fontsize=10, y=1.02)
    save_fig(fig, "fig3_marker_dotplot")


def fig4_networks(adata, top10_by_ct, family_palette, edges_by_ct, markers_by_ct):
    print("\n=== fig4: TF-target bipartite networks (per celltype) ===")
    flagged_edges_log = []

    for ct in CELLTYPES:
        top10 = top10_by_ct[ct]
        top10_gene_set = set(top10["TF"])
        edges = edges_by_ct[ct]
        markers = markers_by_ct[ct]
        marker_edges = edges[edges["target"].isin(markers)]

        # top-N targets per TF, restricted to edges into the marker set
        rows = []
        for _, r in top10.iterrows():
            tf_edges = marker_edges[marker_edges["TF"] == r["TF"]].sort_values(
                "importance", ascending=False).head(TOP_TARGETS_PER_TF)
            rows.append(tf_edges)
        plot_edges = pd.concat(rows, ignore_index=True)
        plot_edges["target_is_top10_tf"] = plot_edges["target"].isin(top10_gene_set)

        n_flagged = int(plot_edges["target_is_top10_tf"].sum())
        if n_flagged:
            for _, r in plot_edges[plot_edges["target_is_top10_tf"]].iterrows():
                sub = adata[adata.obs["celltype"] == ct]
                x1 = np.asarray(sub[:, r["TF"]].X.todense()).flatten()
                x2 = np.asarray(sub[:, r["target"]].X.todense()).flatten()
                corr, _ = pearsonr(x1, x2)
                flagged_edges_log.append({
                    "celltype": ct, "TF": r["TF"], "target": r["target"],
                    "importance": r["importance"], "raw_count_pearson_r": corr,
                })
                print(f"FLAGGED (TF-TF edge): {ct}: {r['TF']} -> {r['target']} "
                      f"(target is also a top-10 TF), importance={r['importance']:.3f}, "
                      f"raw-count Pearson r={corr:.3f}")

        tf_order = top10["TF"].tolist()
        target_order = plot_edges["target"].drop_duplicates().tolist()
        tf_y = {tf: -i for i, tf in enumerate(tf_order)}
        tgt_y = {t: -i * (len(tf_order) / max(len(target_order), 1)) for i, t in enumerate(target_order)}

        fig, ax = plt.subplots(figsize=(9, max(6, 0.55 * len(tf_order))))
        max_imp = plot_edges["importance"].max()

        for _, r in plot_edges.iterrows():
            x0, y0 = 0, tf_y[r["TF"]]
            x1_, y1 = 1, tgt_y[r["target"]]
            lw = 0.5 + 4.0 * (r["importance"] / max_imp)
            if r["target_is_top10_tf"]:
                ax.plot([x0, x1_], [y0, y1], color="crimson", linewidth=lw + 1.0,
                        linestyle="--", zorder=5, alpha=0.9)
            else:
                ax.plot([x0, x1_], [y0, y1], color="grey", linewidth=lw, alpha=0.45, zorder=1)

        fam_series = top10.set_index("TF")["family"]
        for tf, y in tf_y.items():
            fam = fam_series[tf]
            ax.scatter([0], [y], s=180, color=family_palette[fam], edgecolor="black",
                       zorder=6, linewidths=0.8)
            rank = int(top10.set_index("TF").loc[tf, "rank"])
            ax.text(-0.03, y, f"#{rank} {tf} ({fam})", ha="right", va="center", fontsize=6.5)

        for t, y in tgt_y.items():
            is_tf = t in top10_gene_set
            ax.scatter([1], [y], s=60, color="crimson" if is_tf else "steelblue",
                       marker="*" if is_tf else "o", edgecolor="black", zorder=6, linewidths=0.5)
            label = t + ("  [also top-10 TF]" if is_tf else "")
            ax.text(1.03, y, label, ha="left", va="center", fontsize=5.5,
                    color="crimson" if is_tf else "black",
                    fontweight="bold" if is_tf else "normal")

        ax.set_xlim(-0.6, 1.7)
        ax.axis("off")
        ax.set_title(f"{ct}: top-10 TFs and top-{TOP_TARGETS_PER_TF} marker targets each\n"
                     f"edge width = importance; red dashed = target is itself a top-10 TF "
                     f"(co-expression inference limit, see legend note)", fontsize=8.5)

        fam_handles = [mlines.Line2D([], [], marker="o", linestyle="", markersize=9,
                                     markerfacecolor=family_palette[f], markeredgecolor="black",
                                     label=f) for f in sorted(top10["family"].unique())]
        edge_handles = [
            mlines.Line2D([], [], color="grey", linewidth=2, alpha=0.6, label="TF -> marker target"),
            mlines.Line2D([], [], color="crimson", linewidth=2, linestyle="--",
                         label="TF -> target that is also a top-10 TF"),
        ]
        leg1 = ax.legend(handles=fam_handles, title="TF family", loc="upper left",
                         bbox_to_anchor=(-0.55, 1.0), fontsize=6, title_fontsize=7, frameon=False)
        ax.add_artist(leg1)
        ax.legend(handles=edge_handles, loc="lower left", bbox_to_anchor=(-0.55, 0.0),
                 fontsize=6, frameon=False)

        save_fig(fig, f"fig4_networks_{slug(ct)}")

    if flagged_edges_log:
        flag_df = pd.DataFrame(flagged_edges_log)
        flag_df.to_csv(f"{FIGURES}/fig4_flagged_tf_tf_edges.csv", index=False)
        print(f"\nWrote {FIGURES}/fig4_flagged_tf_tf_edges.csv ({len(flag_df)} flagged edge(s))")
    else:
        print("\nNo TF-TF edges encountered in any fig4 network.")


def get_top_markers(adata_norm, ct):
    de = sc.get.rank_genes_groups_df(adata_norm, group=ct)
    de_sorted = de.sort_values("pvals_adj", ascending=True, kind="stable")
    return de_sorted.head(TOP_N_MARKERS)["names"].tolist()


def main():
    import os
    os.makedirs(FIGURES, exist_ok=True)

    print("Loading h5ad...")
    adata = ad.read_h5ad(f"{PROCESSED}/root_singlets.h5ad")
    print(f"h5ad: {adata.shape[0]} cells x {adata.shape[1]} genes")

    top10_by_ct = {}
    edges_by_ct = {}
    for ct in CELLTYPES:
        top10_by_ct[ct] = pd.read_csv(f"{RESULTS}/selected_tfs_{slug(ct)}.csv")
        edges_by_ct[ct] = pd.read_csv(f"{RESULTS}/grnboost2_edges_{slug(ct)}.csv")

    all_families = set()
    for ct in CELLTYPES:
        all_families |= set(top10_by_ct[ct]["family"])
    family_palette = build_family_palette(all_families)
    print(f"\nGlobal TF-family palette covers {len(all_families)} families "
          f"(shared across fig2/fig4): {sorted(all_families)}")

    adata_norm_for_markers = adata.copy()
    sc.pp.normalize_total(adata_norm_for_markers, target_sum=1e4)
    sc.pp.log1p(adata_norm_for_markers)
    sc.tl.rank_genes_groups(adata_norm_for_markers, groupby="celltype", method="wilcoxon", pts=False)
    markers_by_ct = {ct: get_top_markers(adata_norm_for_markers, ct) for ct in CELLTYPES}

    fig1_umap_celltypes(adata)
    fig2_tf_featureplots(adata, top10_by_ct)
    fig3_marker_dotplot(adata)
    fig4_networks(adata, top10_by_ct, family_palette, edges_by_ct, markers_by_ct)

    print("\nAll figures written to figures/ (PNG + PDF, 300 dpi).")


if __name__ == "__main__":
    main()
