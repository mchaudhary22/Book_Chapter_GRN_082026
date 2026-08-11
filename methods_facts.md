# Methods facts — sorghum root atlas GRN pipeline

Consolidated reference for Materials and Methods drafting. Source: `logs/run_log.md`
(full narrative log) and `prompt.md` (pipeline specification). All work done
2026-08-10 on VT ARC TinkerCliffs.

---

## 1. Compute environment

- **Cluster:** VT ARC TinkerCliffs (SLURM). All compute run on compute nodes, never login nodes.
- **R:** `R/4.4.2-gfbf-2024a` + `R-bundle-Bioconductor/3.20-foss-2024a-R-4.4.2` (module system).
  - Seurat 5.3.1, SeuratObject 5.3.0, hdf5r 1.3.11 (pre-installed in module; no personal-library install needed).
  - SeuratDisk: not available in any ARC module; not installed (user decision). Export used MatrixMarket + Python assembly route instead.
- **Python / conda environment `grn`** (Python 3.10.0), built via `mamba create -p ~/.conda/envs/grn`, exported to `environment.yml` (project root, full pinned dependency list, `conda env export --no-builds`):
  - pyscenic 0.12.1, scanpy 1.10.4, anndata 0.11.4, pandas 2.3.3, numpy 1.23.5, matplotlib 3.10.9, seaborn 0.13.2, openpyxl 3.1.5, setuptools 80.10.2
  - arboreto 0.1.6, ctxcore 0.2.0, dask 2024.8.2, distributed 2024.8.2 (transitive dependencies of pyscenic/GRNBoost2)
- **Known dependency fixes applied** (per pipeline spec, all three verified working via synthetic smoke test before production use):
  1. `setuptools<81` pinned (pyscenic 0.12.1 / ctxcore require `pkg_resources`, removed in newer setuptools).
  2. `dask.config.set({"dataframe.query-planning": False})` called before importing arboreto/distributed (arboreto 0.1.6 incompatible with dask's query-planning backend).
  3. `OMP_NUM_THREADS`, `MKL_NUM_THREADS`, `OPENBLAS_NUM_THREADS` set to `1` before importing numpy/sklearn; dask `LocalCluster` constructed explicitly with one thread per worker (`threads_per_worker=1`).
- **Synthetic smoke test:** GRNBoost2 on 50 cells × 20 genes, 5 TFs, seed 42 → **95 edges** (matches expectation exactly), confirming the environment and dependency fixes before production runs.

---

## 2. Input data

- **Seurat object:** `GSE297576_seurat_object.bicolor_root_atlas.RDS`, 1,314,285,536 bytes (confirmed by direct file inspection; note if citing GEO's listing, its rounded "1.2 Gb" display corresponds to this same file — the two are consistent, not a discrepancy in the underlying data).
- **Assays in the RDS:** `RNA` (25,464 features, raw counts) and `SCT` (23,974 features). Object's `DefaultAssay` as loaded is `SCT` — RNA assay selected explicitly for all downstream analysis (raw counts, not SCT-normalized, not `scale.data`).
- **Reductions present in the RDS:** `pca`, `umap` (UMAP embedding used directly from the original object, not recomputed).
- **Doublet/singlet column:** metadata column `doublet` (auto-detected by scanning for the literal value `"Singlet"`; not assumed by name). Values: Doublet = 4,524, Singlet = 14,792, NA = 0.
- **Cell-type column:** `celltype`, 27 levels. (`cluster_label.long` explicitly ignored per spec — leaf-atlas artifact with NAs.)
- **TF list:** `Sbi_RTx430_TF_list.txt`, 1,827 rows, tab-separated (protein ID, TF family), no header. Protein IDs stripped of trailing `.1.p` to yield gene IDs (1:1 mapping, 0 duplicates after stripping).
- **Marker table:** `41477_2025_2047_MOESM3_ESM.xlsx`, sheet `TableS10.sorghum_lit_markers`, 50 rows, columns include `SbicolorRTx430` (gene ID used), `alias` (short gene name), `SbicolorBTx623`.

---

## 3. Pipeline parameters and per-step results

### Step 02 — Seurat → h5ad export (`scripts/02_export_seurat.R`, `scripts/02b_assemble_h5ad.py`)

| Quantity | Value |
|---|---|
| Cells before singlet filter | 19,316 |
| RNA features (before and after filter) | 25,464 |
| Cells after singlet filter | **14,792** |
| Route | MatrixMarket export from R + Python (scanpy/anndata) assembly, not SeuratDisk |
| Output | `data/processed/root_singlets.h5ad`, verified cells × genes orientation (14792, 25464), raw int64 counts (max value 525), `X_umap` carried through in `obsm` |

**Per-cell-type counts (all 27, singlets only):**

| celltype | n | celltype | n | celltype | n |
|---|---|---|---|---|---|
| cortex | 1589 | early-epidermis | 599 | mature-endodermis | 282 |
| pericycle | 1362 | cortical sclerenchyma | 575 | cortical aerenchyma | 281 |
| meristem | 1354 | phloem | 511 | phloem/SE | 264 |
| atrichoblast | 1028 | trichoblast | 504 | elongating-endodermis | 202 |
| lignified-xylem | 818 | LRC | 470 | mature-xylem | 189 |
| dividing cells | 801 | early-cortex | 438 | sclerenchyma | 151 |
| elongating-xylem | 684 | early-stele | 412 | lignified-endodermis | 133 |
| s-phase | 633 | crown root initials | 346 | senescent cortex | 129 |
| procambium | 605 | exodermis | 315 | columella | 117 |

### Step 03 — Regulator preparation (`scripts/03_prepare_inputs.py`)

| Stage | Count |
|---|---|
| TF list (raw) | 1,827 |
| Present in expression matrix | 1,370 |
| Detected in ≥50 nuclei (nonzero raw count) | **1,022** — used as the final regulator set |
| Distinct TF families (post-threshold) | 54 (largest: bHLH 96, MYB 73, bZIP 72, NAC 60, ERF 57, WRKY 52, C2H2 51) |

### Step 04 — GRNBoost2 (`scripts/04_run_grnboost2.py`, `scripts/04_submit.sh`)

Common parameters: GRNBoost2 (arboreto 0.1.6), regulators = the 1,022-gene set above, targets = all genes, **seed = 42**, dask `n_workers=64`, `threads_per_worker=1`, SLURM partition `normal_q`, QOS `tc_normal_short`, account `introtogds`, 1 node / 64 cores / 64 GB per job.

| Cell type | SLURM job ID | Nuclei | Wall time requested | Elapsed | Peak mem (MaxRSS) | Edges | Importance (min–max, mean, median) |
|---|---|---|---|---|---|---|---|
| mature-endodermis (probe) | 6963109 | 282 | 4 h | 00:06:56 (382.0 s) | 16,225,436 KB ≈ 15.5 GB | 1,246,780 | 5.5e-18 – 14.88, mean 0.403, median 0.205 |
| dividing cells | 6963283 | 801 | 8 h | 00:05:29 (290.6 s) | 16,781,100 KB ≈ 16.0 GB | 1,563,314 | 1.4e-18 – 15.20, mean 0.366, median 0.200 |
| exodermis | 6963284 | 315 | 4 h | 00:06:24 (355.0 s) | 16,109,460 KB ≈ 15.4 GB | 1,380,256 | 1.2e-18 – 14.19, mean 0.387, median 0.201 |

All three: 1,022/1,022 regulators present in the subsetted expression data. Runtime is not proportional to nuclei count at this scale (dividing cells, ~3× the probe's nuclei, ran fastest) — dominated by dask/GBM scheduling overhead on the fixed regulator×target grid rather than cell count. Note the two non-probe jobs were submitted concurrently rather than sized from the probe's measured runtime (user decision); wall time was instead scaled from nuclei count (dividing cells got 8 h as a ~3× margin, exodermis kept 4 h).

### Step 05 — TF selection (`scripts/05_select_tfs.py`)

- Marker genes: Wilcoxon differential expression (`sc.tl.rank_genes_groups`, method="wilcoxon", `celltype` vs. rest, computed across all 27 groups in one call) on a `normalize_total(target_sum=1e4)` + `log1p` **copy** of the AnnData (raw counts on disk untouched; normalization applied for the DE test only). Top 50 markers per cell type by ascending adjusted p-value.
- Per TF: GRNBoost2 importance summed across edges into the 50 markers; ranked descending; top 10 retained.
- Outputs: `results/selected_tfs_<celltype>.csv` (columns: TF, summed_importance, rank, family).

**Top-10 TFs per cell type (rank, TF, family, summed importance):**

*mature-endodermis:*
| Rank | TF | Family | Summed importance |
|---|---|---|---|
| 1 | SbiRTX430.09G017700 | MYB | 74.81 |
| 2 | SbiRTX430.09G017600 | MYB | 25.02 |
| 3 | SbiRTX430.04G243200 | ERF | 21.80 |
| 4 | SbiRTX430.02G238700 | HD-ZIP | 20.28 |
| 5 | SbiRTX430.01G471500 | HD-ZIP | 19.72 |
| 6 | SbiRTX430.06G011000 | C2H2 | 18.89 |
| 7 | SbiRTX430.01G110400 | bHLH | 17.16 |
| 8 | SbiRTX430.10G016400 | AP2 | 15.06 |
| 9 | SbiRTX430.04G299700 | Dof | 13.26 |
| 10 | SbiRTX430.09G149400 | NAC | 12.24 |

*dividing cells:*
| Rank | TF | Family | Summed importance |
|---|---|---|---|
| 1 | SbiRTX430.10G018700 | GRAS | 73.84 |
| 2 | SbiRTX430.09G219200 | WRKY | 54.56 |
| 3 | SbiRTX430.09G224300 | C2H2 | 39.35 |
| 4 | SbiRTX430.06G052500 | FAR1 | 35.32 |
| 5 | SbiRTX430.03G428400 | C2H2 | 20.26 |
| 6 | SbiRTX430.09G130300 | AP2 | 19.06 |
| 7 | SbiRTX430.07G054700 | MYB_related | 17.96 |
| 8 | SbiRTX430.09G085600 | MYB | 17.63 |
| 9 | SbiRTX430.01G110400 | bHLH | 15.15 |
| 10 | SbiRTX430.10G044200 | DBB | 11.42 |

*exodermis:*
| Rank | TF | Family | Summed importance |
|---|---|---|---|
| 1 | SbiRTX430.03G038000 | NAC | 74.51 |
| 2 | SbiRTX430.02G190100 | ERF | 27.78 |
| 3 | SbiRTX430.06G011000 | C2H2 | 26.64 |
| 4 | SbiRTX430.04G294200 | MYB_related | 25.09 |
| 5 | SbiRTX430.09G017700 | MYB | 24.56 |
| 6 | SbiRTX430.02G211900 | AP2 | 23.62 |
| 7 | SbiRTX430.04G243200 | ERF | 19.41 |
| 8 | SbiRTX430.08G120400 | WRKY | 17.26 |
| 9 | SbiRTX430.03G310600 | WRKY | 14.40 |
| 10 | SbiRTX430.09G219200 | WRKY | 14.30 |

**Rank-1 vs rank-2 ratio:** mature-endodermis 2.99×, exodermis 2.68×, dividing cells 1.35× — the tight rank-1 agreement across cell types (~74–75) is specific to the top rank; ranks 2–10 diverge substantially between cell types (dividing cells decays much more gradually, reflecting stronger secondary/tertiary regulators), converging again by rank 8–10.

### Step 06 — Figures (`scripts/06_figures.py`)

All at 300 dpi, PNG + PDF, in `figures/`:
- `fig1_umap_celltypes` — UMAP colored by celltype, all 27 types, legend outside plot area.
- `fig2_tf_featureplots` — one combined figure, 3 rows (cell type) × 10 columns (TF rank 1–10) = 30 panels, log-normalized expression, shared UMAP coordinates.
- `fig3_marker_dotplot` — Table S10 markers (46/50 present in the matrix; see §4) across all 27 cell types, `standard_scale="var"`, gene labels from the `alias` column (one true duplicate alias, "SWEET13", disambiguated with a gene-ID suffix).
- `fig4_networks_<celltype>` (×3) — bipartite TF→target networks, top-5 targets per TF (by importance into that cell type's top-50 markers), edge width ∝ importance, TF nodes colored by family from one 13-family palette shared across fig2 and fig4.

---

## 4. Caveats and limitations to state explicitly in Methods

1. **Table S10 marker coverage:** 50 rows → 49 unique gene IDs (1 exact duplicate row) → **46 present** in this Seurat object's RNA assay. The 3 missing (`SbiRTX430.10G175100`, `SbiRTX430.08G154200`, `SbiRTX430.06G155500`) were confirmed absent from **both** the RNA (25,464 ft) and SCT (23,974 ft) assay rownames of the source RDS directly — not an artifact of the export pipeline. Likely explanation: annotation-version differences or an undocumented upstream filter between the published Table S10 and this processed object.

2. **GRNBoost2 importance is not an absolute or probabilistic score.** Read from arboreto 0.1.6 source (`arboreto/core.py`, `arboreto/algo.py`): `grnboost2()` fits one `GradientBoostingRegressor` per target gene (`SGBM_KWARGS`: learning_rate=0.01, n_estimators=5000 max, max_features=0.1, **subsample=0.9**), which enables an out-of-bag early-stopping heuristic (25-round trailing window). Per-edge importance = sklearn's `feature_importances_` (normalized to sum to 1.0 across all 1,022 TF predictors, **separately for each target's own model**) × the number of boosting rounds actually used before early stopping (empirically confirmed on this data: per-target totals range 25–63 rounds, median 27, across all mature-endodermis targets). **Consequence for Methods text:** "summed importance" values are valid for ranking TFs *within the same fixed target set* (as done in Step 05) but are not comparable in absolute magnitude across different marker sets or cell types, since each target's contribution is rescaled by its own convergence depth.

3. **GRNBoost2 exact-tie blocks.** 5.2–13.2% of edges in each cell type's full network (mature-endodermis 13.15%, dividing cells 5.24%, exodermis 11.83%) fall into "tie blocks" — 3+ edges from the same TF sharing a bit-for-bit identical importance value. Root cause, confirmed by direct inspection: these ties occur when multiple target genes have **byte-for-byte identical raw-count expression vectors** (e.g., one verified 19-way tie block where all 19 target genes were detected in exactly the same single nucleus, count=1, zero elsewhere) — identical inputs to independent GBM fits deterministically produce identical `feature_importances_`. This is a mechanical consequence of near-zero/duplicate-pattern expression, not a general GBM artifact, and it does not affect the top-10 TF rankings reported above: none of the 90 top-10-TF-into-marker edges (30 TFs × 3 cell types, plus the ranks-11–30 control of 60 more) derive any importance from tie-block edges, because the top-50 DE marker genes are specifically the most robustly/differentially expressed genes and structurally avoid the low-expression regime where ties concentrate.

4. **TF↔target self-referential edges in fig4.** In all three cell types, at least one of the top-5 target edges for a top-10 TF points to a gene that is *itself* also a top-10 TF for that same cell type (5 such edges total; mature-endodermis has a reciprocal pair). These are genuine GRNBoost2/co-expression-based inferences, not artifacts, and are retained in the figures (flagged with a distinct crimson-dashed style) as an explicit illustration of the method's inherent inability to distinguish direction/causality from co-expression: `SbiRTX430.09G017700` (rank 1 MYB) ↔ `SbiRTX430.09G017600` (rank 2 MYB) in mature-endodermis, raw-count Pearson r=0.562 (p=7.4e-25); full list with importance and correlation values in `figures/fig4_flagged_tf_tf_edges.csv`.

5. **SbiRTX430.09G017700 / SbiRTX430.09G017600 (mature-endodermis rank 1/2 MYBs).** Locus IDs differ by exactly 100 on chromosome 09, consistent with immediately adjacent genes under this genome annotation's ID-increment convention — **not confirmed against a GFF/coordinate file** (none available in the project's raw data), so state as an inference from ID adjacency, not a confirmed tandem-duplicate call if used in the manuscript. Co-expression across the 282 mature-endodermis nuclei: raw-count Pearson r=0.562 (p=7.4e-25), Spearman ρ=0.586 (p=2.4e-27); log-normalized Pearson r=0.487 (p=3.3e-18). Both genes detected in ~83% of nuclei.

6. **Deviations from the original pipeline plan** (both by explicit user decision, recorded in `prompt.md` and `logs/run_log.md`):
   - RDS byte-size correction: 1,288,490,189 (originally stated, an unverified GEO-display estimate) → 1,314,285,536 (actual, confirmed).
   - Step 04 scheduling: the two non-probe cell types were submitted concurrently with the probe rather than sized from its measured runtime; wall time was instead scaled from nuclei count (dividing cells 8 h, exodermis 4 h).

---

*Generated 2026-08-10. Source of truth for anything not captured here: `logs/run_log.md` (full chronological narrative) and `prompt.md` (pipeline specification, kept in sync with all deviations).*
