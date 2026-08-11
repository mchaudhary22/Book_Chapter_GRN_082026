# Chapter source bundle — sorghum root atlas GRN pipeline

Verbatim concatenation of source-of-truth files for book-chapter drafting, assembled 2026-08-10.
Every file below is reproduced in full, unmodified, inside a fenced code block.

---

## 1. methods_facts.md

```markdown
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
```

---

## 2. prompt.md (current version, as executed)

```markdown
# Single-cell GRN inference pipeline — sorghum root atlas (ARC / TinkerCliffs)

## Role

You are executing a reproducible bioinformatics pipeline on Virginia Tech's ARC TinkerCliffs cluster (Linux, SLURM). Write numbered scripts to `scripts/`, checkpoint every stage to disk, and hold no state between stages. Do not modify anything in `data/raw/`.

## Critical cluster rules

- **Never run the GRNBoost2 inference step on a login node.** It must be submitted to SLURM with `sbatch`. Login-node compute is against ARC policy and will be killed.
- Light steps (file inspection, input preparation, figure generation) may run interactively if they take under a few minutes; otherwise submit them too.
- Do not start long-running background processes in the terminal.

## Environment

Do not assume module names or partition names. Discover them (Step 00 below) and report before proceeding.

## Inputs (in `data/raw/`, read-only)

- `GSE297576_seurat_object.bicolor_root_atlas.RDS` — Seurat v5 object, 19,316 nuclei, 25,464 RNA features, gene IDs in `SbiRTX430.*` namespace, 1314285536 bytes (confirmed by `ls -la` in Step 00; the previously stated 1288490189 was an unverified estimate derived from GEO's rounded "1.2 Gb" display, not a real byte count)
- `Sbi_RTx430_TF_list.txt` — 1,827 rows, two tab-separated columns (protein ID like `SbiRTX430.01G002500.1.p`, TF family), no header
- `41477_2025_2047_MOESM3_ESM.xlsx` — sheet `TableS10.sorghum_lit_markers` has marker genes with both `SbicolorRTx430` and `SbicolorBTx623` ID columns

## Fixed parameters

- Cell-type column: `celltype` (27 labels). IGNORE `cluster_label.long` — leaf-atlas artifact with NAs.
- Keep only nuclei where the doublet column equals `Singlet`. **Expect exactly 14,792.**
- Assay: `RNA`, raw counts. Not SCT, not scale.data. Note the object's DefaultAssay is SCT, so bare `nrow(obj)` reports 23,974 — use the RNA assay explicitly and expect 25,464.
- TF detection threshold: present in the matrix AND detected in ≥50 singlet nuclei. **Expect 1,370 present, 1,022 after threshold.**
- Target cell types and expected singlet counts: `exodermis` (315), `dividing cells` (801), `mature-endodermis` (282)
- Random seed: 42 everywhere it can be set

## Known dependency issues — apply these, they are already diagnosed

1. **pySCENIC 0.12.1 requires `setuptools<81`.** Later versions remove `pkg_resources`, which `ctxcore` imports.
2. **arboreto 0.1.6 is incompatible with dask's query-planning backend.** `create_graph()` calls `from_delayed()` on an empty list, raising `TypeError: Must supply at least one delayed object`. Set `dask.config.set({"dataframe.query-planning": False})` **before importing arboreto or distributed.**
3. **Cap BLAS/OpenMP threads before importing numpy or sklearn**, and construct the dask cluster explicitly with one thread per worker. Accepting defaults caused total run failure through worker oversubscription on the previous machine. Set `OMP_NUM_THREADS`, `MKL_NUM_THREADS`, `OPENBLAS_NUM_THREADS` to 1.
4. There is no SSL certificate problem on Linux. Do not port any SSL patch.

---

## Step 00 — Environment reconnaissance (interactive, report before proceeding)

Do not write pipeline code yet. Report:

1. `module avail` output filtered for R, Python, Anaconda/Miniconda — with exact module names and versions
2. Whether an R module provides Seurat 5, or whether Seurat must be installed to a personal library
3. Whether SeuratDisk and hdf5r are available in that R
4. Available SLURM partitions, their time limits, and per-node core and memory limits (`sinfo`, and your account's allocation)
5. Available disk quota in `/projects/intro2gds/Manisha_Book_chapter_082026` and in scratch
6. Whether outbound internet works from a compute node (relevant for conda installs)

**Stop and report. Do not proceed until I confirm.**

## Step 01 — Build the environment

Create a conda environment named `grn` with Python 3.10 and: pyscenic, scanpy, anndata, pandas, numpy, matplotlib, seaborn, openpyxl, and `setuptools<81`. Export to `environment.yml` at the project root with pinned versions and report every resolved version.

Verify with a synthetic smoke test before going further: run GRNBoost2 on a 50-cell, 20-gene random matrix with 5 TFs and confirm it returns 95 edges. Do this on a compute node, not a login node.

**Report and stop.**

## Step 02 — Export Seurat to h5ad (`scripts/02_export_seurat.R`)

Load the RDS, subset to singlets, extract RNA assay counts and full cell metadata, write `data/processed/root_singlets.h5ad`.

Use SeuratDisk if available. If not, write MatrixMarket counts plus metadata and gene names from R, then assemble the h5ad in Python (`scripts/02b_assemble_h5ad.py`). Document which route you used.

**Resolved in Step 00:** SeuratDisk is not available in any ARC module and will not be installed to a personal library (user decision, 2026-08-10). Use the MatrixMarket + Python assembly route.

Report dimensions before and after filtering, and per-cell-type counts for all 27 types. Verify by reading the h5ad back and confirming the matrix orientation is cells × genes.

## Step 03 — Prepare regulators (`scripts/03_prepare_inputs.py`)

Read the h5ad. Read the TF list, strip the trailing `.1.p`, keep the family column. Filter to TFs present in the matrix and detected in ≥50 nuclei.

Write `data/processed/regulators.txt` (one gene ID per line) and `data/processed/tf_families.csv`. Report TF count at each stage and the family breakdown.

**Stop here and let me verify before the expensive step.**

## Step 04 — GRNBoost2 (`scripts/04_run_grnboost2.py` + `scripts/04_submit.sh`)

Write both the Python script and a SLURM batch script. The batch script must request an appropriate partition, wall time, core count, and memory based on what you found in Step 00, and must be submitted with `sbatch`.

**Resolved in Step 00 (user decision, 2026-08-10):** target `normal_q` with QOS `tc_normal_short`. Probe job (`mature-endodermis`) requests 1 node, 64 cores, 64 GB memory, 4 hour wall time; GRNBoost2 dask cluster uses `n_workers=64`, one thread per worker.

**Revised 2026-08-10 (user decision):** rather than waiting for the probe's measured runtime, the remaining two cell types were submitted concurrently with the probe, same resources (1 node, 64 cores, 64 GB, `n_workers=64`, `threads_per_worker=1`, seed 42), wall time scaled from nuclei count instead of measurement: `dividing cells` (801 nuclei, ~3x the probe's 282) got 8 hours; `exodermis` (315 nuclei, close to the probe) kept 4 hours.

For each of the three target cell types: subset the matrix to that cell type's nuclei, run GRNBoost2 with the regulator list as TFs and all genes as candidate targets, seed fixed. Write `results/grnboost2_edges_<celltype>.csv` with columns TF, target, importance.

Support running a single cell type via a command-line argument so a timing probe can be submitted first. Submit `mature-endodermis` alone as the probe, report its runtime and edge count (and peak memory from `sacct`), then size the remaining two.

If one cell type fails, continue with the others and report which failed.

## Step 05 — Select TFs (`scripts/05_select_tfs.py`)

Per cell type: identify marker genes by Wilcoxon differential expression against all other cell types, top 50 by adjusted p-value. For each TF, sum importance across edges into those markers. Rank, take top 10.

Write `results/selected_tfs_<celltype>.csv` with TF, summed importance, rank, family. Report the top 10 per cell type with families.

## Step 06 — Figures (`scripts/06_figures.py`)

Four figures in `figures/`, 300 dpi, PNG and PDF:

1. `fig1_umap_celltypes` — UMAP coloured by `celltype`, all 27 types, readable legend
2. `fig2_tf_featureplots` — UMAP feature plots of the selected TFs, one panel per TF, grouped by cell type
3. `fig3_marker_dotplot` — dotplot of Table S10 root markers (`SbicolorRTx430` column) across all 27 cell types
4. `fig4_networks_<celltype>` — bipartite plot per target cell type: 10 TFs and their top targets, edge width by importance, TF nodes coloured by family

Publication quality: legible fonts, no clipped labels, colourblind-safe palettes.

## Step 07 — Run log

Maintain `logs/run_log.md` throughout: date, module and package versions, SLURM job IDs, resources requested, every parameter, per-step runtimes, checkpoint dimensions, and every deviation from this file with its reason.

---

## Rules

- Stop and ask before deviating from any fixed parameter.
- If an observed count differs from an expected value stated above, stop and report rather than proceeding.
- Never write to `data/raw/`.
- Fail loudly rather than substituting defaults.
- Never run compute on a login node.
- If you change how the pipeline works, update this file so it stays an accurate specification.
```

---

## 3. logs/run_log.md

```markdown
# Run log — sorghum root atlas GRN pipeline

## 2026-08-10 — Step 00: Environment reconnaissance

- Ran on compute node `tc-xe005` (SLURM job 6961514, partition `h200_normal_q`), not a login node.
- Modules discovered:
  - R: `R/4.3.2-gfbf-2023a`, `R/4.4.2-gfbf-2024a`, `R/4.5.2-gfbf-2025b` (default)
  - R bundles: `R-bundle-CRAN/{2023.12-foss-2023a, 2024.11-foss-2024a, 2025.11-foss-2025b(D)}`, `R-bundle-Bioconductor/3.20-foss-2024a-R-4.4.2`
  - Python: `Python/{3.11.3, 3.11.5, 3.12.3, 3.13.1, 3.13.5, 3.14.2(D)}`
  - Conda: `Miniconda3/24.7.1-0`, `Miniconda3/25.11.1-1` (default); `mamba` 2.4.0 available via bundled Miniforge3/25.11.0-1.
- `R/4.4.2-gfbf-2024a` + `R-bundle-Bioconductor/3.20-foss-2024a-R-4.4.2` provides Seurat 5.3.1, SeuratObject 5.3.0, hdf5r 1.3.11 pre-installed. **SeuratDisk not available** in any module.
- SLURM: `normal_q`/`preemptable_q` (CPU, 312 nodes, 96-128 cores, ~252-257 GB/node) is the right target for GRNBoost2 (no GPU needed). Account `introtogds` has QOS `tc_normal_short` (1 day), `tc_normal_base` (7 days), `tc_normal_long` (14 days), etc. ~1.999M core-hours remaining on TinkerCliffs.
- Disk: `/projects/intro2gds` — 8.33 TB / 25 TB used, 5.44M / 10.49M files. Home — 259 GB / 640 GB.
- Outbound internet confirmed working from compute node (pypi.org, repo.anaconda.com, conda.anaconda.org all reachable).
- **Deviation noted and corrected:** `prompt.md` originally stated the RDS file as 1,288,490,189 bytes; actual size is 1,314,285,536 bytes. Confirmed this was an unverified estimate from GEO's rounded "1.2 Gb" display, not a real mismatch. `prompt.md` updated with the correct byte count.

**User decisions recorded (2026-08-10):**
1. Do not install SeuratDisk — use MatrixMarket + Python assembly route for Step 02.
2. Conda env built in default home location (`~/.conda/envs`).
3. Step 04 GRNBoost2 target: `normal_q` partition, QOS `tc_normal_short`.
4. Step 04 probe job (`mature-endodermis`): 1 node, 64 cores, 64 GB memory, 4 hour wall time, dask `n_workers=64`, one thread per worker. Remaining two cell types to be sized from measured probe runtime/peak memory (via `sacct`), not estimated in advance.

## 2026-08-10 — Step 01: Build the environment

- Created conda environment `grn` (Python 3.10) via `mamba create -p ~/.conda/envs/grn -c conda-forge -c bioconda ...`.
  - First attempt failed: `mamba` defaulted to the shared system Miniforge path (`/apps/common/software/Miniforge3/.../envs/`), not writable by this user. Retried with explicit `-p ~/.conda/envs/grn`. Succeeded — 278 packages, no conflicts.
- Resolved versions (see `environment.yml` at project root for full pinned list):
  - python 3.10.0, pyscenic 0.12.1, scanpy 1.10.4, anndata 0.11.4, pandas 2.3.3, numpy 1.23.5, matplotlib 3.10.9, seaborn 0.13.2, openpyxl 3.1.5, setuptools 80.10.2 (satisfies `<81` pin)
  - arboreto 0.1.6, ctxcore 0.2.0, dask 2024.8.2, distributed 2024.8.2 (pulled in transitively by pyscenic)
- Exported `environment.yml` to project root (`conda env export --no-builds`, stripped user-specific `prefix:` line for portability).
- **Synthetic GRNBoost2 smoke test** (`scripts/01_smoke_test_grnboost2.py`): 50 cells × 20 genes, 5 TFs, seed 42, run on compute node `tc-xe005` (not login node).
  - First run failed: `RuntimeError` from Python multiprocessing spawn — script had no `if __name__ == "__main__":` guard, so dask worker processes re-imported and re-executed the module-level `LocalCluster()` call recursively. Fixed by wrapping logic in `main()` guarded by `__main__`. Killed the stuck leftover process tree from the first (broken) run before retrying.
  - Second run: **succeeded — 95 edges**, matching the expected count exactly. Applied all three documented dependency fixes (dask query-planning disabled before importing arboreto/distributed; OMP/MKL/OPENBLAS thread caps set to 1 before numpy import; explicit `LocalCluster(n_workers=4, threads_per_worker=1, processes=True)`).

## 2026-08-10 — Step 02: Export Seurat to h5ad

- Ran on compute node `tc-xe005` (not login node), modules `R/4.4.2-gfbf-2024a` + `R-bundle-Bioconductor/3.20-foss-2024a-R-4.4.2`.
- `scripts/02_export_seurat.R`: RDS load took 11.1 s. Object has assays `RNA`, `SCT` (DefaultAssay as loaded was `SCT`, switched to `RNA` explicitly per fixed parameters). Reductions present: `pca`, `umap`.
- Doublet column auto-detected (not assumed) by scanning metadata for a column containing the value `"Singlet"`: exactly one match, column **`doublet`** (values: Doublet=4524, Singlet=14792, NA=0).
- **Before filtering:** 19,316 cells x 25,464 RNA features (matches expected 25,464).
- **After filtering (Singlet only):** 14,792 cells (matches expected exactly) x 25,464 RNA features.
- `celltype` confirmed 27 levels. All 27 per-cell-type counts recorded in `data/processed/celltype_counts.csv`. Target cell-type checks all passed exactly: exodermis=315, dividing cells=801, mature-endodermis=282.
- RNA layer was a single `counts` layer (no multi-sample split layers needing `JoinLayers`).
- Route used: **MatrixMarket + Python assembly** (SeuratDisk unavailable, per Step 00 / user decision). R wrote `counts.mtx` (25464 x 14792, dgCMatrix), `genes.csv`, `cells.csv`, `metadata.csv`, `umap_embedding.csv` (14792 x 2) to `data/processed/`.
- `scripts/02b_assemble_h5ad.py` assembled `data/processed/root_singlets.h5ad`, transposing to cells x genes.
- **Verified by reading the h5ad back:** shape (14792, 25464), orientation confirmed cells x genes. All 27 celltype counts in the h5ad match the R-side counts exactly. `X` is raw integer counts (dtype int64, max value 525 — not normalized/scaled). `X_umap` present in `obsm`, shape (14792, 2). `var_names` confirmed in `SbiRTX430.*` namespace.
- No deviations from fixed parameters; all expected counts matched on first run — no stop-and-report needed.

## 2026-08-10 — Step 03: Prepare regulators

- `scripts/03_prepare_inputs.py`, run in `grn` conda env.
- TF list: 1,827 rows (matches expected). Stripped trailing `.1.p` from protein IDs to get gene IDs — 0 duplicate gene IDs after stripping (1:1 protein-to-gene mapping in this list).
- Present in RNA matrix: **1,370** (matches expected).
- Detected in >=50 nuclei (nonzero raw count): **1,022** (matches expected).
- 54 distinct TF families represented after thresholding; largest are bHLH (96), MYB (73), bZIP (72), NAC (60), ERF (57).
- Wrote `data/processed/regulators.txt` (1,022 gene IDs) and `data/processed/tf_families.csv` (gene_id, family, n_nuclei_detected).
- No deviations; all three expected counts (1827 -> 1370 -> 1022) matched exactly.

## 2026-08-10 — Step 04: GRNBoost2

- `scripts/04_run_grnboost2.py` + `scripts/04_submit.sh` written. Submission target: `normal_q` partition, QOS `tc_normal_short`, account `introtogds` (per `sacctmgr show user` default account).
- **Probe job (`mature-endodermis`)**: SLURM job **6963109**, submitted via `sbatch scripts/04_submit.sh mature-endodermis 64` (1 node, 64 cores, 64 GB, 4h limit, `n_workers=64`, `threads_per_worker=1`, seed 42). Ran on node `tc077`.
  - **COMPLETED. Elapsed 00:06:56 (script-reported runtime 382.0 s / 6.37 min). Peak memory (MaxRSS) 16,225,436 KB ≈ 15.5 GB** (well under the 64 GB request).
  - 282 cells, 1022/1022 regulators present in data. **1,246,780 edges** written to `results/grnboost2_edges_mature-endodermis.csv`. Importance range 5.5e-18 to 14.9 (mean 0.40, median 0.20).
- **User decision (2026-08-10):** rather than waiting on the probe's measured runtime, submitted the remaining two cell types concurrently, wall time scaled from nuclei count instead of measured: `dividing cells` (801 nuclei, ~3x probe) got 8h; `exodermis` (315 nuclei, close to probe) kept 4h. `prompt.md` updated to record this deviation from the originally planned probe-then-size sequencing.
  - `dividing cells`: SLURM job **6963283** (`sbatch --time=08:00:00 --job-name=grnboost2_dividing_cells scripts/04_submit.sh "dividing cells" 64`)
  - `exodermis`: SLURM job **6963284** (`sbatch --time=04:00:00 --job-name=grnboost2_exodermis scripts/04_submit.sh "exodermis" 64`)
  - **Both COMPLETED.** `dividing cells`: elapsed 00:05:29 (script-reported 290.6 s), MaxRSS 16,781,100 KB ≈ 16.0 GB, 801 cells, 1022/1022 regulators present, **1,563,314 edges** to `results/grnboost2_edges_dividing_cells.csv` (importance range 1.4e-18–15.2, mean 0.366, median 0.200).
  - `exodermis`: elapsed 00:06:24 (script-reported 355.0 s), MaxRSS 16,109,460 KB ≈ 15.4 GB, 315 cells, 1022/1022 regulators present, **1,380,256 edges** to `results/grnboost2_edges_exodermis.csv` (importance range 1.2e-18–14.2, mean 0.387, median 0.201).
  - Note: `dividing cells` (~3x the probe's nuclei count) ran *faster* than the probe (290.6 s vs 382.0 s) — GRNBoost2 runtime here is dominated by regulator/target-set size and dask scheduling overhead more than by cell count at this scale, not a cause for concern.
  - All three cell types now complete. Peak memory across all three jobs stayed at ~15–16 GB, far under the 64 GB requested.

## 2026-08-10 — Step 05: Select TFs

- `scripts/05_select_tfs.py`, run interactively on compute node `tc-xe005` (light step, ~a few minutes).
- Marker genes: Wilcoxon DE (`sc.tl.rank_genes_groups`, `celltype` vs rest, all 27 groups computed once) on a normalized+log1p **copy** of the AnnData (raw counts on disk untouched) — normalization (`normalize_total` target_sum=1e4, `log1p`) applied only for the DE test itself, standard practice for rank-based DE; not a deviation from the "RNA raw counts, not SCT/scale.data" fixed parameter, which governs the assay/slot source, not DE preprocessing. Top 50 markers per cell type taken by ascending adjusted p-value.
- For each TF, summed GRNBoost2 importance across edges into those 50 markers; ranked; top 10 written to `results/selected_tfs_<celltype>.csv` (TF, summed_importance, rank, family).
- **Tie-block characterization** (user-requested, 2026-08-10): exact-tie blocks defined as >=3 edges from the same TF sharing an identical importance value, computed over each cell type's full edge set.
  - mature-endodermis: 163,986 / 1,246,780 edges (13.15%) in tie blocks.
  - dividing cells: 81,963 / 1,563,314 edges (5.24%) in tie blocks.
  - exodermis: 163,224 / 1,380,256 edges (11.83%) in tie blocks.
  - For every cell type, **none** of the top-10 selected TFs derive any of their summed importance from tie-block edges (frac_from_tie_blocks = 0.0000 for all 30 TF entries across the three cell types) — the top-10 rankings are not tie-block artifacts, even though 5-13% of each network's edges genome-wide are exact ties. Full breakdown in `results/tie_block_summary.csv`.
- No stop-and-report triggers — pipeline has no pre-stated expected values for Step 05 outputs.

## 2026-08-10 — Step 05 rank 2-10 comparison (user-requested)

- mature-endodermis: 74.81, 25.02, 21.80, 20.28, 19.72, 18.89, 17.16, 15.06, 13.26, 12.24 (rank 1-10)
- dividing cells: 73.84, 54.56, 39.35, 35.32, 20.26, 19.06, 17.96, 17.63, 15.15, 11.42
- exodermis: 74.51, 27.78, 26.64, 25.09, 24.56, 23.62, 19.41, 17.26, 14.40, 14.30
- Rank-1/rank-2 ratio: mature-endodermis 2.99x, exodermis 2.68x, dividing cells only 1.35x. The tight rank-1 agreement (~74-75) across cell types is specific to the top rank, not a general distributional property — dividing cells has a much more gradual decay (strong secondary/tertiary regulators) while the other two drop sharply after rank 1 then converge with dividing cells by rank 8-10.

## 2026-08-10 — Step 06: Figures

- `scripts/06_figures.py`, run interactively on compute node `tc-xe005`. All 4 figures written to `figures/` as PNG + PDF at 300 dpi.
- **Implementation choices** (spec left these underspecified, documented per "update this file" rule):
  - fig2: one combined figure, 3 rows (celltype) x 10 columns (TF rank 1-10) rather than 30 separate files, per "grouped by cell type."
  - fig4: top 5 targets per TF (by importance into that celltype's own top-50 markers), for legibility with 10 TFs.
  - TF-family colors (13 families across all three top-10 lists) assigned from one global palette shared across fig2 and fig4, so the same family reads as the same color everywhere.
- **Two bugs hit and fixed during this step:**
  1. fig3 `sc.pl.dotplot(..., gene_symbols="marker_alias")` initially passed `var_names` in the gene-ID namespace — when `gene_symbols` is set, scanpy expects `var_names` values to be in the *symbol* namespace instead (it maps symbol -> real var_name internally). Fixed by passing the alias/display-label list as `var_names`. Also found a real duplicate: alias "SWEET13" is reused across two different gene IDs (`SbiRTX430.08G105800`, `SbiRTX430.08G105600`) in Table S10 — disambiguated by appending the gene-ID suffix to the display label for any reused alias, verified unique before plotting.
  2. fig3 initially saved as a blank plot with only the title visible — `sc.pl.dotplot(..., return_fig=True)` returns the `DotPlot` object itself, not a rendered matplotlib figure; `plt.gcf()` grabbed an empty canvas. Fixed by calling `dp.make_figure()` explicitly before accessing `dp.fig`.
  - Both fixes verified by visually inspecting the rendered PNGs (Read tool) after rerun; all 4 figure types confirmed legible with correct content (fig1: 27 celltypes distinguishable with full legend; fig2: 30 panels, no clipped labels; fig3: 46/50 Table S10 markers with sensible cell-type-specific patterns, e.g. xylem markers IRX1/XCP1/CESA7 on xylem cell types; fig4: TF-TF edges render in the expected crimson-dashed style).
  - Table S10 markers: 50 rows -> 49 unique gene IDs (1 exact duplicate row) -> 46 present in the expression matrix. 3 not present: `SbiRTX430.10G175100`, `SbiRTX430.08G154200`, `SbiRTX430.06G155500` (excluded from fig3).
- **User-requested TF-TF edge flag** (generic detection: any fig4 edge whose target is itself one of that celltype's own top-10 TFs), written to `figures/fig4_flagged_tf_tf_edges.csv`. Caught the requested case plus 3 more, not previously identified:
  - mature-endodermis: `SbiRTX430.09G017700 -> SbiRTX430.09G017600` (importance 5.411) **and the reverse** `SbiRTX430.09G017600 -> SbiRTX430.09G017700` (importance 4.585) — both directions present since each is in the other's top-5 marker targets. Raw-count Pearson r=0.562 (matches Step 05 verification).
  - dividing cells: `SbiRTX430.07G054700 -> SbiRTX430.09G224300` (importance 1.731, r=0.257)
  - exodermis: `SbiRTX430.03G310600 -> SbiRTX430.03G038000` (importance 2.131, r=0.373) and `SbiRTX430.09G219200 -> SbiRTX430.03G038000` (importance 1.564, r=0.339)
  - All flagged edges drawn in crimson dashed style with target nodes marked with a star and "[also top-10 TF]" label in fig4; not removed, per instruction.

## 2026-08-10 — Step 06 verification: missing Table S10 markers (user-requested)

- Confirmed directly against the original RDS (not the h5ad export): `SbiRTX430.10G175100`, `SbiRTX430.08G154200`, `SbiRTX430.06G155500` are absent from `rownames(obj[["RNA"]])` (25,464 features) **and** absent from `rownames(obj[["SCT"]])` (23,974 features) in the source Seurat object. No partial/near-string matches found on the locus-ID substring either (ruling out a formatting mismatch). These three markers are genuinely not present anywhere in this Seurat object's feature space — not an artifact of the Step 02 MatrixMarket export or the Step 06 presence filter. Likely explanation: Table S10 (from the published paper) references genes filtered out of this particular processed object, possibly due to annotation-version differences or an upstream low-expression filter not documented in this RDS.

## 2026-08-10 — Step 05 verification (user-requested audit, `scripts/05b_verify_importance_and_ties.py`)

- **Tie-block metric verified real, not a join/key bug.** Control check: ranks 11-30 (60 additional TF entries across 3 cell types) all also show frac_from_tie_blocks=0.0. Manual audit of rank-1 mature-endodermis TF `SbiRTX430.09G017700`: listed all 45 individual edges into the 50 markers — zero duplicate importance values among them (each is a unique float). But the SAME TF has real tie blocks elsewhere: 1,235/7,530 edges (16.4%) in its full network are tied. Example: a 19-way tie at importance=3.156694826102661 — those 19 tied target genes are detected in only 0.4% of nuclei (mean count 0.0035), i.e. essentially undetected/dropout-dominated. By contrast this TF's actual (non-tied) marker edge target is detected in 82.6% of nuclei. **Mechanism identified:** GBM ties arise for near-zero-expression targets where the booster can't discriminate; the top-50 DE markers are specifically the most robustly/differentially expressed genes and structurally avoid this regime, so ties naturally never land on marker edges for well-ranked TFs — the 0.0 values are a real consequence of algorithm + marker-selection interaction, not an artifact of the analysis code.
- **arboreto 0.1.6 importance computation (read from installed `arboreto/core.py`, `arboreto/algo.py`):** `grnboost2()` calls `diy(regressor_type='GBM', regressor_kwargs=SGBM_KWARGS)` where `SGBM_KWARGS = {learning_rate:0.01, n_estimators:5000, max_features:0.1, subsample:0.9}`. Since `subsample<1.0`, `is_oob_heuristic_supported`→True, so `GradientBoostingRegressor.fit()` runs with `EarlyStopMonitor` (OOB-improvement based, 25-round trailing window) and typically halts well before 5000 rounds. `to_feature_importances()` then computes `denormalized_importances = trained_regressor.feature_importances_ * n_estimators` where `n_estimators = len(trained_regressor.estimators_)` (the actual number of boosting rounds used, not the configured max). sklearn's `feature_importances_` is normalized to sum to 1.0 across all ~1,022 TF predictors *for that one target's model*. **Empirically confirmed** on this dataset: summing importance over all TFs for a single target gives near-integer values (e.g. 25.0, 26.0, 27.0, 29.0, 32.0 for the first 15 markers) matching the number of boosting rounds used; range across all 17,806 mature-endodermis targets is 25 (the minimum possible, since `EARLY_STOP_WINDOW_LENGTH=25`) to 63, median 27 — most targets converge in the minimum window.
  - **Methods-text implication:** per-edge importance = (that TF's fractional share of variance-reduction, normalized to sum to 1 across all candidate regulators, in a GBM fit independently for that one target) × (number of boosting rounds used for that target before early stopping, 25-63 in this data, median 27). "Summed importance" for a TF across its top-50-marker edges therefore aggregates fractional shares from up to 50 independently-fit models, each modestly rescaled (~2.5x range) by that target's own convergence depth. It is **not an absolute or probabilistic regulatory-strength score** — TFs are validly comparable to each other *within the same marker set* (same models, same targets), but summed-importance magnitudes should not be compared across different marker sets or cell types without accounting for this per-target scaling.
- **Tie mechanism confirmed exactly** (`scripts/05c_tie_mechanism.py`): the 19-way tie at importance=3.156694826102661 for `SbiRTX430.09G017700` is caused by **byte-for-byte identical target expression vectors**, not general GBM behavior. All 19 tied targets are detected in exactly one nucleus each — and it's the *same single nucleus* (index 145 of 282) for all 19, each with raw count=1 and 0 elsewhere. 19/19 columns are exactly identical (`np.all(X == X[:, [0]])`); the tie-group submatrix has only 19 nonzero entries total across 282×19 cells, all count=1, all in the same row. Since GRNBoost2 fits one independent GBM model per target using the same TF matrix as predictors, **identical y vectors deterministically produce identical fitted models and therefore identical `feature_importances_`** — this is not a coincidental tie or a general property of gradient boosting on similar-but-different data, it is a direct, mechanical consequence of duplicate target columns in the expression matrix (most likely genes that are silent everywhere except one outlier/doublet-like or high-depth nucleus, each detected only there at the minimum measurable count of 1).
- **SbiRTX430.09G017700 vs SbiRTX430.09G017600** (mature-endodermis rank-1/rank-2 MYBs): locus IDs differ by exactly 100 on chromosome 09, consistent with immediately adjacent genes under this annotation's ID-increment convention — **not confirmed via GFF** (no annotation file present in `data/raw`), so this is an inference from ID adjacency only, not a coordinate-verified tandem-duplicate call. Expression correlation across the 282 mature-endodermis nuclei: raw-count Pearson r=0.562 (p=7.4e-25), Spearman rho=0.586 (p=2.4e-27); log-normalized Pearson r=0.487 (p=3.3e-18). Both genes broadly detected (~83% of nuclei) — moderate-strong positive co-expression, consistent with (but not proof of) tandem-duplicate shared regulation.
```

---

## 4. environment.yml

```yaml
name: grn
channels:
  - conda-forge
  - bioconda
dependencies:
  - _openmp_mutex=4.5
  - alsa-lib=1.2.16.1
  - anndata=0.11.4
  - arboreto=0.1.6
  - array-api-compat=1.15.0
  - attrs=26.1.0
  - aws-c-auth=0.9.6
  - aws-c-cal=0.9.13
  - aws-c-common=0.12.6
  - aws-c-compression=0.3.2
  - aws-c-event-stream=0.5.9
  - aws-c-http=0.10.10
  - aws-c-io=0.26.1
  - aws-c-mqtt=0.14.0
  - aws-c-s3=0.11.5
  - aws-c-sdkutils=0.2.4
  - aws-checksums=0.2.10
  - aws-crt-cpp=0.37.3
  - aws-sdk-cpp=1.11.747
  - azure-core-cpp=1.16.2
  - azure-identity-cpp=1.13.3
  - azure-storage-blobs-cpp=12.16.0
  - azure-storage-common-cpp=12.13.0
  - azure-storage-files-datalake-cpp=12.14.0
  - backports.zstd=1.6.0
  - bokeh=3.9.2
  - boltons=26.1.0
  - brotli=1.2.0
  - brotli-bin=1.2.0
  - brotli-python=1.2.0
  - bzip2=1.0.8
  - c-ares=1.34.8
  - ca-certificates=2026.7.22
  - cached-property=1.5.2
  - cached_property=1.5.2
  - cairo=1.18.4
  - click=8.4.2
  - cloudpickle=3.1.2
  - contourpy=1.3.2
  - ctxcore=0.2.0
  - cycler=0.12.1
  - cyrus-sasl=2.1.28
  - cytoolz=1.1.0
  - dask=2024.8.2
  - dask-core=2024.8.2
  - dask-expr=1.1.13
  - dbus=1.16.2
  - dill=0.4.1
  - diptest=0.10.0
  - distributed=2024.8.2
  - double-conversion=3.3.1
  - et_xmlfile=2.0.0
  - exceptiongroup=1.3.1
  - font-ttf-dejavu-sans-mono=2.37
  - font-ttf-inconsolata=3.000
  - font-ttf-source-code-pro=2.038
  - font-ttf-ubuntu=0.83
  - fontconfig=2.18.3
  - fonts-conda-ecosystem=1
  - fonts-conda-forge=1
  - fonttools=4.63.0
  - freetype=2.14.3
  - frozendict=2.4.7
  - fsspec=2026.7.0
  - gflags=2.3.1
  - glog=0.7.1
  - graphite2=1.3.15
  - h2=4.4.1
  - h5py=3.16.0
  - harfbuzz=12.1.0
  - hdf5=1.14.6
  - hpack=4.2.0
  - hyperframe=6.1.0
  - icu=75.1
  - importlib-metadata=9.0.0
  - importlib_metadata=9.0.0
  - interlap=0.2.7
  - jinja2=3.1.6
  - joblib=1.5.3
  - keyutils=1.6.3
  - kiwisolver=1.5.0
  - krb5=1.21.3
  - lcms2=2.19.1
  - ld_impl_linux-64=2.46.1
  - legacy-api-wrap=1.5
  - lerc=4.2.0
  - libabseil=20260107.1
  - libaec=1.1.5
  - libarrow=23.0.1
  - libarrow-acero=23.0.1
  - libarrow-compute=23.0.1
  - libarrow-dataset=23.0.1
  - libarrow-substrait=23.0.1
  - libblas=3.11.0
  - libbrotlicommon=1.2.0
  - libbrotlidec=1.2.0
  - libbrotlienc=1.2.0
  - libcblas=3.11.0
  - libclang-cpp21.1=21.1.8
  - libclang-cpp22.1=22.1.8
  - libclang13=22.1.8
  - libcrc32c=1.1.2
  - libcups=2.3.3
  - libcurl=8.18.0
  - libdeflate=1.25
  - libdrm=2.4.127
  - libedit=3.1.20250104
  - libegl=1.7.0
  - libev=4.33
  - libevent=2.1.12
  - libexpat=2.8.1
  - libffi=3.4.6
  - libfreetype=2.14.3
  - libfreetype6=2.14.3
  - libgcc=16.1.0
  - libgcc-ng=16.1.0
  - libgfortran=16.1.0
  - libgfortran5=16.1.0
  - libgl=1.7.0
  - libglib=2.86.0
  - libglvnd=1.7.0
  - libglx=1.7.0
  - libgomp=16.1.0
  - libgoogle-cloud=2.39.0
  - libgoogle-cloud-storage=2.39.0
  - libgrpc=1.78.1
  - libhwloc=2.13.0
  - libiconv=1.18
  - libjpeg-turbo=3.2.0
  - liblapack=3.11.0
  - libllvm21=21.1.8
  - libllvm22=22.1.8
  - liblzma=5.8.3
  - liblzma-devel=5.8.3
  - libnghttp2=1.68.1
  - libnsl=2.0.1
  - libntlm=1.8
  - libopenblas=0.3.34
  - libopengl=1.7.0
  - libopentelemetry-cpp=1.21.0
  - libopentelemetry-cpp-headers=1.21.0
  - libparquet=23.0.1
  - libpciaccess=0.19
  - libpng=1.6.58
  - libpq=18.1
  - libprotobuf=6.33.5
  - libre2-11=2025.11.05
  - libsqlite=3.53.4
  - libssh2=1.11.1
  - libstdcxx=16.1.0
  - libstdcxx-ng=16.1.0
  - libthrift=0.22.0
  - libtiff=4.7.2
  - libutf8proc=2.11.3
  - libuuid=2.42.2
  - libvulkan-loader=1.4.357.0
  - libwebp-base=1.6.0
  - libxcb=1.17.0
  - libxcrypt=4.4.38
  - libxkbcommon=1.13.2
  - libxml2=2.15.1
  - libxml2-16=2.15.1
  - libxslt=1.1.43
  - libzlib=1.3.2
  - llvmlite=0.48.0
  - locket=1.0.0
  - loompy=3.0.8
  - lz4=4.4.5
  - lz4-c=1.10.0
  - markupsafe=3.0.3
  - matplotlib=3.10.9
  - matplotlib-base=3.10.9
  - msgpack-python=1.2.1
  - multiprocessing_on_dill=3.5.0a4
  - munkres=1.1.4
  - narwhals=2.24.0
  - natsort=8.4.0
  - ncurses=6.6
  - networkx=3.4
  - nlohmann_json=3.12.0
  - nomkl=1.0
  - numba=0.66.0
  - numexpr=2.14.1
  - numpy=1.23.5
  - numpy_groupies=0.11.3
  - openjpeg=2.5.4
  - openldap=2.6.10
  - openpyxl=3.1.5
  - openssl=3.6.3
  - orc=2.3.0
  - packaging=26.3
  - pandas=2.3.3
  - partd=1.4.2
  - patsy=1.0.2
  - pcre2=10.46
  - pillow=12.3.0
  - pip=26.2.1
  - pixman=0.46.4
  - prometheus-cpp=1.3.0
  - psutil=7.2.2
  - pthread-stubs=0.4
  - pyarrow=23.0.1
  - pyarrow-core=23.0.1
  - pynndescent=0.6.0
  - pyparsing=3.3.2
  - pyscenic=0.12.1
  - pyside6=6.9.3
  - pysocks=1.7.1
  - python=3.10.0
  - python-dateutil=2.9.0.post0
  - python-tzdata=2026.3
  - python_abi=3.10
  - pytz=2026.3.post1
  - pyyaml=6.0.3
  - qhull=2020.2
  - qt6-main=6.9.3
  - re2=2025.11.05
  - readline=8.3
  - s2n=1.7.1
  - scanpy=1.10.4
  - scikit-learn=1.7.2
  - scipy=1.15.2
  - seaborn=0.13.2
  - seaborn-base=0.13.2
  - session-info=1.0.0
  - setuptools=80.10.2
  - six=1.17.0
  - snappy=1.2.2
  - sortedcontainers=2.4.0
  - sqlite=3.53.4
  - statsmodels=0.14.6
  - stdlib-list=0.12.0
  - tbb=2023.0.0
  - tblib=3.2.2
  - threadpoolctl=3.6.0
  - tk=8.6.13
  - toolz=1.1.0
  - tornado=6.5.8
  - tqdm=4.70.0
  - typing_extensions=4.16.0
  - tzdata=2026c
  - umap-learn=0.5.12
  - unicodedata2=17.0.1
  - urllib3=2.7.0
  - wayland=1.24.0
  - wheel=0.47.0
  - xcb-util=0.4.1
  - xcb-util-cursor=0.1.6
  - xcb-util-image=0.4.0
  - xcb-util-keysyms=0.4.1
  - xcb-util-renderutil=0.3.10
  - xcb-util-wm=0.4.2
  - xkeyboard-config=2.48
  - xorg-libice=1.1.2
  - xorg-libsm=1.2.6
  - xorg-libx11=1.8.13
  - xorg-libxau=1.0.12
  - xorg-libxcomposite=0.4.7
  - xorg-libxcursor=1.2.3
  - xorg-libxdamage=1.1.6
  - xorg-libxdmcp=1.1.5
  - xorg-libxext=1.3.7
  - xorg-libxfixes=6.0.2
  - xorg-libxi=1.8.3
  - xorg-libxrandr=1.5.5
  - xorg-libxrender=0.9.12
  - xorg-libxtst=1.2.5
  - xorg-libxxf86vm=1.1.7
  - xyzservices=2026.3.0
  - xz=5.8.3
  - xz-gpl-tools=5.8.3
  - xz-tools=5.8.3
  - yaml=0.2.5
  - zict=3.0.0
  - zipp=4.1.0
  - zlib=1.3.2
  - zlib-ng=2.3.3
  - zstd=1.5.7
```

---

## 5. scripts/04_submit.sh

```bash
#!/bin/bash
#SBATCH --job-name=grnboost2
#SBATCH --account=introtogds
#SBATCH --partition=normal_q
#SBATCH --qos=tc_normal_short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=logs/slurm_%x_%j.out
#SBATCH --error=logs/slurm_%x_%j.err

# Usage: sbatch scripts/04_submit.sh <celltype> [n_workers]
# Example (probe): sbatch scripts/04_submit.sh mature-endodermis 64

set -euo pipefail

cd "$SLURM_SUBMIT_DIR"

CELLTYPE="${1:?Usage: sbatch 04_submit.sh <celltype> [n_workers]}"
N_WORKERS="${2:-64}"

module load Miniconda3/25.11.1-1
source activate ~/.conda/envs/grn

echo "Job $SLURM_JOB_ID on $SLURM_NODELIST, celltype=$CELLTYPE, n_workers=$N_WORKERS"

python scripts/04_run_grnboost2.py \
  --celltype "$CELLTYPE" \
  --n-workers "$N_WORKERS" \
  --threads-per-worker 1
```

---

## 6. results/selected_tfs_mature-endodermis.csv

```csv
TF,summed_importance,rank,family
SbiRTX430.09G017700,74.80855472674872,1,MYB
SbiRTX430.09G017600,25.020342526268546,2,MYB
SbiRTX430.04G243200,21.802530534164845,3,ERF
SbiRTX430.02G238700,20.278892844836765,4,HD-ZIP
SbiRTX430.01G471500,19.72186756301535,5,HD-ZIP
SbiRTX430.06G011000,18.894430971822235,6,C2H2
SbiRTX430.01G110400,17.16389829841509,7,bHLH
SbiRTX430.10G016400,15.063471880957001,8,AP2
SbiRTX430.04G299700,13.257940746761227,9,Dof
SbiRTX430.09G149400,12.238974607374423,10,NAC
```

## results/selected_tfs_dividing_cells.csv

```csv
TF,summed_importance,rank,family
SbiRTX430.10G018700,73.84139128317545,1,GRAS
SbiRTX430.09G219200,54.55574911344838,2,WRKY
SbiRTX430.09G224300,39.35428602761884,3,C2H2
SbiRTX430.06G052500,35.32133161907897,4,FAR1
SbiRTX430.03G428400,20.262860025141826,5,C2H2
SbiRTX430.09G130300,19.059925123967158,6,AP2
SbiRTX430.07G054700,17.96262213100104,7,MYB_related
SbiRTX430.09G085600,17.634857036949914,8,MYB
SbiRTX430.01G110400,15.153436219672361,9,bHLH
SbiRTX430.10G044200,11.418339539845666,10,DBB
```

## results/selected_tfs_exodermis.csv

```csv
TF,summed_importance,rank,family
SbiRTX430.03G038000,74.51320347127212,1,NAC
SbiRTX430.02G190100,27.780739864366733,2,ERF
SbiRTX430.06G011000,26.64202242841211,3,C2H2
SbiRTX430.04G294200,25.086774087661475,4,MYB_related
SbiRTX430.09G017700,24.558865086901456,5,MYB
SbiRTX430.02G211900,23.62153414310835,6,AP2
SbiRTX430.04G243200,19.41450923149384,7,ERF
SbiRTX430.08G120400,17.2565523623674,8,WRKY
SbiRTX430.03G310600,14.402172815887267,9,WRKY
SbiRTX430.09G219200,14.302322358952136,10,WRKY
```

---

## 7. results/tie_block_summary.csv

```csv
celltype,n_total_edges,n_tie_edges,frac_tie_edges,n_top10_bulk_from_ties,top10_bulk_TFs
mature-endodermis,1246780,163986,0.13152761513659186,0,
dividing cells,1563314,81963,0.05242900658472962,0,
exodermis,1380256,163224,0.11825632346463265,0,
```

---

## 8. figures/fig4_flagged_tf_tf_edges.csv

```csv
celltype,TF,target,importance,raw_count_pearson_r
mature-endodermis,SbiRTX430.09G017700,SbiRTX430.09G017600,5.411293792873621,0.5617717035698272
mature-endodermis,SbiRTX430.09G017600,SbiRTX430.09G017700,4.585031965942332,0.5617717035698271
dividing cells,SbiRTX430.07G054700,SbiRTX430.09G224300,1.7308574059452282,0.25699147794500643
exodermis,SbiRTX430.03G310600,SbiRTX430.03G038000,2.1308612056384093,0.3727510372859702
exodermis,SbiRTX430.09G219200,SbiRTX430.03G038000,1.563509802264444,0.339094082714446
```

---

## 9. data/processed/celltype_counts.csv

```csv
"Var1","Freq"
"cortex",1589
"pericycle",1362
"meristem",1354
"atrichoblast",1028
"lignified-xylem",818
"dividing cells",801
"elongating-xylem",684
"s-phase",633
"procambium",605
"early-epidermis",599
"cortical sclerenchyma",575
"phloem",511
"trichoblast",504
"LRC",470
"early-cortex",438
"early-stele",412
"crown root initials",346
"exodermis",315
"mature-endodermis",282
"cortical aerenchyma",281
"phloem/SE",264
"elongating-endodermis",202
"mature-xylem",189
"sclerenchyma",151
"lignified-endodermis",133
"senescent cortex",129
"columella",117
```

---

## Figure inventory

All figures generated by `scripts/06_figures.py`, 300 dpi, saved as both PNG and PDF (PDF is vector; pixel dimensions below are for the PNG raster). Visually verified by direct inspection after the two fig3 rendering bugs (see run log Step 06) were fixed.

### fig1_umap_celltypes.png / .pdf
- **Dimensions:** 2672 × 1788 px
- **Displays:** UMAP embedding (from the original Seurat object's stored `umap` reduction, not recomputed) of all 14,792 singlet nuclei, points colored by `celltype`.
- **Elements:** 14,792 points across 27 celltype color groups, one combined legend (27 entries) placed outside the plot area to the right. Palette: `tab20` + `tab20b` combined (27 of the available 40 colors used).

### fig2_tf_featureplots.png / .pdf
- **Dimensions:** 5980 × 1846 px
- **Displays:** Log-normalized (`normalize_total` target_sum=1e4 + `log1p`) expression of each cell type's top-10 selected TF, plotted on the same shared UMAP coordinates used in fig1.
- **Elements:** 30 panels arranged 3 rows (mature-endodermis, dividing cells, exodermis) × 10 columns (TF rank 1–10 within that cell type), each with its own colorbar. Panel titles show rank, gene ID, and TF family.

### fig3_marker_dotplot.png / .pdf
- **Dimensions:** 3471 × 3288 px
- **Displays:** Dotplot of Table S10 literature root markers across all 27 cell types (`sc.pl.dotplot`, `standard_scale="var"`, dot size = fraction of cells expressing, dot color = mean scaled expression), gene labels from the `alias` column with a gene-ID suffix appended for the one reused alias.
- **Elements:** 46 marker genes (columns) × 27 cell types (rows) = 1,242 dot positions, plus a size legend and a color legend.
- **Exact 46 marker genes plotted (display label → gene ID, table order):**

  | Display label | Gene ID |
  |---|---|
  | DOF43/HCA2 | SbiRTX430.04G278400 |
  | IRX1 | SbiRTX430.03G320500 |
  | XCP1 | SbiRTX430.09G009200 |
  | CESA7 | SbiRTX430.02G210900 |
  | Tdy1 | SbiRTX430.09G148300 |
  | SWEET13 (08G105800) | SbiRTX430.08G105800 |
  | SWEET13 (08G105600) | SbiRTX430.08G105600 |
  | PEAR2 | SbiRTX430.09G015100 |
  | SUT1 | SbiRTX430.01G515100 |
  | CCR1 | SbiRTX430.07G152900 |
  | SHR1 | SbiRTX430.02G363200 |
  | NADP-ME | SbiRTX430.03G039400 |
  | PPDK | SbiRTX430.09G139600 |
  | WBC11 | SbiRTX430.06G163900 |
  | STP1 | SbiRTX430.02G005100 |
  | CKX4 | SbiRTX430.03G040200 |
  | LIGULELESS | SbiRTX430.06G264300 |
  | UBC19 | SbiRTX430.07G198700 |
  | ROOT MERISTEM GROWTH FACTOR 1-RELATED(RMGF1-RELATED) | SbiRTX430.02G240800 |
  | PGX1 | SbiRTX430.03G202200 |
  | LAX2 | SbiRTX430.01G462000 |
  | CYSTM | SbiRTX430.02G067200 |
  | APL | SbiRTX430.04G063900 |
  | XCP | SbiRTX430.03G476300 |
  | MYB46 | SbiRTX430.09G230400 |
  | SDT1 | SbiRTX430.02G396700 |
  | WOX11 | SbiRTX430.02G430600 |
  | SCR1h | SbiRTX430.05G015500 |
  | MYB36 | SbiRTX430.07G098600 |
  | LTP_2 | SbiRTX430.03G431300 |
  | S1EXO | SbiRTX430.02G001700 |
  | RBOH | SbiRTX430.08G131100 |
  | BMR(Brown Midrib 30) | SbiRTX430.01G036400 |
  | HHP1 | SbiRTX430.10G217100 |
  | SILAC3 | SbiRTX430.03G368900 |
  | ASFT | SbiRTX430.03G397000 |
  | AS2 | SbiRTX430.01G428000 |
  | glycine-rich cell wall structural protein 2(GRCWSP-2) | SbiRTX430.01G239000 |
  | WRKY53 | SbiRTX430.07G129000 |
  | ARS1 | SbiRTX430.05G175100 |
  | ARS2 | SbiRTX430.08G039500 |
  | GRP4 | SbiRTX430.05G065000 |
  | ARF16 | SbiRTX430.04G231000 |
  | NAC070/BEARSKIN2 | SbiRTX430.06G097600 |
  | SMB/SOMBRERO | SbiRTX430.04G119700 |
  | MEE31 | SbiRTX430.03G094500 |

- **The 4 omitted rows (of the original 50 in Table S10):**

  | Reason | Alias | Gene ID | Note |
  |---|---|---|---|
  | Exact duplicate row | IRX1 | SbiRTX430.03G320500 | Second occurrence of the same gene ID (row 25 of the raw table); collapsed into the single IRX1 entry above, not a distinct marker. |
  | Not present in RNA/SCT assay | FTL12 | SbiRTX430.10G175100 | Confirmed absent from both assays in the source RDS (verified 2026-08-10). |
  | Not present in RNA/SCT assay | OSC1 | SbiRTX430.08G154200 | Confirmed absent from both assays in the source RDS (verified 2026-08-10). |
  | Not present in RNA/SCT assay | NAC74/Dry | SbiRTX430.06G155500 | Confirmed absent from both assays in the source RDS (verified 2026-08-10). |

### fig4_networks_mature-endodermis.png / .pdf
- **Dimensions:** 3290 × 1537 px
- **Displays:** Bipartite network, mature-endodermis top-10 TFs (left column, colored by family, labeled with rank) → their top-5 marker-directed targets each (right column), edge width ∝ GRNBoost2 importance.
- **Elements:** 10 TF nodes, 50 edges drawn, 35 unique target nodes (some targets shared across TFs). 2 edges flagged crimson-dashed (reciprocal TF-TF pair, `SbiRTX430.09G017700` ↔ `SbiRTX430.09G017600`); those 2 target nodes marked with a star and "[also top-10 TF]" label. Two legends: TF family (7 families present in this cell type's top-10) and edge-style key.

### fig4_networks_dividing_cells.png / .pdf
- **Dimensions:** 3290 × 1537 px
- **Displays:** Same layout as above, for dividing cells' top-10 TFs.
- **Elements:** 10 TF nodes, 50 edges drawn, 34 unique target nodes. 1 edge flagged crimson-dashed (`SbiRTX430.07G054700 → SbiRTX430.09G224300`). Family legend covers the 9 families present in this cell type's top-10.

### fig4_networks_exodermis.png / .pdf
- **Dimensions:** 3290 × 1537 px
- **Displays:** Same layout as above, for exodermis' top-10 TFs.
- **Elements:** 10 TF nodes, 50 edges drawn, 27 unique target nodes. 2 edges flagged crimson-dashed (`SbiRTX430.03G310600 → SbiRTX430.03G038000` and `SbiRTX430.09G219200 → SbiRTX430.03G038000`, both targeting the same rank-1 NAC). Family legend covers the 6 families present in this cell type's top-10.

---

## Script inventory

All scripts in `scripts/`, in pipeline order. Runtimes marked "self-timed" come from timers the script itself prints (exact); others are wall-clock approximations from file/log timestamps and may include some non-execution overhead (noted explicitly).

| Script | Description | Runtime |
|---|---|---|
| `01_smoke_test_grnboost2.py` | Runs GRNBoost2 on a synthetic 50-cell × 20-gene random matrix with 5 TFs (seed 42) to verify the conda environment and all three documented dependency fixes before any production use. | Not self-timed; completed within the tool's default synchronous window (well under 3 minutes), dominated by dask `LocalCluster` startup overhead for a trivial matrix. |
| `02_export_seurat.R` | Loads the source Seurat RDS, auto-detects the singlet/doublet column, subsets to the 14,792 singlet nuclei, and exports RNA raw counts (MatrixMarket) plus full cell metadata and the UMAP embedding (CSV) for Python assembly. | **16.95 s, self-timed** (11.12 s RDS load + 5.83 s processing/export/checks). |
| `02b_assemble_h5ad.py` | Reads the MatrixMarket export and assembles it into a cells × genes AnnData `h5ad` file, attaching the UMAP embedding. | Not self-timed; ran synchronously, well under the 5-minute tool timeout (I/O-bound on ~227 MB MatrixMarket input / 220 MB h5ad output). |
| `03_prepare_inputs.py` | Filters the 1,827-row TF list to genes present in the expression matrix and detected in ≥50 nuclei, writing the final 1,022-gene regulator list and per-TF family table. | Not self-timed; ran synchronously in well under a minute. |
| `04_run_grnboost2.py` | Runs GRNBoost2 for one cell type: subsets the h5ad to that cell type's nuclei, uses the 1,022-gene regulator list as TFs and all genes as targets, seed 42, writes `results/grnboost2_edges_<celltype>.csv`. Invoked by `04_submit.sh` inside each SLURM job, one cell type per job via `--celltype`. | **Self-timed per job (printed by the script, confirmed against `sacct` elapsed):** mature-endodermis 382.0 s, dividing cells 290.6 s, exodermis 355.0 s. |
| `04_submit.sh` | SLURM batch wrapper for `04_run_grnboost2.py`: requests 1 node / 64 cores / 64 GB on `normal_q` with QOS `tc_normal_short`, loads the `grn` conda environment, and launches the Python script with the cell type and worker count passed as arguments. | N/A (submission script; job elapsed times are reported under `04_run_grnboost2.py` above). SLURM job IDs: 6963109 (mature-endodermis probe, 4h limit), 6963283 (dividing cells, 8h limit), 6963284 (exodermis, 4h limit). |
| `05_select_tfs.py` | Computes Wilcoxon marker genes (top 50 by adjusted p-value) for each of the 27 cell types in one call, then for each of the 3 target cell types sums GRNBoost2 importance per TF across edges into those 50 markers, ranks, and writes the top-10 TF table plus the tie-block summary. | Not self-timed; background job wall time approximately a few minutes, dominated by the full 27-group Wilcoxon DE computation on 25,464 genes × 14,792 cells. |
| `05b_verify_importance_and_ties.py` | Diagnostic/audit script: re-derives the top-50 markers and tie-block flags to verify `frac_from_tie_blocks=0` is real (ranks-11–30 control across all 3 cell types, full edge-by-edge audit of the mature-endodermis rank-1 TF), and empirically checks arboreto's importance-denormalization mechanism by summing importance per target across all TFs. | Approximately 4 minutes wall time (script finalized 19:42:03, log completed 19:46:17 — includes a full Wilcoxon DE re-run, not a pure execution timer). |
| `05c_tie_mechanism.py` | Drills into one specific 19-way GRNBoost2 tie block to test whether the tied targets have identical or merely proportional expression vectors, confirming the exact mechanism (byte-for-byte identical single-nucleus detection pattern). | Approximately 55 s wall time (script finalized 19:54:29, log completed 19:55:24 — not a pure execution timer, but this script does no expensive recomputation so overhead is minimal). |
| `06_figures.py` | Generates all four publication figures: UMAP colored by cell type, TF feature plots grouped by cell type, Table S10 marker dotplot, and per-cell-type TF→target bipartite networks with the generic TF-TF edge flag; all at 300 dpi PNG+PDF. | Approximately 4.5 minutes wall time for the final successful run (script finalized 20:09:01, log completed 20:13:45 — includes a full Wilcoxon DE re-run for marker identification plus all figure rendering, not a pure execution timer). |

---

*End of bundle. Assembled 2026-08-10 from the project at `/projects/intro2gds/Manisha_Book_chapter_082026`.*
