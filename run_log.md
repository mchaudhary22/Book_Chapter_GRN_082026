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
