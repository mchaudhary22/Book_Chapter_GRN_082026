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
