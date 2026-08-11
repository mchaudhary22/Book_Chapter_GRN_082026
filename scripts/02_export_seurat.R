#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

set.seed(42)

rds_path <- "data/raw/GSE297576_seurat_object.bicolor_root_atlas.RDS"
out_dir <- "data/processed"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

t0 <- Sys.time()
cat("=== Loading RDS ===\n")
obj <- readRDS(rds_path)
cat("Load time:", format(Sys.time() - t0), "\n")

cat("\n=== Object structure ===\n")
cat("Assays:", paste(Assays(obj), collapse = ", "), "\n")
cat("DefaultAssay (as loaded):", DefaultAssay(obj), "\n")
cat("Reductions:", paste(Reductions(obj), collapse = ", "), "\n")
cat("Metadata columns:\n")
print(colnames(obj@meta.data))

DefaultAssay(obj) <- "RNA"
rna_layers <- Layers(obj[["RNA"]])
cat("\nRNA layers:", paste(rna_layers, collapse = ", "), "\n")
if (length(grep("^counts", rna_layers)) > 1) {
  cat("Multiple counts layers detected — joining layers\n")
  obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
}

cat("\n=== BEFORE FILTERING ===\n")
n_cells_before <- ncol(obj)
n_features_before <- nrow(obj[["RNA"]])
cat("Cells:", n_cells_before, "\n")
cat("RNA features:", n_features_before, "(expected 25464)\n")
if (n_features_before != 25464) {
  stop("FATAL: RNA feature count ", n_features_before, " != expected 25464. Stopping.")
}

# --- Identify the doublet column programmatically (do not assume its name) ---
meta <- obj@meta.data
candidate_cols <- c()
for (col in colnames(meta)) {
  vals <- unique(as.character(meta[[col]]))
  if ("Singlet" %in% vals) {
    candidate_cols <- c(candidate_cols, col)
  }
}
cat("\nCandidate doublet columns (contain a 'Singlet' value):\n")
print(candidate_cols)
if (length(candidate_cols) != 1) {
  stop("FATAL: expected exactly one doublet column containing 'Singlet', found ",
       length(candidate_cols), ": ", paste(candidate_cols, collapse = ", "))
}
doublet_col <- candidate_cols[1]
cat("Using doublet column:", doublet_col, "\n")
cat("Value counts:\n")
print(table(meta[[doublet_col]], useNA = "always"))

if (!"celltype" %in% colnames(meta)) {
  stop("FATAL: 'celltype' column not found in metadata.")
}
n_celltypes <- length(unique(meta$celltype))
cat("\ncelltype has", n_celltypes, "levels (expected 27)\n")
if (n_celltypes != 27) {
  stop("FATAL: celltype level count ", n_celltypes, " != expected 27. Stopping.")
}

# --- Subset to singlets ---
singlet_cells <- rownames(meta)[meta[[doublet_col]] == "Singlet"]
cat("\nSinglet cell count:", length(singlet_cells), "(expected 14792)\n")
if (length(singlet_cells) != 14792) {
  stop("FATAL: singlet count ", length(singlet_cells),
       " != expected 14792. Stopping per pipeline rules — do not proceed on mismatch.")
}

obj_sub <- subset(obj, cells = singlet_cells)

cat("\n=== AFTER FILTERING ===\n")
n_cells_after <- ncol(obj_sub)
n_features_after <- nrow(obj_sub[["RNA"]])
cat("Cells:", n_cells_after, "\n")
cat("RNA features:", n_features_after, "\n")

cat("\nPer-cell-type counts (all", length(unique(obj_sub$celltype)), "types):\n")
ct_counts <- sort(table(obj_sub$celltype), decreasing = TRUE)
print(ct_counts)
write.csv(as.data.frame(ct_counts), file.path(out_dir, "celltype_counts.csv"), row.names = FALSE)

targets <- c("exodermis" = 315, "dividing cells" = 801, "mature-endodermis" = 282)
cat("\nTarget cell-type checks:\n")
for (ct in names(targets)) {
  actual <- sum(obj_sub$celltype == ct)
  cat(sprintf("  %s: %d (expected %d)\n", ct, actual, targets[[ct]]))
  if (actual != targets[[ct]]) {
    stop("FATAL: cell type '", ct, "' count ", actual, " != expected ", targets[[ct]], ". Stopping.")
  }
}

# --- Extract raw RNA counts ---
counts <- GetAssayData(obj_sub, assay = "RNA", layer = "counts")
cat("\nCounts matrix dims (features x cells):", paste(dim(counts), collapse = " x "), "\n")
cat("Class:", class(counts)[1], "\n")

cat("\nWriting MatrixMarket + metadata to", out_dir, "...\n")
Matrix::writeMM(counts, file.path(out_dir, "counts.mtx"))
write.csv(data.frame(gene = rownames(counts)), file.path(out_dir, "genes.csv"), row.names = FALSE)
write.csv(data.frame(cell = colnames(counts)), file.path(out_dir, "cells.csv"), row.names = FALSE)
write.csv(obj_sub@meta.data, file.path(out_dir, "metadata.csv"), row.names = TRUE)

if ("umap" %in% Reductions(obj_sub)) {
  umap_emb <- Embeddings(obj_sub, "umap")
  write.csv(umap_emb, file.path(out_dir, "umap_embedding.csv"), row.names = TRUE)
  cat("Wrote UMAP embedding, dims:", paste(dim(umap_emb), collapse = " x "), "\n")
} else {
  cat("NOTE: no 'umap' reduction found on the object (Reductions: ",
      paste(Reductions(obj_sub), collapse = ", "), "). Step 06 UMAP will need to be computed.\n")
}

cat("\n=== SUMMARY ===\n")
cat("Before filtering:", n_cells_before, "cells x", n_features_before, "RNA features\n")
cat("After filtering: ", n_cells_after, "cells x", n_features_after, "RNA features\n")
cat("Total elapsed:", format(Sys.time() - t0), "\n")
cat("Route used: MatrixMarket + Python assembly (SeuratDisk not available per Step 00).\n")
