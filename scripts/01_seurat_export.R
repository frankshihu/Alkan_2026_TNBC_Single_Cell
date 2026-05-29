#!/usr/bin/env Rscript
# =============================================================================
# Script 01: Export Seurat object to AnnData-compatible components
# =============================================================================
# Purpose : Convert the Harmony-integrated Seurat neutrophil object
#           (NP_merged.rds) into flat files that Python/scanpy can read.
# Input   : NP_merged.rds  — Seurat v5 object, RNA assay, Harmony-integrated
# Outputs : counts_matrix.mtx  — log-normalised sparse count matrix (genes x cells)
#           genes.csv           — gene names
#           barcodes.csv        — cell barcodes
#           harmony_embedding.csv — Harmony embedding (cells x 50 dims)
#           umap_embedding.csv    — UMAP coordinates (cells x 2)
#           pca_embedding.csv     — PCA coordinates (cells x 50)
#           metadata.csv          — per-cell metadata
#
# Dependencies : Seurat >= 5.0, Matrix
# Usage        : Rscript 01_seurat_export.R
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})
setwd ("/Volumes/Shi-Lab/Hasan_Standard/")
# ── Parameters ────────────────────────────────────────────────────────────────
INPUT_RDS  <- "./Seurat/Spleen/outputs/NP_merged.rds"   # path to input Seurat object
OUTPUT_DIR <- "./Seurat/Spleen/outputs/export"          # directory for output files

dir.create(OUTPUT_DIR, showWarnings = FALSE)

# ── Load Seurat object ────────────────────────────────────────────────────────
message("Loading Seurat object: ", INPUT_RDS)
seu <- readRDS(INPUT_RDS)
Idents (seu) <- "celltype.l3"
seu <- RenameIdents(seu, "Ccl3hi_GMDSC" = "Ccl3hi_TAN", "GMDSC" = "Cxcr2hi_TAN")
seu$celltype.l3 <- seu@active.ident
seu$celltype.l3 <- factor (seu$celltype.l3, levels = c("Prolif. preNeu (S)", "Prolif. preNeu (G2M)", "Immature Neu",
                                                       "Act. Neu", "Mature Neu", "IFN Neu", "Inflam. Neu",
                                                       "Mature_TAN", "Il1bhi_TAN", "Ccl3hi_TAN", "Cxcr2hi_TAN"))
message("  Cells: ", ncol(seu), "  Genes: ", nrow(seu))
message("  Assays: ", paste(Assays(seu), collapse = ", "))
message("  Reductions: ", paste(Reductions(seu), collapse = ", "))

# ── Add tissue label ──────────────────────────────────────────────────────────
# Samples ending in "T" are tumor; all others are spleen
seu$tissue <- ifelse(grepl("T$", seu$samples), "Tumor", "Spleen")

# ── Extract log-normalised counts (Seurat v5 layer syntax) ───────────────────
message("Extracting log-normalised counts (RNA / data layer)...")
counts_mat <- GetAssayData(seu, assay = "RNA", layer = "data")
message("  Matrix: ", nrow(counts_mat), " genes x ", ncol(counts_mat), " cells")

# ── Extract embeddings ────────────────────────────────────────────────────────
message("Extracting embeddings...")
harmony_emb <- Embeddings(seu, "harmony")   # 50 dims
umap_emb    <- Embeddings(seu, "umap")      # 2 dims
pca_emb     <- Embeddings(seu, "pca")       # 50 dims

# ── Extract metadata ──────────────────────────────────────────────────────────
meta_cols <- c("orig.ident", "samples", "Phase",
               "celltype.l2", "celltype.l3", "clusters",
               "S.Score", "G2M.Score",
               "nCount_RNA", "nFeature_RNA", "percent.mt")
meta <- seu@meta.data[, meta_cols]
meta$tissue <- seu$tissue

# ── Write outputs ─────────────────────────────────────────────────────────────
message("Writing outputs to: ", OUTPUT_DIR)

# Sparse count matrix in Matrix Market format
writeMM(counts_mat, file.path(OUTPUT_DIR, "counts_matrix.mtx"))
write.csv(data.frame(gene    = rownames(counts_mat)),
          file.path(OUTPUT_DIR, "genes.csv"),    row.names = FALSE)
write.csv(data.frame(barcode = colnames(counts_mat)),
          file.path(OUTPUT_DIR, "barcodes.csv"), row.names = FALSE)

# Embeddings
write.csv(harmony_emb, file.path(OUTPUT_DIR, "harmony_embedding.csv"))
write.csv(umap_emb,    file.path(OUTPUT_DIR, "umap_embedding.csv"))
write.csv(pca_emb,     file.path(OUTPUT_DIR, "pca_embedding.csv"))

# Metadata
write.csv(meta, file.path(OUTPUT_DIR, "metadata.csv"))

message("Done. Files written:")
for (f in list.files(OUTPUT_DIR)) {
  sz <- file.size(file.path(OUTPUT_DIR, f))
  message(sprintf("  %-35s  %.1f MB", f, sz / 1e6))
}
