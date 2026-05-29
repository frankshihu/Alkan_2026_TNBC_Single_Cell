#!/usr/bin/env python3
"""
Script 02: Build AnnData object from Seurat export components
=============================================================
Purpose : Assemble the flat files produced by 01_seurat_export.R into a
          single AnnData (.h5ad) object ready for trajectory analysis.
Input   : export/  directory produced by 01_seurat_export.R
Output  : NP_merged.h5ad  — AnnData with expression matrix + embeddings + metadata

Dependencies : anndata, scipy, numpy, pandas
Usage        : python 02_build_anndata.py
"""

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse
import anndata as ad
import os

# ── Parameters ────────────────────────────────────────────────────────────────
EXPORT_DIR  = "export"           # directory from 01_seurat_export.R
OUTPUT_H5AD = "NP_merged.h5ad"  # output AnnData file

# ── Load sparse count matrix ──────────────────────────────────────────────────
print("Loading count matrix...")
X        = scipy.io.mmread(os.path.join(EXPORT_DIR, "counts_matrix.mtx")).T  # cells x genes
X        = scipy.sparse.csr_matrix(X)
genes    = pd.read_csv(os.path.join(EXPORT_DIR, "genes.csv"))["gene"].values
barcodes = pd.read_csv(os.path.join(EXPORT_DIR, "barcodes.csv"))["barcode"].values
print(f"  Matrix shape: {X.shape}  (cells x genes)")

# ── Load metadata ─────────────────────────────────────────────────────────────
print("Loading metadata...")
meta         = pd.read_csv(os.path.join(EXPORT_DIR, "metadata.csv"), index_col=0)
meta.index   = barcodes

# ── Load embeddings ───────────────────────────────────────────────────────────
print("Loading embeddings...")
harmony = pd.read_csv(os.path.join(EXPORT_DIR, "harmony_embedding.csv"), index_col=0)
umap    = pd.read_csv(os.path.join(EXPORT_DIR, "umap_embedding.csv"),    index_col=0)
pca     = pd.read_csv(os.path.join(EXPORT_DIR, "pca_embedding.csv"),     index_col=0)
harmony.index = umap.index = pca.index = barcodes

# ── Build AnnData ─────────────────────────────────────────────────────────────
print("Building AnnData...")
adata            = ad.AnnData(X=X, obs=meta, var=pd.DataFrame(index=genes))
adata.var_names  = genes
adata.obs_names  = barcodes

# Store embeddings
adata.obsm["X_harmony"] = harmony.values
adata.obsm["X_umap"]    = umap.values
adata.obsm["X_pca"]     = pca.values

# Convenience columns used downstream
adata.obs["celltype"] = adata.obs["celltype.l3"].astype(str)
adata.obs["tissue"]   = adata.obs["tissue"].astype(str)
adata.obs["sample"]   = adata.obs["samples"].astype(str)

# ── Save ──────────────────────────────────────────────────────────────────────
print(f"Saving to {OUTPUT_H5AD}...")
adata.write_h5ad(OUTPUT_H5AD)
size_mb = os.path.getsize(OUTPUT_H5AD) / 1e6
print(f"Done. {OUTPUT_H5AD}  ({size_mb:.1f} MB)")
print(adata)
