#!/usr/bin/env python3
"""
Script 03: Neutrophil trajectory analysis — PAGA + DPT pseudotime
==================================================================
Purpose : Run the full trajectory analysis on the combined spleen + tumor
          neutrophil AnnData object:
            1. Build neighbor graph on Harmony embedding (dims 1–40)
            2. Compute diffusion map (15 components)
            3. Compute DPT pseudotime rooted at Prolif. preNeu (S)
            4. Run PAGA cluster connectivity
            5. Export all result tables (CSV)

Input   : NP_merged.h5ad  (from 02_build_anndata.py)
          OR NP_merged_trajectory.h5ad  (pre-computed, skip recomputation)

Outputs : NP_merged_trajectory.h5ad  — AnnData with all trajectory results
          results/cell_metadata_full.csv
          results/pseudotime_summary.csv
          results/paga_connectivity_full.csv
          results/composition_counts.csv
          results/composition_percent.csv
          results/harmony_embedding_40dims.csv

Dependencies : scanpy>=1.9, anndata, numpy, pandas, scipy
Usage        : python 03_trajectory_analysis.py
"""

import scanpy as sc
import numpy as np
import pandas as pd
import scipy.sparse
import os, warnings

warnings.filterwarnings("ignore")
sc.settings.verbosity = 1

# ── Parameters ────────────────────────────────────────────────────────────────
INPUT_H5AD    = "NP_merged.h5ad"           # from 02_build_anndata.py
OUTPUT_H5AD   = "NP_merged_trajectory.h5ad"
OUTPUT_DIR    = "results"
HARMONY_DIMS  = 40     # must match Seurat FindNeighbors dims=1:40
N_NEIGHBORS   = 20     # must match Seurat k.param=20
N_DIFFMAP     = 15     # diffusion map components
N_DPT_DCS     = 10     # DPT uses top N diffusion components
ROOT_CELLTYPE = "Prolif. preNeu (S)"  # trajectory root

# Cell type display order (early → late, spleen then tumor)
CELLTYPE_ORDER = [
    "Prolif. preNeu (S)", "Prolif. preNeu (G2M)", "Immature Neu",
    "Act. Neu", "Mature Neu", "IFN Neu", "Inflam. Neu",
    "Mature_TAN", "Il1bhi_TAN", "Ccl3hi_TAN", "Cxcr2hi_TAN",
]
TUMOR_TYPES  = ["Mature_TAN", "Il1bhi_TAN", "Ccl3hi_TAN", "Cxcr2hi_TAN"]
SPLEEN_TYPES = ["Prolif. preNeu (S)", "Prolif. preNeu (G2M)", "Immature Neu",
                "Act. Neu", "Mature Neu", "IFN Neu", "Inflam. Neu"]

os.makedirs(OUTPUT_DIR, exist_ok=True)

# ── Load data ─────────────────────────────────────────────────────────────────
print(f"Loading {INPUT_H5AD}...")
adata = sc.read_h5ad(INPUT_H5AD)
print(adata)

# Use first HARMONY_DIMS dimensions of Harmony embedding
adata.obsm["X_harmony_40"] = adata.obsm["X_harmony"][:, :HARMONY_DIMS]

# ── Step 1: Neighbor graph ────────────────────────────────────────────────────
print(f"\nStep 1: Building neighbor graph (Harmony dims 1:{HARMONY_DIMS}, "
      f"n_neighbors={N_NEIGHBORS})...")
sc.pp.neighbors(
    adata,
    use_rep="X_harmony_40",
    n_neighbors=N_NEIGHBORS,
    n_pcs=HARMONY_DIMS,
)

# ── Step 2: Diffusion map ─────────────────────────────────────────────────────
print(f"\nStep 2: Computing diffusion map (n_comps={N_DIFFMAP})...")
sc.tl.diffmap(adata, n_comps=N_DIFFMAP)

# ── Step 3: DPT pseudotime ───────────────────────────────────────────────────
print(f"\nStep 3: Computing DPT pseudotime (root = {ROOT_CELLTYPE})...")
# Root: cell with minimum DC1 value within the root cell type
root_mask  = adata.obs["celltype"] == ROOT_CELLTYPE
root_cells = np.where(root_mask)[0]
dc1        = adata.obsm["X_diffmap"][:, 1]   # DC1 (index 1; index 0 is trivial)
adata.uns["iroot"] = int(root_cells[np.argmin(dc1[root_cells])])
print(f"  Root cell index: {adata.uns['iroot']}  "
      f"(DC1={dc1[adata.uns['iroot']]:.4f})")
sc.tl.dpt(adata, n_dcs=N_DPT_DCS)
pt = adata.obs["dpt_pseudotime"]
print(f"  Pseudotime range: [{pt.min():.4f}, {pt.max():.4f}]")

# ── Step 4: PAGA ─────────────────────────────────────────────────────────────
print("\nStep 4: Running PAGA...")
sc.tl.paga(adata, groups="celltype")

# ── Add convenience columns ───────────────────────────────────────────────────
adata.obs["DC1"] = adata.obsm["X_diffmap"][:, 1]
adata.obs["DC2"] = adata.obsm["X_diffmap"][:, 2]

# ── Save processed AnnData ────────────────────────────────────────────────────
print(f"\nSaving processed AnnData to {OUTPUT_H5AD}...")
adata.write_h5ad(OUTPUT_H5AD, compression="gzip")
print(f"  File size: {os.path.getsize(OUTPUT_H5AD)/1e6:.1f} MB")

# ── Export result tables ──────────────────────────────────────────────────────
print("\nExporting result tables...")

# 1. Per-cell metadata + UMAP + pseudotime
meta_out = adata.obs.copy()
meta_out["UMAP_1"] = adata.obsm["X_umap"][:, 0]
meta_out["UMAP_2"] = adata.obsm["X_umap"][:, 1]
meta_out.to_csv(os.path.join(OUTPUT_DIR, "cell_metadata_full.csv"))
print(f"  cell_metadata_full.csv: {meta_out.shape[0]} cells x {meta_out.shape[1]} cols")

# 2. Pseudotime summary per cell type
pt_summary = adata.obs.groupby("celltype")["dpt_pseudotime"].agg(
    n_cells="count", median_PT="median", mean_PT="mean", std_PT="std",
    Q1_PT=lambda x: np.percentile(x, 25),
    Q3_PT=lambda x: np.percentile(x, 75),
    min_PT="min", max_PT="max",
).round(4)
pt_summary["tissue"] = ["Tumor" if c in TUMOR_TYPES else "Spleen"
                         for c in pt_summary.index]
pt_summary = pt_summary.reindex(CELLTYPE_ORDER)
pt_summary.to_csv(os.path.join(OUTPUT_DIR, "pseudotime_summary.csv"))
print(f"  pseudotime_summary.csv")
print(pt_summary[["n_cells", "median_PT", "tissue"]].to_string())

# 3. PAGA connectivity matrix
groups = adata.obs["celltype"].cat.categories.tolist()
conn   = adata.uns["paga"]["connectivities"]
if scipy.sparse.issparse(conn):
    conn = conn.toarray()
paga_df = pd.DataFrame(conn, index=groups, columns=groups)
paga_df.reindex(index=CELLTYPE_ORDER, columns=CELLTYPE_ORDER).round(4).to_csv(
    os.path.join(OUTPUT_DIR, "paga_connectivity_full.csv"))
print(f"  paga_connectivity_full.csv")

# 4. Composition tables
comp     = pd.crosstab(adata.obs["celltype"], adata.obs["sample"])
comp_pct = comp.div(comp.sum(axis=0), axis=1).round(4) * 100
comp.to_csv(os.path.join(OUTPUT_DIR, "composition_counts.csv"))
comp_pct.to_csv(os.path.join(OUTPUT_DIR, "composition_percent.csv"))
print(f"  composition_counts.csv / composition_percent.csv")

# 5. Harmony embedding (40 dims)
harm_df = pd.DataFrame(
    adata.obsm["X_harmony_40"],
    index=adata.obs_names,
    columns=[f"Harmony_{i+1}" for i in range(HARMONY_DIMS)],
)
harm_df.to_csv(os.path.join(OUTPUT_DIR, "harmony_embedding_40dims.csv"))
print(f"  harmony_embedding_40dims.csv")

print("\n=== Trajectory analysis complete ===")
