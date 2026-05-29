# Neutrophil Trajectory Analysis — Reproducibility Package
## Combined Spleen + Tumor Neutrophil PAGA/Pseudotime Analysis

This repository contains all code and data required to reproduce the neutrophil
trajectory analysis from:

> **[Your manuscript title]**  
> [Authors], [Journal], [Year]

---

## Overview

This analysis reconstructs the granulopoiesis trajectory of neutrophils from
spleen (7,376 cells) and tumor (1,803 cells) of 4T1, EMT6, and non-tumor-bearing
mice using PAGA and diffusion pseudotime. The combined dataset contains 9,179
neutrophils across 11 cell types.

**Key findings:**
- Tumor neutrophils (Mature_TAN, Il1bhi_TAN, Ccl3hi_GMDSC, GMDSC) occupy late
  pseudotime positions (median PT = 0.96–0.98), consistent with mature/terminal states
- PAGA reveals distinct lineage relationships: Mature_TAN ↔ Act. Neu;
  Il1bhi_TAN ↔ IFN Neu; GMDSC cluster forms a tight inter-tumor hub
- Both 4T1 and EMT6 models drive identical TAN/GMDSC composition (all FDR = 0.87)

---

## Directory Structure

```
neutrophil_trajectory_analysis/
├── README.md                          # This file
├── requirements.txt                   # Python dependencies
├── r_session_info.txt                 # R session info
│
├── data/
│   ├── NP_merged.rds                  # INPUT: Harmony-integrated Seurat object
│   │                                  # (not included — provide your own)
│   └── NP_merged_trajectory.h5ad     # OUTPUT: Processed AnnData (provided)
│
├── scripts/
│   ├── 01_seurat_export.R             # Step 1: Export Seurat → flat files
│   ├── 02_build_anndata.py            # Step 2: Build AnnData from flat files
│   ├── 03_trajectory_analysis.py      # Step 3: PAGA + pseudotime
│   ├── 04_generate_figures.py         # Step 4: All figures
│   └── 05_statistical_comparisons.py  # Step 5: Statistical tests
│
├── results/                           # CSV result tables (provided)
│   ├── cell_metadata_full.csv
│   ├── pseudotime_summary.csv
│   ├── paga_connectivity_full.csv
│   ├── composition_counts.csv
│   ├── composition_percent.csv
│   ├── harmony_embedding_40dims.csv
│   ├── pseudotime_stats_tumor_vs_spleen.csv
│   ├── composition_4T1_vs_EMT6_tumor.csv
│   └── pseudotime_4T1_vs_EMT6_tumor.csv
│
└── figures/                           # All figures (provided)
    ├── 01_umap_celltype.{png,svg}
    ├── 02_umap_split_tissue.{png,svg}
    ├── 03_umap_pseudotime.{png,svg}
    ├── 04_umap_split_sample.{png,svg}
    ├── 05_paga_graph.{png,svg}
    ├── 06_paga_heatmap.{png,svg}
    ├── 07_pseudotime_violin.{png,svg}
    ├── 08_marker_genes_umap.{png,svg}
    └── 09_comprehensive_summary.{png,svg}
```

---

## Quick Start: Reproduce from Processed Data

If you have `NP_merged_trajectory.h5ad`, you can skip Steps 1–3 and go
directly to figures and statistics:

```bash
# Install Python dependencies
pip install -r requirements.txt

# Generate all figures
python scripts/04_generate_figures.py

# Run statistical comparisons
python scripts/05_statistical_comparisons.py
```

---

## Full Reproduction from Raw Seurat Object

### Prerequisites

**R (≥ 4.4.0):**
```r
install.packages("Seurat")   # v5.5.0
install.packages("Matrix")   # v1.7.5
```

**Python (≥ 3.11):**
```bash
pip install -r requirements.txt
```

### Step-by-step

**Step 1 — Export Seurat object to flat files (R)**
```bash
Rscript scripts/01_seurat_export.R
# Input:  NP_merged.rds
# Output: export/  (counts_matrix.mtx, genes.csv, barcodes.csv,
#                   harmony_embedding.csv, umap_embedding.csv,
#                   pca_embedding.csv, metadata.csv)
```

**Step 2 — Build AnnData (Python)**
```bash
python scripts/02_build_anndata.py
# Input:  export/
# Output: NP_merged.h5ad
```

**Step 3 — Trajectory analysis: PAGA + pseudotime (Python)**
```bash
python scripts/03_trajectory_analysis.py
# Input:  NP_merged.h5ad
# Output: NP_merged_trajectory.h5ad
#         results/  (CSV tables)
```

**Step 4 — Generate all figures (Python)**
```bash
python scripts/04_generate_figures.py
# Input:  NP_merged_trajectory.h5ad
# Output: figures/  (PNG + SVG)
```

**Step 5 — Statistical comparisons (Python)**
```bash
python scripts/05_statistical_comparisons.py
# Input:  NP_merged_trajectory.h5ad
# Output: results/  (statistical test CSVs)
```

---

## Key Parameters

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| Harmony dims used | 1–40 | Matches Seurat `FindNeighbors(dims=1:40)` |
| n_neighbors | 20 | Matches Seurat `k.param=20` |
| Diffusion components | 15 | Sufficient to capture trajectory structure |
| DPT n_dcs | 10 | Top 10 diffusion components for pseudotime |
| Pseudotime root | Prolif. preNeu (S) | Earliest progenitor; cell with min DC1 |
| PAGA grouping | celltype (celltype.l3) | 11 groups: 7 spleen + 4 tumor |

---

## Cell Type Annotations

| Cell type | Tissue | n cells | Description |
|-----------|--------|---------|-------------|
| Prolif. preNeu (S) | Spleen | 467 | S-phase proliferating precursor |
| Prolif. preNeu (G2M) | Spleen | 343 | G2M-phase proliferating precursor |
| Immature Neu | Spleen | 1,264 | Post-mitotic immature neutrophil |
| Act. Neu | Spleen | 561 | Activated neutrophil |
| Mature Neu | Spleen | 1,227 | Mature circulating neutrophil |
| IFN Neu | Spleen | 363 | IFN-stimulated terminal state |
| Inflam. Neu | Spleen | 3,151 | Inflammatory terminal state |
| Mature_TAN | Tumor | 198 | Mature tumor-associated neutrophil |
| Il1bhi_TAN | Tumor | 553 | IL-1β-high inflammatory TAN |
| Ccl3hi_GMDSC | Tumor | 773 | CCL3-high G-MDSC |
| GMDSC | Tumor | 279 | Immunosuppressive G-MDSC |

---

## Software Versions

| Software | Version |
|----------|---------|
| Python | 3.11.14 |
| scanpy | 1.11.4 |
| anndata | 0.12.1 |
| numpy | 2.1.0 |
| pandas | 2.3.3 |
| scipy | 1.15.0 |
| matplotlib | 3.10.9 |
| seaborn | 0.13.2 |
| statsmodels | 0.14.6 |
| R | 4.4.2 |
| Seurat | 5.5.0 |
| Matrix | 1.7.5 |

---

## Statistical Methods

- **PAGA connectivity**: Graph abstraction of cluster-level connectivity
  (Wolf et al., *Genome Biol.* 2019)
- **Diffusion pseudotime (DPT)**: Haghverdi et al., *Nat. Methods* 2016
- **Pseudotime comparisons**: Mann-Whitney U test, FDR correction (Benjamini-Hochberg)
- **Composition comparisons**: Fisher's exact test, FDR correction (Benjamini-Hochberg)
- **Effect size**: Rank-biserial correlation r (r=0.1 small, 0.3 medium, 0.5 large)

---

## Important Caveats

1. **EMT6 tumor cells are sparse** (n=158 total; n=14–69 per cell type).
   Composition comparisons are adequately powered (can detect OR≥1.5–2.0),
   but pseudotime and DE comparisons between 4T1 and EMT6 tumor cells should
   be treated as exploratory.

2. **Tumor cells form a distinct UMAP island** despite Harmony integration,
   reflecting genuine TME-specific transcriptional programs. PAGA connectivity
   (graph-based) is more reliable than UMAP distances for trajectory inference.

3. **No RNA velocity**: Spliced/unspliced counts are not available in this
   object. Pseudotime directionality is inferred from diffusion pseudotime
   rooted at the proliferating precursor population.

---

## Contact

[Your name and contact information]

---

## License

[Your license]
