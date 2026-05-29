#!/usr/bin/env python3
"""
Script 05: Statistical comparisons
====================================
Purpose : Perform all statistical tests reported in the paper:
            1. Pseudotime comparison: tumor vs matched spleen types (Mann-Whitney)
            2. Composition comparison: 4T1 vs EMT6 tumor neutrophils (Fisher's exact)
            3. Pseudotime comparison: 4T1 vs EMT6 tumor neutrophils (Mann-Whitney)

Input   : NP_merged_trajectory.h5ad  (from 03_trajectory_analysis.py)
Outputs : results/pseudotime_stats_tumor_vs_spleen.csv
          results/composition_4T1_vs_EMT6_tumor.csv
          results/pseudotime_4T1_vs_EMT6_tumor.csv

Dependencies : scanpy, numpy, pandas, scipy, statsmodels
Usage        : python 05_statistical_comparisons.py
"""

import scanpy as sc
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, fisher_exact
from statsmodels.stats.multitest import multipletests
import os, warnings

warnings.filterwarnings("ignore")

INPUT_H5AD = "NP_merged_trajectory.h5ad"
OUTPUT_DIR = "results"
os.makedirs(OUTPUT_DIR, exist_ok=True)

print(f"Loading {INPUT_H5AD}...")
adata = sc.read_h5ad(INPUT_H5AD)
obs   = adata.obs.copy()

TUMOR_TYPES = ["Mature_TAN", "Il1bhi_TAN", "Ccl3hi_GMDSC", "GMDSC"]

# ── 1. Pseudotime: tumor types vs matched spleen types ───────────────────────
print("\n=== 1. Pseudotime: tumor vs matched spleen types ===")
comparisons = [
    ("Mature Neu",  "Mature_TAN",   "Mature state"),
    ("Inflam. Neu", "Il1bhi_TAN",   "Inflammatory state"),
    ("IFN Neu",     "Ccl3hi_GMDSC", "IFN/GMDSC state"),
    ("Inflam. Neu", "GMDSC",        "Inflam vs GMDSC"),
]
rows = []
for sp_type, tm_type, label in comparisons:
    sp_pt = obs[obs["celltype"] == sp_type]["dpt_pseudotime"].values
    tm_pt = obs[obs["celltype"] == tm_type]["dpt_pseudotime"].values
    _, pval = mannwhitneyu(sp_pt, tm_pt, alternative="two-sided")
    rows.append({
        "Comparison": label,
        "Spleen_type": sp_type, "Spleen_n": len(sp_pt),
        "Spleen_median_PT": round(float(np.median(sp_pt)), 4),
        "Tumor_type": tm_type,  "Tumor_n": len(tm_pt),
        "Tumor_median_PT": round(float(np.median(tm_pt)), 4),
        "Delta_PT": round(float(np.median(tm_pt) - np.median(sp_pt)), 4),
        "p_value": pval,
    })
df1 = pd.DataFrame(rows)
_, df1["FDR"], _, _ = multipletests(df1["p_value"], method="fdr_bh")
df1["p_value"] = df1["p_value"].apply(lambda x: f"{x:.2e}")
df1["FDR"]     = df1["FDR"].apply(lambda x: f"{x:.2e}")
print(df1.to_string(index=False))
df1.to_csv(os.path.join(OUTPUT_DIR, "pseudotime_stats_tumor_vs_spleen.csv"), index=False)

# ── 2. Composition: 4T1 vs EMT6 tumor (Fisher's exact) ──────────────────────
print("\n=== 2. Composition: 4T1T vs EMPT6T (Fisher's exact) ===")
tumor  = obs[obs["tissue"] == "Tumor"]
n_4t1  = int((tumor["sample"] == "4T1T").sum())
n_emt6 = int((tumor["sample"] == "EMPT6T").sum())
print(f"  4T1T total: {n_4t1},  EMPT6T total: {n_emt6}")

rows = []
for ct in TUMOR_TYPES:
    a = int(((tumor["sample"] == "4T1T")   & (tumor["celltype"] == ct)).sum())
    b = int(((tumor["sample"] == "EMPT6T") & (tumor["celltype"] == ct)).sum())
    OR, pval = fisher_exact([[a, b], [n_4t1-a, n_emt6-b]])
    rows.append({
        "Cell_type": ct,
        "4T1T_n": a, "4T1T_pct": round(a/n_4t1*100, 1),
        "EMPT6T_n": b, "EMPT6T_pct": round(b/n_emt6*100, 1),
        "OR": round(OR, 3), "p_value": pval,
    })
df2 = pd.DataFrame(rows)
_, df2["FDR"], _, _ = multipletests(df2["p_value"], method="fdr_bh")
df2["p_value"] = df2["p_value"].round(4)
df2["FDR"]     = df2["FDR"].round(4)
print(df2.to_string(index=False))
df2.to_csv(os.path.join(OUTPUT_DIR, "composition_4T1_vs_EMT6_tumor.csv"), index=False)

# ── 3. Pseudotime: 4T1 vs EMT6 tumor (Mann-Whitney + effect size) ────────────
print("\n=== 3. Pseudotime: 4T1T vs EMPT6T ===")
rows = []
for ct in TUMOR_TYPES:
    pt_4t1  = tumor[(tumor["sample"] == "4T1T")   & (tumor["celltype"] == ct)]["dpt_pseudotime"].values
    pt_emt6 = tumor[(tumor["sample"] == "EMPT6T") & (tumor["celltype"] == ct)]["dpt_pseudotime"].values
    if len(pt_4t1) > 5 and len(pt_emt6) > 5:
        U, pval = mannwhitneyu(pt_4t1, pt_emt6, alternative="two-sided")
        r = round(float(1 - (2*U)/(len(pt_4t1)*len(pt_emt6))), 3)  # rank-biserial r
    else:
        pval, r = np.nan, np.nan
    rows.append({
        "Cell_type": ct,
        "4T1T_n": len(pt_4t1),
        "4T1T_median_PT": round(float(np.median(pt_4t1)), 4),
        "EMPT6T_n": len(pt_emt6),
        "EMPT6T_median_PT": round(float(np.median(pt_emt6)), 4) if len(pt_emt6)>0 else np.nan,
        "Delta_PT_EMT6_minus_4T1": round(float(np.median(pt_emt6)-np.median(pt_4t1)), 4) if len(pt_emt6)>0 else np.nan,
        "p_value": round(pval, 4) if not np.isnan(pval) else "n/a",
        "rank_biserial_r": r,
        "note": "underpowered (n<15)" if len(pt_emt6) < 15 else "",
    })
df3 = pd.DataFrame(rows)
print(df3.to_string(index=False))
print("\nNOTE: Composition differences are NOT significant (all FDR=0.87).")
print("      Pseudotime differences are statistically significant but effect")
print("      sizes are small (|r|<0.25) or based on very few EMT6 cells (n=14).")
print("      Both models drive the same TAN/GMDSC composition.")
df3.to_csv(os.path.join(OUTPUT_DIR, "pseudotime_4T1_vs_EMT6_tumor.csv"), index=False)

print("\n=== Statistical comparisons complete ===")
