#!/usr/bin/env python3
"""
Script 04: Generate all trajectory figures
==========================================
Purpose : Produce all publication-quality figures from the processed AnnData.
          Figures are saved as both PNG (300 DPI) and SVG (vector, text editable).

Input   : NP_merged_trajectory.h5ad  (from 03_trajectory_analysis.py)
Outputs : figures/
            01_umap_celltype.{png,svg}
            02_umap_split_tissue.{png,svg}
            03_umap_pseudotime.{png,svg}
            04_umap_split_sample.{png,svg}
            05_paga_graph.{png,svg}
            06_paga_heatmap.{png,svg}
            07_pseudotime_violin.{png,svg}
            08_marker_genes_umap.{png,svg}
            09_comprehensive_summary.{png,svg}

Dependencies : scanpy, matplotlib, seaborn, numpy, pandas, scipy
Usage        : python 04_generate_figures.py
"""

import scanpy as sc
import numpy as np
import pandas as pd
import scipy.sparse
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
import seaborn as sns
import os, warnings

warnings.filterwarnings("ignore")
matplotlib.rcParams["font.family"]   = ["Liberation Sans", "Arial", "DejaVu Sans"]
matplotlib.rcParams["svg.fonttype"]  = "none"   # keep text editable in Illustrator

# ── Parameters ────────────────────────────────────────────────────────────────
INPUT_H5AD  = "NP_merged_trajectory.h5ad"
OUTPUT_DIR  = "figures"
DPI         = 300

os.makedirs(OUTPUT_DIR, exist_ok=True)

def savefig(name):
    """Save current figure as PNG + SVG."""
    base = os.path.join(OUTPUT_DIR, name)
    plt.savefig(base + ".png", dpi=DPI, bbox_inches="tight")
    plt.savefig(base + ".svg",          bbox_inches="tight")
    plt.close()
    print(f"  Saved: {name}.png / .svg")

# ── Color palettes ────────────────────────────────────────────────────────────
SPLEEN_COLORS = {
    "Prolif. preNeu (S)":   "#2166AC",
    "Prolif. preNeu (G2M)": "#4393C3",
    "Immature Neu":         "#74ADD1",
    "Act. Neu":             "#ABD9E9",
    "Mature Neu":           "#4DAC26",
    "IFN Neu":              "#7B3294",
    "Inflam. Neu":          "#1A9641",
}
TUMOR_COLORS = {
    "Mature_TAN":    "#D73027",
    "Il1bhi_TAN":    "#F46D43",
    "Ccl3hi_TAN":  "#FDAE61",
    "Cxcr2hi_TAN":         "#FEE08B",
}
ALL_COLORS = {**SPLEEN_COLORS, **TUMOR_COLORS}

CELLTYPE_ORDER = [
    "Prolif. preNeu (S)", "Prolif. preNeu (G2M)", "Immature Neu",
    "Act. Neu", "Mature Neu", "IFN Neu", "Inflam. Neu",
    "Mature_TAN", "Il1bhi_TAN", "Ccl3hi_TAN", "Cxcr2hi_TAN",
]
TUMOR_TYPES = ["Mature_TAN", "Il1bhi_TAN", "Ccl3hi_TAN", "Cxcr2hi_TAN"]

SAMPLE_COLORS = {
    "N_S":    "#636363",
    "4T1S":   "#2166AC",
    "EMPT6S": "#74ADD1",
    "4T1T":   "#D73027",
    "EMPT6T": "#F46D43",
}
SAMPLE_LABELS = {
    "N_S":    "Non-tumor spleen",
    "4T1S":   "4T1 spleen",
    "EMPT6S": "EMT6 spleen",
    "4T1T":   "4T1 tumor",
    "EMPT6T": "EMT6 tumor",
}

# ── Load data ─────────────────────────────────────────────────────────────────
print(f"Loading {INPUT_H5AD}...")
adata = sc.read_h5ad(INPUT_H5AD)
umap  = adata.obsm["X_umap"]
print(adata)

# ── Figure 1: UMAP by cell type ───────────────────────────────────────────────
print("\nGenerating figures...")
fig, ax = plt.subplots(figsize=(8, 6))
colors_all = [ALL_COLORS.get(c, "#AAAAAA") for c in adata.obs["celltype"]]
ax.scatter(umap[:, 0], umap[:, 1], c=colors_all, s=5, alpha=0.6, rasterized=True)
handles = [mpatches.Patch(color=ALL_COLORS[c], label=c) for c in CELLTYPE_ORDER]
ax.legend(handles=handles, fontsize=6, loc="lower left", frameon=True, framealpha=0.8)
ax.set_title("Cell type (combined spleen + tumor)", fontsize=12, fontweight="bold")
ax.axis("off")
plt.tight_layout()
savefig("01_umap_celltype")

# ── Figure 2: UMAP split by tissue ───────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(13, 5))
for ax, tissue in zip(axes, ["Spleen", "Tumor"]):
    mask = adata.obs["tissue"] == tissue
    ax.scatter(umap[:, 0], umap[:, 1], c="#EEEEEE", s=4, rasterized=True, zorder=1)
    fg_colors = [ALL_COLORS.get(c, "#AAAAAA") for c in adata.obs["celltype"][mask]]
    ax.scatter(umap[mask, 0], umap[mask, 1], c=fg_colors, s=8,
               rasterized=True, zorder=2, alpha=0.8)
    ax.set_title(f"{tissue} (n={mask.sum():,})", fontsize=12, fontweight="bold")
    shown = adata.obs["celltype"][mask].unique()
    handles = [mpatches.Patch(color=ALL_COLORS.get(c, "#AAAAAA"), label=c) for c in shown]
    ax.legend(handles=handles, fontsize=6.5, loc="lower left", frameon=True)
    ax.axis("off")
plt.suptitle("Neutrophil UMAP — Spleen vs Tumor", fontsize=13, fontweight="bold")
plt.tight_layout()
savefig("02_umap_split_tissue")

# ── Figure 3: UMAP pseudotime ─────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(13, 5))
sc.pl.umap(adata, color="dpt_pseudotime", ax=axes[0], show=False,
           frameon=False, size=8, color_map="viridis",
           title="Pseudotime (all cells)")
tumor_mask = adata.obs["tissue"] == "Tumor"
sc_bg = axes[1].scatter(umap[:, 0], umap[:, 1],
                         c=adata.obs["dpt_pseudotime"], cmap="viridis",
                         s=4, alpha=0.5, rasterized=True)
axes[1].scatter(umap[tumor_mask, 0], umap[tumor_mask, 1],
                c=adata.obs["dpt_pseudotime"][tumor_mask],
                cmap="viridis", s=20, rasterized=True, zorder=3,
                edgecolors="black", linewidths=0.3)
axes[1].set_title("Pseudotime (tumor cells highlighted)", fontsize=11)
axes[1].axis("off")
plt.colorbar(sc_bg, ax=axes[1], label="Pseudotime", shrink=0.7)
plt.tight_layout()
savefig("03_umap_pseudotime")

# ── Figure 4: UMAP split by sample ───────────────────────────────────────────
samples = ["N_S", "4T1S", "EMPT6S", "4T1T", "EMPT6T"]
fig, axes = plt.subplots(1, 5, figsize=(22, 4))
for ax, samp in zip(axes, samples):
    mask = adata.obs["sample"] == samp
    ax.scatter(umap[:, 0], umap[:, 1], c="#EEEEEE", s=3, rasterized=True, zorder=1)
    ax.scatter(umap[mask, 0], umap[mask, 1], c=SAMPLE_COLORS[samp],
               s=8, rasterized=True, zorder=2, alpha=0.8)
    ax.set_title(f"{SAMPLE_LABELS[samp]}\n(n={mask.sum():,})", fontsize=9, fontweight="bold")
    ax.axis("off")
plt.suptitle("UMAP split by sample", fontsize=13, fontweight="bold")
plt.tight_layout()
savefig("04_umap_split_sample")

# ── Figure 5: PAGA graph (schematic layout) ───────────────────────────────────
SCHEMATIC_POS = {
    "Prolif. preNeu (S)":   (0.0,  0.0),
    "Prolif. preNeu (G2M)": (0.0,  1.5),
    "Immature Neu":         (0.0,  3.0),
    "Act. Neu":             (0.0,  4.5),
    "Mature Neu":           (-1.5, 6.0),
    "IFN Neu":              (-1.5, 7.5),
    "Inflam. Neu":          (-1.5, 9.0),
    "Mature_TAN":           (2.5,  6.0),
    "Il1bhi_TAN":           (4.5,  7.5),
    "Ccl3hi_TAN":         (4.5,  4.5),
    "Cxcr2hi_TAN":                (4.5,  3.0),
}
LABEL_CONFIG = {
    "Prolif. preNeu (S)":   (-0.6,  0.0,  "right"),
    "Prolif. preNeu (G2M)": (-0.6,  1.5,  "right"),
    "Immature Neu":         (-0.6,  3.0,  "right"),
    "Act. Neu":             (-0.6,  4.5,  "right"),
    "Mature Neu":           (-2.1,  6.0,  "right"),
    "IFN Neu":              (-2.1,  7.5,  "right"),
    "Inflam. Neu":          (-2.1,  9.0,  "right"),
    "Mature_TAN":           (3.1,   6.0,  "left"),
    "Il1bhi_TAN":           (5.1,   7.5,  "left"),
    "Ccl3hi_TAN":         (5.1,   4.5,  "left"),
    "Cxcr2hi_TAN":                (5.1,   3.0,  "left"),
}

groups = adata.obs["celltype"].cat.categories.tolist()
conn   = adata.uns["paga"]["connectivities"]
if scipy.sparse.issparse(conn):
    conn = conn.toarray()

fig, ax = plt.subplots(figsize=(12, 11))
threshold = 0.05
for i, g1 in enumerate(groups):
    for j, g2 in enumerate(groups):
        if j <= i:
            continue
        w = conn[i, j]
        if w > threshold:
            p1, p2 = SCHEMATIC_POS[g1], SCHEMATIC_POS[g2]
            color  = "#111111" if w >= 0.5 else "#888888"
            ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color=color,
                    lw=w * 10, alpha=min(0.85, 0.3 + w * 0.7),
                    zorder=2, solid_capstyle="round")
            if w >= 0.5:
                mx, my = (p1[0]+p2[0])/2, (p1[1]+p2[1])/2
                ax.text(mx, my, f"{w:.2f}", fontsize=7, ha="center", va="center",
                        color="white", fontweight="bold", zorder=6,
                        bbox=dict(fc="#333333", alpha=0.75, ec="none",
                                  boxstyle="round,pad=0.15"))
node_r = 0.42
for g in groups:
    p  = SCHEMATIC_POS[g]
    c  = ALL_COLORS.get(g, "#AAAAAA")
    ec = "#CC0000" if g in TUMOR_TYPES else "#333333"
    lw = 2.5 if g in TUMOR_TYPES else 1.2
    ax.add_patch(plt.Circle(p, radius=node_r, color=c, ec=ec, lw=lw, zorder=4))
for g, (lx, ly, ha) in LABEL_CONFIG.items():
    ax.text(lx, ly, g, ha=ha, va="center", fontsize=9.5, fontweight="bold",
            zorder=7, color=ALL_COLORS.get(g, "#000000"),
            bbox=dict(fc="white", alpha=0.85, ec="none", boxstyle="round,pad=0.2"))
ax.text(-2.2, 10.2, "SPLEEN", fontsize=13, fontweight="bold", color="#2166AC", ha="center")
ax.text(4.5,  10.2, "TUMOR",  fontsize=13, fontweight="bold", color="#D73027", ha="center")
ax.axvline(1.3, color="#AAAAAA", lw=1.5, ls="--", alpha=0.5, zorder=1)
legend_handles = [
    mpatches.Patch(color="#2166AC", label="Spleen neutrophils"),
    mpatches.Patch(facecolor="#D73027", edgecolor="#CC0000", lw=2, label="Tumor neutrophils"),
    plt.Line2D([0],[0], color="#111111", lw=4, label="Strong connection (≥0.5)"),
    plt.Line2D([0],[0], color="#888888", lw=1.5, label="Weak connection (0.05–0.5)"),
]
ax.legend(handles=legend_handles, loc="lower right", fontsize=9, frameon=True)
ax.set_title("PAGA connectivity graph — Combined spleen + tumor neutrophils\n"
             "(edge width ∝ connectivity; scores on strong edges ≥0.5)",
             fontsize=12, fontweight="bold")
ax.set_xlim(-3.5, 6.5); ax.set_ylim(-1.0, 11.0)
ax.set_aspect("equal"); ax.axis("off")
plt.tight_layout()
savefig("05_paga_graph")

# ── Figure 6: PAGA connectivity heatmap ──────────────────────────────────────
paga_df = pd.DataFrame(conn, index=groups, columns=groups)
paga_ordered = paga_df.reindex(index=CELLTYPE_ORDER, columns=CELLTYPE_ORDER)
fig, ax = plt.subplots(figsize=(9, 8))
mask_diag = np.eye(len(CELLTYPE_ORDER), dtype=bool)
sns.heatmap(paga_ordered, ax=ax, cmap="YlOrRd", vmin=0, vmax=1,
            mask=mask_diag, annot=True, fmt=".2f", annot_kws={"size": 7},
            linewidths=0.5, linecolor="#CCCCCC",
            cbar_kws={"label": "PAGA connectivity", "shrink": 0.7})
for tick, label in zip(ax.get_xticklabels(), CELLTYPE_ORDER):
    tick.set_color(ALL_COLORS.get(label, "#000000")); tick.set_fontweight("bold"); tick.set_fontsize(8)
for tick, label in zip(ax.get_yticklabels(), CELLTYPE_ORDER):
    tick.set_color(ALL_COLORS.get(label, "#000000")); tick.set_fontweight("bold"); tick.set_fontsize(8)
ax.axhline(7, color="black", lw=2.5, ls="--")
ax.axvline(7, color="black", lw=2.5, ls="--")
ax.text(3.5, -0.6, "Spleen", ha="center", fontsize=10, fontweight="bold", color="#2166AC")
ax.text(9.0, -0.6, "Tumor",  ha="center", fontsize=10, fontweight="bold", color="#D73027")
ax.set_title("PAGA connectivity — Spleen + Tumor neutrophils", fontsize=12, fontweight="bold")
plt.xticks(rotation=45, ha="right"); plt.yticks(rotation=0)
plt.tight_layout()
savefig("06_paga_heatmap")

# ── Figure 7: Pseudotime violin ───────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(13, 5))
plot_data = adata.obs[["celltype", "dpt_pseudotime"]].copy()
plot_data  = plot_data[plot_data["celltype"].isin(CELLTYPE_ORDER)]
sns.violinplot(data=plot_data, x="celltype", y="dpt_pseudotime",
               order=CELLTYPE_ORDER, palette=ALL_COLORS, ax=ax,
               inner="box", cut=0, linewidth=1.2, scale="width")
medians = plot_data.groupby("celltype")["dpt_pseudotime"].median()
for i, ct in enumerate(CELLTYPE_ORDER):
    if ct in medians:
        ax.scatter(i, medians[ct], color="white", s=30, zorder=5,
                   edgecolors="black", lw=1)
ax.axvline(6.5, color="black", lw=2, ls="--", alpha=0.7)
ax.text(3.0, 1.07, "Spleen neutrophils", ha="center", fontsize=10,
        fontweight="bold", color="#2166AC")
ax.text(9.0, 1.07, "Tumor neutrophils",  ha="center", fontsize=10,
        fontweight="bold", color="#D73027")
ax.set_xlabel(""); ax.set_ylabel("Diffusion pseudotime", fontsize=11)
ax.set_title("Pseudotime distribution by cell type", fontsize=12, fontweight="bold")
plt.xticks(rotation=35, ha="right", fontsize=9)
ax.set_ylim(-0.05, 1.12)
ax.spines[["top", "right"]].set_visible(False)
plt.tight_layout()
savefig("07_pseudotime_violin")

# ── Figure 8: Key marker genes on UMAP ───────────────────────────────────────
markers = ["Retnlg", "Wfdc17", "Ngp", "Camp",
           "Il1b", "Ccl3", "Zeb2", "Slfn4",
           "Cd274", "Igf1r", "Slc7a11", "Ifit1"]
present = [g for g in markers if g in adata.var_names]
fig, axes = plt.subplots(3, 4, figsize=(16, 12))
for idx, gene in enumerate(present[:12]):
    ax = axes.flatten()[idx]
    expr = np.array(adata[:, gene].X.todense()).flatten()
    sc_p = ax.scatter(umap[:, 0], umap[:, 1], c=expr, cmap="Reds", s=4,
                      rasterized=True, vmin=0, vmax=np.percentile(expr, 99))
    ax.set_title(gene, fontsize=10, fontweight="bold", fontstyle="italic")
    ax.axis("off")
    plt.colorbar(sc_p, ax=ax, shrink=0.7, pad=0.02)
for idx in range(len(present[:12]), 12):
    axes.flatten()[idx].axis("off")
plt.suptitle("Key marker gene expression on combined UMAP",
             fontsize=13, fontweight="bold")
plt.tight_layout()
savefig("08_marker_genes_umap")

# ── Figure 9: Comprehensive manuscript summary ────────────────────────────────
# Revised layout requested by user:
#   Top row:    A = cell type UMAP; B = split by tissue, spanning two columns
#   Bottom row: C = pseudotime UMAP; D = pseudotime violin; E = PAGA connection map
#   Panel C from the previous version was removed.
#   The figure is optimized for US-letter landscape manuscript placement.

fig = plt.figure(figsize=(11.0, 8.5))  # US letter landscape-compatible
main_gs = gridspec.GridSpec(
    2, 3, figure=fig,
    left=0.055, right=0.985, bottom=0.22, top=0.90,
    hspace=0.34, wspace=0.26,
    width_ratios=[1.05, 1.15, 1.15],
    height_ratios=[1.0, 1.0]
)

PANEL_LABEL_FS = 17
TITLE_FS       = 11.5
AXIS_FS        = 10
TICK_FS        = 8.5
LEGEND_FS      = 7.4
POINT_SIZE     = 3.4


def add_panel_label(ax, label, x=-0.12, y=1.08):
    ax.text(
        x, y, label,
        transform=ax.transAxes,
        fontsize=PANEL_LABEL_FS,
        fontweight="bold",
        va="top", ha="left"
    )


def clean_umap_axis(ax):
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("")
    ax.set_ylabel("")
    for spine in ax.spines.values():
        spine.set_visible(False)


# A: cell type UMAP ------------------------------------------------------------
ax_a = fig.add_subplot(main_gs[0, 0])
ax_a.scatter(
    umap[:, 0], umap[:, 1],
    c=[ALL_COLORS.get(c, "#AAAAAA") for c in adata.obs["celltype"]],
    s=POINT_SIZE, alpha=0.65, rasterized=True, linewidths=0
)
add_panel_label(ax_a, "A")
ax_a.set_title("Cell type", fontsize=TITLE_FS, fontweight="bold", pad=4)
clean_umap_axis(ax_a)


# B: UMAP split by tissue, spanning top middle/right --------------------------
# Use a nested two-panel layout so spleen and tumor have enough width.
gs_b = gridspec.GridSpecFromSubplotSpec(
    1, 2, subplot_spec=main_gs[0, 1:], wspace=0.10
)
ax_b1 = fig.add_subplot(gs_b[0, 0])
ax_b2 = fig.add_subplot(gs_b[0, 1])

for ax, tissue in zip([ax_b1, ax_b2], ["Spleen", "Tumor"]):
    mask = adata.obs["tissue"] == tissue
    ax.scatter(
        umap[:, 0], umap[:, 1],
        c="#E6E6E6", s=2.4, alpha=0.55,
        rasterized=True, linewidths=0, zorder=1
    )
    ax.scatter(
        umap[mask, 0], umap[mask, 1],
        c=[ALL_COLORS.get(c, "#AAAAAA") for c in adata.obs["celltype"][mask]],
        s=4.8, alpha=0.86,
        rasterized=True, linewidths=0, zorder=2
    )
    ax.set_title(f"{tissue}\n(n={mask.sum():,})", fontsize=10.5, fontweight="bold", pad=2)
    clean_umap_axis(ax)

add_panel_label(ax_b1, "B")
ax_b1.text(
    0.0, 1.17, "Split by tissue",
    transform=ax_b1.transAxes,
    fontsize=TITLE_FS, fontweight="bold",
    ha="left", va="bottom"
)


# C: pseudotime UMAP -----------------------------------------------------------
ax_c = fig.add_subplot(main_gs[1, 0])
sc_c = ax_c.scatter(
    umap[:, 0], umap[:, 1],
    c=adata.obs["dpt_pseudotime"], cmap="viridis",
    s=POINT_SIZE, alpha=0.75, rasterized=True, linewidths=0
)
add_panel_label(ax_c, "C")
ax_c.set_title("Pseudotime", fontsize=TITLE_FS, fontweight="bold", pad=4)
clean_umap_axis(ax_c)
cbar_c = plt.colorbar(sc_c, ax=ax_c, fraction=0.046, pad=0.02)
cbar_c.set_label("Pseudotime", fontsize=AXIS_FS)
cbar_c.ax.tick_params(labelsize=TICK_FS)


# D: pseudotime violin ---------------------------------------------------------
ax_d = fig.add_subplot(main_gs[1, 1])
sns.violinplot(
    data=plot_data, x="celltype", y="dpt_pseudotime",
    order=CELLTYPE_ORDER, palette=ALL_COLORS, ax=ax_d,
    inner="box", cut=0, linewidth=1.0, scale="width"
)
for i, ct in enumerate(CELLTYPE_ORDER):
    if ct in medians:
        ax_d.scatter(i, medians[ct], color="white", s=18, zorder=5,
                     edgecolors="black", lw=0.8)

ax_d.axvline(6.5, color="black", lw=1.3, ls="--", alpha=0.65)
add_panel_label(ax_d, "D")
ax_d.set_title("Pseudotime distribution", fontsize=TITLE_FS, fontweight="bold", pad=4)
ax_d.set_xlabel("")
ax_d.set_ylabel("Pseudotime", fontsize=AXIS_FS)
ax_d.tick_params(axis="y", labelsize=TICK_FS)
ax_d.set_ylim(-0.05, 1.08)
ax_d.spines[["top", "right"]].set_visible(False)

# Full names, rotated 90 degrees counter-clockwise, with extra bottom margin.
ax_d.set_xticklabels(CELLTYPE_ORDER, rotation=90, ha="center", va="top", fontsize=7.3)
ax_d.tick_params(axis="x", pad=2)


# E: PAGA connection map -------------------------------------------------------
ax_e = fig.add_subplot(main_gs[1, 2])

# Compact PAGA schematic that fits a single-column panel while keeping labels readable.
MAP_POS = {
    "Prolif. preNeu (S)":   (0.0,  0.0),
    "Prolif. preNeu (G2M)": (0.0,  1.45),
    "Immature Neu":         (0.0,  2.9),
    "Act. Neu":             (0.0,  4.35),
    "Mature Neu":           (-1.35, 5.7),
    "IFN Neu":              (-1.35, 7.1),
    "Inflam. Neu":          (-1.35, 8.5),
    "Mature_TAN":           (2.2,  5.7),
    "Il1bhi_TAN":           (3.9,  7.1),
    "Ccl3hi_TAN":           (3.9,  4.35),
    "Cxcr2hi_TAN":          (3.9,  2.9),
}
MAP_LABELS = {
    "Prolif. preNeu (S)":   "Prolif.\npreNeu (S)",
    "Prolif. preNeu (G2M)": "Prolif.\npreNeu (G2M)",
    "Immature Neu":         "Immature\nNeu",
    "Act. Neu":             "Act.\nNeu",
    "Mature Neu":           "Mature\nNeu",
    "IFN Neu":              "IFN\nNeu",
    "Inflam. Neu":          "Inflam.\nNeu",
    "Mature_TAN":           "Mature\nTAN",
    "Il1bhi_TAN":           "Il1bhi\nTAN",
    "Ccl3hi_TAN":           "Ccl3hi\nTAN",
    "Cxcr2hi_TAN":          "Cxcr2hi\nTAN",
}
MAP_LABEL_CONFIG = {
    "Prolif. preNeu (S)":   (-0.52,  0.0, "right"),
    "Prolif. preNeu (G2M)": (-0.52,  1.45, "right"),
    "Immature Neu":         (-0.52,  2.9, "right"),
    "Act. Neu":             (-0.52,  4.35, "right"),
    "Mature Neu":           (-1.85, 5.7, "right"),
    "IFN Neu":              (-1.85, 7.1, "right"),
    "Inflam. Neu":          (-1.85, 8.5, "right"),
    "Mature_TAN":           (2.72, 5.7, "left"),
    "Il1bhi_TAN":           (4.42, 7.1, "left"),
    "Ccl3hi_TAN":           (4.42, 4.35, "left"),
    "Cxcr2hi_TAN":          (4.42, 2.9, "left"),
}

threshold = 0.05
for i, g1 in enumerate(groups):
    for j, g2 in enumerate(groups):
        if j <= i:
            continue
        w = conn[i, j]
        if w > threshold and g1 in MAP_POS and g2 in MAP_POS:
            p1, p2 = MAP_POS[g1], MAP_POS[g2]
            color = "#111111" if w >= 0.5 else "#888888"
            ax_e.plot(
                [p1[0], p2[0]], [p1[1], p2[1]],
                color=color,
                lw=max(0.6, w * 5.3),
                alpha=min(0.85, 0.25 + w * 0.65),
                zorder=2,
                solid_capstyle="round"
            )

node_r_summary = 0.30
for g in groups:
    if g not in MAP_POS:
        continue
    p_node = MAP_POS[g]
    c_node = ALL_COLORS.get(g, "#AAAAAA")
    ec = "#CC0000" if g in TUMOR_TYPES else "#333333"
    lw = 1.8 if g in TUMOR_TYPES else 1.0
    ax_e.add_patch(
        plt.Circle(p_node, radius=node_r_summary, color=c_node, ec=ec, lw=lw, zorder=4)
    )

for g, (lx, ly, ha) in MAP_LABEL_CONFIG.items():
    ax_e.text(
        lx, ly, MAP_LABELS.get(g, g),
        ha=ha, va="center",
        fontsize=7.5, fontweight="bold",
        color=ALL_COLORS.get(g, "#000000"),
        bbox=dict(fc="white", alpha=0.88, ec="none", boxstyle="round,pad=0.10"),
        zorder=7
    )

ax_e.text(-1.55, 9.35, "SPLEEN", fontsize=9.3, fontweight="bold", color="#2166AC", ha="center")
ax_e.text(3.35, 9.35, "TUMOR", fontsize=9.3, fontweight="bold", color="#D73027", ha="center")
ax_e.axvline(1.15, color="#AAAAAA", lw=1.0, ls="--", alpha=0.5, zorder=1)
add_panel_label(ax_e, "E")
ax_e.set_title("PAGA connection map", fontsize=TITLE_FS, fontweight="bold", pad=4)
ax_e.set_xlim(-2.65, 5.05)
ax_e.set_ylim(-0.65, 9.8)
ax_e.set_aspect("equal")
ax_e.axis("off")


# Global cell-type legend ------------------------------------------------------
legend_handles = [
    mpatches.Patch(color=ALL_COLORS[c], label=c) for c in CELLTYPE_ORDER
]
fig.legend(
    handles=legend_handles,
    loc="lower left",
    bbox_to_anchor=(0.055, 0.025),
    ncol=4,
    fontsize=LEGEND_FS,
    frameon=False,
    columnspacing=1.05,
    handlelength=1.05,
    handletextpad=0.35,
    borderaxespad=0
)

plt.suptitle(
    "Combined spleen + tumor neutrophil trajectory analysis",
    fontsize=14.5, fontweight="bold", y=0.975
)

savefig("09_comprehensive_summary")

print("\n=== All figures generated ===")
print(f"Output directory: {OUTPUT_DIR}/")
for f in sorted(os.listdir(OUTPUT_DIR)):
    sz = os.path.getsize(os.path.join(OUTPUT_DIR, f)) / 1e3
    print(f"  {f:<45} {sz:.0f} KB")
