#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Rhesus SHIV — Embedding on the Final Labeled Object
====================================================
Computes PCA / neighbors / UMAP / Leiden and CLR-normalizes the ADT on the
VDJ-integrated object, so cell types, surface markers, TCR status, and the
recovered SHIV all live on one set of coordinates for plotting.

Reproduces the pipeline's Section 5 recipe so this embedding is consistent with
how the rest of the project was built:
  - HVG: seurat_v3 on the counts layer, exclude SHIV + mito genes, 3000 top
  - PCA 50 comps, neighbors 30 PCs, Leiden res 1.0 (random_state 42)
  - BBKNN batch correction on 'library' (GEM-well level = technical variation)
  - ADT CLR -> obsm['protein_clr']

Runs a before/after BBKNN pass for the batch-mixing QC. Robust to whether X is
already log-normalized (it is, from the correction step) or still raw.

Output: a new embedded object; the plotting script reads it. Labels are NOT
changed here.

Author: Jake Lehle / Kaushal Lab
Date: July 2026
"""

# %% ============================================================
# CELL 1 — CONFIG AND DIALS
# ============================================================
import os
os.environ.setdefault('MPLBACKEND', 'Agg')   # headless: save figures, no GUI/Qt
import matplotlib
matplotlib.use('Agg')
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from scipy import sparse
import warnings
warnings.filterwarnings("ignore")

WORKING_DIR = '/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV'
ADATA_IN  = os.path.join(WORKING_DIR, 'shiv_starsolo', 'merged',
                         'shiv_host_tier1_tcr.h5ad')
ADATA_OUT = os.path.join(WORKING_DIR, 'shiv_starsolo', 'merged',
                         'shiv_host_final_embedded.h5ad')
OUT_DIR   = os.path.join(WORKING_DIR, 'annotation_output', 'embedding')

# --- Dials ---
HVG_FLAVOR       = 'seurat'   # scanpy default (dispersion on log-norm X);
                              # set 'seurat_v3' only to reproduce the old pipeline
N_TOP_HVG        = 3000
N_PCS_PCA        = 50
N_PCS_NEIGHBORS  = 30
LEIDEN_RES       = 1.0
RANDOM_STATE     = 42
BATCH_KEY        = 'library'
COUNTS_LAYER     = 'counts'
DO_PREBBKNN_QC   = True     # compute the pre-correction UMAP for batch QC

plt.rcParams.update({
    'font.size': 28, 'axes.titlesize': 34, 'axes.labelsize': 30,
    'xtick.labelsize': 28, 'ytick.labelsize': 28, 'legend.fontsize': 18,
    'figure.dpi': 100,
})
os.makedirs(OUT_DIR, exist_ok=True)
sc.settings.figdir = OUT_DIR


def savefig(fig, name):
    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(OUT_DIR, f'{name}.{ext}'),
                    bbox_inches='tight', dpi=300)
    print(f"    saved figure: {name}.pdf / .png")


# %% ============================================================
# CELL 2 — LOAD + ENSURE counts LAYER AND LOG-NORM X
# ============================================================
print("=" * 80)
print("CELL 2 — LOAD + EXPRESSION STATE")
print("=" * 80)
adata = sc.read_h5ad(ADATA_IN)
print(f"  {ADATA_IN}")
print(f"  {adata.n_obs:,} cells x {adata.n_vars:,} genes")

xmax = adata.X.max()
head = adata.X[:50].toarray() if sparse.issparse(adata.X) else adata.X[:50]
x_is_raw = float(xmax) > 50 and np.allclose(head, np.round(head))

if COUNTS_LAYER in adata.layers:
    print(f"  counts layer present (raw counts for seurat_v3 HVG)")
    if x_is_raw:
        print("  X looks raw; normalizing X for PCA.")
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    else:
        print("  X already log-normalized; leaving as-is.")
else:
    if x_is_raw:
        print("  no counts layer; X is raw -> stashing counts, normalizing X.")
        adata.layers[COUNTS_LAYER] = adata.X.copy()
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    else:
        print("  WARN no counts layer and X is log-normalized; seurat_v3 needs")
        print("       raw counts, will fall back to flavor='seurat' on X.")


# %% ============================================================
# CELL 3 — HVG (exclude SHIV + mito) + PCA
# ============================================================
print("\n" + "=" * 80)
print("CELL 3 — HVG + PCA")
print("=" * 80)

exclude = pd.Series(False, index=adata.var_names)
for flag in ('shiv', 'mt'):
    if flag in adata.var.columns:
        exclude |= adata.var[flag].astype(bool).values
print(f"  genes excluded from HVG (SHIV/mito): {int(exclude.sum())}")

sub = adata[:, ~exclude.values].copy()
if HVG_FLAVOR == 'seurat_v3':
    if COUNTS_LAYER not in adata.layers:
        raise RuntimeError("seurat_v3 needs raw counts in layers['counts'].")
    sc.pp.highly_variable_genes(sub, n_top_genes=N_TOP_HVG,
                                flavor='seurat_v3', layer=COUNTS_LAYER)
else:
    # scanpy default: dispersion-based on the log-normalized X (no counts needed)
    sc.pp.highly_variable_genes(sub, n_top_genes=N_TOP_HVG, flavor=HVG_FLAVOR)
hvg = sub.var_names[sub.var['highly_variable']].tolist()
del sub

adata.var['highly_variable'] = adata.var_names.isin(hvg)
print(f"  HVGs selected: {adata.var['highly_variable'].sum()} (flavor={HVG_FLAVOR})")

sc.tl.pca(adata, n_comps=N_PCS_PCA, use_highly_variable=True)
print(f"  PCA: {N_PCS_PCA} comps")


# %% ============================================================
# CELL 4 — PRE-BBKNN UMAP (batch QC)
# ============================================================
if DO_PREBBKNN_QC:
    print("\n" + "=" * 80)
    print("CELL 4 — PRE-BBKNN UMAP (batch QC)")
    print("=" * 80)
    sc.pp.neighbors(adata, n_pcs=N_PCS_NEIGHBORS)
    sc.tl.umap(adata)
    adata.obsm['X_umap_prebbknn'] = adata.obsm['X_umap'].copy()

    fig, axes = plt.subplots(1, 4, figsize=(30, 6))
    for ax, key, ttl in zip(
            axes, ['library', 'tissue', 'timepoint', 'tier1_celltype'],
            ['Library', 'Tissue', 'Timepoint', 'Tier 1']):
        sc.pl.umap(adata, color=key, ax=ax, show=False, title=ttl,
                   frameon=False, size=3)
    plt.tight_layout()
    savefig(fig, 'umap_before_bbknn')
    plt.close('all')


# %% ============================================================
# CELL 5 — BBKNN + UMAP + LEIDEN
# ============================================================
print("\n" + "=" * 80)
print("CELL 5 — BBKNN + UMAP + LEIDEN")
print("=" * 80)
try:
    import bbknn
    print(f"  BBKNN on batch_key='{BATCH_KEY}' "
          f"({adata.obs[BATCH_KEY].nunique()} batches)")
    bbknn.bbknn(adata, batch_key=BATCH_KEY, n_pcs=N_PCS_NEIGHBORS)
except ImportError:
    print("  WARN bbknn not installed; using standard neighbors (no correction)")
    sc.pp.neighbors(adata, n_pcs=N_PCS_NEIGHBORS)

sc.tl.umap(adata)
sc.tl.leiden(adata, key_added='clusters', resolution=LEIDEN_RES,
             random_state=RANDOM_STATE)
print(f"  Leiden clusters: {adata.obs['clusters'].nunique()}")

fig, axes = plt.subplots(1, 4, figsize=(30, 6))
for ax, key, ttl in zip(
        axes, ['library', 'clusters', 'tissue', 'tier1_celltype'],
        ['Library', 'Clusters', 'Tissue', 'Tier 1']):
    sc.pl.umap(adata, color=key, ax=ax, show=False, title=ttl,
               frameon=False, size=3)
plt.tight_layout()
savefig(fig, 'umap_after_bbknn')
plt.close('all')


# %% ============================================================
# CELL 6 — CLR-NORMALIZE ADT -> protein_clr
# ============================================================
print("\n" + "=" * 80)
print("CELL 6 — ADT CLR")
print("=" * 80)
if 'protein_counts' in adata.obsm and 'adt_names' in adata.uns:
    raw = np.asarray(adata.obsm['protein_counts'], dtype=float)
    logp = np.log1p(raw)
    gmean = np.exp(np.mean(logp, axis=1, keepdims=True))
    adata.obsm['protein_clr'] = logp - np.log(gmean)
    print(f"  protein_clr: {adata.obsm['protein_clr'].shape} "
          f"({len(adata.uns['adt_names'])} tags)")
    print("  (CD163/CD169 remain near-zero pending BD Ab-seq barcode recovery)")
else:
    print("  WARN no protein_counts / adt_names; skipping ADT CLR")


# %% ============================================================
# CELL 7 — CLUSTER x TIER1 CROSS-TAB (finally available)
# ============================================================
print("\n" + "=" * 80)
print("CELL 7 — CLUSTER x TIER1")
print("=" * 80)
ct = pd.crosstab(adata.obs['clusters'], adata.obs['tier1_celltype'])
ct.to_csv(os.path.join(OUT_DIR, 'cluster_by_tier1.csv'))
print(f"  {adata.obs['clusters'].nunique()} clusters x "
      f"{adata.obs['tier1_celltype'].nunique()} labels; "
      f"saved cluster_by_tier1.csv")
# purity: fraction of each cluster in its dominant label
purity = (ct.max(axis=1) / ct.sum(axis=1)).round(2)
print("  cluster purity (dominant-label fraction) summary:")
print(f"    median {purity.median():.2f}, min {purity.min():.2f}, "
      f"clusters < 0.6 purity: {int((purity < 0.6).sum())}")


# %% ============================================================
# CELL 8 — WRITE
# ============================================================
print("\n" + "=" * 80)
print("CELL 8 — WRITE")
print("=" * 80)
adata.write_h5ad(ADATA_OUT)
print(f"  [SAVED] {ADATA_OUT}")
print(f"  obsm: {list(adata.obsm.keys())}")
print(f"  clusters: {adata.obs['clusters'].nunique()}, "
      f"HVGs: {int(adata.var['highly_variable'].sum())}")
print("  Ready for the plotting script (cell types + surface receptors on UMAP).")
print("=" * 80)
