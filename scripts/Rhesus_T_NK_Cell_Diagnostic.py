#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Rhesus SHIV — Tier 1 T/NK Misclassification Diagnostic (READ-ONLY)
===================================================================
Purpose
-------
Quantify the suspected Tier 1 NK-vs-T misclassification BEFORE changing any
assignment logic. Reads the merged object, recomputes the Tier 1 lineage
scores using the SAME marker sets as the processing pipeline (Section 6 of the
Rhesus SHIV v2 pipeline), and measures how much of the current NK compartment
is actually CD3+ cytotoxic T-cell signal riding in on the shared cytotoxic
effector block (GZMB / GNLY / PRF1 / NKG7) rather than true NK identity
(KLRD1 / KLRC3 / NCAM1 / FCGR3A).

Questions answered (numbers we read together before touching the pipeline):
  1. What fraction of Tier-1 NK calls co-express ANY CD3 (CD3D/CD3E/CD3G)?
     -> the smoking gun.
  2. Do those cells win their NK label on the NK-EXCLUSIVE anchors, or only on
     the SHARED cytotoxic block? -> re-score the NK set as two sub-scores.
  3. How large is the NK-vs-T score margin for the CD3+ NK calls (near-tie the
     gate flips cleanly, vs decisive)?
  4. Where do NK calls sit relative to the Leiden clusters (scattered inside
     CD3+ T clusters = mislabel, vs their own island = real NK)?
  5. CD8 protein (ADT) distribution across NK calls split by CD3 status.
  6. How many 'Unassigned' cells are actually CD3+ (the T cells that fell
     through the argmax / SCORE_THRESHOLD)?
  7. Per-animal NK fraction BEFORE vs AFTER a CD3 gate, elite controller 40707
     flagged. This is the prerequisite for the "naturally NK-high = protected"
     hypothesis: if 40707 stays NK-high after cleaning, the signal is real; if
     it collapses, the story was an annotation artifact.

Design principles carried over: diagnostic before implementation, warn-and-
continue rather than hard-stop, verify mechanism from real data first.

READ-ONLY: nothing is written back to the input .h5ad. Only diagnostic CSVs
and figures are saved to OUT_DIR. No assignment logic is changed here.

There is NO CD3 tag in the 11-marker ADT panel (CD28, CD69, CX3CR1, CD1c,
CD11b, CD206, CD163, CD169, CD68, CD8, CD4), so the CD3 gate is RNA-only for
now; TCR-negativity (scirpy) folds in once VDJ is integrated. CD4/CD8 are
clean TotalSeq-C and are used here for the protein cross-check.

Input : shiv_host_merged.h5ad  (host object + recovered SHIV in obs; produced
        at the end of the STARsolo viral-read-recovery session)
Author: Jake Lehle / Kaushal Lab
Date  : July 2026
"""

# %% ============================================================
# CELL 1 — CONFIGURATION AND DIALS
#   Everything adjustable lives here. Fixed logic is below.
# ============================================================
import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from scipy import sparse
import warnings
warnings.filterwarnings("ignore")

# --- Paths ---
WORKING_DIR = '/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV'
ADATA_IN    = os.path.join(WORKING_DIR, 'shiv_starsolo', 'merged',
                           'shiv_host_merged.h5ad')
OUT_DIR     = os.path.join(WORKING_DIR, 'annotation_output', 'tnk_diagnostic')

# --- Dials ---
SCORE_THRESHOLD   = 0.1        # matches pipeline Tier 1 assignment threshold
CD3_MIN_COUNTS    = 1          # "CD3+" = at least this many CD3 counts (any CD3 gene)
MARGIN_NEARTIE    = 0.05       # |NK - T| below this = near-tie a gate flips cleanly
ELITE_CONTROLLER  = '40707'
RANDOM_STATE      = 0          # deterministic score_genes control sampling

# --- Column-name candidates (script auto-detects which is present) ---
LABEL_COL_CANDIDATES   = ['tier1_celltype']
CLUSTER_COL_CANDIDATES = ['clusters', 'leiden']
ANIMAL_COL             = 'animal_id'
COUNTS_LAYER           = 'counts'   # raw counts layer if present

# --- Marker sets ---
# Tier 1 sets copied verbatim from pipeline Section 6 so the recomputed calls
# reproduce the pipeline's calls (provenance: Rhesus SHIV v2 processing script).
TIER1_MARKERS = {
    'T cells':     ['CD3D', 'CD3E', 'CD3G'],
    'B cells':     ['CD79A', 'MS4A1', 'CD19', 'BANK1', 'PAX5'],
    'NK cells':    ['KLRD1', 'KLRC3', 'GNLY', 'PRF1', 'GZMB', 'NKG7'],
    'Monocytes':   ['S100A8', 'S100A9', 'VCAN', 'CD14', 'LYZ'],
    'Macrophages': ['CD68', 'CD163', 'MRC1', 'MERTK', 'APOC1'],
    'DCs':         ['FLT3', 'CLEC9A', 'CLEC10A', 'IL3RA', 'CLEC4C'],
    'Mast cells':  ['CPA3', 'KIT', 'HPGD', 'ENPP3'],
    'Erythrocytes': ['HBB', 'HBA1', 'HBA2'],
    'Platelets':   ['PPBP', 'PF4', 'GP9'],
}

# The proposed decomposition of the NK set. This is the whole hypothesis:
# the current NK label is driven by SHARED genes that activated CD8 T cells
# also express, not by NK-EXCLUSIVE identity.
CD3_GENES      = ['CD3D', 'CD3E', 'CD3G']
NK_EXCLUSIVE   = ['KLRD1', 'KLRC3', 'NCAM1', 'FCGR3A']   # NCAM1=CD56, FCGR3A=CD16
NK_SHARED      = ['GNLY', 'PRF1', 'GZMB', 'NKG7']        # shared with cytotoxic CD8 T
T_SCORE_GENES  = ['CD3D', 'CD3E', 'CD3G']

# --- Figure style (Jake's standard: 28-34 pt, hex colors, PDF+PNG @300) ---
plt.rcParams.update({
    'font.size': 28,
    'axes.titlesize': 34,
    'axes.labelsize': 30,
    'xtick.labelsize': 28,
    'ytick.labelsize': 28,
    'legend.fontsize': 28,
    'figure.dpi': 100,
})
CD3POS_COLOR = '#D32F2F'   # CD3+ (would flip to T)
CD3NEG_COLOR = '#1976D2'   # CD3- (candidate real NK)
BEFORE_COLOR = '#90A4AE'   # NK fraction before gate
AFTER_COLOR  = '#1976D2'   # NK fraction after gate
ELITE_COLOR  = '#F9A825'   # highlight for 40707

os.makedirs(OUT_DIR, exist_ok=True)


def savefig(fig, name):
    """Save every figure as both PDF and PNG at 300 DPI."""
    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(OUT_DIR, f'{name}.{ext}'),
                    bbox_inches='tight', dpi=300)
    print(f"    saved figure: {name}.pdf / .png")


# %% ============================================================
# CELL 2 — LOAD (READ-ONLY) AND PREFLIGHT INVENTORY
#   Report what is present; warn on anything missing; do not stop.
# ============================================================
print("=" * 80)
print("CELL 2 — LOAD + PREFLIGHT INVENTORY  (read-only)")
print("=" * 80)
print(f"  input: {ADATA_IN}")

adata = sc.read_h5ad(ADATA_IN)
print(f"  loaded: {adata.n_obs:,} cells x {adata.n_vars:,} genes")


def genes_present(genes):
    return [g for g in genes if g in adata.var_names]


def warn_missing(label, wanted, present):
    missing = [g for g in wanted if g not in present]
    if missing:
        print(f"  WARN  {label}: missing {missing} (continuing with {present})")
    else:
        print(f"  OK    {label}: all present {present}")


# Resolve which label / cluster / counts columns actually exist
label_col = next((c for c in LABEL_COL_CANDIDATES if c in adata.obs.columns), None)
cluster_col = next((c for c in CLUSTER_COL_CANDIDATES if c in adata.obs.columns), None)
has_counts = COUNTS_LAYER in adata.layers
has_animal = ANIMAL_COL in adata.obs.columns
has_umap = 'X_umap' in adata.obsm

print(f"\n  label column     : {label_col}")
print(f"  cluster column   : {cluster_col}")
print(f"  counts layer     : {COUNTS_LAYER if has_counts else 'ABSENT (will use X)'}")
print(f"  animal_id in obs : {has_animal}")
print(f"  UMAP in obsm     : {has_umap}")

# ADT inventory: accept CLR-normalized (protein_clr) OR raw (protein_counts).
# adata_qc.h5ad carries raw protein_counts; adata_pp/annotated carry protein_clr.
adt_names = list(adata.uns.get('adt_names', []))
adt_source = None
if 'protein_clr' in adata.obsm and len(adt_names) > 0:
    adt_source = 'protein_clr'      # already CLR-normalized
elif 'protein_counts' in adata.obsm and len(adt_names) > 0:
    adt_source = 'protein_counts'   # raw counts; CLR-normalize on the fly in Cell 6
has_adt = adt_source is not None
cd8_adt_name = None
if has_adt:
    cd8_adt_name = next((n for n in adt_names if 'CD8' in n.upper()), None)
    print(f"  ADT proteins     : {adt_names}")
    print(f"  ADT source       : {adt_source}")
    print(f"  CD8 ADT resolved : {cd8_adt_name}")
else:
    print("  WARN  no ADT (protein_clr / protein_counts + adt_names) found; "
          "CD8-protein cross-check will be skipped.")

# Recovered SHIV inventory (from the STARsolo merge)
for col in ('shiv_umi', 'shiv_pos'):
    print(f"  {col:9s} in obs : {col in adata.obs.columns}")

# Marker presence for the sets we depend on
print("\n  marker-set presence:")
warn_missing('CD3 genes', CD3_GENES, genes_present(CD3_GENES))
warn_missing('NK exclusive', NK_EXCLUSIVE, genes_present(NK_EXCLUSIVE))
warn_missing('NK shared', NK_SHARED, genes_present(NK_SHARED))
for lin, mk in TIER1_MARKERS.items():
    warn_missing(f'Tier1 {lin}', mk, genes_present(mk))


# %% ============================================================
# CELL 3 — PREP EXPRESSION + RECOMPUTE TIER 1 SCORES
#   Ensure X is log-normalized for scoring, recompute the pipeline's
#   Tier 1 lineage scores, and establish the "current" labels.
# ============================================================
print("\n" + "=" * 80)
print("CELL 3 — EXPRESSION PREP + TIER 1 SCORE RECOMPUTE")
print("=" * 80)

# Guard: make sure X is log-normalized before score_genes (in-memory only).
xmax = adata.X.max()
looks_raw = float(xmax) > 50 and np.allclose(
    adata.X[:50].toarray() if sparse.issparse(adata.X) else adata.X[:50],
    np.round(adata.X[:50].toarray() if sparse.issparse(adata.X) else adata.X[:50])
)
if looks_raw:
    print(f"  X looks like raw counts (max={xmax:.0f}); normalizing a working "
          f"copy in memory (source object untouched).")
    if not has_counts:
        # preserve raw counts for CD3 DETECTION before we overwrite X, so Cell 4
        # measures true detection (>=1 raw count) rather than a threshold on
        # log-normalized values.
        adata.layers[COUNTS_LAYER] = adata.X.copy()
        has_counts = True
        print(f"  stashed raw counts into layers['{COUNTS_LAYER}'] for detection")
    else:
        adata.X = adata.layers[COUNTS_LAYER].copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
else:
    print(f"  X appears log-normalized (max={xmax:.2f}); scoring on X as-is.")


def score_set(genes, key, min_genes=2):
    present = genes_present(genes)
    if len(present) >= min_genes:
        sc.tl.score_genes(adata, gene_list=present, score_name=key,
                          random_state=RANDOM_STATE)
    elif len(present) == 1:
        expr = adata[:, present[0]].X
        expr = expr.toarray().flatten() if sparse.issparse(expr) else np.asarray(expr).flatten()
        adata.obs[key] = expr
        print(f"  note {key}: single marker ({present[0]}), using raw expression")
    else:
        adata.obs[key] = 0.0
        print(f"  WARN {key}: no markers present, score set to 0")
    return present


# Recompute the pipeline's Tier 1 lineage scores
tier1_score_cols = {}
for lin, mk in TIER1_MARKERS.items():
    key = f"dx_tier1_{lin.replace(' ', '_')}"
    score_set(mk, key)
    tier1_score_cols[lin] = key

# Recomputed argmax label (reproduces the pipeline call)
lin_names = list(tier1_score_cols.keys())
score_mat = adata.obs[[tier1_score_cols[n] for n in lin_names]].values
best_idx = np.argmax(score_mat, axis=1)
best_score = score_mat[np.arange(len(best_idx)), best_idx]
adata.obs['dx_tier1_recomputed'] = np.where(
    best_score >= SCORE_THRESHOLD,
    [lin_names[i] for i in best_idx], 'Unassigned')

# Choose the working label: prefer the saved pipeline label if present.
if label_col is not None:
    adata.obs['dx_label'] = adata.obs[label_col].astype(str)
    print(f"\n  using saved label column '{label_col}' as the current call.")
    # Report agreement between saved and recomputed, as a sanity check.
    agree = (adata.obs['dx_label'] == adata.obs['dx_tier1_recomputed']).mean()
    print(f"  saved-vs-recomputed Tier1 agreement: {100*agree:.1f}%")
else:
    adata.obs['dx_label'] = adata.obs['dx_tier1_recomputed']
    print("\n  no saved label found; using the recomputed argmax as the "
          "current call (matches pipeline logic).")

print("\n  current Tier 1 distribution:")
print(adata.obs['dx_label'].value_counts().to_string())

# The two NK sub-scores + a clean T score for the margin analysis
score_set(NK_EXCLUSIVE, 'dx_nk_exclusive')
score_set(NK_SHARED,    'dx_nk_shared')
score_set(T_SCORE_GENES, 'dx_t_score')


# %% ============================================================
# CELL 4 — CD3 STATUS + THE SMOKING GUN
#   How many NK calls are CD3+? How many Unassigned are CD3+?
# ============================================================
print("\n" + "=" * 80)
print("CELL 4 — CD3 STATUS  (smoking gun)")
print("=" * 80)

# CD3+ = at least CD3_MIN_COUNTS of any CD3 gene, measured on RAW counts
cd3_present = genes_present(CD3_GENES)
if has_counts:
    cd3_mat = adata[:, cd3_present].layers[COUNTS_LAYER]
else:
    print("  WARN no counts layer; measuring CD3 detection on X instead.")
    cd3_mat = adata[:, cd3_present].X
cd3_mat = cd3_mat.toarray() if sparse.issparse(cd3_mat) else np.asarray(cd3_mat)
adata.obs['dx_cd3_counts'] = cd3_mat.sum(axis=1)
adata.obs['dx_cd3_pos'] = adata.obs['dx_cd3_counts'] >= CD3_MIN_COUNTS

is_nk = adata.obs['dx_label'] == 'NK cells'
is_un = adata.obs['dx_label'] == 'Unassigned'

n_nk = int(is_nk.sum())
n_nk_cd3 = int((is_nk & adata.obs['dx_cd3_pos']).sum())
print(f"  NK-called cells             : {n_nk:,}")
print(f"  ...of which CD3+ (>= {CD3_MIN_COUNTS} CD3 count): "
      f"{n_nk_cd3:,} ({100*n_nk_cd3/max(n_nk,1):.1f}%)  <-- likely mislabeled T")

n_un = int(is_un.sum())
n_un_cd3 = int((is_un & adata.obs['dx_cd3_pos']).sum())
print(f"\n  Unassigned cells            : {n_un:,}")
print(f"  ...of which CD3+             : "
      f"{n_un_cd3:,} ({100*n_un_cd3/max(n_un,1):.1f}%)  <-- T cells that fell through")

# For reference, CD3+ rate in the cells actually called T
is_t = adata.obs['dx_label'] == 'T cells'
if is_t.sum() > 0:
    print(f"\n  (reference) CD3+ rate among 'T cells' calls: "
          f"{100*adata.obs.loc[is_t, 'dx_cd3_pos'].mean():.1f}%")


# %% ============================================================
# CELL 5 — NK SCORE DECOMPOSITION + PREDICTED FLIP + MARGIN
#   Among NK calls: do they win on EXCLUSIVE identity or only SHARED?
#   Predict how many the CD3 gate flips, and how decisive the margin is.
# ============================================================
print("\n" + "=" * 80)
print("CELL 5 — NK DECOMPOSITION + PREDICTED FLIP")
print("=" * 80)

nk = adata.obs.loc[is_nk].copy()

# Does each NK cell score higher on shared cytotoxic than on exclusive NK?
nk['wins_on_shared'] = nk['dx_nk_shared'] > nk['dx_nk_exclusive']
nk['nk_minus_t'] = nk[tier1_score_cols['NK cells']] - nk['dx_t_score']
nk['near_tie'] = nk['nk_minus_t'].abs() < MARGIN_NEARTIE

n_shared = int(nk['wins_on_shared'].sum())
print(f"  NK calls winning higher on SHARED than EXCLUSIVE: "
      f"{n_shared:,}/{len(nk):,} ({100*n_shared/max(len(nk),1):.1f}%)")

# Cross the two axes that matter: CD3 status x what the NK label rests on
xt = pd.crosstab(nk['dx_cd3_pos'].map({True: 'CD3+', False: 'CD3-'}),
                 nk['wins_on_shared'].map({True: 'wins_on_SHARED',
                                           False: 'wins_on_EXCLUSIVE'}))
print("\n  NK calls: CD3 status x NK-label basis")
print(xt.to_string())
print("  (CD3+ & wins_on_SHARED = the classic mislabeled cytotoxic CD8 T cell)")

# Predicted correction under the proposed CD3 hard gate:
#   CD3+ NK  -> reassigned to T (the gate)
#   CD3- NK  -> stays candidate NK (confirmed later by exclusive score + TCR-neg)
nk['predicted'] = np.where(nk['dx_cd3_pos'], 'T cells (flipped)', 'NK (retained)')
print("\n  predicted effect of CD3 hard gate on the NK compartment:")
print(nk['predicted'].value_counts().to_string())

# Margin readout for the flipped set: near-tie flips clean; decisive is a
# stronger claim of a genuinely dual-identity cell worth a closer look.
flip = nk[nk['dx_cd3_pos']]
if len(flip) > 0:
    print(f"\n  CD3+ NK calls (to be flipped): {len(flip):,}")
    print(f"    near-tie (|NK-T| < {MARGIN_NEARTIE}): "
          f"{int(flip['near_tie'].sum()):,} ({100*flip['near_tie'].mean():.1f}%)")
    print(f"    median NK-minus-T margin        : {flip['nk_minus_t'].median():.3f}")

# Figure A: decomposition scatter (exclusive vs shared), colored by CD3 status
try:
    fig, ax = plt.subplots(figsize=(12, 11))
    for status, color in (('CD3-', CD3NEG_COLOR), ('CD3+', CD3POS_COLOR)):
        m = (nk['dx_cd3_pos'].map({True: 'CD3+', False: 'CD3-'}) == status)
        ax.scatter(nk.loc[m, 'dx_nk_exclusive'], nk.loc[m, 'dx_nk_shared'],
                   s=14, alpha=0.5, c=color, label=status, edgecolors='none')
    lims = [min(nk['dx_nk_exclusive'].min(), nk['dx_nk_shared'].min()),
            max(nk['dx_nk_exclusive'].max(), nk['dx_nk_shared'].max())]
    ax.plot(lims, lims, ls='--', lw=2, c='#000000', alpha=0.5)
    ax.set_xlabel('NK-exclusive score\n(KLRD1/KLRC3/NCAM1/FCGR3A)')
    ax.set_ylabel('shared cytotoxic score\n(GNLY/PRF1/GZMB/NKG7)')
    ax.set_title('NK calls: identity vs shared effectors')
    ax.legend(title=None, markerscale=2)
    ax.spines[['top', 'right']].set_visible(False)
    savefig(fig, 'A_nk_score_decomposition')
    plt.close(fig)
except Exception as e:
    print(f"  WARN figure A skipped: {e}")

# Figure B: NK-minus-T margin histogram, split by CD3 status
try:
    fig, ax = plt.subplots(figsize=(12, 9))
    bins = np.linspace(nk['nk_minus_t'].quantile(0.01),
                       nk['nk_minus_t'].quantile(0.99), 40)
    for status, color in (('CD3-', CD3NEG_COLOR), ('CD3+', CD3POS_COLOR)):
        m = (nk['dx_cd3_pos'].map({True: 'CD3+', False: 'CD3-'}) == status)
        ax.hist(nk.loc[m, 'nk_minus_t'], bins=bins, alpha=0.6,
                color=color, label=status)
    ax.axvline(0, ls='--', lw=2, c='#000000', alpha=0.6)
    ax.set_xlabel('NK score minus T score')
    ax.set_ylabel('NK-called cells')
    ax.set_title('How decisively did NK win over T?')
    ax.legend()
    ax.spines[['top', 'right']].set_visible(False)
    savefig(fig, 'B_nk_minus_t_margin')
    plt.close(fig)
except Exception as e:
    print(f"  WARN figure B skipped: {e}")


# %% ============================================================
# CELL 6 — CD8 PROTEIN (ADT) CROSS-CHECK ON NK CALLS
#   CD3+ NK calls should read CD8-protein-high (cytotoxic T); CD3- NK
#   calls should be lower. This is the defensible protein corroboration.
# ============================================================
print("\n" + "=" * 80)
print("CELL 6 — CD8 ADT CROSS-CHECK ON NK CALLS")
print("=" * 80)

if has_adt and cd8_adt_name is not None:
    if adt_source == 'protein_clr':
        adt_mat = np.asarray(adata.obsm['protein_clr'])
    else:
        # CLR-normalize raw protein_counts on the fly (matches pipeline Section 5)
        raw = np.asarray(adata.obsm['protein_counts'], dtype=float)
        logp = np.log1p(raw)
        gmean = np.exp(np.mean(logp, axis=1, keepdims=True))
        adt_mat = logp - np.log(gmean)
        print("  CLR-normalized raw protein_counts on the fly")
    clr = pd.DataFrame(adt_mat, index=adata.obs_names, columns=adt_names)
    adata.obs['dx_cd8_adt'] = clr[cd8_adt_name].values
    grp = adata.obs.loc[is_nk].groupby(
        adata.obs.loc[is_nk, 'dx_cd3_pos'].map({True: 'CD3+', False: 'CD3-'}))['dx_cd8_adt']
    print(f"  CD8 ADT ('{cd8_adt_name}', CLR) among NK calls:")
    print(grp.describe()[['count', 'mean', '50%']].to_string())

    try:
        fig, ax = plt.subplots(figsize=(10, 9))
        data, labels, colors = [], [], []
        for status, color in (('CD3-', CD3NEG_COLOR), ('CD3+', CD3POS_COLOR)):
            m = is_nk & (adata.obs['dx_cd3_pos'].map({True: 'CD3+', False: 'CD3-'}) == status)
            vals = adata.obs.loc[m, 'dx_cd8_adt'].dropna().values
            if len(vals) > 0:
                data.append(vals); labels.append(status); colors.append(color)
        parts = ax.violinplot(data, showmedians=True, showextrema=False)
        for pc, color in zip(parts['bodies'], colors):
            pc.set_facecolor(color); pc.set_alpha(0.6)
        ax.set_xticks(range(1, len(labels) + 1))
        ax.set_xticklabels(labels)
        ax.set_ylabel('CD8 protein (CLR)')
        ax.set_title('CD8 ADT across NK calls')
        ax.spines[['top', 'right']].set_visible(False)
        savefig(fig, 'C_cd8_adt_on_nk_calls')
        plt.close(fig)
    except Exception as e:
        print(f"  WARN figure C skipped: {e}")
else:
    print("  skipped: CD8 ADT not available in this object.")


# %% ============================================================
# CELL 7 — CLUSTER CROSS-TAB
#   Are NK calls scattered inside CD3+ T clusters (mislabel) or a
#   distinct CD3- island (real NK)?
# ============================================================
print("\n" + "=" * 80)
print("CELL 7 — NK CALLS vs LEIDEN CLUSTERS")
print("=" * 80)

if cluster_col is not None:
    per_cluster = adata.obs.groupby(cluster_col).agg(
        n_cells=('dx_label', 'size'),
        pct_cd3_pos=('dx_cd3_pos', lambda s: 100 * s.mean()),
        pct_nk_call=('dx_label', lambda s: 100 * (s == 'NK cells').mean()),
        pct_t_call=('dx_label', lambda s: 100 * (s == 'T cells').mean()),
    ).sort_values('pct_nk_call', ascending=False)
    print(f"  per-cluster summary (cluster col = '{cluster_col}'):")
    print(per_cluster.round(1).to_string())
    per_cluster.round(2).to_csv(os.path.join(OUT_DIR, 'cluster_cd3_nk_summary.csv'))
    print("    saved: cluster_cd3_nk_summary.csv")
    print("  read: clusters that are high pct_nk_call AND high pct_cd3_pos are "
          "mislabeled cytotoxic T; high pct_nk_call with low pct_cd3_pos are "
          "the real NK island.")
else:
    print("  skipped: no cluster column found in obs "
          f"({CLUSTER_COL_CANDIDATES}). Re-run neighbors+leiden upstream if "
          "you want this panel.")


# %% ============================================================
# CELL 8 — PER-ANIMAL NK FRACTION, BEFORE vs AFTER CD3 GATE
#   Prerequisite for the "naturally NK-high = protected" hypothesis.
#   Reported as % of all cells AND % of the (T + NK) compartment, which
#   is the cleaner denominator for "how much of the cytotoxic-lymphocyte
#   pool is truly NK vs mislabeled T". Elite controller 40707 flagged.
# ============================================================
print("\n" + "=" * 80)
print("CELL 8 — PER-ANIMAL NK FRACTION (before vs after CD3 gate)")
print("=" * 80)

if has_animal:
    df = adata.obs[[ANIMAL_COL, 'dx_label', 'dx_cd3_pos']].copy()
    df = df[df[ANIMAL_COL] != 'Unknown']

    # After the CD3 gate: CD3+ NK calls are reassigned to T (in-memory only).
    df['label_after'] = df['dx_label']
    flip_mask = (df['dx_label'] == 'NK cells') & (df['dx_cd3_pos'])
    df.loc[flip_mask, 'label_after'] = 'T cells'

    rows = []
    for animal, g in df.groupby(ANIMAL_COL):
        n_total = len(g)
        n_nk_before = (g['dx_label'] == 'NK cells').sum()
        n_nk_after = (g['label_after'] == 'NK cells').sum()
        # (T + NK) compartment denominators
        tnk_before = ((g['dx_label'] == 'NK cells') | (g['dx_label'] == 'T cells')).sum()
        tnk_after = ((g['label_after'] == 'NK cells') | (g['label_after'] == 'T cells')).sum()
        rows.append({
            'animal': animal,
            'elite_controller': (animal == ELITE_CONTROLLER),
            'n_cells': n_total,
            'NK_pct_allcells_before': 100 * n_nk_before / max(n_total, 1),
            'NK_pct_allcells_after': 100 * n_nk_after / max(n_total, 1),
            'NK_pct_of_TNK_before': 100 * n_nk_before / max(tnk_before, 1),
            'NK_pct_of_TNK_after': 100 * n_nk_after / max(tnk_after, 1),
        })
    per_animal = pd.DataFrame(rows).sort_values('NK_pct_of_TNK_after',
                                                ascending=False)
    per_animal['NK_pct_of_TNK_delta'] = (per_animal['NK_pct_of_TNK_after']
                                         - per_animal['NK_pct_of_TNK_before'])
    print(per_animal.round(2).to_string(index=False))
    per_animal.round(3).to_csv(os.path.join(OUT_DIR, 'per_animal_nk_fraction.csv'),
                               index=False)
    print("    saved: per_animal_nk_fraction.csv")

    # Figure D: per-animal NK % of (T+NK), before vs after, 40707 highlighted
    try:
        fig, ax = plt.subplots(figsize=(14, 9))
        x = np.arange(len(per_animal))
        w = 0.38
        ax.bar(x - w/2, per_animal['NK_pct_of_TNK_before'], width=w,
               color=BEFORE_COLOR, label='before gate')
        ax.bar(x + w/2, per_animal['NK_pct_of_TNK_after'], width=w,
               color=AFTER_COLOR, label='after CD3 gate')
        # highlight the elite controller tick
        labels = []
        for _, r in per_animal.iterrows():
            labels.append(f"{r['animal']}*" if r['elite_controller'] else r['animal'])
        for i, r in enumerate(per_animal.itertuples()):
            if r.elite_controller:
                ax.axvspan(i - 0.5, i + 0.5, color=ELITE_COLOR, alpha=0.15, zorder=0)
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=45, ha='right')
        ax.set_ylabel('NK % of (T + NK)')
        ax.set_title('Per-animal NK fraction\n(* = elite controller 40707)')
        ax.legend()
        ax.spines[['top', 'right']].set_visible(False)
        savefig(fig, 'D_per_animal_nk_fraction')
        plt.close(fig)
    except Exception as e:
        print(f"  WARN figure D skipped: {e}")
else:
    print(f"  skipped: '{ANIMAL_COL}' not in obs.")


# %% ============================================================
# CELL 9 — HEADLINE SUMMARY
#   The few numbers we read together before writing the fix.
# ============================================================
print("\n" + "=" * 80)
print("CELL 9 — HEADLINE SUMMARY")
print("=" * 80)
print(f"  input object            : {os.path.basename(ADATA_IN)}  ({adata.n_obs:,} cells)")
print(f"  NK calls                : {n_nk:,}")
print(f"  NK calls that are CD3+   : {n_nk_cd3:,} "
      f"({100*n_nk_cd3/max(n_nk,1):.1f}%)  -> flip to T under the gate")
print(f"  NK retained (CD3-)       : {n_nk - n_nk_cd3:,}")
print(f"  NK winning only on SHARED: {n_shared:,} "
      f"({100*n_shared/max(n_nk,1):.1f}%)")
print(f"  Unassigned that are CD3+ : {n_un_cd3:,} "
      f"({100*n_un_cd3/max(n_un,1):.1f}%)  -> recovered as T by the gate")
print(f"\n  outputs written to      : {OUT_DIR}")
print("  NOTE: source .h5ad NOT modified. No assignment logic changed. This is")
print("        a read-only diagnostic; the fix comes after we read these numbers.")
print("=" * 80)
