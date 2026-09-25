#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Rhesus SHIV — Tier 1 T/NK Annotation Correction
=================================================
Turns the diagnostic-confirmed mechanism into corrected Tier 1 labels.

What the diagnostic established (on shiv_host_merged.h5ad, 89,844 QC'd cells):
  - 62.2% of Tier-1 NK calls are CD3+ (raw-count detection); the flips are
    decisive (median NK-minus-T margin 0.762), i.e. the shared cytotoxic block
    (GNLY/PRF1/GZMB/NKG7) was overpowering the 3-gene CD3 score in argmax.
  - 39.1% of Unassigned cells are CD3+ T cells that fell through.
  - CD8 ADT does NOT separate CD3+ from CD3- NK calls (both CD8-protein-high),
    because the anti-CD8a tag picks up CD8aa on rhesus NK and CD3-dropout T
    cells alike. So CD3-negative retained NK stay PROVISIONAL pending TCR.

The correction (root-cause, not a patch):
  1. Rebuild the NK Tier-1 score on EXCLUSIVE anchors only (KLRD1/KLRC3/NCAM1,
     andNCR1 if present). The shared cytotoxic block is pulled
     OUT of the discriminator (kept as a separate cytotoxic score for later
     functional work, never for lineage assignment).
  2. CD3 hard gate: a CD3+ cell can NEVER be NK. Any cell that argmaxes to NK
     but is CD3+ is reassigned to T.
  3. Recover CD3+ Unassigned cells to T (the fell-through T cells).
  4. Flag the CD3-negative cells that remain NK as PROVISIONAL (obs
     'nk_provisional'), pending TCR-negativity from the VDJ integration.

Nothing is blessed blindly: BOTH the old-logic and new-logic labels are
computed here, the full transition matrix is printed and saved, the original
labels are preserved in obs, and output goes to a NEW file (the merge is not
overwritten). Read the transition report before treating tier1_celltype as
canonical.

Deferred to Tier 2 (once VDJ is integrated): a called TCR on any of these cells
= T (reassign); TCR-negative + CD3-negative + exclusive-NK = confirmed NK. This
script builds the hook (nk_provisional) but does not make that call.

Run WITH python (sc_pre env). Spyder cells (# %%).

Author: Jake Lehle 
Date: July 2026
"""

# %% ============================================================
# CELL 1 — CONFIGURATION AND DIALS
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
ADATA_OUT   = os.path.join(WORKING_DIR, 'shiv_starsolo', 'merged',
                           'shiv_host_tier1_corrected.h5ad')
OUT_DIR     = os.path.join(WORKING_DIR, 'annotation_output', 'tier1_correction')

# --- Dials ---
SCORE_THRESHOLD  = 0.1        # matches pipeline Tier 1 threshold
CD3_MIN_COUNTS   = 1          # CD3+ = at least this many CD3 counts (raw)
RANDOM_STATE     = 0
ELITE_CONTROLLER = '40707'
ANIMAL_COL       = 'animal_id'
COUNTS_LAYER     = 'counts'

# --- Marker sets ---
# Non-NK lineages: copied verbatim from pipeline Section 6.
TIER1_OTHER = {
    'T cells':      ['CD3D', 'CD3E', 'CD3G'],
    'B cells':      ['CD79A', 'MS4A1', 'CD19', 'BANK1', 'PAX5'],
    'Monocytes':    ['S100A8', 'S100A9', 'VCAN', 'CD14', 'LYZ'],
    'Macrophages':  ['CD68', 'CD163', 'MRC1', 'MERTK', 'APOC1'],
    'DCs':          ['FLT3', 'CLEC9A', 'CLEC10A', 'IL3RA', 'CLEC4C'],
    'Mast cells':   ['CPA3', 'KIT', 'HPGD', 'ENPP3'],
    'Erythrocytes': ['HBB', 'HBA1', 'HBA2'],
    'Platelets':    ['PPBP', 'PF4', 'GP9'],
}

# NK, OLD (6-gene) set: reproduces the pre-correction call for the transition.
NK_OLD = ['KLRD1', 'KLRC3', 'GNLY', 'PRF1', 'GZMB', 'NKG7']

# NK, EXCLUSIVE candidates: assignment uses whichever are present. NCR1 (NKp46)
# and KLRF1 (NKp80) are strong NK-exclusive markers
NK_EXCLUSIVE_CANDIDATES = ['KLRD1', 'KLRC3', 'NCAM1', 
                           'NCR1', 'KLRF1']

# Shared cytotoxic block: kept as an informational score, NOT for assignment.
NK_SHARED = ['GNLY', 'PRF1', 'GZMB', 'NKG7']

CD3_GENES = ['CD3D', 'CD3E', 'CD3G']

# --- Figure style (28-34 pt, hex colors, PDF+PNG @300) ---
plt.rcParams.update({
    'font.size': 28, 'axes.titlesize': 34, 'axes.labelsize': 30,
    'xtick.labelsize': 28, 'ytick.labelsize': 28, 'legend.fontsize': 28,
    'figure.dpi': 100,
})
COL_T       = '#D32F2F'   # T cells
COL_NK      = '#1976D2'   # NK cells
COL_OTHER   = '#90A4AE'   # other / unassigned
BEFORE_COL  = '#90A4AE'
AFTER_COL   = '#1976D2'
ELITE_COL   = '#F9A825'

os.makedirs(OUT_DIR, exist_ok=True)


def savefig(fig, name):
    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(OUT_DIR, f'{name}.{ext}'),
                    bbox_inches='tight', dpi=300)
    print(f"    saved figure: {name}.pdf / .png")


# %% ============================================================
# CELL 2 — LOAD, PRESERVE RAW COUNTS, NORMALIZE FOR SCORING
# ============================================================
print("=" * 80)
print("CELL 2 — LOAD + PREP")
print("=" * 80)
print(f"  input: {ADATA_IN}")

adata = sc.read_h5ad(ADATA_IN)
print(f"  loaded: {adata.n_obs:,} cells x {adata.n_vars:,} genes")


def genes_present(genes):
    return [g for g in genes if g in adata.var_names]


has_counts = COUNTS_LAYER in adata.layers
xmax = adata.X.max()
head = adata.X[:50].toarray() if sparse.issparse(adata.X) else adata.X[:50]
looks_raw = float(xmax) > 50 and np.allclose(head, np.round(head))

if looks_raw:
    print(f"  X looks like raw counts (max={xmax:.0f}); preserving counts and "
          f"normalizing a working copy in memory.")
    if not has_counts:
        adata.layers[COUNTS_LAYER] = adata.X.copy()
        has_counts = True
    else:
        adata.X = adata.layers[COUNTS_LAYER].copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
else:
    print(f"  X appears log-normalized (max={xmax:.2f}); scoring on X as-is.")

# Resolve the exclusive NK set actually available
nk_excl = genes_present(NK_EXCLUSIVE_CANDIDATES)
print(f"\n  NK exclusive markers found: {nk_excl}")
if len(nk_excl) < 2:
    print("  WARN  fewer than 2 NK-exclusive markers present; the exclusive "
          "rescore will be weak. Check for NCR1.")
print(f"  NK shared (informational only): {genes_present(NK_SHARED)}")
print(f"  CD3 genes: {genes_present(CD3_GENES)}")


# %% ============================================================
# CELL 3 — SCORE ALL LINEAGES (NK gets OLD + EXCLUSIVE scores)
# ============================================================
print("\n" + "=" * 80)
print("CELL 3 — LINEAGE SCORING")
print("=" * 80)


def score_set(genes, key, min_genes=2):
    present = genes_present(genes)
    if len(present) >= min_genes:
        sc.tl.score_genes(adata, gene_list=present, score_name=key,
                          random_state=RANDOM_STATE)
    elif len(present) == 1:
        expr = adata[:, present[0]].X
        expr = expr.toarray().flatten() if sparse.issparse(expr) else np.asarray(expr).flatten()
        adata.obs[key] = expr
    else:
        adata.obs[key] = 0.0
        print(f"  WARN {key}: no markers present")
    return present


other_cols = {}
for lin, mk in TIER1_OTHER.items():
    key = f"dx_{lin.replace(' ', '_')}"
    score_set(mk, key)
    other_cols[lin] = key

score_set(NK_OLD, 'dx_nk_old')            # 6-gene, reproduces the old call
score_set(nk_excl, 'dx_nk_excl')          # exclusive anchors, the corrected call
score_set(NK_SHARED, 'dx_cyto_shared')    # informational only
print("  scored: non-NK lineages, dx_nk_old, dx_nk_excl, dx_cyto_shared")


# %% ============================================================
# CELL 4 — CD3 STATUS (raw counts)
# ============================================================
print("\n" + "=" * 80)
print("CELL 4 — CD3 STATUS")
print("=" * 80)

cd3_present = genes_present(CD3_GENES)
src = COUNTS_LAYER if has_counts else None
cd3_mat = (adata[:, cd3_present].layers[COUNTS_LAYER] if src
           else adata[:, cd3_present].X)
cd3_mat = cd3_mat.toarray() if sparse.issparse(cd3_mat) else np.asarray(cd3_mat)
adata.obs['dx_cd3_pos'] = cd3_mat.sum(axis=1) >= CD3_MIN_COUNTS
print(f"  CD3+ cells (>= {CD3_MIN_COUNTS} count on "
      f"{'raw counts' if src else 'X'}): {int(adata.obs['dx_cd3_pos'].sum()):,}")


# %% ============================================================
# CELL 5 — OLD-LOGIC AND NEW-LOGIC LABELS
# ============================================================
print("\n" + "=" * 80)
print("CELL 5 — ASSIGN (old logic vs corrected logic)")
print("=" * 80)


def argmax_label(score_cols_dict):
    names = list(score_cols_dict.keys())
    mat = adata.obs[[score_cols_dict[n] for n in names]].values
    idx = np.argmax(mat, axis=1)
    best = mat[np.arange(len(idx)), idx]
    lab = np.where(best >= SCORE_THRESHOLD,
                   [names[i] for i in idx], 'Unassigned')
    return pd.Series(lab, index=adata.obs_names)


# OLD logic: NK scored on the 6-gene set, no gate (reproduces the diagnostic).
old_cols = dict(other_cols); old_cols['NK cells'] = 'dx_nk_old'
label_old = argmax_label(old_cols)

# NEW logic: NK scored on EXCLUSIVE anchors only.
new_cols = dict(other_cols); new_cols['NK cells'] = 'dx_nk_excl'
label_new = argmax_label(new_cols)

# CD3 hard gate: a CD3+ cell can never be NK -> reassign to T.
gate_nk = (label_new == 'NK cells') & adata.obs['dx_cd3_pos']
label_new[gate_nk] = 'T cells'

# Recover CD3+ Unassigned -> T (the fell-through T cells).
gate_un = (label_new == 'Unassigned') & adata.obs['dx_cd3_pos']
label_new[gate_un] = 'T cells'

adata.obs['tier1_celltype_uncorrected'] = pd.Categorical(label_old)
adata.obs['tier1_celltype'] = pd.Categorical(label_new)

# Provisional NK: retained NK are CD3-negative and rest on the exclusive score;
# CD8 ADT cannot separate CD8aa NK from CD3-dropout T, so mark them pending TCR.
adata.obs['nk_provisional'] = (label_new == 'NK cells')

print(f"  CD3 gate reassigned NK -> T        : {int(gate_nk.sum()):,}")
print(f"  CD3+ Unassigned recovered -> T     : {int(gate_un.sum()):,}")
print(f"  NK retained (CD3-, provisional)    : {int(adata.obs['nk_provisional'].sum()):,}")

print("\n  OLD-logic Tier 1 distribution:")
print(label_old.value_counts().to_string())
print("\n  CORRECTED Tier 1 distribution:")
print(label_new.value_counts().to_string())


# %% ============================================================
# CELL 6 — TRANSITION REPORT (read this before blessing anything)
# ============================================================
print("\n" + "=" * 80)
print("CELL 6 — OLD -> CORRECTED TRANSITION")
print("=" * 80)

trans = pd.crosstab(label_old, label_new,
                    rownames=['old'], colnames=['corrected'])
print(trans.to_string())
trans.to_csv(os.path.join(OUT_DIR, 'tier1_transition_matrix.csv'))
print("    saved: tier1_transition_matrix.csv")

# Focused view: where did the OLD NK calls go?
old_nk_fate = label_new[label_old == 'NK cells'].value_counts()
print("\n  fate of the OLD NK calls:")
print(old_nk_fate.to_string())

try:
    fig, ax = plt.subplots(figsize=(12, 9))
    fates = old_nk_fate.reindex(
        ['T cells', 'NK cells'] +
        [x for x in old_nk_fate.index if x not in ('T cells', 'NK cells')]
    ).dropna()
    colors = [COL_T if f == 'T cells' else COL_NK if f == 'NK cells'
              else COL_OTHER for f in fates.index]
    ax.bar(range(len(fates)), fates.values, color=colors)
    ax.set_xticks(range(len(fates)))
    ax.set_xticklabels(fates.index, rotation=45, ha='right')
    ax.set_ylabel('cells')
    ax.set_title('Where the old NK calls went')
    ax.spines[['top', 'right']].set_visible(False)
    savefig(fig, 'E_old_nk_fate')
    plt.close(fig)
except Exception as e:
    print(f"  WARN figure E skipped: {e}")


# %% ============================================================
# CELL 7 — PER-ANIMAL NK FRACTION (corrected labels), 40707 flagged
# ============================================================
print("\n" + "=" * 80)
print("CELL 7 — PER-ANIMAL NK FRACTION (old vs corrected)")
print("=" * 80)

if ANIMAL_COL in adata.obs.columns:
    df = pd.DataFrame({
        'animal': adata.obs[ANIMAL_COL].astype(str).values,
        'old': label_old.values,
        'new': label_new.values,
    })
    df = df[df['animal'] != 'Unknown']

    rows = []
    for animal, g in df.groupby('animal'):
        nk_o = (g['old'] == 'NK cells').sum()
        nk_n = (g['new'] == 'NK cells').sum()
        tnk_o = ((g['old'] == 'NK cells') | (g['old'] == 'T cells')).sum()
        tnk_n = ((g['new'] == 'NK cells') | (g['new'] == 'T cells')).sum()
        rows.append({
            'animal': animal,
            'elite_controller': animal == ELITE_CONTROLLER,
            'n_cells': len(g),
            'NK_pct_of_TNK_old': 100 * nk_o / max(tnk_o, 1),
            'NK_pct_of_TNK_corrected': 100 * nk_n / max(tnk_n, 1),
        })
    per_animal = pd.DataFrame(rows).sort_values('NK_pct_of_TNK_corrected',
                                                ascending=False)
    per_animal['delta'] = (per_animal['NK_pct_of_TNK_corrected']
                           - per_animal['NK_pct_of_TNK_old'])
    print(per_animal.round(2).to_string(index=False))
    per_animal.round(3).to_csv(
        os.path.join(OUT_DIR, 'per_animal_nk_fraction_corrected.csv'), index=False)
    print("    saved: per_animal_nk_fraction_corrected.csv")
    print("  NOTE: corrected NK is CD3-negative + exclusive-scored but still")
    print("        PROVISIONAL (CD8 ADT cannot split CD8aa NK from dropout T);")
    print("        treat these as an upper bound until TCR-negativity is applied.")

    try:
        fig, ax = plt.subplots(figsize=(14, 9))
        x = np.arange(len(per_animal))
        w = 0.38
        ax.bar(x - w/2, per_animal['NK_pct_of_TNK_old'], width=w,
               color=BEFORE_COL, label='old logic')
        ax.bar(x + w/2, per_animal['NK_pct_of_TNK_corrected'], width=w,
               color=AFTER_COL, label='corrected')
        labels = [f"{r.animal}*" if r.elite_controller else r.animal
                  for r in per_animal.itertuples()]
        for i, r in enumerate(per_animal.itertuples()):
            if r.elite_controller:
                ax.axvspan(i - 0.5, i + 0.5, color=ELITE_COL, alpha=0.15, zorder=0)
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=45, ha='right')
        ax.set_ylabel('NK % of (T + NK)')
        ax.set_title('Per-animal NK fraction\n(* = elite controller 40707)')
        ax.legend()
        ax.spines[['top', 'right']].set_visible(False)
        savefig(fig, 'F_per_animal_nk_corrected')
        plt.close(fig)
    except Exception as e:
        print(f"  WARN figure F skipped: {e}")
else:
    print(f"  skipped: '{ANIMAL_COL}' not in obs.")


# %% ============================================================
# CELL 8 — WRITE CORRECTED OBJECT (new file; merge not overwritten)
# ============================================================
print("\n" + "=" * 80)
print("CELL 8 — WRITE")
print("=" * 80)

adata.write_h5ad(ADATA_OUT)
print(f"  [SAVED] {ADATA_OUT}")
print(f"  new obs columns: tier1_celltype (corrected, canonical), "
      f"tier1_celltype_uncorrected, nk_provisional, dx_* scores")
print(f"  original merge object left untouched at {os.path.basename(ADATA_IN)}")
print("\n  Review the transition matrix (Cell 6) and the per-animal table")
print("  (Cell 7) before treating tier1_celltype as final. TCR-negativity")
print("  from the VDJ integration is the next step to resolve the provisional")
print("  NK cells (Tier 2).")
