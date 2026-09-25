#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Rhesus SHIV - Tier 1 TCR-based T Cell Recovery (step 1 of the recover-first pass)
=================================================================================
Recovers CD3-dropout T cells that fell into non-T Tier-1 labels but carry a
called TCR, so the Tier-2 recluster runs on the COMPLETE T pool. Also recomputes
the strict gamma-delta flag on the recovered pool (CELL 6), because the old flag
was set at VDJ-integration time and never saw the recovered cells.

What the Tier-2 diagnostic established (read-only, on the canonical object):
  - Unassigned pool = 9,832 cells; 0% CD3+ by RNA but 3,493 (35.5%) carry a
    called TCR. These are CD3-dropout T cells Tier 1 missed because it annotated
    on CD3D/CD3E/CD3G RNA, which drops out hardest in effector / cytotoxic T.
  - Additional non-Unassigned flips (Platelet / Erythrocyte / Mast) are cells
    that dropped CD3 and picked up ambient HBB/PPBP; ~56% carry a clonotype-
    definable complete pair, so the majority are real dropout-T, not soup.

The recovery rule (one rule, applied once):
  flip -> 'T cells'  IF  has_tcr == True
                     AND current tier1_celltype is not 'T cells'
                     AND current tier1_celltype is not in the B/myeloid
                         quarantine set {B cells, Monocytes, Macrophages, DCs}

Rationale for the quarantine: B and myeloid lineages carry documented ambient
orphan-TCR rates (B 9.3%, Monocytes 10.7%, Macrophages 14.4%), so a lone TCR
chain there is more likely soup than a real receptor. For Unassigned and for
Erythrocyte / Platelet (whose HBB/PPBP is itself ambient soup), a detected
receptor is the lineage evidence, so those flip. NK-with-TCR is definitionally
T (already handled upstream during VDJ integration) and flips defensively.

Nothing is blessed blindly (same discipline as Rhesus_T_NK_Cell_Correction.py):
the original label and old gamma-delta flag are preserved in obs, the full flip
breakdown, transition matrix, and gamma-delta accounting are printed BEFORE the
write cell, the canonical object is NOT overwritten, and output goes to a NEW
file. In Spyder, run the report cells, read them, then run the write cell. Under
SLURM it runs end-to-end but only ever writes the new file, so it is reversible.

This is step 1 of 3. Step 2 re-runs shiv_tier2_diagnostic.py against the new
object; step 3 (after we read the new tables) writes final subtype labels and
runs the cell-level DEG.

Run WITH python (sc_pre env). Spyder cells (# %%). Headless-safe under SLURM.

Author: Jake Lehle / Kaushal Lab
Date: July 2026
"""

# %% ============================================================
# CELL 1 - CONFIG AND DIALS
# ============================================================
import os
os.environ.setdefault('MPLBACKEND', 'Agg')
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
ADATA_IN    = os.path.join(WORKING_DIR, 'shiv_starsolo', 'merged',
                           'shiv_host_final_embedded.h5ad')
ADATA_OUT   = os.path.join(WORKING_DIR, 'shiv_starsolo', 'merged',
                           'shiv_host_final_tcrTrecovered.h5ad')
OUT_DIR     = os.path.join(WORKING_DIR, 'annotation_output', 'tier1_tcr_recovery')

# ---- behaviour ----
DRY_RUN = False            # True = run every report but skip the write cell

# ---- columns (all present on the canonical object) ----
CELLTYPE_COL   = 'tier1_celltype'
PRE_COL        = 'tier1_celltype_pre_recovery'   # NEW: preserved original label
TCR_COL        = 'has_tcr'
GD_COL         = 'is_gamma_delta'
SUBTYPE_COL    = 'receptor_subtype'              # 'TRG+TRD' = strict gamma-delta
PAIRING_COL    = 'chain_pairing'                 # clonotype-definable vs orphan
ANIMAL_COL     = 'animal_id'
TIME_COL       = 'timepoint'
TISSUE_COL     = 'tissue'
UMAP_KEY       = 'X_umap'
ELITE          = '40707'

# ---- recovery rule ----
TARGET_LABEL   = 'T cells'
QUARANTINE     = {'B cells', 'Monocytes', 'Macrophages', 'DCs'}

# chain_pairing values that mean a complete, clonotype-definable pair
DEFINABLE_PAIRING = {'single pair', 'extra VJ', 'extra VDJ', 'two full chains'}

# strict gamma-delta = called T AND this receptor_subtype
STRICT_GD_SUBTYPE = 'TRG+TRD'

# ---- ordering / colors (all hex per lab convention) ----
CATEGORY_ORDER = ['T cells', 'B cells', 'NK cells', 'Monocytes', 'Macrophages',
                  'DCs', 'Mast cells', 'Erythrocytes', 'Platelets', 'Unassigned']
CELLTYPE_COLORS = {
    'T cells':      '#D32F2F',
    'B cells':      '#1976D2',
    'NK cells':     '#7B1FA2',
    'Monocytes':    '#F9A825',
    'Macrophages':  '#00897B',
    'DCs':          '#6D4C41',
    'Mast cells':   '#455A64',
    'Erythrocytes': '#8D6E63',
    'Platelets':    '#C2185B',
    'Unassigned':   '#BDBDBD',
}
TIME_ORDER   = ['Pre', '21 DPI', 'Necropsy', 'Non-Infected', 'Unknown']
FLIP_COL     = '#D32F2F'
KEEP_COL     = '#90A4AE'
ELITE_COL    = '#F9A825'
UNKNOWN_COL  = '#455A64'

plt.rcParams.update({
    'font.size': 28, 'axes.titlesize': 34, 'axes.labelsize': 30,
    'xtick.labelsize': 24, 'ytick.labelsize': 24, 'legend.fontsize': 18,
    'figure.dpi': 100, 'pdf.fonttype': 42, 'ps.fonttype': 42,
})
os.makedirs(OUT_DIR, exist_ok=True)
sc.settings.figdir = OUT_DIR
pd.set_option('display.width', 200)
pd.set_option('display.max_columns', 40)
pd.set_option('display.max_rows', 100)


def savefig(fig, name):
    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(OUT_DIR, f'{name}.{ext}'),
                    bbox_inches='tight', dpi=300)
    print(f"    saved: {name}.pdf / .png")


def save_csv(df, name, index=True):
    df.to_csv(os.path.join(OUT_DIR, name), index=index)
    print(f"    saved: {name}")


def ordered_present(values, order):
    seen = set(pd.Series(values).dropna().astype(str).unique())
    return [c for c in order if c in seen]


def apply_celltype_colors(adata_obj, key):
    """Fix category order + hex colors so before/after figures match."""
    cats = ordered_present(adata_obj.obs[key], CATEGORY_ORDER)
    adata_obj.obs[key] = pd.Categorical(adata_obj.obs[key].astype(str),
                                        categories=cats, ordered=True)
    adata_obj.uns[f'{key}_colors'] = [CELLTYPE_COLORS.get(c, '#000000')
                                      for c in cats]
    return cats


# %% ============================================================
# CELL 2 - LOAD + PRESERVE ORIGINAL LABEL
# ============================================================
print("=" * 80)
print("CELL 2 - LOAD + PRESERVE")
print("=" * 80)
print(f"  input: {ADATA_IN}")

adata = sc.read_h5ad(ADATA_IN)
print(f"  loaded: {adata.n_obs:,} cells x {adata.n_vars:,} genes")

for col in (CELLTYPE_COL, TCR_COL, ANIMAL_COL):
    if col not in adata.obs.columns:
        raise KeyError(f"required obs column missing: {col}")

# preserve the original canonical Tier 1 label before any edit
adata.obs[PRE_COL] = adata.obs[CELLTYPE_COL].astype(str).values

cur = adata.obs[CELLTYPE_COL].astype(str)
print("\n  Tier 1 distribution (as loaded):")
print(cur.value_counts().to_string())


# %% ============================================================
# CELL 3 - COMPUTE FLIP CANDIDATES + FULL BREAKDOWN (read before write)
# ============================================================
print("\n" + "=" * 80)
print("CELL 3 - FLIP CANDIDATES")
print("=" * 80)

has_tcr = adata.obs[TCR_COL].astype(bool).values
is_T    = (cur == TARGET_LABEL).values
in_quar = cur.isin(QUARANTINE).values
flip_mask = has_tcr & (~is_T) & (~in_quar)

n_flip = int(flip_mask.sum())
print(f"  rule: has_tcr AND not '{TARGET_LABEL}' AND not in "
      f"{sorted(QUARANTINE)}")
print(f"  cells flipping to T: {n_flip:,}")

# breakdown by source lineage
src = cur[flip_mask]
print("\n  flips by SOURCE lineage:")
print(src.value_counts().to_string())

# quarantine sanity: how many TCR+ cells we deliberately did NOT flip
held = has_tcr & (~is_T) & in_quar
print(f"\n  TCR+ cells held in B/myeloid quarantine (NOT flipped): "
      f"{int(held.sum()):,}")
if int(held.sum()):
    print(cur[held].value_counts().to_string())

# source lineage x definable/orphan (trustworthiness of the receptor)
if PAIRING_COL in adata.obs.columns:
    pairing = adata.obs[PAIRING_COL].astype(str)
    definable = pairing.isin(DEFINABLE_PAIRING).values
    ev = np.where(definable, 'clonotype-definable', 'orphan/ambiguous')
    breakdown = pd.crosstab(src.values, ev[flip_mask])
    print("\n  flips by SOURCE lineage x receptor evidence:")
    print(breakdown.to_string())
    save_csv(breakdown, 'tcr_recovery_flip_by_source_evidence.csv')
    n_def = int(definable[flip_mask].sum())
    print(f"\n  of {n_flip:,} flips, {n_def:,} are clonotype-definable "
          f"({100*n_def/max(n_flip,1):.1f}%) -> these add to the clonotype pool.")
else:
    print(f"  WARN '{PAIRING_COL}' not in obs; skipping evidence breakdown.")

# annotate provenance on the object (used by the v2 subset diagnostic)
adata.obs['is_tcr_recovered'] = flip_mask
adata.obs['tcr_recovery_source'] = np.where(flip_mask, adata.obs[PRE_COL], '')


# %% ============================================================
# CELL 4 - PER-ANIMAL / UNKNOWN CHECK + GAMMA-DELTA PREVIEW
# ============================================================
print("\n" + "=" * 80)
print("CELL 4 - PER-ANIMAL CHECK + GAMMA-DELTA PREVIEW")
print("=" * 80)

animal = adata.obs[ANIMAL_COL].astype(str)

# per-animal flip counts (is the recovery skewed to any one animal?)
flip_by_animal = animal[flip_mask].value_counts()
Tbase_by_animal = animal[is_T].value_counts()
per_animal = pd.DataFrame({
    'T_before': Tbase_by_animal,
    'flipped_in': flip_by_animal,
}).fillna(0).astype(int)
per_animal['T_after'] = per_animal['T_before'] + per_animal['flipped_in']
per_animal['pct_growth'] = (100 * per_animal['flipped_in']
                            / per_animal['T_before'].replace(0, np.nan)).round(1)
per_animal['elite'] = per_animal.index == ELITE
per_animal = per_animal.sort_values('flipped_in', ascending=False)
print("  per-animal T recovery (watch skew and 40707):")
print(per_animal.to_string())
save_csv(per_animal, 'tcr_recovery_per_animal.csv')

n_unknown = int((animal[flip_mask] == 'Unknown').sum())
print(f"\n  flipped cells with animal_id == 'Unknown': {n_unknown:,} "
      f"({100*n_unknown/max(n_flip,1):.1f}% of flips)")
print("  NOTE: Unknown-animal cells pad the global UMAP / composition figure but")
print("        contribute nothing to the animal-resolved DEG or the 40707 story.")

# strict gamma-delta PREVIEW (the flag itself is recomputed + written in CELL 6)
if SUBTYPE_COL in adata.obs.columns:
    is_gd_strict = (adata.obs[SUBTYPE_COL].astype(str) == STRICT_GD_SUBTYPE).values
    gd_T_before = int((is_gd_strict & is_T).sum())
    gd_T_after  = int((is_gd_strict & (is_T | flip_mask)).sum())
    print(f"\n  strict paired TRG+TRD among T:  before {gd_T_before:,}  ->  "
          f"after {gd_T_after:,}  (+{gd_T_after - gd_T_before:,})")
if GD_COL in adata.obs.columns:
    is_gd_flag = adata.obs[GD_COL].astype(bool).values
    print(f"  OLD is_gamma_delta flag among T: before "
          f"{int((is_gd_flag & is_T).sum()):,}  ->  after "
          f"{int((is_gd_flag & (is_T | flip_mask)).sum()):,}  "
          f"(stale: flag never saw the flipped cells -> recomputed in CELL 6)")


# %% ============================================================
# CELL 5 - APPLY FLIP + TRANSITION MATRIX
# ============================================================
print("\n" + "=" * 80)
print("CELL 5 - APPLY FLIP")
print("=" * 80)

new_label = cur.copy()
new_label[flip_mask] = TARGET_LABEL
adata.obs[CELLTYPE_COL] = new_label.values  # re-categorized in CELL 7 for colors

trans = pd.crosstab(adata.obs[PRE_COL], adata.obs[CELLTYPE_COL],
                    rownames=['pre_recovery'], colnames=['post_recovery'])
print("  transition matrix (pre_recovery -> post_recovery):")
print(trans.to_string())
save_csv(trans, 'tcr_recovery_transition_matrix.csv')

before_ct = adata.obs[PRE_COL].value_counts()
after_ct  = adata.obs[CELLTYPE_COL].value_counts()
delta = pd.DataFrame({'before': before_ct, 'after': after_ct}).fillna(0).astype(int)
delta['delta'] = delta['after'] - delta['before']
print("\n  Tier 1 distribution before vs after:")
print(delta.to_string())
save_csv(delta, 'tcr_recovery_celltype_counts.csv')

una_before = int(before_ct.get('Unassigned', 0))
una_after  = int(after_ct.get('Unassigned', 0))
print(f"\n  Unassigned: {una_before:,} -> {una_after:,} "
      f"({una_before - una_after:,} recovered)")
print(f"  T cells:    {int(before_ct.get('T cells', 0)):,} -> "
      f"{int(after_ct.get('T cells', 0)):,}")


# %% ============================================================
# CELL 6 - RECOMPUTE is_gamma_delta ON THE RECOVERED POOL
#
# The old flag was set at VDJ-integration time and only ever flagged cells that
# were T THEN (368). Recovery added paired TRG+TRD cells to the T pool, so the
# flag is stale. Recompute from the stable receptor property.
#
# strict gamma-delta T = called T (post-flip) AND receptor_subtype == 'TRG+TRD'.
# receptor_subtype 'TRG+TRD' IS the strict paired set; the copies that remain in
# the B/myeloid quarantine are contamination by our own locked decision and are
# NOT gamma-delta T cells, so we gate on the T label. Accounting must close:
#   (gamma-delta T) + (TRG+TRD held in non-T) == (TRG+TRD total).
# ============================================================
print("\n" + "=" * 80)
print("CELL 6 - RECOMPUTE is_gamma_delta")
print("=" * 80)

if SUBTYPE_COL in adata.obs.columns:
    is_trgtrd = (adata.obs[SUBTYPE_COL].astype(str) == STRICT_GD_SUBTYPE).values
    is_T_now  = (adata.obs[CELLTYPE_COL].astype(str) == TARGET_LABEL).values  # post-flip
    old_flag  = (adata.obs[GD_COL].astype(bool).values if GD_COL in adata.obs.columns
                 else np.zeros(adata.n_obs, bool))
    new_flag  = is_T_now & is_trgtrd

    # preserve the old flag, write the recomputed one
    adata.obs['is_gamma_delta_pre_recovery'] = old_flag
    adata.obs[GD_COL] = new_flag

    n_total = int(is_trgtrd.sum())
    n_gdT   = int(new_flag.sum())
    n_nonT  = int((is_trgtrd & ~is_T_now).sum())
    print(f"  receptor_subtype == '{STRICT_GD_SUBTYPE}' total : {n_total:,}")
    print(f"  strict gamma-delta T (recomputed)       : {n_gdT:,}")
    print(f"  TRG+TRD held in non-T quarantine        : {n_nonT:,}")
    print(f"  accounting: {n_gdT:,} + {n_nonT:,} = {n_gdT + n_nonT:,} "
          f"(should equal {n_total:,})")
    print(f"  old is_gamma_delta flag was True on     : {int(old_flag.sum()):,} "
          f"cells (stale)")

    nonT_src = adata.obs[CELLTYPE_COL].astype(str)[is_trgtrd & ~is_T_now]
    if len(nonT_src):
        print("  non-T TRG+TRD by lineage (expect B/myeloid only):")
        print(nonT_src.value_counts().to_string())

    save_csv(pd.DataFrame({
        'metric': ['TRG+TRD_total', 'gamma_delta_T_recomputed',
                   'TRG+TRD_nonT_quarantine', 'old_flag_true'],
        'count':  [n_total, n_gdT, n_nonT, int(old_flag.sum())],
    }), 'tcr_recovery_gamma_delta_accounting.csv', index=False)
    print("  NOTE: downstream scripts should key gamma-delta on is_gamma_delta")
    print("        (now 413), which equals receptor_subtype=='TRG+TRD' among T.")
else:
    print(f"  WARN '{SUBTYPE_COL}' not in obs; cannot recompute is_gamma_delta.")


# %% ============================================================
# CELL 7 - FIGURES (the ones that visibly change from the flip)
# ============================================================
print("\n" + "=" * 80)
print("CELL 7 - FIGURES")
print("=" * 80)

# fix category order + hex colors on both label columns for matched figures
cats_pre  = apply_celltype_colors(adata, PRE_COL)
cats_post = apply_celltype_colors(adata, CELLTYPE_COL)

# ---- Fig 1: global UMAP, tier1 before vs after ----
if UMAP_KEY in adata.obsm:
    try:
        fig, axes = plt.subplots(1, 2, figsize=(26, 12))
        for ax, key, ttl in zip(
                axes, [PRE_COL, CELLTYPE_COL],
                ['Tier 1 (before recovery)', 'Tier 1 (after TCR recovery)']):
            sc.pl.umap(adata, color=key, ax=ax, show=False, title=ttl,
                       frameon=False, size=10, legend_loc='right margin',
                       legend_fontsize=16)
        plt.tight_layout()
        savefig(fig, 'tcr_recovery_umap_prepost')
        plt.close(fig)
    except Exception as e:
        print(f"  WARN Fig 1 skipped: {e}")
else:
    print(f"  WARN '{UMAP_KEY}' not in obsm; skipping UMAP figure.")

# ---- Fig 2: Tier 1 counts before vs after (grouped bar) ----
try:
    order = ordered_present(list(delta.index), CATEGORY_ORDER)
    d = delta.reindex(order)
    x = np.arange(len(order)); w = 0.4
    fig, ax = plt.subplots(figsize=(18, 10))
    ax.bar(x - w/2, d['before'], width=w, color=KEEP_COL, label='before')
    ax.bar(x + w/2, d['after'],  width=w, color=FLIP_COL, label='after')
    ax.set_xticks(x); ax.set_xticklabels(order, rotation=45, ha='right')
    ax.set_ylabel('cells')
    ax.set_title('Tier 1 cell counts, TCR T recovery')
    ax.legend(); ax.spines[['top', 'right']].set_visible(False)
    plt.tight_layout()
    savefig(fig, 'tcr_recovery_celltype_counts')
    plt.close(fig)
except Exception as e:
    print(f"  WARN Fig 2 skipped: {e}")

# ---- Fig 3: per-animal flipped cells (the quick check, 40707 flagged) ----
try:
    pa = per_animal.drop(index='Unknown', errors='ignore').copy()
    pa = pa.sort_values('flipped_in', ascending=False)
    x = np.arange(len(pa))
    colors = [ELITE_COL if e else KEEP_COL for e in pa['elite']]
    fig, ax = plt.subplots(figsize=(16, 10))
    ax.bar(x, pa['flipped_in'], color=colors)
    if 'Unknown' in per_animal.index:
        ax.bar(len(pa), int(per_animal.loc['Unknown', 'flipped_in']),
               color=UNKNOWN_COL)
    xt = list(pa.index) + (['Unknown'] if 'Unknown' in per_animal.index else [])
    labs = [f"{a}*" if a == ELITE else a for a in xt]
    ax.set_xticks(range(len(xt))); ax.set_xticklabels(labs, rotation=45, ha='right')
    ax.set_ylabel('cells flipped to T')
    ax.set_title('TCR-recovered T cells per animal\n(* = elite controller 40707)')
    ax.spines[['top', 'right']].set_visible(False)
    plt.tight_layout()
    savefig(fig, 'tcr_recovery_per_animal_flips')
    plt.close(fig)
except Exception as e:
    print(f"  WARN Fig 3 skipped: {e}")

# ---- Fig 4: composition by timepoint x tissue, after recovery (stacked prop) ----
try:
    df = adata.obs[[TIME_COL, TISSUE_COL, CELLTYPE_COL]].astype(str).copy()
    df['grp'] = df[TIME_COL] + '\n' + df[TISSUE_COL]
    ct = pd.crosstab(df['grp'], df[CELLTYPE_COL])
    ct = ct.div(ct.sum(axis=1), axis=0)
    col_order = ordered_present(list(ct.columns), CATEGORY_ORDER)
    ct = ct[col_order]
    grp_order = [f"{t}\n{ti}" for t in TIME_ORDER for ti in ('LN', 'PBMC')
                 if f"{t}\n{ti}" in ct.index]
    grp_order += [g for g in ct.index if g not in grp_order]
    ct = ct.reindex(grp_order)
    fig, ax = plt.subplots(figsize=(20, 11))
    bottom = np.zeros(len(ct))
    for c in col_order:
        ax.bar(range(len(ct)), ct[c].values, bottom=bottom,
               color=CELLTYPE_COLORS.get(c, '#000000'), label=c, width=0.8)
        bottom += ct[c].values
    ax.set_xticks(range(len(ct))); ax.set_xticklabels(ct.index)
    ax.set_ylabel('fraction of cells'); ax.set_ylim(0, 1)
    ax.set_title('Tier 1 composition by timepoint x tissue (after recovery)')
    ax.legend(bbox_to_anchor=(1.01, 1), loc='upper left', fontsize=16)
    ax.spines[['top', 'right']].set_visible(False)
    plt.tight_layout()
    savefig(fig, 'tcr_recovery_composition_timepoint_tissue')
    plt.close(fig)
    save_csv(ct.round(4), 'tcr_recovery_composition_timepoint_tissue.csv')
except Exception as e:
    print(f"  WARN Fig 4 skipped: {e}")


# %% ============================================================
# CELL 8 - WRITE (new file; canonical object NOT overwritten)
# ============================================================
print("\n" + "=" * 80)
print("CELL 8 - WRITE")
print("=" * 80)

if DRY_RUN:
    print("  DRY_RUN = True -> nothing written. Set DRY_RUN = False to commit.")
else:
    adata.write_h5ad(ADATA_OUT)
    print(f"  [SAVED] {ADATA_OUT}")
    print(f"  new/updated obs:")
    print(f"    {PRE_COL} (preserved original Tier 1 label)")
    print(f"    is_tcr_recovered, tcr_recovery_source (recovery provenance)")
    print(f"    {GD_COL} (RECOMPUTED strict gamma-delta on recovered pool)")
    print(f"    is_gamma_delta_pre_recovery (preserved old flag)")
    print(f"  {CELLTYPE_COL} now includes {n_flip:,} recovered T cells")
    print(f"  canonical object untouched at {os.path.basename(ADATA_IN)}")
    print(f"""
  NEXT (step 2): re-run shiv_tier2_diagnostic.py with
      ADATA_IN = {os.path.basename(ADATA_OUT)}
      ADATA_OUT = shiv_Tcell_subset_reclustered_v2.h5ad
  Leiden will renumber; we re-map cluster ID -> label once against the new
  tables. The naive/CM/EM/terminal/proliferating axis and the tissue/timepoint
  gradients are driven by tens of thousands of cells and are not erased by
  adding ~{n_flip:,} mostly-CD8-dropout cells, so the DEG design carries over.
""")

print("=" * 80)
print("TIER 1 TCR RECOVERY COMPLETE - canonical object untouched.")
print("=" * 80)
