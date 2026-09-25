#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Rhesus SHIV — Final Object Integrity Check  (READ-ONLY)
========================================================
Confirms shiv_host_tier1_tcr.h5ad carries everything it should after the
metadata-preserving merge, the Tier 1 T/NK correction, and the VDJ integration,
and that the derived columns are internally consistent. Reads only; writes
nothing.

Checks:
  1. Shape, var, obsm (protein_counts), uns (adt_names).
  2. All expected obs columns present, grouped by provenance.
  3. Label + receptor distributions.
  4. Consistency cross-checks:
       - every is_gamma_delta cell is a T cell
       - every nk_tcr_negative cell is an NK cell
       - no NK cell still carries has_tcr (all TCR+ NK were moved to T)
       - paired TRG+TRD cells not labeled T (the 109 gamma-delta the RNA missed)
       - SHIV+ counts preserved (shiv_umi / shiv_pos)
       - animal_id preserved, Unknown fraction sane
  5. Ambient-TCR observation: non-lymphoid lineages carrying has_tcr.

Author: Jake Lehle
Date: July 2026
"""

# %% ============================================================
# CELL 1 — CONFIG
# ============================================================
import os
import numpy as np
import pandas as pd
import scanpy as sc
import warnings
warnings.filterwarnings("ignore")

WORKING_DIR = '/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV'
ADATA_IN = os.path.join(WORKING_DIR, 'shiv_starsolo', 'merged',
                        'shiv_host_tier1_tcr.h5ad')

EXPECTED_OBS = {
    'metadata':  ['condition', 'library', 'tissue', 'animal_id',
                  'timepoint', 'subspecies', 'sample_name'],
    'qc':        ['pct_counts_mt', 'doublet_score', 'predicted_doublet'],
    'shiv':      ['shiv_umi', 'shiv_pos'],
    'tier1':     ['tier1_celltype', 'tier1_celltype_uncorrected',
                  'tier1_celltype_pre_tcr', 'nk_provisional', 'nk_tcr_negative'],
    'tcr':       ['receptor_type', 'receptor_subtype', 'chain_pairing',
                  'has_tcr', 'is_gamma_delta', 't_subtype'],
}
ELITE_CONTROLLER = '40707'


def ok(cond):
    return 'OK  ' if cond else 'FAIL'


# %% ============================================================
# CELL 2 — LOAD + STRUCTURE
# ============================================================
print("=" * 80)
print("CELL 2 — STRUCTURE")
print("=" * 80)
adata = sc.read_h5ad(ADATA_IN)
print(f"  {ADATA_IN}")
print(f"  shape: {adata.n_obs:,} cells x {adata.n_vars:,} genes")

adt_names = list(adata.uns.get('adt_names', []))
has_prot = 'protein_counts' in adata.obsm
print(f"  obsm['protein_counts']: {ok(has_prot)} "
      f"{adata.obsm['protein_counts'].shape if has_prot else '(absent)'}")
print(f"  uns['adt_names']: {ok(len(adt_names) > 0)} ({len(adt_names)} tags)")
if has_prot and adt_names:
    print(f"  ADT n cols == n tags: "
          f"{ok(adata.obsm['protein_counts'].shape[1] == len(adt_names))}")
print(f"  var columns: {list(adata.var.columns)[:8]}"
      f"{' ...' if adata.var.shape[1] > 8 else ''}")


# %% ============================================================
# CELL 3 — OBS COLUMN PRESENCE
# ============================================================
print("\n" + "=" * 80)
print("CELL 3 — OBS COLUMNS")
print("=" * 80)
all_present = True
for group, cols in EXPECTED_OBS.items():
    missing = [c for c in cols if c not in adata.obs.columns]
    all_present &= (len(missing) == 0)
    tag = ok(len(missing) == 0)
    print(f"  {tag} {group:9s}: "
          f"{'all present' if not missing else 'MISSING ' + str(missing)}")
print(f"\n  all expected obs present: {ok(all_present)}")


# %% ============================================================
# CELL 4 — DISTRIBUTIONS
# ============================================================
print("\n" + "=" * 80)
print("CELL 4 — DISTRIBUTIONS")
print("=" * 80)
print("  tier1_celltype (final):")
print(adata.obs['tier1_celltype'].value_counts().to_string())
print("\n  t_subtype:")
print(adata.obs['t_subtype'].value_counts().to_string())
print("\n  receptor_subtype:")
print(adata.obs['receptor_subtype'].value_counts().to_string())


# %% ============================================================
# CELL 5 — CONSISTENCY CROSS-CHECKS
# ============================================================
print("\n" + "=" * 80)
print("CELL 5 — CONSISTENCY")
print("=" * 80)

o = adata.obs

gd = o['is_gamma_delta'].astype(bool)
gd_all_T = bool((o.loc[gd, 'tier1_celltype'] == 'T cells').all()) if gd.any() else True
print(f"  {ok(gd_all_T)} every is_gamma_delta cell is a T cell "
      f"({int(gd.sum())} gamma-delta)")

nktn = o['nk_tcr_negative'].astype(bool)
nktn_all_NK = bool((o.loc[nktn, 'tier1_celltype'] == 'NK cells').all()) if nktn.any() else True
print(f"  {ok(nktn_all_NK)} every nk_tcr_negative cell is an NK cell "
      f"({int(nktn.sum())} NK)")

nk_mask = o['tier1_celltype'] == 'NK cells'
nk_with_tcr = int((nk_mask & o['has_tcr'].astype(bool)).sum())
print(f"  {ok(nk_with_tcr == 0)} no NK cell still carries has_tcr "
      f"(found {nk_with_tcr}; should be 0 after resolution)")

# paired gamma-delta the RNA missed (not labeled T)
paired_gd = o['receptor_subtype'].astype(str) == 'TRG+TRD'
gd_not_T = int((paired_gd & (o['tier1_celltype'] != 'T cells')).sum())
print(f"  NOTE {gd_not_T} paired TRG+TRD cells are NOT labeled T "
      f"(real gamma-delta the RNA argmax missed; candidates for a paired pass)")

# SHIV preserved
n_shiv_pos = int(o['shiv_pos'].sum())
print(f"  SHIV+ cells (shiv_pos): {n_shiv_pos}  "
      f"(max shiv_umi {int(o['shiv_umi'].max())})")
print("  SHIV+ by library:")
print(pd.crosstab(o['library'], o['shiv_pos']).to_string())

# animal_id preserved
n_unknown = int((o['animal_id'].astype(str) == 'Unknown').sum())
print(f"\n  animal_id Unknown: {n_unknown} "
      f"({100*n_unknown/adata.n_obs:.1f}%)")
print(f"  {ok(ELITE_CONTROLLER in set(o['animal_id'].astype(str)))} "
      f"elite controller {ELITE_CONTROLLER} present")


# %% ============================================================
# CELL 6 — AMBIENT-TCR OBSERVATION (not an error)
# ============================================================
print("\n" + "=" * 80)
print("CELL 6 — has_tcr BY LINEAGE (ambient watch)")
print("=" * 80)
ct = pd.crosstab(o['tier1_celltype'], o['has_tcr'].astype(bool))
if True in ct.columns:
    ct['pct_tcr'] = (100 * ct[True] / ct.sum(axis=1)).round(1)
print(ct.to_string())
print("\n  Non-lymphoid lineages (B, Monocytes, Macrophages, ...) carrying a")
print("  'TCR' are almost certainly ambient orphan contigs, not real receptors.")
print("  has_tcr is therefore not clean lineage evidence off the T/NK axis.")


# %% ============================================================
# CELL 7 — VERDICT
# ============================================================
print("\n" + "=" * 80)
print("CELL 7 — VERDICT")
print("=" * 80)
checks = [has_prot, len(adt_names) > 0, all_present, gd_all_T,
          nktn_all_NK, nk_with_tcr == 0,
          ELITE_CONTROLLER in set(o['animal_id'].astype(str))]
print(f"  passed {sum(checks)}/{len(checks)} structural + consistency checks")
print(f"  object: {adata.n_obs:,} cells, {len(adt_names)} ADT tags, "
      f"{n_shiv_pos} SHIV+, {int(gd.sum())} gamma-delta T")
print("  Reminder for plotting: no embedding/clusters yet and ADT is raw")
print("  protein_counts (not CLR); CD163/CD169 blank pending BD Ab-seq recovery.")
print("=" * 80)
