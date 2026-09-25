#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Rhesus SHIV — VDJ (TCR) Integration onto the Tier-1-Corrected Object
=====================================================================
scirpy 0.22.5 workflow (MuData-based; the old merge_with_tcr is gone).

Purpose
-------
1. Attach TCR receptor status to the corrected object.
2. Resolve provisional NK: a provisional NK cell that carries a real productive
   TCR is definitively a T cell -> reassign to T.
3. Call gamma-delta T STRICTLY from a paired TRG+TRD receptor_subtype. The
   inventory showed TRG recovered broadly but TRD sparse, and TRG persists in
   alpha-beta T cells, so "any TRG/TRD" (2,982) overcounts gamma-delta. chain_qc
   classifies alpha-beta cells on their TRB, so the paired TRG+TRD subtype is
   the honest gamma-delta call (a floor, since TRD dropout is real).

What TCR can and cannot do here (from the inventory):
  - ~39% of GEX cells carry a TCR, and C/F map poorly because QC cut those
    pooled control-bearing libraries hard. So a productive TCR reassigns a
    mislabeled provisional NK to T (strong), but TCR-NEGATIVE stays SUPPORTING
    evidence for NK, not proof, because of dropout. TCR-negative provisional NK
    therefore remain provisional.

Barcode join: VDJ obs_names are '<bare>-1'; GEX obs_names are '<bare>-1-<lib>'.
We append the library suffix so they align, then SUBSET AIRR to cells present
in GEX before building the MuData (drops the ~24% of TCR cells whose GEX
partner QC removed, keeping the modalities aligned).

Preserves pre-TCR labels; writes both .h5ad (for the pipeline) and .h5mu (for
the clonotype work next). Read the receptor_subtype distribution and the
provisional-NK-with-TCR count before treating labels as final.

Run WITH python (sc_pre env). Spyder cells (# %%).

Author: Jake Lehle / Kaushal Lab
Date: July 2026
"""

# %% ============================================================
# CELL 1 — CONFIGURATION AND DIALS
# ============================================================
import os
import glob
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

WORKING_DIR    = '/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV'
ANALYSIS_DIR   = os.path.join(WORKING_DIR, 'Analysis')
CORRECTED_H5AD = os.path.join(WORKING_DIR, 'shiv_starsolo', 'merged',
                              'shiv_host_tier1_corrected.h5ad')
OUT_H5AD       = os.path.join(WORKING_DIR, 'shiv_starsolo', 'merged',
                              'shiv_host_tier1_tcr.h5ad')
OUT_H5MU       = os.path.join(WORKING_DIR, 'shiv_starsolo', 'merged',
                              'shiv_host_tcr.h5mu')
OUT_DIR        = os.path.join(WORKING_DIR, 'annotation_output', 'vdj_integration')

LIBRARY_COL = 'library'
ANIMAL_COL  = 'animal_id'
ELITE_CONTROLLER = '40707'

LIB_NAME = {
    'A': 'SHIV_Pre_PBMC',             'B': 'SHIV_Pre_LN',
    'C': 'SHIV_Pre_PBMC_IndianMixed', 'D': 'SHIV_21DPI_Pre_LN_Mixed',
    'E': 'SHIV_21DPI_PBMC',           'F': 'SHIV_21DPI_NonInf_LN_Mixed',
    'G': 'SHIV_Necropsy_PBMC',        'H': 'SHIV_Necropsy_LN',
}
LIBS = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
VDJ_SUBDIR = 'vdj_t'

# Which receptor_subtype values count as a confident gamma-delta call.
GD_SUBTYPES = ['TRG+TRD']          # strict, paired
GD_SUBTYPES_LOOSE = ['TRG+TRD', 'TRD', 'TRG']  # reported for transparency only

plt.rcParams.update({
    'font.size': 28, 'axes.titlesize': 34, 'axes.labelsize': 30,
    'xtick.labelsize': 28, 'ytick.labelsize': 28, 'legend.fontsize': 26,
    'figure.dpi': 100,
})
GD_COLOR = '#F9A825'
ELITE_COL = '#F9A825'
BAR_COLOR = '#1976D2'

os.makedirs(OUT_DIR, exist_ok=True)


def savefig(fig, name):
    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(OUT_DIR, f'{name}.{ext}'),
                    bbox_inches='tight', dpi=300)
    print(f"    saved figure: {name}.pdf / .png")


# %% ============================================================
# CELL 2 — READ VDJ PER LIBRARY, HARMONIZE BARCODES, CONCAT
# ============================================================
print("=" * 80)
print("CELL 2 — READ + HARMONIZE VDJ")
print("=" * 80)

import scirpy as ir
print(f"  scirpy version: {ir.__version__}")


def find_vdj_file(sample_name, fname):
    p = os.path.join(ANALYSIS_DIR, sample_name, 'outs', 'per_sample_outs',
                     sample_name, VDJ_SUBDIR, fname)
    if os.path.exists(p):
        return p
    hits = glob.glob(os.path.join(ANALYSIS_DIR, sample_name, 'outs', '**',
                                  VDJ_SUBDIR, fname), recursive=True)
    return hits[0] if hits else None


airr_list = []
for L in LIBS:
    name = LIB_NAME[L]
    filt = find_vdj_file(name, 'filtered_contig_annotations.csv')
    if not filt:
        print(f"  [WARN] {L} {name}: no filtered_contig_annotations.csv")
        continue
    a = ir.io.read_10x_vdj(filt)
    # append the library suffix so barcodes align with GEX '<bare>-1-<lib>'
    a.obs_names = [f"{bc}-{L}" for bc in a.obs_names]
    a.obs['library'] = L
    airr_list.append(a)
    print(f"  {L} {name:32s}: {a.n_obs:,} VDJ cells")

if not airr_list:
    raise RuntimeError("No VDJ files read.")

airr = ad.concat(airr_list, join='outer', merge='same', index_unique=None)
print(f"\n  combined AIRR cells: {airr.n_obs:,}")


# %% ============================================================
# CELL 3 — LOAD GEX, SUBSET AIRR TO GEX CELLS, BUILD MuData
# ============================================================
print("\n" + "=" * 80)
print("CELL 3 — BUILD MuData")
print("=" * 80)

from mudata import MuData

gex = sc.read_h5ad(CORRECTED_H5AD)
print(f"  GEX (corrected): {gex.n_obs:,} cells")

# Keep only AIRR cells that exist in GEX (drops TCR cells whose GEX partner was
# QC-removed; reported so the loss is explicit).
in_gex = airr.obs_names.isin(gex.obs_names)
print(f"  AIRR cells mapping to GEX: {int(in_gex.sum()):,}/{airr.n_obs:,} "
      f"({100*in_gex.mean():.1f}%)")
airr = airr[in_gex].copy()

mdata = MuData({"gex": gex, "airr": airr})
print(f"  MuData built: gex={gex.n_obs:,}, airr={airr.n_obs:,}")


# %% ============================================================
# CELL 4 — index_chains + chain_qc
# ============================================================
print("\n" + "=" * 80)
print("CELL 4 — index_chains + chain_qc")
print("=" * 80)

ir.pp.index_chains(mdata)
ir.tl.chain_qc(mdata)

airr_obs = mdata["airr"].obs
for col in ['receptor_type', 'receptor_subtype', 'chain_pairing']:
    if col in airr_obs.columns:
        print(f"\n  {col} distribution (TCR cells):")
        print(airr_obs[col].value_counts(dropna=False).to_string())
    else:
        print(f"  WARN chain_qc did not produce '{col}'")


# %% ============================================================
# CELL 5 — PULL TCR STATUS ONTO THE GEX AnnData
# ============================================================
print("\n" + "=" * 80)
print("CELL 5 — ATTACH TCR STATUS TO GEX obs")
print("=" * 80)


def pull(col, default):
    s = airr_obs[col] if col in airr_obs.columns else pd.Series(dtype=object)
    return s.reindex(gex.obs_names).astype(object).where(
        lambda x: x.notna(), default)


gex.obs['receptor_type'] = pull('receptor_type', 'no TCR').values
gex.obs['receptor_subtype'] = pull('receptor_subtype', 'no TCR').values
gex.obs['chain_pairing'] = pull('chain_pairing', 'no TCR').values
gex.obs['has_tcr'] = (gex.obs['receptor_type'] == 'TCR')

n_tcr = int(gex.obs['has_tcr'].sum())
print(f"  cells with a productive TCR (receptor_type == 'TCR'): "
      f"{n_tcr:,} ({100*n_tcr/gex.n_obs:.1f}%)")
print("\n  has_tcr vs current tier1_celltype:")
print(pd.crosstab(gex.obs['tier1_celltype'], gex.obs['has_tcr']).to_string())


# %% ============================================================
# CELL 6 — RESOLVE PROVISIONAL NK WITH A REAL TCR -> T
# ============================================================
print("\n" + "=" * 80)
print("CELL 6 — RESOLVE PROVISIONAL NK")
print("=" * 80)

gex.obs['tier1_celltype_pre_tcr'] = gex.obs['tier1_celltype'].copy()

prov_nk = gex.obs['tier1_celltype'] == 'NK cells'
nk_with_tcr = prov_nk & gex.obs['has_tcr']

print(f"  provisional NK cells            : {int(prov_nk.sum()):,}")
print(f"  ...carrying a productive TCR    : {int(nk_with_tcr.sum()):,} "
      f"-> reassigned to T (definitively T)")
print(f"  ...TCR-negative (stay NK)       : {int((prov_nk & ~gex.obs['has_tcr']).sum()):,}")
print("  (TCR-negative is supporting evidence for NK, not proof; dropout means")
print("   some of these are T cells whose TCR was not recovered.)")

lab = gex.obs['tier1_celltype'].astype(str)
lab[nk_with_tcr.values] = 'T cells'
gex.obs['tier1_celltype'] = pd.Categorical(lab)

# NK confidence flag: retained NK are CD3-negative + exclusive-scored + now
# TCR-negative. Still not proof, but the strongest NK call available pre-clonotype.
gex.obs['nk_tcr_negative'] = (gex.obs['tier1_celltype'] == 'NK cells')
print(f"\n  NK after TCR resolution         : "
      f"{int((gex.obs['tier1_celltype'] == 'NK cells').sum()):,}")


# %% ============================================================
# CELL 7 — GAMMA-DELTA CALL (strict TRG+TRD) + per animal
# ============================================================
print("\n" + "=" * 80)
print("CELL 7 — GAMMA-DELTA")
print("=" * 80)

sub = gex.obs['receptor_subtype'].astype(str)
gex.obs['is_gamma_delta'] = sub.isin(GD_SUBTYPES) & (gex.obs['tier1_celltype'] == 'T cells')
gex.obs['t_subtype'] = np.where(
    gex.obs['is_gamma_delta'], 'gamma-delta T',
    np.where(gex.obs['tier1_celltype'] == 'T cells', 'alpha-beta T', 'not T'))

n_gd_strict = int(gex.obs['is_gamma_delta'].sum())
n_gd_loose = int((sub.isin(GD_SUBTYPES_LOOSE) & (gex.obs['tier1_celltype'] == 'T cells')).sum())
print(f"  gamma-delta T (strict TRG+TRD)     : {n_gd_strict:,}")
print(f"  gamma-delta-ish (incl. lone TRG/TRD): {n_gd_loose:,}  "
      f"(the gap is mostly TRG-only alpha-beta background)")
print("  NOTE: strict count is a floor; TRD dropout means true gamma-delta is")
print("        somewhat higher, but strict is the defensible number for Zoey.")

if ANIMAL_COL in gex.obs.columns:
    gd = gex.obs.loc[gex.obs['is_gamma_delta']]
    per_animal = gd[ANIMAL_COL].astype(str).value_counts()
    print("\n  gamma-delta T per animal (strict):")
    for a in sorted(per_animal.index):
        flag = ' (elite controller)' if a == ELITE_CONTROLLER else ''
        print(f"    {a}: {per_animal[a]}{flag}")

    try:
        fig, ax = plt.subplots(figsize=(13, 9))
        order = per_animal.sort_values(ascending=False)
        colors = [ELITE_COL if a == ELITE_CONTROLLER else BAR_COLOR for a in order.index]
        ax.bar(range(len(order)), order.values, color=colors)
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels([f"{a}*" if a == ELITE_CONTROLLER else a
                            for a in order.index], rotation=45, ha='right')
        ax.set_ylabel('gamma-delta T cells')
        ax.set_title('Gamma-delta T per animal (strict TRG+TRD)\n(* = elite controller)')
        ax.spines[['top', 'right']].set_visible(False)
        savefig(fig, 'I_gamma_delta_per_animal')
        plt.close(fig)
    except Exception as e:
        print(f"  WARN figure I skipped: {e}")


# %% ============================================================
# CELL 8 — WRITE (.h5ad for the pipeline, .h5mu for clonotypes)
# ============================================================
print("\n" + "=" * 80)
print("CELL 8 — WRITE")
print("=" * 80)

# refresh the MuData gex modality with the new obs, then save both
mdata.mod['gex'] = gex
mdata.update()

gex.write_h5ad(OUT_H5AD)
print(f"  [SAVED] {OUT_H5AD}")
print("          new obs: receptor_type, receptor_subtype, chain_pairing,")
print("          has_tcr, tier1_celltype (post-TCR), tier1_celltype_pre_tcr,")
print("          nk_tcr_negative, is_gamma_delta, t_subtype")

try:
    mdata.write(OUT_H5MU)
    print(f"  [SAVED] {OUT_H5MU}  (AIRR retained for clonotype analysis)")
except Exception as e:
    print(f"  WARN could not write .h5mu: {e}")

print("\n  Review the receptor_subtype distribution (Cell 4), the provisional-NK")
print("  resolution (Cell 6), and the strict gamma-delta count (Cell 7) before")
print("  treating tier1_celltype as final. Clonotype expansion (40707 CD8 story)")
print("  is the next step and builds on the .h5mu.")
print("=" * 80)
