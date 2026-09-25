#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Rhesus SHIV — VDJ (TCR) Inventory  (READ-ONLY)
===============================================
Look at what the cellranger multi vdj_t step actually produced, per library,
BEFORE writing the scirpy integration. The rhesus IMGT reference covers all
four TCR loci (TRA/TRB/TRG/TRD), so this is an integration, not a recovery
problem. This inventory confirms that from the OUTPUTS rather than the
reference, and surfaces anything that would trip the integration.

Questions answered (read together before writing the merge):
  1. Where does each library's vdj_t output live, and is it complete?
  2. Recovery: how many cells got a productive TCR contig per library?
  3. Chain coverage in the DATA: counts of productive TRA / TRB / TRG / TRD
     contigs and cells. Are TRG/TRD genuinely recovered (gamma-delta on the
     table for Zoey), or only TRA/TRB?
  4. Barcode join: do VDJ barcodes map onto the corrected object's called
     cells via (library, bare_barcode), and at what yield?
  5. Gamma-delta cells per animal (for the matched cells), via the corrected
     object's animal_id.
  6. Is scirpy installed in sc_pre, and which version (API differs by release)?

READ-ONLY: reads the vdj_t CSVs and the corrected .h5ad, writes only inventory
CSVs and figures to OUT_DIR. Nothing is merged here.

Author: Jake Lehle
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
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

WORKING_DIR   = '/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV'
ANALYSIS_DIR  = os.path.join(WORKING_DIR, 'Analysis')
CORRECTED_H5AD = os.path.join(WORKING_DIR, 'shiv_starsolo', 'merged',
                              'shiv_host_tier1_corrected.h5ad')
OUT_DIR       = os.path.join(WORKING_DIR, 'annotation_output', 'vdj_inventory')

LIBRARY_COL = 'library'
ANIMAL_COL  = 'animal_id'

# letter -> cellranger multi sample name (== per_sample_outs stem)
LIB_NAME = {
    'A': 'SHIV_Pre_PBMC',             'B': 'SHIV_Pre_LN',
    'C': 'SHIV_Pre_PBMC_IndianMixed', 'D': 'SHIV_21DPI_Pre_LN_Mixed',
    'E': 'SHIV_21DPI_PBMC',           'F': 'SHIV_21DPI_NonInf_LN_Mixed',
    'G': 'SHIV_Necropsy_PBMC',        'H': 'SHIV_Necropsy_LN',
}
LIBS = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']

VDJ_SUBDIR = 'vdj_t'   # cellranger multi TCR output folder
CHAINS = ['TRA', 'TRB', 'TRG', 'TRD']

# --- Figure style (28-34 pt, hex, PDF+PNG @300) ---
plt.rcParams.update({
    'font.size': 28, 'axes.titlesize': 34, 'axes.labelsize': 30,
    'xtick.labelsize': 28, 'ytick.labelsize': 28, 'legend.fontsize': 28,
    'figure.dpi': 100,
})
CHAIN_COLORS = {'TRA': '#1976D2', 'TRB': '#D32F2F',
                'TRG': '#F9A825', 'TRD': '#388E3C'}
AB_COLOR = '#1976D2'
GD_COLOR = '#F9A825'

os.makedirs(OUT_DIR, exist_ok=True)


def savefig(fig, name):
    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(OUT_DIR, f'{name}.{ext}'),
                    bbox_inches='tight', dpi=300)
    print(f"    saved figure: {name}.pdf / .png")


def as_bool(series):
    """cellranger writes booleans as True/False or 'true'/'false'; normalize."""
    if series.dtype == bool:
        return series
    return series.astype(str).str.lower().isin(['true', '1', 'yes'])


# %% ============================================================
# CELL 2 — LOCATE per-library vdj_t OUTPUTS
# ============================================================
print("=" * 80)
print("CELL 2 — LOCATE vdj_t OUTPUTS")
print("=" * 80)


def find_vdj_file(sample_name, fname):
    """Locate a vdj_t file under the multi per_sample_outs, with a find fallback."""
    p = os.path.join(ANALYSIS_DIR, sample_name, 'outs', 'per_sample_outs',
                     sample_name, VDJ_SUBDIR, fname)
    if os.path.exists(p):
        return p
    hits = glob.glob(os.path.join(ANALYSIS_DIR, sample_name, 'outs', '**',
                                  VDJ_SUBDIR, fname), recursive=True)
    return hits[0] if hits else None


vdj_paths = {}
for L in LIBS:
    name = LIB_NAME[L]
    filt = find_vdj_file(name, 'filtered_contig_annotations.csv')
    allc = find_vdj_file(name, 'all_contig_annotations.csv')
    airr = find_vdj_file(name, 'airr_rearrangement.tsv')
    vdj_paths[L] = {'filtered': filt, 'all': allc, 'airr': airr}
    status = 'OK' if filt else 'MISSING filtered_contig_annotations.csv'
    print(f"  {L} {name:32s}: {status}")
    if filt:
        print(f"       {filt}")


# %% ============================================================
# CELL 3 — PER-LIBRARY RECOVERY + CHAIN COVERAGE
# ============================================================
print("\n" + "=" * 80)
print("CELL 3 — RECOVERY + CHAIN COVERAGE (from filtered_contig_annotations)")
print("=" * 80)

lib_rows = []
vdj_frames = {}   # keep per-library filtered tables for the join in Cell 4
for L in LIBS:
    filt = vdj_paths[L]['filtered']
    if not filt:
        continue
    df = pd.read_csv(filt)
    vdj_frames[L] = df

    # normalize the flags we rely on
    prod = as_bool(df['productive']) if 'productive' in df.columns else pd.Series(False, index=df.index)
    hi = as_bool(df['high_confidence']) if 'high_confidence' in df.columns else pd.Series(True, index=df.index)
    keep = prod & hi
    dfp = df[keep]

    n_contigs = len(df)
    n_prod = int(keep.sum())
    n_cells = df['barcode'].nunique()
    n_cells_prod = dfp['barcode'].nunique()

    row = {'library': L, 'condition': LIB_NAME[L],
           'contigs': n_contigs, 'productive_contigs': n_prod,
           'cells_any_contig': n_cells, 'cells_productive': n_cells_prod}

    # per-chain productive contig counts and per-cell counts
    for ch in CHAINS:
        ch_mask = dfp['chain'] == ch
        row[f'{ch}_contigs'] = int(ch_mask.sum())
        row[f'{ch}_cells'] = dfp.loc[ch_mask, 'barcode'].nunique()

    # cell-level ab vs gd (descriptive; scirpy chain_qc makes the final call)
    cell_chains = dfp.groupby('barcode')['chain'].apply(set)
    has_ab = cell_chains.apply(lambda s: bool(s & {'TRA', 'TRB'}))
    has_gd = cell_chains.apply(lambda s: bool(s & {'TRG', 'TRD'}))
    row['cells_ab'] = int(has_ab.sum())
    row['cells_gd'] = int(has_gd.sum())
    row['cells_gd_only'] = int((has_gd & ~has_ab).sum())

    lib_rows.append(row)
    print(f"  {L}: {n_cells_prod:,} cells w/ productive TCR "
          f"(TRA={row['TRA_cells']}, TRB={row['TRB_cells']}, "
          f"TRG={row['TRG_cells']}, TRD={row['TRD_cells']}; "
          f"gd-only cells={row['cells_gd_only']})")

vdj_summary = pd.DataFrame(lib_rows).set_index('library')
vdj_summary.to_csv(os.path.join(OUT_DIR, 'vdj_recovery_by_library.csv'))
print("\n  full per-library summary:")
print(vdj_summary.to_string())
print("    saved: vdj_recovery_by_library.csv")

tot_gd = int(vdj_summary['cells_gd'].sum())
tot_gd_only = int(vdj_summary['cells_gd_only'].sum())
print(f"\n  TOTAL gamma-delta cells (any TRG/TRD): {tot_gd:,}  "
      f"(gd-only, no TRA/TRB: {tot_gd_only:,})")
if tot_gd == 0:
    print("  WARN  no TRG/TRD recovered despite the reference covering them; "
          "check the multi config actually enabled gamma-delta or whether the "
          "inner-enrichment primers pulled TRG/TRD.")


# %% ============================================================
# CELL 4 — BARCODE JOIN AGAINST THE CORRECTED OBJECT
# ============================================================
print("\n" + "=" * 80)
print("CELL 4 — BARCODE JOIN vs CORRECTED OBJECT")
print("=" * 80)

join_rows = []
gd_per_animal = {}
if os.path.exists(CORRECTED_H5AD):
    obj = sc.read_h5ad(CORRECTED_H5AD)
    print(f"  corrected object: {obj.n_obs:,} cells")

    obj_lib = (obj.obs[LIBRARY_COL].astype(str).values
               if LIBRARY_COL in obj.obs.columns
               else np.array([b.rsplit('-', 1)[-1] for b in obj.obs_names]))
    obj_bare = np.array([b.split('-')[0] for b in obj.obs_names])
    obj_keys = set(zip(obj_lib, obj_bare))

    # animal lookup per (library, bare_barcode)
    animal_lookup = {}
    if ANIMAL_COL in obj.obs.columns:
        for l, b, a in zip(obj_lib, obj_bare, obj.obs[ANIMAL_COL].astype(str).values):
            animal_lookup[(l, b)] = a

    # per-library: what fraction of GEX cells have TCR, and TCR cells that map
    for L in LIBS:
        if L not in vdj_frames:
            continue
        dfp = vdj_frames[L]
        prod = as_bool(dfp['productive']) if 'productive' in dfp.columns else pd.Series(True, index=dfp.index)
        tcr_bcs = {bc.split('-')[0] for bc in dfp.loc[prod, 'barcode'].unique()}
        matched = {b for b in tcr_bcs if (L, b) in obj_keys}
        n_gex = int((obj_lib == L).sum())
        join_rows.append({
            'library': L,
            'tcr_cells': len(tcr_bcs),
            'tcr_cells_matched': len(matched),
            'gex_cells': n_gex,
            'pct_tcr_mapped': round(100 * len(matched) / max(len(tcr_bcs), 1), 1),
            'pct_gex_with_tcr': round(100 * len(matched) / max(n_gex, 1), 1),
        })
        # gamma-delta cells per animal (matched only)
        cell_chains = dfp[prod].groupby('barcode')['chain'].apply(set)
        for bc, chs in cell_chains.items():
            b = bc.split('-')[0]
            if (L, b) in obj_keys and (chs & {'TRG', 'TRD'}):
                a = animal_lookup.get((L, b), 'Unknown')
                gd_per_animal[a] = gd_per_animal.get(a, 0) + 1

    join = pd.DataFrame(join_rows).set_index('library')
    join.to_csv(os.path.join(OUT_DIR, 'vdj_join_by_library.csv'))
    print("\n  join yield per library:")
    print(join.to_string())
    print("    saved: vdj_join_by_library.csv")

    if gd_per_animal:
        print("\n  gamma-delta cells per animal (matched):")
        for a in sorted(gd_per_animal):
            flag = ' (elite controller)' if a == '40707' else ''
            print(f"    {a}: {gd_per_animal[a]}{flag}")
else:
    print(f"  WARN corrected object not found at {CORRECTED_H5AD}; "
          "run the Tier 1 correction first. Skipping join.")


# %% ============================================================
# CELL 5 — scirpy AVAILABILITY
# ============================================================
print("\n" + "=" * 80)
print("CELL 5 — scirpy CHECK")
print("=" * 80)
try:
    import scirpy as ir
    print(f"  scirpy available: version {ir.__version__}")
    print("  integration path: read_10x_vdj per library -> harmonize barcodes to")
    print("  '<bare>-1-<lib>' -> merge onto corrected obj -> chain_qc for")
    print("  receptor_type (ab vs gd) -> has_tcr / tcr status in obs.")
except Exception as e:
    print(f"  scirpy NOT importable in this env: {e}")
    print("  If missing, install into sc_pre or a dedicated env before integration.")


# %% ============================================================
# CELL 6 — FIGURES
# ============================================================
print("\n" + "=" * 80)
print("CELL 6 — FIGURES")
print("=" * 80)

if lib_rows:
    # Chain composition per library (productive cells per chain)
    try:
        fig, ax = plt.subplots(figsize=(15, 9))
        x = np.arange(len(vdj_summary))
        w = 0.2
        for i, ch in enumerate(CHAINS):
            ax.bar(x + (i - 1.5) * w, vdj_summary[f'{ch}_cells'].values,
                   width=w, color=CHAIN_COLORS[ch], label=ch)
        ax.set_xticks(x)
        ax.set_xticklabels(vdj_summary.index, rotation=0)
        ax.set_xlabel('library')
        ax.set_ylabel('cells with productive chain')
        ax.set_title('TCR chain recovery by library')
        ax.legend(ncol=4)
        ax.spines[['top', 'right']].set_visible(False)
        savefig(fig, 'G_tcr_chain_by_library')
        plt.close(fig)
    except Exception as e:
        print(f"  WARN figure G skipped: {e}")

    # ab vs gd cells per library
    try:
        fig, ax = plt.subplots(figsize=(15, 9))
        x = np.arange(len(vdj_summary))
        w = 0.38
        ax.bar(x - w/2, vdj_summary['cells_ab'].values, width=w,
               color=AB_COLOR, label='alpha/beta')
        ax.bar(x + w/2, vdj_summary['cells_gd'].values, width=w,
               color=GD_COLOR, label='gamma/delta')
        ax.set_xticks(x)
        ax.set_xticklabels(vdj_summary.index)
        ax.set_xlabel('library')
        ax.set_ylabel('cells')
        ax.set_title('alpha/beta vs gamma/delta T cells by library')
        ax.legend()
        ax.spines[['top', 'right']].set_visible(False)
        savefig(fig, 'H_ab_vs_gd_by_library')
        plt.close(fig)
    except Exception as e:
        print(f"  WARN figure H skipped: {e}")


# %% ============================================================
# CELL 7 — HEADLINE SUMMARY
# ============================================================
print("\n" + "=" * 80)
print("CELL 7 — HEADLINE SUMMARY")
print("=" * 80)
if lib_rows:
    print(f"  libraries with vdj_t output : {len(lib_rows)}/{len(LIBS)}")
    print(f"  total cells w/ productive TCR: {int(vdj_summary['cells_productive'].sum()):,}")
    print(f"  total alpha/beta cells       : {int(vdj_summary['cells_ab'].sum()):,}")
    print(f"  total gamma/delta cells      : {tot_gd:,} (gd-only: {tot_gd_only:,})")
    if join_rows:
        j = pd.DataFrame(join_rows)
        print(f"  TCR cells mapping to GEX     : "
              f"{j['tcr_cells_matched'].sum():,}/{j['tcr_cells'].sum():,} "
              f"({100*j['tcr_cells_matched'].sum()/max(j['tcr_cells'].sum(),1):.1f}%)")
print(f"\n  outputs written to: {OUT_DIR}")
print("  Read the recovery, chain coverage, and join yield before we write the")
print("  scirpy integration onto shiv_host_tier1_corrected.h5ad.")
print("=" * 80)
