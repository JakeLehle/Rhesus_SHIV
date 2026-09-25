#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Rhesus SHIV Multi-Omics Processing Pipeline (v2 — Corrected Mapping)
=====================================================================
5' scRNA-seq + CITE-seq + VDJ from Cell Ranger multi outputs.

CRITICAL FIX (May 2026): Library-to-condition mapping corrected.
The original mapping was wrong for 7/8 libraries. Libraries C, D, F
are mixed-condition pools — per-animal metadata is assigned using the
provider's genetic demultiplexing barcode assignments.

Sections:
  0. Setup and configuration
  1. Metadata (corrected library + per-animal mapping)
  2. Read Cell Ranger multi outputs + provider barcode-to-animal lookup
  3. Gene category tagging (mito, ribo, hb, SHIV)
  4. QC (MAD-based + Scrublet doublet removal)
  5. Normalize, embed, batch correct (BBKNN)
  6. Tier 1 annotation (RIRA rhesus markers for PBMC/LN)
  7. Tier 2 annotation (subtypes within lineage)
  8. CITE-seq protein validation + macrophage resolution
  9. SHIV+ cell identification
 10. Biological gene set scoring
 11. Diagnostics, dotplots, proportions
 12. Export

Author: Jake Lehle
Date:   May 2026
"""

# %% ============================================================
# 0. SETUP AND CONFIGURATION
# ============================================================
import os
import glob
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from tqdm import tqdm
from scipy.stats import median_abs_deviation
from scipy import sparse
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings("ignore")

# === ADJUSTABLE PARAMETERS ===
MT_THRESHOLD         = 20
MAD_COUNTS_GENES     = 5
MAD_MITO             = 3
MIN_GENES            = 200
MIN_CELLS_PER_GENE   = 20
EXPECTED_DOUBLET_RATE = 0.08
LEIDEN_RES           = 1.0
N_PCS                = 30
SCORE_THRESHOLD      = 0.1

# === PATHS ===
ANALYSIS_DIR = '/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV/Analysis'
WORKING_DIR  = '/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV'
FIGURES_DIR  = os.path.join(WORKING_DIR, 'figures')
ANNOT_DIR    = os.path.join(WORKING_DIR, 'annotation_output')

# Provider's demultiplexed outputs (for barcode-to-animal mapping)
PROVIDER_DIR = '/master/zwallis/SHARED/10X_ANALYSIS'
SAMPLE_TABLE = os.path.join(WORKING_DIR, 'sample_table.csv')

os.makedirs(FIGURES_DIR, exist_ok=True)
os.makedirs(ANNOT_DIR, exist_ok=True)
os.chdir(WORKING_DIR)
sc.settings.figdir = FIGURES_DIR

# %% ============================================================
# 1. METADATA — Corrected Library + Per-Animal Mapping
#
#    Libraries A-H are PHYSICAL pooling for sequencing.
#    Three libraries (C, D, F) are mixed-condition pools.
#    Per-animal metadata comes from the provider's demultiplexing.
# ============================================================
print("=" * 80)
print("SECTION 1: METADATA (CORRECTED)")
print("=" * 80)

# Library-level metadata (matches corrected cellranger multi output names)
# Tissue is always unambiguous at the library level.
# Timepoint/subspecies are "Mixed" for libraries C, D, F — resolved per-animal.
LIBRARY_META = {
    'SHIV_Pre_PBMC':                {'library': 'A', 'tissue': 'PBMC', 'timepoint': 'Pre',   'subspecies': 'Chinese'},
    'SHIV_Pre_LN':                  {'library': 'B', 'tissue': 'LN',   'timepoint': 'Pre',   'subspecies': 'Chinese'},
    'SHIV_Pre_PBMC_IndianMixed':    {'library': 'C', 'tissue': 'PBMC', 'timepoint': 'Mixed', 'subspecies': 'Mixed'},
    'SHIV_21DPI_Pre_LN_Mixed':      {'library': 'D', 'tissue': 'LN',   'timepoint': 'Mixed', 'subspecies': 'Mixed'},
    'SHIV_21DPI_PBMC':              {'library': 'E', 'tissue': 'PBMC', 'timepoint': '21 DPI','subspecies': 'Chinese'},
    'SHIV_21DPI_NonInf_LN_Mixed':   {'library': 'F', 'tissue': 'LN',   'timepoint': 'Mixed', 'subspecies': 'Mixed'},
    'SHIV_Necropsy_PBMC':           {'library': 'G', 'tissue': 'PBMC', 'timepoint': 'Necropsy','subspecies': 'Chinese'},
    'SHIV_Necropsy_LN':             {'library': 'H', 'tissue': 'LN',   'timepoint': 'Necropsy','subspecies': 'Chinese'},
}

# Per-animal ground truth: (library_letter, animal_id) → metadata
# Verified against provider's sample_table.csv and Zoey's CITEseq Matrix Info
ANIMAL_META = {
    # Library A: Pre PBMC (Chinese)
    ('A', '39272'): {'timepoint': 'Pre',          'subspecies': 'Chinese'},
    ('A', '40702'): {'timepoint': 'Pre',          'subspecies': 'Chinese'},
    ('A', '40707'): {'timepoint': 'Pre',          'subspecies': 'Chinese'},
    ('A', '41861'): {'timepoint': 'Pre',          'subspecies': 'Chinese'},
    # Library B: Pre LN (Chinese)
    ('B', '39272'): {'timepoint': 'Pre',          'subspecies': 'Chinese'},
    ('B', '40702'): {'timepoint': 'Pre',          'subspecies': 'Chinese'},
    ('B', '40707'): {'timepoint': 'Pre',          'subspecies': 'Chinese'},
    ('B', '41861'): {'timepoint': 'Pre',          'subspecies': 'Chinese'},
    # Library C: Pre PBMC — MIXED (Indian + Chinese 41862)
    ('C', '34315'): {'timepoint': 'Non-Infected', 'subspecies': 'Indian'},
    ('C', '41862'): {'timepoint': 'Pre',          'subspecies': 'Chinese'},
    ('C', '41903'): {'timepoint': 'Non-Infected', 'subspecies': 'Indian'},
    # Library D: 21DPI + Pre LN — MIXED
    ('D', '39272'): {'timepoint': '21 DPI',       'subspecies': 'Chinese'},
    ('D', '40702'): {'timepoint': '21 DPI',       'subspecies': 'Chinese'},
    ('D', '40707'): {'timepoint': '21 DPI',       'subspecies': 'Chinese'},
    ('D', '41862'): {'timepoint': 'Pre',          'subspecies': 'Chinese'},
    # Library E: 21DPI PBMC (Chinese)
    ('E', '39272'): {'timepoint': '21 DPI',       'subspecies': 'Chinese'},
    ('E', '40702'): {'timepoint': '21 DPI',       'subspecies': 'Chinese'},
    ('E', '40707'): {'timepoint': '21 DPI',       'subspecies': 'Chinese'},
    ('E', '41861'): {'timepoint': '21 DPI',       'subspecies': 'Chinese'},
    # Library F: 21DPI + NonInf LN — MIXED
    ('F', '34315'): {'timepoint': 'Non-Infected', 'subspecies': 'Indian'},
    ('F', '41861'): {'timepoint': '21 DPI',       'subspecies': 'Chinese'},
    ('F', '41862'): {'timepoint': '21 DPI',       'subspecies': 'Chinese'},
    ('F', '41903'): {'timepoint': 'Non-Infected', 'subspecies': 'Indian'},
    # Library G: Necropsy PBMC (Chinese)
    ('G', '40702'): {'timepoint': 'Necropsy',     'subspecies': 'Chinese'},
    ('G', '40707'): {'timepoint': 'Necropsy',     'subspecies': 'Chinese'},
    ('G', '41861'): {'timepoint': 'Necropsy',     'subspecies': 'Chinese'},
    ('G', '41862'): {'timepoint': 'Necropsy',     'subspecies': 'Chinese'},
    # Library H: Necropsy LN (Chinese)
    ('H', '40702'): {'timepoint': 'Necropsy',     'subspecies': 'Chinese'},
    ('H', '40707'): {'timepoint': 'Necropsy',     'subspecies': 'Chinese'},
    ('H', '41861'): {'timepoint': 'Necropsy',     'subspecies': 'Chinese'},
    ('H', '41862'): {'timepoint': 'Necropsy',     'subspecies': 'Chinese'},
}

# Animal-level constants
ANIMAL_INFO = {
    '34315': {'subspecies': 'Indian',  'infection': 'Non-Infected', 'note': 'control'},
    '41903': {'subspecies': 'Indian',  'infection': 'Non-Infected', 'note': 'control'},
    '39272': {'subspecies': 'Chinese', 'infection': 'Infected',     'note': 'no necropsy'},
    '40702': {'subspecies': 'Chinese', 'infection': 'Infected',     'note': ''},
    '40707': {'subspecies': 'Chinese', 'infection': 'Infected',     'note': 'elite controller'},
    '41861': {'subspecies': 'Chinese', 'infection': 'Infected',     'note': ''},
    '41862': {'subspecies': 'Chinese', 'infection': 'Infected',     'note': 'pooled separately at Pre/21DPI'},
}

print("Library-level metadata:")
for cond, info in LIBRARY_META.items():
    mixed = " [MIXED]" if info['timepoint'] == 'Mixed' else ""
    print(f"  {info['library']}: {cond} — {info['tissue']}, {info['timepoint']}{mixed}")

print(f"\nPer-animal entries: {len(ANIMAL_META)}")
print(f"Animals: {sorted(ANIMAL_INFO.keys())}")

# %% ============================================================
# 2. READ CELL RANGER MULTI OUTPUTS + BARCODE-TO-ANIMAL MAPPING
#
#    Step 1: Read our cellranger multi h5 files (GEX + ADT)
#    Step 2: Load provider's per-animal barcode assignments
#    Step 3: Map barcodes → animal_id → per-animal metadata
# ============================================================
print("\n" + "=" * 80)
print("SECTION 2: DATA INGESTION + BARCODE-TO-ANIMAL MAPPING")
print("=" * 80)

# --- Step 2a: Build barcode-to-animal lookup from provider's outputs ---
print("Building barcode-to-animal lookup from provider's demultiplexed outputs...")

sample_df = pd.read_csv(SAMPLE_TABLE)
print(f"  Provider sample table: {len(sample_df)} entries")

# Map provider run_letter to our library letter
# Provider's GEX run letters match our library letters (A-H)
barcode_to_animal = {}  # {(library_letter, barcode): animal_id}

for _, row in sample_df.iterrows():
    lib = row['run_letter']
    animal = str(row['animal_id'])
    matrix_path = row['matrix_path']

    if pd.isna(matrix_path) or not os.path.exists(str(matrix_path)):
        print(f"    SKIP {animal}_{lib}: matrix not found at {matrix_path}")
        continue

    # Read the provider's barcodes for this animal
    # The h5 file contains the filtered barcodes
    try:
        provider_adata = sc.read_10x_h5(matrix_path, gex_only=True)
        for bc in provider_adata.obs_names:
            # Strip any existing suffix (e.g., -1) for matching
            bc_base = bc.split('-')[0] if '-' in bc else bc
            barcode_to_animal[(lib, bc_base)] = animal
        print(f"    {animal}_{lib}: {provider_adata.n_obs} barcodes loaded")
        del provider_adata
    except Exception as e:
        print(f"    FAIL {animal}_{lib}: {e}")

print(f"  Total barcode-to-animal entries: {len(barcode_to_animal):,}")


# --- Step 2b: Read Cell Ranger multi outputs ---
def find_matrix_h5(condition_name):
    """Locate the filtered feature-barcode matrix h5 for a condition."""
    per_sample_dir = os.path.join(
        ANALYSIS_DIR, condition_name, 'outs', 'per_sample_outs', condition_name
    )
    h5_path = os.path.join(per_sample_dir, 'sample_filtered_feature_bc_matrix.h5')
    if os.path.exists(h5_path):
        return h5_path

    multi_h5 = os.path.join(
        ANALYSIS_DIR, condition_name, 'outs', 'multi', 'count',
        'raw_feature_bc_matrix.h5'
    )
    if os.path.exists(multi_h5):
        return multi_h5
    return None


def read_multimodal_h5(h5_path):
    """Read a Cell Ranger multi h5 file containing GEX + ADT features."""
    adata_full = sc.read_10x_h5(h5_path, gex_only=False)

    if 'feature_types' in adata_full.var.columns:
        ft_col = 'feature_types'
    elif 'feature_type' in adata_full.var.columns:
        ft_col = 'feature_type'
    else:
        return adata_full, None

    gex_mask = adata_full.var[ft_col] == 'Gene Expression'
    adt_mask = adata_full.var[ft_col] == 'Antibody Capture'

    adata_gex = adata_full[:, gex_mask].copy()
    adata_adt = adata_full[:, adt_mask].copy() if adt_mask.sum() > 0 else None

    return adata_gex, adata_adt


def assign_animal_metadata(adata_obj, library_letter, barcode_lookup, animal_meta):
    """Assign per-cell animal_id, timepoint, subspecies using provider barcodes."""
    animal_ids = []
    timepoints = []
    subspecies_list = []

    for bc_full in adata_obj.obs_names:
        # Extract base barcode (before our -Library suffix, and before -1 GEM well suffix)
        # Our barcodes are formatted as: ACGTACGT-1-A (barcode-gemwell-library)
        parts = bc_full.rsplit('-', 1)  # split off library letter
        bc_with_gem = parts[0] if len(parts) > 1 else bc_full
        bc_base = bc_with_gem.split('-')[0]  # strip -1 GEM well

        animal = barcode_lookup.get((library_letter, bc_base), None)

        if animal is not None:
            animal_ids.append(animal)
            meta = animal_meta.get((library_letter, animal), {})
            timepoints.append(meta.get('timepoint', 'Unknown'))
            subspecies_list.append(meta.get('subspecies', 'Unknown'))
        else:
            animal_ids.append('Unknown')
            timepoints.append('Unknown')
            subspecies_list.append('Unknown')

    adata_obj.obs['animal_id'] = animal_ids
    adata_obj.obs['timepoint'] = timepoints
    adata_obj.obs['subspecies'] = subspecies_list

    n_mapped = sum(1 for a in animal_ids if a != 'Unknown')
    n_total = len(animal_ids)
    pct = 100 * n_mapped / n_total if n_total > 0 else 0
    print(f"    Barcode mapping: {n_mapped:,}/{n_total:,} ({pct:.1f}%) assigned to animals")

    if n_mapped < n_total:
        n_unknown = n_total - n_mapped
        print(f"    {n_unknown:,} cells unmatched (ambient/multiplets/threshold diffs)")


# Read and process all libraries
print("\nReading Cell Ranger multi outputs...")
adata_gex_list = []
adata_adt_list = []
sample_summary = []

for condition, info in LIBRARY_META.items():
    h5_path = find_matrix_h5(condition)

    if h5_path is None:
        print(f"  SKIP {condition} (Library {info['library']}): no matrix found")
        sample_summary.append({
            'condition': condition, 'library': info['library'],
            'status': 'MISSING', 'n_cells': 0,
        })
        continue

    print(f"  Reading {condition} (Library {info['library']})...")

    try:
        adata_gex, adata_adt = read_multimodal_h5(h5_path)

        # Attach library-level metadata
        adata_gex.obs['condition'] = condition
        adata_gex.obs['library'] = info['library']
        adata_gex.obs['tissue'] = info['tissue']

        # Make barcodes globally unique: barcode-library
        adata_gex.obs_names = [
            f"{bc}-{info['library']}" for bc in adata_gex.obs_names
        ]

        # Assign per-cell animal metadata using provider barcodes
        assign_animal_metadata(adata_gex, info['library'],
                               barcode_to_animal, ANIMAL_META)

        n_cells = adata_gex.n_obs
        print(f"    GEX: {n_cells:,} cells × {adata_gex.n_vars:,} genes")

        # Report animal breakdown
        animal_counts = adata_gex.obs['animal_id'].value_counts()
        for animal, count in animal_counts.items():
            tp = adata_gex.obs.loc[adata_gex.obs['animal_id'] == animal, 'timepoint'].iloc[0]
            sp = adata_gex.obs.loc[adata_gex.obs['animal_id'] == animal, 'subspecies'].iloc[0]
            note = ANIMAL_INFO.get(animal, {}).get('note', '')
            note_str = f" ({note})" if note else ""
            print(f"      {animal}: {count:,} cells — {tp}, {sp}{note_str}")

        adata_gex_list.append(adata_gex)

        # Handle ADT
        if adata_adt is not None:
            adata_adt.obs_names = [
                f"{bc}-{info['library']}" for bc in adata_adt.obs_names
            ]
            adata_adt_list.append(adata_adt)
            print(f"    ADT: {adata_adt.n_vars} proteins")

        sample_summary.append({
            'condition': condition, 'library': info['library'],
            'status': 'OK', 'n_cells': n_cells,
        })

    except Exception as e:
        print(f"    FAILED: {e}")
        sample_summary.append({
            'condition': condition, 'library': info['library'],
            'status': f'ERROR: {e}', 'n_cells': 0,
        })

# Concatenate GEX
if not adata_gex_list:
    raise RuntimeError("No valid GEX samples found!")

adata = ad.concat(adata_gex_list, join='outer', merge='same')
print(f"\nCombined GEX: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

# Create a clean sample_name: animal_tissue_timepoint
adata.obs['sample_name'] = (
    adata.obs['animal_id'] + '_' +
    adata.obs['tissue'] + '_' +
    adata.obs['timepoint']
)

# Concatenate ADT
if adata_adt_list:
    adata_adt_combined = ad.concat(adata_adt_list, join='outer', merge='same')
    common_bcs = adata.obs_names.intersection(adata_adt_combined.obs_names)
    adt_aligned = adata_adt_combined[common_bcs].copy()
    adt_df = pd.DataFrame(
        adt_aligned.X.toarray() if sparse.issparse(adt_aligned.X) else adt_aligned.X,
        index=adt_aligned.obs_names,
        columns=adt_aligned.var_names,
    )
    adt_df = adt_df.reindex(adata.obs_names, fill_value=0)
    adata.obsm['protein_counts'] = adt_df.values
    adata.uns['adt_names'] = list(adt_df.columns)
    print(f"ADT stored: {adt_df.shape[1]} proteins")

# Summary
print(f"\n=== Ingestion Summary ===")
print(pd.DataFrame(sample_summary).to_string(index=False))
print(f"\nTotal cells: {adata.n_obs:,}")
print(f"\nCells per animal:")
print(adata.obs['animal_id'].value_counts().to_string())
print(f"\nCells per timepoint:")
print(adata.obs['timepoint'].value_counts().to_string())
print(f"\nCells per tissue × timepoint:")
print(pd.crosstab(adata.obs['tissue'], adata.obs['timepoint']))

# %% ============================================================
# 3. GENE CATEGORY TAGGING (mito, ribo, hb, SHIV)
# ============================================================
print("\n" + "=" * 80)
print("SECTION 3: GENE CATEGORY TAGGING")
print("=" * 80)

adata.var['mt'] = adata.var_names.str.startswith('KEG06_')
print(f"  Mitochondrial genes (KEG06_*): {adata.var['mt'].sum()}")

adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL'))
print(f"  Ribosomal genes (RPS/RPL): {adata.var['ribo'].sum()}")

adata.var['hb'] = adata.var_names.str.contains('^HB[^(P)]', case=False)
print(f"  Hemoglobin genes: {adata.var['hb'].sum()}")

SHIV_GENES = ['gag', 'pol', 'vif', 'vpx', 'vpr', 'tat', 'rev', 'vpu', 'env', 'nef']
shiv_in_data = [g for g in SHIV_GENES if g in adata.var_names]
adata.var['shiv'] = adata.var_names.isin(SHIV_GENES)
n_shiv = adata.var['shiv'].sum()
print(f"  SHIV genes: {n_shiv}/{len(SHIV_GENES)} — {shiv_in_data}")

if n_shiv == 0:
    shiv_lower = {g.lower() for g in SHIV_GENES}
    shiv_matches = [g for g in adata.var_names if g.lower() in shiv_lower]
    if shiv_matches:
        adata.var['shiv'] = adata.var_names.isin(shiv_matches)
        n_shiv = adata.var['shiv'].sum()
        print(f"  (Case-insensitive match: {n_shiv}: {shiv_matches})")

sc.pp.calculate_qc_metrics(
    adata, qc_vars=['mt', 'ribo', 'hb', 'shiv'],
    inplace=True, percent_top=[20], log1p=True
)

# SHIV+ summary — now with correct per-animal metadata
if n_shiv > 0:
    shiv_pos = adata.obs['total_counts_shiv'] > 0
    print(f"\n  SHIV+ cells: {shiv_pos.sum():,}/{adata.n_obs:,} ({100*shiv_pos.mean():.2f}%)")

    print("  SHIV+ by timepoint:")
    for tp in sorted(adata.obs['timepoint'].unique()):
        mask = adata.obs['timepoint'] == tp
        n_pos = (shiv_pos & mask).sum()
        n_total = mask.sum()
        pct = 100 * n_pos / n_total if n_total > 0 else 0
        print(f"    {tp}: {n_pos}/{n_total} ({pct:.3f}%)")

    print("  SHIV+ by animal:")
    for animal in sorted(adata.obs['animal_id'].unique()):
        if animal == 'Unknown':
            continue
        mask = adata.obs['animal_id'] == animal
        n_pos = (shiv_pos & mask).sum()
        n_total = mask.sum()
        pct = 100 * n_pos / n_total if n_total > 0 else 0
        note = ANIMAL_INFO.get(animal, {}).get('note', '')
        note_str = f" ({note})" if note else ""
        if n_pos > 0 or animal in ['40707']:
            print(f"    {animal}{note_str}: {n_pos}/{n_total} ({pct:.3f}%)")

# %%  ============================================================
# SECTIONS 4-12: Unchanged from v1 (QC through Export)
# The metadata is now correct at the per-cell level, so all
# downstream groupby operations (by timepoint, animal, tissue)
# will produce correct results.
# ============================================================

# %% ============================================================
# 4. QUALITY CONTROL (MAD-based + Scrublet)
# ============================================================
print("\n" + "=" * 80)
print("SECTION 4: QUALITY CONTROL")
print("=" * 80)

def is_outlier(adata_obj, metric, nmads):
    M = adata_obj.obs[metric]
    median_val = M.median()
    mad = median_abs_deviation(M, nan_policy='omit')
    if mad == 0:
        mad = M.std()
        if mad == 0:
            return pd.Series([False] * len(M), index=M.index)
    outlier_mask = (np.abs(M - median_val) > nmads * mad)
    n_out = outlier_mask.sum()
    print(f"  {metric}: {n_out} outliers ({100*n_out/len(M):.1f}%) "
          f"[median={median_val:.2f}, MAD={mad:.2f}]")
    return outlier_mask

sc.pl.scatter(adata, 'total_counts', 'n_genes_by_counts',
              color='pct_counts_mt', save='_scatter_preQC.pdf')
plt.close('all')

adata.obs['outlier'] = (
    is_outlier(adata, 'log1p_total_counts', MAD_COUNTS_GENES)
    | is_outlier(adata, 'log1p_n_genes_by_counts', MAD_COUNTS_GENES)
    | is_outlier(adata, 'pct_counts_in_top_20_genes', MAD_COUNTS_GENES)
)
print(f"General outlier cells: {adata.obs['outlier'].sum():,}")

adata.obs['mt_outlier'] = (
    is_outlier(adata, 'pct_counts_mt', MAD_MITO)
    | (adata.obs['pct_counts_mt'] > MT_THRESHOLD)
)
print(f"Mito outlier cells: {adata.obs['mt_outlier'].sum():,}")

n_before = adata.n_obs
adata = adata[(~adata.obs['outlier']) & (~adata.obs['mt_outlier'])].copy()
print(f"\nCells removed by outlier/mito filter: {n_before - adata.n_obs:,}")

sc.pp.filter_cells(adata, min_genes=MIN_GENES)
sc.pp.filter_genes(adata, min_cells=MIN_CELLS_PER_GENE)
print(f"After min_genes={MIN_GENES}, min_cells={MIN_CELLS_PER_GENE}:")
print(f"  {adata.n_obs:,} cells × {adata.n_vars:,} genes")

sc.pl.scatter(adata, 'total_counts', 'n_genes_by_counts',
              color='pct_counts_mt', save='_scatter_postQC.pdf')
plt.close('all')

print("\nCells per animal after QC:")
print(adata.obs['animal_id'].value_counts().to_string())
print(f"\nCells per timepoint:")
print(adata.obs['timepoint'].value_counts().to_string())

# %% ============================================================
# 4b. DOUBLET REMOVAL (Scrublet, per-library)
# ============================================================
print("\n--- Doublet Detection (Scrublet, per-library) ---")

try:
    import scrublet as scr

    n_before = adata.n_obs
    doublet_scores = np.zeros(adata.n_obs)
    predicted_doublets = np.zeros(adata.n_obs, dtype=bool)

    # Run per-library (not per-condition) since that's the GEM well level
    for lib in sorted(adata.obs['library'].unique()):
        lib_mask = adata.obs['library'] == lib
        idx = np.where(lib_mask)[0]
        lib_adata = adata[lib_mask].copy()

        if lib_adata.n_obs < 100:
            print(f"  Library {lib}: <100 cells, skip")
            continue

        try:
            scrub = scr.Scrublet(lib_adata.X,
                                 expected_doublet_rate=EXPECTED_DOUBLET_RATE)
            scores, predictions = scrub.scrub_doublets(
                min_counts=2, min_cells=3, verbose=False)
            doublet_scores[idx] = scores
            predicted_doublets[idx] = predictions
            n_d = predictions.sum()
            print(f"  Library {lib}: {n_d} doublets ({100*n_d/len(predictions):.1f}%)")
        except Exception as e:
            print(f"  Library {lib}: Scrublet failed — {e}")

    adata.obs['doublet_score'] = doublet_scores
    adata.obs['predicted_doublet'] = predicted_doublets
    adata = adata[~adata.obs['predicted_doublet']].copy()
    print(f"Doublets removed: {n_before - adata.n_obs:,}")
    print(f"Cells after doublet removal: {adata.n_obs:,}")

except ImportError:
    print("Scrublet not installed — skipping.")
    adata.obs['doublet_score'] = 0.0
    adata.obs['predicted_doublet'] = False

# %% Save QC'd object
adata.write(os.path.join(WORKING_DIR, 'adata_qc.h5ad'))
print(f"\nSaved: adata_qc.h5ad ({adata.n_obs:,} cells)")

# %% ============================================================
# 5. NORMALIZE, EMBED, BATCH CORRECT
# ============================================================
print("\n" + "=" * 80)
print("SECTION 5: NORMALIZATION AND EMBEDDING")
print("=" * 80)

adata.layers['counts'] = adata.X.copy()

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

exclude_mask = adata.var['shiv'] | adata.var['mt']
adata_hvg = adata[:, ~exclude_mask].copy()
sc.pp.highly_variable_genes(adata_hvg, n_top_genes=3000, flavor='seurat_v3',
                            layer='counts')
hvg_genes = adata_hvg.var_names[adata_hvg.var['highly_variable']].tolist()

adata.var['highly_variable'] = False
adata.var.loc[adata.var_names.isin(hvg_genes), 'highly_variable'] = True
print(f"HVGs selected: {adata.var['highly_variable'].sum()}")

sc.tl.pca(adata, n_comps=50, use_highly_variable=True)
sc.pp.neighbors(adata, n_pcs=N_PCS)
sc.tl.umap(adata)
sc.tl.leiden(adata, key_added='clusters', resolution=LEIDEN_RES, random_state=42)
print(f"Clusters (pre-batch correction): {adata.obs['clusters'].nunique()}")

fig, axes = plt.subplots(1, 4, figsize=(26, 5))
sc.pl.umap(adata, color='library',    size=3, ax=axes[0], show=False, title='Library')
sc.pl.umap(adata, color='clusters',   size=3, ax=axes[1], show=False, title='Clusters')
sc.pl.umap(adata, color='tissue',     size=3, ax=axes[2], show=False, title='Tissue')
sc.pl.umap(adata, color='timepoint',  size=3, ax=axes[3], show=False, title='Timepoint')
plt.tight_layout()
fig.savefig(os.path.join(FIGURES_DIR, 'umap_before_batch_correction.pdf'),
            bbox_inches='tight', dpi=150)
plt.close('all')

try:
    import bbknn
    # Batch correct by library (GEM well level — the source of technical variation)
    print(f"\nRunning BBKNN (batch_key='library', {adata.obs['library'].nunique()} batches)...")
    bbknn.bbknn(adata, batch_key='library', n_pcs=N_PCS)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, key_added='clusters', resolution=LEIDEN_RES, random_state=42)
    print(f"Clusters after BBKNN: {adata.obs['clusters'].nunique()}")

    fig, axes = plt.subplots(1, 4, figsize=(26, 5))
    sc.pl.umap(adata, color='library',    size=3, ax=axes[0], show=False, title='Library')
    sc.pl.umap(adata, color='clusters',   size=3, ax=axes[1], show=False, title='Clusters')
    sc.pl.umap(adata, color='tissue',     size=3, ax=axes[2], show=False, title='Tissue')
    sc.pl.umap(adata, color='timepoint',  size=3, ax=axes[3], show=False, title='Timepoint')
    plt.tight_layout()
    fig.savefig(os.path.join(FIGURES_DIR, 'umap_after_batch_correction.pdf'),
                bbox_inches='tight', dpi=150)
    plt.close('all')
except ImportError:
    print("bbknn not installed — skipping batch correction.")

# %% CLR Normalize CITE-seq
if 'protein_counts' in adata.obsm:
    print("\n--- CLR Normalizing CITE-seq ADT ---")
    adt_raw = adata.obsm['protein_counts'].copy()
    adt_log = np.log1p(adt_raw)
    geometric_mean = np.exp(np.mean(adt_log, axis=1, keepdims=True))
    adt_clr = adt_log - np.log(geometric_mean)
    adata.obsm['protein_clr'] = adt_clr
    print(f"  CLR-normalized ADT: {adt_clr.shape}")

adata.write(os.path.join(WORKING_DIR, 'adata_pp.h5ad'))
print(f"\nSaved: adata_pp.h5ad ({adata.n_obs:,} cells)")

# %% ============================================================
# 6. TIER 1: BROAD LINEAGE ANNOTATION
# ============================================================
print("\n" + "=" * 80)
print("SECTION 6: TIER 1 — BROAD LINEAGE ANNOTATION")
print("=" * 80)

TIER1_MARKERS = {
    'T cells':     {'markers': ['CD3D', 'CD3E', 'CD3G']},
    'B cells':     {'markers': ['CD79A', 'MS4A1', 'CD19', 'BANK1', 'PAX5']},
    'NK cells':    {'markers': ['KLRD1', 'KLRC3', 'GNLY', 'PRF1', 'GZMB', 'NKG7']},
    'Monocytes':   {'markers': ['S100A8', 'S100A9', 'VCAN', 'CD14', 'LYZ']},
    'Macrophages': {'markers': ['CD68', 'CD163', 'MRC1', 'MERTK', 'APOC1']},
    'DCs':         {'markers': ['FLT3', 'CLEC9A', 'CLEC10A', 'IL3RA', 'CLEC4C']},
    'Mast cells':  {'markers': ['CPA3', 'KIT', 'HPGD', 'ENPP3']},
    'Erythrocytes':{'markers': ['HBB', 'HBA1', 'HBA2']},
    'Platelets':   {'markers': ['PPBP', 'PF4', 'GP9']},
}

score_cols = {}
for lineage, info in TIER1_MARKERS.items():
    present = [g for g in info['markers'] if g in adata.var_names]
    score_key = f"tier1_{lineage.replace(' ', '_')}"

    if len(present) >= 2:
        sc.tl.score_genes(adata, gene_list=present, score_name=score_key)
        score_cols[lineage] = score_key
        print(f"  {lineage}: {len(present)}/{len(info['markers'])} markers, "
              f"mean={adata.obs[score_key].mean():.3f}")
    elif len(present) == 1:
        expr = adata[:, present[0]].X
        if sparse.issparse(expr): expr = expr.toarray().flatten()
        else: expr = np.array(expr).flatten()
        adata.obs[score_key] = expr
        score_cols[lineage] = score_key
        print(f"  {lineage}: 1 marker ({present[0]})")
    else:
        print(f"  SKIP {lineage}: no markers")

lineage_names = list(score_cols.keys())
lineage_col_list = [score_cols[n] for n in lineage_names]
score_matrix = adata.obs[lineage_col_list].values
best_idx = np.argmax(score_matrix, axis=1)
best_score = np.max(score_matrix, axis=1)

adata.obs['tier1_celltype'] = pd.Categorical(np.where(
    best_score >= SCORE_THRESHOLD,
    [lineage_names[i] for i in best_idx],
    'Unassigned'
))
adata.obs['tier1_score'] = best_score

print(f"\n  Tier 1 assignment:")
print(adata.obs['tier1_celltype'].value_counts().to_string())

# %% ============================================================
# 7. TIER 2: SUBTYPE ANNOTATION
# ============================================================
print("\n" + "=" * 80)
print("SECTION 7: TIER 2 — SUBTYPE ANNOTATION")
print("=" * 80)

TIER2_MARKERS = {
    'T cells': {
        'CD4': ['CD4', 'BCL11B', 'LTB', 'IL7R'],
        'CD8': ['CD8A', 'CD8B', 'GZMK'],
    },
    'NK cells': {
        'NK CD56bright': ['IL7R', 'SELL', 'XCL1'],
        'NK CD56dim':    ['FCGR3A', 'FGFBP2', 'CX3CR1'],
    },
    'B cells': {
        'B naive':  ['MS4A1', 'IGHD', 'TCL1A'],
        'B memory': ['CD27', 'AIM2'],
        'Plasma':   ['JCHAIN', 'MZB1', 'SDC1', 'XBP1'],
    },
    'Monocytes': {
        'Classical':     ['CD14', 'LYZ', 'S100A8', 'S100A9', 'VCAN'],
        'Non-classical': ['FCGR3A', 'CX3CR1', 'ITGAX'],
    },
    'Macrophages': {
        'Tissue-resident':  ['CD68', 'APOC1', 'FABP4', 'C1QA', 'C1QB'],
        'CD169+':           ['SIGLEC1', 'FMNL2'],
        'Inflammatory':     ['CD163', 'CCL2', 'CCL7', 'CXCL3'],
        'IFN-responsive':   ['IDO1', 'CXCL9', 'CXCL10', 'CXCL11', 'CD274'],
    },
    'DCs': {
        'cDC1': ['CLEC9A', 'XCR1', 'BATF3'],
        'cDC2': ['CD1C', 'CLEC10A', 'FCER1A'],
        'pDC':  ['IL3RA', 'CLEC4C', 'TCF4', 'PACSIN1', 'COBLL1'],
    },
}

tier2_score_cols = {}
for lineage, subtypes in TIER2_MARKERS.items():
    tier2_score_cols[lineage] = {}
    print(f"\n  --- {lineage} ---")
    for subtype, markers in subtypes.items():
        present = [g for g in markers if g in adata.var_names]
        score_key = (f"tier2_{lineage.replace(' ', '_')}_"
                     f"{subtype.replace(' ', '_').replace('+', 'pos')}")
        if len(present) >= 2:
            sc.tl.score_genes(adata, gene_list=present, score_name=score_key)
            tier2_score_cols[lineage][subtype] = score_key
            print(f"    {subtype}: {len(present)} markers")
        elif len(present) == 1:
            expr = adata[:, present[0]].X
            if sparse.issparse(expr): expr = expr.toarray().flatten()
            else: expr = np.array(expr).flatten()
            adata.obs[score_key] = expr
            tier2_score_cols[lineage][subtype] = score_key
            print(f"    {subtype}: 1 marker ({present[0]})")
        else:
            print(f"    SKIP {subtype}")

tier2_labels = []
tier2_scores = []
for i in range(adata.n_obs):
    t1 = adata.obs['tier1_celltype'].iloc[i]
    if t1 in tier2_score_cols and tier2_score_cols[t1]:
        subtypes = tier2_score_cols[t1]
        snames = list(subtypes.keys())
        svals = [adata.obs[subtypes[s]].iloc[i] for s in snames]
        best = np.argmax(svals)
        tier2_labels.append(snames[best])
        tier2_scores.append(svals[best])
    elif t1 == 'Unassigned':
        tier2_labels.append('Unassigned')
        tier2_scores.append(0.0)
    else:
        tier2_labels.append(t1)
        tier2_scores.append(adata.obs['tier1_score'].iloc[i])

adata.obs['tier2_celltype'] = pd.Categorical(tier2_labels)
adata.obs['tier2_score'] = tier2_scores
adata.obs['tier2_full'] = [
    f"{t1}: {t2}" if t1 != t2 and t1 != 'Unassigned' else t1
    for t1, t2 in zip(adata.obs['tier1_celltype'], adata.obs['tier2_celltype'])
]

print(f"\n  Tier 2 assignment:")
print(adata.obs['tier2_full'].value_counts().to_string())

# %% ============================================================
# 8. CITE-seq PROTEIN VALIDATION
# ============================================================
print("\n" + "=" * 80)
print("SECTION 8: CITE-seq PROTEIN VALIDATION")
print("=" * 80)

if 'protein_clr' in adata.obsm and 'adt_names' in adata.uns:
    adt_names = adata.uns['adt_names']
    adt_clr = pd.DataFrame(adata.obsm['protein_clr'], index=adata.obs_names, columns=adt_names)
    for protein in adt_names:
        adata.obs[f'ADT_{protein}'] = adt_clr[protein].values
    print(f"  ADT proteins: {adt_names}")

    mac_adts = [c for c in adata.obs.columns if c.startswith('ADT_ADT_CD')
                and any(m in c for m in ['CD68', 'CD163', 'CD169', 'CD206', 'CD11b'])]
    if mac_adts:
        fig, axes = plt.subplots(1, len(mac_adts), figsize=(5*len(mac_adts), 4))
        if len(mac_adts) == 1: axes = [axes]
        for i, col in enumerate(mac_adts):
            sc.pl.umap(adata, color=col, ax=axes[i], show=False, cmap='Reds',
                       frameon=False, title=col.replace('ADT_ADT_', 'Protein: '))
        plt.tight_layout()
        fig.savefig(os.path.join(FIGURES_DIR, 'umap_adt_macrophage_markers.pdf'),
                    bbox_inches='tight', dpi=150)
        plt.close('all')
        print("  Saved: umap_adt_macrophage_markers.pdf")

    # ADT-assisted macrophage rescue
    unassigned_mask = adata.obs['tier1_celltype'] == 'Unassigned'
    if unassigned_mask.sum() > 0 and 'ADT_ADT_CD68' in adata.obs.columns:
        mac_cells = adata.obs['tier1_celltype'] == 'Macrophages'
        cd68_thresh = (adata.obs.loc[mac_cells, 'ADT_ADT_CD68'].quantile(0.25)
                       if mac_cells.sum() > 10 else 0.5)
        rescued = unassigned_mask & (adata.obs['ADT_ADT_CD68'] > cd68_thresh)
        if rescued.sum() > 0:
            print(f"\n  ADT macrophage rescue: {rescued.sum():,} cells")
            adata.obs.loc[rescued, 'tier1_celltype'] = 'Macrophages'
            adata.obs.loc[rescued, 'tier2_full'] = 'Macrophages: ADT-rescued'

    mac_mask = adata.obs['tier1_celltype'] == 'Macrophages'
    print(f"\n  Macrophages total: {mac_mask.sum():,} ({100*mac_mask.mean():.1f}%)")
    for tissue in ['PBMC', 'LN']:
        t_mask = adata.obs['tissue'] == tissue
        n_mac = (mac_mask & t_mask).sum()
        pct = 100 * n_mac / t_mask.sum() if t_mask.sum() > 0 else 0
        print(f"    {tissue}: {n_mac:,} ({pct:.1f}%)")

    # ADT dotplot
    adt_for_plot = ad.AnnData(
        X=adata.obsm['protein_clr'],
        obs=adata.obs[['tier1_celltype']].copy(),
        var=pd.DataFrame(index=adt_names),
    )
    sc.pl.dotplot(adt_for_plot, var_names=adt_names, groupby='tier1_celltype',
                  standard_scale='var', save='_adt_by_tier1_celltype.pdf')
    plt.close('all')
    print("  Saved: dotplot_adt_by_tier1_celltype.pdf")
else:
    print("  No CITE-seq data — skipping")

# %% ============================================================
# 9. SHIV+ CELL IDENTIFICATION
# ============================================================
print("\n" + "=" * 80)
print("SECTION 9: SHIV+ CELL IDENTIFICATION")
print("=" * 80)

if adata.var['shiv'].sum() > 0:
    shiv_genes = adata.var_names[adata.var['shiv']].tolist()
    shiv_expr = adata[:, shiv_genes].layers.get('counts', adata[:, shiv_genes].X)
    if sparse.issparse(shiv_expr): shiv_expr = shiv_expr.toarray()

    adata.obs['shiv_total_counts'] = np.sum(shiv_expr, axis=1)
    adata.obs['shiv_n_genes'] = np.sum(shiv_expr > 0, axis=1)
    adata.obs['shiv_positive'] = adata.obs['shiv_total_counts'] > 0

    n_pos = adata.obs['shiv_positive'].sum()
    print(f"  SHIV+ cells: {n_pos:,}/{adata.n_obs:,} ({100*n_pos/adata.n_obs:.2f}%)")

    print(f"\n  Per-gene detection:")
    for gene in shiv_genes:
        g_expr = adata[:, gene].layers.get('counts', adata[:, gene].X)
        if sparse.issparse(g_expr): g_expr = g_expr.toarray().flatten()
        print(f"    {gene}: {(g_expr > 0).sum():,} cells")

    print(f"\n  SHIV+ by cell type:")
    for ct in adata.obs['tier1_celltype'].cat.categories:
        ct_mask = adata.obs['tier1_celltype'] == ct
        n_pos_ct = (adata.obs['shiv_positive'] & ct_mask).sum()
        if n_pos_ct > 0:
            print(f"    {ct}: {n_pos_ct}/{ct_mask.sum()} ({100*n_pos_ct/ct_mask.sum():.2f}%)")

    print(f"\n  SHIV+ by animal × timepoint:")
    for animal in sorted(adata.obs['animal_id'].unique()):
        if animal == 'Unknown': continue
        for tp in sorted(adata.obs['timepoint'].unique()):
            mask = (adata.obs['animal_id'] == animal) & (adata.obs['timepoint'] == tp)
            if mask.sum() == 0: continue
            n_pos_at = (adata.obs['shiv_positive'] & mask).sum()
            if n_pos_at > 0:
                note = ANIMAL_INFO.get(animal, {}).get('note', '')
                note_str = f" ({note})" if note else ""
                print(f"    {animal}{note_str} @ {tp}: {n_pos_at}/{mask.sum()}")

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    sc.pl.umap(adata, color='shiv_positive', ax=axes[0], show=False,
               title='SHIV+ cells', frameon=False)
    sc.pl.umap(adata, color='shiv_total_counts', ax=axes[1], show=False,
               title='SHIV UMI counts', cmap='Reds', frameon=False)
    sc.pl.umap(adata, color='tier1_celltype', ax=axes[2], show=False,
               title='Cell types', frameon=False)
    plt.tight_layout()
    fig.savefig(os.path.join(FIGURES_DIR, 'umap_shiv_detection.pdf'),
                bbox_inches='tight', dpi=150)
    plt.close('all')
    print("  Saved: umap_shiv_detection.pdf")
else:
    adata.obs['shiv_positive'] = False
    adata.obs['shiv_total_counts'] = 0

# %% ============================================================
# 10. BIOLOGICAL GENE SET SCORING
# ============================================================
print("\n" + "=" * 80)
print("SECTION 10: BIOLOGICAL GENE SET SCORING")
print("=" * 80)

def score_gene_set(adata_obj, gene_list, score_name):
    present = [g for g in gene_list if g in adata_obj.var_names]
    if len(present) < 2:
        adata_obj.obs[score_name] = 0.0
        print(f"  SKIP {score_name}: {len(present)}/{len(gene_list)} present")
        return 0
    sc.tl.score_genes(adata_obj, gene_list=present, score_name=score_name)
    print(f"  {score_name}: {len(present)}/{len(gene_list)} genes")
    return len(present)

score_gene_set(adata, ['PDCD1','LAG3','HAVCR2','TIGIT','CTLA4','TOX','ENTPD1'], 'exhaustion_score')
score_gene_set(adata, ['TBX21','TNF','BHLHE40','IFNG'], 'score_Th1')
score_gene_set(adata, ['RORC','STAT3','RORA','IL17A'], 'score_Th17')
score_gene_set(adata, ['PRF1','GZMB','GNLY','IFNG','FCGR3A','KLRK1','NCR1'], 'NK_activation_score')
score_gene_set(adata, ['CD4','CCR5','CXCR4','APOBEC3G','APOBEC3F','TRIM5','BST2','SAMHD1',
                       'IFITM1','IFITM2','IFITM3','MX2'], 'SHIV_host_factor_score')
# Use MAMU nomenclature for rhesus MHC genes
score_gene_set(adata, ['MAMU-DRA','MAMU-DRB1','CD74','CD80','CD86','CIITA',
                       'TAP1','TAP2','B2M'], 'antigen_presentation_score')
score_gene_set(adata, ['CXCL9','CXCL10','CXCL11','IDO1','CD274',
                       'GBP1','GBP2','GBP5','STAT1','IRF1','IRF7'], 'IFN_responsive_score')

# %% ============================================================
# 11. DIAGNOSTICS
# ============================================================
print("\n" + "=" * 80)
print("SECTION 11: DIAGNOSTICS")
print("=" * 80)

fig, axes = plt.subplots(1, 4, figsize=(26, 5))
sc.pl.umap(adata, color='tier1_celltype', ax=axes[0], show=False, title='Tier 1', frameon=False)
sc.pl.umap(adata, color='tissue',         ax=axes[1], show=False, title='Tissue', frameon=False)
sc.pl.umap(adata, color='timepoint',      ax=axes[2], show=False, title='Timepoint', frameon=False)
sc.pl.umap(adata, color='animal_id',      ax=axes[3], show=False, title='Animal', frameon=False)
plt.tight_layout()
fig.savefig(os.path.join(ANNOT_DIR, 'umap_tier1_overview.pdf'), bbox_inches='tight', dpi=300)
plt.close('all')
print("  Saved: umap_tier1_overview.pdf")

fig, ax = plt.subplots(figsize=(12, 8))
sc.pl.umap(adata, color='tier2_full', ax=ax, show=False, title='Tier 2: All Subtypes',
           frameon=False, legend_loc='right margin', legend_fontsize=7)
plt.tight_layout()
fig.savefig(os.path.join(ANNOT_DIR, 'umap_tier2_full.pdf'), bbox_inches='tight', dpi=300)
plt.close('all')
print("  Saved: umap_tier2_full.pdf")

# Tier 1 dotplot
TIER1_DOTPLOT = {
    'T cells': ['CD3D','CD3E'], 'B cells': ['CD79A','MS4A1','CD19'],
    'NK cells': ['KLRD1','GNLY','NKG7','GZMB'], 'Monocytes': ['S100A8','S100A9','VCAN','CD14'],
    'Macrophages': ['CD68','CD163','MRC1','APOC1'], 'DCs': ['FLT3','CLEC9A','CD1C','IL3RA'],
}
t1_groups = {ct: [g for g in genes if g in adata.var_names] for ct, genes in TIER1_DOTPLOT.items()}
t1_groups = {k: v for k, v in t1_groups.items() if v}
t1_order = [t for t in t1_groups if t in adata.obs['tier1_celltype'].values]

if t1_groups and t1_order:
    sc.pl.dotplot(adata[adata.obs['tier1_celltype'].isin(t1_order)],
                  var_names=t1_groups, groupby='tier1_celltype',
                  categories_order=t1_order, standard_scale='var',
                  save='_tier1_markers.pdf')
    plt.close('all')
    print("  Saved: dotplot_tier1_markers.pdf")

# Annotation summary
print("\n--- Annotation Summary ---")
for lineage in TIER1_MARKERS:
    n = (adata.obs['tier1_celltype'] == lineage).sum()
    pct = 100 * n / adata.n_obs
    flag = " *** LOW" if lineage == 'Macrophages' and pct < 0.5 else ""
    print(f"  {lineage}: {n:>8,} ({pct:5.1f}%){flag}")

n_un = (adata.obs['tier1_celltype'] == 'Unassigned').sum()
print(f"  Unassigned: {n_un:,} ({100*n_un/adata.n_obs:.1f}%)")

# Proportions by animal × timepoint × tissue
props = pd.crosstab(
    [adata.obs['animal_id'], adata.obs['tissue'], adata.obs['timepoint']],
    adata.obs['tier1_celltype'], normalize='index'
) * 100
props.to_csv(os.path.join(ANNOT_DIR, 'tier1_proportions_by_animal.csv'))
print(f"  Proportions saved: tier1_proportions_by_animal.csv")

# %% ============================================================
# 12. EXPORT
# ============================================================
print("\n" + "=" * 80)
print("SECTION 12: EXPORT")
print("=" * 80)

adata.write(os.path.join(WORKING_DIR, 'adata_annotated.h5ad'))
print(f"  Saved: adata_annotated.h5ad")
print(f"  Shape: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
print(f"  Animals: {adata.obs['animal_id'].nunique()}")
print(f"  Tier 1: {adata.obs['tier1_celltype'].nunique()} lineages")
print(f"  Tier 2: {adata.obs['tier2_full'].nunique()} subtypes")
print(f"  SHIV+ cells: {adata.obs['shiv_positive'].sum():,}")

tier_summary = adata.obs['tier2_full'].value_counts().to_frame('n_cells')
tier_summary['pct'] = (tier_summary['n_cells'] / adata.n_obs * 100).round(1)
tier_summary.to_csv(os.path.join(ANNOT_DIR, 'tier2_counts.csv'))

print(f"\n{'=' * 80}")
print("PIPELINE COMPLETE")
print(f"{'=' * 80}")
print(f"  Next steps:")
print(f"    - VDJ integration (scirpy)")
print(f"    - SHIV barcode extraction from BAMs")
print(f"    - Elite controller deep dive (40707)")
print(f"    - Viral dynamics: Pre → 21 DPI → Necropsy by cell type × animal")
