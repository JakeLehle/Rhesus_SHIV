#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
04 — Ingest the v11 Cell Ranger outputs, assign metadata, run QC
=================================================================

Derived from scanpy_analysis_script.py sections 0 to 4. Those are the only
sections of that 1,023-line script that feed the production object: sections 5
to 12 re-embed and re-annotate, and the real chain discards all of it by
branching off adata_qc.h5ad and redoing both in shiv_embeding.py and the
Tier 1 correction scripts.

Reads   pipeline_results/03_cellranger_multi/SHIV_<COND>/outs/...
Writes  pipeline_results/04_ingest_qc/
            adata_prefilter.h5ad    EVERY called cell, before any filtering
            adata_qc.h5ad           after QC, the object the chain continues on
            qc_dropped_ledger.csv   every dropped barcode and why
            qc_summary.csv, shiv_whitelist.csv


THE SHIV TAG IS A HARD STOP, NOT WARN-AND-CONTINUE
---------------------------------------------------
The original tagged the virus like this:

    SHIV_GENES = ['gag','pol','vif','vpx','vpr','tat','rev','vpu','env','nef']
    adata.var['shiv'] = adata.var_names.isin(SHIV_GENES)
    if n_shiv == 0:                      # case-insensitive fallback
        shiv_lower = {g.lower() for g in SHIV_GENES}
        ...

Under the v11 collapsed GTF none of those ten names exist, so n_shiv is 0 and
the fallback fires. It then matches host **TAT**, tyrosine aminotransferase,
which check_all_libraries_shiv.py already documented as colliding with viral
*tat*. The run would NOT error. It would report per-cell, per-animal and
per-timepoint SHIV UMI that are actually tyrosine aminotransferase expression,
and the numbers would look entirely plausible.

So: exact, case-sensitive match on the single collapsed feature, no fallback,
and the script STOPS if it is absent. Warn-and-continue is the convention on
this project and this is a deliberate exception, because continuing produces
a confidently wrong answer rather than a missing one.


QC, AND WHAT CHANGED ON 2026-11-02
-----------------------------------
    mitochondrial   pct_counts_mt > 20%            fixed threshold, no MAD
    counts/genes    5-MAD, two-sided, global       unchanged
    complexity      5-MAD on pct_counts_in_top_20  unchanged
    floor           min_genes 200, min_cells 20    unchanged
    doublets        Scrublet per library           unchanged

The mitochondrial rule used to be `MAD(pct_counts_mt, 3) | (> 20%)`, where
the MAD half made the ceiling irrelevant. Measured by 04b and 04c:

  * global mito median 1.77%, MAD 1.25, so the 3-MAD bound sat at 5.51%
  * the two Indian control animals have a median of 7.58%, ABOVE that bound,
    so most of that population was an outlier by construction
  * 72% of the Indian controls were removed against 6.5% of the Chinese
  * of the Indian cells removed, 18,418 went to the MAD rule and only 2,253
    actually exceeded 20%
  * the survivors were not even a cleaner population: 4x the mito of the
    Chinese survivors and MORE genes, 724 against 666

04b Cell 6 also established this is an animal effect rather than a library
one. Inside libraries C and F, which carry Indian and Chinese animals in the
same GEM well, the Chinese cells lose ~16% and the Indian cells lose 54% and
91%.

A fixed 20% is the field-standard rule, states something absolute about cell
quality, and cannot be defeated by a sample that is uniformly degraded, which
a purely relative rule can. README.md and SESSION_CONTEXT.md describe QC as
per-sample one-sided MAD with MT off; neither has ever matched the code, and
both should now be corrected to the table above.

One wrinkle left alone deliberately: the two-sided bound on
log1p_total_counts works out to 139-8,436 UMI, so the richest cells are
removed along with the debris. That is 751 cells and it is where a cell
carrying both host and viral transcripts sits. The viral whitelist below
protects those specific barcodes, which is the narrower fix.


TWO ADDITIONS, BOTH VISIBLE RATHER THAN SILENT
-----------------------------------------------
1. VIRAL WHITELIST. Any barcode carrying at least MIN_WHITELIST_UMI viral UMI
   is exempt from outlier and doublet removal, and is flagged with
   `qc_would_have_dropped` so the exemption is auditable rather than hidden.
   A cell carrying host plus viral transcripts can read as a doublet to
   Scrublet, which makes that a false positive in the filter rather than a
   real exclusion.
2. DROPPED-CELL LEDGER. adata_prefilter.h5ad holds every called cell before
   any filtering, and qc_dropped_ledger.csv records each dropped barcode with
   its metrics and the reason. SESSION_CONTEXT.md currently tells you to load
   adata_qc.h5ad to recover a QC-dropped SHIV+ cell, which has never been
   possible: the old script wrote that file AFTER the outlier filter, the gene
   and cell minimums and doublet removal. Nothing on disk held the dropped
   droplets. Now something does.


Environment:  conda activate sc_pre
Author: Jake Lehle
Date:   October 2026
"""

# %% ============================================================
# CELL 1 — CONFIGURATION (every tunable lives here)
# ============================================================
import os
import glob

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy.stats import median_abs_deviation
from scipy import sparse

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import warnings
warnings.filterwarnings('ignore')

try:
    from shiv_config import WORKDIR                                 # noqa: F401
except Exception:                                                   # noqa: BLE001
    WORKDIR = '/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV'

# --- inputs -----------------------------------------------------------------
# The v11 run. Point this at Analysis/ only to reproduce the v10 object.
ANALYSIS_DIR = os.path.join(WORKDIR, 'pipeline_results', '03_cellranger_multi')
PROVIDER_DIR = '/master/zwallis/SHARED/10X_ANALYSIS'
SAMPLE_TABLE = os.path.join(WORKDIR, 'sample_table.csv')

# --- the viral feature, v11 collapsed GTF ----------------------------------
SHIV_GENE = 'SHIVAD8EO'              # gene_name in the collapsed GTF
SHIV_GENE_ID = 'SHIVAD8EO_provirus'  # gene_id, used only as a fallback lookup
HOST_COLLISION = ['TAT']             # tyrosine aminotransferase, NOT viral
REQUIRE_SHIV_FEATURE = True          # hard stop if absent; see the header

# --- QC ---------------------------------------------------------------------
# MITOCHONDRIAL: one absolute threshold, no MAD. Decided 2026-11-02 from the
# sweep in scripts/04c_mito_threshold_sweep.py.
#
# The previous rule was `MAD(pct_counts_mt, 3) | (pct_counts_mt > 20)`, and
# the MAD half made the ceiling irrelevant. 04b measured the global mito
# median at 1.77% with MAD 1.25, putting the 3-MAD bound at 5.51%, while the
# two Indian control animals have a median of 7.58%. More than half of that
# population was an outlier by construction, and 72% of it was removed, for
# 6.5% of the Chinese animals. Of the Indian cells removed, 18,418 went to the
# MAD rule and only 2,253 exceeded 20%.
#
# The filter also failed at its own job: the surviving Indian cells still
# carried 4x the mito of the surviving Chinese cells AND had more genes
# (724 vs 666). It cost the control arm its power without making the groups
# comparable.
#
# 04b Cell 6 showed this is an animal effect, not a library one: inside
# libraries C and F, which carry both, the Chinese cells lose ~16% while the
# Indian cells lose 54% and 91%.
#
# A fixed 20% is the field-standard rule, needs no explanation in Methods, and
# unlike a purely relative rule it cannot be defeated by a sample that is
# uniformly degraded. Methods sentence: "Cells with more than 20%
# mitochondrial reads were removed."
MT_THRESHOLD = 20

# COUNTS / GENES / COMPLEXITY: unchanged, two-sided MAD over the whole object.
# One known wrinkle, left alone deliberately rather than changed in the same
# pass: two-sided on log1p_total_counts works out to 139-8,436 UMI, so the
# richest cells are removed as well as the debris. That is 751 cells, and it
# is where a cell carrying both host and viral transcripts would sit. The
# viral whitelist below already protects those specific barcodes.
MAD_COUNTS_GENES = 5
MIN_GENES = 200
MIN_CELLS_PER_GENE = 20
EXPECTED_DOUBLET_RATE = 0.08

# Scrublet's automatic threshold is failing on this dataset: six of eight
# libraries return 0.0% and library G returns 5.9%, for 1,134 total against an
# expected ~8%. That is the known failure mode where the threshold lands in the
# tail of a unimodal simulated-doublet distribution.
#
# Removal is therefore OFF. Deleting 1,121 cells from library G and ~0 from the
# other seven is not doublet removal, it is a batch artifact written into the
# object before the embedding runs. Any threshold we substituted here would be
# invented rather than measured.
#
# doublet_score and predicted_doublet still travel with the object, and 05
# handles doublets at annotation by flagging subclusters on marker
# co-expression (CD3D with CD79A, CD3D with LYZ), which is evidence rather
# than a guessed cut.
REMOVE_DOUBLETS = False

# --- viral whitelist --------------------------------------------------------
WHITELIST_VIRAL_CELLS = True
MIN_WHITELIST_UMI = 1

# Ingestion is deterministic and happens before any QC, so a rerun can reuse
# the prefilter object instead of re-reading 8 matrices and 31 provider
# matrices. Set False to re-ingest from scratch.
REUSE_PREFILTER = True

# --- output -----------------------------------------------------------------
OUT_DIR = os.path.join(WORKDIR, 'pipeline_results', '04_ingest_qc')
PREFILTER_SHARED = os.path.join(OUT_DIR, 'adata_prefilter.h5ad')
FIG_DIR = os.path.join(OUT_DIR, 'figures')
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)
sc.settings.figdir = FIG_DIR

FIG_DPI = 300
FS_LABEL, FS_TICK = 30, 28
HEX_KEEP, HEX_DROP, HEX_VIRAL = '#2E7D32', '#C62828', '#1565C0'

pd.set_option('display.width', 220)


def hr(t):
    print('\n' + '=' * 94, flush=True)
    print(t, flush=True)
    print('=' * 94, flush=True)


def save_csv(df, name, index=False):
    df.to_csv(os.path.join(OUT_DIR, name), index=index, lineterminator='\n')
    print(f'  wrote {name}  ({df.shape[0]:,} rows)', flush=True)


# %% ============================================================
# CELL 2 — METADATA (corrected library and per-animal mapping)
# ============================================================
hr('CELL 2 — metadata')

# Libraries are SEQUENCING POOLS, not conditions. Tissue is unambiguous at the
# library level; timepoint and subspecies are Mixed for C, D and F and are
# resolved per animal from the provider's genetic demultiplexing.
LIBRARY_META = {
    'SHIV_Pre_PBMC':              {'library': 'A', 'tissue': 'PBMC', 'timepoint': 'Pre',      'subspecies': 'Chinese'},
    'SHIV_Pre_LN':                {'library': 'B', 'tissue': 'LN',   'timepoint': 'Pre',      'subspecies': 'Chinese'},
    'SHIV_Pre_PBMC_IndianMixed':  {'library': 'C', 'tissue': 'PBMC', 'timepoint': 'Mixed',    'subspecies': 'Mixed'},
    'SHIV_21DPI_Pre_LN_Mixed':    {'library': 'D', 'tissue': 'LN',   'timepoint': 'Mixed',    'subspecies': 'Mixed'},
    'SHIV_21DPI_PBMC':            {'library': 'E', 'tissue': 'PBMC', 'timepoint': '21 DPI',   'subspecies': 'Chinese'},
    'SHIV_21DPI_NonInf_LN_Mixed': {'library': 'F', 'tissue': 'LN',   'timepoint': 'Mixed',    'subspecies': 'Mixed'},
    'SHIV_Necropsy_PBMC':         {'library': 'G', 'tissue': 'PBMC', 'timepoint': 'Necropsy', 'subspecies': 'Chinese'},
    'SHIV_Necropsy_LN':           {'library': 'H', 'tissue': 'LN',   'timepoint': 'Necropsy', 'subspecies': 'Chinese'},
}

ANIMAL_META = {}
for _lib, _animals in {
    'A': [('39272', 'Pre', 'Chinese'), ('40702', 'Pre', 'Chinese'),
          ('40707', 'Pre', 'Chinese'), ('41861', 'Pre', 'Chinese')],
    'B': [('39272', 'Pre', 'Chinese'), ('40702', 'Pre', 'Chinese'),
          ('40707', 'Pre', 'Chinese'), ('41861', 'Pre', 'Chinese')],
    'C': [('34315', 'Non-Infected', 'Indian'), ('41862', 'Pre', 'Chinese'),
          ('41903', 'Non-Infected', 'Indian')],
    'D': [('39272', '21 DPI', 'Chinese'), ('40702', '21 DPI', 'Chinese'),
          ('40707', '21 DPI', 'Chinese'), ('41862', 'Pre', 'Chinese')],
    'E': [('39272', '21 DPI', 'Chinese'), ('40702', '21 DPI', 'Chinese'),
          ('40707', '21 DPI', 'Chinese'), ('41861', '21 DPI', 'Chinese')],
    'F': [('34315', 'Non-Infected', 'Indian'), ('41861', '21 DPI', 'Chinese'),
          ('41862', '21 DPI', 'Chinese'), ('41903', 'Non-Infected', 'Indian')],
    'G': [('40702', 'Necropsy', 'Chinese'), ('40707', 'Necropsy', 'Chinese'),
          ('41861', 'Necropsy', 'Chinese'), ('41862', 'Necropsy', 'Chinese')],
    'H': [('40702', 'Necropsy', 'Chinese'), ('40707', 'Necropsy', 'Chinese'),
          ('41861', 'Necropsy', 'Chinese'), ('41862', 'Necropsy', 'Chinese')],
}.items():
    for _a, _tp, _ss in _animals:
        ANIMAL_META[(_lib, _a)] = {'timepoint': _tp, 'subspecies': _ss}

ANIMAL_INFO = {
    '34315': {'subspecies': 'Indian',  'infection': 'Non-Infected', 'note': 'control'},
    '41903': {'subspecies': 'Indian',  'infection': 'Non-Infected', 'note': 'control'},
    '39272': {'subspecies': 'Chinese', 'infection': 'Infected',     'note': 'no necropsy'},
    '40702': {'subspecies': 'Chinese', 'infection': 'Infected',     'note': ''},
    '40707': {'subspecies': 'Chinese', 'infection': 'Infected',     'note': 'elite controller'},
    '41861': {'subspecies': 'Chinese', 'infection': 'Infected',     'note': ''},
    '41862': {'subspecies': 'Chinese', 'infection': 'Infected',     'note': 'pooled separately at Pre/21DPI'},
}

for cond, info in LIBRARY_META.items():
    mixed = ' [MIXED]' if info['timepoint'] == 'Mixed' else ''
    print(f"  {info['library']}: {cond:<28} {info['tissue']:<5} "
          f"{info['timepoint']}{mixed}")
print(f'\n  per-animal entries: {len(ANIMAL_META)}   '
      f'animals: {sorted(ANIMAL_INFO)}')


# %% ============================================================
# CELL 3 — BARCODE TO ANIMAL, FROM THE PROVIDER'S DEMULTIPLEXING
# ============================================================
hr('CELL 3 — ingestion')

REUSED = False
if REUSE_PREFILTER and os.path.exists(PREFILTER_SHARED):
    print(f'  reusing {PREFILTER_SHARED}')
    print('  Ingestion is deterministic and happens before any QC, so the')
    print('  prefilter object is mode-independent and is loaded rather than')
    print('  rebuilt. Set REUSE_PREFILTER = False to re-ingest from scratch.')
    adata = sc.read_h5ad(PREFILTER_SHARED)
    print(f'  loaded {adata.n_obs:,} cells x {adata.n_vars:,} genes')
    REUSED = True

barcode_to_animal = {}
if not REUSED:
    if not os.path.exists(SAMPLE_TABLE):
        print(f'  WARN  sample table not found: {SAMPLE_TABLE}')
        print('        every cell will carry animal_id "Unknown".')
    else:
        sample_df = pd.read_csv(SAMPLE_TABLE)
        print(f'  sample table: {len(sample_df)} entries')
        for _, row in sample_df.iterrows():
            lib, animal = row['run_letter'], str(row['animal_id'])
            mp = row['matrix_path']
            if pd.isna(mp) or not os.path.exists(str(mp)):
                print(f'    SKIP {animal}_{lib}: matrix not found')
                continue
            try:
                prov = sc.read_10x_h5(mp, gex_only=True)
                for bc in prov.obs_names:
                    barcode_to_animal[(lib, bc.split('-')[0])] = animal
                print(f'    {animal}_{lib}: {prov.n_obs:,} barcodes')
                del prov
            except Exception as e:                                      # noqa: BLE001
                print(f'    FAIL {animal}_{lib}: {type(e).__name__}: {e}')
    print(f'  total barcode-to-animal entries: {len(barcode_to_animal):,}')


    hr('CELL 4 — read Cell Ranger multi outputs')


    def find_matrix_h5(cond):
        p = os.path.join(ANALYSIS_DIR, cond, 'outs', 'per_sample_outs', cond,
                         'sample_filtered_feature_bc_matrix.h5')
        if os.path.exists(p):
            return p
        hits = sorted(glob.glob(os.path.join(
            ANALYSIS_DIR, cond, 'outs', '*', 'sample_filtered_feature_bc_matrix.h5')))
        return hits[0] if hits else None


    def split_modalities(h5):
        full = sc.read_10x_h5(h5, gex_only=False)
        full.var_names_make_unique()
        col = next((c for c in ('feature_types', 'feature_type')
                    if c in full.var.columns), None)
        if col is None:
            return full, None
        gex = full[:, full.var[col] == 'Gene Expression'].copy()
        adt_mask = full.var[col] == 'Antibody Capture'
        adt = full[:, adt_mask].copy() if adt_mask.sum() else None
        return gex, adt


    def assign_animals(a, lib):
        ids, tps, sss = [], [], []
        for bc_full in a.obs_names:
            base = bc_full.rsplit('-', 1)[0].split('-')[0]
            animal = barcode_to_animal.get((lib, base))
            if animal is None:
                ids.append('Unknown'); tps.append('Unknown'); sss.append('Unknown')
            else:
                m = ANIMAL_META.get((lib, animal), {})
                ids.append(animal)
                tps.append(m.get('timepoint', 'Unknown'))
                sss.append(m.get('subspecies', 'Unknown'))
        a.obs['animal_id'] = ids
        a.obs['timepoint'] = tps
        a.obs['subspecies'] = sss
        n = sum(1 for x in ids if x != 'Unknown')
        print(f'    demux: {n:,}/{len(ids):,} ({100 * n / max(1, len(ids)):.1f}%) assigned')


    gex_list, adt_list, summary = [], [], []
    for cond, info in LIBRARY_META.items():
        lib = info['library']
        h5 = find_matrix_h5(cond)
        if h5 is None:
            print(f'  SKIP {cond} (library {lib}): no matrix under {ANALYSIS_DIR}')
            summary.append({'condition': cond, 'library': lib,
                            'status': 'MISSING', 'n_cells': 0})
            continue
        print(f'  {lib}: {cond}')
        try:
            gex, adt = split_modalities(h5)
            gex.obs['condition'] = cond
            gex.obs['library'] = lib
            gex.obs['tissue'] = info['tissue']
            assign_animals(gex, lib)
            gex.obs_names = [f'{bc}-{lib}' for bc in gex.obs_names]
            gex_list.append(gex)
            print(f'    GEX: {gex.n_obs:,} cells x {gex.n_vars:,} genes')
            if adt is not None:
                adt.obs_names = [f'{bc}-{lib}' for bc in adt.obs_names]
                adt_list.append(adt)
                print(f'    ADT: {adt.n_vars} proteins')
            summary.append({'condition': cond, 'library': lib,
                            'status': 'OK', 'n_cells': int(gex.n_obs)})
        except Exception as e:                                          # noqa: BLE001
            print(f'    FAILED: {type(e).__name__}: {e}')
            summary.append({'condition': cond, 'library': lib,
                            'status': f'ERROR: {e}', 'n_cells': 0})

    if not gex_list:
        raise SystemExit(f'No matrices found under {ANALYSIS_DIR}. '
                         'Has 03_cellranger_multi.sh finished?')

    adata = ad.concat(gex_list, join='outer', merge='same')
    adata.obs['sample_name'] = (adata.obs['animal_id'].astype(str) + '_'
                                + adata.obs['tissue'].astype(str) + '_'
                                + adata.obs['timepoint'].astype(str))
    print(f'\n  combined GEX: {adata.n_obs:,} cells x {adata.n_vars:,} genes')

    if adt_list:
        adt_all = ad.concat(adt_list, join='outer', merge='same')
        common = adata.obs_names.intersection(adt_all.obs_names)
        al = adt_all[common].copy()
        df = pd.DataFrame(al.X.toarray() if sparse.issparse(al.X) else al.X,
                          index=al.obs_names, columns=al.var_names)
        df = df.reindex(adata.obs_names, fill_value=0)
        adata.obsm['protein_counts'] = df.values
        adata.uns['adt_names'] = list(df.columns)
        print(f'  ADT stored: {df.shape[1]} proteins')

    sumdf = pd.DataFrame(summary)
    print()
    print(sumdf.to_string(index=False))
    save_csv(sumdf, 'ingestion_summary.csv')

    n_libs = int((sumdf.status == 'OK').sum())
    if n_libs < 8:
        print(f'\n  WARN  only {n_libs}/8 libraries ingested. Continuing, but the')
        print('        object is incomplete and downstream counts will not match')
        print('        the published 89,844.')


    hr('CELL 5 — gene categories')

    adata.var['mt'] = adata.var_names.str.startswith('KEG06_')   # rhesus mito
    adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL'))
    adata.var['hb'] = adata.var_names.str.contains('^HB[^(P)]', case=False)
    print(f"  mitochondrial (KEG06_*) : {int(adata.var['mt'].sum())}")
    print(f"  ribosomal (RPS/RPL)     : {int(adata.var['ribo'].sum())}")
    print(f"  hemoglobin              : {int(adata.var['hb'].sum())}")

    # EXACT, CASE-SENSITIVE, NO FALLBACK. See the header for why.
    adata.var['shiv'] = adata.var_names == SHIV_GENE
    n_shiv = int(adata.var['shiv'].sum())

    if n_shiv == 0 and SHIV_GENE_ID and 'gene_ids' in adata.var.columns:
        by_id = adata.var['gene_ids'].astype(str) == SHIV_GENE_ID
        if by_id.sum():
            adata.var['shiv'] = by_id.values
            n_shiv = int(by_id.sum())
            print(f"  NOTE  matched the virus by gene_id '{SHIV_GENE_ID}' because "
                  f"gene_name '{SHIV_GENE}' was absent")

    print(f"  SHIV features           : {n_shiv}  "
          f"({', '.join(adata.var_names[adata.var['shiv']])})")

    present_collisions = [g for g in HOST_COLLISION if g in adata.var_names]
    if present_collisions:
        print(f"  host look-alikes in the matrix, NOT tagged as viral: "
              f"{present_collisions}")
        for g in present_collisions:
            assert not adata.var.loc[g, 'shiv'], f'host {g} was tagged as viral'
        print('  verified: none of them is tagged shiv')

    if n_shiv == 0:
        msg = (f"\n  STOP: no viral feature found. Expected gene_name "
               f"'{SHIV_GENE}'.\n"
               f"        This matrix was probably built on the v10 reference, "
               f"whose\n"
               f"        viral features are the ten ORF names. Continuing would "
               f"report\n"
               f"        zero viral signal everywhere, or, with the old "
               f"case-insensitive\n"
               f"        fallback, would report host TAT as virus. Neither fails "
               f"loudly,\n"
               f"        which is why this is a hard stop.")
        if REQUIRE_SHIV_FEATURE:
            raise SystemExit(msg)
        print(msg + '\n        REQUIRE_SHIV_FEATURE is False, continuing anyway.')

    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo', 'hb', 'shiv'],
                               inplace=True, percent_top=[20], log1p=True)

    adata.obs['shiv_umi_cellranger'] = np.asarray(
        adata[:, adata.var['shiv']].X.sum(axis=1)).ravel() if n_shiv else 0.0
    n_viral_cells = int((adata.obs['shiv_umi_cellranger'] > 0).sum())
    print(f"\n  cells with >=1 viral UMI: {n_viral_cells:,} / {adata.n_obs:,}")
    print(f"  total viral UMI         : "
          f"{adata.obs['shiv_umi_cellranger'].sum():,.0f}")
    for tp in sorted(adata.obs['timepoint'].astype(str).unique()):
        m = adata.obs['timepoint'].astype(str) == tp
        print(f'    {tp:<14} {int((adata.obs["shiv_umi_cellranger"][m] > 0).sum()):>6} '
              f'viral cells / {int(m.sum()):>7,}')


# %% ============================================================
# CELL 5b — VERIFY THE VIRAL TAG (runs on BOTH paths)
# ============================================================
# This must run whether the object was just built or reloaded from disk.
# Cells 4 and 5 sit under `if not REUSED:`, so on the reuse path the viral
# tagging block never executes. That left n_shiv undefined (the NameError at
# Cell 7) and, more seriously, SKIPPED THE HOST TAT GUARD, which is the one
# check in this script that exists to stop a confidently wrong answer.
# Everything below reads state off the object itself rather than off variables
# set earlier, so it is correct on both paths.
hr('CELL 5b — verify the viral tag')

required_obs = ['library', 'animal_id', 'timepoint', 'tissue',
                'pct_counts_mt', 'total_counts', 'n_genes_by_counts',
                'log1p_total_counts', 'log1p_n_genes_by_counts',
                'pct_counts_in_top_20_genes']
missing_obs = [c for c in required_obs if c not in adata.obs.columns]
if missing_obs:
    raise SystemExit(
        f'  STOP: the object is missing obs columns {missing_obs}.\n'
        '        If it was reloaded, it predates this script. Set\n'
        '        REUSE_PREFILTER = False and re-ingest.')

if 'shiv' not in adata.var.columns:
    raise SystemExit(
        '  STOP: adata.var has no "shiv" column.\n'
        '        A reused prefilter object was written before the viral tag\n'
        '        existed. Set REUSE_PREFILTER = False and re-ingest.')

n_shiv = int(adata.var['shiv'].sum())
viral_names = list(adata.var_names[adata.var['shiv']])
print(f'  SHIV features tagged : {n_shiv}'
      + (f'  ({", ".join(viral_names)})' if viral_names else ''))

# The reason this script has a hard stop at all. Under the v11 collapsed GTF
# none of the ten ORF names exist, and the old case-insensitive fallback would
# match host TAT, tyrosine aminotransferase, reporting it as viral UMI without
# raising anything.
present_collisions = [g for g in HOST_COLLISION if g in adata.var_names]
if present_collisions:
    print(f'  host look-alikes in the matrix : {present_collisions}')
    for g in present_collisions:
        if bool(adata.var.loc[g, 'shiv']):
            raise SystemExit(
                f'  STOP: host gene {g} is tagged as viral. That is the exact\n'
                '        failure this guard exists to catch.')
    print('  verified: none of them is tagged shiv')

if n_shiv == 0 and REQUIRE_SHIV_FEATURE:
    raise SystemExit(
        f'  STOP: no viral feature found. Expected gene_name "{SHIV_GENE}".\n'
        '        This matrix was probably built on the v10 reference, whose\n'
        '        viral features are the ten ORF names. Continuing would report\n'
        '        zero viral signal everywhere and would not fail loudly.')

if 'shiv_umi_cellranger' not in adata.obs.columns:
    adata.obs['shiv_umi_cellranger'] = (
        np.asarray(adata[:, adata.var['shiv']].X.sum(axis=1)).ravel()
        if n_shiv else 0.0)
    print('  recomputed shiv_umi_cellranger from the matrix')

_vc = int((adata.obs['shiv_umi_cellranger'] > 0).sum())
print(f'  cells with >=1 viral UMI : {_vc:,} / {adata.n_obs:,}')
print(f'  total viral UMI          : '
      f'{adata.obs["shiv_umi_cellranger"].sum():,.0f}')


# %% ============================================================
# CELL 6 — SAVE THE PRE-FILTER OBJECT (so dropped cells are recoverable)
# ============================================================
hr('CELL 6 — pre-filter object')

prefilter_path = os.path.join(OUT_DIR, 'adata_prefilter.h5ad')
if REUSED and os.path.abspath(prefilter_path) == os.path.abspath(PREFILTER_SHARED):
    print('  prefilter was reused from this same path; not rewriting it.')
else:
    adata.write(prefilter_path)
    print(f'  wrote adata_prefilter.h5ad  ({adata.n_obs:,} cells)')
print('  This is what SESSION_CONTEXT.md has always claimed adata_qc.h5ad was.')
print('  It is the only object from which a QC-dropped barcode can be recovered.')

prefilter_obs = adata.obs.copy()


# %% ============================================================
# CELL 7 — QC FILTER (frozen) WITH THE VIRAL WHITELIST
# ============================================================
hr('CELL 7 — QC')


def is_outlier(a, metric, nmads):
    """Two-sided MAD over the whole object.

    Unchanged from scanpy_analysis_script.py. Applies to counts, genes and
    the top-20% complexity metric only; the mitochondrial filter is a fixed
    threshold and does not go through here.
    """
    M = a.obs[metric]
    med = M.median()
    mad = median_abs_deviation(M, nan_policy='omit')
    if mad == 0:
        mad = M.std()
        if mad == 0:
            return pd.Series(False, index=M.index)
    mask = np.abs(M - med) > nmads * mad
    lo, hi = med - nmads * mad, med + nmads * mad
    print(f'    {metric}: {int(mask.sum()):,} outliers '
          f'({100 * mask.mean():.1f}%) [median={med:.2f}, MAD={mad:.2f}, '
          f'bounds {lo:.2f} to {hi:.2f}]')
    return mask


print('  outlier detection (two-sided, global, 5 MAD):')
adata.obs['outlier'] = (
    is_outlier(adata, 'log1p_total_counts', MAD_COUNTS_GENES)
    | is_outlier(adata, 'log1p_n_genes_by_counts', MAD_COUNTS_GENES)
    | is_outlier(adata, 'pct_counts_in_top_20_genes', MAD_COUNTS_GENES))
# MITOCHONDRIAL FILTER: a single absolute threshold, no MAD. See the header
# of this script for the measurement that drove the change.
adata.obs['mt_outlier'] = adata.obs['pct_counts_mt'] > MT_THRESHOLD
print(f'\n  mito filter: pct_counts_mt > {MT_THRESHOLD}%  -> '
      f'{int(adata.obs["mt_outlier"].sum()):,} cells '
      f'({100 * adata.obs["mt_outlier"].mean():.1f}%)')
for _g, _m in [('Indian controls', adata.obs['animal_id'].astype(str)
                .isin(['34315', '41903'])),
               ('Chinese infected', adata.obs['animal_id'].astype(str)
                .isin(['39272', '40702', '40707', '41861', '41862']))]:
    if _m.any():
        print(f'    {_g:<18} {int(adata.obs.loc[_m, "mt_outlier"].sum()):>6,} '
              f'of {int(_m.sum()):>7,} '
              f'({100 * adata.obs.loc[_m, "mt_outlier"].mean():.1f}%)')

# --- the whitelist ----------------------------------------------------------
adata.obs['shiv_whitelist'] = (
    (adata.obs['shiv_umi_cellranger'] >= MIN_WHITELIST_UMI)
    if WHITELIST_VIRAL_CELLS else False)
would_drop = adata.obs['outlier'] | adata.obs['mt_outlier']
adata.obs['qc_would_have_dropped'] = would_drop & adata.obs['shiv_whitelist']

n_wl = int(adata.obs['shiv_whitelist'].sum())
n_rescued = int(adata.obs['qc_would_have_dropped'].sum())
print(f'\n  viral whitelist: {n_wl:,} barcodes carry >={MIN_WHITELIST_UMI} viral UMI')
print(f'  of those, {n_rescued} would have been dropped by QC and are kept,')
print('  flagged qc_would_have_dropped so the exemption is auditable.')

drop_mask = would_drop & ~adata.obs['shiv_whitelist']

# --- the ledger -------------------------------------------------------------
ledger = prefilter_obs.loc[would_drop, [
    'library', 'animal_id', 'timepoint', 'tissue', 'total_counts',
    'n_genes_by_counts', 'pct_counts_mt', 'pct_counts_in_top_20_genes',
    'shiv_umi_cellranger']].copy()
ledger['reason'] = np.where(
    adata.obs.loc[would_drop, 'mt_outlier'], 'mt_outlier', 'count_gene_outlier')
ledger['whitelisted_kept'] = adata.obs.loc[would_drop, 'shiv_whitelist'].values
ledger.index.name = 'barcode'
save_csv(ledger.reset_index(), 'qc_dropped_ledger.csv')

n_before = adata.n_obs
adata = adata[~drop_mask].copy()
print(f'\n  removed by outlier/mito filter: {n_before - adata.n_obs:,}')

sc.pp.filter_cells(adata, min_genes=MIN_GENES)
sc.pp.filter_genes(adata, min_cells=MIN_CELLS_PER_GENE)
print(f'  after min_genes={MIN_GENES}, min_cells={MIN_CELLS_PER_GENE}: '
      f'{adata.n_obs:,} cells x {adata.n_vars:,} genes')

if n_shiv and SHIV_GENE not in adata.var_names:
    print(f'\n  WARN  the viral feature was removed by '
          f'min_cells={MIN_CELLS_PER_GENE},')
    print('        because it is detected in fewer cells than that. Viral')
    print('        counts are preserved in obs as shiv_umi_cellranger, which')
    print('        was computed before gene filtering.')


# %% ============================================================
# CELL 8 — DOUBLET SCORING (per library) WITH THE SAME WHITELIST
# Scoring only by default. REMOVE_DOUBLETS in Cell 1 controls removal.
# ============================================================
hr('CELL 8 — Scrublet')

try:
    import scrublet as scr
    n_before = adata.n_obs
    scores = np.zeros(adata.n_obs)
    preds = np.zeros(adata.n_obs, dtype=bool)
    for lib in sorted(adata.obs['library'].astype(str).unique()):
        m = (adata.obs['library'].astype(str) == lib).values
        idx = np.where(m)[0]
        if m.sum() < 100:
            print(f'  library {lib}: <100 cells, skipped')
            continue
        try:
            s, p = scr.Scrublet(
                adata[m].X, expected_doublet_rate=EXPECTED_DOUBLET_RATE
            ).scrub_doublets(min_counts=2, min_cells=3, verbose=False)
            scores[idx], preds[idx] = s, p
            print(f'  library {lib}: {int(p.sum()):,} doublets '
                  f'({100 * p.mean():.1f}%)')
        except Exception as e:                                      # noqa: BLE001
            print(f'  library {lib}: Scrublet failed — {type(e).__name__}: {e}')

    adata.obs['doublet_score'] = scores
    adata.obs['predicted_doublet'] = preds

    # A cell carrying host AND viral transcripts can read as a doublet. That
    # makes it a false positive in the filter, not a real exclusion, so the
    # whitelist applies here too.
    wl = adata.obs['shiv_whitelist'].values.astype(bool)
    adata.obs['doublet_whitelisted'] = preds & wl
    n_dwl = int(adata.obs['doublet_whitelisted'].sum())
    if n_dwl:
        print(f'  {n_dwl} viral barcode(s) called doublet and KEPT '
              f'(flagged doublet_whitelisted)')
    if REMOVE_DOUBLETS:
        adata = adata[~(preds & ~wl)].copy()
        print(f'  doublets removed: {n_before - adata.n_obs:,}')
    else:
        print(f'  REMOVE_DOUBLETS = False, so none were removed.')
        print(f'  {int(preds.sum()):,} cells carry predicted_doublet = True '
              f'and stay in the object.')
        print('  Scrublet is misbehaving on this dataset (see Cell 1), so a')
        print('  threshold here would be invented. 05 flags doublet-like')
        print('  subclusters on marker co-expression instead.')
except ImportError:
    print('  scrublet not installed; skipped')
    adata.obs['doublet_score'] = 0.0
    adata.obs['predicted_doublet'] = False
    adata.obs['doublet_whitelisted'] = False

print(f'\n  FINAL: {adata.n_obs:,} cells x {adata.n_vars:,} genes')


# %% ============================================================
# CELL 9 — SAVE AND SUMMARISE
# ============================================================
hr('CELL 9 — save')

qc_path = os.path.join(OUT_DIR, 'adata_qc.h5ad')
adata.write(qc_path)
print(f'  wrote adata_qc.h5ad  ({adata.n_obs:,} cells)')

surv = adata.obs['shiv_umi_cellranger'] > 0
wl_df = adata.obs.loc[surv, ['library', 'animal_id', 'timepoint', 'tissue',
                             'shiv_umi_cellranger', 'shiv_whitelist',
                             'qc_would_have_dropped', 'doublet_whitelisted']]
wl_df.index.name = 'barcode'
save_csv(wl_df.reset_index(), 'shiv_whitelist.csv')

rows = [{'stage': 'cells loaded (pre-filter)', 'n_cells': len(prefilter_obs)},
        {'stage': 'after QC + doublets', 'n_cells': int(adata.n_obs)},
        {'stage': 'viral cells pre-filter',
         'n_cells': int((prefilter_obs['shiv_umi_cellranger'] > 0).sum())},
        {'stage': 'viral cells surviving', 'n_cells': int(surv.sum())}]
qsum = pd.DataFrame(rows)
print()
print(qsum.to_string(index=False))
save_csv(qsum, 'qc_summary.csv')

print('\n  cells per animal:')
print(adata.obs['animal_id'].value_counts().to_string())
print('\n  cells per timepoint:')
print(adata.obs['timepoint'].value_counts().to_string())


hr('DONE — 04_ingest_qc.py')
print(f'  outputs: {OUT_DIR}')
print(f'    adata_prefilter.h5ad   every called cell, nothing filtered')
print(f'    adata_qc.h5ad          the object stage 05 continues from')
print(f'    qc_dropped_ledger.csv  every dropped barcode and why')
print('\n  Compare total cells against the v10 object (89,844 after QC). A')
print('  large shift means the reference change moved the QC thresholds and')
print('  the Tier 2 cluster fingerprints will need re-adjudication.')
print('\n  Next: 05_shiv_merge.py')
