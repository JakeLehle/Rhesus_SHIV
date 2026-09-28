#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Rhesus SHIV - Tier 2 Myeloid Annotation Write
==============================================
Writes final Tier-2 myeloid labels onto the purity-gated myeloid subset from
shiv_tier2_myeloid_diagnostic_v2.py, carries the two state scores and the QC
context as per-cell attributes rather than labels, and propagates the labels back
onto the global object WITHOUT touching any non-myeloid cell.

Zoey's constraint (2026-09-25) is enforced structurally: these labels only ever
subdivide the Tier-1 myeloid population. CELL 6 refuses to write onto a cell
whose tier1_celltype is not in MYELOID_TIER1, and asserts it afterwards.

LABEL MAP, locked after reading the v2 diagnostic output (9,698 gated cells,
15 post-gate subclusters). Every entry has its basis recorded so we do not
relitigate it:

  5  (890)  Non-classical / CD163+ Mono-Mac   z 2.64 margin 2.43; CFD, FCGR3,
                                              MS4A7, LST1, LILRB1, CD68, ITGB2
  11 (189)  Tissue Macrophage (C1Q+/CD206+)   z 3.69 margin 2.31; 99.5% LN;
                                              C1QA/B/C, APOE, APOC1, SELENOP,
                                              DNASE1L3, LIPA, CD68, GPNMB
  12 (163)  pDC                               z 3.72 margin 4.02; TCF4, LILRA4,
                                              IRF8, JCHAIN, BLNK, CCDC50, SCT
  13 (149)  Migratory / Mature DC             z 3.72 margin 1.90; 100% LN;
                                              LAMP3, FSCN1, CD83, BIRC3, LY75,
                                              FLT3, CTSH
  9  (368)  Antigen Presenting DC             OVERRIDE. argmax said cDC1 at
                                              z 3.03, but cDC1 RAW is 0.350 and
                                              no cluster in the compartment
                                              exceeds 0.350 for cDC1, while this
                                              cluster's class II raw is 1.915.
                                              Markers are CD74 + MAMU-DRA/DQA1/
                                              DQB1/DRB1/DRB5/DPA and contain NO
                                              CLEC9A/XCR1/BATF3/CADM1. The z
                                              ranking inverted because cDC1 is
                                              flat near zero everywhere, so a
                                              small bump z-scores large. See
                                              RAW_FLOOR note below.
  2  (1068) Classical Monocyte                z 0.85 margin 0.65; LYZ, VCAN,
                                              CD14, CD163, FN1, PLBD1
  0  (1469) Classical Monocyte  [QC-FLAGGED]  z 1.08 margin 0.78 BUT 79.9%
  4  (970)  Classical Monocyte  [QC-FLAGGED]  animal 34315 / 89.5% animal 34315,
                                              82.7% / 90.1% Non-Infected, mito z
                                              1.98 / 1.62, five KEG06_ genes in
                                              the top 15, plus CAMP. Genuine
                                              classical monocytes but confounded
                                              with one Indian control animal and
                                              a stress signature. Labelled so
                                              they are not lost, qc_exclude=True
                                              so they stay out of DEG.
  1  (1222) Monocyte, MHC-II-high             no call (margin 0.11). CD74 +
                                              MAMU-DRA/DPA/DRB1, state_APC 0.58,
                                              IEG z 1.92. Lineage is monocyte;
                                              the competing signal is a state.
                                              Labelled at the lineage level with
                                              the state carried separately.
  7  (627)  Monocyte, IFN-stimulated          no call. IFI27, IFI6, state_IFN
                                              1.82, 67.3% 21 DPI, 62% library E.
                                              Same treatment as cluster 1.
  8  (472)  Unresolved myeloid  [QC-FLAGGED]  no call; 66.7% animal 41861, IEG
                                              z 1.64, MHC-II + BTG2/ZFP36.
                                              Single-animal artifact-suspect.
  3  (972)  Unresolved myeloid                no call; nothing above z -0.30,
  10 (236)  Unresolved myeloid                ribosomal z -2.87 (cluster 3).
                                              Kept, unlabelled at subtype level.
  6  (771)  EXCLUDE - T/NK contamination      SKAP1, NKG7, GZMB, GNLY, CCL5,
                                              GZMH, CD3E, CD247, CAMK4;
                                              z_Contam_T 3.29, z_Contam_NK 3.61.
                                              Survived the v2 gate because the
                                              gate's secondary condition requires
                                              contamination to beat the myeloid
                                              score, and this cluster has a high
                                              myeloid score too (1.148) - the
                                              signature of a myeloid-lymphocyte
                                              DOUBLET. The AND-veto is wrong for
                                              that case; see GATE_FIX note.
  14 (132)  EXCLUDE - platelets               PPBP, NRGN, PF4V1, CAVIN2, CD9,
                                              LOC703451. Called Classical
                                              Monocyte at margin 0.78 because the
                                              v2 contamination panels covered
                                              T/B/NK but NOT platelets or
                                              erythrocytes, even though Tier 1
                                              reports 482 platelets and 300
                                              erythrocytes.

THREE DIAGNOSTIC GAPS this script works around. All three need fixing in a v3 of
the diagnostic, but none requires a rerun before annotating, because each is
resolvable by cluster identity here:

  RAW_FLOOR  - the call gate checks a z floor and a margin but not a RAW floor.
               A signature that is flat near zero across all clusters produces a
               large z from a trivial bump. Cluster 9 is the case. v3 should
               require the raw signature mean to clear a floor as well.
  GATE_FIX   - the cluster purity gate ANDs the MAD-outlier test with
               "contamination beats myeloid". Doublets satisfy the first and fail
               the second. The outlier test should stand alone, with the myeloid
               comparison reported rather than used as a veto.
  CONTAM_SET - add platelet (PPBP, PF4, GP9, NRGN, CAVIN2) and erythrocyte
               (HBB, HBA1, HBA2, ALAS2) panels to CONTAM_PANELS.

QC-excluded cells are FLAGGED, never deleted. The DEG script respects qc_exclude.
Both input objects are left untouched; output goes to new files.

Author: Jake Lehle / Kaushal Lab
Date: September 2026
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
import warnings
warnings.filterwarnings("ignore")

WORKING_DIR = '/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV'
MERGED      = os.path.join(WORKING_DIR, 'shiv_starsolo', 'merged')

SUBSET_IN   = os.path.join(MERGED, 'shiv_Myeloid_subset_reclustered_v2.h5ad')
GLOBAL_IN   = os.path.join(MERGED, 'shiv_host_tier2_annotated.h5ad')
SUBSET_OUT  = os.path.join(MERGED, 'shiv_Myeloid_annotated.h5ad')
GLOBAL_OUT  = os.path.join(MERGED, 'shiv_host_tier2_myeloid_annotated.h5ad')
OUT_DIR     = os.path.join(WORKING_DIR, 'annotation_output', 'tier2_myeloid_annotation')

DRY_RUN = False        # True = print everything, write nothing

CLUSTER_KEY   = 'sub_leiden'
TIER1_COL     = 'tier1_celltype'
SUBTYPE_COL   = 'tier2_myeloid_subtype'     # NEW on subset and global
GLOBAL_T2COL  = 'tier2_celltype'            # existing T-cell-aware column
ANIMAL_COL    = 'animal_id'
TIME_COL      = 'timepoint'
TISSUE_COL    = 'tissue'
MYELOID_TIER1 = ['Monocytes', 'Macrophages', 'DCs']

# ---- fingerprint guard ----------------------------------------------------
# The label map keys on post-gate cluster IDs. Reclustering is seeded
# (random_state=0) so it reproduces, but if any diagnostic dial changed the IDs
# would shift silently and every label would be wrong. Refuse to proceed unless
# the cluster sizes match what the map was read from.
EXPECTED_SIZES = {
    '0': 1469, '1': 1222, '2': 1068, '3': 972, '4': 970, '5': 890, '6': 771,
    '7': 627, '8': 472, '9': 368, '10': 236, '11': 189, '12': 163, '13': 149,
    '14': 132,
}
ENFORCE_FINGERPRINT = True
FINGERPRINT_TOLERANCE = 0      # exact match required; set >0 to allow drift

# ---- the label map --------------------------------------------------------
LABEL_MAP = {
    '0':  'Classical Monocyte',
    '1':  'Monocyte, MHC-II-high',
    '2':  'Classical Monocyte',
    '3':  'Unresolved myeloid',
    '4':  'Classical Monocyte',
    '5':  'Non-classical / CD163+ Mono-Mac',
    '6':  'EXCLUDE: T/NK contamination',
    '7':  'Monocyte, IFN-stimulated',
    '8':  'Unresolved myeloid',
    '9':  'Antigen Presenting DC',
    '10': 'Unresolved myeloid',
    '11': 'Tissue Macrophage (C1Q+/CD206+)',
    '12': 'pDC',
    '13': 'Migratory / Mature DC',
    '14': 'EXCLUDE: Platelets',
}

# clusters flagged qc_exclude=True, with the reason recorded per cell
QC_FLAG = {
    '0':  'mito_high_single_animal_34315',
    '4':  'mito_high_single_animal_34315',
    '8':  'single_animal_41861_IEG_high',
    '6':  'lymphocyte_contamination_or_doublet',
    '14': 'platelet_contamination',
}
# clusters that are not myeloid at all and must never enter a myeloid analysis
NON_MYELOID_CLUSTERS = ['6', '14']

# state score columns carried from the diagnostic (per-cell, never labels)
STATE_COLS = ['sig_Antigen_Presenting_Myeloid', 'sig_IFN_Stimulated_Myeloid']
QC_SCORE_COLS = ['qc_QC_mitochondrial', 'qc_QC_ribosomal', 'qc_QC_immediate_early',
                 'qc_Contam_T', 'qc_Contam_B', 'qc_Contam_NK']

FIG_DPI = 300
plt.rcParams.update({
    'font.size': 28, 'axes.titlesize': 34, 'axes.labelsize': 30,
    'xtick.labelsize': 24, 'ytick.labelsize': 24, 'legend.fontsize': 18,
    'pdf.fonttype': 42, 'ps.fonttype': 42,
})
os.makedirs(OUT_DIR, exist_ok=True)
sc.settings.figdir = OUT_DIR
pd.set_option('display.width', 220)
pd.set_option('display.max_columns', 60)
pd.set_option('display.max_rows', 300)


def hr(t):
    print('\n' + '=' * 90)
    print(t)
    print('=' * 90)


def save_csv(df, name, index=True):
    p = os.path.join(OUT_DIR, name)
    df.to_csv(p, index=index)
    print(f'    saved: {name}  ({len(df)} rows)')


# %% ============================================================
# CELL 2 - LOAD + FINGERPRINT GUARD
# ============================================================
hr('CELL 2 - load subset and verify cluster fingerprint')

if not os.path.exists(SUBSET_IN):
    raise FileNotFoundError(f'myeloid subset not found: {SUBSET_IN}\n'
                            'Run shiv_tier2_myeloid_diagnostic_v2.py first.')
sub = sc.read_h5ad(SUBSET_IN)
print(f'  loaded: {sub.n_obs:,} cells x {sub.n_vars:,} genes')

if CLUSTER_KEY not in sub.obs.columns:
    raise KeyError(f'{CLUSTER_KEY} not in obs. Available: {list(sub.obs.columns)}')
sub.obs[CLUSTER_KEY] = sub.obs[CLUSTER_KEY].astype(str)
sizes = sub.obs[CLUSTER_KEY].value_counts().to_dict()

print('\n  cluster sizes, observed vs expected:')
mismatch = []
allk = sorted(set(sizes) | set(EXPECTED_SIZES), key=lambda x: int(x))
for k in allk:
    obs_n, exp_n = sizes.get(k, 0), EXPECTED_SIZES.get(k, 0)
    flag = '' if abs(obs_n - exp_n) <= FINGERPRINT_TOLERANCE else '   <-- MISMATCH'
    if flag:
        mismatch.append((k, obs_n, exp_n))
    print(f'    {k:>3s}: observed {obs_n:6,}   expected {exp_n:6,}{flag}')

if mismatch:
    msg = (f'\n  *** FINGERPRINT MISMATCH on {len(mismatch)} cluster(s).\n'
           '  The label map was read from a specific diagnostic run. If the\n'
           '  clustering has shifted, applying this map would assign wrong labels\n'
           '  to every cell. Re-read the diagnostic output and rebuild LABEL_MAP,\n'
           '  or set ENFORCE_FINGERPRINT = False only if you have confirmed the\n'
           '  cluster identities yourself.')
    print(msg)
    if ENFORCE_FINGERPRINT:
        raise RuntimeError('Cluster fingerprint does not match LABEL_MAP. Refusing '
                           'to write labels.')
else:
    print('\n  fingerprint OK, cluster identities match the label map')

unmapped = [k for k in sizes if k not in LABEL_MAP]
if unmapped:
    raise RuntimeError(f'clusters with no LABEL_MAP entry: {unmapped}')


# %% ============================================================
# CELL 3 - APPLY LABELS, QC FLAGS, STATE ATTRIBUTES
# ============================================================
hr('CELL 3 - apply labels and per-cell attributes')

sub.obs[SUBTYPE_COL] = pd.Categorical(
    sub.obs[CLUSTER_KEY].map(LABEL_MAP).astype(str))

sub.obs['myeloid_qc_reason'] = sub.obs[CLUSTER_KEY].map(QC_FLAG).fillna('').astype(str)
sub.obs['qc_exclude_myeloid'] = sub.obs['myeloid_qc_reason'] != ''
sub.obs['is_myeloid'] = ~sub.obs[CLUSTER_KEY].isin(NON_MYELOID_CLUSTERS)

print('  Tier-2 myeloid subtype composition:')
comp = sub.obs[SUBTYPE_COL].value_counts()
print(comp.to_string())

print(f'\n  qc_exclude_myeloid: {int(sub.obs["qc_exclude_myeloid"].sum()):,} cells')
print(sub.obs.loc[sub.obs['qc_exclude_myeloid'], 'myeloid_qc_reason']
      .value_counts().to_string())
print(f'\n  not myeloid at all (is_myeloid False): '
      f'{int((~sub.obs["is_myeloid"]).sum()):,} cells '
      f'(clusters {NON_MYELOID_CLUSTERS})')

# state scores stay per-cell, never labels. Renamed to make that explicit.
for c in STATE_COLS:
    if c in sub.obs.columns:
        new = 'state_' + c.replace('sig_', '')
        sub.obs[new] = sub.obs[c].astype(float)
        print(f'  carried state attribute: {c} -> {new}')
    else:
        print(f'  WARN: state column absent: {c}')
for c in QC_SCORE_COLS:
    if c in sub.obs.columns:
        sub.obs[c.replace('qc_', 'score_')] = sub.obs[c].astype(float)

# the analysis-ready set: myeloid, not QC-excluded, subtype resolved
sub.obs['myeloid_analysis_ready'] = (
    sub.obs['is_myeloid'] & ~sub.obs['qc_exclude_myeloid'] &
    ~sub.obs[SUBTYPE_COL].astype(str).str.startswith(('EXCLUDE', 'Unresolved')))
print(f'\n  ANALYSIS-READY myeloid cells: '
      f'{int(sub.obs["myeloid_analysis_ready"].sum()):,} of {sub.n_obs:,}')
print(sub.obs.loc[sub.obs['myeloid_analysis_ready'], SUBTYPE_COL]
      .value_counts().to_string())


# %% ============================================================
# CELL 4 - COMPOSITION TABLES FOR THE DEG STEP
# ============================================================
hr('CELL 4 - subtype x animal x timepoint x tissue (DEG feasibility)')

ready = sub.obs[sub.obs['myeloid_analysis_ready']]
for cols, name in [([SUBTYPE_COL, TISSUE_COL], 'tissue'),
                   ([SUBTYPE_COL, TIME_COL], 'timepoint'),
                   ([SUBTYPE_COL, ANIMAL_COL], 'animal')]:
    if all(c in ready.columns for c in cols):
        ct = pd.crosstab(ready[cols[0]], ready[cols[1]].astype(str))
        print(f'\n  subtype x {name}:')
        print(ct.to_string())
        save_csv(ct, f'myeloid_subtype_by_{name}.csv')

if all(c in ready.columns for c in [SUBTYPE_COL, TISSUE_COL, TIME_COL, ANIMAL_COL]):
    full = (ready.groupby([SUBTYPE_COL, TISSUE_COL, TIME_COL, ANIMAL_COL],
                          observed=True).size().rename('n_cells').reset_index())
    save_csv(full, 'myeloid_subtype_full_breakdown.csv', index=False)
    # which subtypes can actually support a per-animal pseudobulk contrast
    feas = (full[full['n_cells'] >= 10]
            .groupby([SUBTYPE_COL, TISSUE_COL, TIME_COL], observed=True)[ANIMAL_COL]
            .nunique().rename('n_animals_ge10cells').reset_index())
    print('\n  animals with >=10 cells per subtype x tissue x timepoint')
    print('  (pseudobulk needs >=2 on BOTH sides of a contrast):')
    print(feas.to_string(index=False))
    save_csv(feas, 'myeloid_pseudobulk_feasibility.csv', index=False)
    thin = feas[feas['n_animals_ge10cells'] < 2]
    if len(thin):
        print(f'\n  *** {len(thin)} subtype x tissue x timepoint cell(s) have <2 '
              'animals with >=10 cells. Those cannot support pseudobulk and will')
        print('  appear in the DEG manifest as wilcoxon-only / pseudoreplicated.')


# %% ============================================================
# CELL 5 - FIGURES
# ============================================================
hr('CELL 5 - figures')

try:
    for color in [SUBTYPE_COL, 'qc_exclude_myeloid', 'myeloid_analysis_ready']:
        if color in sub.obs.columns:
            sc.pl.umap(sub, color=color, show=False,
                       save=f'_myeloid_annotated_{color}.pdf', legend_fontsize=14)
    print('  saved UMAP panels')
except Exception as e:
    print(f'  WARN: UMAP panels failed ({e})')

try:
    fig, ax = plt.subplots(figsize=(18, 11))
    d = comp[comp.index.map(lambda x: not str(x).startswith('EXCLUDE'))]
    ax.barh(range(len(d)), d.values, color='#2C6E9B', edgecolor='white', linewidth=2)
    ax.set_yticks(range(len(d)))
    ax.set_yticklabels(d.index, fontsize=22)
    ax.invert_yaxis()
    ax.set_xlabel('cells', fontsize=30)
    ax.set_title('Tier-2 myeloid subtypes', fontsize=34)
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)
    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(OUT_DIR, f'myeloid_subtype_counts.{ext}'),
                    dpi=FIG_DPI, bbox_inches='tight')
    plt.close(fig)
    print('  saved myeloid_subtype_counts.pdf / .png')
except Exception as e:
    print(f'  WARN: bar figure failed ({e})')


# %% ============================================================
# CELL 6 - PROPAGATE TO GLOBAL OBJECT (myeloid cells ONLY)
# ============================================================
hr('CELL 6 - propagate onto the global object')

if not os.path.exists(GLOBAL_IN):
    print(f'  WARN: global object not found ({GLOBAL_IN}); subset still written')
    glob_ok = False
else:
    glob_ok = True
    adata = sc.read_h5ad(GLOBAL_IN)
    print(f'  global: {adata.n_obs:,} cells')

    # Zoey's constraint, enforced: only Tier-1 myeloid cells may receive a
    # Tier-2 myeloid label. Anything else keeps whatever it already had.
    myeloid_mask = adata.obs[TIER1_COL].astype(str).isin(MYELOID_TIER1).values
    print(f'  Tier-1 myeloid in global: {int(myeloid_mask.sum()):,}')

    adata.obs[SUBTYPE_COL] = 'not myeloid'
    adata.obs['myeloid_qc_reason'] = ''
    adata.obs['qc_exclude_myeloid'] = False
    adata.obs['myeloid_analysis_ready'] = False

    shared = adata.obs_names.intersection(sub.obs_names)
    print(f'  obs_name overlap subset -> global: {len(shared):,} of {sub.n_obs:,}')
    if len(shared) == 0:
        raise RuntimeError('Zero obs_name overlap between subset and global. The '
                           'join would silently write nothing. Check that both '
                           'objects came from the same pipeline pass.')
    if len(shared) < sub.n_obs:
        print(f'  WARN: {sub.n_obs - len(shared):,} subset cells are absent from '
              'the global object and cannot be propagated')

    # Propagate, then RESTORE DTYPES explicitly. The .astype(object) round trip
    # leaves boolean columns holding Python bools in an object-dtype Series, and
    # h5py refuses those with "Can't implicitly convert non-string objects to
    # strings". That crashed the global write in testing AFTER the subset file had
    # already been written, which is the worst possible failure point because it
    # leaves a half-finished state on disk. Cast back before writing anything.
    STR_COLS  = [SUBTYPE_COL, 'myeloid_qc_reason']
    BOOL_COLS = ['qc_exclude_myeloid', 'myeloid_analysis_ready']
    for col in STR_COLS + BOOL_COLS:
        vals = adata.obs[col].astype(object).copy()
        vals.loc[shared] = sub.obs.loc[shared, col].values
        adata.obs[col] = vals
    for col in STR_COLS:
        adata.obs[col] = adata.obs[col].astype(str)
    for col in BOOL_COLS:
        adata.obs[col] = adata.obs[col].fillna(False).astype(bool)
    for c in [c for c in sub.obs.columns if c.startswith(('state_', 'score_'))]:
        v = pd.Series(np.nan, index=adata.obs_names, dtype=float)
        v.loc[shared] = sub.obs.loc[shared, c].astype(float).values
        adata.obs[c] = v

    # assert the constraint actually held
    wrote_to_non_myeloid = (
        (adata.obs[SUBTYPE_COL].astype(str) != 'not myeloid').values & ~myeloid_mask)
    n_bad = int(wrote_to_non_myeloid.sum())
    if n_bad:
        raise RuntimeError(f'{n_bad} non-myeloid cell(s) received a Tier-2 myeloid '
                           'label. That violates the Tier-1 constraint; aborting '
                           'before write.')
    print('  constraint verified: no non-myeloid cell carries a myeloid subtype')

    adata.obs[SUBTYPE_COL] = pd.Categorical(adata.obs[SUBTYPE_COL].astype(str))
    print('\n  global Tier-2 myeloid composition:')
    print(adata.obs[SUBTYPE_COL].value_counts().to_string())


# %% ============================================================
# CELL 7 - WRITE
# ============================================================
hr('CELL 7 - write')

# Validate BOTH objects before writing EITHER. The failure mode this prevents is
# the subset writing successfully and the global then crashing on an h5py dtype
# error, which leaves a half-finished annotation on disk that looks complete.
def obs_write_problems(ad_obj, tag):
    bad = []
    for c in ad_obj.obs.columns:
        s = ad_obj.obs[c]
        if s.dtype != object:
            continue
        kinds = set(type(v).__name__ for v in s.dropna().unique()[:200])
        if kinds - {'str'}:
            bad.append((tag, c, sorted(kinds)))
    return bad


problems = obs_write_problems(sub, 'subset')
if glob_ok:
    problems += obs_write_problems(adata, 'global')
if problems:
    print('  *** object-dtype obs columns that h5py cannot write:')
    for tag, c, kinds in problems:
        print(f'    {tag}.obs["{c}"] holds {kinds}')
    raise RuntimeError('Fix these dtypes before writing. Nothing was written, so '
                       'no half-finished state is on disk.')
print('  pre-write dtype check passed for all objects')

if DRY_RUN:
    print('  DRY_RUN is True; nothing written')
else:
    sub.write(SUBSET_OUT)
    print(f'  wrote {SUBSET_OUT}')
    if glob_ok:
        adata.write(GLOBAL_OUT)
        print(f'  wrote {GLOBAL_OUT}')

n_ready = int(sub.obs['myeloid_analysis_ready'].sum())
print(f"""
  SUMMARY
  -------
  myeloid subset cells       : {sub.n_obs:,}
  labelled subtypes          : {sub.obs[SUBTYPE_COL].nunique()}
  excluded as non-myeloid    : {int((~sub.obs['is_myeloid']).sum()):,}  (T/NK, platelets)
  QC-flagged (kept, flagged) : {int(sub.obs['qc_exclude_myeloid'].sum()):,}
  ANALYSIS-READY for DEG     : {n_ready:,}

  NEXT: shiv_deg_comparisons.py picks these up automatically. Its TIER2_COLS
  already lists '{SUBTYPE_COL}' first, and it honours qc_exclude. Set
  QC_COL = 'qc_exclude_myeloid' there if you want the myeloid flags respected in
  addition to the existing global qc_exclude, and run the MANIFEST_ONLY pass
  first to see which myeloid contrasts survive the cell-count gates.
""")
print('=' * 90)
print('TIER-2 MYELOID ANNOTATION COMPLETE')
print('=' * 90)
