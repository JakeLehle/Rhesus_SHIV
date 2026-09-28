#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Rhesus SHIV — Comprehensive DEG Comparison Suite
=================================================
Runs Zoey's full DEG manifest across tissues, timepoints, animals, subspecies and
Tier-2 cell populations, writes one DEG table plus top-50 up/down per contrast,
and emits the threshold summary that decides the GSEA cutoff.

BUILT TO BE COMPREHENSIVE (Jake's call, 2026-09): every population is run in all
three tissue scopes (PBMC, LN, pooled) even where Zoey's document was silent on
tissue, plus the cross-tissue contrasts, plus a batch-controlled variant of the
subspecies contrast. Extra contrasts are labelled in the manifest so the ones she
explicitly asked for can be pulled out with a filter.

SIX REQUESTED CONTRASTS DO NOT EXIST IN THE DATA
------------------------------------------------
Built from ANIMAL_META in scanpy_analysis_script.py:

  tissue  39272            40702/40707/41861   41862
  PBMC    Pre, 21 DPI      all three           Pre, Necropsy
  LN      Pre, 21 DPI      all three           all three

39272 has no necropsy at all (flagged 'no necropsy' in ANIMAL_INFO). 41862 has no
21 DPI PBMC, because its 21 DPI sample sits in library F, which is lymph node. So:

  PBMC 39272  21 DPI vs Necropsy   IMPOSSIBLE
  PBMC 39272  Pre    vs Necropsy   IMPOSSIBLE
  PBMC 41862  Pre    vs 21 DPI     IMPOSSIBLE
  PBMC 41862  21 DPI vs Necropsy   IMPOSSIBLE
  LN   39272  21 DPI vs Necropsy   IMPOSSIBLE
  LN   39272  Pre    vs Necropsy   IMPOSSIBLE

These are emitted in the manifest with runnable=False and a reason, never
silently dropped, so the gap is visible to Zoey rather than looking like an
analysis that quietly returned fewer rows.

The grouped contrasts are also UNBALANCED and not paired end to end: PBMC is
Pre n=5, 21 DPI n=4, Necropsy n=4 with a DIFFERENT animal missing at each end.
LN is Pre 5, 21 DPI 5, Necropsy 4. The manifest records n per side so this is
visible in every output row.

TWO EVIDENCE LAYERS, AND WHY THE DISTINCTION MATTERS FOR THE THRESHOLD
-----------------------------------------------------------------------
  pseudobulk  per-animal summed counts -> logCPM -> paired or Welch t-test.
              This is the credible layer (standing lab convention). Needs
              >= MIN_ANIMALS_PER_SIDE animals on BOTH sides.
  wilcoxon    cell-level rank_genes_groups. Fast, sensitive, and PSEUDOREPLICATED
              whenever a side has fewer than 2 animals.

Every individual-animal contrast is n=1 vs n=1 by construction, so it has NO
replication and only the wilcoxon layer can run. Those rows carry
pseudoreplicated=True in the output table itself, not just in conversation,
because 24 such contrasts will otherwise read like 24 independent findings.

This matters directly for the threshold question. If raw-p and FDR counts are
pooled across grouped and individual contrasts, the pseudoreplicated ones will
dominate the totals (cell-level n is in the thousands, so p values are tiny) and
push the chosen cutoff too permissive. The threshold table is therefore SPLIT by
evidence layer, and the recommendation is read off the pseudobulk layer.

Read-only w.r.t. the input object. Spyder cells (# %%). Headless under SLURM.

Author: Jake Lehle / Kaushal Lab
Date: September 2026
"""

# %% ============================================================
# CELL 1 — CONFIG AND DIALS
# ============================================================
import os
os.environ.setdefault('MPLBACKEND', 'Agg')
import matplotlib
matplotlib.use('Agg')
import itertools
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from scipy import sparse, stats
import warnings
warnings.filterwarnings("ignore")

WORKING_DIR = '/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV'
MERGED      = os.path.join(WORKING_DIR, 'shiv_starsolo', 'merged')

# Prefer the object carrying BOTH Tier-2 T and Tier-2 myeloid labels.
ADATA_CANDIDATES = [
    os.path.join(MERGED, 'shiv_host_tier2_myeloid_annotated.h5ad'),
    os.path.join(MERGED, 'shiv_host_tier2_annotated.h5ad'),
    os.path.join(MERGED, 'shiv_host_final_tcrTrecovered.h5ad'),
]

OUT_DIR   = os.path.join(WORKING_DIR, 'annotation_output', 'deg_comparisons')
TAB_DIR   = os.path.join(OUT_DIR, 'tables')
ORA_DIR   = os.path.join(OUT_DIR, 'ora_gene_lists')
TOP_DIR   = os.path.join(OUT_DIR, 'top50')

# ---- two-stage running ----
# 176 of the 288 runnable contrasts are individual-animal, n=1 vs n=1, and every
# one is pseudoreplicated by construction. They are 61% of the runnable set and
# they will dominate the wall clock while contributing nothing to the threshold
# decision, which must be read off the pseudobulk layer.
#
# So: first pass with CREDIBLE_ONLY = True runs the 97 pseudobulk-capable
# contrasts, which is what the GSEA threshold is chosen from and what goes to
# Zoey. Review DEG_counts_primary_table.csv and DEG_lfc_sweep.csv, then set it
# False for the complete run including the pseudoreplicated ones.
CREDIBLE_ONLY = True

# ---- RUN THIS FIRST ----
# True  = build and write the manifest, print the plan, run NO tests. Review the
#         manifest with Zoey before spending the compute.
# False = full run.
MANIFEST_ONLY = True

# ---- columns ----
ANIMAL_COL = 'animal_id'
TIME_COL   = 'timepoint'
TISSUE_COL = 'tissue'
LIB_COL    = 'library'
SUBSP_COL  = 'subspecies'
TIER1_COL  = 'tier1_celltype'
# Tier-2 label columns to harvest populations from, in priority order. A cell's
# population is the first non-null, non-placeholder value across these.
TIER2_COLS = ['tier2_myeloid_subtype', 'tier2_celltype', 'tier2_subtype']

# ---- placeholder handling (FIXED 2026-09-26) ----
# The myeloid annotation writes 'not myeloid' onto every non-myeloid cell, which
# is 80,146 of 89,844. With the original hard-coded placeholder list, that value
# was treated as a real label, so:
#   1. tier2_myeloid_subtype claimed 100% of cells and tier2_celltype NEVER got
#      to fill anything, which silently DELETED all four T cell subsets Zoey
#      asked for (Naive/Resting CD4+, CD4+ CM, CD8+ EM, CD8+ Terminal Effector).
#   2. 'not myeloid' itself became an 80k-cell "population" and would have had a
#      full ~62-contrast set run on it, pooling T, B, NK, platelets and
#      erythrocytes into one meaningless comparison, at the highest cost of
#      anything in the run.
#   3. 'EXCLUDE: T/NK contamination', 'EXCLUDE: Platelets' and 'Unresolved
#      myeloid' would each have got a full contrast set too.
# SKIP_VALUES pass through to the next column in TIER2_COLS. DENY_PREFIXES mark a
# cell unlabelled so it is dropped from every population.
SKIP_VALUES   = {'nan', 'None', '', 'NA', 'Unassigned', 'not T', 'not myeloid'}
DENY_PREFIXES = ('EXCLUDE:', 'Unresolved')

# QC columns, all ORed. The myeloid annotation adds its own flag, so both are
# honoured rather than only the global one.
QC_COLS = ['qc_exclude', 'qc_exclude_myeloid']
COUNTS_LAYER = 'counts'

# ---- cohort ----
INFECTED_ANIMALS = ['39272', '40702', '40707', '41861', '41862']   # Chinese
CONTROL_ANIMALS  = ['34315', '41903']                              # Indian
TIMEPOINTS       = ['Pre', '21 DPI', 'Necropsy']
TP_PAIRS         = [('Pre', '21 DPI'), ('21 DPI', 'Necropsy'), ('Pre', 'Necropsy')]
CONTROL_TP       = 'Non-Infected'
TISSUES          = ['PBMC', 'LN']

# ---- gates ----
MIN_CELLS_PER_SIDE   = 30     # below this a contrast is skipped with a reason
MIN_ANIMALS_PER_SIDE = 2      # pseudobulk requires this on BOTH sides
MIN_CELLS_PER_ANIMAL = 10     # an animal contributes a pseudobulk sample above this
MIN_GENE_CELLS       = 10     # gene must be detected in this many cells to be tested
DUPLICATE_TOLERANCE  = 0.10   # pooled vs single-tissue cell counts within 10% = same comparison

N_TOP                = 50     # top N up and down

# ---- populations ----
RUN_ALL_CELLS   = True        # the 'ALL' population (every cell in scope)
MIN_POP_CELLS   = 100         # a Tier-2 label must have this many cells overall
# Auto-discovery was right while the annotation was in flux. Now that both
# Tier-2 passes are locked, NAME the populations. It makes the hand-off auditable
# and it stops a stray label value from quietly acquiring a 62-contrast set.
# Set to None to go back to auto-discovery.
EXPLICIT_POPULATIONS = [
    # --- myeloid Tier 2 (shiv_tier2_myeloid_annotation_write.py) ---
    'Classical Monocyte',
    'Non-classical / CD163+ Mono-Mac',
    'Tissue Macrophage (C1Q+/CD206+)',
    'Antigen Presenting DC',
    'pDC',
    'Migratory / Mature DC',
    'Monocyte, MHC-II-high',
    'Monocyte, IFN-stimulated',
    # --- T cell Tier 2 (shiv_tier2_annotation_write.py) ---
    # Zoey's four requested subsets. Confirm these strings against
    # adata.obs['tier2_celltype'].value_counts() on the first run; the label-map
    # docstring in that script uses these names.
    'CD8 effector memory',
    'CD8 terminal effector',
    'CD4 central memory',
    'Naive/resting T',
]

# ---- tissue scopes ----
TISSUE_SCOPES = ['PBMC', 'LN', 'Both']

# ---- extras beyond Zoey's document ----
INCLUDE_CROSS_TISSUE      = True   # PBMC vs LN at each timepoint
INCLUDE_SPECIES_CROSS_TIS = True   # Chinese PBMC vs Indian LN and the reverse
INCLUDE_SPECIES_WITHIN_LIB = True  # batch-controlled subspecies contrast, below

# The subspecies contrast is NOT fully confounded with batch, which is worth
# exploiting. Library C carries both Indian animals AND Chinese 41862 in PBMC;
# library F carries both Indian animals AND Chinese 41861/41862 in LN. So a
# within-library contrast exists. Using all five Chinese pulls four of them from
# libraries A and B and puts batch straight back in. We run BOTH and report them
# side by side rather than picking silently.
#
# IMPORTANT ASYMMETRY, found in testing. The Chinese animals inside each mixed
# library are NOT at the same timepoint:
#   library C (PBMC) holds Chinese 41862 at Pre      -> Pre vs Non-Infected
#   library F (LN)   holds Chinese 41861/41862 at 21 DPI, NOT Pre
# So the batch-controlled LN contrast cannot be "Chinese Pre vs Indian" at all.
# Asking for Pre in library F returns zero cells. The nearest available LN
# contrast is Chinese 21 DPI vs Indian Non-Infected, which answers a different
# question (infected vs uninfected, not subspecies) and is labelled as such.
SPECIES_WITHIN_LIB = {
    'PBMC': {'library': 'C', 'chinese_tp': 'Pre',
             'note': 'batch-controlled; Chinese side is 41862 only, n=1'},
    'LN':   {'library': 'F', 'chinese_tp': '21 DPI',
             'note': 'batch-controlled BUT Chinese side is 21 DPI not Pre, so '
                     'this is infected-vs-uninfected, NOT a clean subspecies test'},
}

# ---- thresholds under evaluation ----
# Zoey's requested output (2026-09-26): one table of DEG counts at
# |log2FC| > PRIMARY_LFC, with a raw-p column and an FDR column, so the GSEA
# filter can be chosen from data rather than convention.
PRIMARY_ALPHA = 0.05
PRIMARY_LFC   = 1.5
# The other levels are still computed and written, because |log2FC| > 1.5 is a
# hard cut and at n=2-4 animals per side it may return near zero. If it does, the
# neighbouring levels are what tell you whether the threshold or the power is the
# limiting factor. They cost nothing to carry.
ALPHA_LEVELS = [0.05, 0.10, 0.25]
LFC_LEVELS   = [0.0, 0.5, 1.0, 1.5, 2.0]

# Over-representation analysis (enrichKEGG / enrichGO / DAVID / Enrichr) needs a
# GENE LIST, not a ranked vector, and directional pathways need the up and down
# lists tested SEPARATELY. So the counts below are split by direction, and a
# contrast is only enrichable if BOTH directions clear a floor.
MIN_GENES_FOR_ORA = 10        # per direction; below this an ORA is not worth running

FIG_DPI = 300
HEX = {'pb': '#2C6E9B', 'wx': '#C4453C', 'grid': '#D8D8D8', 'acc': '#E8A33D'}
plt.rcParams.update({
    'font.size': 28, 'axes.titlesize': 34, 'axes.labelsize': 30,
    'xtick.labelsize': 24, 'ytick.labelsize': 24, 'legend.fontsize': 20,
    'pdf.fonttype': 42, 'ps.fonttype': 42,
})
for d in (OUT_DIR, TAB_DIR, TOP_DIR, ORA_DIR):
    os.makedirs(d, exist_ok=True)
pd.set_option('display.width', 240)
pd.set_option('display.max_columns', 80)
pd.set_option('display.max_rows', 500)


def hr(t):
    print('\n' + '=' * 92)
    print(t)
    print('=' * 92)


def slug(s):
    out = str(s)
    for a, b in [(' ', '_'), ('/', '-'), ('+', 'pos'), ('(', ''), (')', ''),
                 (',', ''), ('.', ''), (':', ''), ('*', '')]:
        out = out.replace(a, b)
    return out


def bh_fdr(p):
    """Benjamini-Hochberg. NaN-safe, returns array aligned to input."""
    p = np.asarray(p, dtype=float)
    out = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if not ok.any():
        return out
    pv = p[ok]
    n = pv.size
    order = np.argsort(pv)
    ranked = pv[order]
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    res = np.empty(n)
    res[order] = q
    out[ok] = res
    return out


# %% ============================================================
# CELL 2 — LOAD AND VALIDATE
# ============================================================
hr('CELL 2 — load object and validate metadata')

ADATA_IN = next((c for c in ADATA_CANDIDATES if os.path.exists(c)), None)
if ADATA_IN is None:
    raise FileNotFoundError('No input object. Looked for:\n  ' +
                            '\n  '.join(ADATA_CANDIDATES))
print(f'  input: {os.path.basename(ADATA_IN)}')

adata = sc.read_h5ad(ADATA_IN)
print(f'  loaded: {adata.n_obs:,} cells x {adata.n_vars:,} genes')

for c in [ANIMAL_COL, TIME_COL, TISSUE_COL]:
    if c not in adata.obs.columns:
        raise KeyError(f'required obs column missing: {c}')
for c in [LIB_COL, SUBSP_COL]:
    if c not in adata.obs.columns:
        print(f'  WARN: {c} absent; contrasts needing it will be skipped')

qc_present = [c for c in QC_COLS if c in adata.obs.columns]
if qc_present:
    before = adata.n_obs
    drop = np.zeros(adata.n_obs, dtype=bool)
    for c in qc_present:
        d = adata.obs[c].fillna(False).astype(bool).values
        print(f'  {c}: {int(d.sum()):,} flagged')
        drop |= d
    adata = adata[~drop].copy()
    print(f'  QC filter ({" OR ".join(qc_present)}): {before:,} -> {adata.n_obs:,} cells')
else:
    print(f'  NOTE: none of {QC_COLS} present; no QC filtering applied')

# ---- verify .X is lognormalized before the cell-level engine trusts it ----
# The wilcoxon engine reads .X directly. If .X holds raw counts the test still
# runs and still returns p values, but the log fold changes are computed off the
# wrong scale and are silently wrong. Caught in testing against a counts-in-X
# object, so we check rather than assume.
DE_EXPR_LAYER = None
try:
    _probe = adata.X[:2000]
    _probe = _probe.toarray() if sparse.issparse(_probe) else np.asarray(_probe)
    _looks_like_counts = bool(np.all(np.equal(np.mod(_probe, 1), 0))) and _probe.max() > 30
except Exception as e:
    print(f'  WARN: could not probe .X ({e}); assuming lognormalized')
    _looks_like_counts = False

if _looks_like_counts:
    print('\n  *** .X looks like RAW COUNTS (all-integer, max > 30).')
    print('      The cell-level wilcoxon layer needs lognormalized values or its')
    print('      log2FC column is wrong. Building a lognorm layer for DE use.')
    print('      The input object on disk is NOT modified.')
    src_layer = COUNTS_LAYER if COUNTS_LAYER in adata.layers else None
    tmpX = adata.layers[src_layer].copy() if src_layer else adata.X.copy()
    _t = sc.AnnData(X=tmpX)
    sc.pp.normalize_total(_t, target_sum=1e4)
    sc.pp.log1p(_t)
    adata.layers['_lognorm_for_de'] = _t.X
    DE_EXPR_LAYER = '_lognorm_for_de'
    del _t, tmpX
    print('      built layer "_lognorm_for_de"')
else:
    print('  .X passes the lognorm sanity check (non-integer values present)')

obs = adata.obs
obs[ANIMAL_COL] = obs[ANIMAL_COL].astype(str)
obs[TIME_COL] = obs[TIME_COL].astype(str)
obs[TISSUE_COL] = obs[TISSUE_COL].astype(str)

print('\n  animal x tissue x timepoint cell counts:')
av = (obs.groupby([TISSUE_COL, ANIMAL_COL, TIME_COL]).size()
      .rename('n_cells').reset_index())
print(av.to_string(index=False))
av.to_csv(os.path.join(OUT_DIR, 'cell_counts_by_animal_tissue_timepoint.csv'),
          index=False)

# ---- population column ----
pop = pd.Series(index=obs.index, dtype=object)
used_cols = []
for c in TIER2_COLS:
    if c not in obs.columns:
        continue
    used_cols.append(c)
    v = obs[c].astype(str)
    skip = v.isin(SKIP_VALUES)
    deny = v.str.startswith(DENY_PREFIXES)
    filled_before = int(pop.notna().sum())
    # denied values are recorded as a sentinel so a later column cannot revive
    # them, then stripped below
    cand = v.where(~skip)
    cand = cand.mask(deny, '__DENIED__')
    pop = pop.where(pop.notna(), cand)
    gained = int(pop.notna().sum()) - filled_before
    n_skip, n_deny = int(skip.sum()), int(deny.sum())
    print(f'  {c}: filled {gained:,} cells  (skipped {n_skip:,} placeholder, '
          f'denied {n_deny:,})')
    if n_skip:
        print(f'      skipped values: '
              f'{sorted(set(v[skip].unique()) & set(SKIP_VALUES))}')
    if n_deny:
        print(f'      denied values: {sorted(set(v[deny].unique()))}')
pop = pop.replace('__DENIED__', np.nan)
if not used_cols:
    print(f'  WARN: none of {TIER2_COLS} present; falling back to {TIER1_COL}')
    if TIER1_COL in obs.columns:
        pop = obs[TIER1_COL].astype(str)
adata.obs['_population'] = pop.fillna('Unlabelled')
print(f'\n  population column built from: {used_cols}')
counts = adata.obs['_population'].value_counts()
print(counts.to_string())

if EXPLICIT_POPULATIONS:
    POPULATIONS = [p for p in EXPLICIT_POPULATIONS if p in set(counts.index)]
    missing_pop = [p for p in EXPLICIT_POPULATIONS if p not in set(counts.index)]
    if missing_pop:
        print(f'  WARN: requested populations absent: {missing_pop}')
else:
    POPULATIONS = [p for p, n in counts.items()
                   if n >= MIN_POP_CELLS and p not in ('Unlabelled',)]
if RUN_ALL_CELLS:
    POPULATIONS = ['ALL'] + POPULATIONS
print(f'\n  populations to test ({len(POPULATIONS)}): {POPULATIONS}')


# %% ============================================================
# CELL 3 — BUILD THE CONTRAST MANIFEST
# ============================================================
hr('CELL 3 — build contrast manifest')


def cells_for(population, tissue_scope, animals=None, timepoints=None,
              library=None, subspecies=None):
    """Boolean mask over adata for one side of a contrast."""
    m = np.ones(adata.n_obs, dtype=bool)
    if population != 'ALL':
        m &= (adata.obs['_population'].astype(str) == population).values
    if tissue_scope != 'Both':
        m &= (adata.obs[TISSUE_COL] == tissue_scope).values
    if animals is not None:
        m &= adata.obs[ANIMAL_COL].isin(animals).values
    if timepoints is not None:
        m &= adata.obs[TIME_COL].isin(timepoints).values
    if library is not None and LIB_COL in adata.obs.columns:
        m &= (adata.obs[LIB_COL].astype(str) == library).values
    if subspecies is not None and SUBSP_COL in adata.obs.columns:
        m &= (adata.obs[SUBSP_COL].astype(str) == subspecies).values
    return m


def side_spec(**kw):
    return kw


manifest = []


def add(population, scope, family, name, a_spec, b_spec, requested_by_zoey,
        note=''):
    manifest.append({
        'contrast_id': None, 'population': population, 'tissue_scope': scope,
        'family': family, 'contrast': name,
        'side_a': a_spec, 'side_b': b_spec,
        'requested': requested_by_zoey, 'note': note,
    })


for population in POPULATIONS:
    for scope in TISSUE_SCOPES:
        # --- grouped timepoint pairs (Chinese infected) ---
        for x, y in TP_PAIRS:
            add(population, scope, 'grouped_timepoint', f'{x}_vs_{y}',
                side_spec(animals=INFECTED_ANIMALS, timepoints=[x]),
                side_spec(animals=INFECTED_ANIMALS, timepoints=[y]),
                True)
        # --- per-animal timepoint pairs ---
        for a in INFECTED_ANIMALS:
            for x, y in TP_PAIRS:
                add(population, scope, 'individual_timepoint',
                    f'{a}_{x}_vs_{y}',
                    side_spec(animals=[a], timepoints=[x]),
                    side_spec(animals=[a], timepoints=[y]),
                    True)
        # --- subspecies, all animals (batch-confounded) ---
        add(population, scope, 'subspecies',
            'ChinesePre_vs_IndianNonInf',
            side_spec(animals=INFECTED_ANIMALS, timepoints=['Pre']),
            side_spec(animals=CONTROL_ANIMALS, timepoints=[CONTROL_TP]),
            True, note='batch-confounded: Chinese mostly lib A/B, Indian lib C/F')
        # --- subspecies, within-library (batch-controlled) ---
        if INCLUDE_SPECIES_WITHIN_LIB and scope in SPECIES_WITHIN_LIB:
            cfg = SPECIES_WITHIN_LIB[scope]
            lib, ctp = cfg['library'], cfg['chinese_tp']
            add(population, scope, 'subspecies_within_library',
                f'Chinese{slug(ctp)}_vs_IndianNonInf_lib{lib}',
                side_spec(animals=INFECTED_ANIMALS, timepoints=[ctp], library=lib),
                side_spec(animals=CONTROL_ANIMALS, timepoints=[CONTROL_TP], library=lib),
                False, note=cfg['note'])

    # --- cross-tissue contrasts (scope is inherently Both) ---
    if INCLUDE_CROSS_TISSUE:
        for tp in TIMEPOINTS:
            add(population, 'Both', 'cross_tissue', f'PBMC_vs_LN_{slug(tp)}',
                side_spec(animals=INFECTED_ANIMALS, timepoints=[tp]),
                side_spec(animals=INFECTED_ANIMALS, timepoints=[tp]),
                True, note='side A restricted to PBMC, side B to LN')
    if INCLUDE_SPECIES_CROSS_TIS:
        add(population, 'Both', 'cross_tissue_subspecies',
            'ChinesePre_PBMC_vs_Indian_LN',
            side_spec(animals=INFECTED_ANIMALS, timepoints=['Pre']),
            side_spec(animals=CONTROL_ANIMALS, timepoints=[CONTROL_TP]),
            True, note='confounds tissue AND subspecies; interpret with care')
        add(population, 'Both', 'cross_tissue_subspecies',
            'ChinesePre_LN_vs_Indian_PBMC',
            side_spec(animals=INFECTED_ANIMALS, timepoints=['Pre']),
            side_spec(animals=CONTROL_ANIMALS, timepoints=[CONTROL_TP]),
            True, note='confounds tissue AND subspecies; interpret with care')

# cross-tissue families need an explicit per-side tissue override
for row in manifest:
    if row['family'] == 'cross_tissue':
        row['side_a'] = {**row['side_a'], '_tissue': 'PBMC'}
        row['side_b'] = {**row['side_b'], '_tissue': 'LN'}
    elif row['family'] == 'cross_tissue_subspecies':
        if row['contrast'].startswith('ChinesePre_PBMC'):
            row['side_a'] = {**row['side_a'], '_tissue': 'PBMC'}
            row['side_b'] = {**row['side_b'], '_tissue': 'LN'}
        else:
            row['side_a'] = {**row['side_a'], '_tissue': 'LN'}
            row['side_b'] = {**row['side_b'], '_tissue': 'PBMC'}

for i, row in enumerate(manifest):
    row['contrast_id'] = f'C{i:05d}'

print(f'  contrasts constructed: {len(manifest)}')


def resolve_side(population, scope, spec):
    s = dict(spec)
    tissue = s.pop('_tissue', None) or scope
    return cells_for(population, tissue, **s), tissue


# ---- annotate runnability ----
rows = []
for row in manifest:
    ma, ta = resolve_side(row['population'], row['tissue_scope'], row['side_a'])
    mb, tb = resolve_side(row['population'], row['tissue_scope'], row['side_b'])
    na, nb = int(ma.sum()), int(mb.sum())
    aa = sorted(set(adata.obs.loc[ma, ANIMAL_COL])) if na else []
    ab = sorted(set(adata.obs.loc[mb, ANIMAL_COL])) if nb else []

    runnable, reason = True, ''
    if na == 0 or nb == 0:
        runnable = False
        empty = 'A' if na == 0 else ('B' if nb == 0 else 'both')
        reason = f'no cells on side {empty} (this combination does not exist)'
    elif na < MIN_CELLS_PER_SIDE or nb < MIN_CELLS_PER_SIDE:
        runnable = False
        reason = f'below MIN_CELLS_PER_SIDE ({na} vs {nb} < {MIN_CELLS_PER_SIDE})'

    n_anim_a = len([a for a in aa
                    if (adata.obs.loc[ma, ANIMAL_COL] == a).sum() >= MIN_CELLS_PER_ANIMAL])
    n_anim_b = len([a for a in ab
                    if (adata.obs.loc[mb, ANIMAL_COL] == a).sum() >= MIN_CELLS_PER_ANIMAL])
    pb_ok = runnable and n_anim_a >= MIN_ANIMALS_PER_SIDE and n_anim_b >= MIN_ANIMALS_PER_SIDE
    pseudorep = not pb_ok

    rows.append({**{k: v for k, v in row.items() if k not in ('side_a', 'side_b')},
                 'tissue_a': ta, 'tissue_b': tb,
                 'n_cells_a': na, 'n_cells_b': nb,
                 'animals_a': ';'.join(aa), 'animals_b': ';'.join(ab),
                 'n_animals_a': n_anim_a, 'n_animals_b': n_anim_b,
                 'runnable': runnable, 'skip_reason': reason,
                 'pseudobulk_possible': pb_ok,
                 'pseudoreplicated': pseudorep,
                 # power, made visible per row: a pseudobulk t-test on 2 vs 3
                 # animals has 3 residual df, and FDR across ~12,000 genes at
                 # that df is almost always empty. This column is why.
                 'min_animals_per_side': min(n_anim_a, n_anim_b),
                 'pseudobulk_resid_df': (n_anim_a + n_anim_b - 2) if pb_ok else 0})

man_df = pd.DataFrame(rows)

# ---- near-duplicate detection (added 2026-09-26) ----
# The pooled "Both" scope is a genuinely different comparison only when the
# population actually occupies both tissues. For a population that lives in one
# compartment - Classical Monocyte is 1,067 PBMC cells against 1 in lymph node -
# the Both contrast is the single-tissue contrast run again under another name.
# In the first credible list that was 16 of 51 Both contrasts, so 97 credible
# contrasts were really 81 distinct ones. Left unflagged they are counted twice
# in every median in the threshold table, which biases the GSEA decision toward
# whichever populations happen to be tissue-restricted.
man_df['near_duplicate_of'] = ''
for (popn, con), g in man_df[man_df['runnable']].groupby(['population', 'contrast']):
    both = g[g['tissue_scope'] == 'Both']
    single = g[g['tissue_scope'] != 'Both']
    for bi, b in both.iterrows():
        for _, sg in single.iterrows():
            da = abs(b['n_cells_a'] - sg['n_cells_a']) / max(b['n_cells_a'], 1)
            db = abs(b['n_cells_b'] - sg['n_cells_b']) / max(b['n_cells_b'], 1)
            if da < DUPLICATE_TOLERANCE and db < DUPLICATE_TOLERANCE:
                man_df.at[bi, 'near_duplicate_of'] = sg['contrast_id']
                break

man_df.to_csv(os.path.join(OUT_DIR, 'contrast_manifest.csv'), index=False)
print(f'  wrote contrast_manifest.csv')
n_dup = int((man_df['near_duplicate_of'] != '').sum())
if n_dup:
    print(f'\n  NEAR-DUPLICATES: {n_dup} pooled-tissue contrast(s) reproduce a '
          'single-tissue contrast')
    print('  of the same population, because that population lives in one tissue.')
    print('  They still run (the pooling is not wrong), but they are flagged in')
    print('  near_duplicate_of and excluded from the threshold medians.')
    dd = man_df[man_df['near_duplicate_of'] != '']
    print(dd.groupby('population').size().rename('n_duplicate_contrasts')
          .to_string())

print('\n  manifest summary by family (pseudobulk/pseudoreplicated counted among')
print('  RUNNABLE contrasts only, so the columns add up):')
_run = man_df[man_df['runnable']]
summ = (man_df.groupby('family').agg(total=('contrast_id', 'size'),
                                     runnable=('runnable', 'sum')))
summ = summ.join(_run.groupby('family').agg(
    pseudobulk=('pseudobulk_possible', 'sum'),
    pseudoreplicated=('pseudoreplicated', 'sum'))).fillna(0).astype(int).reset_index()
print(summ.to_string(index=False))

# REPORTING FIX (2026-09-26). This previously de-duplicated across populations
# and labelled the result "population-independent", which was wrong and actively
# misleading: a contrast appeared if ANY ONE population lacked cells, so
# Classical Monocyte (1 LN cell) made the LN timepoint contrasts look impossible
# and Tissue Macrophage (1 PBMC cell) did the same for PBMC, when both are fine
# for most populations. Report two separate things instead.
impossible = man_df[(~man_df['runnable']) &
                    (man_df['skip_reason'].str.contains('does not exist'))]
n_pop = man_df['population'].nunique()
if len(impossible):
    print(f'\n  {len(impossible)} contrast x population combination(s) skipped '
          'because that combination has no cells.')

    # (a) contrasts impossible for EVERY population -> a genuine data gap
    per = (impossible.groupby(['tissue_scope', 'family', 'contrast'])
           .size().rename('n_populations_affected').reset_index())
    universal = per[per['n_populations_affected'] >= n_pop]
    print(f'\n  (a) impossible for ALL {n_pop} populations '
          f'(a real gap in the study design): {len(universal)}')
    if len(universal):
        print(universal[['tissue_scope', 'family', 'contrast']]
              .sort_values(['tissue_scope', 'contrast']).to_string(index=False))
    else:
        print('      none')

    # (b) population-specific losses -> a cell-count problem for that population
    partial = per[per['n_populations_affected'] < n_pop]
    print(f'\n  (b) lost for SOME populations only (cell-count limited), '
          f'top 20 by breadth:')
    print(partial.sort_values('n_populations_affected', ascending=False)
          .head(20).to_string(index=False))

    # (c) which populations lose the most, so the thin ones are obvious
    bypop = (impossible.groupby('population').size()
             .rename('contrasts_lost').reset_index()
             .merge(man_df.groupby('population').size().rename('contrasts_total')
                    .reset_index(), on='population'))
    bypop['pct_lost'] = (100 * bypop['contrasts_lost'] /
                         bypop['contrasts_total']).round(1)
    print('\n  (c) contrasts lost per population:')
    print(bypop.sort_values('pct_lost', ascending=False).to_string(index=False))
    bypop.to_csv(os.path.join(OUT_DIR, 'contrasts_lost_by_population.csv'),
                 index=False)

# the runnable set, which is what actually goes to Zoey
run_df = man_df[man_df['runnable']].copy()
cred = run_df[run_df['pseudobulk_possible']]
print(f'\n  RUNNABLE: {len(run_df)} contrasts, of which {len(cred)} support '
      'pseudobulk (the credible layer)')
print('\n  credible (pseudobulk-capable) contrasts per population x family:')
if len(cred):
    print(pd.crosstab(cred['population'], cred['family']).to_string())
    cred[['contrast_id', 'population', 'tissue_scope', 'family', 'contrast',
          'n_cells_a', 'n_cells_b', 'n_animals_a', 'n_animals_b']].to_csv(
        os.path.join(OUT_DIR, 'credible_contrasts.csv'), index=False)
    print('\n  wrote credible_contrasts.csv  <- this is the list for Zoey')

if MANIFEST_ONLY:
    print('\n' + '*' * 92)
    print('  MANIFEST_ONLY is True. No tests were run.')
    print('  Review contrast_manifest.csv, then set MANIFEST_ONLY = False.')
    print('*' * 92)
    raise SystemExit(0)


# %% ============================================================
# CELL 4 — DE ENGINES
# ============================================================
hr('CELL 4 — DE engines')


def _counts_matrix(mask):
    sl = adata[mask]
    if COUNTS_LAYER in adata.layers:
        X = sl.layers[COUNTS_LAYER]
    else:
        X = sl.X
    return X


def pseudobulk_de(mask_a, mask_b, genes_keep):
    """Per-animal summed counts -> logCPM -> paired or Welch t-test.

    Paired when the animal sets match exactly (the timepoint contrasts, where the
    same animals appear on both sides). Otherwise Welch. Which one ran is
    returned so it lands in the output table.
    """
    def build(mask):
        an = adata.obs.loc[mask, ANIMAL_COL].values
        X = _counts_matrix(mask)
        X = X.tocsr() if sparse.issparse(X) else sparse.csr_matrix(X)
        out = {}
        for a in np.unique(an):
            sel = (an == a)
            if sel.sum() < MIN_CELLS_PER_ANIMAL:
                continue
            out[a] = np.asarray(X[sel].sum(axis=0)).ravel()
        return out

    A, B = build(mask_a), build(mask_b)
    if len(A) < MIN_ANIMALS_PER_SIDE or len(B) < MIN_ANIMALS_PER_SIDE:
        return None, 'insufficient animals'

    shared = sorted(set(A) & set(B))
    paired = (sorted(A) == sorted(B)) and len(shared) >= MIN_ANIMALS_PER_SIDE

    def logcpm(d, keys):
        M = np.vstack([d[k] for k in keys]).astype(float)
        tot = M.sum(axis=1, keepdims=True)
        tot[tot == 0] = 1.0
        return np.log2(M / tot * 1e6 + 1.0)

    if paired:
        ka = kb = shared
        La, Lb = logcpm(A, ka), logcpm(B, kb)
        t, p = stats.ttest_rel(La[:, genes_keep], Lb[:, genes_keep], axis=0)
        method = f'pseudobulk_paired_t(n={len(shared)})'
    else:
        ka, kb = sorted(A), sorted(B)
        La, Lb = logcpm(A, ka), logcpm(B, kb)
        t, p = stats.ttest_ind(La[:, genes_keep], Lb[:, genes_keep],
                               axis=0, equal_var=False)
        method = f'pseudobulk_welch_t(nA={len(ka)},nB={len(kb)})'

    lfc = La[:, genes_keep].mean(axis=0) - Lb[:, genes_keep].mean(axis=0)

    # NOTE ON A REVERTED CHANGE (2026-09-26). A variance floor was added here to
    # stop suspected n=2 false positives, then removed. A single-seed simulation
    # had suggested plain Welch emits FDR-significant genes under a pure null at
    # 2 vs 4 animals. A properly controlled test across 40 independent seeds
    # showed the opposite: plain Welch is CORRECTLY CALIBRATED (mean 0.0 false
    # positives), and the variance floor INTRODUCED them (mean 7.5, max 32 at
    # 3 vs 4 animals), because flooring the variances makes them more equal and
    # so inflates the Welch-Satterthwaite degrees of freedom. Do not re-add a
    # variance floor without testing its null calibration over many seeds.
    # Proper moderation needs limma/eBayes or pydeseq2, not a percentile clamp.

    df = pd.DataFrame({
        'gene': adata.var_names.values[genes_keep],
        'log2FC': lfc, 'stat': t, 'pval': p,
    })
    df['padj'] = bh_fdr(df['pval'].values)
    return df, method


def wilcoxon_de(mask_a, mask_b, genes_keep):
    """Cell-level rank_genes_groups. Pseudoreplicated when a side has <2 animals.

    FIX (2026-09-26): this used to test ALL genes while pseudobulk tested only
    the detected ones, so the two layers were FDR-corrected over different
    universes. On the same contrast pseudobulk tested 4,471 genes and wilcoxon
    17,965, a 4x difference in multiple-testing burden, and BH scales its
    threshold by n. That penalised the cell-level layer worst on the small
    myeloid populations, which are exactly the ones where it is the only usable
    layer. Both engines now receive genes_keep.
    """
    m = mask_a | mask_b
    tmp = adata[m][:, genes_keep].copy()
    if DE_EXPR_LAYER is not None:
        tmp.X = tmp.layers[DE_EXPR_LAYER]
    lab = np.where(mask_a[m], 'A', 'B')
    tmp.obs['_grp'] = pd.Categorical(lab, categories=['A', 'B'])
    try:
        sc.tl.rank_genes_groups(tmp, '_grp', groups=['A'], reference='B',
                                method='wilcoxon', use_raw=False)
    except Exception as e:
        return None, f'wilcoxon failed: {e}'
    r = tmp.uns['rank_genes_groups']
    df = pd.DataFrame({
        'gene': [x[0] for x in r['names']],
        'log2FC': [x[0] for x in r['logfoldchanges']],
        'stat': [x[0] for x in r['scores']],
        'pval': [x[0] for x in r['pvals']],
        'padj': [x[0] for x in r['pvals_adj']],
    })
    return df, 'wilcoxon_cell_level'


# %% ============================================================
# CELL 5 — RUN
# ============================================================
hr('CELL 5 — run contrasts')

spec_by_id = {r['contrast_id']: r for r in manifest}
results = []
n_run = 0

n_skipped_credible = 0
for _, row in man_df.iterrows():
    cid = row['contrast_id']
    if not row['runnable']:
        continue
    if CREDIBLE_ONLY and not row['pseudobulk_possible']:
        n_skipped_credible += 1
        continue
    spec = spec_by_id[cid]
    ma, _ = resolve_side(row['population'], row['tissue_scope'], spec['side_a'])
    mb, _ = resolve_side(row['population'], row['tissue_scope'], spec['side_b'])

    # gene filter: detected in enough cells across the two sides
    sl = adata[ma | mb]
    X = sl.layers[COUNTS_LAYER] if COUNTS_LAYER in adata.layers else sl.X
    det = np.asarray((X > 0).sum(axis=0)).ravel()
    genes_keep = np.where(det >= MIN_GENE_CELLS)[0]
    if genes_keep.size < 50:
        continue

    stem = f"{cid}_{slug(row['population'])}_{row['tissue_scope']}_{slug(row['contrast'])}"

    for layer in ('pseudobulk', 'wilcoxon'):
        if layer == 'pseudobulk':
            if not row['pseudobulk_possible']:
                continue
            df, method = pseudobulk_de(ma, mb, genes_keep)
        else:
            df, method = wilcoxon_de(ma, mb, genes_keep)
        if df is None:
            print(f'  {cid} {layer}: {method}')
            continue

        df = df.dropna(subset=['pval']).copy()
        df['contrast_id'] = cid
        df['population'] = row['population']
        df['tissue_scope'] = row['tissue_scope']
        df['contrast'] = row['contrast']
        df['evidence_layer'] = layer
        df['method'] = method
        df['pseudoreplicated'] = (layer == 'wilcoxon') and bool(row['pseudoreplicated'])
        df = df.sort_values('pval')
        df.to_csv(os.path.join(TAB_DIR, f'{stem}_{layer}.csv'), index=False)

        up = df[df['log2FC'] > 0].nsmallest(N_TOP, 'pval')
        dn = df[df['log2FC'] < 0].nsmallest(N_TOP, 'pval')
        top = pd.concat([up.assign(direction='up'), dn.assign(direction='down')])
        top.to_csv(os.path.join(TOP_DIR, f'{stem}_{layer}_top{N_TOP}.csv'),
                   index=False)

        # ---- ORA-ready gene lists at the primary threshold ----
        # One file per direction plus the background universe, which is this
        # contrast's own tested gene set. Using all 17,965 genes as background
        # when only 4,471 were testable inflates every enrichment p value.
        for gate, col in [('rawP', 'pval'), ('FDR', 'padj')]:
            sel = df[(df[col] < PRIMARY_ALPHA) &
                     (df['log2FC'].abs() >= PRIMARY_LFC)]
            for direction, sub in [('up', sel[sel['log2FC'] > 0]),
                                   ('down', sel[sel['log2FC'] < 0])]:
                if len(sub) >= MIN_GENES_FOR_ORA:
                    sub[['gene', 'log2FC', 'pval', 'padj']].to_csv(
                        os.path.join(ORA_DIR,
                                     f'{stem}_{layer}_{gate}_{direction}.csv'),
                        index=False)
        pd.DataFrame({'gene': df['gene']}).to_csv(
            os.path.join(ORA_DIR, f'{stem}_{layer}_background.csv'), index=False)

        rec = {'contrast_id': cid, 'population': row['population'],
               'tissue_scope': row['tissue_scope'], 'family': row['family'],
               'contrast': row['contrast'], 'evidence_layer': layer,
               'method': method,
               'pseudoreplicated': df['pseudoreplicated'].iloc[0] if len(df) else False,
               'n_cells_a': row['n_cells_a'], 'n_cells_b': row['n_cells_b'],
               'n_animals_a': row['n_animals_a'], 'n_animals_b': row['n_animals_b'],
               'min_animals_per_side': row['min_animals_per_side'],
               'pseudobulk_resid_df': row['pseudobulk_resid_df'],
               'near_duplicate_of': row['near_duplicate_of'],
               'n_genes_tested': len(df)}
        for a in ALPHA_LEVELS:
            for l in LFC_LEVELS:
                tag = f'a{a}_lfc{l}'
                rec[f'raw_{tag}'] = int(((df['pval'] < a) &
                                         (df['log2FC'].abs() >= l)).sum())
                rec[f'fdr_{tag}'] = int(((df['padj'] < a) &
                                         (df['log2FC'].abs() >= l)).sum())
        # directional counts at the primary threshold, for over-representation
        pt = f'a{PRIMARY_ALPHA}_lfc{PRIMARY_LFC}'
        for gate, col in [('raw', 'pval'), ('fdr', 'padj')]:
            sel = (df[col] < PRIMARY_ALPHA) & (df['log2FC'].abs() >= PRIMARY_LFC)
            rec[f'{gate}_{pt}_up'] = int((sel & (df['log2FC'] > 0)).sum())
            rec[f'{gate}_{pt}_down'] = int((sel & (df['log2FC'] < 0)).sum())
            rec[f'{gate}_{pt}_ORA_ready'] = bool(
                (sel & (df['log2FC'] > 0)).sum() >= MIN_GENES_FOR_ORA and
                (sel & (df['log2FC'] < 0)).sum() >= MIN_GENES_FOR_ORA)
        # the ORA background universe is THIS contrast's tested gene set, not all
        # genes in the object. Written out so the enrichment uses the right one.
        rec['ora_background_n'] = len(df)
        results.append(rec)

    n_run += 1
    if n_run % 25 == 0:
        print(f'  ... {n_run} contrasts done')

res_df = pd.DataFrame(results)
res_df.to_csv(os.path.join(OUT_DIR, 'deg_counts_by_contrast.csv'), index=False)
print(f'\n  ran {n_run} contrasts, {len(res_df)} contrast x layer results')
if CREDIBLE_ONLY:
    print(f'  CREDIBLE_ONLY skipped {n_skipped_credible} pseudoreplicated-only '
          'contrast(s). Set CREDIBLE_ONLY = False for the complete run.')


# %% ============================================================
# CELL 6 — THRESHOLD SUMMARY (split by evidence layer)
# ============================================================
hr('CELL 6 — threshold summary for GSEA cutoff selection')

if len(res_df):
    # ---------------------------------------------------------------
    # THE PRIMARY TABLE (Zoey's requested output)
    # One row per contrast x evidence layer. DEG count at |log2FC| > PRIMARY_LFC,
    # split into a raw-p column and an FDR column, so the GSEA filter is chosen
    # from the data.
    # ---------------------------------------------------------------
    ptag = f'a{PRIMARY_ALPHA}_lfc{PRIMARY_LFC}'
    raw_c, fdr_c = f'raw_{ptag}', f'fdr_{ptag}'
    if raw_c not in res_df.columns:
        raise RuntimeError(f'primary threshold columns missing ({raw_c}). '
                           f'Check PRIMARY_LFC={PRIMARY_LFC} is in LFC_LEVELS.')

    primary = res_df[['population', 'tissue_scope', 'family', 'contrast',
                      'evidence_layer', 'method', 'pseudoreplicated',
                      'n_cells_a', 'n_cells_b', 'n_animals_a', 'n_animals_b',
                      'min_animals_per_side', 'pseudobulk_resid_df',
                      'near_duplicate_of',
                      'n_genes_tested', raw_c, fdr_c]].copy()
    primary = primary.rename(columns={
        raw_c: f'n_DEG_rawP{PRIMARY_ALPHA}_absLFC{PRIMARY_LFC}',
        fdr_c: f'n_DEG_FDR{PRIMARY_ALPHA}_absLFC{PRIMARY_LFC}'})
    primary = primary.sort_values(['evidence_layer', 'population', 'tissue_scope',
                                  'family', 'contrast'])
    primary.to_csv(os.path.join(OUT_DIR, 'DEG_counts_primary_table.csv'), index=False)
    print(f'\n  PRIMARY TABLE  |log2FC| > {PRIMARY_LFC}, alpha {PRIMARY_ALPHA}')
    print(f'  wrote DEG_counts_primary_table.csv  ({len(primary)} rows)')
    print('\n  pseudobulk layer, non-pseudoreplicated (the credible rows):')
    pshow = primary[(primary['evidence_layer'] == 'pseudobulk') &
                    (~primary['pseudoreplicated'])]
    if len(pshow):
        print(pshow.drop(columns=['evidence_layer', 'pseudoreplicated',
                                  'method']).to_string(index=False))
    else:
        print('    none')

    # collapsed view: one row per population, both columns, credible layer only
    if len(pshow):
        rc = f'n_DEG_rawP{PRIMARY_ALPHA}_absLFC{PRIMARY_LFC}'
        fc = f'n_DEG_FDR{PRIMARY_ALPHA}_absLFC{PRIMARY_LFC}'
        bypop = (pshow.groupby('population')
                 .agg(contrasts=('contrast', 'size'),
                      median_raw=(rc, 'median'), max_raw=(rc, 'max'),
                      median_fdr=(fc, 'median'), max_fdr=(fc, 'max'))
                 .round(1).reset_index())
        print(f'\n  per population, |log2FC| > {PRIMARY_LFC}, pseudobulk only:')
        print(bypop.to_string(index=False))
        bypop.to_csv(os.path.join(OUT_DIR, 'DEG_counts_by_population.csv'),
                     index=False)

    # ---------------------------------------------------------------
    # the other LFC levels, for deciding whether 1.5 is the limiting factor
    # ---------------------------------------------------------------
    cols = [c for c in res_df.columns if c.startswith(('raw_', 'fdr_'))]
    thr = (res_df.groupby(['evidence_layer', 'pseudoreplicated'])[cols]
           .agg(['median', 'mean', 'max']).round(1))
    print('\n  ALL thresholds, by evidence layer (context for the choice above):')
    print(thr.to_string())
    thr.to_csv(os.path.join(OUT_DIR, 'threshold_summary_by_layer.csv'))

    # LFC sweep on the credible layer: is 1.5 too strict for this n?
    # duplicates excluded so a tissue-restricted population is not counted twice
    credr = res_df[(res_df['evidence_layer'] == 'pseudobulk') &
                   (~res_df['pseudoreplicated']) &
                   (res_df['near_duplicate_of'] == '')]
    n_dup_dropped = int(((res_df['evidence_layer'] == 'pseudobulk') &
                         (res_df['near_duplicate_of'] != '')).sum())
    if n_dup_dropped:
        print(f'\n  ({n_dup_dropped} near-duplicate pseudobulk row(s) excluded from '
              'the medians below)')
    if len(credr):
        sweep = []
        for l in LFC_LEVELS:
            t = f'a{PRIMARY_ALPHA}_lfc{l}'
            if f'raw_{t}' in credr.columns:
                sweep.append({'abs_log2FC': l,
                              'median_raw_p<0.05': credr[f'raw_{t}'].median(),
                              'median_FDR<0.05': credr[f'fdr_{t}'].median(),
                              'contrasts_with_0_at_FDR': int((credr[f'fdr_{t}'] == 0).sum()),
                              'n_contrasts': len(credr)})
        sw = pd.DataFrame(sweep)
        print('\n  LFC SWEEP, pseudobulk layer only:')
        print(sw.to_string(index=False))
        sw.to_csv(os.path.join(OUT_DIR, 'DEG_lfc_sweep.csv'), index=False)
            # DEG yield against residual df: the honest read on whether a zero is
        # biology or power
        if 'pseudobulk_resid_df' in credr.columns:
            pw = (credr.groupby('pseudobulk_resid_df')
                  .agg(contrasts=('contrast', 'size'),
                       median_raw=(f'raw_a{PRIMARY_ALPHA}_lfc{PRIMARY_LFC}', 'median'),
                       median_fdr=(f'fdr_a{PRIMARY_ALPHA}_lfc{PRIMARY_LFC}', 'median'))
                  .reset_index())
            print('\n  DEG yield by residual df (power), at the primary threshold:')
            print(pw.to_string(index=False))
            pw.to_csv(os.path.join(OUT_DIR, 'DEG_yield_by_power.csv'), index=False)

        prim = sw[sw['abs_log2FC'] == PRIMARY_LFC]
        if len(prim) and float(prim['median_FDR<0.05'].iloc[0]) == 0:
            print(f'\n  *** At |log2FC| > {PRIMARY_LFC} the median FDR-significant')
            print('      count is ZERO on the credible layer. Read the sweep above')
            print('      before committing to this filter for GSEA: if lower LFC')
            print('      levels are also zero at FDR, the limit is animal n, not')
            print('      the fold-change cut, and the honest GSEA input is a')
            print('      ranked list on raw p with no hard threshold.')

    byfam = (res_df.groupby(['family', 'evidence_layer'])[cols]
             .median().round(1).reset_index())
    print('\n  median DEG count by family x layer:')
    print(byfam.to_string(index=False))
    byfam.to_csv(os.path.join(OUT_DIR, 'threshold_summary_by_family.csv'),
                 index=False)

    credible = res_df[(res_df['evidence_layer'] == 'pseudobulk') &
                      (~res_df['pseudoreplicated'])]
    print('\n  RECOMMENDATION BASIS (pseudobulk, non-pseudoreplicated only):')
    if len(credible):
        for a in ALPHA_LEVELS:
            for l in LFC_LEVELS:
                tag = f'a{a}_lfc{l}'
                print(f'    raw p<{a}, |lfc|>={l}: median {credible[f"raw_{tag}"].median():.0f} '
                      f'DEGs | FDR<{a}: median {credible[f"fdr_{tag}"].median():.0f}')
        # ---- ORA readiness: the number that actually decides the threshold ----
        pt = f'a{PRIMARY_ALPHA}_lfc{PRIMARY_LFC}'
        print(f'\n  ORA READINESS at |log2FC|>{PRIMARY_LFC}, alpha {PRIMARY_ALPHA}')
        print(f'  (KEGG/GO over-representation needs >={MIN_GENES_FOR_ORA} genes in BOTH')
        print('   directions, since up and down are enriched separately)')
        for layer in ('pseudobulk', 'wilcoxon'):
            sub = res_df[(res_df['evidence_layer'] == layer) &
                         (~res_df['pseudoreplicated']) &
                         (res_df['near_duplicate_of'] == '')]
            if not len(sub):
                continue
            for gate in ('fdr', 'raw'):
                c = f'{gate}_{pt}_ORA_ready'
                if c in sub.columns:
                    n = int(sub[c].sum())
                    print(f'    {layer:11s} {gate.upper():4s}: {n:3d} of {len(sub)} '
                          f'contrasts enrichable  ({100 * n / len(sub):4.1f}%)')
        print('\n  If FDR leaves most contrasts unenrichable, the choice is between')
        print('  reporting those as "no significant genes" (a result, and honest) or')
        print('  moving to raw p with the fold-change floor and saying so in methods.')
        print('  Relaxing FDR to 0.10 or 0.25 is the middle option; both are in the')
        print('  threshold tables above.')
    else:
        print('    none available; every contrast is pseudoreplicated at this n.')
        print('    That is itself the finding: report the threshold off the')
        print('    wilcoxon layer only with the pseudoreplication stated.')

    try:
        fig, ax = plt.subplots(figsize=(18, 10))
        for layer, color in [('pseudobulk', HEX['pb']), ('wilcoxon', HEX['wx'])]:
            d = res_df[res_df['evidence_layer'] == layer]
            if not len(d):
                continue
            ax.hist(np.log10(d['fdr_a0.05_lfc0.0'] + 1), bins=40, alpha=0.65,
                    label=f'{layer} (n={len(d)})', color=color, edgecolor='white')
        ax.set_xlabel('log10(DEGs at FDR<0.05 + 1)', fontsize=30)
        ax.set_ylabel('contrasts', fontsize=30)
        ax.set_title('DEG yield by evidence layer', fontsize=34)
        ax.legend(fontsize=22, frameon=False)
        for sp in ('top', 'right'):
            ax.spines[sp].set_visible(False)
        for ext in ('pdf', 'png'):
            fig.savefig(os.path.join(OUT_DIR, f'deg_yield_by_layer.{ext}'),
                        dpi=FIG_DPI, bbox_inches='tight')
        plt.close(fig)
        print('\n  saved deg_yield_by_layer.pdf / .png')
    except Exception as e:
        print(f'  WARN: figure failed ({e})')

print(f'\n  outputs: {OUT_DIR}')
print('=' * 92)
print('DEG COMPARISON SUITE COMPLETE — input object untouched.')
print('=' * 92)
