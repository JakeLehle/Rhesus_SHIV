#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
04b — Why QC removes 73% of the control animals (READ-ONLY)
============================================================

04 reported this split:

    timepoint        pre-QC    post-QC    loss
    21 DPI           37,325     32,257    -13.6%
    Necropsy         37,429     33,375    -10.8%
    Pre              16,410     12,994    -20.8%
    Non-Infected     28,918      7,753    -73.2%

The two Indian non-infected controls, 34315 and 41903, lose roughly five times
the fraction the infected timepoints lose. The mito filter flags 29,216 cells
overall, 22.9% of everything, and appears to be landing on them specifically.

This matters beyond tidiness. Every surviving control cell is a heavily
selected subset, so an infected-versus-control contrast compares a lightly
filtered population against a severely filtered one. That is a referee
question, and it is cheap to characterise now and expensive to discover after
annotation is built on top of it.

It also connects to an open question already in the rebuild plan: why myeloid
clusters 0 and 4 are mito-high in animal 34315 only, currently handled by
flagging rather than by understanding. This script tests whether that is
confined to two clusters or is the whole control arm.

THE DECISIVE TEST (CELL 6)
--------------------------
Animal and library are almost confounded here, because the Indian controls
appear only in C and F. But not entirely: C also carries Chinese 41862, and F
also carries Chinese 41861 and 41862. So within one GEM well, with identical
handling downstream of the animal, Indian and Chinese cells can be compared
directly. If the Indian animals lose more INSIDE the same library, the cause
is upstream of library prep, in the sample. If the whole library loses
equally, it is the library.

WHICH RULE FIRES (CELL 5)
-------------------------
The mito filter is two rules OR'd together, and they have different meanings:

    pct_counts_mt > MT_THRESHOLD        a fixed 20% cutoff
    is_outlier(pct_counts_mt, 3 MAD)    relative to the GLOBAL median

The second is computed over the whole object, so a population with genuinely
higher mitochondrial content is penalised for differing from the majority
rather than for being poor quality. Which of the two is doing the work is the
difference between "these cells are degraded" and "these cells are unlike the
other animals". Reported separately.

QC CONSTANTS ARE COPIED FROM 04 AND MUST STAY IN STEP. Cell 3 recomputes the
flags from adata_prefilter.h5ad and asserts the totals reproduce what 04
printed; if they do not, the two scripts have drifted and everything below is
measuring something else.

Inputs (read-only)
  pipeline_results/04_ingest_qc/adata_prefilter.h5ad
  pipeline_results/04_ingest_qc/adata_qc.h5ad          (cross-check only)

Outputs
  qc_loss_by_animal.csv, qc_loss_by_library.csv, qc_loss_within_library.csv
  mito_distribution_by_animal.csv, qc_rule_attribution.csv
  F04b_qc_loss.pdf / .png

Environment:  conda activate sc_pre
Author: Jake Lehle
Date:   November 2026
"""

# %% ============================================================
# CELL 1 — CONFIGURATION (must match 04)
# ============================================================
import os

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import median_abs_deviation

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import warnings
warnings.filterwarnings('ignore')

try:
    from shiv_config import WORKDIR                                 # noqa: F401
except Exception:                                                   # noqa: BLE001
    WORKDIR = '/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV'

IN_DIR = os.path.join(WORKDIR, 'pipeline_results', '04_ingest_qc')
PREFILTER = os.path.join(IN_DIR, 'adata_prefilter.h5ad')
QC_OBJECT = os.path.join(IN_DIR, 'adata_qc.h5ad')

# --- copied verbatim from 04. If 04 changes, change these together. --------
MT_THRESHOLD = 20
MAD_COUNTS_GENES = 5
MAD_MITO = 3
MIN_GENES = 200

# 04 printed these; Cell 3 asserts it reproduces them
EXPECT_OUTLIER_TOTAL = 34537 + 10        # drops + whitelisted rescues
EXPECT_PREFILTER_CELLS = 127417

INDIAN = ['34315', '41903']
CHINESE = ['39272', '40702', '40707', '41861', '41862']

OUT_DIR = os.path.join(WORKDIR, 'pipeline_results', '04b_qc_loss')
os.makedirs(OUT_DIR, exist_ok=True)
pd.set_option('display.width', 240)
pd.set_option('display.max_columns', 60)

FIG_DPI = 300
FS_TITLE, FS_LABEL, FS_TICK, FS_LEG = 34, 30, 28, 28
HEX_INDIAN, HEX_CHINESE = '#C62828', '#2E7D32'
HEX_THRESH, HEX_MAD, HEX_GRID = '#1565C0', '#F9A825', '#E0E0E0'


def hr(t):
    print('\n' + '=' * 94, flush=True)
    print(t, flush=True)
    print('=' * 94, flush=True)


def save_csv(df, name, index=False):
    df.to_csv(os.path.join(OUT_DIR, name), index=index, lineterminator='\n')
    print(f'  wrote {name}  ({df.shape[0]:,} rows)', flush=True)


# %% ============================================================
# CELL 2 — LOAD
# ============================================================
hr('CELL 2 — load the pre-filter object')

if not os.path.exists(PREFILTER):
    raise SystemExit(f'not found: {PREFILTER}\nRun 04_ingest_qc.py first.')

a = sc.read_h5ad(PREFILTER)
obs = a.obs.copy()
del a                                 # only obs is needed; free the matrix
print(f'  pre-filter cells : {len(obs):,}')
print(f'  columns          : {len(obs.columns)}')

need = ['library', 'animal_id', 'timepoint', 'subspecies', 'tissue',
        'log1p_total_counts', 'log1p_n_genes_by_counts',
        'pct_counts_in_top_20_genes', 'pct_counts_mt', 'n_genes_by_counts',
        'total_counts', 'shiv_umi_cellranger']
missing = [c for c in need if c not in obs.columns]
if missing:
    raise SystemExit(f'pre-filter object is missing columns: {missing}')

obs['animal_id'] = obs['animal_id'].astype(str)
obs['library'] = obs['library'].astype(str)
obs['grp'] = np.where(obs.animal_id.isin(INDIAN), 'Indian control',
                      np.where(obs.animal_id.isin(CHINESE), 'Chinese infected',
                               'Unassigned'))

if os.path.exists(QC_OBJECT):
    q = sc.read_h5ad(QC_OBJECT)
    n_surv = q.n_obs
    surviving = set(q.obs_names)
    del q
    print(f'  survived QC      : {n_surv:,}')
else:
    surviving, n_surv = None, None
    print('  WARN  adata_qc.h5ad not found; the survival cross-check is skipped')


# %% ============================================================
# CELL 3 — RECOMPUTE THE QC FLAGS AND VERIFY THEY MATCH 04
# ============================================================
hr('CELL 3 — recompute QC flags, assert they reproduce 04')


def is_outlier(series, nmads):
    """Identical to 04: two-sided, global median and MAD."""
    med = series.median()
    mad = median_abs_deviation(series, nan_policy='omit')
    if mad == 0:
        mad = series.std()
        if mad == 0:
            return pd.Series(False, index=series.index), med, 0.0
    return (np.abs(series - med) > nmads * mad), med, mad


o1, m1, d1 = is_outlier(obs['log1p_total_counts'], MAD_COUNTS_GENES)
o2, m2, d2 = is_outlier(obs['log1p_n_genes_by_counts'], MAD_COUNTS_GENES)
o3, m3, d3 = is_outlier(obs['pct_counts_in_top_20_genes'], MAD_COUNTS_GENES)
omt, mmt, dmt = is_outlier(obs['pct_counts_mt'], MAD_MITO)

obs['f_counts'] = o1
obs['f_genes'] = o2
obs['f_top20'] = o3
obs['f_mt_mad'] = omt
obs['f_mt_hard'] = obs['pct_counts_mt'] > MT_THRESHOLD
obs['outlier'] = o1 | o2 | o3
obs['mt_outlier'] = omt | obs['f_mt_hard']
obs['would_drop'] = obs['outlier'] | obs['mt_outlier']

print(f'  log1p_total_counts       {int(o1.sum()):>7,}  '
      f'[median {m1:.2f}, MAD {d1:.2f}]')
print(f'  log1p_n_genes_by_counts  {int(o2.sum()):>7,}  '
      f'[median {m2:.2f}, MAD {d2:.2f}]')
print(f'  pct_counts_in_top_20     {int(o3.sum()):>7,}  '
      f'[median {m3:.2f}, MAD {d3:.2f}]')
print(f'  pct_counts_mt (MAD)      {int(omt.sum()):>7,}  '
      f'[median {mmt:.2f}, MAD {dmt:.2f}]')
print(f'  pct_counts_mt (>{MT_THRESHOLD}%)       '
      f'{int(obs["f_mt_hard"].sum()):>7,}')
print(f'  mt_outlier (either)      {int(obs["mt_outlier"].sum()):>7,}')
print(f'  would_drop (any rule)    {int(obs["would_drop"].sum()):>7,}')

ok = True
if len(obs) != EXPECT_PREFILTER_CELLS:
    print(f'\n  WARN  pre-filter cells {len(obs):,} != expected '
          f'{EXPECT_PREFILTER_CELLS:,}')
    ok = False
if abs(int(obs['would_drop'].sum()) - EXPECT_OUTLIER_TOTAL) > 5:
    print(f'\n  WARN  would_drop {int(obs["would_drop"].sum()):,} != expected '
          f'~{EXPECT_OUTLIER_TOTAL:,}')
    ok = False
if ok:
    print('\n  reproduces 04. The attribution below is measuring the same '
          'filter.')
else:
    print('\n  *** 04 and 04b have drifted. Reconcile the constants in Cell 1 '
          'before\n      trusting anything below.')


# %% ============================================================
# CELL 4 — LOSS BY ANIMAL AND BY LIBRARY
# ============================================================
hr('CELL 4 — who is being removed')


def loss_table(by):
    g = obs.groupby(by, observed=True)
    t = pd.DataFrame({
        'n_cells': g.size(),
        'drop_any': g['would_drop'].sum(),
        'drop_mt': g['mt_outlier'].sum(),
        'drop_mt_hard': g['f_mt_hard'].sum(),
        'drop_mt_mad_only': g.apply(
            lambda d: int((d['f_mt_mad'] & ~d['f_mt_hard']).sum())),
        'drop_counts_genes': g.apply(
            lambda d: int((d['outlier'] & ~d['mt_outlier']).sum())),
        'median_pct_mt': g['pct_counts_mt'].median(),
        'median_genes': g['n_genes_by_counts'].median(),
        'median_counts': g['total_counts'].median(),
    })
    t['pct_lost'] = (100 * t['drop_any'] / t['n_cells']).round(1)
    t['pct_lost_to_mt'] = (100 * t['drop_mt'] / t['n_cells']).round(1)
    return t.reset_index()


by_animal = loss_table('animal_id')
by_animal['group'] = np.where(by_animal.animal_id.isin(INDIAN), 'Indian ctrl',
                              np.where(by_animal.animal_id.isin(CHINESE),
                                       'Chinese inf', 'Unassigned'))
by_animal = by_animal.sort_values('pct_lost', ascending=False)
print('  BY ANIMAL')
print(by_animal[['animal_id', 'group', 'n_cells', 'pct_lost',
                 'pct_lost_to_mt', 'median_pct_mt', 'median_genes']]
      .to_string(index=False))
save_csv(by_animal, 'qc_loss_by_animal.csv')

by_lib = loss_table('library').sort_values('library')
print('\n  BY LIBRARY')
print(by_lib[['library', 'n_cells', 'pct_lost', 'pct_lost_to_mt',
              'median_pct_mt', 'median_genes']].to_string(index=False))
save_csv(by_lib, 'qc_loss_by_library.csv')


# %% ============================================================
# CELL 5 — WHICH RULE IS DOING THE WORK
# ============================================================
hr('CELL 5 — fixed 20% cutoff, or relative-to-global MAD')

rows = []
for grp in ['Indian control', 'Chinese infected', 'Unassigned']:
    m = obs['grp'] == grp
    if not m.any():
        continue
    sub = obs[m]
    rows.append({
        'group': grp,
        'n_cells': int(m.sum()),
        'median_pct_mt': round(sub['pct_counts_mt'].median(), 2),
        'mean_pct_mt': round(sub['pct_counts_mt'].mean(), 2),
        'p90_pct_mt': round(sub['pct_counts_mt'].quantile(0.90), 2),
        f'n_over_{MT_THRESHOLD}pct': int(sub['f_mt_hard'].sum()),
        f'pct_over_{MT_THRESHOLD}pct': round(100 * sub['f_mt_hard'].mean(), 1),
        'n_mad_only': int((sub['f_mt_mad'] & ~sub['f_mt_hard']).sum()),
        'pct_mad_only': round(
            100 * (sub['f_mt_mad'] & ~sub['f_mt_hard']).mean(), 1),
    })
attr = pd.DataFrame(rows)
print(attr.to_string(index=False))
save_csv(attr, 'qc_rule_attribution.csv')

print(f'\n  The global mito median is {mmt:.2f}% with MAD {dmt:.2f}, so the')
print(f'  3-MAD bound sits at {mmt + MAD_MITO * dmt:.2f}%, well BELOW the '
      f'{MT_THRESHOLD}% cutoff.')
print('  That means the MAD rule, not the fixed cutoff, is the binding')
print('  constraint, and it is defined relative to the majority population.')
print('  A group with genuinely higher mitochondrial content is then removed')
print('  for being unlike the others rather than for being degraded.')

mito_by_animal = (obs.groupby('animal_id', observed=True)['pct_counts_mt']
                  .describe(percentiles=[0.25, 0.5, 0.75, 0.9]).round(2)
                  .reset_index())
save_csv(mito_by_animal, 'mito_distribution_by_animal.csv')


# %% ============================================================
# CELL 6 — THE DECISIVE TEST: SAME LIBRARY, DIFFERENT ANIMAL
# ============================================================
hr('CELL 6 — animal effect or library effect')

print('  Libraries C and F each carry Indian controls AND a Chinese animal in')
print('  the same GEM well. Within one library, everything downstream of the')
print('  animal is identical, so a difference there is upstream of library')
print('  prep. A library that loses uniformly points at the library instead.\n')

within = []
for lib in sorted(obs['library'].unique()):
    sub = obs[obs.library == lib]
    if sub['grp'].nunique() < 2:
        continue
    for grp in ['Indian control', 'Chinese infected', 'Unassigned']:
        s = sub[sub.grp == grp]
        if not len(s):
            continue
        within.append({
            'library': lib, 'group': grp, 'n_cells': len(s),
            'pct_lost': round(100 * s['would_drop'].mean(), 1),
            'pct_lost_to_mt': round(100 * s['mt_outlier'].mean(), 1),
            'median_pct_mt': round(s['pct_counts_mt'].median(), 2),
            'median_genes': int(s['n_genes_by_counts'].median()),
        })
wdf = pd.DataFrame(within)
if len(wdf):
    print(wdf.to_string(index=False))
    save_csv(wdf, 'qc_loss_within_library.csv')

    print('\n  within-library contrast, Indian minus Chinese:')
    for lib in sorted(wdf.library.unique()):
        s = wdf[wdf.library == lib]
        i = s[s.group == 'Indian control']
        c = s[s.group == 'Chinese infected']
        if len(i) and len(c):
            d_lost = i.iloc[0]['pct_lost'] - c.iloc[0]['pct_lost']
            d_mt = i.iloc[0]['median_pct_mt'] - c.iloc[0]['median_pct_mt']
            print(f'    library {lib}: loss {i.iloc[0]["pct_lost"]:5.1f}% vs '
                  f'{c.iloc[0]["pct_lost"]:5.1f}%  (diff {d_lost:+5.1f} pts) | '
                  f'median mito {i.iloc[0]["median_pct_mt"]:5.2f}% vs '
                  f'{c.iloc[0]["median_pct_mt"]:5.2f}%  (diff {d_mt:+5.2f})')
    print('\n  A large positive difference means the ANIMAL, not the library.')
    print('  Near zero means the library, and the Indian animals only look')
    print('  worse because they sit in the two libraries that did badly.')
else:
    print('  No library carries more than one group; the test cannot run.')


# %% ============================================================
# CELL 7 — WHAT SURVIVES, AND IS IT A BIASED SUBSET
# ============================================================
hr('CELL 7 — the surviving control population')

obs['survived'] = (~obs['would_drop']) & (obs['n_genes_by_counts'] >= MIN_GENES)
rows = []
for grp in ['Indian control', 'Chinese infected']:
    for state in [True, False]:
        s = obs[(obs.grp == grp) & (obs.survived == state)]
        if not len(s):
            continue
        rows.append({
            'group': grp, 'survived': state, 'n_cells': len(s),
            'median_genes': int(s['n_genes_by_counts'].median()),
            'median_counts': int(s['total_counts'].median()),
            'median_pct_mt': round(s['pct_counts_mt'].median(), 2),
        })
surv_df = pd.DataFrame(rows)
print(surv_df.to_string(index=False))
save_csv(surv_df, 'surviving_population.csv')

print('\n  Compare median_genes between the surviving Indian and surviving')
print('  Chinese rows. If the survivors are similar, the filter removed a')
print('  distinct poor-quality population and the remainder is comparable.')
print('  If the Indian survivors still differ, the filter selected a tail and')
print('  any infected-versus-control contrast is partly a QC artefact.')

if surviving is not None:
    recomputed = set(obs.index[obs['survived']])
    overlap = len(recomputed & surviving)
    print(f'\n  cross-check against adata_qc.h5ad: {overlap:,} of '
          f'{len(surviving):,} survivors reproduced')
    print('  (the gap is doublet removal and min_cells, which run after this)')


# %% ============================================================
# CELL 8 — FIGURE
# ============================================================
hr('CELL 8 — figure')

order = ([x for x in by_animal.animal_id if x in INDIAN]
         + [x for x in by_animal.animal_id if x in CHINESE]
         + [x for x in by_animal.animal_id
            if x not in INDIAN and x not in CHINESE])
fig, axes = plt.subplots(1, 2, figsize=(30, 12),
                         gridspec_kw={'width_ratios': [1.4, 1]})

ax = axes[0]
data = [obs.loc[obs.animal_id == an, 'pct_counts_mt'].values for an in order]
cols = [HEX_INDIAN if an in INDIAN else HEX_CHINESE if an in CHINESE
        else '#9E9E9E' for an in order]
parts = ax.violinplot(data, showextrema=False, widths=0.85)
for b, c in zip(parts['bodies'], cols):
    b.set_facecolor(c)
    b.set_alpha(0.75)
for i, d in enumerate(data, start=1):
    ax.plot([i - 0.2, i + 0.2], [np.median(d)] * 2, color='#212121', lw=4)
ax.axhline(MT_THRESHOLD, color=HEX_THRESH, lw=4, ls='--',
           label=f'fixed cutoff {MT_THRESHOLD}%')
ax.axhline(mmt + MAD_MITO * dmt, color=HEX_MAD, lw=4, ls=':',
           label=f'{MAD_MITO}-MAD bound {mmt + MAD_MITO * dmt:.1f}%')
ax.set_xticks(range(1, len(order) + 1))
ax.set_xticklabels(order, fontsize=FS_TICK, rotation=45, ha='right')
ax.set_ylabel('percent mitochondrial counts', fontsize=FS_LABEL)
ax.set_title('Mitochondrial content per animal', fontsize=FS_TITLE, pad=16)
ax.tick_params(labelsize=FS_TICK)
ax.legend(fontsize=FS_LEG, frameon=False, loc='upper right')
ax.grid(axis='y', color=HEX_GRID, lw=1)
ax.set_axisbelow(True)
ax.set_ylim(0, min(60, float(obs['pct_counts_mt'].quantile(0.999))))

ax = axes[1]
bars = by_animal.set_index('animal_id').loc[order]
x = np.arange(len(order))
ax.bar(x, bars['pct_lost'], color=[c for c in cols], alpha=0.5,
       label='lost to any rule')
ax.bar(x, bars['pct_lost_to_mt'], color=[c for c in cols],
       label='lost to the mito rule')
ax.set_xticks(x)
ax.set_xticklabels(order, fontsize=FS_TICK, rotation=45, ha='right')
ax.set_ylabel('percent of cells removed', fontsize=FS_LABEL)
ax.set_title('QC loss per animal', fontsize=FS_TITLE, pad=16)
ax.tick_params(labelsize=FS_TICK)
ax.legend(fontsize=FS_LEG, frameon=False, loc='upper left')
ax.grid(axis='y', color=HEX_GRID, lw=1)
ax.set_axisbelow(True)

for ext in ('pdf', 'png'):
    fig.savefig(os.path.join(OUT_DIR, f'F04b_qc_loss.{ext}'), dpi=FIG_DPI,
                bbox_inches='tight')
plt.close(fig)
print('  wrote F04b_qc_loss.pdf / .png')

hr('DONE — 04b_qc_loss_diagnostic.py')
print(f'  outputs: {OUT_DIR}')
print('  nothing was modified.')
print('\n  Read Cell 6 first. It decides whether this is the animals or the')
print('  libraries, and that is what determines whether the control arm can')
print('  be used as it stands.')
