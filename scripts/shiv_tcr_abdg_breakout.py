#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Rhesus SHIV - TCR alpha-beta vs gamma-delta breakout (per animal x timepoint)
=============================================================================
Zoey asked to break the per-animal alpha-beta / gamma-delta split out BY
TIMEPOINT, to test whether the elite controller (40707) holds a steady
gamma-delta fraction while the other animals spike then fall as cells exhaust.
Second grouping: Chinese (infected) vs Indian (non-infected control) animals.

Definitions (consistent with the locked Tier 1/2 accounting):
  - gamma-delta  = is_gamma_delta  (strict paired TRG+TRD; 413 T cells total)
  - alpha-beta   = receptor_subtype == 'TRA+TRB'
  - receptor-defined = ab OR gd  (ambiguous / multichain / no-TCR excluded)
  - primary readout = pct_gd = gd / (ab + gd)   [of receptor-defined T]
    (using the receptor-defined denominator removes the timepoint-varying TCR
     recovery rate as a confound; pct of all T is reported too, for reference)

Run on biological T only (qc_exclude == False) so the QC doublet/mito bins do
not contribute spurious receptors. Read-only; writes tables + figures.

Caveats printed in the output:
  - Indian = 2 non-infected controls, small n; Chinese-vs-Indian is confounded
    with infection status. Report, do not over-read.
  - Small (ab+gd) denominators are flagged unreliable (< MIN_RECEPTOR).

Author: Jake Lehle
Date: July 2026
"""

# %% ============================================================
# CELL 1 - CONFIG
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
ADATA_IN    = os.path.join(MERGED, 'shiv_Tcell_annotated.h5ad')
OUT_DIR     = os.path.join(WORKING_DIR, 'annotation_output', 'tcr_abgd_breakout')

QC_COL       = 'qc_exclude'
ANIMAL_COL   = 'animal_id'
TIME_COL     = 'timepoint'
TISSUE_COL   = 'tissue'
GD_COL       = 'is_gamma_delta'
SUBTYPE_RCPT = 'receptor_subtype'
AB_VALUE     = 'TRA+TRB'

CHINESE = ['39272', '40702', '40707', '41861', '41862']
INDIAN  = ['34315', '41903']
ELITE   = '40707'
TIME_ORDER_CH = ['Pre', '21 DPI', 'Necropsy']     # longitudinal Chinese
MIN_RECEPTOR  = 20                                 # below this, % is unreliable

# colors (hex)
ELITE_COL = '#F9A825'
COMP_COL  = '#90A4AE'
CH_COL    = '#D32F2F'
IN_COL    = '#1976D2'
AB_COL    = '#90A4AE'
GD_COL_HX = '#D32F2F'

plt.rcParams.update({
    'font.size': 28, 'axes.titlesize': 34, 'axes.labelsize': 30,
    'xtick.labelsize': 24, 'ytick.labelsize': 24, 'legend.fontsize': 20,
    'figure.dpi': 100, 'pdf.fonttype': 42, 'ps.fonttype': 42,
})
os.makedirs(OUT_DIR, exist_ok=True)
pd.set_option('display.width', 200); pd.set_option('display.max_rows', 200)


def savefig(fig, name):
    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(OUT_DIR, f'{name}.{ext}'), bbox_inches='tight', dpi=300)
    print(f"    saved: {name}.pdf / .png")


def save_csv(df, name, index=True):
    df.to_csv(os.path.join(OUT_DIR, name), index=index)
    print(f"    saved: {name}")


# %% ============================================================
# CELL 2 - LOAD + DEFINE ab / gd
# ============================================================
print("=" * 80)
print("CELL 2 - LOAD + DEFINE ab / gd")
print("=" * 80)

adata = sc.read_h5ad(ADATA_IN)
bio = adata[~adata.obs[QC_COL].astype(bool)].copy()
print(f"  biological T: {bio.n_obs:,}")

o = bio.obs
o['_animal'] = o[ANIMAL_COL].astype(str)
o['_subspecies'] = np.where(o['_animal'].isin(CHINESE), 'Chinese',
                     np.where(o['_animal'].isin(INDIAN), 'Indian', 'Other'))
o['_gd'] = o[GD_COL].astype(bool)
o['_ab'] = (o[SUBTYPE_RCPT].astype(str) == AB_VALUE)
o['_recep'] = o['_ab'] | o['_gd']

print(f"  alpha-beta (TRA+TRB): {int(o['_ab'].sum()):,}")
print(f"  gamma-delta (strict): {int(o['_gd'].sum()):,}")
print(f"  receptor-defined T  : {int(o['_recep'].sum()):,} "
      f"({100*o['_recep'].mean():.1f}% of biological T)")


def breakout(df, group_cols):
    """Return counts + pct_gd per group. pct_gd of receptor-defined; also of all T."""
    g = df.groupby(group_cols, observed=True)
    tbl = pd.DataFrame({
        'n_T':      g.size(),
        'n_ab':     g['_ab'].sum().astype(int),
        'n_gd':     g['_gd'].sum().astype(int),
        'n_recep':  g['_recep'].sum().astype(int),
    })
    tbl['pct_gd_of_recep'] = (100 * tbl['n_gd'] / tbl['n_recep'].replace(0, np.nan)).round(2)
    tbl['pct_gd_of_T']     = (100 * tbl['n_gd'] / tbl['n_T'].replace(0, np.nan)).round(2)
    tbl['reliable']        = tbl['n_recep'] >= MIN_RECEPTOR
    return tbl


# %% ============================================================
# CELL 3 - PER ANIMAL x TIMEPOINT (pooled tissue) - the headline
# ============================================================
print("\n" + "=" * 80)
print("CELL 3 - PER ANIMAL x TIMEPOINT (pooled tissue)")
print("=" * 80)

at = breakout(o[o['_animal'].isin(CHINESE + INDIAN)], ['_animal', TIME_COL])
print(at.to_string())
save_csv(at, 'abgd_per_animal_timepoint.csv')
print(f"\n  NOTE: rows with n_recep < {MIN_RECEPTOR} are flagged reliable=False; "
      f"read their % with caution.")


# %% ============================================================
# CELL 4 - PER ANIMAL x TIMEPOINT x TISSUE (supplementary)
# ============================================================
print("\n" + "=" * 80)
print("CELL 4 - PER ANIMAL x TIMEPOINT x TISSUE")
print("=" * 80)
att = breakout(o[o['_animal'].isin(CHINESE + INDIAN)], ['_animal', TIME_COL, TISSUE_COL])
print(att.to_string())
save_csv(att, 'abgd_per_animal_timepoint_tissue.csv')


# %% ============================================================
# CELL 5 - CHINESE vs INDIAN (confounded with infection; small Indian n)
# ============================================================
print("\n" + "=" * 80)
print("CELL 5 - CHINESE vs INDIAN")
print("=" * 80)
sp_overall = breakout(o[o['_subspecies'].isin(['Chinese', 'Indian'])], ['_subspecies'])
print("  overall:")
print(sp_overall.to_string())
save_csv(sp_overall, 'abgd_chinese_vs_indian_overall.csv')

sp_tp = breakout(o[o['_subspecies'].isin(['Chinese', 'Indian'])], ['_subspecies', TIME_COL])
print("\n  by timepoint (Indian only exist at Non-Infected):")
print(sp_tp.to_string())
save_csv(sp_tp, 'abgd_chinese_vs_indian_timepoint.csv')
print("\n  CAVEAT: Indian = 2 non-infected controls (small n); Chinese vs Indian")
print("  is confounded with infection status. Descriptive, not a test.")


# %% ============================================================
# CELL 6 - FIGURES
# ============================================================
print("\n" + "=" * 80)
print("CELL 6 - FIGURES")
print("=" * 80)

# Fig 1: per-animal gd-fraction trajectory across timepoints (steady vs spike)
try:
    piv = at['pct_gd_of_recep'].unstack(TIME_COL)
    rel = at['reliable'].unstack(TIME_COL)
    piv = piv.reindex(index=[a for a in CHINESE if a in piv.index],
                      columns=[t for t in TIME_ORDER_CH if t in piv.columns])
    rel = rel.reindex_like(piv)
    x = np.arange(piv.shape[1])
    fig, ax = plt.subplots(figsize=(14, 10))
    for animal in piv.index:
        y = piv.loc[animal].values.astype(float)
        is_elite = animal == ELITE
        ax.plot(x, y, marker='o', markersize=16 if is_elite else 11,
                linewidth=4 if is_elite else 2.5,
                color=ELITE_COL if is_elite else COMP_COL,
                label=f"{animal}*" if is_elite else animal, zorder=3 if is_elite else 2)
        # hollow the unreliable points
        for xi, r in zip(x, rel.loc[animal].values):
            if not bool(r):
                ax.plot(xi, piv.loc[animal].values[xi], marker='o', markersize=18,
                        markerfacecolor='white',
                        markeredgecolor=ELITE_COL if is_elite else COMP_COL,
                        markeredgewidth=2, zorder=4)
    ax.set_xticks(x); ax.set_xticklabels(piv.columns)
    ax.set_ylabel('% gamma-delta (of ab+gd)')
    ax.set_title('gamma-delta fraction over infection\n(* = elite 40707; hollow = n<%d)' % MIN_RECEPTOR)
    ax.legend(title='animal'); ax.spines[['top', 'right']].set_visible(False)
    plt.tight_layout()
    savefig(fig, 'abgd_trajectory_per_animal')
    plt.close(fig)
except Exception as e:
    print(f"  WARN Fig 1 skipped: {e}")

# Fig 2: Chinese vs Indian gd fraction (overall)
try:
    fig, ax = plt.subplots(figsize=(9, 9))
    d = sp_overall.reindex(['Chinese', 'Indian'])
    ax.bar(range(len(d)), d['pct_gd_of_recep'].values, color=[CH_COL, IN_COL])
    ax.set_xticks(range(len(d))); ax.set_xticklabels(d.index)
    ax.set_ylabel('% gamma-delta (of ab+gd)')
    ax.set_title('gamma-delta fraction: Chinese vs Indian')
    ax.spines[['top', 'right']].set_visible(False)
    for i, (v, n) in enumerate(zip(d['pct_gd_of_recep'], d['n_recep'])):
        ax.text(i, v, f"n={int(n)}", ha='center', va='bottom', fontsize=18)
    plt.tight_layout()
    savefig(fig, 'abgd_chinese_vs_indian')
    plt.close(fig)
except Exception as e:
    print(f"  WARN Fig 2 skipped: {e}")

# Fig 3: absolute ab / gd stacked per animal x timepoint (denominator sizes)
try:
    sub = at.reset_index()
    sub = sub[sub['_animal'].isin(CHINESE)]
    sub['xt'] = sub['_animal'] + '\n' + sub[TIME_COL].astype(str)
    order = [f"{a}\n{t}" for a in CHINESE for t in TIME_ORDER_CH
             if f"{a}\n{t}" in set(sub['xt'])]
    sub = sub.set_index('xt').reindex(order)
    x = np.arange(len(sub))
    fig, ax = plt.subplots(figsize=(20, 10))
    ax.bar(x, sub['n_ab'], color=AB_COL, label='alpha-beta')
    ax.bar(x, sub['n_gd'], bottom=sub['n_ab'], color=GD_COL_HX, label='gamma-delta')
    ax.set_xticks(x); ax.set_xticklabels(sub.index, rotation=45, ha='right')
    ax.set_ylabel('receptor-defined T cells')
    ax.set_title('alpha-beta / gamma-delta counts per animal x timepoint')
    ax.legend(); ax.spines[['top', 'right']].set_visible(False)
    plt.tight_layout()
    savefig(fig, 'abgd_counts_per_animal_timepoint')
    plt.close(fig)
except Exception as e:
    print(f"  WARN Fig 3 skipped: {e}")


def _plot_trajectory(ax, piv, rel):
    """Per-animal gamma-delta trajectory on one axis (elite gold, hollow=unreliable)."""
    x = np.arange(piv.shape[1])
    for animal in piv.index:
        y = piv.loc[animal].values.astype(float)
        is_elite = animal == ELITE
        ax.plot(x, y, marker='o', markersize=16 if is_elite else 11,
                linewidth=4 if is_elite else 2.5,
                color=ELITE_COL if is_elite else COMP_COL,
                label=f"{animal}*" if is_elite else animal,
                zorder=3 if is_elite else 2)
        for xi, r in zip(x, rel.loc[animal].values):
            if not bool(r):
                ax.plot(xi, piv.loc[animal].values[xi], marker='o', markersize=18,
                        markerfacecolor='white',
                        markeredgecolor=ELITE_COL if is_elite else COMP_COL,
                        markeredgewidth=2, zorder=4)
    ax.set_xticks(x); ax.set_xticklabels(piv.columns)
    ax.spines[['top', 'right']].set_visible(False)


# Fig 4: tissue-split per-animal trajectory (PBMC | LN), shared y for comparison
try:
    fig, axes = plt.subplots(1, 2, figsize=(26, 11), sharey=True)
    for ax, tissue in zip(axes, ['PBMC', 'LN']):
        sub_t = att.xs(tissue, level=TISSUE_COL)          # index = (_animal, timepoint)
        piv = sub_t['pct_gd_of_recep'].unstack(TIME_COL)
        rel = sub_t['reliable'].unstack(TIME_COL)
        piv = piv.reindex(index=[a for a in CHINESE if a in piv.index],
                          columns=[t for t in TIME_ORDER_CH if t in piv.columns])
        rel = rel.reindex_like(piv).fillna(False)
        _plot_trajectory(ax, piv, rel)
        ax.set_title(tissue)
    axes[0].set_ylabel('% gamma-delta (of ab+gd)')
    axes[0].legend(title='animal', fontsize=16)
    fig.suptitle('gamma-delta fraction over infection, split by tissue\n'
                 f'(* = elite 40707; hollow = n<{MIN_RECEPTOR})', fontsize=32)
    plt.tight_layout()
    savefig(fig, 'abgd_trajectory_per_animal_by_tissue')
    plt.close(fig)
except Exception as e:
    print(f"  WARN Fig 4 skipped: {e}")


# %% ============================================================
# CELL 7 - SUMMARY
# ============================================================
print("\n" + "=" * 80)
print("CELL 7 - SUMMARY")
print("=" * 80)
print(f"""
  Outputs in {OUT_DIR}:
   - abgd_per_animal_timepoint.csv (+ trajectory figure)  <- the steady-vs-spike
     read: does 40707 (gold) hold a flatter gamma-delta fraction while others
     rise then fall? Read hollow points (n<{MIN_RECEPTOR}) with caution.
   - abgd_per_animal_timepoint_tissue.csv (+ by-tissue trajectory figure)  <- the
     honest view: PBMC | LN panels. Check whether the pooled spike-and-fall
     survives within a tissue or is a PBMC/LN sampling shift across timepoints.
   - abgd_chinese_vs_indian_{{overall,timepoint}}.csv (+ figure)  <- confounded
     with infection; Indian = 2 small non-infected controls. Descriptive only.

  Denominator note: pct is of receptor-defined T (ab+gd), which removes the
  timepoint-varying TCR recovery rate; pct_gd_of_T is reported alongside.
""")
print("=" * 80)
print("TCR ab/gd BREAKOUT COMPLETE.")
print("=" * 80)
