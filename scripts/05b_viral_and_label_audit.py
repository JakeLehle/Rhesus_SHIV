#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
05b — Is the viral signal real, and are the Tier 1 labels real?
================================================================
READ-ONLY. Runs off adata_annotated.h5ad in a couple of minutes. Writes
nothing into any path a pipeline stage owns.

WHY THIS EXISTS
---------------
05 reported 27 of 59 viral cells as B cells. HIV and SHIV bind B cells via
complement receptors and DC-SIGN, and B cells can carry virions on their
surface as a trans-infection reservoir, but productive infection of B cells is
not supported: they lack CD4. So either the viral call is wrong, the B label is
wrong, or both.

Three things in the 05 output already argue that the viral calls are largely
ambient rather than tropic, and this script tests each one properly:

  1. Two of the 59 viral cells are Erythrocytes and one is a Mast cell. An
     erythrocyte has no nucleus and no CD4. A viral UMI there cannot be
     infection, so it is proof that ambient viral RNA is present in this
     dataset at a rate that produces impossible calls.

  2. Under a null of "viral cells are drawn in proportion to lineage
     abundance", which is what ambient soup does, B cells come out at a
     nominal p = 0.013 (27 observed against 18.4 expected). Erythrocytes come
     out at p = 0.012. The impossible category is as significant as the
     interesting one, and neither survives correction for the ten lineages
     tested (Bonferroni 0.12). So the B enrichment is not established.

  3. Mean viral burden is 3.0 UMI per cell. A productively infected cell
     transcribes virus heavily and should carry a substantial fraction of its
     transcriptome as viral reads. 1 to 3 UMI out of several thousand is the
     ambient range. 05 printed the total but never the per-cell burden, which
     is the measurement that distinguishes the two.

THE DECISIVE TEST IS THE CONTROL ARM, AND 05 NEVER RAN IT
----------------------------------------------------------
26,026 cells are from Non-Infected Indian control animals. If viral UMI appear
there, the rate is a direct measurement of ambient plus index hopping, because
those animals have no virus to detect. CELL 3 is that test. Everything else in
this script is supporting.

IT ALSO AUDITS THE TIER 1 LABELS, BECAUSE TWO DEFECTS ARE IN THE 05 CODE
-------------------------------------------------------------------------
Found by reading 05, not inferred from the output:

  DEFECT A  The CD3 hard gate is applied to only two labels:

                lab2 = np.where((lab2 == 'NK cells')   & cd3_pos, 'T cells', lab2)
                lab2 = np.where((lab2 == 'Unassigned') & cd3_pos, 'T cells', lab2)

            A CD3-positive cell that argmaxed to B cells, Monocytes,
            Macrophages, DCs, Mast, Erythrocytes or Platelets keeps that label
            and is never examined. The gate exists because the scores are
            untrustworthy; restricting it to the two labels where the v10
            collision happened to be known has no justification. CELL 5
            measures what that asymmetry costs.

  DEFECT B  Nothing requires a cell to express any marker of the lineage it
            was assigned. Assignment is argmax over sc.tl.score_genes with a
            floor on the SCORE:

                lab = np.array(list(TIER1_MARKERS))[S.argmax(axis=1)]
                lab = np.where(S.max(axis=1) < SCORE_THRESHOLD, 'Unassigned', lab)

            score_genes subtracts a control-gene mean, so in a sparse cell both
            terms sit near zero and noise picks the winner. A cell with zero
            counts of CD79A, MS4A1, CD19, BANK1 and PAX5 can still be labelled
            B cells. CELL 4 measures how often that happens, per lineage.

            Note that 31% B cells is NOT by itself evidence of a bug: four of
            the eight libraries are lymph node, where B cells legitimately run
            40 to 60% of lymphocytes. CELL 4 settles it on marker detection
            instead of on plausibility, and CELL 6 splits it by tissue.

Inputs   pipeline_results/05_embed_annotate/adata_annotated.h5ad
         pipeline_results/04_ingest_qc/adata_prefilter.h5ad   (optional)
Outputs  pipeline_results/05b_viral_label_audit/
         viral_cells_annotated.csv      the 59 cells, one row each, read by eye
         viral_ambient_test.csv         observed vs expected per lineage
         viral_by_infection_status.csv  THE control test
         tier1_marker_support.csv       per label, who expresses what
         tier1_cd3_by_label.csv         what Defect A costs
         F05b_*.pdf / .png

Environment:  conda activate sc_pre
Submit with:  sbatch run_05b_viral_and_label_audit.sh   (or just run it)
Author: Jake Lehle
Date:   November 2026
"""

# %% ============================================================
# CELL 1 — CONFIGURATION
# ============================================================
import os
os.environ.setdefault('MPLBACKEND', 'Agg')

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse, stats

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import warnings
warnings.filterwarnings('ignore')

WORKDIR = '/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV'
ANNOT = os.path.join(WORKDIR, 'pipeline_results', '05_embed_annotate',
                     'adata_annotated.h5ad')
PREFILTER = os.path.join(WORKDIR, 'pipeline_results', '04_ingest_qc',
                         'adata_prefilter.h5ad')
OUT_DIR = os.path.join(WORKDIR, 'pipeline_results', '05b_viral_label_audit')

COUNTS_LAYER = 'counts'
SHIV_COL = 'shiv_umi_cellranger'
ANIMAL_COL = 'animal_id'
TIME_COL = 'timepoint'
TISSUE_COL = 'tissue'

# Indian animals are the uninfected controls. This is the control arm.
INDIAN_CONTROL = ['34315', '41903']
CONTROL_TIMEPOINT = 'Non-Infected'

# Tier 1 panels, copied verbatim from 05 Cell 1 so the audit tests the same
# gene sets that made the calls.
TIER1_MARKERS = {
    'T cells':      ['CD3D', 'CD3E', 'CD3G'],
    'B cells':      ['CD79A', 'MS4A1', 'CD19', 'BANK1', 'PAX5'],
    'NK cells':     ['KLRD1', 'KLRC3', 'GNLY', 'PRF1', 'GZMB', 'NKG7'],
    'Monocytes':    ['S100A8', 'S100A9', 'VCAN', 'CD14', 'LYZ'],
    'Macrophages':  ['CD68', 'CD163', 'MRC1', 'MERTK', 'APOC1'],
    'DCs':          ['FLT3', 'CLEC9A', 'CLEC10A', 'IL3RA', 'CLEC4C'],
    'Mast cells':   ['CPA3', 'KIT', 'HPGD', 'ENPP3'],
    'Erythrocytes': ['HBB', 'HBA1', 'HBA2'],
    'Platelets':    ['PPBP', 'PF4', 'GP9'],
}
CD3_GENES = ['CD3D', 'CD3E', 'CD3G']
POS_MIN = 1               # raw counts for "detected"

# SHIV AD8-EO is CCR5-tropic and CD4-dependent. A cell that never had CD4
# cannot be productively infected, so CD4 and CCR5 are an orthogonal check on
# any cell claimed as infected.
ENTRY_GENES = ['CD4', 'CCR5', 'CXCR4']

# Lineages that cannot be productively infected under any reading. A viral UMI
# here is a direct measurement of the ambient rate.
IMPOSSIBLE_LINEAGES = ['Erythrocytes', 'Platelets', 'Mast cells']

# A productively infected cell should carry a real viral fraction, not one
# stray molecule. These are reporting bands, not a classifier.
BURDEN_BANDS = [1, 2, 3, 5, 10, 25]
VIRAL_FRAC_BANDS = [1e-4, 1e-3, 1e-2, 1e-1]

FIG_DPI = 300
plt.rcParams.update({
    'font.size': 28, 'axes.titlesize': 34, 'axes.labelsize': 30,
    'xtick.labelsize': 28, 'ytick.labelsize': 28, 'legend.fontsize': 20,
    'figure.dpi': 100, 'pdf.fonttype': 42, 'ps.fonttype': 42,
})
HEX_CTRL, HEX_INF, HEX_GRID = '#C62828', '#2E7D32', '#E0E0E0'
HEX_BAR, HEX_IMPOSSIBLE = '#1976D2', '#F9A825'

os.makedirs(OUT_DIR, exist_ok=True)
pd.set_option('display.width', 250)
pd.set_option('display.max_columns', 60)
pd.set_option('display.max_rows', 200)


def hr(t):
    print('\n' + '=' * 94, flush=True)
    print(t, flush=True)
    print('=' * 94, flush=True)


def save_csv(df, name, index=True):
    df.to_csv(os.path.join(OUT_DIR, name), index=index, lineterminator='\n')
    print(f'  wrote {name}  ({df.shape[0]:,} rows)', flush=True)


def savefig(fig, stem):
    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(OUT_DIR, f'{stem}.{ext}'), dpi=FIG_DPI,
                    bbox_inches='tight')
    plt.close(fig)
    print(f'  wrote {stem}.pdf / .png', flush=True)


def present(genes, a):
    return [g for g in genes if g in a.var_names]


def counts_of(a, genes):
    """Per-cell raw counts for each gene, as a DataFrame. Empty if none."""
    g = present(genes, a)
    if not g:
        return pd.DataFrame(index=a.obs_names)
    m = a[:, g].layers[COUNTS_LAYER] if COUNTS_LAYER in a.layers else a[:, g].X
    m = m.toarray() if sparse.issparse(m) else np.asarray(m)
    return pd.DataFrame(m, index=a.obs_names, columns=g)


# %% ============================================================
# CELL 2 — LOAD
# ============================================================
hr('CELL 2 — load')

if not os.path.exists(ANNOT):
    raise SystemExit(f'not found: {ANNOT}\nRun 05_embed_annotate.py first.')

a = sc.read_h5ad(ANNOT)
print(f'  {a.n_obs:,} cells x {a.n_vars:,} genes')

if COUNTS_LAYER not in a.layers:
    print(f'  WARN  layers["{COUNTS_LAYER}"] absent; marker detection will')
    print('        fall back to normalized X, which still answers "detected'
          ' or not"')
if SHIV_COL not in a.obs.columns:
    raise SystemExit(f'  STOP: obs["{SHIV_COL}"] absent. Nothing to audit.')

a.obs[ANIMAL_COL] = a.obs[ANIMAL_COL].astype(str)
a.obs['tier1_celltype'] = a.obs['tier1_celltype'].astype(str)
a.obs['is_control'] = (a.obs[ANIMAL_COL].isin(INDIAN_CONTROL)
                       | (a.obs[TIME_COL].astype(str) == CONTROL_TIMEPOINT))
v = a.obs[SHIV_COL].values.astype(float)
a.obs['shiv_pos'] = v > 0

n_v = int(a.obs['shiv_pos'].sum())
print(f'  viral cells (>=1 UMI) : {n_v}')
print(f'  total viral UMI       : {v.sum():,.0f}')
print(f'  mean burden           : {v[v > 0].mean():.2f} UMI per viral cell')
print(f'  control-arm cells     : {int(a.obs["is_control"].sum()):,}')

if 'total_counts' in a.obs.columns:
    tot = a.obs['total_counts'].values.astype(float)
else:
    X = a.layers[COUNTS_LAYER] if COUNTS_LAYER in a.layers else a.X
    tot = np.asarray(X.sum(axis=1)).ravel()
a.obs['_total'] = tot
a.obs['viral_frac'] = np.where(tot > 0, v / tot, 0.0)


# %% ============================================================
# CELL 3 — THE CONTROL TEST. Everything else is supporting.
# ============================================================
hr('CELL 3 — viral UMI in the uninfected control arm')

print('  The Indian animals were never infected. Any viral UMI there is')
print('  ambient, index hopping, or a mapping artifact. The rate in the')
print('  control arm is the floor below which nothing in the infected arm')
print('  can be called real.')
print()

grp = np.where(a.obs['is_control'], 'control (Non-Infected)', 'infected')
tab = pd.DataFrame({
    'n_cells': pd.Series(grp).value_counts(),
    'viral_cells': pd.Series(grp)[a.obs['shiv_pos'].values].value_counts(),
    'viral_UMI': pd.Series(v).groupby(grp).sum(),
    'total_UMI': pd.Series(tot).groupby(grp).sum(),
}).fillna(0)
tab['pct_cells_viral'] = (100 * tab.viral_cells / tab.n_cells).round(4)
tab['viral_UMI_per_1M'] = (1e6 * tab.viral_UMI / tab.total_UMI).round(3)
print(tab.to_string())
save_csv(tab, 'viral_by_infection_status.csv')

ctrl_rate = float(tab.loc[tab.index.str.startswith('control'),
                          'pct_cells_viral'].iloc[0]) \
    if tab.index.str.startswith('control').any() else 0.0
inf_rate = float(tab.loc['infected', 'pct_cells_viral']) \
    if 'infected' in tab.index else 0.0

print()
if ctrl_rate == 0:
    print('  ZERO viral cells in the control arm.')
    print('  That is the good outcome: it rules out index hopping and a')
    print('  mapping artifact as the whole explanation, because both would')
    print('  deposit reads in the control animals too. It does NOT rule out')
    print('  within-library ambient, since soup is per-GEM-well and the')
    print('  control animals sit in libraries C and F with infected animals.')
    print('  CELL 4 handles that case.')
else:
    print(f'  {ctrl_rate:.4f}% of control cells carry viral UMI, against')
    print(f'  {inf_rate:.4f}% of infected cells.')
    if inf_rate > 0:
        print(f'  ratio infected:control = {inf_rate / ctrl_rate:.1f}x')
    print('  A ratio near 1 means the infected-arm signal IS the background.')
    print('  Subtract the control rate before claiming any infected count.')

# Per library, because ambient is per GEM well and libraries C/F are mixed.
print('\n  per library (C and F carry both subspecies in one well):')
lib = pd.DataFrame({
    'n_cells': a.obs.groupby('library', observed=True).size(),
    'n_control': a.obs.groupby('library', observed=True)['is_control'].sum(),
    'viral_cells': a.obs.groupby('library', observed=True)['shiv_pos'].sum(),
    'viral_UMI': a.obs.groupby('library', observed=True)[SHIV_COL].sum(),
})
ctl = a.obs[a.obs['is_control']]
lib['viral_cells_in_controls'] = (
    ctl.groupby('library', observed=True)['shiv_pos'].sum()
    .reindex(lib.index, fill_value=0))
print(lib.to_string())
save_csv(lib, 'viral_by_library.csv')
print('\n  viral_cells_in_controls is the within-well ambient measurement.')
print('  Non-zero in a mixed library means soup is reaching control cells in')
print('  the same well, which is the mechanism that would also put virus in')
print('  B cells and erythrocytes.')


# %% ============================================================
# CELL 4 — CELLS THAT CANNOT BE INFECTED
# ============================================================
hr('CELL 4 — the impossible lineages')

print('  Erythrocytes have no nucleus and no CD4. Platelets have no nucleus.')
print('  Mast cells are not a documented productive target. A viral UMI in')
print('  any of them is a direct measurement of the ambient rate, with no')
print('  interpretation required.')
print()

imp = a.obs['tier1_celltype'].isin(IMPOSSIBLE_LINEAGES)
rows = []
for lin in IMPOSSIBLE_LINEAGES + ['ALL IMPOSSIBLE', 'everything else']:
    if lin == 'ALL IMPOSSIBLE':
        m = imp.values
    elif lin == 'everything else':
        m = ~imp.values
    else:
        m = (a.obs['tier1_celltype'] == lin).values
    if m.sum() == 0:
        continue
    rows.append({
        'lineage': lin, 'n_cells': int(m.sum()),
        'viral_cells': int(a.obs['shiv_pos'].values[m].sum()),
        'pct_viral': round(100 * a.obs['shiv_pos'].values[m].mean(), 4),
        'viral_UMI': float(v[m].sum()),
        'viral_UMI_per_1M': round(1e6 * v[m].sum() / max(tot[m].sum(), 1), 3),
    })
impt = pd.DataFrame(rows)
print(impt.to_string(index=False))
save_csv(impt, 'viral_in_impossible_lineages.csv', index=False)

n_imp = int(a.obs['shiv_pos'].values[imp.values].sum())
print()
if n_imp:
    print(f'  {n_imp} viral cell(s) sit in lineages that cannot be infected.')
    print(f'  Of {n_v} viral calls, that is {100 * n_imp / n_v:.1f}% known-false.')
    print('  Ambient viral RNA is therefore present in this dataset at a')
    print('  measurable rate. The question for B cells is no longer whether')
    print('  ambient exists, it is whether B cells carry MORE than ambient.')
else:
    print('  No viral UMI in the impossible lineages. That is a point in')
    print('  favour of the calls being real, though with so few cells in')
    print('  those lineages the test has little power.')


# %% ============================================================
# CELL 5 — ARE THE VIRAL CELLS ENRICHED BEYOND LINEAGE ABUNDANCE?
# ============================================================
hr('CELL 5 — ambient null: viral cells drawn in proportion to abundance')

print('  Ambient soup contaminates cells in proportion to how many of them')
print('  there are. If the viral calls are tropic, some lineage must carry')
print('  MORE than its share. This tests each lineage against that null.')
print()

lin_n = a.obs['tier1_celltype'].value_counts()
vir_n = a.obs.loc[a.obs['shiv_pos'], 'tier1_celltype'].value_counts()
rows = []
for lin, n in lin_n.items():
    p = n / a.n_obs
    o = int(vir_n.get(lin, 0))
    e = p * n_v
    rows.append({
        'lineage': lin, 'n_cells': int(n), 'pct_of_object': round(100 * p, 2),
        'expected': round(e, 2), 'observed': o,
        'fold': round(o / e, 2) if e > 0 else np.nan,
        'p_greater': round(
            stats.binomtest(o, n_v, p, alternative='greater').pvalue, 4),
    })
amb = pd.DataFrame(rows).sort_values('n_cells', ascending=False)
amb['p_bonferroni'] = (amb.p_greater * len(amb)).clip(upper=1.0).round(3)
amb['cannot_be_infected'] = amb.lineage.isin(IMPOSSIBLE_LINEAGES)
print(amb.to_string(index=False))
save_csv(amb, 'viral_ambient_test.csv', index=False)

keep = amb[amb.expected >= 1]
chi = float(((keep.observed - keep.expected) ** 2 / keep.expected).sum())
df_ = max(len(keep) - 1, 1)
print(f'\n  chi-square over bins with expected >=1: {chi:.2f} on {df_} df, '
      f'p = {1 - stats.chi2.cdf(chi, df_):.4f}')

sig = amb[amb.p_bonferroni < 0.05]
print()
if len(sig) == 0:
    print('  NO lineage is enriched beyond its abundance after correcting for')
    print('  the number of lineages tested. With 59 viral cells there is not')
    print('  enough signal to establish tropism either way, and the honest')
    print('  statement is that the data do not support a lineage claim.')
else:
    print('  enriched after correction:')
    print(sig[['lineage', 'observed', 'expected', 'fold',
               'p_bonferroni', 'cannot_be_infected']].to_string(index=False))
    if sig.cannot_be_infected.any():
        print('\n  NOTE a lineage that CANNOT be infected is among the')
        print('  enriched ones. Whatever mechanism produces that signal')
        print('  produces the others, so none of them is evidence of tropism.')


# %% ============================================================
# CELL 6 — PER-CELL VIRAL BURDEN
# ============================================================
hr('CELL 6 — viral burden per cell')

print('  A productively infected cell transcribes virus heavily, so viral')
print('  reads should be a real fraction of its transcriptome. One or two')
print('  stray molecules out of several thousand is the ambient range. 05')
print('  printed the total and never this.')
print()

vp = v[a.obs['shiv_pos'].values]
fp = a.obs.loc[a.obs['shiv_pos'], 'viral_frac'].values
tp = tot[a.obs['shiv_pos'].values]

print(f'  viral UMI per viral cell: min {vp.min():.0f}  median '
      f'{np.median(vp):.0f}  mean {vp.mean():.2f}  max {vp.max():.0f}')
print(f'  viral fraction of transcriptome: median {np.median(fp):.2e}  '
      f'max {fp.max():.2e}')
print(f'  total UMI in those cells: median {np.median(tp):,.0f}')

bands = pd.DataFrame({
    'min_viral_UMI': BURDEN_BANDS,
    'n_cells': [int((vp >= b).sum()) for b in BURDEN_BANDS],
})
bands['pct_of_viral_cells'] = (100 * bands.n_cells / n_v).round(1)
print('\n  cells at or above each burden:')
print(bands.to_string(index=False))
save_csv(bands, 'viral_burden_bands.csv', index=False)

fbands = pd.DataFrame({
    'min_viral_fraction': VIRAL_FRAC_BANDS,
    'n_cells': [int((fp >= b).sum()) for b in VIRAL_FRAC_BANDS],
})
print('\n  cells at or above each viral fraction of their transcriptome:')
print(fbands.to_string(index=False))
save_csv(fbands, 'viral_fraction_bands.csv', index=False)

print('\n  For scale: in published scRNA-seq of productively infected cells,')
print('  viral reads commonly run 1% or more of the transcriptome. Count how')
print('  many cells clear 1e-2 above. If that number is zero, nothing in')
print('  this dataset looks productively infected, in any lineage.')

fig, axes = plt.subplots(1, 2, figsize=(26, 11))
axes[0].hist(vp, bins=np.arange(0.5, max(vp.max(), 2) + 1.5, 1),
             color=HEX_BAR, edgecolor='white')
axes[0].set_xlabel('viral UMI per cell')
axes[0].set_ylabel('cells')
axes[0].set_title('Viral burden', fontsize=34)
axes[1].scatter(tp, np.maximum(fp, 1e-6), s=180, c=HEX_BAR, alpha=0.75,
                edgecolor='white', linewidth=1.5)
axes[1].set_xscale('log')
axes[1].set_yscale('log')
axes[1].axhline(1e-2, color=HEX_CTRL, lw=3, ls='--')
axes[1].text(tp.min(), 1.3e-2, '1% of transcriptome', color=HEX_CTRL,
             fontsize=24)
axes[1].set_xlabel('total UMI in the cell')
axes[1].set_ylabel('viral fraction')
axes[1].set_title('Burden against depth', fontsize=34)
for ax in axes:
    ax.grid(color=HEX_GRID, lw=1)
    ax.set_axisbelow(True)
savefig(fig, 'F05b_viral_burden')


# %% ============================================================
# CELL 7 — THE 59 CELLS, ONE ROW EACH
# ============================================================
hr('CELL 7 — every viral cell, with the evidence against it')

mk = counts_of(a, sorted({g for gs in TIER1_MARKERS.values() for g in gs}
                         | set(ENTRY_GENES)))
cd3 = present(CD3_GENES, a)
b_mk = present(TIER1_MARKERS['B cells'], a)

cols = {}
for c in ('library', ANIMAL_COL, TIME_COL, TISSUE_COL, SHIV_COL, 'viral_frac',
          '_total', 'tier1_celltype', 'tier2_subtype', 'clusters',
          'sub_leiden', 'has_tcr', 'doublet_score', 'predicted_doublet',
          'qc_exclude', 'dx_cd3_pos'):
    if c in a.obs.columns:
        cols[c] = a.obs[c]
tbl = pd.DataFrame(cols)[a.obs['shiv_pos'].values].copy()
tbl = tbl.rename(columns={'_total': 'total_UMI', SHIV_COL: 'viral_UMI'})

if cd3:
    tbl['CD3_counts'] = mk.loc[tbl.index, cd3].sum(axis=1)
if b_mk:
    tbl['Bmarker_counts'] = mk.loc[tbl.index, b_mk].sum(axis=1)
    tbl['Bmarkers_detected'] = (mk.loc[tbl.index, b_mk] >= POS_MIN).sum(axis=1)
for g in present(ENTRY_GENES, a):
    tbl[g] = mk.loc[tbl.index, g]
tbl['is_control'] = a.obs.loc[tbl.index, 'is_control']
tbl.index.name = 'barcode'
tbl = tbl.sort_values('viral_UMI', ascending=False)

print(tbl.to_string())
save_csv(tbl.reset_index(), 'viral_cells_annotated.csv', index=False)

print('\n  HOW TO READ THIS')
print('  viral_UMI 1 or 2 with viral_frac near 1e-4 is ambient.')
print('  A B-labelled row with Bmarkers_detected = 0 was never shown to be a')
print('  B cell (Defect B). A B-labelled row with high CD3_counts is Defect A.')
print('  CD4 = 0 in a cell claimed as infected is an orthogonal problem:')
print('  SHIV AD8-EO is CD4-dependent.')

if b_mk:
    bv = tbl[tbl.tier1_celltype == 'B cells']
    if len(bv):
        print(f'\n  of the {len(bv)} B-labelled viral cells:')
        print(f'    express >=1 B marker : '
              f'{int((bv.Bmarkers_detected >= 1).sum())}')
        print(f'    express ZERO B marker: '
              f'{int((bv.Bmarkers_detected == 0).sum())}')
        if 'CD3_counts' in bv:
            print(f'    CD3-positive         : '
                  f'{int((bv.CD3_counts >= POS_MIN).sum())}')
        if 'CD4' in bv:
            print(f'    CD4-positive         : '
                  f'{int((bv.CD4 >= POS_MIN).sum())}')
        print(f'    median viral UMI     : {bv.viral_UMI.median():.0f}')
        print(f'    median viral fraction: {bv.viral_frac.median():.2e}')


# %% ============================================================
# CELL 8 — DEFECT B: DOES EACH LABEL EXPRESS ITS OWN MARKERS?
# ============================================================
hr('CELL 8 — marker support per Tier 1 label')

print('  For every label, the fraction of its cells that actually detect each')
print('  of its markers. A label whose cells mostly express none of its own')
print('  markers was assigned on score noise, not on evidence.')
print()

rows = []
for lin, genes in TIER1_MARKERS.items():
    m = (a.obs['tier1_celltype'] == lin).values
    if m.sum() == 0:
        continue
    g = present(genes, a)
    sub = mk.loc[a.obs_names[m], g] if g else pd.DataFrame(index=a.obs_names[m])
    det = (sub >= POS_MIN) if g else pd.DataFrame(index=sub.index)
    r = {'label': lin, 'n_cells': int(m.sum()),
         'n_markers': len(g),
         'pct_zero_markers': round(
             100 * (det.sum(axis=1) == 0).mean(), 1) if g else np.nan,
         'median_markers_detected': float(
             det.sum(axis=1).median()) if g else np.nan}
    for gene in genes:
        r[f'pct_{gene}'] = (round(100 * det[gene].mean(), 1)
                            if gene in det.columns else np.nan)
    rows.append(r)
sup = pd.DataFrame(rows)
front = ['label', 'n_cells', 'n_markers', 'pct_zero_markers',
         'median_markers_detected']
print(sup[front].to_string(index=False))
print()
print(sup[['label'] + [c for c in sup.columns
                       if c.startswith('pct_') and c != 'pct_zero_markers']
          ].to_string(index=False))
save_csv(sup, 'tier1_marker_support.csv', index=False)

print('\n  pct_zero_markers is the number that matters. A label where a large')
print('  share of cells detect NONE of its markers is not a cell type, it is')
print('  the argmax of a noisy score. NaN in a pct_ column means that gene is')
print('  absent from the matrix entirely.')

worst = sup.dropna(subset=['pct_zero_markers']).nlargest(3, 'pct_zero_markers')
print('\n  worst three:')
print(worst[front].to_string(index=False))


# %% ============================================================
# CELL 9 — DEFECT A: WHAT THE ASYMMETRIC CD3 GATE COSTS
# ============================================================
hr('CELL 9 — CD3-positive cells sitting in non-T labels')

if not cd3:
    print('  no CD3 gene present; cannot run.')
else:
    cd3c = mk[cd3].sum(axis=1)
    a.obs['_cd3'] = cd3c.reindex(a.obs_names).values
    rows = []
    for lin in lin_n.index:
        m = (a.obs['tier1_celltype'] == lin).values
        c = a.obs['_cd3'].values[m] >= POS_MIN
        gated = lin in ('NK cells', 'Unassigned', 'T cells')
        rows.append({
            'label': lin, 'n_cells': int(m.sum()),
            'cd3_positive': int(c.sum()),
            'pct_cd3_positive': round(100 * c.mean(), 1),
            'gate_was_applied': gated,
        })
    g = pd.DataFrame(rows).sort_values('cd3_positive', ascending=False)
    print(g.to_string(index=False))
    save_csv(g, 'tier1_cd3_by_label.csv', index=False)

    ungated = g[~g.gate_was_applied]
    n_missed = int(ungated.cd3_positive.sum())
    print(f'\n  CD3-positive cells in labels the gate never touched: '
          f'{n_missed:,}')
    print('  Each one is a cell where the score said "not T" and the CD3')
    print('  counts said "T", and 05 kept the score because of which label')
    print('  it happened to be. That is Defect A.')
    print()
    print('  CAVEAT, and it matters: CD3_MIN_COUNTS = 1 is a lenient gate.')
    print('  Ambient CD3 mRNA in a T-cell-rich well will make genuinely')
    print('  non-T cells CD3-positive at one count. So this number is an')
    print('  upper bound on the mislabelling, not a count of it. Raising the')
    print('  threshold to 2 or 3 is the way to separate the two, and that is')
    print('  a decision rather than a fix.')
    for thr in (1, 2, 3, 5):
        n = int(((a.obs['_cd3'].values >= thr)
                 & ~a.obs['tier1_celltype'].isin(
                     ['T cells', 'NK cells', 'Unassigned']).values).sum())
        print(f'    at CD3 >= {thr}: {n:,} CD3-positive cells in ungated labels')


# %% ============================================================
# CELL 10 — IS 31% B CELLS PLAUSIBLE? SPLIT IT BY TISSUE.
# ============================================================
hr('CELL 10 — lineage composition by tissue')

print('  Four of the eight libraries are lymph node, where B cells run 40 to')
print('  60% of lymphocytes. So 31% B overall is not by itself a bug. If the')
print('  B fraction is high in LN and low in PBMC, the composition is real')
print('  and Defect B is a labelling precision problem rather than a mass')
print('  mislabelling. If B is high in BOTH, something is wrong.')
print()

if TISSUE_COL in a.obs.columns:
    comp = pd.crosstab(a.obs['tier1_celltype'], a.obs[TISSUE_COL],
                       normalize='columns').mul(100).round(1)
    print(comp.to_string())
    save_csv(comp, 'tier1_composition_by_tissue.csv')

    if 'B cells' in comp.index:
        print('\n  B cells by tissue:')
        for t in comp.columns:
            print(f'    {t:<8} {comp.loc["B cells", t]:.1f}%')
        print('  Expected from the literature: PBMC 5-15%, LN 40-60%.')

    comp_lib = pd.crosstab(a.obs['tier1_celltype'], a.obs['library'],
                           normalize='columns').mul(100).round(1)
    print('\n  by library:')
    print(comp_lib.to_string())
    save_csv(comp_lib, 'tier1_composition_by_library.csv')
else:
    print(f'  no {TISSUE_COL} column; skipped')


# %% ============================================================
# CELL 11 — WERE ANY VIRAL BARCODES LOST BEFORE 05?
# ============================================================
hr('CELL 11 — viral barcodes in the prefilter object')

if not os.path.exists(PREFILTER):
    print(f'  {PREFILTER} not found; skipped')
else:
    pre = sc.read_h5ad(PREFILTER)
    if SHIV_COL in pre.obs.columns:
        pv = pre.obs[SHIV_COL].values.astype(float)
        print(f'  prefilter: {pre.n_obs:,} cells, '
              f'{int((pv > 0).sum())} viral cells, {pv.sum():,.0f} UMI')
        print(f'  annotated: {a.n_obs:,} cells, {n_v} viral cells, '
              f'{v.sum():,.0f} UMI')
        lost = int((pv > 0).sum()) - n_v
        print(f'  viral barcodes not in the annotated object: {lost}')
        if lost:
            keep = pre.obs_names[pv > 0].difference(a.obs_names)
            cc = ['library', ANIMAL_COL, TIME_COL, SHIV_COL,
                  'pct_counts_mt', 'total_counts']
            cc = [c for c in cc if c in pre.obs.columns]
            print('\n  the dropped ones:')
            print(pre.obs.loc[keep, cc].to_string())
            save_csv(pre.obs.loc[keep, cc], 'viral_barcodes_dropped.csv')
            print('\n  These were removed by the gene filter, not by QC (the')
            print('  whitelist exempts viral barcodes from QC). Recoverable')
            print('  from adata_prefilter.h5ad if 06 needs them.')
    else:
        print(f'  prefilter has no {SHIV_COL}; skipped')
    del pre


# %% ============================================================
# CELL 12 — WHAT THIS ADDS UP TO
# ============================================================
hr('DONE — 05b_viral_and_label_audit.py')
print(f'  outputs: {OUT_DIR}')
print('  nothing was modified.')
print()
print('  Read in this order:')
print('    1. viral_by_infection_status.csv  the control arm. If the rate')
print('       there is not far below the infected arm, the signal IS the')
print('       background and no lineage claim survives.')
print('    2. viral_in_impossible_lineages.csv  the ambient rate, measured')
print('       in cells that cannot be infected.')
print('    3. viral_fraction_bands.csv  how many cells clear 1% of their')
print('       transcriptome. If zero, nothing looks productively infected.')
print('    4. viral_ambient_test.csv  whether any lineage beats its own')
print('       abundance after correction.')
print('    5. tier1_marker_support.csv  whether the B label means anything.')
print('    6. tier1_cd3_by_label.csv  what the asymmetric gate costs.')
print()
print('  The B-cell question splits in two and they have different answers:')
print('    - Is the viral UMI in those cells real infection? Almost')
print('      certainly not, and the same is true of most of the T cells.')
print('      59 cells at 3 UMI each is an ambient-level signal.')
print('    - Is the B label on those cells correct? CELL 8 and CELL 10')
print('      answer that, and it is a separate problem worth fixing on its')
print('      own merits regardless of the virus.')
