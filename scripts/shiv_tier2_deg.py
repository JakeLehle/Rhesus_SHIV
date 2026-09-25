#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Rhesus SHIV - Tier 2 T Cell Cell-Level DEG (step 3b of 3)
=========================================================
Cell-level differential expression on the annotated T subset, built for the
n=5-infected reality: every contrast is WITHIN a fixed subtype and WITHIN the
infected Chinese animals, and the elite-controller contrast (n=1) is reported as
a described case with a per-animal consistency panel, NOT a significance-tested
group comparison.

Why cell-level and not pseudobulk: five infected animals is too few for a
DESeq2/edgeR animal-level test. Cell-level Wilcoxon is used for RANKING only.
Cells within an animal are not independent, so every cell-level p-value here is
pseudoreplicated and inflated; the EVIDENCE layer is the per-animal consistency
panel (does the effect hold across all four comparator animals, or is it driven
by one). Reported effects lead with effect size + consistency, not p-value.

Comparison group discipline (confounds): infection is confounded with subspecies
and sample quality (all infected are Chinese; both controls are Indian and tiny).
So the elite comparison is 40707 vs the FOUR other infected Chinese only. Indian
controls and Unknown-animal cells are excluded from every contrast here.

Robustness: percent-expressing is computed directly from the matrix via a dict
(frac_expr_dict), not scanpy's version-variable pts columns and not pandas .map
(which can silently mis-align). The gene picker relaxes the pct filter if the
strict pass is empty, so a contrast never silently collapses.

Contrasts (all cell-level, all on qc_exclude == False):
  1. ELITE (centerpiece): 40707 vs {39272,40702,41861,41862} within CD8 effector
     memory, ranked, then per-animal consistency across the 4 comparators.
  2. ELITE (secondary): same design within CD4 central memory.
  3. TRAJECTORY: Pre vs Necropsy within CD8 effector memory (39272 has no
     necropsy; the per-animal check drops it automatically).
  4. TISSUE / reservoir: LN vs PBMC within CD8 effector memory, plus a direct
     per-animal CXCR5 readout in LN cells (the follicular-exclusion probe).
  5. CONTINUOUS PROGRAMS: Th1/Th17/Tfh/Treg/Cytotoxic/Proliferating/Activation
     and gamma-delta fraction, per animal within subtype.
  6. TERMINAL EFFECTOR: compositional depletion of 40707, reported not tested.

Read-only w.r.t. objects; writes DEG CSVs + figures. Spyder cells (# %%).

Author: Jake Lehle / Kaushal Lab
Date: July 2026
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
from scipy import sparse
import warnings
warnings.filterwarnings("ignore")

WORKING_DIR = '/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV'
MERGED      = os.path.join(WORKING_DIR, 'shiv_starsolo', 'merged')
ADATA_IN    = os.path.join(MERGED, 'shiv_Tcell_annotated.h5ad')
OUT_DIR     = os.path.join(WORKING_DIR, 'annotation_output', 'tier2_DEG')

# ---- columns ----
SUBTYPE_COL = 'tier2_subtype'
QC_COL      = 'qc_exclude'
ANIMAL_COL  = 'animal_id'
TIME_COL    = 'timepoint'
TISSUE_COL  = 'tissue'
GD_COL      = 'is_gamma_delta'

# ---- groups ----
INFECTED_CHINESE = ['39272', '40702', '40707', '41861', '41862']
ELITE            = '40707'
COMPARATORS      = ['39272', '40702', '41861', '41862']

# ---- subtype targets ----
CD8_EM   = 'CD8 effector memory'
CD4_CM   = 'CD4 central memory'
CD8_TERM = 'CD8 terminal effector'
NAIVE    = 'Naive/resting T'

# ---- DEG thresholds (cell-level = ranking; effect size + consistency = evidence) ----
DEG_METHOD   = 'wilcoxon'
N_TOP        = 25          # top up / top down to carry into consistency panel
LFC_MIN      = 0.5         # |log2FC| filter (log2 0.5 = 1.41x)
PCT_MIN      = 0.10        # min fraction expressing in the higher group
FDR_MAX      = 0.05        # BH-adjusted p cutoff (cell-level, noise floor)
MIN_CELLS    = 50          # skip a per-animal group below this in a contrast
CAND_HEAD    = 250         # candidate genes per direction to compute pct on

# ---- continuous programs (per-animal, not clusters) ----
SIG_PROGRAMS = ['sig_Th1', 'sig_Th17', 'sig_Tfh', 'sig_Treg', 'sig_Cytotoxic',
                'sig_Proliferating', 'sig_Activation']
CXCR5_GENE   = 'CXCR5'

# ---- colors (hex) ----
ELITE_COL = '#F9A825'
COMP_COL  = '#90A4AE'
UP_COL    = '#D32F2F'
DN_COL    = '#1976D2'
ANIMAL_ORDER = ['39272', '40702', ELITE, '41861', '41862']

plt.rcParams.update({
    'font.size': 28, 'axes.titlesize': 34, 'axes.labelsize': 30,
    'xtick.labelsize': 22, 'ytick.labelsize': 20, 'legend.fontsize': 18,
    'figure.dpi': 100, 'pdf.fonttype': 42, 'ps.fonttype': 42,
})
os.makedirs(OUT_DIR, exist_ok=True)
sc.settings.figdir = OUT_DIR
pd.set_option('display.width', 220)
pd.set_option('display.max_columns', 40)
pd.set_option('display.max_rows', 120)


def savefig(fig, name):
    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(OUT_DIR, f'{name}.{ext}'),
                    bbox_inches='tight', dpi=300)
    print(f"    saved: {name}.pdf / .png")


def save_csv(df, name, index=True):
    df.to_csv(os.path.join(OUT_DIR, name), index=index)
    print(f"    saved: {name}")


def frac_expr_dict(ad, genes, mask):
    """Fraction of masked cells expressing each gene (>0), as a dict. Robust to
    np.matrix, empty masks, and duplicate/missing gene names."""
    g = [x for x in dict.fromkeys(genes) if x in ad.var_names]  # dedupe, keep order
    if not g or int(np.sum(mask)) == 0:
        return {x: 0.0 for x in genes}
    X = ad[mask, :][:, g].X
    X = X.toarray() if sparse.issparse(X) else np.asarray(X)
    fr = np.asarray((np.asarray(X) > 0).mean(axis=0)).ravel()
    return {gene: float(fr[i]) for i, gene in enumerate(g)}


def per_animal_expr(ad, genes, animals):
    """Mean log-norm expression per animal (rows=animals, cols=genes)."""
    g = [x for x in dict.fromkeys(genes) if x in ad.var_names]
    if not g:
        return pd.DataFrame(index=animals)
    X = ad[:, g].X
    X = X.toarray() if sparse.issparse(X) else np.asarray(X)
    df = pd.DataFrame(np.asarray(X), columns=g, index=ad.obs_names)
    df['_a'] = ad.obs[ANIMAL_COL].astype(str).values
    return df.groupby('_a')[g].mean().reindex(animals)


def consistency_panel(ad, genes, direction, elite=ELITE, comparators=COMPARATORS):
    """Mean expr per animal + consistency count. Empty-safe (keeps columns)."""
    animals = [elite] + comparators
    cols = animals + ['consistency', 'n_comp', 'direction']
    g = [x for x in genes if x in ad.var_names]
    if not g:
        return pd.DataFrame(columns=cols)
    m = per_animal_expr(ad, g, animals)          # animals x genes
    e = m.loc[elite]
    comp = m.loc[comparators]
    cons = (e > comp).sum(axis=0) if direction == 'up' else (e < comp).sum(axis=0)
    out = m.T.copy()                             # genes x animals
    out['consistency'] = cons.reindex(out.index).values
    out['n_comp'] = len(comparators)
    out['direction'] = direction
    return out


def run_elite_deg(ad_subtype, tag):
    """Cell-level Wilcoxon: elite vs pooled comparators, + consistency panel."""
    a = ad_subtype[ad_subtype.obs[ANIMAL_COL].astype(str).isin(INFECTED_CHINESE)].copy()
    a.obs['_grp'] = pd.Categorical(
        np.where(a.obs[ANIMAL_COL].astype(str) == ELITE, 'elite', 'other'))
    n_e = int((a.obs['_grp'] == 'elite').sum())
    n_o = int((a.obs['_grp'] == 'other').sum())
    print(f"  [{tag}] elite {n_e:,} vs other-4 {n_o:,} cells")
    if n_e < MIN_CELLS or n_o < MIN_CELLS:
        print(f"  [{tag}] SKIP: a group is below {MIN_CELLS} cells.")
        return None

    sc.tl.rank_genes_groups(a, groupby='_grp', groups=['elite'], reference='other',
                            method=DEG_METHOD, use_raw=False)
    df = sc.get.rank_genes_groups_df(a, group='elite').rename(
        columns={'names': 'gene', 'logfoldchanges': 'log2fc'})
    df = df.dropna(subset=['log2fc'])
    df['pvals_adj'] = df['pvals_adj'].fillna(1.0)
    save_csv(df, f'deg_elite_{tag}_full.csv', index=False)

    elite_mask = (a.obs['_grp'] == 'elite').values
    other_mask = ~elite_mask

    # funnel: effect size -> FDR -> fraction expressing (in the higher group)
    lfc = df[df['log2fc'].abs() >= LFC_MIN]
    sig = lfc[lfc['pvals_adj'] < FDR_MAX].copy()
    pe = frac_expr_dict(a, sig['gene'].tolist(), elite_mask)
    po = frac_expr_dict(a, sig['gene'].tolist(), other_mask)
    sig['pct_elite'] = sig['gene'].map(pe).astype(float)
    sig['pct_other'] = sig['gene'].map(po).astype(float)
    up = sig[(sig['log2fc'] >= LFC_MIN) & (sig['pct_elite'] >= PCT_MIN)] \
        .sort_values('log2fc', ascending=False)
    dn = sig[(sig['log2fc'] <= -LFC_MIN) & (sig['pct_other'] >= PCT_MIN)] \
        .sort_values('log2fc')
    print(f"  [{tag}] funnel: all {len(df):,} -> |log2fc|>={LFC_MIN} {len(lfc):,} "
          f"-> FDR<{FDR_MAX} {len(sig):,} -> pct>={PCT_MIN} up {len(up)} / down {len(dn)}")
    save_csv(pd.concat([up, dn]).round(4), f'deg_elite_{tag}_significant.csv', index=False)
    up = up.head(N_TOP); dn = dn.head(N_TOP)
    if len(up) + len(dn) == 0:
        print(f"  [{tag}] no genes pass |log2fc|>={LFC_MIN}, FDR<{FDR_MAX}, pct>={PCT_MIN}.")
        return None

    con = pd.concat([consistency_panel(a, up['gene'].tolist(), 'up'),
                     consistency_panel(a, dn['gene'].tolist(), 'down')], axis=0)
    if con.empty or 'direction' not in con.columns:
        print(f"  [{tag}] no genes to report after filtering.")
        return None

    eff = sig.set_index('gene')
    con['log2fc']    = eff['log2fc'].reindex(con.index).values
    con['pvals_adj'] = eff['pvals_adj'].reindex(con.index).values
    con['pct_elite'] = eff['pct_elite'].reindex(con.index).values
    con['pct_other'] = eff['pct_other'].reindex(con.index).values
    con = con.sort_values(['direction', 'consistency', 'log2fc'],
                          ascending=[True, False, False])

    print(f"  [{tag}] top consistent genes (consistency = k / {len(COMPARATORS)}, "
          f"FDR<{FDR_MAX}):")
    show = con[[ELITE] + COMPARATORS + ['consistency', 'log2fc', 'pct_elite',
                                        'pct_other', 'pvals_adj']].round(4)
    print(show.to_string())
    save_csv(con.round(4), f'deg_elite_{tag}_consistency.csv')
    return con


# %% ============================================================
# CELL 2 - LOAD + FILTER TO BIOLOGICAL T
# ============================================================
print("=" * 80)
print("CELL 2 - LOAD + FILTER")
print("=" * 80)

sub = sc.read_h5ad(ADATA_IN)
print(f"  loaded: {sub.n_obs:,} cells; X.max()={float(sub.X.max()):.2f} "
      f"({'log-norm' if sub.X.max() < 50 else 'raw'})")

bio = sub[~sub.obs[QC_COL].astype(bool)].copy()
print(f"  biological T (qc_exclude False): {bio.n_obs:,}")
print("  by subtype:")
print(bio.obs[SUBTYPE_COL].value_counts().to_string())
print(f"\n  infected-Chinese biological T: "
      f"{int(bio.obs[ANIMAL_COL].astype(str).isin(INFECTED_CHINESE).sum()):,}")


# %% ============================================================
# CELL 3 - CONTRAST 1: ELITE within CD8 EFFECTOR MEMORY (centerpiece)
# ============================================================
print("\n" + "=" * 80)
print("CELL 3 - ELITE DEG within CD8 EFFECTOR MEMORY")
print("=" * 80)
cd8em = bio[bio.obs[SUBTYPE_COL] == CD8_EM].copy()
con_cd8em = run_elite_deg(cd8em, 'CD8em')


# %% ============================================================
# CELL 4 - CONTRAST 1b: ELITE within CD4 CENTRAL MEMORY (secondary)
# ============================================================
print("\n" + "=" * 80)
print("CELL 4 - ELITE DEG within CD4 CENTRAL MEMORY")
print("=" * 80)
cd4cm = bio[bio.obs[SUBTYPE_COL] == CD4_CM].copy()
con_cd4cm = run_elite_deg(cd4cm, 'CD4cm')


# %% ============================================================
# CELL 5 - CONTRAST 2: TIMEPOINT (Pre vs Necropsy) within CD8 EM
# ============================================================
print("\n" + "=" * 80)
print("CELL 5 - TRAJECTORY: Pre vs Necropsy within CD8 EM")
print("=" * 80)
tj = cd8em[cd8em.obs[TIME_COL].astype(str).isin(['Pre', 'Necropsy'])
           & cd8em.obs[ANIMAL_COL].astype(str).isin(INFECTED_CHINESE)].copy()
print("  cells per timepoint x animal (watch 39272 = no necropsy):")
print(pd.crosstab(tj.obs[ANIMAL_COL].astype(str), tj.obs[TIME_COL].astype(str)).to_string())

n_pre = int((tj.obs[TIME_COL] == 'Pre').sum())
n_nec = int((tj.obs[TIME_COL] == 'Necropsy').sum())
if n_pre >= MIN_CELLS and n_nec >= MIN_CELLS:
    sc.tl.rank_genes_groups(tj, groupby=TIME_COL, groups=['Necropsy'],
                            reference='Pre', method=DEG_METHOD, use_raw=False)
    tjdf = sc.get.rank_genes_groups_df(tj, group='Necropsy').rename(
        columns={'names': 'gene', 'logfoldchanges': 'log2fc'}).dropna(subset=['log2fc'])
    save_csv(tjdf, 'deg_trajectory_CD8em_necropsy_vs_pre_full.csv', index=False)

    nec_mask = (tj.obs[TIME_COL] == 'Necropsy').values
    pre_mask = (tj.obs[TIME_COL] == 'Pre').values
    tjdf['pvals_adj'] = tjdf['pvals_adj'].fillna(1.0)
    sig = tjdf[(tjdf['log2fc'].abs() >= LFC_MIN) & (tjdf['pvals_adj'] < FDR_MAX)]
    up_c = sig[sig['log2fc'] >= LFC_MIN].sort_values('log2fc', ascending=False)
    dn_c = sig[sig['log2fc'] <= -LFC_MIN].sort_values('log2fc')
    pct_nec = frac_expr_dict(tj, up_c['gene'].tolist(), nec_mask)
    pct_pre = frac_expr_dict(tj, dn_c['gene'].tolist(), pre_mask)
    tj_up = up_c[up_c['gene'].map(pct_nec).fillna(0) >= PCT_MIN].head(N_TOP)
    tj_dn = dn_c[dn_c['gene'].map(pct_pre).fillna(0) >= PCT_MIN].head(N_TOP)
    print(f"  trajectory funnel: |log2fc|>={LFC_MIN} & FDR<{FDR_MAX} & pct>={PCT_MIN} "
          f"-> necropsy-up {len(tj_up)} / necropsy-down {len(tj_dn)}")

    both = [a for a in INFECTED_CHINESE
            if {'Pre', 'Necropsy'} <= set(tj.obs.loc[tj.obs[ANIMAL_COL] == a, TIME_COL])]
    print(f"  animals with both timepoints (per-animal check set): {both}")
    genes = tj_up['gene'].tolist() + tj_dn['gene'].tolist()
    rows = []
    for a in both:
        aa = tj[tj.obs[ANIMAL_COL] == a]
        for tp in ('Pre', 'Necropsy'):
            m = per_animal_expr(aa[aa.obs[TIME_COL] == tp], genes, [a])
            if not m.empty:
                r = m.iloc[0].to_dict(); r['animal'] = a; r['timepoint'] = tp
                rows.append(r)
    if rows:
        pa = pd.DataFrame(rows).set_index(['animal', 'timepoint'])
        save_csv(pa.round(4), 'deg_trajectory_CD8em_per_animal.csv')
        print("  per-animal Pre/Necropsy means saved for the top trajectory genes.")
    print("  top necropsy-up genes:", tj_up['gene'].head(12).tolist())
    print("  top necropsy-down genes:", tj_dn['gene'].head(12).tolist())
else:
    print(f"  SKIP: Pre {n_pre}, Necropsy {n_nec} (need >= {MIN_CELLS} each).")


# %% ============================================================
# CELL 6 - CONTRAST 3: TISSUE (LN vs PBMC) within CD8 EM + CXCR5 readout
# ============================================================
print("\n" + "=" * 80)
print("CELL 6 - TISSUE: LN vs PBMC within CD8 EM + CXCR5 reservoir probe")
print("=" * 80)
ts = cd8em[cd8em.obs[ANIMAL_COL].astype(str).isin(INFECTED_CHINESE)].copy()
print("  cells per tissue:")
print(ts.obs[TISSUE_COL].value_counts().to_string())
n_ln = int((ts.obs[TISSUE_COL] == 'LN').sum())
n_pb = int((ts.obs[TISSUE_COL] == 'PBMC').sum())
if n_ln >= MIN_CELLS and n_pb >= MIN_CELLS:
    sc.tl.rank_genes_groups(ts, groupby=TISSUE_COL, groups=['LN'], reference='PBMC',
                            method=DEG_METHOD, use_raw=False)
    tsdf = sc.get.rank_genes_groups_df(ts, group='LN').rename(
        columns={'names': 'gene', 'logfoldchanges': 'log2fc'}).dropna(subset=['log2fc'])
    tsdf['pvals_adj'] = tsdf['pvals_adj'].fillna(1.0)
    save_csv(tsdf, 'deg_tissue_CD8em_LN_vs_PBMC_full.csv', index=False)
    ts_sig = tsdf[(tsdf['log2fc'].abs() >= LFC_MIN) & (tsdf['pvals_adj'] < FDR_MAX)]
    ln_mask = (ts.obs[TISSUE_COL] == 'LN').values
    pb_mask = (ts.obs[TISSUE_COL] == 'PBMC').values
    up_c = ts_sig[ts_sig['log2fc'] >= LFC_MIN]
    dn_c = ts_sig[ts_sig['log2fc'] <= -LFC_MIN]
    pe = frac_expr_dict(ts, up_c['gene'].tolist(), ln_mask)
    po = frac_expr_dict(ts, dn_c['gene'].tolist(), pb_mask)
    ts_up = up_c[up_c['gene'].map(pe).fillna(0) >= PCT_MIN].sort_values('log2fc', ascending=False)
    ts_dn = dn_c[dn_c['gene'].map(po).fillna(0) >= PCT_MIN].sort_values('log2fc')
    save_csv(pd.concat([ts_up, ts_dn]).round(4), 'deg_tissue_CD8em_significant.csv', index=False)
    print(f"  tissue funnel: |log2fc|>={LFC_MIN} & FDR<{FDR_MAX} & pct>={PCT_MIN} "
          f"-> LN-up {len(ts_up)} / PBMC-up {len(ts_dn)}")
    print("  top LN-up genes:", ts_up['gene'].head(12).tolist())
else:
    print(f"  SKIP tissue DEG: LN {n_ln}, PBMC {n_pb}.")

# CXCR5 per-animal readout in LN cells, by subtype (direct, not a cluster test)
print("\n  CXCR5 in LN cells, per animal x subtype (fraction+ / mean expr):")
ln_bio = bio[(bio.obs[TISSUE_COL] == 'LN')
             & bio.obs[ANIMAL_COL].astype(str).isin(INFECTED_CHINESE)].copy()
if CXCR5_GENE in ln_bio.var_names:
    x = ln_bio[:, CXCR5_GENE].X
    x = x.toarray().ravel() if sparse.issparse(x) else np.asarray(x).ravel()
    ln_df = pd.DataFrame({
        'animal': ln_bio.obs[ANIMAL_COL].astype(str).values,
        'subtype': ln_bio.obs[SUBTYPE_COL].astype(str).values,
        'cxcr5': x, 'pos': (x > 0).astype(int)})
    cx = ln_df.groupby(['subtype', 'animal']).agg(
        n=('cxcr5', 'size'), pct_pos=('pos', 'mean'), mean_expr=('cxcr5', 'mean'))
    cx['pct_pos'] = (100 * cx['pct_pos']).round(1)
    cx['mean_expr'] = cx['mean_expr'].round(3)
    print(cx.to_string())
    save_csv(cx, 'cxcr5_LN_per_animal_by_subtype.csv')
else:
    print(f"  WARN {CXCR5_GENE} not in var_names.")


# %% ============================================================
# CELL 7 - CONTRAST 5: CONTINUOUS PROGRAMS per animal (within subtype)
# ============================================================
print("\n" + "=" * 80)
print("CELL 7 - CONTINUOUS PROGRAMS per animal")
print("=" * 80)
progs = [c for c in SIG_PROGRAMS if c in bio.obs.columns]
for stype in (CD8_EM, CD4_CM):
    a = bio[(bio.obs[SUBTYPE_COL] == stype)
            & bio.obs[ANIMAL_COL].astype(str).isin(INFECTED_CHINESE)]
    if a.n_obs == 0:
        continue
    tbl = a.obs.groupby(ANIMAL_COL, observed=True)[progs].mean()
    if GD_COL in a.obs.columns:
        tbl['gamma_delta_frac'] = a.obs.groupby(ANIMAL_COL, observed=True)[GD_COL].mean()
    tbl = tbl.reindex([x for x in ANIMAL_ORDER if x in tbl.index])
    tbl.columns = [c.replace('sig_', '') for c in tbl.columns]
    print(f"\n  [{stype}] mean program score per animal (40707 = elite):")
    print(tbl.round(3).to_string())
    save_csv(tbl.round(4), f'programs_per_animal_{stype.replace(" ", "_")}.csv')


# %% ============================================================
# CELL 8 - CONTRAST 6: TERMINAL EFFECTOR compositional depletion (reported)
# ============================================================
print("\n" + "=" * 80)
print("CELL 8 - TERMINAL EFFECTOR compositional depletion")
print("=" * 80)
bio_ic = bio[bio.obs[ANIMAL_COL].astype(str).isin(INFECTED_CHINESE)]
comp = pd.crosstab(bio_ic.obs[ANIMAL_COL].astype(str), bio_ic.obs[SUBTYPE_COL])
term_frac = (comp.get(CD8_TERM, pd.Series(0, index=comp.index))
             / comp.sum(axis=1)).reindex([x for x in ANIMAL_ORDER if x in comp.index])
print("  terminal-effector fraction of T per infected animal:")
print((100 * term_frac).round(2).to_string())
term_n = int(comp.loc[ELITE, CD8_TERM]) if (ELITE in comp.index and CD8_TERM in comp.columns) else 0
print(f"\n  40707 = {100*term_frac.get(ELITE, np.nan):.2f}%  vs "
      f"comparator mean {100*term_frac.drop(index=ELITE, errors='ignore').mean():.2f}% "
      f"-> depletion (n={term_n} cells, too few to test per-cell state).")
save_csv(term_frac.round(4).to_frame('terminal_frac'),
         'terminal_effector_fraction_per_animal.csv')


# %% ============================================================
# CELL 9 - FIGURES
# ============================================================
print("\n" + "=" * 80)
print("CELL 9 - FIGURES")
print("=" * 80)


def consistency_heatmap(con, tag, title):
    if con is None or con.empty:
        print(f"  [{tag}] no consistency panel to plot.")
        return
    animals = [c for c in ANIMAL_ORDER if c in con.columns]
    M = con[animals].astype(float)
    Z = M.sub(M.mean(axis=1), axis=0).div(M.std(axis=1).replace(0, np.nan), axis=0).fillna(0)
    fig, ax = plt.subplots(figsize=(1.1 * len(animals) + 6, 0.42 * len(Z) + 4))
    im = ax.imshow(Z.values, cmap='RdBu_r', vmin=-1.8, vmax=1.8, aspect='auto')
    ax.set_xticks(range(len(animals)))
    ax.set_xticklabels([f"{a}*" if a == ELITE else a for a in animals])
    ax.set_yticks(range(len(Z)))
    ax.set_yticklabels([f"{g}  ({int(c)}/{int(n)})"
                        for g, c, n in zip(Z.index, con['consistency'], con['n_comp'])])
    if ELITE in animals:
        ecol = animals.index(ELITE)
        ax.add_patch(plt.Rectangle((ecol - 0.5, -0.5), 1, len(Z),
                                   fill=False, edgecolor='#F9A825', lw=4))
    ax.set_title(title)
    ax.set_xlabel('animal (* = elite)'); ax.set_ylabel('gene (consistency)')
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label='z (per gene)')
    plt.tight_layout()
    savefig(fig, f'deg_elite_{tag}_consistency_heatmap')
    plt.close(fig)


consistency_heatmap(con_cd8em, 'CD8em',
                    '40707 vs 4 comparators: CD8 EM\n(consistency in parentheses)')
consistency_heatmap(con_cd4cm, 'CD4cm',
                    '40707 vs 4 comparators: CD4 CM\n(consistency in parentheses)')

# CXCR5 per-animal bar in LN (CD8 EM + CD4 CM)
try:
    cxp = pd.read_csv(os.path.join(OUT_DIR, 'cxcr5_LN_per_animal_by_subtype.csv'))
    fig, ax = plt.subplots(figsize=(15, 9))
    stypes = [CD8_EM, CD4_CM]
    width = 0.38
    for j, st in enumerate(stypes):
        d = cxp[cxp['subtype'] == st].set_index('animal').reindex(
            [a for a in ANIMAL_ORDER if a in cxp['animal'].values])
        x = np.arange(len(d))
        ax.bar(x + (j - 0.5) * width, d['pct_pos'].values, width=width,
               label=st, color=[UP_COL, DN_COL][j], alpha=0.85)
    ax.set_xticks(np.arange(len(ANIMAL_ORDER)))
    ax.set_xticklabels([f"{a}*" if a == ELITE else a for a in ANIMAL_ORDER],
                       rotation=45, ha='right')
    ax.set_ylabel('% CXCR5+ (LN cells)')
    ax.set_title('CXCR5 in LN by animal and subtype\n(* = elite controller)')
    ax.legend(); ax.spines[['top', 'right']].set_visible(False)
    plt.tight_layout()
    savefig(fig, 'cxcr5_LN_per_animal')
    plt.close(fig)
except Exception as e:
    print(f"  WARN CXCR5 figure skipped: {e}")

# terminal-effector fraction bar
try:
    fig, ax = plt.subplots(figsize=(12, 9))
    tf = (100 * term_frac).reindex([a for a in ANIMAL_ORDER if a in term_frac.index])
    colors = [ELITE_COL if a == ELITE else COMP_COL for a in tf.index]
    ax.bar(range(len(tf)), tf.values, color=colors)
    ax.set_xticks(range(len(tf)))
    ax.set_xticklabels([f"{a}*" if a == ELITE else a for a in tf.index],
                       rotation=45, ha='right')
    ax.set_ylabel('% terminal effector of T')
    ax.set_title('CD8 terminal effector depletion in 40707\n(* = elite controller)')
    ax.spines[['top', 'right']].set_visible(False)
    plt.tight_layout()
    savefig(fig, 'terminal_effector_fraction')
    plt.close(fig)
except Exception as e:
    print(f"  WARN terminal figure skipped: {e}")


# %% ============================================================
# CELL 10 - SUMMARY
# ============================================================
print("\n" + "=" * 80)
print("CELL 10 - SUMMARY")
print("=" * 80)
print(f"""
  Outputs in {OUT_DIR}:
   - deg_elite_CD8em_significant.csv = ALL genes passing |log2fc|>=0.5, FDR<0.05,
     pct>=0.10; the _consistency.csv is the top set with the per-animal panel.
   - deg_elite_CD8em_consistency.csv (+ heatmap)  <- CENTERPIECE. Read the
     consistency column: {len(COMPARATORS)}/{len(COMPARATORS)} = 40707 is on the same side vs EVERY
     comparator (credible); lower = driven by a subset. Effect size (log2fc) and
     pct back it; the cell-level p is ranking-only.
   - deg_elite_CD4cm_consistency.csv (+ heatmap)  <- secondary elite contrast.
   - deg_trajectory_CD8em_*  Pre vs necropsy, with the per-animal check (39272
     has no necropsy, so it is out of the paired set).
   - deg_tissue_CD8em_LN_vs_PBMC_full.csv + cxcr5_LN_per_animal_by_subtype.csv
     (+ bar)  <- reservoir / follicular-exclusion probe, direct per-animal CXCR5.
   - programs_per_animal_*  Th1/Th17/Tfh/Treg/Cytotoxic/Prolif/Activation +
     gamma-delta fraction per animal (continuous, not clusters).
   - terminal_effector_fraction_per_animal.csv (+ bar)  <- compositional only.

  FRAMING for July 30: elite findings are a DESCRIBED CASE (n=1) supported by
  per-animal consistency, not group significance tests. Lead with effect size +
  consistency. Terminal effector is depletion, reported not tested.
""")
print("=" * 80)
print("TIER 2 CELL-LEVEL DEG COMPLETE (step 3b).")
print("=" * 80)
