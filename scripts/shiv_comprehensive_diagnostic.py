#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Rhesus SHIV — Comprehensive Read-Only Diagnostic (Object Level)
===============================================================
Reads the canonical embedded object and prints everything we need to make
decisions on four horizon items, WITHOUT changing a single label or coordinate.
Nothing is written back to the object; only CSV summaries are exported.

Modules
-------
  0  Object sanity      : columns / obsm / layers / uns / var flags present?
  1  Paired-receptor 109: the TRG+TRD-paired cells NOT labeled T. Prior label,
                          cluster, CD3 status, B/myeloid doublet flags. Settles
                          the open paired-pass decision (flip clean, quarantine
                          suspicious).
  2  Clonotype feasibility: clonotype-definable (paired) T cells per animal /
                          timepoint, 40707 broken out CD4 vs CD8 (RNA and ADT).
                          Is the 40707 CD8-expansion story powered?
  3  Plotting inputs    : all 11 ADT tags in protein_clr (9 clean vs CD163/CD169
                          blank, verified empirically); RIRA RNA overlay markers
                          present in var_names?
  4  Low-purity clusters: the clusters below the purity cutoff, their dominant
                          labels, marker means, ADT profile, data-driven top
                          markers. Scopes the Tier 2 recluster.

Design
------
  - Read-only. No obs/var/obsm/label is modified and the object is NOT re-saved.
  - Print-first + CSV export (no figures; figures live in the plotting script).
  - Warn-and-continue: every module checks for its columns/genes and skips with
    a WARN if something is missing rather than hard-stopping.
  - All tunable knobs live in Cell 1.
  - Spyder cells (# %%). Runs headless under SLURM (sc_pre env) or interactively.

Author: Jake Lehle / Kaushal Lab
Date: July 2026
"""

# %% ============================================================
# CELL 1 — CONFIG AND DIALS
# ============================================================
import os
os.environ.setdefault('MPLBACKEND', 'Agg')   # headless-safe even though no plots
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
import warnings
warnings.filterwarnings("ignore")

# --- Paths ---
WORKING_DIR = '/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV'
ADATA_IN    = os.path.join(WORKING_DIR, 'shiv_starsolo', 'merged',
                           'shiv_host_final_embedded.h5ad')
OUT_DIR     = os.path.join(WORKING_DIR, 'annotation_output',
                           'comprehensive_diagnostic')

# --- Key obs / structure columns (presence-checked before use) ---
COL_CELLTYPE     = 'tier1_celltype'          # corrected, post-TCR, canonical
COL_CELLTYPE_PRE = 'tier1_celltype_pre_tcr'  # pre-TCR snapshot (VDJ script)
COL_SUBTYPE      = 'receptor_subtype'
COL_PAIRING      = 'chain_pairing'
COL_HASTCR       = 'has_tcr'
COL_CLUSTERS     = 'clusters'
COL_CD3POS       = 'dx_cd3_pos'              # raw-count CD3+ flag (tier1 script)
ANIMAL_COL       = 'animal_id'
TISSUE_COL       = 'tissue'
TIME_COL         = 'timepoint'
COUNTS_LAYER     = 'counts'
ELITE            = '40707'

# --- Gamma-delta paired subtype string (auto-detected if this exact one absent) ---
GD_SUBTYPE = 'TRG+TRD'

# --- chain_pairing categories that mean a complete VJ+VDJ pair (clonotype-definable) ---
DEFINABLE_PAIRINGS = ['single pair', 'extra VJ', 'extra VDJ', 'two full chains']

# --- Thresholds ---
PURITY_CUTOFF  = 0.6     # clusters below this are flagged for Tier 2
CD3_MIN_COUNTS = 1       # CD3+ = >= this many summed CD3 counts (raw)
POS_MIN_COUNTS = 1       # gene "detected" = >= this many raw counts

# --- Marker panels (Mmul10 symbols; presence-checked, warn-and-continue) ---
CD3_GENES     = ['CD3D', 'CD3E', 'CD3G']
B_GENES       = ['CD79A', 'CD79B', 'MS4A1', 'BANK1', 'PAX5']
MYELOID_GENES = ['LYZ', 'CD14', 'S100A8', 'S100A9', 'VCAN', 'FCGR3']
CD8_GENES     = ['CD8A', 'CD8B']
CD4_GENES     = ['CD4']

# RIRA (Mahyari et al.) overlay panel: plotting-readiness + Tier 2 scoping
RIRA_PANEL = {
    'T':       ['CD3D', 'CD3E', 'CD3G', 'CD4', 'CD8A', 'CD8B',
                'IL7R', 'CCR7', 'SELL', 'GZMK'],
    'NK':      ['KLRD1', 'KLRC3', 'NCAM1', 'NCR1', 'KLRF1', 'GNLY', 'NKG7'],
    'B':       ['CD79A', 'CD79B', 'MS4A1', 'CD19', 'BANK1'],
    'Plasma':  ['JCHAIN', 'MZB1', 'XBP1', 'PRDM1'],
    'Myeloid': ['LYZ', 'CD14', 'FCGR3', 'CD68', 'CD163', 'C1QA'],
    'DC':      ['FLT3', 'CLEC9A', 'IRF8', 'IL3RA', 'CLEC4C'],
    'Prolif':  ['MKI67', 'TOP2A'],
}

# ADT tags to profile per low-purity cluster (whatever subset is present)
ADT_PROFILE_TAGS = ['CD4', 'CD8', 'CD11b', 'CD68', 'CD206', 'CD28', 'CD69',
                    'CX3CR1', 'CD1c', 'CD163', 'CD169']

RUN_RANK_GENES = True    # data-driven top markers for low-purity clusters
N_TOP_MARKERS  = 15

os.makedirs(OUT_DIR, exist_ok=True)
pd.set_option('display.width', 170)
pd.set_option('display.max_columns', 50)
pd.set_option('display.max_rows', 120)


def hr(title):
    print("\n" + "=" * 80)
    print(title)
    print("=" * 80)


def save_csv(df, name):
    p = os.path.join(OUT_DIR, name)
    df.to_csv(p)
    print(f"    saved: {name}")


# %% ============================================================
# CELL 2 — LOAD + MODULE 0 (object sanity)
# ============================================================
hr("CELL 2 — LOAD + MODULE 0 (object sanity)")

adata = sc.read_h5ad(ADATA_IN)
print(f"  {ADATA_IN}")
print(f"  {adata.n_obs:,} cells x {adata.n_vars:,} genes")


def present(genes):
    return [g for g in genes if g in adata.var_names]


def gene_vec(gene, layer=None):
    sub = adata[:, gene]
    X = sub.layers[layer] if layer else sub.X
    return X.toarray().ravel() if sparse.issparse(X) else np.asarray(X).ravel()


def set_sum(genes, layer=None):
    genes = present(genes)
    if not genes:
        return np.zeros(adata.n_obs)
    sub = adata[:, genes]
    M = sub.layers[layer] if layer else sub.X
    M = M.toarray() if sparse.issparse(M) else np.asarray(M)
    return M.sum(axis=1)


def resolve_adt(tag, names):
    """Map a short ADT name ('CD8') to the full stored tag ('CD8_TotalSeqC').
    Exact match first, then token-before-underscore, then prefix (CD28 ->
    CD28.2_TotalSeqC). Returns None if nothing matches (warn-and-continue)."""
    if tag in names:
        return tag
    short_map = {}
    for full in names:
        short_map.setdefault(full.split('_')[0], full)
    if tag in short_map:
        return short_map[tag]
    for full in names:
        if full.startswith(tag):
            return full
    return None


have_counts = COUNTS_LAYER in adata.layers
print(f"  counts layer present: {have_counts}")

# X state (expect log-normalized)
xmax = float(adata.X.max())
print(f"  X.max() = {xmax:.2f} ({'looks raw' if xmax > 50 else 'looks log-normalized'})")

# obs columns
print("\n  --- obs columns present? ---")
key_cols = [COL_CELLTYPE, COL_CELLTYPE_PRE, COL_SUBTYPE, COL_PAIRING, COL_HASTCR,
            COL_CLUSTERS, COL_CD3POS, ANIMAL_COL, TISSUE_COL, TIME_COL,
            'is_gamma_delta', 'nk_tcr_negative', 'nk_provisional',
            'shiv_umi', 'shiv_pos']
for c in key_cols:
    print(f"    {'OK ' if c in adata.obs.columns else 'MISSING':<8}{c}")

# obsm / layers / uns
print(f"\n  obsm keys : {list(adata.obsm.keys())}")
print(f"  layers    : {list(adata.layers.keys())}")
print(f"  uns keys  : {list(adata.uns.keys())}")
if 'adt_names' in adata.uns:
    adt_names = list(adata.uns['adt_names'])
    print(f"  adt_names ({len(adt_names)}): {adt_names}")
else:
    adt_names = []
    print("  WARN: uns['adt_names'] missing")

# var flags
print("\n  --- var flags ---")
for f in ['highly_variable', 'shiv', 'mt']:
    if f in adata.var.columns:
        print(f"    {f:<16}: {int(adata.var[f].astype(bool).sum()):,} genes")
    else:
        print(f"    {f:<16}: MISSING")

# headline categorical distributions
for c in [COL_CELLTYPE, COL_SUBTYPE, COL_PAIRING]:
    if c in adata.obs.columns:
        print(f"\n  {c} distribution:")
        print(adata.obs[c].value_counts(dropna=False).to_string())

if COL_HASTCR in adata.obs.columns:
    n_tcr = int(adata.obs[COL_HASTCR].astype(bool).sum())
    print(f"\n  has_tcr: {n_tcr:,} ({100*n_tcr/adata.n_obs:.1f}%)")

# SHIV sanity
if 'shiv_pos' in adata.obs.columns:
    print(f"  SHIV+ cells (shiv_pos): {int(adata.obs['shiv_pos'].astype(bool).sum()):,}")
if 'shiv_umi' in adata.obs.columns:
    su = pd.to_numeric(adata.obs['shiv_umi'], errors='coerce').fillna(0)
    print(f"  SHIV UMI: max {int(su.max())}, cells >=2 UMI {int((su >= 2).sum()):,}")


# %% ============================================================
# CELL 3 — MODULE 1 (paired-receptor 109: flip vs quarantine)
# ============================================================
hr("CELL 3 — MODULE 1 (paired TRG+TRD not labeled T)")

if COL_SUBTYPE not in adata.obs.columns or COL_CELLTYPE not in adata.obs.columns:
    print("  WARN: receptor_subtype or tier1_celltype missing; skipping Module 1.")
else:
    subtypes = adata.obs[COL_SUBTYPE].astype(str)
    print("  receptor_subtype values present:")
    print(subtypes.value_counts(dropna=False).to_string())

    # resolve the gamma-delta paired string
    gd_str = GD_SUBTYPE
    if gd_str not in set(subtypes):
        cand = [v for v in subtypes.unique()
                if ('TRG' in v and 'TRD' in v) or 'gamma' in v.lower()]
        if cand:
            gd_str = cand[0]
            print(f"\n  NOTE: '{GD_SUBTYPE}' not found; using detected '{gd_str}'")
        else:
            print(f"\n  WARN: no paired gamma-delta subtype found; skipping Module 1.")
            gd_str = None

    if gd_str is not None:
        is_t = adata.obs[COL_CELLTYPE].astype(str) == 'T cells'
        gd_paired = subtypes == gd_str
        the109 = gd_paired & ~is_t

        print(f"\n  total {gd_str} paired cells      : {int(gd_paired.sum()):,}")
        print(f"  ...already labeled T (gd T)    : {int((gd_paired & is_t).sum()):,}")
        print(f"  ...NOT labeled T (the target)  : {int(the109.sum()):,}")

        idx = the109.values
        sub = adata[idx]

        print("\n  prior tier1_celltype of the non-T paired cells:")
        print(sub.obs[COL_CELLTYPE].astype(str).value_counts().to_string())

        if COL_CLUSTERS in adata.obs.columns:
            print("\n  Leiden cluster membership:")
            print(sub.obs[COL_CLUSTERS].astype(str).value_counts().to_string())

        # CD3 status
        if COL_CD3POS in adata.obs.columns:
            cd3_pos = adata.obs[COL_CD3POS].astype(bool).values
            cd3_src = "dx_cd3_pos"
        else:
            cd3_sum = set_sum(CD3_GENES, layer=COUNTS_LAYER if have_counts else None)
            cd3_pos = cd3_sum >= CD3_MIN_COUNTS
            cd3_src = f"CD3 counts >= {CD3_MIN_COUNTS}"
        print(f"\n  CD3+ among the paired non-T cells ({cd3_src}): "
              f"{int(cd3_pos[idx].sum())}/{int(idx.sum())} "
              f"({100*cd3_pos[idx].mean():.1f}%)")

        # conflict flags: B / myeloid / erythrocyte-platelet markers on raw counts
        layer = COUNTS_LAYER if have_counts else None
        b_sum   = set_sum(['CD79A', 'MS4A1'], layer=layer)
        mye_sum = set_sum(['LYZ', 'CD14'],    layer=layer)
        ery_sum = set_sum(['HBB', 'PPBP'],    layer=layer)   # erythrocyte + platelet
        b_flag   = b_sum   >= POS_MIN_COUNTS
        mye_flag = mye_sum >= POS_MIN_COUNTS
        ery_flag = ery_sum >= POS_MIN_COUNTS
        print(f"  B-marker+ (CD79A/MS4A1) among them  : "
              f"{int(b_flag[idx].sum())} ({100*b_flag[idx].mean():.1f}%)")
        print(f"  myeloid+ (LYZ/CD14) among them      : "
              f"{int(mye_flag[idx].sum())} ({100*mye_flag[idx].mean():.1f}%)")
        print(f"  ery/platelet+ (HBB/PPBP) among them : "
              f"{int(ery_flag[idx].sum())} ({100*ery_flag[idx].mean():.1f}%)")

        # QUARANTINE on B / myeloid conflict ONLY. Ery/platelet (HBB/PPBP) is
        # ambient soup, not a lineage doublet, so it stays visible in the CSV
        # but does NOT disqualify a cell that carries a paired productive
        # TRG+TRD receptor. The paired receptor IS the lineage evidence.
        quarantine = b_flag | mye_flag
        flip = ~quarantine
        flip_confident   = flip &  cd3_pos     # real receptor + CD3 present
        flip_provisional = flip & ~cd3_pos     # real receptor, CD3 dropout -> provisional gd

        n_q    = int((the109.values & quarantine).sum())
        n_conf = int((the109.values & flip_confident).sum())
        n_prov = int((the109.values & flip_provisional).sum())
        n_ery  = int((the109.values & ery_flag).sum())
        print(f"\n  CANDIDATE CALL on the {int(idx.sum())} paired non-T cells")
        print(f"  (quarantine on B/myeloid only; ery/platelet visible, not disqualifying):")
        print(f"    B/myeloid+ (doublet signature)        -> QUARANTINE            : {n_q}")
        print(f"    clean, CD3+  (receptor + CD3)          -> FLIP to T, confident  : {n_conf}")
        print(f"    clean, CD3-  (receptor, CD3 dropout)   -> FLIP to T, provisional: {n_prov}")
        print(f"    (ery/platelet+ among the above, informational only)  : {n_ery}")

        n_gd_floor = int((gd_paired & is_t).sum())
        print(f"\n  strict gamma-delta floor (already T)      : {n_gd_floor}")
        print(f"  + paired-pass recovered (confident+prov)  : {n_conf + n_prov}")
        print(f"  = gamma-delta after paired pass           : {n_gd_floor + n_conf + n_prov}")

        # export per-cell detail
        det = pd.DataFrame({
            'barcode':   sub.obs_names,
            'tier1':     sub.obs[COL_CELLTYPE].astype(str).values,
            'cluster':   sub.obs[COL_CLUSTERS].astype(str).values if COL_CLUSTERS in adata.obs else '',
            'cd3_pos':   cd3_pos[idx],
            'b_flag':    b_flag[idx],
            'myeloid_flag': mye_flag[idx],
            'ery_plt_flag': ery_flag[idx],
            'animal':    sub.obs[ANIMAL_COL].astype(str).values if ANIMAL_COL in adata.obs else '',
            'tissue':    sub.obs[TISSUE_COL].astype(str).values if TISSUE_COL in adata.obs else '',
            'timepoint': sub.obs[TIME_COL].astype(str).values if TIME_COL in adata.obs else '',
        })
        det['call'] = np.where(
            det['b_flag'] | det['myeloid_flag'], 'quarantine',
            np.where(det['cd3_pos'], 'flip_confident', 'flip_provisional'))
        save_csv(det.set_index('barcode'), 'module1_paired_nonT_cells.csv')

    # ambient TCR on non-T lineages (the looser orphan-contig issue)
    if COL_HASTCR in adata.obs.columns:
        print("\n  ambient has_tcr rate on non-T lineages (orphan contigs):")
        ct = adata.obs[COL_CELLTYPE].astype(str)
        for lin in ['B cells', 'Monocytes', 'Macrophages']:
            m = ct == lin
            if m.sum() > 0:
                r = adata.obs.loc[m, COL_HASTCR].astype(bool).mean()
                print(f"    {lin:<14}: {100*r:.1f}%  (n={int(m.sum()):,})")


# %% ============================================================
# CELL 4 — MODULE 2 (clonotype feasibility: 40707 CD8 story)
# ============================================================
hr("CELL 4 — MODULE 2 (clonotype-definable T cells; 40707 CD8 vs CD4)")

if COL_PAIRING not in adata.obs.columns:
    print("  WARN: chain_pairing missing; skipping Module 2.")
else:
    pairing = adata.obs[COL_PAIRING].astype(str)
    print("  chain_pairing values present:")
    print(pairing.value_counts(dropna=False).to_string())

    present_def = [p for p in DEFINABLE_PAIRINGS if p in set(pairing)]
    missing_def = [p for p in DEFINABLE_PAIRINGS if p not in set(pairing)]
    if missing_def:
        print(f"\n  NOTE: definable categories not seen (ignored): {missing_def}")
    unclassified = sorted(set(pairing) - set(DEFINABLE_PAIRINGS)
                          - {'no TCR', 'nan', 'None'})
    if unclassified:
        print(f"  NOTE: pairings present but NOT counted as definable: {unclassified}")

    definable = pairing.isin(present_def).values
    is_t = (adata.obs[COL_CELLTYPE].astype(str) == 'T cells').values
    def_t = definable & is_t
    print(f"\n  clonotype-definable cells (complete pair): {int(definable.sum()):,}")
    print(f"  ...that are T cells                       : {int(def_t.sum()):,}")

    # CD4/CD8 split among definable T cells — RNA (raw counts) and ADT side by side
    layer = COUNTS_LAYER if have_counts else None
    cd8_rna = set_sum(CD8_GENES, layer=layer) >= POS_MIN_COUNTS
    cd4_rna = set_sum(CD4_GENES, layer=layer) >= POS_MIN_COUNTS

    def rna_call(i):
        if cd8_rna[i] and not cd4_rna[i]:
            return 'CD8'
        if cd4_rna[i] and not cd8_rna[i]:
            return 'CD4'
        if cd4_rna[i] and cd8_rna[i]:
            return 'DP'
        return 'DN'

    # ADT-based CD8 call: RELATIVE (CD8 CLR > CD4 CLR per cell). Above-zero is
    # useless here: the CD8 tag reads positive in ~97% of cells (broad
    # background), so a >0 test just re-counts def_t. The honest surface call is
    # CD8 vs CD4 head-to-head, which is also what the real annotation split uses.
    cd8_adt = cd4_adt = None
    if 'protein_clr' in adata.obsm and adt_names:
        clr = np.asarray(adata.obsm['protein_clr'], dtype=float)
        cd8_tag = resolve_adt('CD8', adt_names)
        cd4_tag = resolve_adt('CD4', adt_names)
        print(f"  ADT resolved: CD8 -> {cd8_tag}, CD4 -> {cd4_tag}")
        if cd8_tag and cd4_tag:
            cd8_clr = clr[:, adt_names.index(cd8_tag)]
            cd4_clr = clr[:, adt_names.index(cd4_tag)]
            cd8_adt = cd8_clr > cd4_clr
            cd4_adt = cd4_clr > cd8_clr
        else:
            print("  WARN: need both CD8 and CD4 tags for the relative ADT call.")

    # per-animal x timepoint definable-T and definable-CD8-T
    if ANIMAL_COL in adata.obs.columns and TIME_COL in adata.obs.columns:
        rna_lab = np.array([rna_call(i) for i in range(adata.n_obs)])
        df = pd.DataFrame({
            'animal':    adata.obs[ANIMAL_COL].astype(str).values,
            'timepoint': adata.obs[TIME_COL].astype(str).values,
            'tissue':    adata.obs[TISSUE_COL].astype(str).values if TISSUE_COL in adata.obs else '',
            'def_t':     def_t,
            'cd8_rna':   def_t & (rna_lab == 'CD8'),
            'cd4_rna':   def_t & (rna_lab == 'CD4'),
            'cd8_adt':   def_t & (cd8_adt if cd8_adt is not None else np.zeros(adata.n_obs, bool)),
        })

        summ = (df.groupby(['animal', 'timepoint'])
                  [['def_t', 'cd8_rna', 'cd4_rna', 'cd8_adt']].sum()
                  .astype(int))
        print("\n  definable T cells per animal x timepoint "
              "(def_t | CD8_rna | CD4_rna | CD8_adt):")
        print(summ.to_string())
        save_csv(summ, 'module2_definable_T_per_animal_timepoint.csv')

        print(f"\n  --- elite controller {ELITE}, by timepoint x tissue ---")
        ec = df[df['animal'] == ELITE]
        if len(ec):
            ectab = (ec.groupby(['timepoint', 'tissue'])
                       [['def_t', 'cd8_rna', 'cd4_rna', 'cd8_adt']].sum().astype(int))
            print(ectab.to_string())
            save_csv(ectab, 'module2_elite_40707_definable.csv')
            n_ec_cd8 = int(ec['cd8_rna'].sum())
            print(f"\n  {ELITE} definable CD8 T (RNA) total: {n_ec_cd8}")
            print("  Read: enough definable CD8 T to look for clonal expansion, or")
            print("        too few to power it? This is a feasibility count, not a")
            print("        clonotype result (define_clonotypes runs on the .h5mu next).")
        else:
            print(f"  WARN: no cells for animal {ELITE} found.")
    else:
        print("  WARN: animal_id/timepoint missing; per-group breakdown skipped.")


# %% ============================================================
# CELL 5 — MODULE 3 (plotting inputs: ADT + RNA overlay markers)
# ============================================================
hr("CELL 5 — MODULE 3 (plotting-readiness: ADT tags + RNA markers)")

if 'protein_clr' in adata.obsm and adt_names:
    clr = np.asarray(adata.obsm['protein_clr'], dtype=float)
    rows = []
    for j, tag in enumerate(adt_names):
        v = clr[:, j]
        rows.append({
            'tag': tag,
            'median_clr': np.median(v),
            'iqr_clr': np.percentile(v, 75) - np.percentile(v, 25),
            'frac_gt0': float((v > 0).mean()),
            'max_clr': v.max(),
        })
    adt_df = pd.DataFrame(rows).set_index('tag').round(3)
    # blank = essentially no positive cells; CLR never rises above the geo-mean.
    # (IQR is a poor criterion for CLR: a blank tag still has spread around a
    #  strongly negative center, which is why the old IQR test missed CD163/CD169.)
    adt_df['likely_blank'] = (adt_df['frac_gt0'] < 0.01) | (adt_df['max_clr'] <= 0)
    print("  ADT CLR summary (blank tags expected: CD163, CD169):")
    print(adt_df.to_string())
    save_csv(adt_df, 'module3_adt_clr_summary.csv')
    blanks = adt_df.index[adt_df['likely_blank']].tolist()
    print(f"\n  flagged likely-blank tags: {blanks}")
    print(f"  usable tags for plotting  : "
          f"{[t for t in adt_names if t not in blanks]}")
else:
    print("  WARN: protein_clr / adt_names missing; ADT check skipped.")

print("\n  --- RIRA RNA overlay markers present in var_names? ---")
overlay_rows = []
for lineage, genes in RIRA_PANEL.items():
    found = present(genes)
    missing = [g for g in genes if g not in found]
    print(f"    {lineage:<8}: {len(found)}/{len(genes)} present"
          + (f"  MISSING {missing}" if missing else ""))
    for g in genes:
        overlay_rows.append({'lineage': lineage, 'gene': g,
                             'present': g in adata.var_names})
save_csv(pd.DataFrame(overlay_rows).set_index('gene'),
         'module3_rira_marker_presence.csv')


# %% ============================================================
# CELL 6 — MODULE 4 (low-purity clusters: Tier 2 scoping)
# ============================================================
hr("CELL 6 — MODULE 4 (low-purity clusters)")

if COL_CLUSTERS not in adata.obs.columns:
    print("  WARN: clusters missing; skipping Module 4.")
else:
    ct = pd.crosstab(adata.obs[COL_CLUSTERS], adata.obs[COL_CELLTYPE])
    size = ct.sum(axis=1)
    dom = ct.idxmax(axis=1)
    purity = (ct.max(axis=1) / size).round(3)
    # second-most label
    def second(row):
        s = row.sort_values(ascending=False)
        return s.index[1] if len(s) > 1 and s.iloc[1] > 0 else ''
    second_lab = ct.apply(second, axis=1)

    clus_df = pd.DataFrame({
        'n_cells': size, 'dominant': dom, 'purity': purity,
        'second': second_lab,
    }).sort_values('purity')
    print("  per-cluster purity (sorted ascending):")
    print(clus_df.to_string())
    save_csv(clus_df, 'module4_cluster_purity.csv')

    low = clus_df.index[clus_df['purity'] < PURITY_CUTOFF].tolist()
    print(f"\n  clusters below purity {PURITY_CUTOFF}: {low}")

    if low:
        # RIRA panel mean expression (log-norm X) per low-purity cluster, by lineage
        print("\n  --- mean RIRA-panel expression per low-purity cluster (log-norm X) ---")
        clusters = adata.obs[COL_CLUSTERS].astype(str)
        panel_rows = []
        for lineage, genes in RIRA_PANEL.items():
            g = present(genes)
            if not g:
                continue
            M = adata[:, g].X
            M = M.toarray() if sparse.issparse(M) else np.asarray(M)
            score = M.mean(axis=1)   # simple mean over present panel genes
            for cl in low:
                m = (clusters == str(cl)).values
                panel_rows.append({'cluster': str(cl), 'lineage': lineage,
                                   'mean_expr': round(float(score[m].mean()), 3)})
        panel = pd.DataFrame(panel_rows).pivot(index='cluster',
                                               columns='lineage', values='mean_expr')
        print(panel.to_string())
        save_csv(panel, 'module4_lowpurity_rira_scores.csv')

        # ADT profile per low-purity cluster (resolve short -> full tag; keep
        # the short label as the display column)
        if 'protein_clr' in adata.obsm and adt_names:
            tag_pairs = [(t, resolve_adt(t, adt_names)) for t in ADT_PROFILE_TAGS]
            tag_pairs = [(s, f) for s, f in tag_pairs if f is not None]
            clr = np.asarray(adata.obsm['protein_clr'], dtype=float)
            adt_rows = []
            for cl in low:
                m = (clusters == str(cl)).values
                row = {'cluster': str(cl)}
                for short, full in tag_pairs:
                    row[short] = round(float(clr[m, adt_names.index(full)].mean()), 2)
                adt_rows.append(row)
            adt_prof = pd.DataFrame(adt_rows).set_index('cluster')
            print("\n  --- mean ADT CLR per low-purity cluster ---")
            print(adt_prof.to_string())
            save_csv(adt_prof, 'module4_lowpurity_adt_profile.csv')

        # data-driven top markers for the low-purity clusters
        if RUN_RANK_GENES:
            try:
                sc.tl.rank_genes_groups(adata, groupby=COL_CLUSTERS, groups=[str(c) for c in low],
                                        method='wilcoxon', n_genes=N_TOP_MARKERS,
                                        use_raw=False)
                print(f"\n  --- top {N_TOP_MARKERS} markers per low-purity cluster (wilcoxon) ---")
                names = adata.uns['rank_genes_groups']['names']
                topdf = pd.DataFrame({str(c): [names[str(c)][i] for i in range(N_TOP_MARKERS)]
                                      for c in low})
                print(topdf.to_string(index=False))
                save_csv(topdf, 'module4_lowpurity_top_markers.csv')
            except Exception as e:
                print(f"  WARN rank_genes_groups skipped: {e}")


# %% ============================================================
# CELL 7 — DECISION SUMMARY
# ============================================================
hr("CELL 7 — DECISION SUMMARY (what to read, what to decide)")
print("""
  MODULE 1 (paired pass): read module1_paired_nonT_cells.csv 'call' column.
    -> quarantine the 'quarantine' rows (B/myeloid doublets); flip the
       'flip_confident' rows to T; flip 'flip_provisional' to T but tag them
       provisional gamma-delta (real paired receptor, CD3 dropout). ery_plt_flag
       is informational, not a disqualifier.

  MODULE 2 (clonotype): read module2_elite_40707_definable.csv.
    -> if 40707 has a workable pool of definable CD8 T at 21DPI / necropsy, the
       CD8-expansion story is powered and define_clonotypes on the .h5mu is
       worth building. If the count is tiny, we learn that here, not later.

  MODULE 3 (plotting): module3_adt_clr_summary.csv should flag CD163/CD169 as
    likely_blank and leave 9 usable tags. Any missing RIRA markers are noted in
    module3_rira_marker_presence.csv for the overlay panel.

  MODULE 4 (Tier 2): module4_cluster_purity.csv lists the low-purity clusters;
    the RIRA-score and ADT-profile tables suggest whether each is a doublet,
    a transitional state, or a real subtype to split.

  Nothing above was written back to the object. Read the tables, then we decide.
""")
print("=" * 80)
print("DIAGNOSTIC COMPLETE — read-only. Object unchanged.")
print("=" * 80)
