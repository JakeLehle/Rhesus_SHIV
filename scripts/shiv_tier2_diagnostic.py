#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Rhesus SHIV — Tier 2 Reclustering + Signature Diagnostic (per lineage)
======================================================================
Breaks one Tier-1 lineage out of the canonical object, reclusters it on its own
variable genes with library batch correction, scores functional subtype
signatures, and reports subcluster-level score profiles, provisional subtypes,
data-driven markers, and animal / timepoint / tissue / infection enrichment.

This is the PREP diagnostic for Tier 2 annotation and DEG, not the annotation
itself. It does NOT bless final subtype names onto cells and does NOT run DEG.
We read the subcluster-by-signature table together, decide the annotation
mapping, then a second pass writes labels and runs DEG.

Built from the session notes:
  - Zoey's GZMB concern: the cytotoxic block (GZMB/GNLY/PRF1/NKG7) is scored as
    its OWN signature and we confirm the cytotoxic subclusters are CD3+ CD8 T,
    i.e. the signal is CD8 biology, not misassigned NK.
  - CXCR5 / reservoir hypothesis: a Tfh panel (CXCR5/PDCD1/BCL6/ICOS) is scored
    so we can see whether 40707 diverges there.
  - Transitional-vs-memory B caution: transitional is scored SEPARATELY from
    memory so we do not call something memory that makes no antibody.
  - "Unassigned that should be T": when LINEAGE == 'T cells', profile the
    Unassigned pool (CD3 / TCR status) to see how many are cryptic T cells.
  - Animal enrichment per subcluster for ALL animals, plus timepoint / tissue /
    infection, to find animal-skewed subclusters for DEG.

Read-only w.r.t. the canonical object. Writes a NEW subset object
(shiv_<lineage>_subset_reclustered.h5ad) + CSVs + figures. Spyder cells (# %%).
Runs headless under SLURM (sc_pre env).

Author: Jake Lehle
Date: July 2026
"""

# %% ============================================================
# CELL 1 — CONFIG AND DIALS
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
ADATA_IN    = os.path.join(WORKING_DIR, 'shiv_starsolo', 'merged',
                           'shiv_host_final_embedded.h5ad')

# ---- which lineage to break out (run T cells first) ----
LINEAGE = 'T cells'         # 'T cells' | 'B cells' | 'Monocytes'
TAG = {'T cells': 'Tcell', 'B cells': 'Bcell', 'Monocytes': 'Mono'}.get(LINEAGE, 'lineage')

OUT_DIR   = os.path.join(WORKING_DIR, 'annotation_output', f'tier2_{TAG}')
ADATA_OUT = os.path.join(WORKING_DIR, 'shiv_starsolo', 'merged',
                         f'shiv_{TAG}_subset_reclustered.h5ad')

# ---- recluster dials ----
N_HVG        = 2000
N_PCS        = 50
N_NEIGH_PCS  = 30
SUB_RES      = 1.0          # sub-Leiden resolution (0.6 for coarser)
BATCH_KEY    = 'library'    # BBKNN on GEM-well technical variation (animal is biology)
RANDOM_STATE = 0
COUNTS_LAYER = 'counts'
DO_PREBBKNN  = True         # compute pre-correction UMAP (raw animal separation)

# ---- columns ----
CELLTYPE   = 'tier1_celltype'
ANIMAL_COL = 'animal_id'
TIME_COL   = 'timepoint'
TISSUE_COL = 'tissue'
ELITE      = '40707'
INFECTED_ANIMALS = ['39272', '40702', '40707', '41861', '41862']
CONTROL_ANIMALS  = ['34315', '41903']
EXCLUDE_ANIMALS  = ['Unknown']

# ---- functional signature panels (present-checked; warn-and-continue) ----
PANELS = {
    'T cells': {
        'Naive':         ['CCR7', 'SELL', 'TCF7', 'LEF1', 'IL7R'],
        'CentralMem':    ['CCR7', 'SELL', 'IL7R', 'CD28', 'GPR183'],
        'EffectorMem':   ['GZMK', 'CCL5', 'CST7', 'KLRG1'],
        'Cytotoxic':     ['GZMB', 'GNLY', 'PRF1', 'NKG7', 'GZMH', 'FGFBP2'],  # the GZMB block
        'Th1':           ['TBX21', 'IFNG', 'CXCR3', 'STAT1'],
        'Th17':          ['RORC', 'IL17A', 'CCR6', 'IL23R', 'KLRB1'],
        'Tfh':           ['CXCR5', 'PDCD1', 'BCL6', 'ICOS', 'IL21'],
        'Treg':          ['FOXP3', 'IL2RA', 'CTLA4', 'IKZF2'],
        'Proliferating': ['MKI67', 'TOP2A'],
    },
    'B cells': {
        'Naive':         ['TCL1A', 'IGHD', 'FCER2', 'IL4R'],
        'Memory':        ['CD27', 'TNFRSF13B', 'AIM2', 'TXNIP'],
        'Transitional':  ['MME', 'VPREB3', 'SOX4', 'CD24'],
        'Plasma':        ['MZB1', 'XBP1', 'PRDM1', 'JCHAIN', 'SDC1'],
        'GerminalCenter':['BCL6', 'AICDA', 'RGS13', 'MEF2B'],
        'Activated':     ['CD69', 'CD83', 'NR4A1'],
    },
    'Monocytes': {
        'Classical':     ['CD14', 'S100A8', 'S100A9', 'VCAN'],
        'NonClassical':  ['FCGR3', 'CDKN1C', 'LILRB2', 'MS4A7'],
        'Intermediate':  ['CD74', 'FCGR3', 'HLA-DRA'],
        'Macrophage':    ['CD68', 'CD163', 'MRC1', 'MERTK', 'APOC1', 'C1QA'],
        'DC':            ['FLT3', 'CLEC9A', 'CLEC10A', 'IL3RA'],
    },
}
PANEL = PANELS.get(LINEAGE, {})

# T-cell CD4/CD8 split markers
CD8_GENES = ['CD8A', 'CD8B']
CD4_GENES = ['CD4']
CD3_GENES = ['CD3D', 'CD3E', 'CD3G']
POS_MIN   = 1

N_TOP_MARKERS = 15
ENRICH_CLIP   = 2.0

plt.rcParams.update({
    'font.size': 28, 'axes.titlesize': 34, 'axes.labelsize': 30,
    'xtick.labelsize': 24, 'ytick.labelsize': 24, 'legend.fontsize': 18,
    'figure.dpi': 100, 'pdf.fonttype': 42, 'ps.fonttype': 42,
})
os.makedirs(OUT_DIR, exist_ok=True)
sc.settings.figdir = OUT_DIR
pd.set_option('display.width', 200)
pd.set_option('display.max_columns', 60)
pd.set_option('display.max_rows', 200)


def savefig(fig, name):
    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(OUT_DIR, f'{name}.{ext}'),
                    bbox_inches='tight', dpi=300)
    print(f"    saved: {name}.pdf / .png")


def save_csv(df, name):
    df.to_csv(os.path.join(OUT_DIR, name))
    print(f"    saved: {name}")


def log2_enrichment(cluster_series, group_series):
    """log2( P(group|cluster) / P(group) ), NaN where group absent from cluster."""
    ct = pd.crosstab(cluster_series, group_series)
    p_obs = ct.div(ct.sum(axis=1), axis=0)
    p_base = ct.sum(axis=0) / ct.values.sum()
    e = np.log2(p_obs.div(p_base, axis=1)).replace([np.inf, -np.inf], np.nan)
    return ct, e


# %% ============================================================
# CELL 2 — LOAD + SUBSET TO LINEAGE
# ============================================================
print("=" * 80)
print(f"CELL 2 — LOAD + SUBSET TO '{LINEAGE}'")
print("=" * 80)

adata = sc.read_h5ad(ADATA_IN)
print(f"  full object: {adata.n_obs:,} cells")
adt_names = list(adata.uns.get('adt_names', []))

sub = adata[adata.obs[CELLTYPE].astype(str) == LINEAGE].copy()
print(f"  {LINEAGE}: {sub.n_obs:,} cells")


def present(a, genes):
    return [g for g in genes if g in a.var_names]


def gset_sum(a, genes, layer=None):
    genes = present(a, genes)
    if not genes:
        return np.zeros(a.n_obs)
    m = a[:, genes].layers[layer] if layer else a[:, genes].X
    m = m.toarray() if sparse.issparse(m) else np.asarray(m)
    return m.sum(axis=1)


have_counts = COUNTS_LAYER in sub.layers
print(f"  counts layer: {have_counts}; X.max()={float(sub.X.max()):.2f} "
      f"({'log-norm' if sub.X.max() < 50 else 'raw'})")
print(f"  animals: {sub.obs[ANIMAL_COL].astype(str).value_counts().to_dict()}")


# %% ============================================================
# CELL 3 — RECLUSTER (HVG, PCA, pre-BBKNN UMAP, BBKNN, sub-Leiden)
# ============================================================
print("\n" + "=" * 80)
print("CELL 3 — RECLUSTER")
print("=" * 80)

sc.pp.highly_variable_genes(sub, n_top_genes=N_HVG, flavor='seurat')
print(f"  HVG (seurat, within-lineage): {int(sub.var['highly_variable'].sum())}")
sc.tl.pca(sub, n_comps=N_PCS, use_highly_variable=True)

if DO_PREBBKNN:
    sc.pp.neighbors(sub, n_pcs=N_NEIGH_PCS)
    sc.tl.umap(sub)
    sub.obsm['X_umap_pre'] = sub.obsm['X_umap'].copy()
    print("  pre-BBKNN UMAP computed (raw animal separation)")

try:
    import bbknn
    bbknn.bbknn(sub, batch_key=BATCH_KEY, n_pcs=N_NEIGH_PCS)
    print(f"  BBKNN on '{BATCH_KEY}' ({sub.obs[BATCH_KEY].nunique()} batches)")
except ImportError:
    print("  WARN bbknn missing; standard neighbors (no correction)")
    sc.pp.neighbors(sub, n_pcs=N_NEIGH_PCS)

sc.tl.umap(sub)
sc.tl.leiden(sub, key_added='sub_leiden', resolution=SUB_RES,
             random_state=RANDOM_STATE)
n_sub = sub.obs['sub_leiden'].nunique()
print(f"  sub-Leiden (res {SUB_RES}): {n_sub} subclusters")
print(sub.obs['sub_leiden'].value_counts().sort_index().to_string())


# %% ============================================================
# CELL 4 — SCORE FUNCTIONAL SIGNATURES
# ============================================================
print("\n" + "=" * 80)
print("CELL 4 — SIGNATURE SCORES")
print("=" * 80)

score_cols = []
for name, genes in PANEL.items():
    g = present(sub, genes)
    miss = [x for x in genes if x not in g]
    if len(g) >= 1:
        sc.tl.score_genes(sub, gene_list=g, score_name=f'sig_{name}',
                          random_state=RANDOM_STATE)
        score_cols.append(f'sig_{name}')
        print(f"  {name:<14}: scored on {len(g)}/{len(genes)}"
              + (f"  MISSING {miss}" if miss else ""))
    else:
        print(f"  {name:<14}: WARN no markers present, skipped")


# %% ============================================================
# CELL 5 — SUBCLUSTER x SIGNATURE + PROVISIONAL SUBTYPE + TOP MARKERS
# ============================================================
print("\n" + "=" * 80)
print("CELL 5 — SUBCLUSTER PROFILES")
print("=" * 80)

# mean score per subcluster
mean_scores = sub.obs.groupby('sub_leiden')[score_cols].mean()
mean_scores.columns = [c.replace('sig_', '') for c in mean_scores.columns]
print("  mean signature score per subcluster:")
print(mean_scores.round(3).to_string())
save_csv(mean_scores.round(4), f'tier2_{TAG}_subcluster_signature_means.csv')

# provisional subtype: z-score each signature across cells, argmax per cell, mode per subcluster
Z = sub.obs[score_cols].copy()
Z = (Z - Z.mean()) / (Z.std() + 1e-9)
sub.obs['provisional_subtype'] = (Z.idxmax(axis=1)
                                  .str.replace('sig_', '', regex=False).values)
prov = (sub.obs.groupby('sub_leiden')['provisional_subtype']
        .agg(lambda s: s.value_counts().index[0]))
prov_frac = (sub.obs.groupby('sub_leiden')['provisional_subtype']
             .agg(lambda s: round(s.value_counts(normalize=True).iloc[0], 2)))
prov_tbl = pd.DataFrame({'provisional_subtype': prov, 'dominant_frac': prov_frac,
                         'n_cells': sub.obs['sub_leiden'].value_counts()})
print("\n  provisional subtype per subcluster (argmax of z-scored signatures):")
print(prov_tbl.to_string())
save_csv(prov_tbl, f'tier2_{TAG}_provisional_subtypes.csv')

# data-driven top markers per subcluster
try:
    sc.tl.rank_genes_groups(sub, groupby='sub_leiden', method='wilcoxon',
                            n_genes=N_TOP_MARKERS, use_raw=False)
    names = sub.uns['rank_genes_groups']['names']
    topdf = pd.DataFrame({cl: [names[cl][i] for i in range(N_TOP_MARKERS)]
                          for cl in sub.obs['sub_leiden'].cat.categories})
    print(f"\n  top {N_TOP_MARKERS} markers per subcluster:")
    print(topdf.to_string(index=False))
    save_csv(topdf, f'tier2_{TAG}_top_markers.csv')
except Exception as e:
    print(f"  WARN rank_genes_groups skipped: {e}")


# %% ============================================================
# CELL 6 — CD4/CD8 SPLIT (T cells only): RNA + relative ADT
# ============================================================
if LINEAGE == 'T cells':
    print("\n" + "=" * 80)
    print("CELL 6 — CD4/CD8 SPLIT (RNA + relative ADT)")
    print("=" * 80)
    layer = COUNTS_LAYER if have_counts else None
    cd8_rna = gset_sum(sub, CD8_GENES, layer) >= POS_MIN
    cd4_rna = gset_sum(sub, CD4_GENES, layer) >= POS_MIN

    def rna_call(i):
        if cd8_rna[i] and not cd4_rna[i]: return 'CD8'
        if cd4_rna[i] and not cd8_rna[i]: return 'CD4'
        if cd4_rna[i] and cd8_rna[i]:     return 'DP'
        return 'DN'
    sub.obs['cd48_rna'] = [rna_call(i) for i in range(sub.n_obs)]

    # relative ADT: CD8 CLR > CD4 CLR
    if 'protein_clr' in sub.obsm and adt_names:
        clr = np.asarray(sub.obsm['protein_clr'], dtype=float)
        def find(tag):
            if tag in adt_names: return adt_names.index(tag)
            for i, n in enumerate(adt_names):
                if n.split('_')[0] == tag or n.startswith(tag): return i
            return None
        i8, i4 = find('CD8'), find('CD4')
        if i8 is not None and i4 is not None:
            sub.obs['cd48_adt'] = np.where(clr[:, i8] > clr[:, i4], 'CD8', 'CD4')
        else:
            print("  WARN CD8/CD4 ADT tags not resolved")
    tab = pd.crosstab(sub.obs['sub_leiden'], sub.obs['cd48_rna'])
    if 'cd48_adt' in sub.obs:
        tab['adt_CD8'] = sub.obs.groupby('sub_leiden')['cd48_adt'].apply(
            lambda s: int((s == 'CD8').sum()))
    print("  CD4/CD8 per subcluster (RNA counts; adt_CD8 = relative-ADT CD8):")
    print(tab.to_string())
    save_csv(tab, f'tier2_{TAG}_cd4cd8_split.csv')


# %% ============================================================
# CELL 7 — ENRICHMENT: animal / timepoint / tissue / infection
# ============================================================
print("\n" + "=" * 80)
print("CELL 7 — SUBCLUSTER ENRICHMENT")
print("=" * 80)

sub.obs['_animal'] = sub.obs[ANIMAL_COL].astype(str)
sub.obs['infection'] = np.where(sub.obs['_animal'].isin(INFECTED_ANIMALS), 'infected',
                        np.where(sub.obs['_animal'].isin(CONTROL_ANIMALS), 'control', 'other'))

keep = ~sub.obs['_animal'].isin(EXCLUDE_ANIMALS)
_, e_animal = log2_enrichment(sub.obs.loc[keep, 'sub_leiden'],
                              sub.obs.loc[keep, '_animal'])
order = [a for a in INFECTED_ANIMALS if a in e_animal.columns] + \
        [a for a in CONTROL_ANIMALS if a in e_animal.columns]
order += [a for a in e_animal.columns if a not in order]
e_animal = e_animal[order]
print("  animal enrichment per subcluster (log2 obs/exp; ALL animals):")
print(e_animal.round(2).to_string())
save_csv(e_animal.round(3), f'tier2_{TAG}_animal_enrichment.csv')

# flag animal-skewed subclusters (any animal > 1 log2 = >2x its baseline)
print("\n  animal-skewed subclusters (max |log2| >= 1.0):")
for cl in e_animal.index:
    row = e_animal.loc[cl].dropna()
    if len(row) and row.abs().max() >= 1.0:
        top = row.idxmax()
        flag = ' <ELITE>' if top == ELITE else ''
        print(f"    subcluster {cl}: {top} at {row[top]:.2f}{flag}  "
              f"(provisional {prov.get(cl, '?')})")

for col, label in [(TIME_COL, 'timepoint'), (TISSUE_COL, 'tissue'),
                   ('infection', 'infection')]:
    if col in sub.obs.columns:
        _, e = log2_enrichment(sub.obs['sub_leiden'], sub.obs[col].astype(str))
        save_csv(e.round(3), f'tier2_{TAG}_{label}_enrichment.csv')
        print(f"  saved {label} enrichment.")


# %% ============================================================
# CELL 8 — LINEAGE-SPECIFIC CHECKS
# ============================================================
print("\n" + "=" * 80)
print("CELL 8 — LINEAGE-SPECIFIC CHECKS")
print("=" * 80)

if LINEAGE == 'T cells':
    # Zoey's GZMB concern: cytotoxic-high subclusters should be CD3+ CD8 T
    if 'sig_Cytotoxic' in score_cols:
        cyto_by = sub.obs.groupby('sub_leiden')['sig_Cytotoxic'].mean()
        hi = cyto_by[cyto_by > cyto_by.quantile(0.75)].index.tolist()
        layer = COUNTS_LAYER if have_counts else None
        cd3 = gset_sum(sub, CD3_GENES, layer) >= POS_MIN
        print("  GZMB/cytotoxic check (cytotoxic-high subclusters should be CD3+ CD8 T):")
        for cl in hi:
            m = (sub.obs['sub_leiden'] == cl).values
            cd3f = 100 * cd3[m].mean()
            cd8f = 100 * (np.array(sub.obs['cd48_rna'])[m] == 'CD8').mean() \
                if 'cd48_rna' in sub.obs else np.nan
            print(f"    subcluster {cl}: cytotoxic {cyto_by[cl]:.2f}, "
                  f"CD3+ {cd3f:.0f}%, CD8(RNA) {cd8f:.0f}%")
        print("  => high CD3+/CD8 here = the cytotoxic signal is CD8 T biology,")
        print("     supporting the Tier 1 decision to pull GZMB out of NK scoring.")

    # cryptic T among Unassigned (from the note: why did T cells fall through?)
    una = adata[adata.obs[CELLTYPE].astype(str) == 'Unassigned']
    if una.n_obs:
        layer_u = COUNTS_LAYER if COUNTS_LAYER in una.layers else None
        cd3u = gset_sum(una, CD3_GENES, layer_u) >= POS_MIN
        tcru = una.obs['has_tcr'].astype(bool).values if 'has_tcr' in una.obs else np.zeros(una.n_obs, bool)
        print(f"\n  Unassigned pool: {una.n_obs:,} cells")
        print(f"    CD3+ (cryptic T by RNA)     : {int(cd3u.sum()):,} ({100*cd3u.mean():.1f}%)")
        print(f"    has TCR (cryptic T by VDJ)  : {int(tcru.sum()):,} ({100*tcru.mean():.1f}%)")
        print(f"    CD3+ OR TCR (likely T)      : {int((cd3u | tcru).sum()):,} "
              f"({100*(cd3u | tcru).mean():.1f}%)")
        print("  => these are candidates for recovery to T in the annotation pass.")

elif LINEAGE == 'B cells':
    # transitional-vs-memory caution: memory-scored subclusters w/o plasma/secretory
    if 'sig_Memory' in score_cols:
        print("  B memory-vs-transitional caution (memory call needs surface/secretory support):")
        for cl in mean_scores.index:
            mem = mean_scores.loc[cl].get('Memory', np.nan)
            trans = mean_scores.loc[cl].get('Transitional', np.nan)
            plasma = mean_scores.loc[cl].get('Plasma', np.nan)
            print(f"    subcluster {cl}: Memory {mem:.2f}, Transitional {trans:.2f}, "
                  f"Plasma {plasma:.2f}, provisional {prov.get(cl, '?')}")
        print("  => a 'Memory'-provisional subcluster with high Transitional and near-zero")
        print("     Plasma is likely early/transitional, not antibody-producing memory.")
else:
    print(f"  (no extra lineage-specific check for {LINEAGE})")


# %% ============================================================
# CELL 9 — FIGURES
# ============================================================
print("\n" + "=" * 80)
print("CELL 9 — FIGURES")
print("=" * 80)

# UMAP panels (post-BBKNN): subclusters, animal, timepoint, provisional subtype
fig, axes = plt.subplots(2, 2, figsize=(24, 20))
for ax, key, ttl in zip(axes.ravel(),
                        ['sub_leiden', '_animal', TIME_COL, 'provisional_subtype'],
                        ['Sub-Leiden', 'Animal', 'Timepoint', 'Provisional subtype']):
    if key in sub.obs.columns:
        sc.pl.umap(sub, color=key, ax=ax, show=False, title=ttl, frameon=False,
                   size=10, legend_loc='on data' if key in ('sub_leiden',) else 'right margin',
                   legend_fontsize=14)
plt.tight_layout()
savefig(fig, f'tier2_{TAG}_umap_panels')
plt.close(fig)

# pre vs post BBKNN by animal (technical vs biological separation)
if DO_PREBBKNN and 'X_umap_pre' in sub.obsm:
    fig, axes = plt.subplots(1, 2, figsize=(24, 11))
    for ax, key, ttl in zip(axes, ['X_umap_pre', 'X_umap'],
                            ['Animal (pre-BBKNN)', 'Animal (post-BBKNN)']):
        sub.obsm['X_umap_tmp'] = sub.obsm[key]
        sc.pl.embedding(sub, basis='X_umap_tmp', color='_animal', ax=ax,
                        show=False, title=ttl, frameon=False, size=10,
                        legend_fontsize=16)
    del sub.obsm['X_umap_tmp']
    plt.tight_layout()
    savefig(fig, f'tier2_{TAG}_prepost_bbknn_animal')
    plt.close(fig)

# signature score heatmap (subcluster x signature, z across subclusters)
if score_cols:
    Zc = (mean_scores - mean_scores.mean()) / (mean_scores.std() + 1e-9)
    fig, ax = plt.subplots(figsize=(1.1 * Zc.shape[1] + 5, 0.55 * Zc.shape[0] + 4))
    im = ax.imshow(Zc.values, cmap='RdBu_r', vmin=-2, vmax=2, aspect='auto')
    ax.set_xticks(range(Zc.shape[1])); ax.set_xticklabels(Zc.columns, rotation=45, ha='right')
    ax.set_yticks(range(Zc.shape[0])); ax.set_yticklabels(Zc.index)
    ax.set_xlabel('signature'); ax.set_ylabel('subcluster')
    ax.set_title(f'{LINEAGE} subcluster signature profile (z across subclusters)')
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    plt.tight_layout()
    savefig(fig, f'tier2_{TAG}_signature_heatmap')
    plt.close(fig)

# animal enrichment heatmap
fig, ax = plt.subplots(figsize=(1.3 * e_animal.shape[1] + 5, 0.55 * e_animal.shape[0] + 4))
cmap = plt.cm.RdBu_r.copy(); cmap.set_bad('#E0E0E0')
im = ax.imshow(e_animal.values, cmap=cmap, vmin=-ENRICH_CLIP, vmax=ENRICH_CLIP, aspect='auto')
ax.set_xticks(range(e_animal.shape[1]))
ax.set_xticklabels([f"{a}*" if a == ELITE else a for a in e_animal.columns], rotation=45, ha='right')
ax.set_yticks(range(e_animal.shape[0])); ax.set_yticklabels(e_animal.index)
ax.set_xlabel('animal'); ax.set_ylabel('subcluster')
ax.set_title(f'{LINEAGE} animal enrichment (log2 obs/exp); * = elite controller')
plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
plt.tight_layout()
savefig(fig, f'tier2_{TAG}_animal_enrichment_heatmap')
plt.close(fig)


# %% ============================================================
# CELL 10 — WRITE SUBSET OBJECT + SUMMARY
# ============================================================
print("\n" + "=" * 80)
print("CELL 10 — WRITE SUBSET OBJECT")
print("=" * 80)
sub.write_h5ad(ADATA_OUT)
print(f"  [SAVED] {ADATA_OUT}")
print(f"  new obs: sub_leiden, provisional_subtype, "
      f"{'cd48_rna/cd48_adt, ' if LINEAGE == 'T cells' else ''}infection, sig_* scores")
print(f"  canonical object untouched at {os.path.basename(ADATA_IN)}")
print(f"""
  READ (to plan Tier 2 annotation + DEG next):
   - tier2_{TAG}_subcluster_signature_means.csv + the signature heatmap: which
     functional subtype each subcluster is, and whether Th1/Th17/Tfh separate.
   - tier2_{TAG}_provisional_subtypes.csv: argmax call per subcluster (a draft,
     not final; we map names together).
   - tier2_{TAG}_animal_enrichment.csv + heatmap: animal-skewed subclusters (ALL
     animals) = the DEG candidates. Watch 40707.
   - {'tier2_'+TAG+'_cd4cd8_split.csv + the GZMB/cytotoxic check above.' if LINEAGE=='T cells' else 'B/mono lineage-specific check above.'}

  NOT DONE here (next pass): final subtype labels written to cells, and the DEG
  on the animal-skewed / elite-controller subcluster. This diagnostic gathers
  what we need to decide those.
""")
print("=" * 80)
print(f"TIER 2 DIAGNOSTIC COMPLETE ({LINEAGE}) — canonical object untouched.")
print("=" * 80)
