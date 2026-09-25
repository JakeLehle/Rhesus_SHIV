#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Rhesus SHIV - Tier 2 T Cell Annotation Write (step 3a of 3)
===========================================================
Writes final Tier-2 T-cell subtype labels onto the v2 recovered T subset, adds
the per-cell attributes we decided to keep OUT of the labels (CD4/CD8 surface
call, proliferation gate, CD4 activation score), propagates the fine T labels
back onto the global object for a Tier-2 UMAP, and saves both to NEW files.

Label map (v2 sub_leiden -> subtype), locked after re-verifying every cluster
against the v2 diagnostic tables (signature means, ADT split, markers, animal /
tissue / timepoint enrichment, and the cluster-9 doublet adjudication):

  0, 5  -> CD8 effector memory      (CD3-high, 90-95% ADT-CD8, NKG7/GZMB/CCL5;
                                     5 carries GZMK; both PBMC-enriched)
  10    -> CD8 terminal effector    (cytotoxic 2.38, GNLY/PRF1/HOPX/TYROBP,
                                     ~100% ADT-CD8, LN-excluded, necropsy-high)
  3, 4  -> CD4 central memory       (77-82% CD4 by ADT; 4 activation-leaning,
                                     captured as a per-cell score, not a label)
  1,2,6,7 -> Naive/resting T        (7 clean naive LEF1/BACH2; 1/2/6 ribosomal-
                                     structured; CD4/CD8 kept as an ADT attribute)
  8     -> QC: high-mito            (KEG06 top-to-bottom, 100% Indian controls)
  9     -> QC: T-B doublet          (dblt 0.200 vs 0.05-0.11 field; 79% co-express
                                     CD79A/CD79B/MS4A1 with 89% CD3+)

Decisions carried as per-cell attributes rather than subtype names:
  - CD4/CD8 surface call: cd48_adt (CD8 CLR > CD4 CLR per cell). The naive pool
    stays ONE label carrying this attribute (to be split with Zoey later).
  - Proliferation: is_proliferating gate (MKI67/TOP2A) + the continuous
    sig_Proliferating score. The clean proliferating cluster dissolved on the
    recovered pool, so proliferation is per-cell, reported per animal in 3b.
  - CD4 activation: sig_Activation on immediate-early genes. NOTE this is partly
    a tissue-dissociation artifact, so it is a soft score, never a label.

QC-excluded cells are FLAGGED (qc_exclude), not deleted. Biological analysis in
3b filters qc_exclude == False. Both input objects are left untouched; output
goes to new files. Spyder cells (# %%). Headless-safe (sc_pre env).

This is step 3a. Step 3b (cell-level DEG) is written after we review this output.

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
SUBSET_IN   = os.path.join(MERGED, 'shiv_Tcell_subset_reclustered_v2.h5ad')
GLOBAL_IN   = os.path.join(MERGED, 'shiv_host_final_tcrTrecovered.h5ad')
SUBSET_OUT  = os.path.join(MERGED, 'shiv_Tcell_annotated.h5ad')
GLOBAL_OUT  = os.path.join(MERGED, 'shiv_host_tier2_annotated.h5ad')
OUT_DIR     = os.path.join(WORKING_DIR, 'annotation_output', 'tier2_annotation')

DRY_RUN = False

# ---- columns ----
CLUSTER_KEY  = 'sub_leiden'
TIER1_COL    = 'tier1_celltype'
SUBTYPE_COL  = 'tier2_subtype'      # NEW on subset
GLOBAL_T2COL = 'tier2_celltype'     # NEW on global (fine T + coarse others)
ANIMAL_COL   = 'animal_id'
TIME_COL     = 'timepoint'
TISSUE_COL   = 'tissue'
COUNTS_LAYER = 'counts'
UMAP_KEY     = 'X_umap'
ELITE        = '40707'

# ---- the locked label map (v2 numbering) ----
SUBTYPE_MAP = {
    '0':  'CD8 effector memory',
    '5':  'CD8 effector memory',
    '10': 'CD8 terminal effector',
    '3':  'CD4 central memory',
    '4':  'CD4 central memory',
    '1':  'Naive/resting T',
    '2':  'Naive/resting T',
    '6':  'Naive/resting T',
    '7':  'Naive/resting T',
    '8':  'QC: high-mito',
    '9':  'QC: T-B doublet',
}
QC_LABELS = {'QC: high-mito', 'QC: T-B doublet'}

# ---- per-cell attribute gene sets ----
PROLIF_GENES     = ['MKI67', 'TOP2A']                 # proliferation gate
ACTIVATION_GENES = ['FOS', 'FOSB', 'JUN', 'JUNB', 'JUND',
                    'EGR1', 'NR4A1', 'DUSP1']          # immediate-early (soft)
CD8_GENES = ['CD8A', 'CD8B']
CD4_GENES = ['CD4']
POS_MIN   = 1

# ADT tags to overlay on the T-subset UMAP as an annotation cross-check
ADT_OVERLAY = ['CD8', 'CD4', 'CX3CR1', 'CD28', 'CD69']

# ---- ordering / colors (all hex) ----
SUBTYPE_ORDER = ['Naive/resting T', 'CD4 central memory', 'CD8 effector memory',
                 'CD8 terminal effector', 'QC: high-mito', 'QC: T-B doublet']
TIER2_COLORS = {
    'Naive/resting T':       '#64B5F6',
    'CD4 central memory':    '#43A047',
    'CD8 effector memory':   '#FB8C00',
    'CD8 terminal effector': '#D32F2F',
    'QC: high-mito':         '#BDBDBD',
    'QC: T-B doublet':       '#9E9E9E',
}
# coarse non-T lineages for the global Tier-2 UMAP
GLOBAL_EXTRA = {
    'B cells': '#1976D2', 'NK cells': '#7B1FA2', 'Monocytes': '#F9A825',
    'Macrophages': '#00897B', 'DCs': '#6D4C41', 'Mast cells': '#455A64',
    'Erythrocytes': '#8D6E63', 'Platelets': '#C2185B', 'Unassigned': '#E0E0E0',
}
GLOBAL_COLORS = {**TIER2_COLORS, **GLOBAL_EXTRA}
GLOBAL_ORDER = (['Naive/resting T', 'CD4 central memory', 'CD8 effector memory',
                 'CD8 terminal effector']
                + ['B cells', 'NK cells', 'Monocytes', 'Macrophages', 'DCs',
                   'Mast cells', 'Erythrocytes', 'Platelets', 'Unassigned']
                + ['QC: high-mito', 'QC: T-B doublet'])
TIME_ORDER = ['Pre', '21 DPI', 'Necropsy', 'Non-Infected', 'Unknown']

plt.rcParams.update({
    'font.size': 28, 'axes.titlesize': 34, 'axes.labelsize': 30,
    'xtick.labelsize': 24, 'ytick.labelsize': 24, 'legend.fontsize': 18,
    'figure.dpi': 100, 'pdf.fonttype': 42, 'ps.fonttype': 42,
})
os.makedirs(OUT_DIR, exist_ok=True)
sc.settings.figdir = OUT_DIR
pd.set_option('display.width', 200)
pd.set_option('display.max_columns', 40)
pd.set_option('display.max_rows', 100)


def savefig(fig, name):
    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(OUT_DIR, f'{name}.{ext}'),
                    bbox_inches='tight', dpi=300)
    print(f"    saved: {name}.pdf / .png")


def save_csv(df, name, index=True):
    df.to_csv(os.path.join(OUT_DIR, name), index=index)
    print(f"    saved: {name}")


def present(a, genes):
    return [g for g in genes if g in a.var_names]


def gset_sum(a, genes, layer=None):
    g = present(a, genes)
    if not g:
        return np.zeros(a.n_obs)
    m = a[:, g].layers[layer] if layer else a[:, g].X
    m = m.toarray() if sparse.issparse(m) else np.asarray(m)
    return m.sum(axis=1)


def ordered_present(values, order):
    seen = set(pd.Series(values).dropna().astype(str).unique())
    return [c for c in order if c in seen]


def apply_colors(a, key, order, cmap):
    cats = ordered_present(a.obs[key], order)
    a.obs[key] = pd.Categorical(a.obs[key].astype(str), categories=cats, ordered=True)
    a.uns[f'{key}_colors'] = [cmap.get(c, '#000000') for c in cats]
    return cats


def resolve_adt(adt_names, tag):
    if tag in adt_names:
        return adt_names.index(tag)
    for i, n in enumerate(adt_names):
        if n.split('_')[0] == tag:
            return i
    for i, n in enumerate(adt_names):
        if n.startswith(tag):
            return i
    return None


# %% ============================================================
# CELL 2 - LOAD SUBSET + WRITE SUBTYPE LABELS
# ============================================================
print("=" * 80)
print("CELL 2 - LOAD SUBSET + WRITE SUBTYPE LABELS")
print("=" * 80)

sub = sc.read_h5ad(SUBSET_IN)
print(f"  subset: {sub.n_obs:,} cells x {sub.n_vars:,} genes")

cl = sub.obs[CLUSTER_KEY].astype(str)
unmapped = sorted(set(cl.unique()) - set(SUBTYPE_MAP))
if unmapped:
    raise ValueError(f"sub_leiden values with no label in SUBTYPE_MAP: {unmapped}")

sub.obs[SUBTYPE_COL] = cl.map(SUBTYPE_MAP).values
sub.obs['qc_exclude'] = sub.obs[SUBTYPE_COL].isin(QC_LABELS).values

print("  subtype label counts:")
print(sub.obs[SUBTYPE_COL].value_counts().to_string())
print(f"\n  qc_exclude (clusters 8 + 9): {int(sub.obs['qc_exclude'].sum()):,} cells")
print(f"  biological T cells (qc_exclude == False): "
      f"{int((~sub.obs['qc_exclude']).sum()):,}")


# %% ============================================================
# CELL 3 - PER-CELL ATTRIBUTES (surface call, prolif gate, activation score)
# ============================================================
print("\n" + "=" * 80)
print("CELL 3 - PER-CELL ATTRIBUTES")
print("=" * 80)
have_counts = COUNTS_LAYER in sub.layers
layer = COUNTS_LAYER if have_counts else None

# CD4/CD8 surface call (carry the diagnostic's cd48_adt; recompute if absent)
if 'cd48_adt' not in sub.obs.columns:
    adt_names = list(sub.uns.get('adt_names', []))
    if 'protein_clr' in sub.obsm and adt_names:
        clr = np.asarray(sub.obsm['protein_clr'], dtype=float)
        i8, i4 = resolve_adt(adt_names, 'CD8'), resolve_adt(adt_names, 'CD4')
        if i8 is not None and i4 is not None:
            sub.obs['cd48_adt'] = np.where(clr[:, i8] > clr[:, i4], 'CD8', 'CD4')
if 'cd48_adt' in sub.obs.columns:
    print("  surface CD4/CD8 call (cd48_adt, CD8 CLR > CD4 CLR):")
    print(pd.crosstab(sub.obs[SUBTYPE_COL], sub.obs['cd48_adt']).to_string())

# proliferation gate
prolif_hits = present(sub, PROLIF_GENES)
sub.obs['is_proliferating'] = gset_sum(sub, PROLIF_GENES, layer) >= POS_MIN
print(f"\n  proliferation gate on {prolif_hits}: "
      f"{int(sub.obs['is_proliferating'].sum()):,} cells "
      f"({100*sub.obs['is_proliferating'].mean():.1f}%)")

# CD4 activation score (immediate-early; soft, partly dissociation artifact)
act_hits = present(sub, ACTIVATION_GENES)
if act_hits:
    sc.tl.score_genes(sub, gene_list=act_hits, score_name='sig_Activation',
                      random_state=0)
    print(f"  sig_Activation scored on {len(act_hits)}/{len(ACTIVATION_GENES)}: "
          f"{act_hits}")
    miss = [g for g in ACTIVATION_GENES if g not in act_hits]
    if miss:
        print(f"    (missing, skipped: {miss})")
else:
    sub.obs['sig_Activation'] = 0.0
    print("  WARN no activation genes present; sig_Activation = 0")

# quick within-CD4-CM activation read (4 was the activation-leaning cluster)
cm = sub.obs[SUBTYPE_COL] == 'CD4 central memory'
if cm.any():
    print("\n  mean sig_Activation within CD4 CM by source sub_leiden:")
    print(sub.obs.loc[cm].groupby(CLUSTER_KEY, observed=True)['sig_Activation']
          .mean().round(3).to_string())


# %% ============================================================
# CELL 4 - LABEL SUMMARY + PER-ANIMAL COMPOSITION (read before propagate)
# ============================================================
print("\n" + "=" * 80)
print("CELL 4 - LABEL SUMMARY + PER-ANIMAL COMPOSITION")
print("=" * 80)

# subtype x key-signature sanity (do the labels hold on the scores?)
sig_cols = [c for c in sub.obs.columns if c.startswith('sig_')]
if sig_cols:
    sanity = sub.obs.groupby(SUBTYPE_COL, observed=True)[sig_cols].mean().round(3)
    sanity.columns = [c.replace('sig_', '') for c in sanity.columns]
    print("  mean signature score per subtype (label sanity):")
    print(sanity.to_string())
    save_csv(sanity, 'tier2_subtype_signature_sanity.csv')

# per-animal subtype composition (BIOLOGICAL cells only, drop Unknown animal)
bio = (~sub.obs['qc_exclude']) & (sub.obs[ANIMAL_COL].astype(str) != 'Unknown')
comp = pd.crosstab(sub.obs.loc[bio, ANIMAL_COL].astype(str),
                   sub.obs.loc[bio, SUBTYPE_COL])
comp_frac = comp.div(comp.sum(axis=1), axis=0).round(4)
print("\n  per-animal subtype counts (biological T only):")
print(comp.to_string())
print("\n  per-animal subtype fraction:")
print(comp_frac.to_string())
save_csv(comp, 'tier2_per_animal_subtype_counts.csv')
save_csv(comp_frac, 'tier2_per_animal_subtype_fraction.csv')

# flag where 40707 sits (sets up the 3b DEG framing, not a test)
if ELITE in comp_frac.index:
    print(f"\n  elite controller {ELITE} subtype fractions vs cohort mean:")
    cohort = comp_frac.drop(index=ELITE).mean()
    cmp = pd.DataFrame({'40707': comp_frac.loc[ELITE], 'cohort_mean': cohort.round(4),
                        'delta': (comp_frac.loc[ELITE] - cohort).round(4)})
    print(cmp.to_string())
    save_csv(cmp, 'tier2_elite_vs_cohort_composition.csv')


# %% ============================================================
# CELL 5 - PROPAGATE FINE T LABELS ONTO THE GLOBAL OBJECT
# ============================================================
print("\n" + "=" * 80)
print("CELL 5 - PROPAGATE TO GLOBAL")
print("=" * 80)

glob = sc.read_h5ad(GLOBAL_IN)
print(f"  global: {glob.n_obs:,} cells")

# start from tier1, overwrite T cells with the fine subtype by barcode
glob.obs[GLOBAL_T2COL] = glob.obs[TIER1_COL].astype(str).values
common = glob.obs_names.intersection(sub.obs_names)
print(f"  barcodes matched subset<->global: {len(common):,} "
      f"(subset has {sub.n_obs:,})")
if len(common) != sub.n_obs:
    print("  WARN not all subset barcodes found in global; check obs_names.")

glob.obs.loc[common, GLOBAL_T2COL] = sub.obs.loc[common, SUBTYPE_COL].values
# carry qc_exclude for T cells (False elsewhere)
glob.obs['qc_exclude'] = False
glob.obs.loc[common, 'qc_exclude'] = sub.obs.loc[common, 'qc_exclude'].values

# any generic 'T cells' left means a subset barcode failed to map
leftover_T = int((glob.obs[GLOBAL_T2COL] == 'T cells').sum())
print(f"  generic 'T cells' remaining on global (should be 0): {leftover_T:,}")
print("\n  global Tier-2 distribution:")
print(glob.obs[GLOBAL_T2COL].value_counts().to_string())


# %% ============================================================
# CELL 6 - FIGURES
# ============================================================
print("\n" + "=" * 80)
print("CELL 6 - FIGURES")
print("=" * 80)

# ---- Fig 1: global Tier-2 UMAP (fine T + coarse others) ----
if UMAP_KEY in glob.obsm:
    try:
        apply_colors(glob, GLOBAL_T2COL, GLOBAL_ORDER, GLOBAL_COLORS)
        fig, ax = plt.subplots(figsize=(16, 13))
        sc.pl.umap(glob, color=GLOBAL_T2COL, ax=ax, show=False, frameon=False,
                   size=10, title='Global Tier 2 (fine T + coarse lineages)',
                   legend_loc='right margin', legend_fontsize=16)
        plt.tight_layout()
        savefig(fig, 'tier2_global_umap')
        plt.close(fig)
    except Exception as e:
        print(f"  WARN Fig 1 skipped: {e}")

# ---- Fig 2: T-subset UMAP by subtype ----
if UMAP_KEY in sub.obsm:
    try:
        apply_colors(sub, SUBTYPE_COL, SUBTYPE_ORDER, TIER2_COLORS)
        fig, ax = plt.subplots(figsize=(15, 13))
        sc.pl.umap(sub, color=SUBTYPE_COL, ax=ax, show=False, frameon=False,
                   size=16, title='T cell subtypes (recovered pool)',
                   legend_loc='right margin', legend_fontsize=18)
        plt.tight_layout()
        savefig(fig, 'tier2_Tsubset_umap_subtype')
        plt.close(fig)
    except Exception as e:
        print(f"  WARN Fig 2 skipped: {e}")

# ---- Fig 3: ADT overlay on T-subset UMAP (annotation cross-check) ----
adt_names = list(sub.uns.get('adt_names', []))
if UMAP_KEY in sub.obsm and 'protein_clr' in sub.obsm and adt_names:
    try:
        clr = np.asarray(sub.obsm['protein_clr'], dtype=float)
        keys = []
        for tag in ADT_OVERLAY:
            i = resolve_adt(adt_names, tag)
            if i is not None:
                k = f'adt_{tag}'
                sub.obs[k] = clr[:, i]
                keys.append((k, adt_names[i]))
        n = len(keys)
        ncol = min(3, n); nrow = int(np.ceil(n / ncol))
        fig, axes = plt.subplots(nrow, ncol, figsize=(9 * ncol, 8 * nrow))
        axes = np.atleast_1d(axes).ravel()
        for ax, (k, full) in zip(axes, keys):
            sc.pl.umap(sub, color=k, ax=ax, show=False, frameon=False, size=10,
                       title=full.replace('_TotalSeqC', ''), cmap='viridis',
                       colorbar_loc='right')
        for ax in axes[n:]:
            ax.axis('off')
        plt.tight_layout()
        savefig(fig, 'tier2_Tsubset_umap_ADT_overlay')
        plt.close(fig)
        for k, _ in keys:
            del sub.obs[k]
    except Exception as e:
        print(f"  WARN Fig 3 skipped: {e}")

# ---- Fig 4: per-animal subtype composition (biological T, stacked) ----
try:
    cols = ordered_present(list(comp_frac.columns), SUBTYPE_ORDER)
    cf = comp_frac[cols]
    an_order = [a for a in ['39272', '40702', ELITE, '41861', '41862',
                            '34315', '41903'] if a in cf.index]
    an_order += [a for a in cf.index if a not in an_order]
    cf = cf.reindex(an_order)
    fig, ax = plt.subplots(figsize=(16, 10))
    bottom = np.zeros(len(cf))
    for c in cols:
        ax.bar(range(len(cf)), cf[c].values, bottom=bottom,
               color=TIER2_COLORS.get(c, '#000000'), label=c, width=0.8)
        bottom += cf[c].values
    labs = [f"{a}*" if a == ELITE else a for a in cf.index]
    ax.set_xticks(range(len(cf))); ax.set_xticklabels(labs, rotation=45, ha='right')
    ax.set_ylabel('fraction of T cells'); ax.set_ylim(0, 1)
    ax.set_title('T subtype composition per animal\n(* = elite controller 40707)')
    ax.legend(bbox_to_anchor=(1.01, 1), loc='upper left', fontsize=16)
    ax.spines[['top', 'right']].set_visible(False)
    plt.tight_layout()
    savefig(fig, 'tier2_per_animal_subtype_composition')
    plt.close(fig)
except Exception as e:
    print(f"  WARN Fig 4 skipped: {e}")


# %% ============================================================
# CELL 7 - WRITE BOTH OBJECTS (new files; inputs untouched)
# ============================================================
print("\n" + "=" * 80)
print("CELL 7 - WRITE")
print("=" * 80)

if DRY_RUN:
    print("  DRY_RUN = True -> nothing written.")
else:
    sub.write_h5ad(SUBSET_OUT)
    print(f"  [SAVED] {SUBSET_OUT}")
    print(f"    new obs: {SUBTYPE_COL}, qc_exclude, is_proliferating, "
          f"sig_Activation (cd48_adt carried)")
    glob.write_h5ad(GLOBAL_OUT)
    print(f"  [SAVED] {GLOBAL_OUT}")
    print(f"    new obs: {GLOBAL_T2COL} (fine T + coarse lineages), qc_exclude")
    print(f"  inputs untouched: {os.path.basename(SUBSET_IN)}, "
          f"{os.path.basename(GLOBAL_IN)}")
    print(f"""
  REVIEW before step 3b (DEG):
   - tier2_Tsubset_umap_ADT_overlay: CD8 protein on CD8 EM/terminal, CD4 on
     CD4 CM, CX3CR1 on terminal effector. If those track, the labels hold.
   - tier2_per_animal_subtype_fraction.csv + composition figure: sane splits,
     and where 40707 sits (tier2_elite_vs_cohort_composition.csv).
   - tier2_subtype_signature_sanity.csv: each subtype tops its own signature.

  THEN 3b runs cell-level DEG on qc_exclude == False:
   - 40707 vs other 4 infected within pooled CD8 EM (0+5), + per-animal panel
   - Pre vs necropsy within CD8 EM (per-animal caution: 39272 has no necropsy)
   - LN vs PBMC within the CD8 pool, + per-animal CXCR5 readout in LN cells
   - terminal-effector = compositional depletion (reported, not a DEG)
""")

print("=" * 80)
print("TIER 2 ANNOTATION WRITE COMPLETE (step 3a) - inputs untouched.")
print("=" * 80)
