#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Rhesus SHIV - Tier 2 Animal Enrichment + Tier 2 UMAPs (T cells AND myeloid)
============================================================================
Supersedes SHIV_tier2_myeloid_animal_enrichment.py. Two changes from that one:

  1. BOTH COMPARTMENTS, ON NAMED SUBTYPES. The T cell enrichment heatmap Zoey
     already has was made before the T labels were locked, so its rows are bare
     sub_leiden IDs. They are locked now (shiv_tier2_annotation_write.py), so
     this regenerates the T figure on 'CD8 effector memory' etc. and produces
     the myeloid one the same way. Same maths, same scale, same layout, so the
     two sit side by side and the T figure finally matches the myeloid one.

  2. TIER 2 UMAPS, three embeddings. Global, T subset, myeloid subset, each
     coloured by Tier 2 label and by subcluster.

PLOTTING. scatter_cat / color_map / cat_legend are ported from
shiv_plot_figures_updated_T_cell_corrected.py CELL 2 with ONE change: that
version closes over a module-level XY, which cannot work here because three
different embeddings are plotted in one run, so xy is a parameter now. Point
size, shuffle-before-draw, rasterised scatter, white-stroked on-data labels,
hidden spines, 300 DPI PDF + PNG are all unchanged.

THE GLOBAL TIER 2 LABEL IS COMPOSED, and this is the one judgement call worth
knowing about. No single obs column holds the full Tier 2 picture:
  tier2_celltype          fine T subtypes + COARSE Tier 1 for everything else
  tier2_myeloid_subtype   fine myeloid subtypes, 'not myeloid' for every T cell
So the global UMAP takes tier2_celltype and overwrites the myeloid cells with
their myeloid subtype. The composition happens IN MEMORY only; no object on
disk is modified. If you ever want that column persisted it belongs in an
annotation-write script, not a plotting script.

ENRICHMENT MATHS, unchanged from shiv_tier2_diagnostic_v2.py CELL 7:
log2( P(animal | row) / P(animal) ), NaN where the animal contributed no cells,
painted gray. The denominator is the COMPARTMENT on each figure, so a T log2 and
a myeloid log2 are not on a shared denominator and must not be compared directly.

Read-only. Spyder cells (# %%). Headless-safe (sc_pre env).

Author: Jake Lehle
Date: October 2026
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
import anndata as ad
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D
import warnings
warnings.filterwarnings("ignore")

# adjustText drives the non-overlapping on-plot labels. Warn-and-continue: if it
# is missing the UMAPs still render, falling back to plain centroid labels that
# can collide. Install into sc_pre with:  pip install adjustText
try:
    from adjustText import adjust_text
    HAVE_ADJUST_TEXT = True
except ImportError:
    HAVE_ADJUST_TEXT = False

WORKING_DIR = '/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV'
MERGED      = os.path.join(WORKING_DIR, 'shiv_starsolo', 'merged')
OUT_DIR     = os.path.join(WORKING_DIR, 'annotation_output', 'tier2_enrichment_umaps')

# global object: prefer the one carrying BOTH annotation passes
GLOBAL_CANDIDATES = [
    os.path.join(MERGED, 'shiv_host_tier2_myeloid_annotated.h5ad'),
    os.path.join(MERGED, 'shiv_host_tier2_annotated.h5ad'),
]

# ---- cohort ----
ANIMAL_COL = 'animal_id'
TISSUE_COL = 'tissue'
TIME_COL   = 'timepoint'
UMAP_KEY   = 'X_umap'
CLUSTER_KEY = 'sub_leiden'
ELITE      = '40707'
INFECTED_ANIMALS = ['39272', '40702', '40707', '41861', '41862']
CONTROL_ANIMALS  = ['34315', '41903']
EXCLUDE_ANIMALS  = ['Unknown']

# ---- scale: shared by both figures. Changing it changes BOTH, on purpose ----
ENRICH_CLIP     = 2.0
ANNOTATE_VALUES = True
ANNOT_FONTSIZE  = 14
HEX_MISSING     = '#E0E0E0'    # gray for NaN (animal contributed no cells)

# ---- what goes on the enrichment figures ----
DROP_EXCLUDE_ROWS = True    # 'EXCLUDE: ...' rows are not that lineage
DROP_QC_ROWS      = False   # QC rows kept but suffixed and sorted to the bottom
QC_SUFFIX         = '  [QC]'
MIN_CELLS_PER_ROW = 20

# ---- UMAP style, from shiv_plot_figures_updated_T_cell_corrected.py ----
PT_SIZE   = 10.0
PT_SIZE_SUBSET = 16.0        # subsets are smaller, so points can be bigger
LABEL_FS  = 22
LEGEND_FS = 18

# ---- on-plot labels (the gen_mpl_labels approach) ----
# Labels sit in the embedding at each group's median position and are then
# pushed apart by adjustText, with a leader line back to the centroid when a
# label has to travel. This replaces the side legend on every UMAP.
LABEL_ON_PLOT      = True
SHOW_SIDE_LEGEND   = False   # True draws BOTH, for a crowded global panel
LABELS_REQUIRE_ADJUSTTEXT = False
# False: no adjustText -> labels still drawn, at raw centroids, with a warning.
# True : no adjustText -> no labels at all, which is what your version does.
LABEL_FONTSIZE     = 22      # your snippet used 14; 22 matches the 28-34 house
                             # scale as closely as a 20-category UMAP allows
EXCLUDE_FROM_LABELS = ('None', 'Unknown', 'Unassigned', 'nan')
LABEL_ARROW_COLOR  = '#000000'
# white halo then a thin black keyline, so a label stays readable over any
# colour of point. Verbatim from your snippet.
COMBINED_EFFECTS = [
    pe.withStroke(linewidth=6, foreground='white'),
    pe.withStroke(linewidth=1, foreground='black'),
    pe.Normal(),
]
FALLBACK_PALETTE = ['#4C72B0', '#DD8452', '#55A868', '#C44E52', '#8172B3',
                    '#937860', '#DA8BC3', '#8C8C8C', '#CCB974', '#64B5CD',
                    '#E377C2', '#7F7F7F', '#BCBD22', '#17BECF']

# ---- Tier 2 colours (all hex) ----
# T palette lifted verbatim from shiv_tier2_annotation_write.py so the T cells
# keep the colours they already have in every figure Zoey has seen.
T_COLORS = {
    'Naive/resting T':       '#64B5F6',
    'CD4 central memory':    '#43A047',
    'CD8 effector memory':   '#FB8C00',
    'CD8 terminal effector': '#D32F2F',
    'QC: high-mito':         '#BDBDBD',
    'QC: T-B doublet':       '#9E9E9E',
}
# Myeloid had no per-subtype palette, only the generic HEX family dict in the
# diagnostic. Built here on that family logic so colour carries lineage:
# monocyte blues, macrophage red, DC greens, IFN orange, APC purple, gray junk.
MYELOID_COLORS = {
    'Classical Monocyte':              '#2C6E9B',
    'Monocyte, MHC-II-high':           '#5B9BD5',
    'Monocyte, IFN-stimulated':        '#E8A33D',
    'Non-classical / CD163+ Mono-Mac': '#1F4E79',
    'Tissue Macrophage (C1Q+/CD206+)': '#C4453C',
    'Antigen Presenting DC':           '#7B5EA7',
    'pDC':                             '#8FBF4D',
    'Migratory / Mature DC':           '#5E8C3A',
    'Unresolved myeloid':              '#BDBDBD',
    'EXCLUDE: T/NK contamination':     '#9E9E9E',
    'EXCLUDE: Platelets':              '#8A8A8A',
}
# coarse Tier 1 lineages that stay coarse on the global figure
GLOBAL_EXTRA = {
    'B cells': '#1976D2', 'NK cells': '#7B1FA2', 'Mast cells': '#455A64',
    'Erythrocytes': '#8D6E63', 'Platelets': '#C2185B', 'Unassigned': '#E0E0E0',
    'Monocytes': '#F9A825', 'Macrophages': '#00897B', 'DCs': '#6D4C41',
}

# ---- the two compartments ----
COMPARTMENTS = [
    {
        'name':      'T cells',
        'tag':       'Tcell',
        'path':      os.path.join(MERGED, 'shiv_Tcell_annotated.h5ad'),
        'row_key':   'tier2_subtype',
        'qc_col':    'qc_exclude',
        'colors':    T_COLORS,
        'axis_label': 'T cell subtype',
        # QC rows here are mito-high and T-B doublet, both animal-skewed
        'qc_note':   'clusters 8 and 9: high-mito (100% Indian controls) and '
                     'T-B doublet. Both are animal-skewed by construction.',
    },
    {
        'name':      'Myeloid',
        'tag':       'Myeloid',
        'path':      os.path.join(MERGED, 'shiv_Myeloid_annotated.h5ad'),
        'row_key':   'tier2_myeloid_subtype',
        'qc_col':    'qc_exclude_myeloid',
        'colors':    MYELOID_COLORS,
        'axis_label': 'myeloid subtype',
        'qc_note':   'clusters 0 and 4 are ~80-90% animal 34315, cluster 8 is '
                     '67% animal 41861. Their large values ARE the artifact.',
    },
]

# ---- global Tier 2 composition ----
GLOBAL_T_COL       = 'tier2_celltype'          # fine T + coarse everything else
GLOBAL_MYELOID_COL = 'tier2_myeloid_subtype'   # fine myeloid, 'not myeloid' elsewhere
NOT_MYELOID_TOKEN  = 'not myeloid'

plt.rcParams.update({
    'font.size': 28, 'axes.titlesize': 34, 'axes.labelsize': 30,
    'xtick.labelsize': 24, 'ytick.labelsize': 24, 'legend.fontsize': LEGEND_FS,
    'figure.dpi': 100, 'pdf.fonttype': 42, 'ps.fonttype': 42,
})
os.makedirs(OUT_DIR, exist_ok=True)
pd.set_option('display.width', 220)
pd.set_option('display.max_columns', 60)
pd.set_option('display.max_rows', 300)


def hr(t):
    print('\n' + '=' * 90)
    print(t)
    print('=' * 90)


def savefig(fig, name):
    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(OUT_DIR, f'{name}.{ext}'),
                    bbox_inches='tight', dpi=300)
    print(f'    saved: {name}.pdf / .png')


# %% ============================================================
# CELL 2 - HELPERS (enrichment maths + the ported plotting functions)
# ============================================================
hr('CELL 2 - helpers')


def log2_enrichment(row_series, group_series):
    """log2( P(group|row) / P(group) ), NaN where group absent from row.

    Verbatim from shiv_tier2_diagnostic_v2.py CELL 1 so both compartments use
    the same statistic the T figure was always built on. A zero cell gives
    log2(0) = -inf, sent to NaN and painted gray: 'this animal contributed no
    cells here', which is a different claim from 'depleted here'.
    """
    ct = pd.crosstab(row_series, group_series)
    p_obs = ct.div(ct.sum(axis=1), axis=0)
    p_base = ct.sum(axis=0) / ct.values.sum()
    e = np.log2(p_obs.div(p_base, axis=1)).replace([np.inf, -np.inf], np.nan)
    return ct, e


def color_map(obj, col, override=None):
    """Reuse uns['{col}_colors'] if present and length-matched, else override,
    else the fallback palette. Ported from the plotting script; the override arg
    is new so the locked Tier 2 palettes win over whatever scanpy last stored."""
    s = obj.obs[col].astype('category')
    cats = list(s.cat.categories)
    if override:
        missing = [c for c in cats if c not in override]
        cmap = {c: override.get(
            c, FALLBACK_PALETTE[i % len(FALLBACK_PALETTE)])
            for i, c in enumerate(cats)}
        if missing:
            print(f'    WARN: no colour defined for {missing}; using fallback')
        return cmap, cats
    stored = list(obj.uns.get(f'{col}_colors', []))
    if len(stored) >= len(cats):
        return {c: stored[i] for i, c in enumerate(cats)}, cats
    return {c: FALLBACK_PALETTE[i % len(FALLBACK_PALETTE)]
            for i, c in enumerate(cats)}, cats


def gen_mpl_labels(xy, series, ax, cmap=None, exclude=EXCLUDE_FROM_LABELS,
                   adjust_kwargs=None, text_kwargs=None, color_by_group=True):
    """Non-overlapping group labels placed in the embedding itself.

    Your gen_mpl_labels. Four deliberate differences from the version you sent,
    all of them stated here so nothing is a silent edit:

    1. SIGNATURE is (xy, series) not (adata, groupby). It has to be: this script
       labels three embeddings in one run and the global one is coloured by a
       composed column that exists only in memory, so there is no single
       adata.obsm['X_umap'] to read. Colours arrive as a dict rather than being
       looked up in adata.uns[groupby + '_colors'] for the same reason, and
       because the locked Tier 2 palettes are explicit dicts, not stored in uns.

    2. NO text_from_points. Yours defaults adjust_kwargs to
       {'text_from_points': False}, which was an adjustText 0.x parameter. It no
       longer exists in 1.x: it is swallowed into **kwargs, forwarded to the
       annotate call and ignored, so it reads like an active setting while doing
       nothing. The 1.x equivalent is avoid_self, which is what the call sites
       pass. Checked against adjustText 1.4.0.

    3. MISSING adjustText still labels. Yours warns and returns, so the figure
       comes out with no labels at all, which is a confusing way to fail right
       after asking for on-plot labels. This places them at raw centroids and
       says so. Set LABELS_REQUIRE_ADJUSTTEXT = True for your behaviour.

    4. EMPTY GROUPS are skipped. A categorical carrying an unused category
       yields an empty group, np.median of which is nan, and a label then lands
       at nan. Your caller avoids it with remove_unused_categories(); this does
       both, because the subsets here genuinely do carry unused categories.

    Returns the list of Text objects, or [] if nothing was labelled.
    """
    adjust_kwargs = dict(adjust_kwargs or {})
    text_kwargs = dict(text_kwargs or {})
    s = pd.Series(np.asarray(series).astype(str))

    medians = {}
    for g, idx in s.groupby(s).groups.items():
        if g in exclude:
            continue
        pos = np.asarray(idx, dtype=int)
        if len(pos) == 0:
            continue
        med = np.median(xy[pos, :2], axis=0)
        if not np.all(np.isfinite(med)):
            print(f'    WARN: {g} has no finite median position; not labelled')
            continue
        medians[g] = med
    if not medians:
        return []

    if not HAVE_ADJUST_TEXT and LABELS_REQUIRE_ADJUSTTEXT:
        print('    WARN: adjustText not installed; labels skipped entirely '
              '(LABELS_REQUIRE_ADJUSTTEXT is True).  pip install adjustText')
        return []

    # text colour defaults to black and is overridden per group, as yours does
    texts = []
    for g, (x, y) in medians.items():
        kw = dict(text_kwargs)
        kw['color'] = (cmap.get(g, '#000000')
                       if (color_by_group and cmap is not None) else '#000000')
        texts.append(ax.text(x=x, y=y, s=str(g), **kw))

    if HAVE_ADJUST_TEXT:
        adjust_text(texts, ax=ax, **adjust_kwargs)
    else:
        print('    WARN: adjustText not installed; labels placed at raw '
              'centroids and may overlap.  pip install adjustText')
    return texts


def scatter_cat(ax, xy, series, cmap, title, on_data=False, shuffle=True,
                pt_size=PT_SIZE, label_fontsize=LABEL_FONTSIZE):
    """Categorical UMAP scatter. Ported from
    shiv_plot_figures_updated_T_cell_corrected.py CELL 2; xy is a parameter
    instead of a module-level global so one run can plot three embeddings.

    on_data now routes through gen_mpl_labels, so labels land in the embedding
    and are pushed apart rather than being dropped at raw centroids where they
    collided (clusters 7/8/9 on the T subset did exactly that)."""
    # drop unused categories, as your caller does with
    # cat.remove_unused_categories(): an unused level otherwise reserves a
    # legend entry and a nan-positioned label
    s = series.astype('category').cat.remove_unused_categories()
    colors = s.map(cmap).values
    order = np.arange(len(s))
    if shuffle:
        np.random.default_rng(0).shuffle(order)
    ax.scatter(xy[order, 0], xy[order, 1], c=colors[order], s=pt_size,
               linewidths=0, rasterized=True)
    ax.set_title(title)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    if on_data:
        gen_mpl_labels(
            xy, s, ax, cmap=cmap,
            adjust_kwargs=dict(
                # avoid_self is the 1.x name for what text_from_points did in
                # 0.x: let a label sit on its own cluster rather than being
                # repelled by the points it is labelling.
                avoid_self=False,
                # shrinkA pulls the leader line back off the text box.
                # adjustText warns that arrows can strike through text without
                # it on a non-FancyArrowPatch transform, which is what a plain
                # axes gives us.
                arrowprops=dict(arrowstyle='-', color=LABEL_ARROW_COLOR,
                                lw=1.0, shrinkA=4)),
            text_kwargs=dict(fontsize=label_fontsize, fontweight='bold',
                             ha='center', va='center',
                             path_effects=COMBINED_EFFECTS),
            color_by_group=True)


def cat_legend(ax, cmap, ncol=1, fs=LEGEND_FS):
    """Ported unchanged."""
    handles = [Line2D([0], [0], marker='o', linestyle='', markersize=11,
                      markerfacecolor=cmap[c], markeredgecolor='none',
                      label=str(c))
               for c in cmap]
    ax.legend(handles=handles, loc='center left', bbox_to_anchor=(1.0, 0.5),
              frameon=False, fontsize=fs, ncol=ncol, handletextpad=0.3)


def is_qc_row(r):
    """A row is QC either because it carries the suffix or because its locked
    label already says so ('QC: high-mito'). Used for sorting QC rows to the
    bottom and for the do-not-read-as-biology note."""
    return QC_SUFFIX.strip() in r or str(r).startswith('QC')


def order_animals(cols):
    o = [a for a in INFECTED_ANIMALS if a in cols] + \
        [a for a in CONTROL_ANIMALS if a in cols]
    return o + [a for a in cols if a not in o]


print('  log2_enrichment  : from shiv_tier2_diagnostic_v2.py')
print('  scatter_cat      : from shiv_plot_figures_updated_T_cell_corrected.py')
print('                     (xy parameterised; three embeddings in one run)')


# %% ============================================================
# CELL 3 - ENRICHMENT HEATMAPS, ONE PER COMPARTMENT
# ============================================================
hr('CELL 3 - Tier 2 animal enrichment, both compartments')

loaded = {}
enrich_results = {}

for comp in COMPARTMENTS:
    name, tag = comp['name'], comp['tag']
    print(f'\n{"-" * 90}\n{name}\n{"-" * 90}')

    if not os.path.exists(comp['path']):
        print(f'  WARN: {os.path.basename(comp["path"])} not found; '
              f'{name} skipped')
        continue
    sb = ad.read_h5ad(comp['path'])
    loaded[name] = sb
    print(f'  loaded: {os.path.basename(comp["path"])}  ({sb.n_obs:,} cells)')

    row_key = comp['row_key']
    if row_key not in sb.obs.columns:
        print(f'  WARN: {row_key} absent; falling back to {CLUSTER_KEY}')
        row_key = CLUSTER_KEY
    if row_key not in sb.obs.columns or ANIMAL_COL not in sb.obs.columns:
        print(f'  WARN: required columns missing; {name} skipped')
        continue

    obs = sb.obs.copy()
    obs['_row'] = obs[row_key].astype(str)
    obs['_animal'] = obs[ANIMAL_COL].astype(str)
    n0 = len(obs)

    obs = obs[~obs['_animal'].isin(EXCLUDE_ANIMALS)]
    print(f'  excluded animals {EXCLUDE_ANIMALS}: {n0 - len(obs):,} cells')

    if DROP_EXCLUDE_ROWS:
        is_excl = obs['_row'].str.startswith('EXCLUDE')
        if is_excl.any():
            print('  dropped non-lineage EXCLUDE rows: '
                  f'{obs.loc[is_excl, "_row"].value_counts().to_dict()}')
        obs = obs[~is_excl]

    qc_col = comp['qc_col']
    if qc_col in obs.columns:
        qc = obs[qc_col].astype(bool)
        print(f'  {qc_col}: {int(qc.sum()):,} cells')
        if DROP_QC_ROWS:
            obs = obs[~qc]
            print('    DROP_QC_ROWS True -> removed')
        else:
            # Suffix PER CELL, not per label. A label can span both a clean and
            # a QC-flagged cluster: Classical Monocyte is cluster 2 (clean) plus
            # clusters 0 and 4 (79.9% / 89.5% animal 34315), and Unresolved
            # myeloid is clusters 3 and 10 (clean) plus cluster 8 (66.7% animal
            # 41861). Suffixing the whole label would tar the clean cells with
            # the confounded ones and hide the comparison that actually matters,
            # which is the clean row sitting flat while the QC row spikes on its
            # one animal. Splitting keeps both visible as separate rows.
            if qc.any():
                split = sorted(set(obs.loc[qc, '_row']) & set(obs.loc[~qc, '_row']))
                # a label that already says QC does not need the suffix twice
                needs = qc.values & ~obs['_row'].str.startswith('QC').values
                obs['_row'] = np.where(needs, obs['_row'] + QC_SUFFIX, obs['_row'])
                print(f'    kept and split, suffixed "{QC_SUFFIX.strip()}": '
                      f'{sorted(obs.loc[qc.values, "_row"].unique())}')
                if split:
                    print(f'    labels split across clean and QC rows: {split}')
                print(f'    {comp["qc_note"]}')
    else:
        print(f'  WARN: {qc_col} absent; no QC flagging applied')

    counts = obs['_row'].value_counts()
    thin = counts[counts < MIN_CELLS_PER_ROW]
    if len(thin):
        print(f'  rows under {MIN_CELLS_PER_ROW} cells, dropped: '
              f'{thin.to_dict()}')
        obs = obs[~obs['_row'].isin(thin.index)]

    ct, e = log2_enrichment(obs['_row'], obs['_animal'])
    cols = order_animals(list(e.columns))
    e, ct = e[cols], ct[cols]
    row_order = sorted(e.index, key=lambda r: (is_qc_row(r), r))
    e, ct = e.loc[row_order], ct.loc[row_order]
    enrich_results[name] = (e, ct)

    print(f'\n  {name} animal enrichment (log2 obs/exp):')
    print(e.round(2).to_string())
    stem = f'tier2_{tag}_animal_enrichment'
    e.round(3).to_csv(os.path.join(OUT_DIR, f'{stem}.csv'))
    ct.to_csv(os.path.join(OUT_DIR, f'{stem}_counts.csv'))
    print(f'    saved: {stem}.csv / {stem}_counts.csv')

    # --- gray cells, explained for the caption ---
    gray = [(r, a) for r in e.index for a in e.columns
            if np.isnan(e.loc[r, a])]
    if not gray:
        print('\n  no gray cells: every animal contributed to every row')
    else:
        print(f'\n  {len(gray)} gray cell(s), n = 0 not low:')
        has_t = TISSUE_COL in obs.columns
        if has_t:
            row_tis = (obs.groupby(['_row', TISSUE_COL], observed=True).size()
                       .unstack(fill_value=0))
            an_tis = (obs.groupby(['_animal', TISSUE_COL], observed=True).size()
                      .unstack(fill_value=0))
        for r, a in gray:
            reason = ''
            if has_t:
                tis = [t for t in row_tis.columns if row_tis.loc[r, t] > 0]
                miss = [t for t in tis
                        if a in an_tis.index and an_tis.loc[a, t] == 0]
                if miss and len(miss) == len(tis):
                    reason = (f'  <- sampling: row is {"/".join(tis)}-restricted'
                              f' and {a} has no {"/".join(miss)} cells here')
            print(f'    {r:<42s} x {a}: 0 of {int(ct[a].sum()):,}{reason}')

    # --- the figure ---
    fig, ax = plt.subplots(figsize=(1.6 * e.shape[1] + 8, 0.75 * e.shape[0] + 4))
    cmap = plt.cm.RdBu_r.copy()
    cmap.set_bad(HEX_MISSING)
    im = ax.imshow(e.values, cmap=cmap, vmin=-ENRICH_CLIP, vmax=ENRICH_CLIP,
                   aspect='auto')
    ax.set_xticks(range(e.shape[1]))
    ax.set_xticklabels([f'{a}*' if a == ELITE else a for a in e.columns],
                       rotation=45, ha='right')
    ax.set_yticks(range(e.shape[0]))
    ax.set_yticklabels([f'{r}  (n={int(ct.loc[r].sum()):,})' for r in e.index])
    ax.set_xlabel('animal')
    ax.set_ylabel(comp['axis_label'])
    ax.set_title(f'{name} animal enrichment (log2 obs/exp); * = elite controller')
    if ANNOTATE_VALUES:
        for i in range(e.shape[0]):
            for j in range(e.shape[1]):
                v = e.values[i, j]
                if not np.isnan(v):
                    ax.text(j, i, f'{v:.1f}', ha='center', va='center',
                            fontsize=ANNOT_FONTSIZE,
                            color='white' if abs(v) > ENRICH_CLIP * 0.6
                            else 'black')
    cb = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label('log2 enrichment', fontsize=22)
    cb.ax.tick_params(labelsize=18)
    plt.tight_layout()
    savefig(fig, f'{stem}_heatmap')
    plt.close(fig)

    # --- what it says ---
    print(f'\n  {name} rows with any animal at |log2| >= 1.0:')
    for r in e.index:
        row = e.loc[r].dropna()
        if len(row) and row.abs().max() >= 1.0:
            top = row.idxmax()
            flag = '  <ELITE>' if top == ELITE else ''
            note = '  [QC-confounded]' if is_qc_row(r) else ''
            print(f'    {r:<42s} top {top} at {row[top]:+.2f}{flag}{note}')
    if ELITE in e.columns:
        print(f'\n  elite controller {ELITE}, most to least enriched:')
        for r, v in e[ELITE].dropna().sort_values(ascending=False).items():
            print(f'    {v:+.2f}  {r}')


# %% ============================================================
# CELL 4 - GLOBAL TIER 2 UMAP (composed label)
# ============================================================
hr('CELL 4 - global Tier 2 UMAP')

gpath = next((p for p in GLOBAL_CANDIDATES if os.path.exists(p)), None)
glob = None
if gpath is None:
    print('  WARN: no global object found; global UMAP skipped')
    print('  looked for:\n    ' + '\n    '.join(GLOBAL_CANDIDATES))
else:
    glob = ad.read_h5ad(gpath)
    print(f'  loaded: {os.path.basename(gpath)}  ({glob.n_obs:,} cells)')

    if UMAP_KEY not in glob.obsm:
        print(f'  WARN: {UMAP_KEY} not in obsm; global UMAP skipped')
        glob = None

if glob is not None:
    # compose: fine T labels from tier2_celltype, fine myeloid overwritten in.
    # In memory only. The object on disk is not touched.
    if GLOBAL_T_COL in glob.obs.columns:
        composed = glob.obs[GLOBAL_T_COL].astype(str).copy()
        src = GLOBAL_T_COL
    else:
        composed = glob.obs['tier1_celltype'].astype(str).copy()
        src = 'tier1_celltype'
        print(f'  WARN: {GLOBAL_T_COL} absent; composing from tier1_celltype')

    n_my = 0
    if GLOBAL_MYELOID_COL in glob.obs.columns:
        my = glob.obs[GLOBAL_MYELOID_COL].astype(str)
        take = (my != NOT_MYELOID_TOKEN) & my.notna() & (my != 'nan')
        composed.loc[take] = my.loc[take]
        n_my = int(take.sum())
        print(f'  composed: {src} + {GLOBAL_MYELOID_COL} '
              f'({n_my:,} myeloid cells refined)')
    else:
        print(f'  WARN: {GLOBAL_MYELOID_COL} absent; myeloid stays coarse')

    glob.obs['_tier2_composed'] = pd.Categorical(composed)
    print('\n  global Tier 2 composition:')
    print(glob.obs['_tier2_composed'].value_counts().to_string())

    palette = {**GLOBAL_EXTRA, **T_COLORS, **MYELOID_COLORS}
    g_cmap, g_cats = color_map(glob, '_tier2_composed', override=palette)
    GXY = np.asarray(glob.obsm[UMAP_KEY])

    fig, ax = plt.subplots(figsize=(18, 16))
    scatter_cat(ax, GXY, glob.obs['_tier2_composed'], g_cmap,
                'Tier 2 cell types, global embedding', pt_size=PT_SIZE,
                on_data=LABEL_ON_PLOT)
    if SHOW_SIDE_LEGEND or not LABEL_ON_PLOT:
        cat_legend(ax, g_cmap, ncol=1)
    plt.tight_layout()
    savefig(fig, 'tier2_global_umap_celltypes')
    plt.close(fig)

    # the same embedding by global Leiden cluster, for the cluster-to-label read
    if 'clusters' in glob.obs.columns:
        c_cmap, _ = color_map(glob, 'clusters')
        fig, ax = plt.subplots(figsize=(16, 14))
        scatter_cat(ax, GXY, glob.obs['clusters'], c_cmap,
                    'Global Leiden clusters', on_data=True, pt_size=PT_SIZE)
        plt.tight_layout()
        savefig(fig, 'tier2_global_umap_clusters')
        plt.close(fig)

        xt = pd.crosstab(glob.obs['clusters'], glob.obs['_tier2_composed'])
        xt.to_csv(os.path.join(OUT_DIR, 'global_cluster_by_tier2.csv'))
        print('    saved: global_cluster_by_tier2.csv')


# %% ============================================================
# CELL 5 - SUBSET UMAPS: subcluster AND Tier 2 label, per compartment
# ============================================================
hr('CELL 5 - subset UMAPs')

for comp in COMPARTMENTS:
    name, tag = comp['name'], comp['tag']
    sb = loaded.get(name)
    if sb is None:
        print(f'  {name}: object not loaded; skipped')
        continue
    if UMAP_KEY not in sb.obsm:
        print(f'  WARN: {name} has no {UMAP_KEY}; subset UMAP skipped')
        continue

    XY = np.asarray(sb.obsm[UMAP_KEY])
    row_key = comp['row_key'] if comp['row_key'] in sb.obs.columns else CLUSTER_KEY
    have_cl = CLUSTER_KEY in sb.obs.columns
    ncol = 2 if have_cl else 1

    fig, axes = plt.subplots(1, ncol, figsize=(16 * ncol, 13))
    axes = np.atleast_1d(axes)

    lab_cmap, _ = color_map(sb, row_key, override=comp['colors'])
    scatter_cat(axes[0], XY, sb.obs[row_key], lab_cmap,
                f'{name} Tier 2 subtypes', pt_size=PT_SIZE_SUBSET,
                on_data=LABEL_ON_PLOT)

    if have_cl:
        cl_cmap, _ = color_map(sb, CLUSTER_KEY)
        scatter_cat(axes[1], XY, sb.obs[CLUSTER_KEY].astype(str), cl_cmap,
                    f'{name} subclusters', on_data=True,
                    pt_size=PT_SIZE_SUBSET)

    # With labels in the embedding the legend is redundant. If you do turn it
    # back on, note cat_legend anchors to (1.0, 0.5) of the axis it is given,
    # which on a two-panel figure is the GAP BETWEEN the panels, so it has to
    # go on the LAST axis to land outside the figure.
    if SHOW_SIDE_LEGEND or not LABEL_ON_PLOT:
        cat_legend(axes[-1], lab_cmap, ncol=1)

    plt.tight_layout()
    savefig(fig, f'tier2_{tag}_umap_subtype_and_subcluster')
    plt.close(fig)

    # the breakout that makes the two panels readable together
    if have_cl:
        xt = pd.crosstab(sb.obs[CLUSTER_KEY].astype(str),
                         sb.obs[row_key].astype(str))
        xt.to_csv(os.path.join(OUT_DIR, f'tier2_{tag}_cluster_by_subtype.csv'))
        print(f'    saved: tier2_{tag}_cluster_by_subtype.csv')
        print(f'  {name} subcluster -> subtype:')
        for cl in sorted(xt.index, key=lambda x: int(x) if x.isdigit() else 1e9):
            lab = xt.loc[cl].idxmax()
            print(f'    {cl:>3s} -> {lab}  (n={int(xt.loc[cl].sum()):,})')


# %% ============================================================
# CELL 6 - SUMMARY
# ============================================================
hr('CELL 6 - summary')

print(f"""  FIGURES in {OUT_DIR}

  tier2_Tcell_animal_enrichment_heatmap      T subtypes x animal
  tier2_Myeloid_animal_enrichment_heatmap    myeloid subtypes x animal
  tier2_global_umap_celltypes                composed Tier 2, global embedding
  tier2_global_umap_clusters                 global Leiden, labels on data
  tier2_Tcell_umap_subtype_and_subcluster    T subset, both colourings
  tier2_Myeloid_umap_subtype_and_subcluster  myeloid subset, both colourings

  CAVEATS for the captions
  ------------------------
  - The two enrichment figures share a scale (+/-{ENRICH_CLIP:.1f} log2) but NOT a
    denominator. Each is normalised within its own compartment, so a T value and
    a myeloid value are not comparable as numbers. Say so in the legend.
  - Gray = n of 0, not depletion. A depleted-but-present combination is blue.
  - The global Tier 2 label is composed in memory from two columns because no
    single column holds it. Nothing on disk was modified.
  - BBKNN corrected on library, not animal, so surviving animal structure is
    biological, but it is still confounded with tissue and timepoint sampling.
    These are scouting figures, suggesting where to look, not tests.

  LABELS: on-plot via gen_mpl_labels + adjustText, coloured by group, leader
  lines where a label had to travel. adjustText present: {HAVE_ADJUST_TEXT}.
  {'' if HAVE_ADJUST_TEXT else 'INSTALL IT (pip install adjustText) or labels sit at raw centroids and collide.'}
  Side legends are off (SHOW_SIDE_LEGEND); set it True to draw both.
""")
print('=' * 90)
print('TIER 2 ENRICHMENT + UMAPS COMPLETE - no object was modified.')
print('=' * 90)
