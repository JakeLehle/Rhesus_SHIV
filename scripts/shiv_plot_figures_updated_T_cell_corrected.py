#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Rhesus SHIV — Figure Set on the Final Embedded Object (Read-Only)
=================================================================
Publication/committee figures for Zoey off shiv_host_final_embedded.h5ad. Pure
plotting pass: reads the BBKNN-corrected UMAP, the tier1 labels, the Leiden
clusters, and the CLR ADT already on the object, and writes figures. The object
is NOT modified.

Figures
-------
  1  cell types + Leiden clusters on the UMAP (labels on data)
  2  the 9 usable surface markers on the UMAP (CLR feature plots, 3x3)
  3  metadata on the UMAP (tissue / timepoint / animal / library)
  4  cell-type composition, stacked bars, by tissue x timepoint

Choices baked in
----------------
  - Cell-type colors reuse uns['tier1_celltype_colors'] for consistency with
    prior figures; categoricals without stored colors get a colorblind-safe
    fallback palette.
  - CD163/CD169 are left OUT of the surface panel (they read blank pending BD
    Ab-seq recovery); noted in the caption rather than shown as empty panels.
  - The CD8 tag is anti-CD8a and cannot separate CD8aa NK from CD3-dropout T,
    so its panel is titled "CD8 (anti-CD8a)"; it is a surface map, not a
    CD8 T-cell map, and nothing here implies the CD4/CD8 split is settled.
  - Composition is the infected trajectory Pre -> 21 DPI -> Necropsy, split
    PBMC vs LN.

Figures saved as BOTH .pdf and .png at 300 DPI, 28-34pt text, hex colors.
Spyder cells (# %%). Runs headless under SLURM (sc_pre env).

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
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D
import warnings
warnings.filterwarnings("ignore")

WORKING_DIR = '/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV'
MERGED      = os.path.join(WORKING_DIR, 'shiv_starsolo', 'merged')
ADATA_IN    = os.path.join(MERGED, 'shiv_Tcell_annotated.h5ad')


OUT_DIR     = os.path.join(WORKING_DIR, 'annotation_output', 'figures')

UMAP_KEY  = 'X_umap'
CELLTYPE  = 'tier1_celltype'
CLUSTERS  = 'clusters'

# The 9 usable ADT tags (CD163/CD169 excluded: blank pending BD Ab-seq)
ADT_USABLE = ['CD28.2_TotalSeqC', 'CD69_TotalSeqC', 'CX3CR1_TotalSeqC',
              'CD1c_TotalSeqC', 'CD11b_TotalSeqC', 'CD206_TotalSeqC',
              'CD68_TotalSeqC', 'CD8_TotalSeqC', 'CD4_TotalSeqC']
ADT_DISPLAY = {
    'CD28.2_TotalSeqC': 'CD28',  'CD69_TotalSeqC': 'CD69',
    'CX3CR1_TotalSeqC': 'CX3CR1', 'CD1c_TotalSeqC': 'CD1c',
    'CD11b_TotalSeqC': 'CD11b',  'CD206_TotalSeqC': 'CD206',
    'CD68_TotalSeqC': 'CD68',    'CD8_TotalSeqC': 'CD8 (anti-CD8a)',
    'CD4_TotalSeqC': 'CD4',
}

# Composition axes
TIMEPOINT_ORDER = ['Pre', '21 DPI', 'Necropsy']   # infected trajectory
TISSUES         = ['PBMC', 'LN']

# TCR / animal figures
SUBTYPE_COL      = 'receptor_subtype'
ALPHA_BETA       = 'TRA+TRB'
ELITE            = '40707'
INFECTED_ANIMALS = ['39272', '40702', '40707', '41861', '41862']
CONTROL_ANIMALS  = ['34315', '41903']
EXCLUDE_ANIMALS  = ['Unknown']
ELITE_COLOR      = '#F9A825'
INFECTED_COLOR   = '#1976D2'
CONTROL_COLOR    = '#90A4AE'
ENRICH_CLIP      = 2.0            # symmetric log2 enrichment scale for the heatmap

# Style
CMAP_ADT   = 'magma_r'          # continuous surface-marker colormap
PT_SIZE    = 10.0
ADT_PMIN, ADT_PMAX = 2, 98    # per-panel percentile scaling for CLR
LABEL_FS   = 22               # on-data label font
LEGEND_FS  = 18               # dense-legend font
FALLBACK_PALETTE = ['#4C72B0', '#DD8452', '#55A868', '#C44E52', '#8172B3',
                    '#937860', '#DA8BC3', '#8C8C8C', '#CCB974', '#64B5CD',
                    '#E377C2', '#7F7F7F', '#BCBD22', '#17BECF']

plt.rcParams.update({
    'font.size': 28, 'axes.titlesize': 34, 'axes.labelsize': 30,
    'xtick.labelsize': 28, 'ytick.labelsize': 28, 'legend.fontsize': LEGEND_FS,
    'figure.dpi': 100, 'pdf.fonttype': 42, 'ps.fonttype': 42,
})
os.makedirs(OUT_DIR, exist_ok=True)


def savefig(fig, name):
    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(OUT_DIR, f'{name}.{ext}'),
                    bbox_inches='tight', dpi=300)
    print(f"    saved: {name}.pdf / .png")


# %% ============================================================
# CELL 2 — LOAD + COLOR MAPS
# ============================================================
print("=" * 80)
print("CELL 2 — LOAD + COLOR MAPS")
print("=" * 80)

adata = sc.read_h5ad(ADATA_IN)
print(f"  {adata.n_obs:,} cells x {adata.n_vars:,} genes")

if UMAP_KEY not in adata.obsm:
    raise RuntimeError(f"{UMAP_KEY} not in obsm; nothing to plot.")
XY = np.asarray(adata.obsm[UMAP_KEY])

adt_names = list(adata.uns.get('adt_names', []))
have_adt = ('protein_clr' in adata.obsm) and len(adt_names) > 0
CLR = np.asarray(adata.obsm['protein_clr'], dtype=float) if have_adt else None
if not have_adt:
    print("  WARN: protein_clr / adt_names missing; surface-marker figure skipped.")


def color_map(col):
    """Reuse uns['{col}_colors'] if present and length-matched, else fallback."""
    s = adata.obs[col].astype('category')
    cats = list(s.cat.categories)
    key = f'{col}_colors'
    stored = list(adata.uns.get(key, []))
    if len(stored) >= len(cats):
        return {c: stored[i] for i, c in enumerate(cats)}, cats
    return {c: FALLBACK_PALETTE[i % len(FALLBACK_PALETTE)]
            for i, c in enumerate(cats)}, cats


def scatter_cat(ax, series, cmap, title, on_data=False, shuffle=True):
    s = series.astype('category')
    colors = s.map(cmap).values
    order = np.arange(len(s))
    if shuffle:
        np.random.default_rng(0).shuffle(order)
    ax.scatter(XY[order, 0], XY[order, 1], c=colors[order], s=PT_SIZE,
               linewidths=0, rasterized=True)
    ax.set_title(title)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    if on_data:
        for c in s.cat.categories:
            m = (s == c).values
            if m.sum() == 0:
                continue
            cx, cy = np.median(XY[m, 0]), np.median(XY[m, 1])
            t = ax.text(cx, cy, str(c), fontsize=LABEL_FS, fontweight='bold',
                        ha='center', va='center')
            t.set_path_effects([pe.withStroke(linewidth=3, foreground='white')])


def scatter_cont(ax, vals, title):
    order = np.argsort(vals)   # high values on top
    vmin, vmax = np.percentile(vals, [ADT_PMIN, ADT_PMAX])
    sm = ax.scatter(XY[order, 0], XY[order, 1], c=vals[order], s=PT_SIZE,
                    cmap=CMAP_ADT, vmin=vmin, vmax=vmax, linewidths=0,
                    rasterized=True)
    ax.set_title(title)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    cb = plt.colorbar(sm, ax=ax, fraction=0.046, pad=0.04)
    cb.ax.tick_params(labelsize=20)
    cb.outline.set_visible(False)


def cat_legend(ax, cmap, ncol=1, fs=LEGEND_FS):
    handles = [Line2D([0], [0], marker='o', linestyle='', markersize=11,
                      markerfacecolor=cmap[c], markeredgecolor='none', label=str(c))
               for c in cmap]
    ax.legend(handles=handles, loc='center left', bbox_to_anchor=(1.0, 0.5),
              frameon=False, fontsize=fs, ncol=ncol, handletextpad=0.3)


ct_cmap, ct_cats = color_map(CELLTYPE)
print(f"  cell types ({len(ct_cats)}): {ct_cats}")
print(f"  using stored colors: {'tier1_celltype_colors' in adata.uns}")


# %% ============================================================
# CELL 3 — FIGURE 1 (cell types + Leiden clusters)
# ============================================================
print("\n" + "=" * 80)
print("CELL 3 — FIGURE 1 (cell types + clusters)")
print("=" * 80)

fig, axes = plt.subplots(1, 2, figsize=(26, 12))
scatter_cat(axes[0], adata.obs[CELLTYPE], ct_cmap,
            'Cell types (Tier 1)', on_data=True)
cl_cmap, _ = color_map(CLUSTERS)
scatter_cat(axes[1], adata.obs[CLUSTERS], cl_cmap,
            f'Leiden clusters (n={adata.obs[CLUSTERS].nunique()})', on_data=True)
plt.tight_layout()
savefig(fig, 'fig1_celltypes_clusters_umap')
plt.close(fig)


# %% ============================================================
# CELL 4 — FIGURE 2 (9 usable surface markers, CLR)
# ============================================================
print("\n" + "=" * 80)
print("CELL 4 — FIGURE 2 (surface markers, CLR)")
print("=" * 80)

if have_adt:
    present_tags = [t for t in ADT_USABLE if t in adt_names]
    missing = [t for t in ADT_USABLE if t not in adt_names]
    if missing:
        print(f"  WARN: tags not found, skipped: {missing}")
    fig, axes = plt.subplots(3, 3, figsize=(26, 24))
    for ax, tag in zip(axes.ravel(), present_tags):
        vals = CLR[:, adt_names.index(tag)]
        scatter_cont(ax, vals, ADT_DISPLAY.get(tag, tag.split('_')[0]))
    for ax in axes.ravel()[len(present_tags):]:
        ax.axis('off')
    fig.suptitle('Surface markers (CLR)   |   CD163, CD169 pending BD Ab-seq recovery',
                 fontsize=30, y=0.99)
    plt.tight_layout(rect=[0, 0, 1, 0.98])
    savefig(fig, 'fig2_surface_markers_umap')
    plt.close(fig)
else:
    print("  skipped (no ADT).")


# %% ============================================================
# CELL 5 — FIGURE 3 (metadata: tissue / timepoint / animal / library)
# ============================================================
print("\n" + "=" * 80)
print("CELL 5 — FIGURE 3 (metadata)")
print("=" * 80)

meta_cols = [c for c in ['tissue', 'timepoint', 'animal_id', 'library']
             if c in adata.obs.columns]
fig, axes = plt.subplots(2, 2, figsize=(24, 20))
for ax, col in zip(axes.ravel(), meta_cols):
    cmap, cats = color_map(col)
    scatter_cat(ax, adata.obs[col], cmap, col, on_data=False)
    ncol = 2 if len(cats) > 8 else 1
    cat_legend(ax, cmap, ncol=ncol)
for ax in axes.ravel()[len(meta_cols):]:
    ax.axis('off')
plt.tight_layout()
savefig(fig, 'fig3_metadata_umap')
plt.close(fig)


# %% ============================================================
# CELL 6 — FIGURE 4 (composition by tissue x timepoint)
# ============================================================
print("\n" + "=" * 80)
print("CELL 6 — FIGURE 4 (composition by tissue x timepoint)")
print("=" * 80)

if 'tissue' in adata.obs.columns and 'timepoint' in adata.obs.columns:
    obs = adata.obs
    tissues = [t for t in TISSUES if t in obs['tissue'].unique().tolist()]
    fig, axes = plt.subplots(1, len(tissues), figsize=(11 * len(tissues), 10),
                             squeeze=False)
    axes = axes.ravel()
    comp_tables = {}
    for ax, tissue in zip(axes, tissues):
        sub = obs[(obs['tissue'] == tissue) &
                  (obs['timepoint'].isin(TIMEPOINT_ORDER))]
        ct = pd.crosstab(sub['timepoint'], sub[CELLTYPE], normalize='index')
        ct = ct.reindex(index=TIMEPOINT_ORDER).fillna(0.0)
        ct = ct.reindex(columns=ct_cats, fill_value=0.0)
        counts = pd.crosstab(sub['timepoint'], sub[CELLTYPE]).reindex(
            index=TIMEPOINT_ORDER).fillna(0).sum(axis=1).astype(int)
        comp_tables[tissue] = ct

        bottom = np.zeros(len(ct))
        x = np.arange(len(ct))
        for celltype in ct_cats:
            ax.bar(x, ct[celltype].values, bottom=bottom, width=0.7,
                   color=ct_cmap[celltype], label=celltype)
            bottom += ct[celltype].values
        ax.set_xticks(x)
        ax.set_xticklabels([f"{tp}\n(n={counts.get(tp, 0):,})"
                            for tp in TIMEPOINT_ORDER])
        ax.set_title(tissue)
        ax.set_ylabel('cell-type proportion')
        ax.set_ylim(0, 1)
        for sp in ('top', 'right'):
            ax.spines[sp].set_visible(False)

    handles = [Line2D([0], [0], marker='s', linestyle='', markersize=14,
                      markerfacecolor=ct_cmap[c], markeredgecolor='none', label=c)
               for c in ct_cats]
    fig.legend(handles=handles, loc='center left', bbox_to_anchor=(1.0, 0.5),
               frameon=False, fontsize=LEGEND_FS)
    fig.suptitle('Cell-type composition across the infected trajectory',
                 fontsize=30, y=1.02)
    plt.tight_layout()
    savefig(fig, 'fig4_composition_tissue_timepoint')
    plt.close(fig)

    # export the proportions for Zoey
    out = pd.concat({k: v for k, v in comp_tables.items()}, axis=0)
    out.index.names = ['tissue', 'timepoint']
    out.round(4).to_csv(os.path.join(OUT_DIR, 'fig4_composition_proportions.csv'))
    print("    saved: fig4_composition_proportions.csv")
else:
    print("  WARN: tissue/timepoint missing; composition figure skipped.")


# %% ============================================================
# CELL 7 — FIGURE 5 (alpha-beta T cells by animal: raw vs relative)
# ============================================================
print("\n" + "=" * 80)
print("CELL 7 — FIGURE 5 (alpha-beta T by animal)")
print("=" * 80)

if SUBTYPE_COL in adata.obs.columns and 'animal_id' in adata.obs.columns:
    is_ab = ((adata.obs[SUBTYPE_COL].astype(str) == ALPHA_BETA) &
             (adata.obs[CELLTYPE].astype(str) == 'T cells')).values
    animal = adata.obs['animal_id'].astype(str).values
    df = pd.DataFrame({'animal': animal, 'ab': is_ab})
    df = df[~df['animal'].isin(EXCLUDE_ANIMALS)]
    grp = df.groupby('animal').agg(total=('ab', 'size'), ab=('ab', 'sum'))
    grp['ab_prop'] = grp['ab'] / grp['total']

    def a_color(a):
        if a == ELITE:
            return ELITE_COLOR
        if a in CONTROL_ANIMALS:
            return CONTROL_COLOR
        return INFECTED_COLOR

    def bar_panel(ax, series, ylabel, title):
        s = series.sort_values(ascending=False)
        colors = [a_color(a) for a in s.index]
        ax.bar(range(len(s)), s.values, color=colors)
        ax.set_xticks(range(len(s)))
        ax.set_xticklabels([f"{a}*" if a == ELITE else a for a in s.index],
                           rotation=45, ha='right')
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        for sp in ('top', 'right'):
            ax.spines[sp].set_visible(False)
        rank = list(s.index).index(ELITE) + 1 if ELITE in s.index else None
        return rank

    fig, axes = plt.subplots(1, 2, figsize=(24, 10))
    r_raw = bar_panel(axes[0], grp['ab'], 'alpha-beta T cells',
                      'Alpha-beta T cells per animal (raw)')
    r_rel = bar_panel(axes[1], grp['ab_prop'],
                      'alpha-beta T / animal cells',
                      "Alpha-beta T as fraction of animal's cells (relative)")
    handles = [Line2D([0], [0], marker='s', linestyle='', markersize=14,
                      markerfacecolor=c, markeredgecolor='none', label=lab)
               for c, lab in [(ELITE_COLOR, 'elite controller 40707'),
                              (INFECTED_COLOR, 'infected'),
                              (CONTROL_COLOR, 'non-infected control')]]
    fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, 1.06),
               frameon=False, ncol=3, fontsize=LEGEND_FS)
    plt.tight_layout()
    savefig(fig, 'fig5_alphabeta_T_by_animal')
    plt.close(fig)

    grp.round(4).to_csv(os.path.join(OUT_DIR, 'fig5_alphabeta_T_by_animal.csv'))
    print("    saved: fig5_alphabeta_T_by_animal.csv")
    print(f"  elite controller {ELITE} rank: {r_raw} of {len(grp)} (raw), "
          f"{r_rel} of {len(grp)} (relative)")
    print("  (relative = alpha-beta T normalized to that animal's total captured")
    print("   cells; aggregates across tissue/timepoint. Split by timepoint on request.)")
else:
    print("  WARN: receptor_subtype / animal_id missing; Figure 5 skipped.")


# %% ============================================================
# CELL 8 — FIGURE 6 (animal enrichment per cluster; Tier 2 scouting)
# ============================================================
print("\n" + "=" * 80)
print("CELL 8 — FIGURE 6 (animal enrichment per cluster)")
print("=" * 80)
print("  NOTE: global 18-cluster view, NOT Tier 2. BBKNN corrected library")
print("        (technical) batch, not animal (biology), so animal enrichment")
print("        here is biological signal and marks where Tier 2 should dig.")

if CLUSTERS in adata.obs.columns and 'animal_id' in adata.obs.columns:
    a = adata.obs['animal_id'].astype(str)
    keep = ~a.isin(EXCLUDE_ANIMALS)
    ct = pd.crosstab(adata.obs.loc[keep, CLUSTERS], a[keep])
    # order animals: infected first (elite flagged), then controls
    order = [x for x in INFECTED_ANIMALS if x in ct.columns] + \
            [x for x in CONTROL_ANIMALS if x in ct.columns] + \
            [x for x in ct.columns if x not in INFECTED_ANIMALS + CONTROL_ANIMALS]
    ct = ct[order]
    # observed P(animal|cluster) vs baseline P(animal); log2 enrichment
    p_obs = ct.div(ct.sum(axis=1), axis=0)
    p_base = ct.sum(axis=0) / ct.values.sum()
    enrich = np.log2(p_obs.div(p_base, axis=1))
    enrich = enrich.replace([np.inf, -np.inf], np.nan)  # animal absent -> blank
    enrich.round(3).to_csv(os.path.join(OUT_DIR, 'fig6_cluster_animal_enrichment.csv'))

    fig, ax = plt.subplots(figsize=(1.4 * len(order) + 6, 0.6 * len(enrich) + 4))
    cmap = plt.cm.RdBu_r.copy()
    cmap.set_bad('#E0E0E0')
    im = ax.imshow(enrich.values, cmap=cmap, vmin=-ENRICH_CLIP, vmax=ENRICH_CLIP,
                   aspect='auto')
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([f"{a_}*" if a_ == ELITE else a_ for a_ in order],
                       rotation=45, ha='right')
    ax.set_yticks(range(len(enrich)))
    ax.set_yticklabels(enrich.index)
    ax.set_xlabel('animal')
    ax.set_ylabel('cluster')
    ax.set_title('Animal enrichment per cluster  (log2 obs/expected)\n'
                 'global clustering, Tier 2 scouting; * = elite controller')
    # annotate values
    for i in range(enrich.shape[0]):
        for j in range(enrich.shape[1]):
            v = enrich.values[i, j]
            if not np.isnan(v):
                ax.text(j, i, f"{v:.1f}", ha='center', va='center',
                        fontsize=14,
                        color='white' if abs(v) > ENRICH_CLIP * 0.6 else 'black')
    cb = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label('log2 enrichment', fontsize=22)
    cb.ax.tick_params(labelsize=18)
    plt.tight_layout()
    savefig(fig, 'fig6_cluster_animal_enrichment')
    plt.close(fig)
    print("    saved: fig6_cluster_animal_enrichment.csv")

    # flag the clusters most skewed toward the elite controller
    if ELITE in enrich.columns:
        top = enrich[ELITE].dropna().sort_values(ascending=False).head(3)
        print(f"  clusters most enriched for {ELITE} (log2 obs/exp):")
        for cl, v in top.items():
            dom = adata.obs.loc[adata.obs[CLUSTERS] == cl, CELLTYPE].mode()
            dom = dom.iloc[0] if len(dom) else '?'
            print(f"    cluster {cl}: {v:.2f}  (dominant tier1: {dom})")
else:
    print("  WARN: clusters / animal_id missing; Figure 6 skipped.")


# %% ============================================================
# CELL 9 — DONE
# ============================================================
print("\n" + "=" * 80)
print("DONE — figures in", OUT_DIR)
print("=" * 80)
print("""  fig1  cell types + Leiden clusters (labels on data)
  fig2  9 surface markers, CLR (CD163/CD169 excluded, blank pending BD Ab-seq)
  fig3  metadata: tissue / timepoint / animal / library
  fig4  composition by tissue x timepoint (+ proportions CSV)
  fig5  alpha-beta T by animal, raw vs relative (+ CSV)
  fig6  animal enrichment per cluster, Tier 2 scouting (+ CSV)

  Caveats for the caption: the CD8 panel is anti-CD8a surface signal (cannot
  split CD8aa NK from CD3-dropout T); CD4/CD8 lineage is a Tier 2 question, not
  settled here. fig6 is the global 18-cluster view, not Tier 2 annotation.
  Read-only: the object was not modified.""")
