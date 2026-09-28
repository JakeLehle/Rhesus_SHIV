#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Rhesus SHIV - enrichment figures: KEGG bubble, GO term network, volcano
========================================================================
Python equivalents of the three R figures Jake uses, built from files already
on disk. Recomputes no statistics.

WHAT EACH FIGURE IS, AND WHAT IT REPLACES
------------------------------------------
1. KEGG BUBBLE          <- clusterProfiler::dotplot(kegg_enrich)
   One panel per direction, up and down side by side on a shared term axis.
   x = fold enrichment, y = pathway, bubble area = overlap (the Count in
   clusterProfiler), colour = -log10 corrected p. Up is #C4453C, down is
   #2C6E9B, matching the palette the DEG and enrichment scripts already use.

   One deliberate change from dotplot(). clusterProfiler puts GeneRatio
   (k/N) on x. That is not comparable between an up list of 19 genes and a
   down list of 301, because N differs, so the same biology lands at a very
   different x. Fold enrichment, (k/N)/(n/M), divides that out and IS
   comparable across the two panels. Both columns are written to the
   accompanying CSV so you can plot GeneRatio instead if Zoey prefers the
   familiar look.

2. GO TERM NETWORK      <- emapplot(pairwise_termsim(go_enrich))
   This is the one you said made GO reviewable. emapplot's default similarity
   is the JACCARD INDEX of the two terms' hit-gene sets, and that is exactly
   what is reproduced here: nodes are enriched GO terms, an edge is drawn when
   JC >= JACCARD_MIN, edge width scales with JC, node area scales with overlap
   and node colour is -log10 corrected p.

   This is only possible because the hit-gene lists validated on your run.
   Every term's genes column was confirmed against g:Profiler's own
   intersection_size, so the Jaccard values are computed on real gene sets.
   Had that check failed, this figure could not be built honestly and the
   script would have refused rather than drawn an empty-looking network.

   Connected components are the point of the figure: a component is a block of
   GO terms describing one underlying process through shared genes, so 40 terms
   collapse into a handful of stories. Components are labelled and their sizes
   are reported in the console and in the CSV.

3. VOLCANO              <- create_enhanced_volcano()
   Per contrast, from the DEG table rather than the enrichment output.
   x = log2FC, y = -log10(FDR), cut lines at the agreed |log2FC| > 1.5 and
   FDR < 0.05, up and down counts annotated, same three colours.

DEPENDENCIES
------------
numpy, scipy, pandas, matplotlib only. networkx is used for the layout when it
is importable, and there is a self-contained Fruchterman-Reingold fallback in
CELL 3 so the network figure still works if sc_pre does not have it. Nothing
here needs installing.

OUTPUT
------
  04_figures/
    kegg_bubble/     {contrast}_KEGG_bubble.pdf/.png   + _KEGG_bubble.csv
    go_network/      {contrast}_{ont}_network.pdf/.png + _network_nodes.csv
                                                       + _network_edges.csv
    volcano/         {contrast}_{layer}_volcano.pdf/.png
    figure_manifest.csv

All text 28-34pt, saved as both PDF and PNG at 300 DPI.

Author: Jake Lehle
Date: September 2026
"""

# %% ============================================================
# CELL 1 - CONFIG
# ============================================================
import os
import re
import sys
import glob
import math
import textwrap
from collections import defaultdict

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import Normalize, LinearSegmentedColormap
from matplotlib.cm import ScalarMappable

WORKING_DIR = '/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV'
ANN_DIR = os.path.join(WORKING_DIR, 'annotation_output')
ENRICH_DIR = os.path.join(ANN_DIR, 'enrichment')
SUM_DIR = os.path.join(ENRICH_DIR, '03_summary')
FIG_DIR = os.path.join(ENRICH_DIR, '04_figures')
DEG_TAB_DIR = os.path.join(ANN_DIR, 'deg_comparisons', 'tables')

BUB_DIR = os.path.join(FIG_DIR, 'kegg_bubble')
NET_DIR = os.path.join(FIG_DIR, 'go_network')
VOL_DIR = os.path.join(FIG_DIR, 'volcano')

SIG_CSV = os.path.join(SUM_DIR, 'enrichment_significant.csv')

# ---- thresholds, kept identical to the upstream scripts ----
ENRICH_ALPHA = 0.05      # corrected p for a term to be shown
PRIMARY_ALPHA = 0.05     # volcano FDR line
PRIMARY_LFC = 1.5        # volcano fold-change lines

# ---- figure dials ----
TOP_N_BUBBLE = 15        # pathways per direction, by corrected p
TOP_N_NETWORK = 40       # terms per ontology fed to the network
JACCARD_MIN = 0.20       # edge cutoff, emapplot's default neighbourhood
MIN_TERMS_FOR_NETWORK = 5
LABEL_WRAP = 34          # characters before wrapping a term name
NET_LABEL_TOP = 12       # label only the most significant N nodes
VOLCANO_LAYER = 'pseudobulk'

ENGINE = 'gprofiler'     # which engine's results to plot
GO_SOURCES = ['GO:BP', 'GO:MF', 'GO:CC']
KEGG_SOURCE = 'KEGG'

FIG_DPI = 300
HEX = {
    'up': '#C4453C', 'down': '#2C6E9B', 'ns': '#B8B8B8',
    'grid': '#D8D8D8', 'edge': '#9AA5AD', 'text': '#222222',
}
# sequential ramps for -log10 p, one per direction, light -> saturated
CMAP = {
    'up': LinearSegmentedColormap.from_list('up', ['#F6D9D6', '#C4453C', '#7E2019']),
    'down': LinearSegmentedColormap.from_list('down', ['#D4E3EE', '#2C6E9B', '#17405C']),
}

plt.rcParams.update({
    'font.size': 28, 'axes.titlesize': 34, 'axes.labelsize': 30,
    'xtick.labelsize': 24, 'ytick.labelsize': 24, 'legend.fontsize': 22,
    'figure.titlesize': 34,
    'pdf.fonttype': 42, 'ps.fonttype': 42,
    'savefig.bbox': 'tight', 'figure.facecolor': 'white',
    'axes.facecolor': 'white',
})

for d in (FIG_DIR, BUB_DIR, NET_DIR, VOL_DIR):
    os.makedirs(d, exist_ok=True)

pd.set_option('display.width', 220)
pd.set_option('display.max_columns', 40)

_report = []
_manifest = []


def log(msg=''):
    print(msg, flush=True)
    _report.append(str(msg))


def hr(title, ch='='):
    log('')
    log(ch * 90)
    log(f'{title}')
    log(ch * 90)


def pretty(stem):
    """Contrast stem -> something readable on a figure title."""
    s = re.sub(r'^C\d+_', '', str(stem))
    s = s.replace('_pseudobulk', '').replace('_wilcoxon', '')
    s = s.replace('_', ' ').replace(' vs ', ' vs ')
    return s.strip()


def save(fig, outdir, name, kind, contrast):
    """PDF and PNG, 300 DPI, per the project convention."""
    for ext in ('pdf', 'png'):
        p = os.path.join(outdir, f'{name}.{ext}')
        fig.savefig(p, dpi=FIG_DPI, facecolor='white')
    plt.close(fig)
    _manifest.append({'kind': kind, 'contrast': contrast, 'name': name,
                      'pdf': os.path.join(outdir, f'{name}.pdf'),
                      'png': os.path.join(outdir, f'{name}.png')})
    log(f'    saved {name}.pdf / .png')


def wrap(s, n=LABEL_WRAP):
    return '\n'.join(textwrap.wrap(str(s), n)) if len(str(s)) > n else str(s)


# %% ============================================================
# CELL 2 - LOAD
# ============================================================
hr('CELL 2 - load enrichment results')

if not os.path.exists(SIG_CSV):
    sys.exit(f'[FATAL] {SIG_CSV} not found. Run SHIV_enrichment_run.py first.')

sig = pd.read_csv(SIG_CSV)
log(f'  significant terms : {len(sig):,}')
log(f'  columns           : {list(sig.columns)}')

need = {'contrast_id', 'stem', 'direction', 'source', 'term_id', 'term_name',
        'p_adj', 'overlap', 'term_size', 'query_size', 'background_size'}
missing = need - set(sig.columns)
if missing:
    sys.exit(f'[FATAL] enrichment_significant.csv is missing {sorted(missing)}')

if 'engine' in sig.columns:
    before = len(sig)
    sig = sig[sig['engine'] == ENGINE].copy()
    log(f'  engine == {ENGINE}  : {len(sig):,} of {before:,}')

sig = sig[sig['p_adj'] < ENRICH_ALPHA].copy()
sig['nlp'] = -np.log10(sig['p_adj'].clip(lower=1e-300))

# Fold enrichment: (k/N) / (n/M). Comparable between an up list and a down
# list of very different size, which GeneRatio is not.
sig['gene_ratio'] = sig['overlap'] / sig['query_size'].replace(0, np.nan)
sig['bg_ratio'] = sig['term_size'] / sig['background_size'].replace(0, np.nan)
sig['fold_enrichment'] = sig['gene_ratio'] / sig['bg_ratio']

has_genes = 'genes' in sig.columns and sig['genes'].notna().any() \
    and (sig['genes'].astype(str).str.len() > 0).any()
log(f'  hit-gene lists     : {"present" if has_genes else "ABSENT"}')

log('\n  terms by source and direction:')
log(pd.crosstab(sig['source'], sig['direction']).to_string())

contrasts = sorted(sig['stem'].unique())
log(f'\n  contrasts with at least one significant term: {len(contrasts)}')
for c in contrasts:
    sub = sig[sig['stem'] == c]
    log(f'    {c}  ({len(sub)} terms)')


# %% ============================================================
# CELL 3 - GRAPH LAYOUT (networkx if present, else self-contained)
# ============================================================
hr('CELL 3 - layout backend')

try:
    import networkx as nx
    _HAS_NX = True
    log(f'  networkx {nx.__version__} available, using spring_layout')
except Exception:
    _HAS_NX = False
    log('  networkx not available, using the built-in Fruchterman-Reingold')
    log('  fallback (numpy only). Layout differs cosmetically; the graph,')
    log('  the edges and the components are identical either way.')


def fr_layout(n_nodes, edges, weights, seed=0, iters=260):
    """Fruchterman-Reingold, enough for a 40-node term map. Deterministic.

    Repulsion between every pair, attraction along edges scaled by similarity,
    linearly cooling step size. Only needed when networkx is absent."""
    rng = np.random.default_rng(seed)
    pos = rng.normal(0, 1, (n_nodes, 2))
    if n_nodes == 1:
        return np.zeros((1, 2))
    k = 1.0 / math.sqrt(n_nodes)
    t = 0.12
    E = np.array(edges, dtype=int) if len(edges) else np.zeros((0, 2), int)
    W = np.array(weights, dtype=float) if len(weights) else np.zeros(0)
    for _ in range(iters):
        d = pos[:, None, :] - pos[None, :, :]
        dist = np.linalg.norm(d, axis=-1)
        np.fill_diagonal(dist, np.inf)
        rep = (k ** 2 / dist ** 2)[:, :, None] * d
        disp = np.nan_to_num(rep).sum(axis=1)
        if len(E):
            dv = pos[E[:, 0]] - pos[E[:, 1]]
            dl = np.linalg.norm(dv, axis=1)
            dl[dl == 0] = 1e-9
            att = ((dl / k)[:, None] * dv / dl[:, None]) * W[:, None]
            np.add.at(disp, E[:, 0], -att)
            np.add.at(disp, E[:, 1], att)
        dl = np.linalg.norm(disp, axis=1)
        dl[dl == 0] = 1e-9
        pos += (disp / dl[:, None]) * np.minimum(dl, t)[:, None]
        t *= 0.985
    span = pos.max(axis=0) - pos.min(axis=0)
    span[span == 0] = 1.0
    return (pos - pos.min(axis=0)) / span * 2 - 1


def components(n_nodes, edges):
    """Connected components by BFS. No dependency."""
    adj = defaultdict(list)
    for a, b in edges:
        adj[a].append(b)
        adj[b].append(a)
    comp = [-1] * n_nodes
    c = 0
    for s in range(n_nodes):
        if comp[s] != -1:
            continue
        stack, comp[s] = [s], c
        while stack:
            u = stack.pop()
            for v in adj[u]:
                if comp[v] == -1:
                    comp[v] = c
                    stack.append(v)
        c += 1
    return np.array(comp), c


# --- self-test: the layout and component code must be right before use ---
_e = [(0, 1), (1, 2), (3, 4)]          # two components: {0,1,2} and {3,4}, 5 alone
_c, _n = components(6, _e)
assert _n == 3, f'component count wrong: {_n}'
assert _c[0] == _c[1] == _c[2] and _c[3] == _c[4] and _c[5] not in (_c[0], _c[3])
_p = fr_layout(6, _e, [1.0] * 3)
assert _p.shape == (6, 2) and np.isfinite(_p).all(), 'layout produced bad coords'
assert np.linalg.norm(_p[0] - _p[1]) < np.linalg.norm(_p[0] - _p[5]), \
    'connected nodes should end up closer than unconnected ones'
log('  [self-test] components and layout PASS')


def jaccard_matrix(gene_sets):
    """Pairwise Jaccard index, which is what emapplot's pairwise_termsim uses
    by default. |A and B| / |A or B|."""
    n = len(gene_sets)
    J = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            a, b = gene_sets[i], gene_sets[j]
            if not a or not b:
                continue
            u = len(a | b)
            if u:
                J[i, j] = J[j, i] = len(a & b) / u
    return J


# hand-checkable: {1,2,3} vs {2,3,4} -> 2 shared, 4 union -> 0.5
_J = jaccard_matrix([{1, 2, 3}, {2, 3, 4}, {9}])
assert abs(_J[0, 1] - 0.5) < 1e-12, f'jaccard wrong: {_J[0,1]}'
assert _J[0, 2] == 0.0
log('  [self-test] jaccard PASS')


# --- readability machinery ------------------------------------------------
# A raw spring layout of enriched GO terms is close to unusable: terms that
# share genes pile on top of each other, and the isolated terms, of which
# there are always many, get flung into a ring around the edge that reads as
# structure when it is nothing of the sort. Three passes fix that.

def separate(P, r, iters=420, pad=7.0):
    """Push overlapping nodes apart. P and r are both in POINTS, so the
    result respects the drawn marker sizes rather than arbitrary data units."""
    P = P.astype(float).copy()
    n = len(P)
    if n < 2:
        return P
    for _ in range(iters):
        moved = False
        d = P[:, None, :] - P[None, :, :]
        dist = np.linalg.norm(d, axis=-1)
        need = r[:, None] + r[None, :] + pad
        np.fill_diagonal(dist, np.inf)
        bad = dist < need
        if not bad.any():
            break
        for i, j in zip(*np.where(np.triu(bad, 1))):
            v = P[i] - P[j]
            dl = np.linalg.norm(v)
            if dl < 1e-9:
                v = np.random.default_rng(i * 97 + j).normal(0, 1, 2)
                dl = np.linalg.norm(v)
            push = (need[i, j] - dl) / 2.0
            P[i] += v / dl * push
            P[j] -= v / dl * push
            moved = True
        if not moved:
            break
    return P


def grid_positions(n, x0, x1, y0, y1, per_row):
    """Tidy grid for the isolated terms, kept out of the main graph area."""
    if n == 0:
        return np.zeros((0, 2))
    rows = int(math.ceil(n / per_row))
    xs = np.linspace(x0, x1, per_row)
    ys = np.linspace(y1, y0, rows) if rows > 1 else np.array([(y0 + y1) / 2])
    out = [(xs[i % per_row], ys[i // per_row]) for i in range(n)]
    return np.array(out)


def place_labels(idx, P, r, texts, fs, ax_w, ax_h, y_floor=0.0):
    """Greedy non-overlapping label placement, all in POINTS.

    For each label, try eight positions around its node and keep the one that
    collides least with the nodes and with the labels already placed. Returns
    offsets in points, which is what annotate(textcoords='offset points')
    wants, so nothing has to be converted back."""
    boxes = []
    offs = {}
    for i in idx:
        lines = str(texts[i]).split('\n')
        w = max(len(l) for l in lines) * fs * 0.58
        h = len(lines) * fs * 1.28
        best, best_cost = None, np.inf
        for ang in np.deg2rad([90, 45, 135, 0, 180, 300, 240, 270]):
            for mult in (1.0, 1.5, 2.0, 2.6, 3.3, 4.1):
                rad = r[i] + 10 + (h / 2) * mult
                cx = P[i, 0] + math.cos(ang) * rad
                cy = P[i, 1] + math.sin(ang) * rad
                bx = (cx - w / 2, cy - h / 2, cx + w / 2, cy + h / 2)
                cost = 0.0
                # off-canvas is the worst outcome, weight it heavily
                if bx[0] < 0 or bx[1] < 0 or bx[2] > ax_w or bx[3] > ax_h:
                    cost += 900
                # never let a cluster label drop into the isolated-terms strip
                if bx[1] < y_floor:
                    cost += 900
                for ob in boxes:
                    ox = max(0, min(bx[2], ob[2]) - max(bx[0], ob[0]))
                    oy = max(0, min(bx[3], ob[3]) - max(bx[1], ob[1]))
                    cost += ox * oy / 12.0
                for j in range(len(P)):
                    nx_, ny_ = P[j]
                    ddx = max(abs(nx_ - cx) - w / 2, 0)
                    ddy = max(abs(ny_ - cy) - h / 2, 0)
                    if math.hypot(ddx, ddy) < r[j]:
                        cost += 60
                cost += 0.22 * rad          # prefer staying near the node
                if cost < best_cost:
                    best_cost, best = cost, (cx, cy, bx)
        cx, cy, bx = best
        boxes.append(bx)
        offs[i] = (cx - P[i, 0], cy - P[i, 1])
    return offs


_P = separate(np.array([[0.0, 0.0], [1.0, 0.0]]), np.array([20.0, 20.0]))
assert np.linalg.norm(_P[0] - _P[1]) >= 40.0, 'separate() did not resolve overlap'
_g = grid_positions(5, 0, 100, 0, 50, 3)
assert _g.shape == (5, 2) and len(np.unique(_g, axis=0)) == 5, 'grid collided'
log('  [self-test] node separation and isolate grid PASS')


# %% ============================================================
# CELL 4 - KEGG BUBBLE PLOTS
# ============================================================
hr('CELL 4 - KEGG bubble plots')

kegg = sig[sig['source'] == KEGG_SOURCE].copy()
log(f'  significant KEGG terms: {len(kegg)} across '
    f'{kegg["stem"].nunique()} contrast(s)')

if kegg.empty:
    log('  none, skipping')
else:
    for stem, sub in kegg.groupby('stem'):
        cid = str(sub['contrast_id'].iloc[0])
        log(f'\n  {stem}')
        panels = []
        for direction in ('up', 'down'):
            d = sub[sub['direction'] == direction] \
                .nsmallest(TOP_N_BUBBLE, 'p_adj') \
                .sort_values('fold_enrichment')
            panels.append((direction, d))
            log(f'    {direction:4s}: {len(d)} pathway(s) plotted '
                f'of {int((sub["direction"] == direction).sum())}')

        if all(len(d) == 0 for _, d in panels):
            continue

        n_rows = max(1, max(len(d) for _, d in panels))
        # One bubble scale across BOTH panels, so a pathway with 12 genes looks
        # the same size whichever direction it came from, and one shared size
        # legend is then honest. Per-panel scaling would silently make a
        # 3-gene down pathway as big as a 40-gene up one.
        omax = int(max(1, sub['overlap'].max()))

        def bubble_area(counts):
            return 200 + 2800 * (np.asarray(counts, float) / omax)

        # narrow panels that carry little, instead of giving an empty
        # direction half the canvas
        ratios = [max(1.0, len(d)) for _, d in panels]
        fig, axes = plt.subplots(
            1, 2, figsize=(30, max(11, 1.15 * n_rows + 5.5)),
            gridspec_kw={'width_ratios': ratios})

        for ax, (direction, d) in zip(axes, panels):
            colour = HEX[direction]
            n_tot = int((sub['direction'] == direction).sum())
            if len(d) == 0:
                ax.text(0.5, 0.5, f'no significant\n{direction} pathways',
                        ha='center', va='center', fontsize=28,
                        color=HEX['ns'], transform=ax.transAxes)
                ax.set_xticks([])
                ax.set_yticks([])
                for sp in ax.spines.values():
                    sp.set_visible(False)
                ax.set_title(f'{direction.capitalize()}regulated',
                             fontsize=32, color=colour, pad=16)
                continue

            y = np.arange(len(d))
            # A single term gives a degenerate colour range, which matplotlib
            # renders as an unreadable offset like "1e-9+3.26". Anchor the
            # scale at 0 instead so the bar still means something.
            lo, hi = float(d['nlp'].min()), float(d['nlp'].max())
            if hi - lo < 1e-9:
                lo = 0.0
                hi = max(hi, 1e-9)
            sc = ax.scatter(d['fold_enrichment'], y, s=bubble_area(d['overlap']),
                            c=d['nlp'], cmap=CMAP[direction],
                            norm=Normalize(vmin=lo, vmax=hi),
                            edgecolors='#333333', linewidths=1.2, zorder=3)
            ax.set_yticks(y)
            ax.set_yticklabels([wrap(t) for t in d['term_name']], fontsize=22)
            ax.set_xlabel('Fold enrichment', fontsize=30)
            ax.set_title(f'{direction.capitalize()}regulated  (n = {n_tot})',
                         fontsize=32, color=colour, pad=16)
            ax.grid(axis='x', color=HEX['grid'], lw=1.0, zorder=0)
            ax.set_axisbelow(True)
            ax.axvline(1.0, color='#777777', ls='--', lw=1.6, zorder=1)
            ax.margins(x=0.22, y=0.10 if len(d) > 2 else 0.35)
            for sp in ('top', 'right'):
                ax.spines[sp].set_visible(False)

            cb = fig.colorbar(sc, ax=ax, pad=0.02, fraction=0.045)
            cb.set_label(r'$-\log_{10}$ adj. $p$', fontsize=24)
            cb.ax.tick_params(labelsize=20)
            cb.formatter.set_useOffset(False)
            cb.update_ticks()

        # shared size legend, outside the panels so it never sits on data
        steps = sorted({max(1, omax // 4), max(1, omax // 2), omax})
        handles = [Line2D([], [], marker='o', ls='', markerfacecolor='white',
                          markeredgecolor='#333333',
                          markersize=np.sqrt(bubble_area(sv)) / 3.0,
                          label=str(sv)) for sv in steps]
        fig.legend(handles=handles, title='Genes in pathway',
                   loc='lower center', bbox_to_anchor=(0.5, -0.13), ncol=len(steps),
                   frameon=True, fontsize=22, title_fontsize=24,
                   handletextpad=1.4, columnspacing=3.5, borderpad=1.0)

        fig.suptitle(f'KEGG pathways  |  {pretty(stem)}', fontsize=34, y=1.005)
        fig.tight_layout()
        save(fig, BUB_DIR, f'{cid}_{pretty(stem).replace(" ", "_")}_KEGG_bubble',
             'kegg_bubble', stem)

        out = sub[['contrast_id', 'stem', 'direction', 'term_id', 'term_name',
                   'p_adj', 'overlap', 'term_size', 'query_size',
                   'background_size', 'gene_ratio', 'fold_enrichment']]
        out.sort_values(['direction', 'p_adj']).to_csv(
            os.path.join(BUB_DIR,
                         f'{cid}_{pretty(stem).replace(" ", "_")}_KEGG_bubble.csv'),
            index=False)


# %% ============================================================
# CELL 5 - GO TERM SIMILARITY NETWORKS
# ============================================================
hr('CELL 5 - GO term networks (emapplot equivalent)')

if not has_genes:
    log('  [SKIP] no hit-gene lists in enrichment_significant.csv, so term')
    log('  similarity cannot be computed. This figure is built on real gene')
    log('  overlap; it will not be faked from p values.')
else:
    go = sig[sig['source'].isin(GO_SOURCES)].copy()
    log(f'  significant GO terms: {len(go)}')

    for (stem, source, direction), sub in go.groupby(
            ['stem', 'source', 'direction']):
        if len(sub) < MIN_TERMS_FOR_NETWORK:
            continue
        cid = str(sub['contrast_id'].iloc[0])
        ont = source.replace('GO:', '')
        sub = sub.nsmallest(TOP_N_NETWORK, 'p_adj').reset_index(drop=True)

        gsets = [set(str(g).split(';')) - {''} for g in sub['genes']]
        if sum(len(g) for g in gsets) == 0:
            continue

        J = jaccard_matrix(gsets)
        iu = np.triu_indices(len(sub), 1)
        edges = [(int(a), int(b)) for a, b in zip(*iu) if J[a, b] >= JACCARD_MIN]
        weights = [float(J[a, b]) for a, b in edges]

        comp, n_comp = components(len(sub), edges)
        sizes = pd.Series(comp).value_counts()
        n_iso = int((sizes == 1).sum())

        log(f'\n  {stem}  {source}  {direction}')
        log(f'    terms {len(sub)}, edges {len(edges)} at JC >= {JACCARD_MIN}, '
            f'components {n_comp} ({n_iso} isolated)')
        if len(sizes[sizes > 1]):
            log(f'    clustered component sizes: '
                f'{sorted(sizes[sizes > 1].tolist(), reverse=True)}')

        # --- split connected terms from isolated ones ---
        # Isolated terms share no genes with anything else here. Laying them
        # out together with the graph throws them into a ring that looks like
        # structure; they belong in their own strip, clearly separated.
        csize = np.array([sizes[c] for c in comp])
        conn = np.where(csize > 1)[0]
        iso = np.where(csize == 1)[0]

        FIG_W, FIG_H = 24.0, 19.0
        fig = plt.figure(figsize=(FIG_W, FIG_H))
        ax = fig.add_axes([0.03, 0.04, 0.80, 0.84])
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        AX_W = 0.80 * FIG_W * 72.0          # axes size in POINTS
        AX_H = 0.84 * FIG_H * 72.0

        area = (400 + 5200 * (sub['overlap'] /
                              max(1, sub['overlap'].max()))).to_numpy()
        rad = np.sqrt(area / math.pi)        # marker radius in points

        # reserve the bottom strip for isolates when there are any
        iso_rows = 0 if len(iso) == 0 else math.ceil(len(iso) / min(9, len(iso)))
        iso_frac = 0.0 if len(iso) == 0 else min(0.36, 0.085 * iso_rows + 0.085)
        main_y0 = iso_frac

        P = np.zeros((len(sub), 2))
        if len(conn):
            sub_e = [(int(np.where(conn == a)[0][0]),
                      int(np.where(conn == b)[0][0])) for a, b in edges]
            if _HAS_NX:
                G = nx.Graph()
                G.add_nodes_from(range(len(conn)))
                for (a, b), w in zip(sub_e, weights):
                    G.add_edge(a, b, weight=w)
                p = nx.spring_layout(G, weight='weight', seed=1,
                                     k=2.4 / math.sqrt(max(1, len(conn))),
                                     iterations=400)
                lay = np.array([p[i] for i in range(len(conn))])
            else:
                lay = fr_layout(len(conn), sub_e, weights, seed=1)
            span = lay.max(axis=0) - lay.min(axis=0)
            span[span == 0] = 1.0
            lay = (lay - lay.min(axis=0)) / span
            lay[:, 0] = 0.06 + lay[:, 0] * 0.88
            lay[:, 1] = main_y0 + 0.06 + lay[:, 1] * (0.94 - main_y0 - 0.08)
            P[conn] = lay

        if len(iso):
            per_row = min(9, len(iso))
            P[iso] = grid_positions(len(iso), 0.07, 0.93,
                                    0.055, max(0.065, iso_frac - 0.085),
                                    per_row)

        # de-overlap, in points, then back to axes fraction
        Ppt = P * np.array([AX_W, AX_H])
        if len(conn):
            moved = separate(Ppt[conn], rad[conn])
            Ppt[conn] = moved
        P = Ppt / np.array([AX_W, AX_H])

        for (a, b), w in zip(edges, weights):
            ax.plot([P[a, 0], P[b, 0]], [P[a, 1], P[b, 1]],
                    color=HEX['edge'], lw=1.6 + 9.0 * w,
                    alpha=0.40 + 0.45 * w, zorder=1, solid_capstyle='round')

        norm = Normalize(vmin=sub['nlp'].min(),
                         vmax=max(sub['nlp'].max(), sub['nlp'].min() + 1e-9))
        sc = ax.scatter(P[:, 0], P[:, 1], s=area, c=sub['nlp'],
                        cmap=CMAP[direction], norm=norm,
                        edgecolors='#2B2B2B', linewidths=1.8, zorder=3)

        if len(iso):
            ax.axhline(iso_frac, color=HEX['grid'], lw=2.0, ls='--', zorder=0)
            ax.text(0.5, iso_frac - 0.014, f'{len(iso)} term(s) sharing no '
                    f'genes with any other at JC >= {JACCARD_MIN}',
                    ha='center', va='top', fontsize=20, color='#666666',
                    transform=ax.transAxes)

        # Label the most significant terms, preferring ones inside a cluster,
        # since an isolated term explains itself and a cluster is what the
        # reader actually needs named.
        order = sub['nlp'].to_numpy().argsort()[::-1]
        pref = [i for i in order if csize[i] > 1] + \
               [i for i in order if csize[i] == 1]
        lab_idx = pref[:min(NET_LABEL_TOP, len(sub))]
        texts = [wrap(t, 24) for t in sub['term_name']]
        offs = place_labels(lab_idx, Ppt, rad, texts, 19, AX_W, AX_H,
                            y_floor=iso_frac * AX_H + 14)
        for i in lab_idx:
            ax.annotate(texts[i], (P[i, 0], P[i, 1]),
                        textcoords='offset points', xytext=offs[i],
                        ha='center', va='center', fontsize=19,
                        color=HEX['text'], zorder=5,
                        bbox=dict(boxstyle='round,pad=0.30', fc='white',
                                  ec='#BBBBBB', alpha=0.90),
                        arrowprops=dict(arrowstyle='-', color='#AAAAAA',
                                        lw=1.2, shrinkA=0, shrinkB=2))

        clustered = sorted(sizes[sizes > 1].tolist(), reverse=True)
        ax.set_title(f'GO {ont} term network  |  {pretty(stem)}  |  '
                     f'{direction}regulated\n'
                     f'{len(sub)} terms, edge = Jaccard >= {JACCARD_MIN}, '
                     f'{len(clustered)} cluster(s) '
                     f'{clustered if clustered else ""}',
                     fontsize=29, pad=24)
        ax.set_xticks([])
        ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(False)

        cax = fig.add_axes([0.855, 0.10, 0.022, 0.62])
        cb = fig.colorbar(sc, cax=cax)
        cb.set_label(r'$-\log_{10}$ adj. $p$', fontsize=26)
        cb.ax.tick_params(labelsize=22)

        omax = int(sub['overlap'].max())
        steps = sorted({max(1, omax // 4), max(1, omax // 2), omax})
        handles = [Line2D([], [], marker='o', ls='', markerfacecolor='white',
                          markeredgecolor='#2B2B2B',
                          markersize=np.sqrt(400 + 5200 * (s / max(1, omax))) / 3.2,
                          label=str(s)) for s in steps]
        fig.legend(handles=handles, title='Genes in term',
                   loc='upper right', bbox_to_anchor=(0.995, 0.93),
                   frameon=True, fontsize=21, title_fontsize=23,
                   labelspacing=2.0, borderpad=1.0)

        base = f'{cid}_{pretty(stem).replace(" ", "_")}_{ont}_{direction}_network'
        save(fig, NET_DIR, base, 'go_network', stem)

        nodes_out = sub[['term_id', 'term_name', 'p_adj', 'overlap',
                         'term_size', 'genes']].copy()
        nodes_out['component'] = comp
        nodes_out['component_size'] = [int(sizes[c]) for c in comp]
        nodes_out.sort_values(['component_size', 'component', 'p_adj'],
                              ascending=[False, True, True]).to_csv(
            os.path.join(NET_DIR, f'{base}_nodes.csv'), index=False)
        pd.DataFrame({
            'term_a': [sub.loc[a, 'term_id'] for a, _ in edges],
            'term_a_name': [sub.loc[a, 'term_name'] for a, _ in edges],
            'term_b': [sub.loc[b, 'term_id'] for _, b in edges],
            'term_b_name': [sub.loc[b, 'term_name'] for _, b in edges],
            'jaccard': weights,
        }).sort_values('jaccard', ascending=False).to_csv(
            os.path.join(NET_DIR, f'{base}_edges.csv'), index=False)


# %% ============================================================
# CELL 6 - VOLCANO PLOTS
# ============================================================
hr('CELL 6 - volcano plots')

# BUG FIXED 2026-09-28. The first version looked for "{stem}.csv" and
# "{stem}_{layer}.csv" and missed every table.
#
# The reason: enrichment_significant.csv does NOT carry the full table stem.
# The prep script parses the ORA filename with
#     ^(C\d+)_(.+)_{LAYER}_{GATE}_(up|down)\.csv$
# so the contrast id goes into its own column and `stem` is only the middle,
# e.g. contrast_id "C00061" + stem "ALL_Both_PBMC_vs_LN_Necropsy". The DEG
# table on disk is the three parts joined back up:
#     C00061_ALL_Both_PBMC_vs_LN_Necropsy_pseudobulk.csv
# I rebuilt the name from `stem` alone, so the C-prefix was missing and all
# nine lookups failed. Same class of mistake as the entrez bug: an identifier
# reassembled by assumption instead of from the column that holds it.
#
# Now the id comes from contrast_id, several spellings are tried, and if
# nothing resolves the script PRINTS what is actually in the directory rather
# than reporting nine identical misses that do not say why.

# 'enriched' = contrasts with at least one significant term (what the figures
#              above cover)
# 'listed'   = every contrast that produced an ORA gene list at this gate and
#              layer, so contrasts that passed the DEG threshold but enriched
#              to nothing still get a volcano. This is the default: an empty
#              enrichment is worth seeing the volcano for.
# 'all'      = every DEG table at this layer (about 97 figures)
VOLCANO_SCOPE = 'listed'
ORA_LIST_DIR = os.path.join(ANN_DIR, 'deg_comparisons', 'ora_gene_lists')

if not os.path.isdir(DEG_TAB_DIR):
    log(f'  [SKIP] {DEG_TAB_DIR} not found')
else:
    available = sorted(f for f in os.listdir(DEG_TAB_DIR) if f.endswith('.csv'))

    targets = []          # (contrast_id, stem)
    if VOLCANO_SCOPE == 'enriched':
        targets = sorted(set(zip(sig['contrast_id'].astype(str), sig['stem'])))
    elif VOLCANO_SCOPE == 'listed' and os.path.isdir(ORA_LIST_DIR):
        pat = re.compile(rf'^(C\d+)_(.+)_{VOLCANO_LAYER}_(?:FDR|rawP)_'
                         r'(?:up|down)\.csv$')
        seen = set()
        for fn in sorted(os.listdir(ORA_LIST_DIR)):
            m = pat.match(fn)
            if m and m.groups() not in seen:
                seen.add(m.groups())
                targets.append((m.group(1), m.group(2)))
        targets = sorted(set(targets))
    if VOLCANO_SCOPE == 'all' or not targets:
        pat = re.compile(rf'^(C\d+)_(.+)_{VOLCANO_LAYER}\.csv$')
        targets = sorted({(m.group(1), m.group(2)) for m in
                          (pat.match(f) for f in available) if m})

    log(f'  scope = {VOLCANO_SCOPE}, layer = {VOLCANO_LAYER}, '
        f'{len(targets)} contrast(s)')

    n_ok = n_miss = 0
    for cid, stem in targets:
        # rebuild the table name from the id AND the stem, not the stem alone
        cands = [f'{cid}_{stem}_{VOLCANO_LAYER}.csv',
                 f'{cid}_{stem}.csv',
                 f'{stem}_{VOLCANO_LAYER}.csv']
        path = next((os.path.join(DEG_TAB_DIR, c) for c in cands
                     if os.path.exists(os.path.join(DEG_TAB_DIR, c))), None)
        if path is None:
            hits = [f for f in available
                    if f.startswith(cid + '_') and VOLCANO_LAYER in f]
            if hits:
                path = os.path.join(DEG_TAB_DIR, hits[0])
        if path is None:
            n_miss += 1
            if n_miss == 1:
                log(f'    [MISS] {cid} / {stem}')
                log(f'           tried: {cands}')
                log(f'           directory holds e.g.: {available[:3]}')
            continue

        df = pd.read_csv(path)
        need_cols = {'log2FC', 'pval', 'padj'}
        if not need_cols <= set(df.columns):
            log(f'    [SKIP] {cid}: table lacks {sorted(need_cols - set(df.columns))}')
            continue
        df = df.dropna(subset=['pval', 'padj', 'log2FC'])
        if df.empty:
            continue

        df['nlq'] = -np.log10(df['padj'].clip(lower=1e-300))
        up = (df['padj'] < PRIMARY_ALPHA) & (df['log2FC'] > PRIMARY_LFC)
        dn = (df['padj'] < PRIMARY_ALPHA) & (df['log2FC'] < -PRIMARY_LFC)
        ns = ~(up | dn)

        fig, ax = plt.subplots(figsize=(16, 15))
        ax.scatter(df.loc[ns, 'log2FC'], df.loc[ns, 'nlq'], s=34,
                   c=HEX['ns'], alpha=0.45, linewidths=0, zorder=1)
        ax.scatter(df.loc[dn, 'log2FC'], df.loc[dn, 'nlq'], s=62,
                   c=HEX['down'], alpha=0.90, linewidths=0, zorder=2)
        ax.scatter(df.loc[up, 'log2FC'], df.loc[up, 'nlq'], s=62,
                   c=HEX['up'], alpha=0.90, linewidths=0, zorder=2)

        ax.axvline(-PRIMARY_LFC, ls='--', color='#444444', lw=1.8)
        ax.axvline(PRIMARY_LFC, ls='--', color='#444444', lw=1.8)
        ax.axhline(-np.log10(PRIMARY_ALPHA), ls='--', color='#444444', lw=1.8)

        ax.set_xlabel(r'$\log_2$(fold change)', fontsize=32)
        ax.set_ylabel(r'$-\log_{10}$(FDR)', fontsize=32)
        ax.set_title(f'{cid}  {pretty(stem)}\n{VOLCANO_LAYER}, '
                     rf'|$\log_2$FC| > {PRIMARY_LFC}, FDR < {PRIMARY_ALPHA}',
                     fontsize=30, pad=20)
        ax.grid(color=HEX['grid'], lw=0.9, zorder=0)
        ax.set_axisbelow(True)
        for sp in ('top', 'right'):
            ax.spines[sp].set_visible(False)

        ax.legend(handles=[
            Line2D([], [], marker='o', ls='', color=HEX['up'], markersize=18,
                   label=f'Up  ({int(up.sum())})'),
            Line2D([], [], marker='o', ls='', color=HEX['down'], markersize=18,
                   label=f'Down  ({int(dn.sum())})'),
            Line2D([], [], marker='o', ls='', color=HEX['ns'], markersize=18,
                   label=f'Not significant  ({int(ns.sum())})'),
        ], loc='upper left', frameon=True, fontsize=24, borderpad=0.9)

        fig.tight_layout()
        save(fig, VOL_DIR,
             f'{cid}_{stem}_{VOLCANO_LAYER}_volcano', 'volcano', stem)
        log(f'    {cid} {pretty(stem)}: up {int(up.sum())}, '
            f'down {int(dn.sum())}')
        n_ok += 1

    log(f'\n  volcanoes written: {n_ok}, missing tables: {n_miss}')
    if n_ok == 0 and targets:
        log('  [ERROR] no DEG table resolved for any contrast. The directory')
        log(f'  listing above shows the real filenames in {DEG_TAB_DIR}.')



# %% ============================================================
# CELL 7 - MANIFEST
# ============================================================
hr('CELL 7 - done')

man = pd.DataFrame(_manifest)
if len(man):
    man.to_csv(os.path.join(FIG_DIR, 'figure_manifest.csv'), index=False)
    log('\n  figures written:')
    log(man['kind'].value_counts().to_string())
    log(f'\n  total: {len(man)} figure(s), each as PDF and PNG at {FIG_DPI} DPI')
else:
    log('  no figures were produced')

log(f'\n  kegg bubble : {BUB_DIR}')
log(f'  go network  : {NET_DIR}')
log(f'  volcano     : {VOL_DIR}')

rp = os.path.join(ENRICH_DIR, 'logs', 'enrichment_figures_report.txt')
os.makedirs(os.path.dirname(rp), exist_ok=True)
with open(rp, 'w') as fh:
    fh.write('\n'.join(_report))
print(f'\n  report -> {rp}')
