#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Rhesus SHIV - Tier 2 Myeloid Diagnostic v2 (purity-gated, collinearity-aware)
==============================================================================
v2 exists because the v1 run (2026-09-25, 11,794 cells, 16 subclusters) was NOT
annotatable. Four separate problems, all addressed here. Read before rerunning.

PROBLEM 1 - 62% of the "myeloid" subset was not myeloid biology.
  Reading every v1 subcluster against its own top markers:
    cluster 6  (703)   CD3E/CD247/NKG7/GNLY/GZMB/SKAP1     -> T/NK
    cluster 7  (671)   MS4A1/CD79A/CD79B/BANK1/EBF1/FCRL1  -> B cells
    cluster 9  (429)   SKAP1/CAMK4/PRKCH/ARHGAP15          -> T cells
    cluster 12 (337)   NKG7/GZMB/GNLY/CCL5/CST7            -> NK/CTL
    cluster 0  (1615)  90% animal 34315, 3x KEG06_ mito + CAMP in top15
    cluster 8  (659)   RPS/RPL top-to-bottom               -> ambient
    cluster 11 (339)   95% animal 41861, BTG2/ZFP36        -> dissociation stress
    clusters 1,4 (2590) housekeeping-flat, nothing above |z| 1.1
  Tier-1 winner-take-all scoring leaks lymphocytes into the myeloid call. FIX:
  a hard lineage-purity gate (CELL 3b) that excludes CD3 / CD79 / NK-positive
  cells BEFORE subclustering, mirroring the CD3 hard gate already established
  for the T lineage, plus mito / ribosomal / IEG QC signatures so artifact
  clusters self-identify instead of being labelled.

PROBLEM 2 - the panels are collinear; some cannot separate even in principle.
  Across the 16 v1 subclusters, mean pairwise signature r = 0.41 and PC1 alone
  explained 49% of the variance. Worst offenders:
    Tissue Mac     vs CD206+ Mac  r = 0.99  (share C1QA,C1QB,C1QC,LGMN,MSR1)
    Intermediate   vs CD163+ Mac  r = 0.95  (share CTSD,LILRB1,MS4A7)
    Classical Mono vs Inflam Mac  r = 0.95  (share CTSD,S100A8,S100A9)
  After dropouts, Non-classical Monocytes resolved to MS4A7/LST1/LILRB1 and ALL
  THREE are also in CD163+ Macrophage, making it a strict subset of another
  panel. FIX: CELL 6b computes and prints the collinearity report every run so
  a pair that cannot be separated is stated rather than silently resolved by
  argmax. This is a panel-design finding for Zoey; code cannot fix it.

PROBLEM 3 - the v1 provisional call was argmax with no floor.
  Clusters 8, 9, 12 and 15 were called CD206+ Macrophage at z_top of -0.28,
  -0.36, -0.33 and -0.27: the argmax of an all-negative row. That is 1,573
  cells that would have entered the annotation as macrophages on no evidence.
  FIX: MIN_Z_FOR_CALL floor plus the margin gate. Anything failing either is
  reported as "no call" and stays unlabelled.

PROBLEM 4 - two real populations are absent from the panel, one symbol set wrong.
    v1 cluster 14 (164)  TCF4/LILRA4/IRF8/IRF7/GZMB/JCHAIN/BLNK -> pDC,
                         mislabelled "IFN-Stimulated Myeloid"
    v1 cluster 15 (148)  LAMP3/FSCN1/CD83/BIRC3/LY75/FLT3       -> LAMP3+
                         mature / migratory DC, mislabelled "CD206+ Macrophage"
  FIX: pDC and Migratory DC panels added, marked as additions beyond Zoey list.
    Antigen Presenting Myeloid scored 0/4 because it asked for HLA-. v1 cluster
  10 top markers are CD74/MAMU-DRA/MAMU-DQA1/MAMU-DQB1/MAMU-DRB1/MAMU-DRB5/
  MAMU-DPA, so the population exists and the symbols were the problem. FIX:
  USE_RHESUS_ALIASES adds the MAMU- forms. Zoey original strings are kept in
  the file so every substitution is visible and reversible in one dial.
    v1 cluster 5 top markers include bare FCGR3, settling FCGR3A. Same mechanism.

Still open: FCN1, SPP1 and IFITM3 all missed. The v1 resolver only does prefix
matching so it cannot find an ortholog filed under a LOC identifier. CELL 2 now
also searches every .var column for the symbol before giving up.

Read-only w.r.t. the canonical object. Writes NO labels. Spyder cells (# %%).

Author: Jake Lehle / Kaushal Lab
Date: September 2026
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
MERGED      = os.path.join(WORKING_DIR, 'shiv_starsolo', 'merged')

# Input: the Tier-2 annotated object if it exists (carries the fine T labels),
# else the TCR-recovered object. Myeloid annotation does not depend on the T
# labels, but chaining off the annotated object keeps one object moving forward.
ADATA_CANDIDATES = [
    os.path.join(MERGED, 'shiv_host_tier2_annotated.h5ad'),
    os.path.join(MERGED, 'shiv_host_final_tcrTrecovered.h5ad'),
]

TAG       = 'Myeloid'
VERSION   = 'v2'
OUT_DIR   = os.path.join(WORKING_DIR, 'annotation_output', f'tier2_{TAG}_{VERSION}')
ADATA_OUT = os.path.join(MERGED, f'shiv_{TAG}_subset_reclustered_{VERSION}.h5ad')

# ---- which Tier-1 lineages make up the myeloid compartment ----
# Unlike the T cell diagnostic this is a UNION of three Tier-1 labels, because
# Zoey's subsets cross the monocyte / macrophage / DC boundary (Monocyte-like DC
# and the macrophage subsets in particular).
MYELOID_TIER1 = ['Monocytes', 'Macrophages', 'DCs']
PROFILE_UNASSIGNED = True    # report how many Unassigned look myeloid

# ---- recluster dials ----
N_HVG        = 2000
N_PCS        = 50
N_NEIGH_PCS  = 30
SUB_RES      = 1.0
BATCH_KEY    = 'library'     # GEM-well technical variation; animal is biology
RANDOM_STATE = 0
COUNTS_LAYER = 'counts'
DO_PREBBKNN  = True

# ---- columns ----
CELLTYPE   = 'tier1_celltype'
ANIMAL_COL = 'animal_id'
TIME_COL   = 'timepoint'
TISSUE_COL = 'tissue'
LIB_COL    = 'library'
ELITE      = '40707'
INFECTED_ANIMALS = ['39272', '40702', '40707', '41861', '41862']
CONTROL_ANIMALS  = ['34315', '41903']
EXCLUDE_ANIMALS  = ['Unknown']

# Zoey's original cDC2 panel included HLA-DRA and HLA-DPA1. Setting this True
# restores them exactly as she sent them, at the cost of making cDC2 and
# Antigen Presenting Myeloid inseparable. Default False; CELL 6b reports the
# consequence either way so the choice stays visible.
CDC2_KEEP_CLASS_II = False

# ---- Zoey's Tier-2 myeloid panels, VERBATIM apart from her own 2026-09-25
# ---- corrections (IFITM8->IFITM2, APO1->APOC1) and the cDC2 class II move ----
# Do not "fix" symbols here. The resolver reports; we edit after she confirms.
MYELOID_PANELS = {
    'Classical Monocytes':       ['LYZ', 'S100A8', 'S100A9', 'S100A12', 'CTSD',
                                  'CTSS', 'FCN1', 'VCAN', 'CTSD'],
    'Intermediate Monocytes':    ['LYZ', 'CTSS', 'LILRB1', 'FCGR3A', 'CTSD', 'MS4A7'],
    # IFITM8 -> IFITM2 confirmed by Zoey 2026-09-25. IFITM3 also missed in the
    # v1 run (only IFITM10 present); left as requested pending a LOC lookup.
    'Non-classical Monocytes':   ['FCGR3A', 'MS4A7', 'LST1', 'IFITM3', 'IFITM2', 'LILRB1'],
    # Class II REMOVED from cDC2 (see CDC2_KEEP_CLASS_II below). With the MAMU-
    # aliases in place, cDC2 and Antigen Presenting Myeloid correlated at
    # r = 0.985 sharing 5 MHC-II genes, and the test-run cluster that should
    # have been cDC2 came back NO CALL as "tied with Antigen_Presenting_Myeloid".
    # CD1C / FCER1A / CLEC10A / CD1E are what distinguish a cDC2 from a generic
    # MHC-II-high cell, so cDC2 keeps those and the APC panel owns class II.
    'cDC2':                      (['CD1C', 'FCER1A', 'CLEC10A', 'CD1D', 'CST3', 'CD1E']
                                  + (['HLA-DRA', 'HLA-DPA1'] if CDC2_KEEP_CLASS_II else [])),
    'Monocyte-like DC':          ['CD1C', 'FCGR3A', 'CD14', 'CD163', 'CD36',
                                  'CD1E', 'LILRB1'],
    'Tissue Macrophage':         ['CD68', 'APOC1', 'APOE', 'LGMN', 'CTSB',
                                  'MSR1', 'C1QC', 'C1QB', 'C1QA'],
    'CD163+ Macrophage':         ['CD163', 'LST1', 'MS4A7', 'LILRB1', 'CTSD',
                                  'FCER1G', 'TYROBP'],
    'CD206+ Macrophage':         ['MRC1', 'MSR1', 'C1QC', 'C1QB', 'C1QA',
                                  'APOC1', 'LGMN'],   # APO1 -> APOC1, Zoey 2026-09-25
    'Inflammatory Macrophage':   ['S100A8', 'S100A9', 'IL1B', 'CTSD', 'CTSB',
                                  'SPP1', 'LGALS3', 'LILRB1'],
    'IFN-Stimulated Myeloid':    ['ISG15', 'IFIT1', 'IFIT2', 'IFIT3', 'MX1',
                                  'OAS1', 'OASL', 'IFI6'],
    'Antigen Presenting Myeloid': ['HLA-DRA', 'HLA-DPA1', 'HLA-DPB1', 'HLA-DQA1'],
}

# ---- panels ADDED beyond Zoey's list (v2) ----
# v1 cluster 14 and 15 are real, well-described myeloid populations that her
# eleven-subset list has no entry for, so the scorer had to put them somewhere
# and put them somewhere wrong. Scored separately and clearly marked as ours.
EXTRA_PANELS = {
    'pDC':                       ['TCF4', 'LILRA4', 'IRF8', 'IRF7', 'JCHAIN',
                                  'BLNK', 'CCDC50', 'PLD4', 'SCT'],
    'Migratory/Mature DC':       ['LAMP3', 'FSCN1', 'CD83', 'BIRC3', 'LY75',
                                  'CCR7', 'MARCKSL1', 'CTSH'],
    'cDC1':                      ['CLEC9A', 'XCR1', 'BATF3', 'CADM1'],
}
INCLUDE_EXTRA_PANELS = True

# ---- QC signatures, so artifact clusters announce themselves ----
# v1 clusters 0, 8 and 11 were mito-high, ribosomal and IEG-high respectively
# and all three still received a biological provisional call. Scoring these
# explicitly means a cluster that is really a QC artifact shows up as one.
QC_PANELS = {
    'QC_mitochondrial': None,     # filled at runtime from KEG06_ prefix
    'QC_ribosomal':     None,     # filled at runtime from RPS/RPL prefix
    'QC_immediate_early': ['FOS', 'FOSB', 'JUN', 'JUNB', 'EGR1', 'ZFP36',
                           'BTG2', 'DUSP1', 'NR4A1'],
}
# Lymphocyte panels used ONLY to detect contamination, never to label myeloid.
CONTAM_PANELS = {
    'Contam_T':  ['CD3D', 'CD3E', 'CD3G', 'CD247', 'SKAP1', 'CAMK4', 'TRAC'],
    'Contam_B':  ['CD79A', 'CD79B', 'MS4A1', 'BANK1', 'EBF1', 'FCRL1'],
    'Contam_NK': ['NKG7', 'GNLY', 'KLRD1', 'PRF1', 'GZMH', 'CST7'],
}

# ---- rhesus symbol aliases (v2) ----
# Every entry below is a substitution of Zoey's string, applied ONLY when
# USE_RHESUS_ALIASES is True, and always reported in the resolution table so
# the change is visible. Evidence for each is in the module docstring.
USE_RHESUS_ALIASES = True
RHESUS_ALIASES = {
    'HLA-DRA':  ['MAMU-DRA', 'MAMU-DRB1', 'MAMU-DRB5'],
    'HLA-DPA1': ['MAMU-DPA', 'MAMU-DPB'],
    'HLA-DPB1': ['MAMU-DPB', 'MAMU-DPA'],
    'HLA-DQA1': ['MAMU-DQA1', 'MAMU-DQB1'],
    'FCGR3A':   ['FCGR3'],
}

# ---- provisional call gates (v2) ----
# v1 called four clusters CD206+ Macrophage off an all-negative row. A call now
# needs BOTH a positive z floor and separation from the runner-up.
MIN_Z_FOR_CALL   = 0.75
MIN_MARGIN_FOR_CALL = 0.50

# ---- lineage purity gate (v2) ----
# Applied BEFORE subclustering. A cell with lymphocyte transcripts above these
# counts is not used to define myeloid structure. Excluded cells are FLAGGED and
# reported, never deleted from the source object.
# The gate operates at CLUSTER level by default. Two per-cell designs were tried
# and both failed in testing; CELL 3b documents exactly how and why.
PURITY_GATE_CLUSTER = True    # judge whole subclusters (the design that works)
PURITY_GATE_PERCELL = False   # per-cell gating; off, see CELL 3b CRITERION NOTE
PERCELL_MARGIN      = 0.10    # only used if PURITY_GATE_PERCELL is True
CLUSTER_CONTAM_MARGIN = 0.0   # secondary: median contam must exceed median myeloid by this
CLUSTER_CONTAM_NMAD   = 5.0   # primary: median + N*MAD outlier across clusters
RECLUSTER_AFTER_GATE  = True  # recompute HVG/PCA/neighbors/leiden on the clean set
MYELOID_ANCHORS  = ['LYZ', 'CST3', 'TYROBP', 'FCER1G', 'CD68', 'AIF1', 'PSAP']

# ---- TWO-LEVEL TIER 2 (Zoey 2026-09-25) ----
# Her instruction: these markers exist only to SUBDIVIDE the Tier-1 myeloid
# population, they must never be assigned to a non-myeloid cluster, a call has
# to clear a confidence threshold, and where classifications overlap we should
# combine them into groups the data can actually separate.
#
# The merge structure below is not a guess, it is read off the correlation
# matrix from the real v1 run (11,794 cells, 16 subclusters). Panels that
# correlated above ~0.8 across subclusters are grouped, because no clustering
# can separate signatures that are measuring largely the same genes:
#
#   Classical Mono  ~ Inflammatory Mac   r=0.95   share CTSD,S100A8,S100A9
#   Classical Mono  ~ Monocyte-like DC   r=0.84
#   Mono-like DC    ~ Inflammatory Mac   r=0.83
#   Intermediate    ~ CD163+ Mac         r=0.95   share CTSD,LILRB1,MS4A7
#   Non-classical   ~ CD163+ Mac         r=0.82   share ALL 3 resolved genes
#   Intermediate    ~ Non-classical      r=0.79
#   Tissue Mac      ~ CD206+ Mac         r=0.99   share C1QA,C1QB,C1QC,LGMN,MSR1
#
# The call is made at the GROUP level, where the evidence supports it. The fine
# panels are still scored and reported underneath as sub-evidence, so nothing
# Zoey asked for is lost and a group can be split later if a discriminating
# marker arrives.
TIER2_GROUPS = {
    'Classical / Inflammatory Monocyte': ['Classical Monocytes',
                                          'Inflammatory Macrophage',
                                          'Monocyte-like DC'],
    'Non-classical / CD163+ Mono-Mac':   ['Intermediate Monocytes',
                                          'Non-classical Monocytes',
                                          'CD163+ Macrophage'],
    'Tissue Macrophage (C1Q+/CD206+)':   ['Tissue Macrophage',
                                          'CD206+ Macrophage'],
    'cDC2':                              ['cDC2'],
    'Antigen Presenting Myeloid':        ['Antigen Presenting Myeloid'],
    'IFN-Stimulated Myeloid':            ['IFN-Stimulated Myeloid'],
    'pDC':                               ['pDC'],
    'Migratory/Mature DC':               ['Migratory/Mature DC'],
    'cDC1':                              ['cDC1'],
}
CALL_AT_GROUP_LEVEL = True   # False reverts to calling on the 14 fine panels

# ---- STATES, not lineages ----
# Two of Zoey's eleven describe a CELL STATE rather than a lineage, and treating
# them as competing labels creates ties that no threshold can fix:
#
#   Antigen Presenting Myeloid is MHC-II-high. But cDC2s ARE professional
#   antigen-presenting cells, as are activated monocytes and monocyte-derived
#   DC. In testing, the cluster that should have been cDC2 came back NO CALL
#   with cDC2 at z=2.41 and Antigen Presenting Myeloid at z=2.51, margin 0.10.
#   Both calls are correct, which is exactly the problem: one is the lineage and
#   the other is what that lineage is doing.
#
#   IFN-Stimulated Myeloid is the same shape. Interferon stimulation is a state
#   any myeloid lineage can enter, and in the v1 real run it is what caused the
#   pDC cluster to be mislabelled, since pDCs constitutively carry IRF7/GZMB.
#
# This project already has the right pattern for this. The T cell annotation
# carries proliferation and CD4 activation as PER-CELL SCORES rather than
# subtype names, for precisely this reason (see shiv_tier2_annotation_write.py).
# So these two are scored, reported per subcluster and per cell, and excluded
# from the lineage call. A cluster can then be "cDC2, MHC-II-high" or
# "Classical Monocyte, IFN-stimulated", which is more informative than either
# label alone and cannot tie.
STATE_PANELS = ['Antigen Presenting Myeloid', 'IFN-Stimulated Myeloid']
TREAT_STATES_SEPARATELY = True

MIN_GENES_TO_SCORE = 2       # sc.tl.score_genes needs a set; 1 gene = raw expr
N_TOP_MARKERS = 15
ENRICH_CLIP   = 2.0
SINGLE_ANIMAL_FLAG = 0.60    # a cluster this dominated by one animal is flagged

# ---- figure conventions ----
FIG_DPI = 300
HEX = {
    'mono':   '#2C6E9B', 'mac':    '#C4453C', 'dc':     '#8FBF4D',
    'ifn':    '#E8A33D', 'apc':    '#7B5EA7', 'other':  '#8A8A8A',
    'lo':     '#3B6EA5', 'mid':    '#F2F2F2', 'hi':     '#B5322C',
}

plt.rcParams.update({
    'font.size': 28, 'axes.titlesize': 34, 'axes.labelsize': 30,
    'xtick.labelsize': 24, 'ytick.labelsize': 24, 'legend.fontsize': 18,
    'figure.dpi': 100, 'pdf.fonttype': 42, 'ps.fonttype': 42,
})
os.makedirs(OUT_DIR, exist_ok=True)
sc.settings.figdir = OUT_DIR
pd.set_option('display.width', 220)
pd.set_option('display.max_columns', 60)
pd.set_option('display.max_rows', 400)


def hr(t):
    print('\n' + '=' * 90)
    print(t)
    print('=' * 90)


def save_csv(df, name, index=True):
    p = os.path.join(OUT_DIR, name)
    df.to_csv(p, index=index)
    print(f'    saved: {name}  ({len(df)} rows)')


def save_fig(fig, stem):
    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(OUT_DIR, f'{stem}.{ext}'), dpi=FIG_DPI,
                    bbox_inches='tight')
    print(f'    saved: {stem}.pdf / .png')


# %% ============================================================
# CELL 2 — LOAD + GENE SYMBOL RESOLVER (report, do not substitute)
# ============================================================
hr('CELL 2 — load object + resolve panel gene symbols')

ADATA_IN = None
for c in ADATA_CANDIDATES:
    if os.path.exists(c):
        ADATA_IN = c
        break
if ADATA_IN is None:
    raise FileNotFoundError(
        'No input object found. Looked for:\n  ' + '\n  '.join(ADATA_CANDIDATES))
print(f'  input: {os.path.basename(ADATA_IN)}')

adata = sc.read_h5ad(ADATA_IN)
print(f'  loaded: {adata.n_obs:,} cells x {adata.n_vars:,} genes')

var_set = set(adata.var_names)
var_upper = {v.upper(): v for v in adata.var_names}


def suggest_alternatives(g, limit=6):
    """Plausible var_names for a symbol that did not resolve.

    Deliberately suggestion-only. Rhesus renames MHC to MAMU-, and paralog
    suffixes drift (FCGR3A vs FCGR3), so we surface candidates for Zoey to
    confirm rather than silently picking one.
    """
    g_up = g.upper()
    cands = []

    # MHC: HLA-DRA -> MAMU-DRA and friends
    if g_up.startswith('HLA-'):
        stem = g_up[4:]
        cands += [v for v in adata.var_names if v.upper() == f'MAMU-{stem}']
        cands += [v for v in adata.var_names
                  if v.upper().startswith('MAMU-') and stem[:2] in v.upper()]

    # paralog suffix drift: FCGR3A -> FCGR3 / FCGR3B ; APO1 -> APOC1 / APOE
    base = g_up.rstrip('ABCDEFG0123456789')
    if len(base) >= 3:
        cands += [v for v in adata.var_names if v.upper().startswith(base)]

    # family members: IFITM8 -> IFITM1/2/3
    if any(ch.isdigit() for ch in g_up):
        fam = g_up.rstrip('0123456789')
        if len(fam) >= 4:
            cands += [v for v in adata.var_names if v.upper().startswith(fam)]

    # v2: search every .var column for the symbol. Mmul10 files a number of
    # orthologs under LOC identifiers with the readable symbol only in a
    # secondary column, which prefix matching alone can never find. This is the
    # path for FCN1 / SPP1 / IFITM3, which all missed in v1.
    for col in adata.var.columns:
        try:
            ser = adata.var[col].astype(str)
        except Exception:
            continue
        hit = adata.var_names[ser.str.upper() == g_up]
        if len(hit):
            cands += [f'{h}[{col}]' for h in hit[:2]]

    seen, out = set(), []
    for c in cands:
        if c not in seen and c.upper() != g_up:
            seen.add(c)
            out.append(c)
    return out[:limit]


def apply_aliases(panel, genes):
    """Substitute rhesus symbols for the ones that do not exist in Mmul10.

    Returns (genes_after, notes). Every substitution is recorded so the
    resolution CSV shows exactly what was swapped and can be reverted by
    setting USE_RHESUS_ALIASES = False.
    """
    if not USE_RHESUS_ALIASES:
        return list(genes), []
    out, notes = [], []
    for g in genes:
        if g in var_set or g not in RHESUS_ALIASES:
            out.append(g)
            continue
        repl = [a for a in RHESUS_ALIASES[g] if a in var_set]
        if repl:
            out.extend(repl)
            notes.append(f'{g}->{"+".join(repl)}')
        else:
            out.append(g)
    return out, notes


resolve_rows = []
PANELS_RESOLVED = {}
empty_panels = []

ALL_PANELS = dict(MYELOID_PANELS)
if INCLUDE_EXTRA_PANELS:
    ALL_PANELS.update(EXTRA_PANELS)

for panel, genes in ALL_PANELS.items():
    genes, alias_notes = apply_aliases(panel, genes)
    if alias_notes:
        print(f'  {panel:28s} aliases applied: {", ".join(alias_notes)}')
        for nte in alias_notes:
            resolve_rows.append({'panel': panel, 'requested': nte.split('->')[0],
                                 'status': 'ALIASED', 'suggestions': nte})
    uniq = list(dict.fromkeys(genes))          # de-dup, keep order (CTSD twice)
    dups = len(genes) - len(uniq)
    found, missing = [], []
    for g in uniq:
        if g in var_set:
            found.append(g)
        elif g.upper() in var_upper:
            found.append(var_upper[g.upper()])   # case-only difference, safe
        else:
            missing.append(g)
            resolve_rows.append({
                'panel': panel, 'requested': g, 'status': 'MISSING',
                'suggestions': ';'.join(suggest_alternatives(g)),
            })
    for g in found:
        resolve_rows.append({'panel': panel, 'requested': g,
                             'status': 'found', 'suggestions': ''})
    PANELS_RESOLVED[panel] = found
    flag = ''
    if not found:
        flag = '   *** EMPTY — CANNOT BE SCORED ***'
        empty_panels.append(panel)
    elif len(found) == 1:
        flag = '   (single gene: scored as raw expression, not a signature)'
    print(f'  {panel:28s} {len(found)}/{len(uniq)} resolved'
          f'{f" [{dups} dup dropped]" if dups else ""}{flag}')
    if missing:
        print(f'      missing: {", ".join(missing)}')

# QC panels: mito and ribosomal are prefix-derived from the actual var_names,
# because rhesus mito genes are KEG06_ locus tags rather than MT-.
QC_PANELS['QC_mitochondrial'] = [v for v in adata.var_names
                                 if v.startswith('KEG06_') or v.startswith('MT-')]
QC_PANELS['QC_ribosomal'] = [v for v in adata.var_names
                             if v.startswith(('RPS', 'RPL'))]
print(f'\n  QC panels: mito {len(QC_PANELS["QC_mitochondrial"])} genes, '
      f'ribosomal {len(QC_PANELS["QC_ribosomal"])} genes')
QC_RESOLVED = {k: [g for g in (v or []) if g in var_set] for k, v in QC_PANELS.items()}
CONTAM_RESOLVED = {k: [g for g in v if g in var_set] for k, v in CONTAM_PANELS.items()}
for k, v in CONTAM_RESOLVED.items():
    print(f'  {k}: {len(v)}/{len(CONTAM_PANELS[k])} resolved -> {v}')

res_df = pd.DataFrame(resolve_rows)
save_csv(res_df, 'myeloid_gene_resolution.csv', index=False)

if empty_panels:
    print('\n' + '!' * 90)
    print('  EMPTY SIGNATURES — these subsets CANNOT be annotated from this run:')
    for p in empty_panels:
        print(f'    - {p}   (requested: {", ".join(MYELOID_PANELS[p])})')
    print('  Not scored, not labelled, and excluded from the heatmap. This is a')
    print('  missing-input problem, NOT a biological absence. Confirm the symbols')
    print('  with Zoey, edit MYELOID_PANELS, rerun. Everything else proceeds.')
    print('!' * 90)

miss_df = res_df[res_df['status'] == 'MISSING']
if len(miss_df):
    print('\n  unresolved symbols with candidates (for the Zoey email):')
    print(miss_df[['panel', 'requested', 'suggestions']].to_string(index=False))


# %% ============================================================
# CELL 3 — SUBSET THE MYELOID COMPARTMENT
# ============================================================
hr('CELL 3 — subset myeloid compartment')

if CELLTYPE not in adata.obs.columns:
    raise KeyError(f'{CELLTYPE} not in obs. Available: {list(adata.obs.columns)[:40]}')

print('  Tier-1 composition of the full object:')
print(adata.obs[CELLTYPE].value_counts().to_string())

present_lineages = [l for l in MYELOID_TIER1 if l in set(adata.obs[CELLTYPE].astype(str))]
absent = [l for l in MYELOID_TIER1 if l not in present_lineages]
if absent:
    print(f'\n  WARN: Tier-1 label(s) not present, skipping: {absent}')
if not present_lineages:
    raise ValueError('None of the myeloid Tier-1 labels are present. Check CELLTYPE.')

mask = adata.obs[CELLTYPE].astype(str).isin(present_lineages).values
if EXCLUDE_ANIMALS and ANIMAL_COL in adata.obs.columns:
    mask &= ~adata.obs[ANIMAL_COL].astype(str).isin(EXCLUDE_ANIMALS).values

sub = adata[mask].copy()
print(f'\n  myeloid subset: {sub.n_obs:,} cells from {present_lineages}')
print(sub.obs[CELLTYPE].value_counts().to_string())

if PROFILE_UNASSIGNED and 'Unassigned' in set(adata.obs[CELLTYPE].astype(str)):
    una = adata[adata.obs[CELLTYPE].astype(str) == 'Unassigned']
    probe = [g for g in ['LYZ', 'CD68', 'S100A8', 'CST3', 'TYROBP', 'FCER1G']
             if g in var_set]
    if probe and una.n_obs:
        X = una[:, probe].X
        X = X.toarray() if sparse.issparse(X) else np.asarray(X)
        frac = float((X > 0).any(axis=1).mean())
        print(f'\n  Unassigned pool: {una.n_obs:,} cells, '
              f'{100 * frac:.1f}% express >=1 of {probe}')
        print('  (reported only; nothing is reassigned here)')


# %% ============================================================
# CELL 3b — LINEAGE PURITY GATE  (v2, the main fix)
# ============================================================
hr('CELL 3b — lineage purity gate')

# v1 let 2,140 cells of frank T, B and NK lineage define myeloid subclusters
# (clusters 6, 7, 9, 12). Tier-1 winner-take-all scoring is what lets them in:
# a B cell with some ambient LYZ can outscore its own lineage. The fix is the
# same one already established for the T lineage, a hard marker gate, applied
# BEFORE the structure is computed so contaminants cannot define a cluster.
#
# Nothing is deleted. Excluded cells are flagged and reported, and the counts
# are written out so the size of the leak is on record.

# CRITERION NOTE - two failed designs before this one, both caught in testing.
#
# ATTEMPT 1, per-cell absolute counts: exclude a cell with >=2 nonzero
# lymphocyte markers. Against a synthetic object of known composition this
# removed 91.7% of cells including every genuine myeloid block, because with ~19
# marker genes and ordinary ambient background, two being nonzero happens by
# chance in nearly every cell. Absolute nonzero counts are not evidence.
#
# ATTEMPT 2, per-cell relative score: exclude when a lymphocyte signature scores
# above the cell's own myeloid anchor score. Better, but still removed 2,089
# cells where only 1,250 were contaminants, AND it has a fatal biological flaw:
# pDCs and migratory DCs express little LYZ / CD68 / AIF1, so a myeloid-anchor
# requirement excludes precisely the two populations v1 got wrong and that we
# added panels for. A gate that removes the cells you are trying to rescue is
# worse than no gate.
#
# WHAT WE ACTUALLY OBSERVED decides the design. In v1 the contamination was not
# scattered cells, it was whole subclusters: 6 (T/NK), 7 (B), 9 (T), 12 (NK/CTL),
# each a coherent block of hundreds of cells. So the gate belongs at CLUSTER
# level, where the evidence is a median over hundreds of cells rather than one
# noisy profile. CELL 3b only SCORES; CELL 4b clusters, judges whole clusters,
# and optionally reclusters on what survives. Per-cell gating stays available
# but off, because it is the design that failed twice.

def _score(ad, genes, name):
    genes = [g for g in genes if g in ad.var_names]
    if len(genes) < 2:
        ad.obs[name] = 0.0
        return np.zeros(ad.n_obs)
    sc.tl.score_genes(ad, gene_list=genes, score_name=name, use_raw=False)
    return ad.obs[name].values


# score on lognorm values, never on the counts layer
_t = sub.copy()
if COUNTS_LAYER in _t.layers:
    _t.X = _t.layers[COUNTS_LAYER].copy()
    sc.pp.normalize_total(_t, target_sum=1e4)
    sc.pp.log1p(_t)

contam_scores = {}
for _k, _g in CONTAM_RESOLVED.items():
    contam_scores[_k] = _score(_t, _g, f'score_{_k}')
    sub.obs[f'score_{_k}'] = contam_scores[_k]
sub.obs['score_myeloid'] = _score(_t, MYELOID_ANCHORS, 'score_myeloid')
del _t

_cm = np.vstack(list(contam_scores.values()))
sub.obs['contam_best'] = _cm.max(axis=0)
sub.obs['contam_which'] = np.array(list(contam_scores))[_cm.argmax(axis=0)]

print('  per-cell scores computed (no exclusion here):')
print(sub.obs[['score_myeloid', 'contam_best'] +
              [f'score_{k}' for k in CONTAM_RESOLVED]].describe().round(3).to_string())

if PURITY_GATE_PERCELL:
    keep = ~((sub.obs['contam_best'] > sub.obs['score_myeloid'] + PERCELL_MARGIN) &
             (sub.obs['contam_best'] > 0)).values
    print(f'\n  PURITY_GATE_PERCELL is ON: dropping {int((~keep).sum()):,} cells')
    print('  (off by default; see CRITERION NOTE for why)')
    sub_all = sub.copy()
    sub = sub[keep].copy()
else:
    sub_all = sub.copy()
    print('\n  per-cell gating OFF (default). Cluster-level gate runs in CELL 4b.')


# %% ============================================================
# CELL 4 — RECLUSTER
# ============================================================
hr('CELL 4 — recluster myeloid subset')

if COUNTS_LAYER in sub.layers:
    sub.X = sub.layers[COUNTS_LAYER].copy()
    sc.pp.normalize_total(sub, target_sum=1e4)
    sc.pp.log1p(sub)
    print(f'  renormalized from layer "{COUNTS_LAYER}"')
else:
    print(f'  WARN: layer "{COUNTS_LAYER}" absent; using .X as-is (assumed lognorm)')

sub.raw = sub
# Keep a pristine FULL-GENE copy. CELL 4b may recluster on the clean set, and
# rebuilding from an object that has already been cut down to HVGs would make
# .raw HVG-only, so every marker outside the HVG set silently vanishes from the
# signature scoring. Caught in testing, where CTSD, LILRB1, C1QA and the MAMU-
# genes were all dropped that way.
sub_fullgenes = sub.copy()
sc.pp.highly_variable_genes(sub, n_top_genes=N_HVG, batch_key=BATCH_KEY
                            if BATCH_KEY in sub.obs.columns else None)
sub = sub[:, sub.var.highly_variable].copy()
sc.pp.scale(sub, max_value=10)
sc.tl.pca(sub, n_comps=min(N_PCS, sub.n_vars - 1), svd_solver='arpack',
          random_state=RANDOM_STATE)

if DO_PREBBKNN:
    sc.pp.neighbors(sub, n_neighbors=15, n_pcs=N_NEIGH_PCS, random_state=RANDOM_STATE)
    sc.tl.umap(sub, random_state=RANDOM_STATE)
    sub.obsm['X_umap_prebbknn'] = sub.obsm['X_umap'].copy()
    print('  pre-correction UMAP stored as X_umap_prebbknn')

if BATCH_KEY in sub.obs.columns and sub.obs[BATCH_KEY].nunique() > 1:
    try:
        import bbknn
        bbknn.bbknn(sub, batch_key=BATCH_KEY, n_pcs=N_NEIGH_PCS)
        print(f'  BBKNN on "{BATCH_KEY}" ({sub.obs[BATCH_KEY].nunique()} batches)')
    except Exception as e:
        print(f'  WARN: BBKNN failed ({e}); falling back to standard neighbors')
        sc.pp.neighbors(sub, n_neighbors=15, n_pcs=N_NEIGH_PCS,
                        random_state=RANDOM_STATE)
else:
    sc.pp.neighbors(sub, n_neighbors=15, n_pcs=N_NEIGH_PCS, random_state=RANDOM_STATE)

sc.tl.umap(sub, random_state=RANDOM_STATE)
sc.tl.leiden(sub, resolution=SUB_RES, key_added='sub_leiden',
             random_state=RANDOM_STATE)
print(f'\n  sub_leiden clusters: {sub.obs["sub_leiden"].nunique()}')
print(sub.obs['sub_leiden'].value_counts().sort_index().to_string())


# %% ============================================================
# CELL 4b — CLUSTER-LEVEL PURITY GATE (v2, the design that works)
# ============================================================
hr('CELL 4b — cluster-level purity gate')

# A whole subcluster whose median lymphocyte score beats its median myeloid
# score is a lymphocyte cluster, however it got into the Tier-1 myeloid call.
# Judging at cluster level uses hundreds of cells per decision instead of one
# noisy profile, and it does not punish pDCs or migratory DCs for expressing
# little LYZ, because those form their own clusters with LOW T/B/NK scores.

def cluster_purity_table(ad):
    cols = ['score_myeloid', 'contam_best'] + [f'score_{k}' for k in CONTAM_RESOLVED]
    cols = [c for c in cols if c in ad.obs.columns]
    t = ad.obs.groupby('sub_leiden')[cols].median()
    t.insert(0, 'n_cells', ad.obs['sub_leiden'].value_counts().sort_index())
    lin_cols = [f'score_{k}' for k in CONTAM_RESOLVED if f'score_{k}' in t.columns]
    t['contam_lineage'] = t[lin_cols].idxmax(axis=1).str.replace('score_Contam_', '')

    # PRIMARY TEST: is this cluster an OUTLIER for a lymphocyte lineage relative
    # to the other clusters? A plain "contam > myeloid" comparison fires on any
    # cluster with a low myeloid score even when no lymphocyte program is
    # present at all. Caught in testing, where a mitochondria-high artifact
    # cluster (all lineage scores within 0.03 of zero) was called a T cell
    # cluster purely because its myeloid score was negative. A QC artifact is a
    # QC artifact, not contamination, and the QC columns in CELL 6 are where it
    # belongs.
    #
    # Median + N*MAD across clusters is the same robust-outlier convention used
    # for per-sample QC elsewhere in this project, and it needs no absolute
    # threshold tuned to one dataset.
    outlier = pd.Series(False, index=t.index)
    detail = pd.Series('', index=t.index)
    for c in lin_cols:
        v = t[c].astype(float)
        med = v.median()
        mad = float(np.median(np.abs(v - med)))
        if mad <= 0:
            mad = float(v.std(ddof=0)) or 1e-9
        cut = med + CLUSTER_CONTAM_NMAD * mad
        hit = (v > cut) & (v > 0)
        detail[hit] = detail[hit] + c.replace('score_Contam_', '') + ';'
        outlier |= hit
    t['contam_outlier_lineages'] = detail.str.rstrip(';')

    # SECONDARY CONFIRMATION: the lymphocyte program must also beat the cluster's
    # own myeloid score. Both must hold.
    t['is_contaminant'] = outlier & (
        t['contam_best'] > t['score_myeloid'] + CLUSTER_CONTAM_MARGIN)
    return t


pur = cluster_purity_table(sub)
print(pur.round(3).to_string())
save_csv(pur.round(4), 'myeloid_cluster_purity.csv')

bad_clusters = pur.index[pur['is_contaminant']].tolist()
n_bad = int(pur.loc[pur['is_contaminant'], 'n_cells'].sum())
if bad_clusters:
    print(f'\n  *** {len(bad_clusters)} subcluster(s) are lymphocyte, not myeloid: '
          f'{bad_clusters}  ({n_bad:,} cells, '
          f'{100 * n_bad / sub.n_obs:.1f}% of the compartment)')
    print(pur.loc[bad_clusters, ['n_cells', 'contam_lineage', 'score_myeloid',
                                 'contam_best']].round(3).to_string())
    print('\n  This is Tier-1 leakage. Worth fixing at Tier 1 as well, since the')
    print('  same cells are miscounted in every downstream composition number.')
else:
    print('\n  No subcluster is contaminant-dominated. Tier-1 myeloid call is clean.')

sub.obs['cluster_is_contaminant'] = sub.obs['sub_leiden'].map(
    pur['is_contaminant']).astype(bool)

if PURITY_GATE_CLUSTER and bad_clusters:
    sub_all = sub.copy()
    sub = sub[~sub.obs['cluster_is_contaminant'].values].copy()
    print(f'\n  retained {sub.n_obs:,} of {sub_all.n_obs:,} cells')
    if RECLUSTER_AFTER_GATE:
        print('  reclustering on the clean set (contaminants distorted the HVGs)')
        # rebuild from the FULL-GENE copy, not the HVG-subset object
        keep_idx = sub.obs_names
        sub = sub_fullgenes[keep_idx].copy()
        if COUNTS_LAYER in sub.layers:
            sub.X = sub.layers[COUNTS_LAYER].copy()
            sc.pp.normalize_total(sub, target_sum=1e4)
            sc.pp.log1p(sub)
        sub.raw = sub
        sub_fullgenes = sub.copy()
        sc.pp.highly_variable_genes(sub, n_top_genes=N_HVG,
                                    batch_key=BATCH_KEY if BATCH_KEY in sub.obs.columns else None)
        sub = sub[:, sub.var.highly_variable].copy()
        sc.pp.scale(sub, max_value=10)
        sc.tl.pca(sub, n_comps=min(N_PCS, sub.n_vars - 1), svd_solver='arpack',
                  random_state=RANDOM_STATE)
        if BATCH_KEY in sub.obs.columns and sub.obs[BATCH_KEY].nunique() > 1:
            try:
                import bbknn
                bbknn.bbknn(sub, batch_key=BATCH_KEY, n_pcs=N_NEIGH_PCS)
            except Exception as e:
                print(f'  WARN: BBKNN failed ({e}); standard neighbors')
                sc.pp.neighbors(sub, n_neighbors=15, n_pcs=N_NEIGH_PCS,
                                random_state=RANDOM_STATE)
        else:
            sc.pp.neighbors(sub, n_neighbors=15, n_pcs=N_NEIGH_PCS,
                            random_state=RANDOM_STATE)
        sc.tl.umap(sub, random_state=RANDOM_STATE)
        sc.tl.leiden(sub, resolution=SUB_RES, key_added='sub_leiden',
                     random_state=RANDOM_STATE)
        print(f'  post-gate clusters: {sub.obs["sub_leiden"].nunique()}')
        print(sub.obs['sub_leiden'].value_counts().sort_index().to_string())
        pur2 = cluster_purity_table(sub)
        still = pur2.index[pur2['is_contaminant']].tolist()
        print(f'\n  residual contaminant clusters after recluster: '
              f'{still if still else "none"}')
        save_csv(pur2.round(4), 'myeloid_cluster_purity_postgate.csv')


# %% ============================================================
# CELL 5 — SCORE SIGNATURES
# ============================================================
hr('CELL 5 — score myeloid signatures')

qc_cols = {}
for panel, genes in {**QC_RESOLVED, **CONTAM_RESOLVED}.items():
    if len(genes) < MIN_GENES_TO_SCORE:
        print(f'  SKIP {panel}: {len(genes)} gene(s)')
        continue
    key = 'qc_' + panel.replace(' ', '_')
    sc.tl.score_genes(sub, gene_list=genes, score_name=key, use_raw=True)
    qc_cols[panel] = key
print(f'  QC / contamination signatures scored: {len(qc_cols)}')

sig_cols = {}
for panel, genes in PANELS_RESOLVED.items():
    key = 'sig_' + panel.replace(' ', '_').replace('+', 'pos').replace('-', '_')
    if len(genes) == 0:
        print(f'  SKIP {panel}: no resolved genes')
        continue
    if len(genes) >= MIN_GENES_TO_SCORE:
        sc.tl.score_genes(sub, gene_list=genes, score_name=key, use_raw=True)
    else:
        src = sub.raw[:, genes[0]].X
        src = src.toarray().flatten() if sparse.issparse(src) else np.asarray(src).flatten()
        sub.obs[key] = src
    sig_cols[panel] = key
    print(f'  {panel:28s} -> {key}  ({len(genes)} genes)')

if not sig_cols:
    raise RuntimeError('No signature could be scored. Fix gene symbols first.')


# %% ============================================================
# CELL 6 — SUBCLUSTER x SIGNATURE TABLE (the one we read together)
# ============================================================
hr('CELL 6 — subcluster x signature profile')

prof = sub.obs.groupby('sub_leiden')[list(sig_cols.values())].mean()
prof.columns = [c.replace('sig_', '') for c in prof.columns]
prof.insert(0, 'n_cells', sub.obs['sub_leiden'].value_counts().sort_index())

# z-score each signature ACROSS clusters so panels with different dynamic range
# are comparable. Raw means are saved alongside so nothing is hidden.
zs = prof.drop(columns=['n_cells'])
zs = (zs - zs.mean()) / zs.std(ddof=0).replace(0, np.nan)

print('\n  RAW signature means per subcluster:')
print(prof.round(3).to_string())
print('\n  Z-SCORED across subclusters (this is the one to read):')
print(zs.round(2).to_string())

# ---- collapse fine panels into groups BEFORE calling (Zoey 2026-09-25) ----
# A group's score is the mean of its member panels' z-scores across subclusters.
# Singleton groups are unchanged. Members that failed to score (empty panel) are
# skipped, and a group with no scoring member is dropped with a warning rather
# than silently becoming a zero column.
zs_fine = zs.copy()
if CALL_AT_GROUP_LEVEL:
    def col_of(panel):
        return panel.replace(' ', '_').replace('+', 'pos').replace('-', '_')

    gcols, gmissing, scols = {}, {}, {}
    for gname, members in TIER2_GROUPS.items():
        if TREAT_STATES_SEPARATELY and all(m in STATE_PANELS for m in members):
            have_s = [col_of(m) for m in members if col_of(m) in zs_fine.columns]
            if have_s:
                scols[gname] = zs_fine[have_s].mean(axis=1)
            continue
        have = [col_of(m) for m in members if col_of(m) in zs_fine.columns]
        lack = [m for m in members if col_of(m) not in zs_fine.columns]
        if not have:
            gmissing[gname] = members
            continue
        gcols[gname] = zs_fine[have].mean(axis=1)
        if lack:
            print(f'  {gname}: built from {len(have)} of {len(members)} members '
                  f'(missing: {lack})')
    if gmissing:
        print('  WARN: group(s) dropped, no member scored:')
        for g, m in gmissing.items():
            print(f'    {g}  <- {m}')
    if not gcols:
        raise RuntimeError('No Tier-2 group could be built. Check panel scoring.')
    zs = pd.DataFrame(gcols)
    zs_state = pd.DataFrame(scols) if scols else pd.DataFrame(index=zs.index)
    print(f'\n  calling at GROUP level: {zs.shape[1]} lineage groups from '
          f'{zs_fine.shape[1]} fine panels')
    if len(zs_state.columns):
        print(f'  held out as STATES (not lineage calls): '
              f'{list(zs_state.columns)}')
        print('\n  STATE z across subclusters (read alongside the call):')
        print(zs_state.round(2).to_string())
        save_csv(zs_state.round(4), 'myeloid_subcluster_state_z.csv')
    print('\n  GROUP-level z across subclusters:')
    print(zs.round(2).to_string())
    save_csv(zs.round(4), 'myeloid_subcluster_group_z.csv')
    print('\n  fine-panel z retained as sub-evidence in '
          'myeloid_subcluster_signature_z.csv')
else:
    print('  calling at FINE panel level (CALL_AT_GROUP_LEVEL is False)')

prov = zs.idxmax(axis=1)
margin = zs.apply(lambda r: (r.nlargest(2).iloc[0] - r.nlargest(2).iloc[1])
                  if r.notna().sum() >= 2 else np.nan, axis=1)
ztop = zs.max(axis=1)

# v2 GATE. v1 used bare argmax, which called four clusters CD206+ Macrophage off
# an all-negative row (z_top -0.28 to -0.27). A call now needs a positive z floor
# AND separation from the runner-up; failing either it is "no call", and no call
# is a perfectly good answer for a cluster these panels cannot describe.
fail_z = ztop < MIN_Z_FOR_CALL
fail_m = margin < MIN_MARGIN_FOR_CALL
reason = np.where(fail_z & fail_m, f'z<{MIN_Z_FOR_CALL} and margin<{MIN_MARGIN_FOR_CALL}',
                  np.where(fail_z, f'z_top<{MIN_Z_FOR_CALL} (no signature is elevated)',
                           np.where(fail_m, f'margin<{MIN_MARGIN_FOR_CALL} (tied with runner-up)', '')))
call = pd.DataFrame({
    'n_cells': prof['n_cells'],
    'provisional': np.where(fail_z | fail_m, 'NO CALL', prov),
    'argmax_was': prov,
    'z_top': ztop.round(2),
    'margin_to_2nd': margin.round(2),
    'runner_up': zs.apply(lambda r: r.nlargest(2).index[-1]
                          if r.notna().sum() >= 2 else '', axis=1),
    'no_call_reason': reason,
})

# state scores sit next to the lineage call, never competing with it
try:
    if CALL_AT_GROUP_LEVEL and len(zs_state.columns):
        for c in zs_state.columns:
            call['state_' + c] = zs_state[c].round(2)
except NameError:
    pass

# QC and contamination context per subcluster, so an artifact cluster is visible
# next to its call rather than discovered later in the marker table.
if qc_cols:
    qprof = sub.obs.groupby('sub_leiden')[list(qc_cols.values())].mean()
    qz = (qprof - qprof.mean()) / qprof.std(ddof=0).replace(0, np.nan)
    qz.columns = [c.replace('qc_', 'z_') for c in qz.columns]
    call = call.join(qz.round(2))

if ANIMAL_COL in sub.obs.columns:
    afrac = pd.crosstab(sub.obs['sub_leiden'], sub.obs[ANIMAL_COL].astype(str))
    afrac = afrac.div(afrac.sum(axis=1), axis=0)
    call['top_animal'] = afrac.idxmax(axis=1)
    call['top_animal_frac'] = afrac.max(axis=1).round(3)
    call['single_animal_flag'] = call['top_animal_frac'] >= SINGLE_ANIMAL_FLAG

print('\n  PROVISIONAL CALL (gated)')
print(f'  a call requires z_top >= {MIN_Z_FOR_CALL} AND margin >= {MIN_MARGIN_FOR_CALL}')
print(call.to_string())

n_nocall = int((call['provisional'] == 'NO CALL').sum())
print(f'\n  {n_nocall} of {len(call)} subclusters have NO CALL.')
if n_nocall > len(call) / 2:
    print('  More than half the compartment is undescribed by these panels.')
    print('  That is a panel-coverage problem, not a clustering problem. Read')
    print('  the marker table for the no-call clusters before relabelling.')
if 'single_animal_flag' in call.columns and call['single_animal_flag'].any():
    flagged = call[call['single_animal_flag']]
    print(f'\n  {len(flagged)} subcluster(s) dominated by a single animal '
          f'(>= {SINGLE_ANIMAL_FLAG:.0%}):')
    print(flagged[['n_cells', 'top_animal', 'top_animal_frac',
                   'provisional']].to_string())
    print('  Treat these as batch or per-sample artifacts until shown otherwise.')

save_csv(prof.round(4), 'myeloid_subcluster_signature_means.csv')
save_csv(zs.round(4), 'myeloid_subcluster_signature_z.csv')
save_csv(call, 'myeloid_subcluster_provisional_call.csv')


# %% ============================================================
# CELL 6b — PANEL COLLINEARITY (v2)
# ============================================================
hr('CELL 6b — panel collinearity: which subsets CANNOT be separated')

# In v1 this was the hidden reason the calls looked arbitrary. Tissue Macrophage
# and CD206+ Macrophage correlated at r = 0.99 because CD206+ resolved to six
# genes of which five are also in Tissue Macrophage. No amount of clustering
# separates two panels that are measuring the same genes, so we state it.
# computed on the FINE panels, since the whole point is to justify the grouping
corr = zs_fine.corr()
print('  fine-panel correlation across subclusters:')
print(corr.round(2).to_string())
save_csv(corr.round(4), 'myeloid_panel_correlation.csv')

pairs = []
names = list(PANELS_RESOLVED)
for i in range(len(names)):
    for j in range(i + 1, len(names)):
        a, b = names[i], names[j]
        if a not in corr.columns.map(lambda c: c) and a not in corr.columns:
            continue
        ga, gb = set(PANELS_RESOLVED[a]), set(PANELS_RESOLVED[b])
        if not ga or not gb:
            continue
        shared = ga & gb
        ca = a.replace(' ', '_').replace('+', 'pos').replace('-', '_')
        cb = b.replace(' ', '_').replace('+', 'pos').replace('-', '_')
        r = corr.loc[ca, cb] if (ca in corr.index and cb in corr.columns) else np.nan
        pairs.append({'panel_a': a, 'panel_b': b, 'r': round(float(r), 3)
                      if np.isfinite(r) else np.nan,
                      'shared_genes': len(shared),
                      'frac_of_smaller': round(len(shared) / min(len(ga), len(gb)), 2),
                      'genes': ','.join(sorted(shared))})
pair_df = pd.DataFrame(pairs).sort_values('r', ascending=False)
print('\n  most collinear panel pairs:')
print(pair_df.head(10).to_string(index=False))
save_csv(pair_df, 'myeloid_panel_collinearity.csv', index=False)

hard = pair_df[(pair_df['r'] > 0.9) | (pair_df['frac_of_smaller'] >= 0.8)]
if len(hard):
    print('\n  *** PANELS THAT CANNOT BE SEPARATED AS WRITTEN:')
    for _, r in hard.iterrows():
        print(f'    {r["panel_a"]} <-> {r["panel_b"]}  r={r["r"]}  '
              f'shared {r["shared_genes"]} ({r["frac_of_smaller"]:.0%} of the '
              f'smaller panel): {r["genes"]}')
    print('  A subcluster assigned to either one of these is arbitrary between')
    print('  them. This needs discriminating markers from Zoey, not more code.')


# %% ============================================================
# CELL 7 — DATA-DRIVEN MARKERS PER SUBCLUSTER
# ============================================================
hr('CELL 7 — data-driven markers per subcluster')

try:
    sc.tl.rank_genes_groups(sub, 'sub_leiden', method='wilcoxon', use_raw=True)
    mk = pd.DataFrame(sub.uns['rank_genes_groups']['names']).head(N_TOP_MARKERS)
    print(mk.to_string())
    save_csv(mk, 'myeloid_subcluster_top_markers.csv', index=False)
except Exception as e:
    print(f'  WARN: rank_genes_groups failed ({e})')


# %% ============================================================
# CELL 8 — ANIMAL / TIMEPOINT / TISSUE ENRICHMENT
# ============================================================
hr('CELL 8 — composition and enrichment per subcluster')

for col, label in [(ANIMAL_COL, 'animal'), (TIME_COL, 'timepoint'),
                   (TISSUE_COL, 'tissue'), (LIB_COL, 'library')]:
    if col not in sub.obs.columns:
        print(f'  WARN: {col} absent; skipping {label} enrichment')
        continue
    ct = pd.crosstab(sub.obs['sub_leiden'], sub.obs[col].astype(str))
    frac = ct.div(ct.sum(axis=1), axis=0)
    overall = sub.obs[col].astype(str).value_counts(normalize=True)
    enr = frac.div(overall.reindex(frac.columns).values, axis=1)
    enr = enr.clip(upper=ENRICH_CLIP)
    print(f'\n  {label} fraction per subcluster:')
    print(frac.round(3).to_string())
    print(f'\n  {label} enrichment (obs/exp, clipped at {ENRICH_CLIP}):')
    print(enr.round(2).to_string())
    save_csv(frac.round(4), f'myeloid_subcluster_{label}_fraction.csv')
    save_csv(enr.round(3), f'myeloid_subcluster_{label}_enrichment.csv')

if ANIMAL_COL in sub.obs.columns:
    el = sub.obs[sub.obs[ANIMAL_COL].astype(str) == ELITE]
    if len(el):
        print(f'\n  elite controller {ELITE}: {len(el):,} myeloid cells, '
              f'subcluster spread:')
        print((el['sub_leiden'].value_counts(normalize=True).sort_index()
               .round(3)).to_string())


# %% ============================================================
# CELL 9 — FIGURES
# ============================================================
hr('CELL 9 — figures')

try:
    from matplotlib.colors import LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list(
        'shiv_div', [HEX['lo'], HEX['mid'], HEX['hi']])
    fig, ax = plt.subplots(figsize=(max(14, 1.4 * zs.shape[1]),
                                    max(10, 0.85 * zs.shape[0])))
    vmax = float(np.nanmax(np.abs(zs.values))) if np.isfinite(zs.values).any() else 1.0
    im = ax.imshow(zs.values, cmap=cmap, vmin=-vmax, vmax=vmax, aspect='auto')
    ax.set_xticks(range(zs.shape[1]))
    ax.set_xticklabels(zs.columns, rotation=45, ha='right', fontsize=22)
    ax.set_yticks(range(zs.shape[0]))
    ax.set_yticklabels([f'{i}  (n={int(prof["n_cells"].iloc[k]):,})'
                        for k, i in enumerate(zs.index)], fontsize=22)
    ax.set_title('Myeloid subcluster x signature (z across subclusters)',
                 fontsize=34, pad=20)
    cb = fig.colorbar(im, ax=ax, shrink=0.8)
    cb.ax.tick_params(labelsize=22)
    save_fig(fig, 'myeloid_subcluster_signature_heatmap')
    plt.close(fig)
except Exception as e:
    print(f'  WARN: heatmap failed ({e})')

try:
    for color in ['sub_leiden', CELLTYPE, TISSUE_COL, TIME_COL]:
        if color not in sub.obs.columns:
            continue
        sc.pl.umap(sub, color=color, show=False,
                   save=f'_myeloid_{color}.pdf', legend_fontsize=16)
    print('  saved UMAP panels to figdir')
except Exception as e:
    print(f'  WARN: UMAP panels failed ({e})')


# %% ============================================================
# CELL 10 — WRITE SUBSET OBJECT + SUMMARY
# ============================================================
hr('CELL 10 — write subset object')

sub.write(ADATA_OUT)
print(f'  wrote {ADATA_OUT}')

print(f"""
  SUMMARY
  -------
  myeloid cells (gated) : {sub.n_obs:,}   (pre-gate {sub_all.n_obs:,})
  contaminant clusters  : {len(bad_clusters)} ({n_bad:,} cells)
  subclusters           : {sub.obs['sub_leiden'].nunique()}
  signatures scored     : {len(sig_cols)} of {len(ALL_PANELS)}
  signatures EMPTY      : {len(empty_panels)}  {empty_panels if empty_panels else ''}
  unresolved symbols    : {len(miss_df)}  (see myeloid_gene_resolution.csv)
  subclusters NO CALL   : {n_nocall} of {len(call)}
  inseparable panel prs : {len(hard)}  (see myeloid_panel_collinearity.csv)

  NEXT, in this order:
   1. myeloid_panel_collinearity.csv  - any pair listed as inseparable cannot be
      annotated apart no matter what the clustering does. That is a question for
      Zoey (discriminating markers), not a code fix.
   2. myeloid_subcluster_provisional_call.csv - read no_call_reason, the QC z
      columns and single_animal_flag ALONGSIDE the call. A cluster with high
      z_QC_mitochondrial or single_animal_flag True is an artifact even when a
      signature scores well on it.
   3. myeloid_purity_gate_percell.csv - confirm the gate removed lymphocytes and
      not genuine myeloid cells, using the Tier-1 crosstab in CELL 3b.
  Only then does a label map go into shiv_tier2_myeloid_annotation_write.py.
  NOTHING is labelled by this script.
""")
print('=' * 90)
print('MYELOID TIER-2 DIAGNOSTIC COMPLETE — no labels written.')
print('=' * 90)
