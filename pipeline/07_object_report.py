#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
07 — Object report: what is in adata_annotated.h5ad and what it can answer
===========================================================================

    in   pipeline_results/05_embed_annotate/adata_annotated.h5ad
         pipeline_results/06_viral_barcode/barcode_inventory.csv   (optional)
    out  pipeline_results/07_object_report/
         07_object_report.md     <- the report, readable on its own
         plus one CSV per table

WHAT THIS IS, AND WHAT IT IS NOT
---------------------------------
This is a census of the finished object, written to be read by a person who
did not build it. It answers three questions:

  1. What is actually in here?  Every obs column, its type, how much of it is
     populated, and which stage wrote it.
  2. How many cells, broken down every way that matters: library, animal,
     timepoint, tissue, Tier 1 lineage, Tier 2 subtype, and the crossings of
     those.
  3. Which contrasts the design can actually support, and which look available
     but are confounded.

**There is no plotting in this script and no matplotlib import.** That is
deliberate, not an omission. The figure set is a separate decision and will be
made once the plots are agreed; automating a figure nobody has chosen yet
produces work that gets thrown away. This report is the input to that
conversation.

It also writes no h5ad and modifies nothing. It opens the object backed and
reads obs, var and metadata only, so it cannot corrupt the thing it audits and
it does not need the memory to hold the matrix.

WHY THE "ANSWERABLE QUESTIONS" SECTION EXISTS
----------------------------------------------
A count table tells you what was sequenced. It does not tell you what can be
compared, and the two differ here for reasons that are easy to miss:

  - Libraries are MULTIPLEXED. Several animals share one 10x channel, so
    library is not animal and a per-library number is not a per-animal number.
  - Composition differs between arms. Library C is dominated by monocytes, so
    an infected-versus-control contrast is confounded by composition before
    any biology enters.
  - 05 writes a placeholder animal for cells that failed demultiplexing. Those
    cells are real data with no animal, and they are reported separately
    rather than folded into a total or silently dropped.

So the last section states, per contrast, whether the cells exist, and names
the confound where there is one. Nothing in it is a result; it is a map of what
a result could be built from.

Environment:  conda activate sc_pre
Submit with:  sbatch run_07_object_report.sh
Author: Jake Lehle
Date:   November 2026
"""

# %% ============================================================
# CELL 1 — CONFIGURATION (every tunable lives here)
# ============================================================
import os
import hashlib
import re
import sys
import textwrap
from datetime import datetime, timezone

import numpy as np
import pandas as pd

import warnings
warnings.filterwarnings('ignore')

WORKDIR = '/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV'
ANNOT_H5AD = os.path.join(WORKDIR, 'pipeline_results', '05_embed_annotate',
                          'adata_annotated.h5ad')
PREFILTER_H5AD = os.path.join(WORKDIR, 'pipeline_results', '04_ingest_qc',
                              'adata_prefilter.h5ad')
BARCODE_INV = os.path.join(WORKDIR, 'pipeline_results', '06_viral_barcode',
                           'barcode_inventory.csv')
OUT_DIR = os.path.join(WORKDIR, 'pipeline_results', '07_object_report')
REPORT_MD = os.path.join(OUT_DIR, '07_object_report.md')

# --- the factors this report breaks cells down by -------------------------
# Discovered against the object; a missing one is reported, never fatal.
# `subspecies` was missing from the first version and it should never have
# been: Chinese versus Indian rhesus is a major biological variable here, and
# it turns out to be perfectly confounded with infection status (see the
# separability check in CELL 12).
FACTORS = ['library', 'animal_id', 'subspecies', 'timepoint', 'tissue',
           'tier1_celltype', 'tier2_subtype']
TIER_COLS = ['tier1_celltype', 'tier2_subtype']

# --- placeholders ---------------------------------------------------------
# 05 writes 'Unknown' for a cell that failed demultiplexing. It is a real cell
# with a missing label, which is neither an extra animal nor a NaN, so it is
# counted and reported in its own right.
PLACEHOLDERS = {'unknown', 'unassigned', 'unlinked', 'none', 'nan', 'na',
                'n/a', '', 'nd', 'unresolved'}

# --- study design, for the answerable-questions section --------------------
CONTROL_ANIMALS = ['34315', '41903']     # uninfected Indian controls
MIN_CELLS_PER_GROUP = 50                 # below this a group is too thin
MIN_CELLS_PER_CELLTYPE = 200             # per lineage, for a per-lineage test
TIME_ORDER = ['Pre', '21 DPI', 'Necropsy', 'Non-Infected', 'Unknown']

# --- which stage wrote which column, by pattern ---------------------------
# First match wins, so order matters: specific before general.
COLUMN_ORIGIN = [
    ('04 ingest/QC',  r'^(library|sample|animal_id|timepoint|tissue|origin|'
                      r'species|sex|batch|hashtag|donor)'),
    ('04 ingest/QC',  r'^(n_genes|n_counts|total_counts|n_genes_by_counts|'
                      r'log1p_|pct_counts_|mt$|ribo$|qc_|outlier)'),
    ('04 doublets',   r'(doublet|scrublet)'),
    ('05 embedding',  r'^(leiden|louvain|cluster)'),
    ('05 subcluster', r'^(sub_|subcluster|subleiden)'),
    ('05 Tier 1',     r'^(tier1_|lineage)'),
    ('05 Tier 2',     r'^(tier2_|subtype)'),
    ('05 signatures', r'(^|_)sig_|^dx_|_score$'),
    ('05 annotation', r'^(margin|top_|dominant_frac|marker_support|'
                      r'n_markers|alt_marker|require_)'),
    ('05 viral',      r'^(shiv_|viral_)'),
    ('05 VDJ/TCR',    r'^(has_tcr|has_ir|clonotype|IR_|chain|receptor|'
                      r'TR[ABGD]|v_call|j_call|d_call|cdr3|is_gd|gd_)'),
    ('05 CITE-seq',   r'^(adt_|ADT_|prot_|cite_|CITE)'),
    ('03 CellRanger', r'^(cr_|cellranger)'),
    # Added after the first real run put 12 columns in 'unclassified'. Each of
    # these is a real pipeline output; leaving them unclassified made the
    # census look like it had leftovers in it when it did not.
    ('04 ingest/QC',  r'^(subspecies|condition)$'),
    ('05 Tier 1',     r'^(nk_provisional|cd3_counts|dx_cd3)'),
    ('05 Tier 2',     r'^(provisional_subtype|prolif_counts|is_proliferating)$'),
    ('05 CD4/CD8',    r'^cd48_'),
    ('05 VDJ/TCR',    r'^(is_gamma_delta|tcr_recovered)$'),
    ('04 ingest/QC',  r'^mt_outlier$'),
]

# --- printing -------------------------------------------------------------
MAX_LEVELS_PRINT = 40        # levels of a categorical to print in full
MAX_COLS_WIDE = 14           # columns before a crosstab is printed long
FLOAT_FMT = '{:,.3f}'

os.makedirs(OUT_DIR, exist_ok=True)
pd.set_option('display.width', 250)
pd.set_option('display.max_columns', 80)
pd.set_option('display.max_rows', 400)
pd.set_option('display.float_format', FLOAT_FMT.format)

REPORT = []          # accumulates the markdown report


def say(s='', report=True):
    """Print to the log and, by default, into the report."""
    print(s, flush=True)
    if report:
        REPORT.append(s)


def head(title, level=2):
    bar = '=' * 94
    print('\n' + bar, flush=True)
    print(title, flush=True)
    print(bar, flush=True)
    REPORT.append('')
    REPORT.append('#' * level + ' ' + title)
    REPORT.append('')


def kv(lines):
    """Print key-value lines and put them in the report as ONE fenced block.
    Appending them individually makes markdown run them together into a
    single paragraph, because a 2-space indent is not a code block."""
    lines = [l for l in lines if l is not None]
    for l in lines:
        print('  ' + l, flush=True)
    REPORT.append('```')
    REPORT.extend(lines)
    REPORT.append('```')
    REPORT.append('')


def sub(title):
    """A sub-section separator: a rule in the log, a real heading in the md."""
    print(f'\n  --- {title} ---', flush=True)
    REPORT.append('')
    REPORT.append('### ' + title)
    REPORT.append('')


def note(body, report=True):
    """A paragraph of explanation, wrapped for the log, kept raw for the md."""
    for line in textwrap.wrap(body, 88):
        print('  ' + line, flush=True)
    if report:
        REPORT.append(body)
        REPORT.append('')


def emit(df, csv_name=None, index=True, note_text=None):
    """Print a table, write it as CSV, and put it in the report verbatim."""
    if df is None or (hasattr(df, 'empty') and df.empty):
        say('  (no rows)')
        return
    txt = df.to_string(index=index)
    print(textwrap.indent(txt, '  '), flush=True)
    REPORT.append('```')
    REPORT.append(txt)
    REPORT.append('```')
    REPORT.append('')
    if csv_name:
        path = os.path.join(OUT_DIR, csv_name)
        df.to_csv(path, index=index, lineterminator='\n')
        print(f'  -> {csv_name}', flush=True)
    if note_text:
        note(note_text)


def is_placeholder(s):
    return s.astype(str).str.strip().str.lower().isin(PLACEHOLDERS)


def real_levels(s):
    v = s.dropna().astype(str)
    return sorted(set(v[~v.str.strip().str.lower().isin(PLACEHOLDERS)]))


def classify(col):
    for origin, pat in COLUMN_ORIGIN:
        if re.search(pat, col):
            return origin
    return 'unclassified'


def pct(n, d):
    return f'{100.0 * n / d:.1f}%' if d else 'n/a'


# %% ============================================================
# CELL 2 — STRUCTURAL INVENTORY
# ============================================================
head('1. Object structure', 2)

if not os.path.exists(ANNOT_H5AD):
    sys.exit(f'  STOP: {ANNOT_H5AD} not found. Run 05 first.')

import anndata as ad

A = ad.read_h5ad(ANNOT_H5AD, backed='r')
n_obs, n_var = A.shape
stamp = datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%MZ')
size_gb = os.path.getsize(ANNOT_H5AD) / 1e9

REPORT.insert(0, f'# Rhesus_SHIV object report\n')
REPORT.insert(1, f'`adata_annotated.h5ad` — {n_obs:,} cells x {n_var:,} '
                 f'features. Generated {stamp} by `07_object_report.py`.\n')
REPORT.insert(2, 'No figures: the plot set is a separate decision. '
                 'This report is the input to it.\n')

kv([f'file      : {ANNOT_H5AD}',
    f'size      : {size_gb:.2f} GB',
    f'modified  : {datetime.fromtimestamp(os.path.getmtime(ANNOT_H5AD)).isoformat(sep=" ", timespec="seconds")}',
    f'shape     : {n_obs:,} cells x {n_var:,} features',
    f'obs cols  : {len(A.obs.columns)}',
    f'var cols  : {len(A.var.columns)}'])

def real_keys(mapping):
    """anndata returns a None key in .layers for X itself when the object is
    opened backed. It is not a layer and reporting it as one is noise, so
    only genuine string keys get through."""
    try:
        return [k for k in mapping.keys() if isinstance(k, str)]
    except Exception:
        return []


rows = []
for slot, keys in (('obsm', real_keys(A.obsm)),
                   ('varm', real_keys(A.varm)),
                   ('layers', real_keys(A.layers)),
                   ('obsp', real_keys(A.obsp)),
                   ('uns', real_keys(A.uns))):
    for k in keys:
        shape = ''
        try:
            if slot == 'obsm':
                shape = str(A.obsm[k].shape)
            elif slot == 'varm':
                shape = str(A.varm[k].shape)
            elif slot == 'layers':
                shape = str(A.layers[k].shape)
            elif slot == 'obsp':
                shape = str(A.obsp[k].shape)
            else:
                v = A.uns[k]
                shape = (f'{type(v).__name__}'
                         + (f' len={len(v)}' if hasattr(v, '__len__')
                            and not isinstance(v, str) else ''))
        except Exception as e:
            shape = f'<unreadable: {type(e).__name__}>'
        rows.append({'slot': slot, 'key': k, 'shape_or_type': shape})
emit(pd.DataFrame(rows), 'object_slots.csv', index=False)

note('`layers` is where the usable matrices live. `counts` holds raw integer '
     'counts and `X` holds the log-normalized values, so any new analysis '
     'that needs counts must take the layer and not X. `obsm` holds the '
     'embeddings: a pre-BBKNN PCA and UMAP alongside the corrected ones, '
     'which is what makes the batch-correction check in 05 possible after '
     'the fact.')

has_raw = A.raw is not None
kv([f'.raw present: {has_raw}'
    + ('  (an earlier, wider feature set is preserved here)' if has_raw else '')])

OBS = A.obs.copy()
VAR = A.var.copy()
layer_keys = real_keys(A.layers)
obsm_keys = real_keys(A.obsm)
A.file.close()
say('\n  obs and var read into memory; the matrix was never loaded.')


# %% ============================================================
# CELL 3 — COLUMN CENSUS
# ============================================================
head('2. Column census: every obs column, what it is, how complete', 2)

note('This is the table to read before writing any new script against this '
     'object. `populated` is the fraction of cells with a real value, '
     'counting the placeholder strings 05 writes for failed assignments as '
     'NOT populated, because a column that is 100% non-null and 40% '
     '"Unknown" is not a usable column at face value.')

rows = []
for c in OBS.columns:
    s = OBS[c]
    dt = str(s.dtype)
    n_na = int(s.isna().sum())
    if dt in ('object', 'category', 'bool') or str(dt).startswith('string'):
        # NOT s.fillna('nan'): on a categorical column pandas refuses a value
        # that is not already a category. astype(str) renders NaN as the
        # string 'nan', which is in PLACEHOLDERS, so subtracting n_na leaves
        # the genuine placeholder labels 05 wrote.
        ph = max(int(is_placeholder(s).sum()) - n_na, 0)
        nun = int(s.dropna().astype(str).nunique())
        example = '; '.join(map(str, real_levels(s)[:3]))
    else:
        ph = 0
        nun = int(s.dropna().nunique())
        v = pd.to_numeric(s, errors='coerce').dropna()
        example = (f'min {v.min():,.3g} / med {v.median():,.3g} / '
                   f'max {v.max():,.3g}') if len(v) else ''
    good = n_obs - n_na - ph
    rows.append({
        'column': c, 'dtype': dt, 'written_by': classify(c),
        'n_levels_or_range': nun, 'n_missing': n_na, 'n_placeholder': ph,
        'populated': pct(good, n_obs), 'example_or_range': example[:60],
    })
CENSUS = pd.DataFrame(rows).sort_values(['written_by', 'column'])
emit(CENSUS, 'column_census.csv', index=False)

say()
emit(CENSUS.groupby('written_by').size().rename('n_columns').to_frame(),
     'columns_by_stage.csv')

unc = CENSUS[CENSUS.written_by == 'unclassified']
if len(unc):
    note(f'{len(unc)} column(s) did not match any known stage pattern. They '
         'are listed here rather than hidden, because an unrecognized column '
         'is either something the pipeline wrote and this report does not '
         'know about, or something left over from an older run.')
    emit(unc[['column', 'dtype', 'n_levels_or_range', 'populated']],
         'columns_unclassified.csv', index=False)

thin = CENSUS[CENSUS.populated.str.rstrip('%').astype(float) < 95]
if len(thin):
    note('Columns below 95% populated. These are the ones that will quietly '
         'drop cells from a groupby or a filter, so check them before using '
         'one as a key.')
    emit(thin[['column', 'written_by', 'n_missing', 'n_placeholder',
               'populated']], 'columns_incomplete.csv', index=False)


# %% ============================================================
# CELL 4 — FEATURE (var) CENSUS
# ============================================================
head('3. Features', 2)

say(f'  {n_var:,} features')
rows = []
for c in VAR.columns:
    s = VAR[c]
    if s.dtype == bool or str(s.dtype) in ('object', 'category'):
        rows.append({'var_column': c, 'dtype': str(s.dtype),
                     'n_true_or_levels': int(s.sum()) if s.dtype == bool
                     else int(s.astype(str).nunique())})
    else:
        rows.append({'var_column': c, 'dtype': str(s.dtype),
                     'n_true_or_levels': int(s.notna().sum())})
emit(pd.DataFrame(rows), 'var_census.csv', index=False)

# the feature classes that matter for any downstream analysis
classes = {
    'viral (SHIV)':        r'^SHIVAD8EO',
    'mitochondrial (KEG06_)': r'^KEG06_',
    'ribosomal':           r'^(RP[SL]|MRP[SL])',
    'MHC class I (MAMU)':  r'^MAMU-(A|A3|AG|B|E|F)$',
    'MHC class II (MAMU)': r'^MAMU-(D[MOPQR])',
    'other MAMU-':         r'^MAMU-',
    'LOC (unannotated)':   r'^LOC\d+',
    'ENSMMUG':             r'^ENSMMUG',
}
rows = []
for name, pat in classes.items():
    m = VAR.index.str.match(pat)
    rows.append({'feature_class': name, 'n_features': int(np.sum(m)),
                 'pct_of_object': pct(int(np.sum(m)), n_var)})
emit(pd.DataFrame(rows), 'feature_classes.csv', index=False)

note('The LOC fraction is the practical ceiling on marker-panel work in '
     'Mmul10: a gene that exists only as a LOC identifier cannot be matched '
     'by symbol, which is why 05 carries a GENE_ALIASES map and why some '
     'human panel genes (CD27, CD24, IGHD, IL17A) could not be recovered at '
     'all. Any new panel should be resolved against this object before it is '
     'trusted, not after.')

if 'highly_variable' in VAR.columns:
    say(f'\n  highly_variable: {int(VAR.highly_variable.sum()):,} features')
for flag in ('shiv', 'mt'):
    if flag in VAR.columns:
        say(f'  var["{flag}"] flags {int(VAR[flag].astype(bool).sum()):,} '
            'features')


# %% ============================================================
# CELL 5 — ONE-WAY CELL COUNTS
# ============================================================
head('4. Cell counts, one factor at a time', 2)

present = [f for f in FACTORS if f in OBS.columns]
missing = [f for f in FACTORS if f not in OBS.columns]
if missing:
    note(f'WARNING: expected factor(s) {missing} are not in obs. Every '
         'breakdown below that would have used them is skipped.')

for f in present:
    sub(f'{f}')
    vc = OBS[f].astype(str).value_counts(dropna=False)
    d = vc.rename('cells').to_frame()
    d['pct'] = (100.0 * d.cells / n_obs).round(2)
    d['placeholder'] = is_placeholder(pd.Series(d.index)).values
    if len(d) > MAX_LEVELS_PRINT:
        note(f'{len(d)} levels; showing the top {MAX_LEVELS_PRINT}. The CSV '
             'has all of them.')
        emit(d.head(MAX_LEVELS_PRINT), f'counts_{f}.csv')
        d.to_csv(os.path.join(OUT_DIR, f'counts_{f}.csv'), lineterminator='\n')
    else:
        emit(d, f'counts_{f}.csv')
    n_ph = int(d.loc[d.placeholder, 'cells'].sum())
    if n_ph:
        note(f'{n_ph:,} cells ({pct(n_ph, n_obs)}) carry a placeholder value '
             f'for `{f}`. They are real cells with an unresolved label. They '
             'are not dropped anywhere in this report, and any analysis '
             f'keyed on `{f}` has to decide about them explicitly.')


# %% ============================================================
# CELL 6 — THE DESIGN: WHAT WAS SEQUENCED AGAINST WHAT
# ============================================================
head('5. The design', 2)

note('These crossings are the study design as the data actually realises it, '
     'which is the thing to check before proposing any comparison.')

if {'library', 'animal_id'} <= set(OBS.columns):
    sub('library x animal (the multiplexing)')
    ct = pd.crosstab(OBS.library.astype(str), OBS.animal_id.astype(str),
                     margins=True, margins_name='TOTAL')
    emit(ct, 'design_library_by_animal.csv')
    per_lib = {lib: real_levels(g.animal_id)
               for lib, g in OBS.groupby('library', observed=True)}
    multi = {k: v for k, v in per_lib.items() if len(v) > 1}
    note(f'{len(multi)} of {len(per_lib)} libraries contain more than one '
         'animal. **Library is not animal.** Two consequences that matter '
         'well beyond this table: a per-library number is not a per-animal '
         'number, and ambient RNA pools across the whole suspension before '
         'droplet formation, so in a multiplexed library the ambient '
         'transcripts in a droplet come from every animal in the pool rather '
         'than from the animal whose cell is in that droplet. That is why '
         '06 cannot attribute an ambient viral molecule to an animal by its '
         'cell barcode.')

for a, b in (('animal_id', 'timepoint'), ('animal_id', 'tissue'),
             ('timepoint', 'tissue'), ('library', 'timepoint'),
             ('library', 'tissue')):
    if {a, b} <= set(OBS.columns):
        sub(f'{a} x {b}')
        ct = pd.crosstab(OBS[a].astype(str), OBS[b].astype(str),
                         margins=True, margins_name='TOTAL')
        emit(ct, f'design_{a}_by_{b}.csv')

if {'animal_id', 'timepoint', 'tissue'} <= set(OBS.columns):
    sub('animal x timepoint x tissue (the full design cell)')
    g = (OBS.groupby([OBS.animal_id.astype(str), OBS.timepoint.astype(str),
                      OBS.tissue.astype(str)], observed=True)
         .size().rename('cells').reset_index())
    g.columns = ['animal_id', 'timepoint', 'tissue', 'cells']
    g = g.sort_values(['animal_id', 'timepoint', 'tissue'])
    emit(g, 'design_cells.csv', index=False)
    note('Every row here is a design cell with data in it. A comparison that '
         'needs a design cell absent from this table cannot be made, however '
         'reasonable it sounds.')


# %% ============================================================
# CELL 7 — COMPOSITION BY TIER 1
# ============================================================
head('6. Composition: Tier 1 lineage against every factor', 2)

if 'tier1_celltype' in OBS.columns:
    for f in [x for x in ('library', 'animal_id', 'timepoint', 'tissue')
              if x in OBS.columns]:
        sub(f'tier1_celltype x {f}, cell counts')
        ct = pd.crosstab(OBS.tier1_celltype.astype(str), OBS[f].astype(str),
                         margins=True, margins_name='TOTAL')
        emit(ct, f'composition_tier1_by_{f}.csv')
        sub(f'tier1_celltype x {f}, column percent')
        pctab = pd.crosstab(OBS.tier1_celltype.astype(str),
                            OBS[f].astype(str), normalize='columns') * 100
        emit(pctab.round(2), f'composition_tier1_by_{f}_pct.csv')

    # the composition warning, computed rather than asserted
    if 'library' in OBS.columns:
        p = pd.crosstab(OBS.tier1_celltype.astype(str),
                        OBS.library.astype(str), normalize='columns') * 100
        spread = (p.max(axis=1) - p.min(axis=1)).sort_values(ascending=False)
        sub('lineage share: spread across libraries '
            '(max - min, percentage points)')
        d = pd.DataFrame({'min_pct': p.min(axis=1).round(2),
                          'max_pct': p.max(axis=1).round(2),
                          'spread_pp': spread.round(2),
                          'max_in_library': p.idxmax(axis=1),
                          'min_in_library': p.idxmin(axis=1)})
        emit(d.sort_values('spread_pp', ascending=False),
             'composition_spread_by_library.csv')
        worst = d.spread_pp.idxmax()
        note(f'The widest spread is `{worst}`, which ranges from '
             f'{d.loc[worst, "min_pct"]:.1f}% in library '
             f'{d.loc[worst, "min_in_library"]} to '
             f'{d.loc[worst, "max_pct"]:.1f}% in library '
             f'{d.loc[worst, "max_in_library"]}. **Composition differs this '
             'much between libraries, so any contrast across libraries is '
             'confounded by composition before biology enters.** A per-cell '
             'differential test on pooled cells will read that compositional '
             'difference as a biological one. Report composition alongside '
             'every cross-library contrast, and prefer a within-lineage test '
             'or a per-animal pseudobulk over a pooled per-cell test.')
else:
    note('`tier1_celltype` is absent, so no composition tables were built.')


# %% ============================================================
# CELL 8 — TIER 2 BREAKDOWN
# ============================================================
head('7. Tier 2 subtypes', 2)

if 'tier2_subtype' in OBS.columns:
    if 'tier1_celltype' in OBS.columns:
        sub('tier2_subtype within tier1_celltype')
        d = (OBS.groupby([OBS.tier1_celltype.astype(str),
                          OBS.tier2_subtype.astype(str)], observed=True)
             .size().rename('cells').reset_index())
        d.columns = ['tier1_celltype', 'tier2_subtype', 'cells']
        tot = d.groupby('tier1_celltype').cells.transform('sum')
        d['pct_of_lineage'] = (100.0 * d.cells / tot).round(2)
        emit(d.sort_values(['tier1_celltype', 'cells'],
                           ascending=[True, False]),
             'tier2_within_tier1.csv', index=False)

    unres = OBS.tier2_subtype.astype(str)
    n_un = int(is_placeholder(unres).sum())
    note(f'{n_un:,} cells ({pct(n_un, n_obs)}) have no resolved Tier 2 '
         'subtype. Tier 2 was pinned by agreement rather than finished, so '
         'this number is a known state of the object and not a defect to '
         'chase. A Tier 2 label is a subcluster-mean assignment with a '
         'reported margin, so treat a subtype with a small margin as '
         'provisional.')

    for f in [x for x in ('timepoint', 'tissue', 'animal_id')
              if x in OBS.columns]:
        sub(f'tier2_subtype x {f}')
        ct = pd.crosstab(OBS.tier2_subtype.astype(str), OBS[f].astype(str),
                         margins=True, margins_name='TOTAL')
        emit(ct, f'tier2_by_{f}.csv')

    # margin / confidence columns, if 05 wrote them
    mcols = [c for c in OBS.columns
             if re.search(r'margin|dominant_frac|top_family_range', c)]
    if mcols:
        sub('Tier 2 confidence columns')
        emit(OBS[mcols].describe().T.round(4), 'tier2_confidence_stats.csv')
        note('`dominant_frac` is a HOMOGENEITY measure, not a confidence '
             'measure: a subcluster where 64% of cells vote the same way is '
             'homogeneous, which says nothing about whether the winning '
             'label is right. Read it with the margin, never instead of it.')
else:
    note('`tier2_subtype` is absent.')


# %% ============================================================
# CELL 9 — QC AND INTEGRITY
# ============================================================
head('8. Integrity checks: duplicated columns and compositional outliers', 2)

# ---- DUPLICATED COLUMNS ------------------------------------------------
# Two measurements of different things cannot produce identical values. When
# they do, one column was written twice under two names, and anybody who uses
# the second one gets a silently wrong answer. The first real run found
# ADT_CD163_BDAbSeq and ADT_CD169_BDAbSeq identical on every value.
sub('columns with byte-identical contents')
note('Candidates are grouped by an ORDER-PRESERVING hash and then each group '
     'is confirmed by exact comparison, so this table has no false positives. '
     'The first version of this check summed the per-element hashes, which '
     'for a boolean column collapses to a function of the True count alone: '
     'any two bool columns with the same number of True values collided even '
     'when the True cells were entirely different. That produced four '
     'spurious groups on the 2026-11-04 run, one of which wrongly implicated '
     '`shiv_confident`.')


def col_sig(s):
    """A hash of the column's values IN ORDER. Summing per-element hashes
    discards position, which is what caused the false positives."""
    try:
        h = pd.util.hash_pandas_object(
            s.astype(str) if str(s.dtype) in ('object', 'category') else s,
            index=False)
        return hashlib.sha1(h.values.tobytes()).hexdigest()
    except Exception:
        return None


def truly_equal(a, b):
    """Exact comparison, NaN-aware, so a hash collision cannot survive."""
    x, y = OBS[a], OBS[b]
    if str(x.dtype) in ('object', 'category') or str(y.dtype) in ('object',
                                                                  'category'):
        return x.astype(str).equals(y.astype(str))
    xv, yv = pd.to_numeric(x, errors='coerce'), pd.to_numeric(y, errors='coerce')
    return bool(np.array_equal(xv.values, yv.values, equal_nan=True))


sig = {}
for c in OBS.columns:
    k = col_sig(OBS[c])
    if k is not None:
        sig.setdefault(k, []).append(c)

dups = []
for group in (v for v in sig.values() if len(v) > 1):
    # verify pairwise against the first member; keep only what really matches
    confirmed = [group[0]] + [c for c in group[1:]
                              if truly_equal(group[0], c)]
    if len(confirmed) > 1:
        dups.append(confirmed)
if dups:
    # Not every duplicate is a defect. Two names for one quantity is
    # redundancy; two names for two DIFFERENT measurements is a bug. Saying
    # which is which is the difference between an actionable table and an
    # alarming one.
    BENIGN = [
        (r'^n_genes', r'^n_genes', 'scanpy writes n_genes and '
         'n_genes_by_counts as the same quantity'),
        (r'_shiv$', r'^shiv_umi', 'the same viral UMI count under two names'),
        (r'^dx_NK_cells$', r'^dx_nk_old$', '05 keeps the superseded NK panel '
         'beside the current one; identical means the panel did not change'),
        (r'^shiv_whitelist$', r'^shiv_detected$', 'detected is defined as '
         'being on the whitelist, so identity is by construction'),
        (r'^tcr_recovered$', r'^marker_support_exempt$',
         '05 assigns both from the same boolean (lines 1524 and 1538), so '
         'they are one flag under two names'),
        (r'^outlier$', r'^qc_would_have_dropped$',
         'survivorship: 04 drops (outlier | mt_outlier) & ~shiv_whitelist, '
         'so every outlier left in this object is whitelisted and both '
         'columns collapse onto the rescued cells'),
    ]

    def judge(g):
        for pa, pb, why in BENIGN:
            if any(re.search(pa, x) for x in g) and any(re.search(pb, x)
                                                        for x in g):
                return 'expected alias', why
        if all(c.startswith('ADT_') for c in g):
            return 'BUG', ('distinct antibodies cannot produce identical '
                           'values; one column was written twice')
        return 'REVIEW', 'no known reason for these to be identical'

    rows = []
    for g in dups:
        v, why = judge(g)
        rows.append({'identical_group': ' == '.join(g), 'n_columns': len(g),
                     'dtype': str(OBS[g[0]].dtype), 'verdict': v,
                     'reason': why})
    D = pd.DataFrame(rows).sort_values('verdict')
    emit(D, 'duplicate_columns.csv', index=False)
    bad = D[D.verdict != 'expected alias']
    if len(bad):
        note(f'**{len(bad)} group(s) need attention.** Two measurements of '
             'different things cannot agree on every one of '
             f'{n_obs:,} cells, so one column was written twice under two '
             'names and anyone reaching for the second gets a silently wrong '
             'answer. Nothing else in the pipeline would catch this.')
    else:
        note('Every group is an expected alias: two names for one quantity, '
             'not two measurements disagreeing. Nothing to fix.')
    for g in dups:
        if all(c.startswith('ADT_') for c in g):
            lo = OBS[g[0]]
            note(f'`{"` and `".join(g)}` are ADT columns and are also '
                 f'entirely non-detected (max {lo.max():.3f} after CLR), so '
                 'they are both duplicated and uninformative. Check the '
                 'feature reference and the ADT ingestion in 04 before '
                 'treating either as a measurement.')

    # The outlier case is benign but it leaves a trap, so say so explicitly
    # rather than letting the 'expected alias' label imply there is nothing
    # to know.
    if 'outlier' in OBS.columns and any(
            {'outlier', 'qc_would_have_dropped'} <= set(g) for g in dups):
        n_out = int(OBS['outlier'].astype(bool).sum())
        note(f'**`outlier` does not mean what its name says in this object.** '
             f'04 computes it from a 5-MAD test on three QC metrics, which '
             'flags a large fraction of cells, and then drops every outlier '
             'that is not on the viral whitelist. So the '
             f'{n_out:,} cells still carrying `outlier = True` here are not '
             '"the QC outliers": they are the whitelisted cells that QC would '
             'have removed and that were deliberately rescued. The general '
             'outlier population is gone from this object and lives only in '
             '`adata_prefilter.h5ad`. Anyone treating this column as a QC '
             f'flag on the annotated object will be reasoning about {n_out:,} '
             'rescued cells while believing they are looking at all outliers.')
else:
    say('  none. No two obs columns hold identical values.')

# ---- CONSTANT COLUMNS --------------------------------------------------
sub('columns with only one value')
const = [c for c in OBS.columns if OBS[c].astype(str).nunique() <= 1]
if const:
    emit(pd.DataFrame({'column': const,
                       'the_only_value': [str(OBS[c].iloc[0]) for c in const]}),
         'constant_columns.csv', index=False)
    note('A column with one value carries no information. Either the logic '
         'that writes it never fired, or it is a leftover. `doublet_whitelisted` '
         'being all-False means the doublet whitelist never matched anything, '
         'which is worth knowing before relying on it.')
else:
    say('  none.')

# ---- COMPOSITIONAL OUTLIERS BY ANIMAL ----------------------------------
# The composition spread in section 6 is reported per library, which hides
# whether a library is odd or one ANIMAL inside it is. That distinction
# decides whether the cause is batch or biology, so compute it.
if {'animal_id', 'tier1_celltype'} <= set(OBS.columns):
    sub('per-animal lineage share, and the outliers in it')
    p = (pd.crosstab(OBS.tier1_celltype.astype(str), OBS.animal_id.astype(str),
                     normalize='columns') * 100)
    real_cols = [c for c in p.columns if c.lower() not in PLACEHOLDERS]
    p = p[real_cols]
    med = p.median(axis=1)
    rows = []
    for lin in p.index:
        for an in p.columns:
            v, m = p.loc[lin, an], med[lin]
            if m > 0.5 and v > 3 * m and v > 10:
                rows.append({'tier1_celltype': lin, 'animal_id': an,
                             'this_animal_pct': round(v, 2),
                             'median_animal_pct': round(m, 2),
                             'fold_over_median': round(v / m, 1)})
    if rows:
        emit(pd.DataFrame(rows).sort_values('fold_over_median',
                                            ascending=False),
             'composition_animal_outliers.csv', index=False)
        note('**A single animal, not a library, drives the composition '
             'spread.** That distinction matters: a library-level difference '
             'points at batch or sample handling, whereas one animal '
             'differing this far from its peers is either real biology in '
             'that animal, a sample-handling problem specific to it, or '
             'mis-demultiplexing. It cannot be corrected by batch '
             'correction, and it means the animals involved cannot be pooled '
             'into an arm without the composition being reported alongside.')
    else:
        say('  no animal exceeds 3x the median share of any lineage.')

head('8b. QC state and data integrity', 2)

qc_rows = []
for c in ('qc_exclude', 'predicted_doublet', 'doublet', 'scrublet_call'):
    if c in OBS.columns:
        s = OBS[c]
        try:
            n = int(s.astype(bool).sum())
        except Exception:
            n = int((s.astype(str).str.lower()
                     .isin(['true', '1', 'doublet', 'yes'])).sum())
        qc_rows.append({'flag': c, 'n_cells': n, 'pct': pct(n, n_obs)})
if qc_rows:
    emit(pd.DataFrame(qc_rows), 'qc_flags.csv', index=False)
    note('These cells are FLAGGED and still present. 04 carries doublets '
         'forward rather than deleting them, after the original hardcoded '
         'removal was found to take 1,134 cells of which 1,121 came from '
         'library G alone, which is a batch artifact rather than a biological '
         'doublet rate. Any analysis must filter on `qc_exclude == False` '
         'itself; nothing upstream has done it.')

numcols = [c for c in ('n_genes_by_counts', 'n_genes', 'total_counts',
                       'pct_counts_mt', 'pct_counts_ribo')
           if c in OBS.columns]
if numcols:
    sub('QC distributions, whole object')
    emit(OBS[numcols].describe().T.round(3), 'qc_distributions.csv')
    if 'library' in OBS.columns:
        sub('median QC per library')
        emit(OBS.groupby(OBS.library.astype(str))[numcols].median().round(3),
             'qc_median_by_library.csv')
        note('Depth differs between libraries, and depth drives detection. '
             'A lineage or a transcript that looks absent from a shallow '
             'library may simply be under-sampled there, so a presence or '
             'absence claim across libraries needs the depth column beside '
             'it. This is the same confound behind the necropsy observation '
             'in 06: virus present in LN and absent from PBMC is '
             'inseparable from the depth difference between those libraries.')

# duplicate and index integrity
dup = int(OBS.index.duplicated().sum())
say(f'\n  duplicated cell barcodes: {dup}')
if dup:
    note('WARNING: duplicated obs index entries. Any merge or map keyed on '
         'the barcode will silently multiply rows. Fix before proceeding.')
say(f'  barcode format example : {OBS.index[0]}')
note("04 writes the barcode as `<bare>-1-<library>`, while the CB tag in a "
     "BAM is `<bare>-1`. Any join from alignment data to this object has to "
     "add the library suffix, which is what 06 does.")


# %% ============================================================
# CELL 10 — THE VIRAL COLUMNS
# ============================================================
head('9. Viral columns, and what they can carry', 2)

vcols = [c for c in OBS.columns if re.match(r'^(shiv_|viral_)', c)]
if vcols:
    emit(OBS[vcols].describe(include='all').T, 'viral_columns_stats.csv')
    for c in vcols:
        if OBS[c].dtype == bool or OBS[c].astype(str).nunique() <= 10:
            sub(f'{c}')
            emit(OBS[c].astype(str).value_counts().rename('cells').to_frame(),
                 f'viral_{c}_counts.csv')

    if {'shiv_call', 'tier1_celltype'} <= set(OBS.columns):
        sub('shiv_call x tier1_celltype')
        emit(pd.crosstab(OBS.tier1_celltype.astype(str),
                         OBS.shiv_call.astype(str), margins=True,
                         margins_name='TOTAL'),
             'viral_call_by_tier1.csv')

    note('**Read this before using any viral column.** 05b established that '
         'per-cell viral attribution below `shiv_confident` is not '
         'publishable: the strongest B cell carried less viral signal than '
         'the strongest erythrocyte, and erythrocytes cannot be infected. '
         'The signal at low UMI counts is ambient RNA, which deposits roughly '
         'a fixed number of molecules per droplet and therefore inflates the '
         'apparent rate in damaged or near-empty droplets. Use the confident '
         'column only, and treat a cell-level viral call as a statement '
         'about the droplet rather than about the cell.')
    note('06 sharpened this rather than softening it. The one viral barcode '
         'in that stage that could be tied to a cell at all was in an '
         'erythrocyte. So the barcode arm yields a per-animal and per-library '
         'inventory, and there is no cell-level viral signal in this dataset '
         'to build a tropism claim on.')
else:
    note('No viral columns found in obs.')

if os.path.exists(BARCODE_INV):
    sub('06 barcode inventory')
    emit(pd.read_csv(BARCODE_INV), 'barcode_inventory_copy.csv', index=False)
    note('Carried through from 06 so this report stands alone. `animal_evidence` '
         'records how the animal was established: `cell` (the barcode shares a '
         'droplet with a cell assigned to that animal), `library` (the library '
         'holds exactly one animal), or `catalog` (the core is on that '
         "animal's plasma list). For an ambient molecule in a multiplexed "
         'library the catalog is the stronger of those, not the weaker, '
         'because droplet-sharing says little about the source animal.')
else:
    note(f'06 inventory not found at {BARCODE_INV}; run 06 for the barcode '
         'section.')


# %% ============================================================
# CELL 11 — VDJ / TCR AND CITE-seq COVERAGE
# ============================================================
head('10. VDJ/TCR and CITE-seq coverage', 2)

# This used a narrower pattern than COLUMN_ORIGIN, so the first real run
# listed has_tcr and chain_pairing but silently omitted receptor_subtype,
# is_gamma_delta and tcr_recovered. One pattern, used in both places.
TCR_PAT = (r'^(has_tcr|has_ir|clonotype|clone_id|IR_|is_gd|gd_|chain|'
           r'receptor|is_gamma_delta|tcr_recovered|cdr3|[vdj]_call)')
tcols = [c for c in OBS.columns if re.match(TCR_PAT, c)]
if tcols:
    kv([f'TCR-related columns ({len(tcols)}): ' + ', '.join(tcols)])
    for c in tcols:
        if OBS[c].dtype == bool:
            n = int(OBS[c].sum())
            say(f'    {c:<28} {n:,} cells ({pct(n, n_obs)})')
    if {'has_tcr', 'tier1_celltype'} <= set(OBS.columns):
        sub('has_tcr x tier1_celltype')
        ct = pd.crosstab(OBS.tier1_celltype.astype(str),
                         OBS.has_tcr.astype(str), margins=True,
                         margins_name='TOTAL')
        emit(ct, 'tcr_by_tier1.csv')
        note('A TCR in a B cell or a monocyte is ambient or a doublet, not a '
             'receptor, which is why 05 quarantines TCR-based T-cell recovery '
             'away from those lineages. The off-lineage rate in this table '
             'is a usable estimate of the ambient contamination level in the '
             'VDJ library.')
    cl = [c for c in OBS.columns if re.search(r'clonotype|clone_id', c)]
    cl = [c for c in cl if OBS[c].astype(str).nunique() > 20]
    if not cl:
        note('**There is no clonotype ID column in this object.** '
             '`chain_pairing` and `receptor_subtype` describe each cell\'s '
             'chain configuration, which is not the same thing: neither one '
             'groups cells that share a receptor. So clonal expansion, '
             'clonal overlap between timepoints and clonotype-to-phenotype '
             'links cannot be computed from this object as it stands. 05 '
             'runs scirpy far enough to produce the chain columns but does '
             'not write a clonotype definition, so recovering this needs a '
             '`scirpy.tl.define_clonotypes` step added to 05 and a re-run. '
             'Worth deciding on before any T-cell analysis is planned, '
             'because clonal structure is usually the point of adding VDJ.')
    for c in cl:
        s = OBS[c].dropna().astype(str)
        s = s[~s.str.lower().isin(PLACEHOLDERS)]
        if len(s):
            vc = s.value_counts()
            kv([f'{c}: {s.nunique():,} distinct clonotypes over {len(s):,} '
                f'cells',
                f'  largest clonotype {vc.iloc[0]:,} cells; '
                f'{int((vc > 1).sum()):,} clonotypes expanded (>1 cell)'])
else:
    note('No VDJ/TCR columns found in obs.')

acols = [c for c in OBS.columns if re.match(r'^(adt_|ADT_|prot_|cite_)', c)]
if acols:
    say(f'\n  CITE-seq (ADT) columns: {len(acols)}')
    emit(OBS[acols].describe().T.round(3), 'adt_stats.csv')
else:
    note('No ADT columns in obs. If CITE-seq lives in a layer or in .obsm '
         'instead, the object-structure table in section 1 shows where.')


# %% ============================================================
# CELL 12 — WHAT CAN THIS OBJECT ANSWER
# ============================================================
head('11. What this object can and cannot answer', 2)

note('Derived from the design tables above, not from an assumption about the '
     'study. A contrast appears as available only if the cells exist on both '
     'sides at or above the thresholds in CELL 1 '
     f'({MIN_CELLS_PER_GROUP} cells per group, '
     f'{MIN_CELLS_PER_CELLTYPE} per lineage).')

verdicts = []


def verdict(question, available, detail, confound=''):
    verdicts.append({'question': question,
                     'status': 'available' if available else 'NOT available',
                     'detail': detail, 'confound': confound})


if {'animal_id', 'timepoint'} <= set(OBS.columns):
    d = (OBS.assign(a=OBS.animal_id.astype(str), t=OBS.timepoint.astype(str))
         .query('a.str.lower() not in @PLACEHOLDERS', engine='python'))
    tab = d.groupby(['a', 't'], observed=True).size().unstack(fill_value=0)
    tab = tab.loc[:, [c for c in TIME_ORDER if c in tab.columns]
                  + [c for c in tab.columns if c not in TIME_ORDER]]
    sub('cells per animal per timepoint (placeholders excluded)')
    emit(tab, 'answerable_animal_by_timepoint.csv')

    longit = {a: [t for t in r.index if r[t] >= MIN_CELLS_PER_GROUP]
              for a, r in tab.iterrows()}
    multi = {a: ts for a, ts in longit.items() if len(ts) > 1}
    verdict('Within-animal change over time',
            len(multi) > 0,
            f'{len(multi)} of {len(longit)} animals have >= '
            f'{MIN_CELLS_PER_GROUP} cells at more than one timepoint: '
            + '; '.join(f'{a} ({len(ts)} tps)' for a, ts in multi.items()),
            'Paired within animal, so this is the strongest design available '
            'here. Still confounded by library depth and composition if the '
            'timepoints sit in different libraries, which they do.')

    infected = [a for a in tab.index if a not in CONTROL_ANIMALS]
    ctrl = [a for a in tab.index if a in CONTROL_ANIMALS]

    # --- is the infected/control split separable from library? ---------
    # This is the question that decides whether an infected-versus-control
    # contrast is interpretable at all. If the two arms sit in disjoint sets
    # of libraries then arm and batch are the same variable and nothing can
    # separate them. If some library holds BOTH arms, that library supports a
    # within-batch comparison, which is far stronger than the pooled one.
    within = []
    if 'library' in OBS.columns and ctrl and infected:
        sub('is the infected/control split separable from library?')
        lt = pd.crosstab(
            OBS.library.astype(str),
            OBS.animal_id.astype(str).map(
                lambda a: 'control' if a in CONTROL_ANIMALS
                else ('placeholder' if a.lower() in PLACEHOLDERS
                      else 'infected')))
        for col in ('infected', 'control'):
            if col not in lt.columns:
                lt[col] = 0
        emit(lt, 'confound_arm_by_library.csv')
        within = [lib for lib in lt.index
                  if lt.loc[lib, 'infected'] >= MIN_CELLS_PER_GROUP
                  and lt.loc[lib, 'control'] >= MIN_CELLS_PER_GROUP]
        if within:
            note(f'**Librar(ies) {within} contain both arms**, so a '
                 'within-library infected-versus-control comparison is '
                 'possible there and is free of the batch confound that '
                 'defeats the pooled version. This is the cleanest contrast '
                 'in the dataset and it is worth building on before anything '
                 'pooled. Check the timepoint composition inside that '
                 'library before using it, because sharing a library does '
                 'not mean sharing a timepoint.')
        else:
            note('**No library contains both arms above threshold.** Arm and '
                 'library are therefore the same variable, and no analysis '
                 'can separate an infection effect from a batch effect. Any '
                 'infected-versus-control result would be uninterpretable.')

    # --- are the controls timepoint-matched? ---------------------------
    ctrl_tps, inf_tps = set(), set()
    for a in ctrl:
        ctrl_tps |= {t for t in tab.columns if tab.loc[a, t] >= MIN_CELLS_PER_GROUP}
    for a in infected:
        inf_tps |= {t for t in tab.columns if tab.loc[a, t] >= MIN_CELLS_PER_GROUP}
    shared_tps = sorted(ctrl_tps & inf_tps)

    conf = (f'n={len(ctrl)} controls is the binding constraint: a '
            f'{len(ctrl)}-animal control arm cannot separate '
            'animal-to-animal variation from an infection effect. ')
    if not shared_tps:
        conf += (f'The controls also share NO timepoint with the infected '
                 f'animals (controls at {sorted(ctrl_tps)}, infected at '
                 f'{sorted(inf_tps)}), so the comparison is not '
                 'timepoint-matched and an infection effect cannot be '
                 'separated from a sampling-occasion effect. ')
    else:
        conf += f'Shared timepoint(s) for a matched comparison: {shared_tps}. '
    conf += ('Composition also differs sharply between the arms (section 6 '
             'and the per-animal outlier table in section 8), so a pooled '
             'per-cell test will read composition as biology.')
    if within:
        conf += (f' Mitigation: librar(ies) {within} hold both arms and '
                 'support a within-batch comparison.')
    verdict('Infected versus uninfected control',
            len(ctrl) > 0 and len(infected) > 0,
            f'{len(infected)} infected animal(s), {len(ctrl)} control(s): '
            f'infected {infected}, control {ctrl}',
            conf)

    # --- any factor perfectly confounded with the arm? -----------------
    for f in ('subspecies', 'tissue'):
        if f not in OBS.columns:
            continue
        m = OBS.assign(
            arm=OBS.animal_id.astype(str).map(
                lambda a: 'control' if a in CONTROL_ANIMALS
                else ('placeholder' if a.lower() in PLACEHOLDERS
                      else 'infected')),
            k=OBS[f].astype(str))
        m = m[(m.arm != 'placeholder') & (~m.k.str.lower().isin(PLACEHOLDERS))]
        ct = pd.crosstab(m.k, m.arm)
        if ct.empty or ct.shape[1] < 2:
            continue
        sub(f'is {f} separable from infection status?')
        emit(ct, f'confound_{f}_by_arm.csv')
        perfect = bool((ct == 0).any(axis=1).all())
        if perfect:
            note(f'**`{f}` is perfectly confounded with infection status.** '
                 f'Every level of `{f}` appears in exactly one arm, so an '
                 f'infected-versus-control difference and a `{f}` difference '
                 'are the same comparison and no analysis of this dataset '
                 'can tell them apart. Any such result has to be reported '
                 f'with that stated, not as an infection effect.')
            verdict(f'Infection effect separable from {f}', False,
                    f'No. Every `{f}` level sits in exactly one arm.',
                    'Structural, not a power problem. It cannot be fixed by '
                    'adding cells.')
        else:
            note(f'`{f}` appears in both arms, so the two are at least '
                 'partially separable.')

if {'tissue', 'animal_id'} <= set(OBS.columns):
    d = OBS.assign(a=OBS.animal_id.astype(str), s=OBS.tissue.astype(str))
    tab = d.groupby(['a', 's'], observed=True).size().unstack(fill_value=0)
    paired = [a for a, r in tab.iterrows()
              if (r >= MIN_CELLS_PER_GROUP).sum() > 1
              and a.lower() not in PLACEHOLDERS]
    sub('cells per animal per tissue')
    emit(tab, 'answerable_animal_by_tissue.csv')
    verdict('PBMC versus lymph node within the same animal',
            len(paired) > 0,
            f'{len(paired)} animal(s) have both tissues above threshold: '
            f'{paired}',
            'Tissue and library are largely collinear here, so a tissue '
            'effect and a library effect cannot be fully separated.')

if 'tier1_celltype' in OBS.columns:
    vc = OBS.tier1_celltype.astype(str).value_counts()
    big = vc[vc >= MIN_CELLS_PER_CELLTYPE]
    small = vc[vc < MIN_CELLS_PER_CELLTYPE]
    sub('lineages large enough for a within-lineage analysis')
    emit(pd.DataFrame({'cells': big}), 'answerable_lineages.csv')
    verdict('Within-lineage differential expression',
            len(big) > 0,
            f'{len(big)} lineage(s) at or above {MIN_CELLS_PER_CELLTYPE} '
            f'cells: {list(big.index)}',
            f'{len(small)} lineage(s) are below it and should not carry a '
            f'differential claim: {list(small.index)}')
    if 'animal_id' in OBS.columns:
        pb = (OBS.assign(a=OBS.animal_id.astype(str),
                         c=OBS.tier1_celltype.astype(str))
              .groupby(['c', 'a'], observed=True).size().unstack(fill_value=0))
        ok = (pb >= MIN_CELLS_PER_GROUP).sum(axis=1)
        sub('animals per lineage with enough cells for pseudobulk')
        emit(pd.DataFrame({'animals_above_threshold': ok.sort_values(
            ascending=False)}), 'answerable_pseudobulk.csv')
        note('Per-animal pseudobulk is the analysis this design supports '
             'best: it treats the animal as the replicate, which is the unit '
             'the experiment randomised, and it is immune to the '
             'compositional confound that defeats a pooled per-cell test. '
             'A lineage needs enough animals in this column on both sides of '
             'whatever contrast is being made.')

verdict('Cell-level viral tropism', False,
        'No. 05b and 06 both land in the same place: low-UMI viral signal is '
        'ambient, and the one 06 barcode attributable to a cell was in an '
        'erythrocyte.',
        'Not fixable by analysis. It needs targeted enrichment.')
verdict('Viral barcode diversity or lineage tracing', False,
        '3 barcodes recovered across the whole experiment, 1 attributable to '
        'an animal by cell.',
        'Coverage, not method: the cassette is 34bp of a 10,296bp contig. '
        'Targeted VpxF1 amplicon enrichment is the only route.')

V = pd.DataFrame(verdicts)
# Printed as blocks, not as a table. The detail and confound fields are
# sentences, and a DataFrame with sentence-length cells is 400 characters wide
# and unreadable in a log or a markdown file.
head('Verdicts', 3)
for r in V.itertuples():
    tag = 'AVAILABLE    ' if r.status == 'available' else 'NOT AVAILABLE'
    print(f'\n  [{tag}] {r.question}', flush=True)
    REPORT.append(f'**[{r.status.upper()}] {r.question}**')
    REPORT.append('')
    for label, body in (('', r.detail), ('confound: ', r.confound)):
        if not body:
            continue
        for i, line in enumerate(textwrap.wrap(label + body, 84)):
            print('      ' + line, flush=True)
        REPORT.append(f'- {label}{body}')
    REPORT.append('')
V.to_csv(os.path.join(OUT_DIR, 'answerable_questions.csv'), index=False,
         lineterminator='\n')
print('\n  -> answerable_questions.csv', flush=True)


# %% ============================================================
# CELL 13 — THE SHORTLIST
# ============================================================
head('12. The columns that matter', 2)

note('Everything above is the full census. This is the subset to reach for '
     'first, with what each is good for and what it cannot do.')

SHORT = [
    ('library', 'the 10x channel; the batch unit and BATCH_KEY in 05',
     'not an animal and not a treatment arm; several animals share one'),
    ('animal_id', 'the experimental replicate; the right unit for pseudobulk',
     'carries a placeholder level for cells that failed demultiplexing'),
    ('timepoint', 'Pre / 21 DPI / Necropsy', 'largely collinear with library'),
    ('tissue', 'PBMC or LN', 'largely collinear with library'),
    ('tier1_celltype', 'the lineage call; stable and the basis of everything',
     'Monocytes still contain some T and NK contamination that only a Tier 1 '
     'change would remove'),
    ('tier2_subtype', 'subtype within lineage, by subcluster-mean signature',
     'pinned by agreement rather than finished; read the margin with it'),
    ('qc_exclude', 'the QC filter you must apply yourself',
     'nothing upstream applied it; doublets are flagged, not removed'),
    ('outlier', 'NOT a general QC flag on this object',
     'survivorship makes it the whitelisted cells QC would have dropped and '
     'that were rescued; the real outlier population is only in '
     'adata_prefilter.h5ad'),
    ('counts (layer)', 'raw integer counts for any new normalization or DE',
     'X is log-normalized, so do not feed X to a count-based model'),
    ('X_pca / X_umap and the pre-BBKNN pair (obsm)',
     'corrected and uncorrected embeddings, so correction can be audited',
     'UMAP distances are not quantitative; never measure on them'),
    # Corrected after the first real run. The level is 'confident', not
    # 'shiv_confident', and there IS a bool column of that name, so the
    # original entry was wrong twice over.
    ('shiv_confident', 'the only viral evidence safe to use; a bool column, '
     "and the matching level of shiv_call is 'confident'",
     'every weaker level is ambient, and even this one is a statement about '
     'the droplet rather than about the cell'),
    ('has_tcr', 'whether a cell has a recovered TCR',
     'off-lineage TCRs are ambient or doublets, so read it against the '
     'lineage, never alone'),
    ('chain_pairing / receptor_subtype',
     'scirpy chain configuration and TRA+TRB versus TRG+TRD',
     'NOT a clonotype. There is no clonotype ID column in this object, so '
     'clonal expansion cannot be analysed without re-running the scirpy '
     'clonotype step'),
]
rows = [{'item': c, 'use_for': u, 'limitation': l} for c, u, l in SHORT]
SL = pd.DataFrame(rows)


def locate(item):
    """Say WHERE a shortlist item lives, rather than guessing a column name.
    An item can be an obs column, a layer, an obsm key, or a level of a
    categorical column, and reporting 'present: False' for a level that does
    exist as a level would be misleading."""
    name = item.split(' ')[0].strip("'")
    # The '==' form is checked FIRST. "shiv_call == 'shiv_confident'" starts
    # with a name that IS an obs column, so a column test placed earlier would
    # return 'obs column' and never count the level, which is the whole point
    # of listing it this way.
    if '==' in item:
        col, _, lvl = (p.strip().strip("'") for p in item.partition('=='))
        if col in OBS.columns:
            n = int((OBS[col].astype(str) == lvl).sum())
            return (f'level of obs["{col}"], {n:,} cells' if n
                    else f'level of obs["{col}"] NOT PRESENT')
        return f'obs["{col}"] missing'
    if name in OBS.columns:
        return 'obs column'
    for slot, keys in (('layer', layer_keys), ('obsm', obsm_keys)):
        if any(name.split('/')[0].strip() in k for k in keys):
            return slot
    if '/' in item:   # an "A / B" pair
        parts = [p.strip() for p in item.split('/')]
        hits = [p for p in parts if p in OBS.columns]
        if hits:
            return 'obs column(s): ' + ', '.join(hits)
    return 'NOT FOUND — check before using'


SL['found_in'] = [locate(i) for i in SL.item]
emit(SL[['item', 'found_in', 'use_for', 'limitation']], 'column_shortlist.csv',
     index=False)
missing_short = SL[SL['found_in'].str.contains('NOT')]
if len(missing_short):
    note(f'{len(missing_short)} shortlist item(s) were not found in this '
         'object. Either 05 named them differently or they were never '
         'written, and either way a script that assumes them will fail.')


# %% ============================================================
# CELL 14 — WRITE THE REPORT
# ============================================================
head('13. Done', 2)

REPORT.append('')
REPORT.append('---')
REPORT.append('')
REPORT.append('Generated by `07_object_report.py`. No figures by design: the '
              'plot set is chosen separately and automated once agreed.')

with open(REPORT_MD, 'w') as fh:
    fh.write('\n'.join(REPORT) + '\n')

csvs = sorted(f for f in os.listdir(OUT_DIR) if f.endswith('.csv'))
print(f'  wrote {os.path.basename(REPORT_MD)} '
      f'({len(REPORT):,} lines) and {len(csvs)} CSV(s)', flush=True)
print(f'  outputs: {OUT_DIR}', flush=True)
print('\n  Read 07_object_report.md top to bottom. Sections 5, 6 and 11 are '
      'the ones that decide what is worth plotting.', flush=True)
print('\n  No figures were produced. That is intentional.', flush=True)
