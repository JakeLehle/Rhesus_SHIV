#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Rhesus SHIV - ORA gene list backfill
=====================================
Rebuilds annotation_output/deg_comparisons/ora_gene_lists/ from the per-contrast
DEG tables that are already on disk. Runs in seconds. Recomputes NOTHING.

WHY THIS EXISTS
---------------
ora_gene_lists/ is empty because the ORA-list writer never executed. It lives
inside the per-contrast test loop in shiv_deg_comparisons.py, and the last run
of that script was the MANIFEST_ONLY = True pass, which builds the manifest and
exits before the loop. CELL 1 creates every output directory up front, which is
why the folder exists but has nothing in it. The tables/ and top50/ content is
from the earlier full run, which predates the ORA block.

So there are two ways forward and this is the cheap one:
  (a) set MANIFEST_ONLY = False and re-run the whole DEG script. Correct, but it
      recomputes every pseudobulk and every cell-level Wilcoxon across 89k cells
      for hours to regenerate tables that already exist unchanged.
  (b) this script. The ORA block is a pure function of the per-contrast table:
      threshold on padj or pval, gate on |log2FC|, split by sign, write. The
      tables on disk already carry gene, log2FC, stat, pval, padj, already
      dropna'd on pval and already sorted, which is exactly the state the
      in-pipeline block sees. Same input, same code, same output.

Use (a) if you ever change a statistic. Use (b) when only the lists are missing.

WHAT IT WRITES  (identical names to the in-pipeline block, so the prep script
                 finds them with no change to GATE or LAYER)
    ora_gene_lists/{stem}_{layer}_{FDR|rawP}_{up|down}.csv
    ora_gene_lists/{stem}_{layer}_background.csv

VERIFICATION
------------
CELL 2 re-derives the thresholds from shiv_deg_comparisons.py rather than
hard-coding them, so the two can never silently drift. CELL 5 cross-checks the
directional counts against deg_counts_by_contrast.csv wherever that file carries
the matching columns, and reports any contrast where they disagree. If the
counts file is from the older run and has no such columns it says so instead of
pretending it checked.

Read-only w.r.t. tables/ and top50/. Spyder cells (# %%).

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

import numpy as np
import pandas as pd

WORKING_DIR = '/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV'
OUT_DIR = os.path.join(WORKING_DIR, 'annotation_output', 'deg_comparisons')
TAB_DIR = os.path.join(OUT_DIR, 'tables')
ORA_DIR = os.path.join(OUT_DIR, 'ora_gene_lists')
COUNTS_CSV = os.path.join(OUT_DIR, 'deg_counts_by_contrast.csv')

# The DEG script is the single source of truth for these. CELL 2 reads them out
# of it; the values here are only the fallback if that file is not beside us.
DEG_SCRIPT_CANDIDATES = [
    os.path.join(WORKING_DIR, 'SHIV_DEG_comparison.py'),
    os.path.join(WORKING_DIR, 'shiv_deg_comparisons.py'),
]
FALLBACK = {'PRIMARY_ALPHA': 0.05, 'PRIMARY_LFC': 1.5, 'MIN_GENES_FOR_ORA': 10}

OVERWRITE = True      # False = refuse to touch a non-empty ora_gene_lists

pd.set_option('display.width', 220)
pd.set_option('display.max_columns', 40)
pd.set_option('display.max_rows', 300)

_report = []


def log(msg=''):
    print(msg, flush=True)
    _report.append(str(msg))


def hr(title, ch='='):
    log('')
    log(ch * 90)
    log(f'{title}')
    log(ch * 90)


# %% ============================================================
# CELL 2 - recover the thresholds from the DEG script
# ============================================================
hr('CELL 2 - thresholds')


def scrape_constants(path, names):
    """Pull `NAME = <number>` out of the DEG script. Regex, not import, because
    importing it would execute the whole pipeline."""
    out = {}
    try:
        with open(path) as fh:
            src = fh.read()
    except OSError:
        return out
    for n in names:
        m = re.search(rf'^{n}\s*=\s*([0-9.]+)', src, re.M)
        if m:
            out[n] = float(m.group(1))
    return out


PARAMS = dict(FALLBACK)
src_used = 'fallback defaults in this script'
for cand in DEG_SCRIPT_CANDIDATES:
    if os.path.exists(cand):
        got = scrape_constants(cand, list(FALLBACK))
        if len(got) == len(FALLBACK):
            PARAMS.update(got)
            src_used = os.path.basename(cand)
            break

PRIMARY_ALPHA = float(PARAMS['PRIMARY_ALPHA'])
PRIMARY_LFC = float(PARAMS['PRIMARY_LFC'])
MIN_GENES_FOR_ORA = int(PARAMS['MIN_GENES_FOR_ORA'])

log(f'  source              : {src_used}')
log(f'  PRIMARY_ALPHA       : {PRIMARY_ALPHA}')
log(f'  PRIMARY_LFC         : {PRIMARY_LFC}')
log(f'  MIN_GENES_FOR_ORA   : {MIN_GENES_FOR_ORA}')

if src_used.startswith('fallback'):
    log('  [WARN] could not read the DEG script, so these are this script\'s own')
    log('         defaults. Confirm they match before trusting the lists.')

for d in DEG_SCRIPT_CANDIDATES:
    if os.path.exists(d):
        with open(d) as fh:
            mo = re.search(r'^MANIFEST_ONLY\s*=\s*(\w+)', fh.read(), re.M)
        if mo:
            log(f'  (for the record, {os.path.basename(d)} currently has '
                f'MANIFEST_ONLY = {mo.group(1)})')
        break


# %% ============================================================
# CELL 3 - inventory the tables
# ============================================================
hr('CELL 3 - input tables')

if not os.path.isdir(TAB_DIR):
    raise RuntimeError(f'No tables directory at {TAB_DIR}. '
                       'The DEG run has not produced per-contrast tables.')

tables = sorted(glob.glob(os.path.join(TAB_DIR, '*.csv')))
log(f'  tables dir          : {TAB_DIR}')
log(f'  csv files found     : {len(tables)}')

if not tables:
    raise RuntimeError('tables/ is empty. Run the DEG script with '
                       'MANIFEST_ONLY = False first.')

LAYER_RE = re.compile(r'^(?P<stem>.+)_(?P<layer>pseudobulk|wilcoxon)\.csv$')

parsed, unparsed = [], []
for t in tables:
    m = LAYER_RE.match(os.path.basename(t))
    if m:
        parsed.append((t, m.group('stem'), m.group('layer')))
    else:
        unparsed.append(os.path.basename(t))

log(f'  parsed              : {len(parsed)}')
if unparsed:
    log(f'  [WARN] {len(unparsed)} file(s) did not match '
        f'<stem>_<layer>.csv and are skipped:')
    for u in unparsed[:10]:
        log(f'           {u}')

by_layer = pd.Series([p[2] for p in parsed]).value_counts()
for k, v in by_layer.items():
    log(f'    {k:12s} {v:4d} contrast table(s)')

n_stems = len({p[1] for p in parsed})
log(f'  distinct contrasts  : {n_stems}')

os.makedirs(ORA_DIR, exist_ok=True)
existing = glob.glob(os.path.join(ORA_DIR, '*.csv'))
log(f'  ora_gene_lists      : {ORA_DIR}')
log(f'  already present     : {len(existing)} file(s)')
if existing and not OVERWRITE:
    raise RuntimeError('ora_gene_lists is not empty and OVERWRITE is False.')


# %% ============================================================
# CELL 4 - write the lists
# ============================================================
hr('CELL 4 - write ORA gene lists')

REQUIRED = ['gene', 'log2FC', 'pval', 'padj']
GATES = [('rawP', 'pval'), ('FDR', 'padj')]

rows = []
n_written = n_bg = 0
bad_cols = []

for path, stem, layer in parsed:
    df = pd.read_csv(path)

    missing = [c for c in REQUIRED if c not in df.columns]
    if missing:
        bad_cols.append((os.path.basename(path), missing))
        continue

    # Mirror the in-pipeline state exactly. The tables were written post-dropna
    # and post-sort, but re-applying dropna is free and makes this script safe
    # against a table written by any earlier revision.
    df = df.dropna(subset=['pval']).copy()

    rec = {'stem': stem, 'layer': layer, 'n_genes_tested': len(df)}

    lfc_ok = df['log2FC'].abs() >= PRIMARY_LFC
    for gate, col in GATES:
        sel = df[(df[col] < PRIMARY_ALPHA) & lfc_ok]
        up = sel[sel['log2FC'] > 0]
        dn = sel[sel['log2FC'] < 0]
        rec[f'{gate}_up'] = len(up)
        rec[f'{gate}_down'] = len(dn)
        rec[f'{gate}_ORA_ready'] = bool(len(up) >= MIN_GENES_FOR_ORA and
                                        len(dn) >= MIN_GENES_FOR_ORA)
        for direction, sub in [('up', up), ('down', dn)]:
            if len(sub) >= MIN_GENES_FOR_ORA:
                sub[REQUIRED].to_csv(
                    os.path.join(ORA_DIR,
                                 f'{stem}_{layer}_{gate}_{direction}.csv'),
                    index=False)
                n_written += 1

    # The background is this contrast's own tested gene set. Always written,
    # even when no direction clears the floor, because the enrichment needs it
    # and a missing background is what silently turns an ORA genome-wide.
    pd.DataFrame({'gene': df['gene']}).to_csv(
        os.path.join(ORA_DIR, f'{stem}_{layer}_background.csv'), index=False)
    n_bg += 1

    rows.append(rec)

if bad_cols:
    log(f'  [ERROR] {len(bad_cols)} table(s) missing required columns:')
    for f, m in bad_cols[:10]:
        log(f'            {f}: missing {m}')

inv = pd.DataFrame(rows)
log(f'  gene list files written : {n_written}')
log(f'  background files written: {n_bg}')
log(f'  total in ora_gene_lists : '
    f'{len(glob.glob(os.path.join(ORA_DIR, "*.csv")))}')

if len(inv):
    hr('enrichable contrasts (both directions clear the floor)', '-')
    for gate, _ in GATES:
        sub = inv.groupby('layer')[f'{gate}_ORA_ready'].agg(['sum', 'count'])
        sub.columns = ['enrichable', 'total']
        log(f'\n  gate = {gate}')
        log(sub.to_string())

    hr('the FDR / pseudobulk set, which is what the prep script consumes', '-')
    pb = inv[(inv['layer'] == 'pseudobulk') & inv['FDR_ORA_ready']]
    log(f'  {len(pb)} contrast(s) enrichable\n')
    if len(pb):
        show = pb[['stem', 'FDR_up', 'FDR_down', 'n_genes_tested']] \
            .sort_values('FDR_up', ascending=False)
        log(show.to_string(index=False))
    else:
        log('  NONE. Check PRIMARY_LFC and PRIMARY_ALPHA above, and look at')
        log('  the rawP column: if rawP is healthy and FDR is empty, that is')
        log('  the multiple-testing burden, not a bug in this script.')

    inv.to_csv(os.path.join(OUT_DIR, 'ora_list_inventory.csv'), index=False)
    log(f'\n  wrote {os.path.join(OUT_DIR, "ora_list_inventory.csv")}')


# %% ============================================================
# CELL 5 - cross-check against the counts the DEG run recorded
# ============================================================
hr('CELL 5 - cross-check vs deg_counts_by_contrast.csv')

if not os.path.exists(COUNTS_CSV):
    log('  [SKIP] deg_counts_by_contrast.csv not found. Nothing to check against.')
elif not len(inv):
    log('  [SKIP] no inventory rows to check.')
else:
    cnt = pd.read_csv(COUNTS_CSV)
    pt = f'a{PRIMARY_ALPHA}_lfc{PRIMARY_LFC}'
    want = {'raw': ('rawP', f'raw_{pt}'), 'fdr': ('FDR', f'fdr_{pt}')}

    have = {k: v for k, v in want.items()
            if f'{v[1]}_up' in cnt.columns and f'{v[1]}_down' in cnt.columns}

    if not have:
        log(f'  [SKIP] that file carries no *_{pt}_up / *_{pt}_down columns.')
        log('         It is from the run that predates the ORA block, so there')
        log('         is nothing to compare. This is expected, not a failure.')
        log('         Available count columns, for reference:')
        cc = [c for c in cnt.columns if c.startswith(('raw_a', 'fdr_a'))]
        log(f'           {cc[:8]}{" ..." if len(cc) > 8 else ""}')
    else:
        # The counts file keys on contrast_id + evidence_layer; the stem starts
        # with the contrast_id, so recover it from there.
        inv2 = inv.copy()
        inv2['contrast_id'] = inv2['stem'].str.split('_').str[0]
        key = ['contrast_id', 'evidence_layer']
        mg = inv2.merge(cnt, left_on=['contrast_id', 'layer'], right_on=key,
                        how='inner', suffixes=('', '_ref'))
        log(f'  matched {len(mg)} of {len(inv2)} contrast x layer rows')

        any_bad = False
        for gate_key, (gate, pref) in have.items():
            for d in ('up', 'down'):
                a = mg[f'{gate}_{d}'].astype(int)
                b = mg[f'{pref}_{d}'].astype(int)
                bad = mg[a != b]
                if len(bad):
                    any_bad = True
                    log(f'  [MISMATCH] {gate} {d}: {len(bad)} row(s) disagree')
                    cols = ['contrast_id', 'layer', f'{gate}_{d}', f'{pref}_{d}']
                    log(bad[cols].head(10).to_string(index=False))
                else:
                    log(f'  [OK] {gate} {d}: all {len(mg)} rows match')
        if not any_bad and len(mg):
            log('\n  The backfilled lists reproduce the counts the DEG run')
            log('  recorded. The lists are the same ones the pipeline would')
            log('  have written.')


# %% ============================================================
# CELL 6 - report
# ============================================================
hr('CELL 6 - done')
log('  Next: python SHIV_enrichment_prep.py')
log('  It looks for *_FDR_up.csv / *_FDR_down.csv plus *_background.csv')
log('  under ora_gene_lists, which now exist.')

rp = os.path.join(OUT_DIR, 'ora_list_backfill_report.txt')
with open(rp, 'w') as fh:
    fh.write('\n'.join(_report))
print(f'\n  report -> {rp}')
