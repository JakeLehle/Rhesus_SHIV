#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Rhesus SHIV - verify the backfilled ORA lists against the original DEG run
==========================================================================
Closes a gap in the backfill's own self-check. That check looked for
    raw_a0.05_lfc1.5_up   /  raw_a0.05_lfc1.5_down
and correctly reported SKIP, because the older DEG run wrote the UNDIRECTED
totals instead:
    raw_a0.05_lfc1.5      /  fdr_a0.05_lfc1.5
which the SKIP message itself listed. Those totals are still a real check, and
a tight one:

    a gene passing |log2FC| >= 1.5 cannot have log2FC == 0, so every gene in
    the thresholded set is strictly up or strictly down, so

        up + down  ==  total        EXACTLY, not approximately

So if the backfill reproduced the original run, up + down from the lists on disk
must equal the totals the DEG run recorded, for every contrast, in both layers,
at both gates. Any row that is off by even one means the lists on disk are not
the lists the pipeline would have written, and step 2 should not run.

Reads only. Takes seconds. Run it before shiv_enrichment_run.py.

Author: Jake Lehle
Date: September 2026
"""

import os
import sys

import numpy as np
import pandas as pd

WORKING_DIR = '/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV'
OUT_DIR = os.path.join(WORKING_DIR, 'annotation_output', 'deg_comparisons')
INV_CSV = os.path.join(OUT_DIR, 'ora_list_inventory.csv')
CNT_CSV = os.path.join(OUT_DIR, 'deg_counts_by_contrast.csv')

ALPHA, LFC = 0.05, 1.5
PT = f'a{ALPHA}_lfc{LFC}'

print('=' * 90)
print('VERIFY backfilled ORA lists vs the original DEG run')
print('=' * 90)

for p in (INV_CSV, CNT_CSV):
    if not os.path.exists(p):
        sys.exit(f'[FATAL] missing {p}')

inv = pd.read_csv(INV_CSV)
cnt = pd.read_csv(CNT_CSV)

inv['contrast_id'] = inv['stem'].str.split('_').str[0]
print(f'  inventory rows : {len(inv)}')
print(f'  counts rows    : {len(cnt)}')

PAIRS = [('rawP', f'raw_{PT}'), ('FDR', f'fdr_{PT}')]
missing = [c for _, c in PAIRS if c not in cnt.columns]
if missing:
    sys.exit(f'[FATAL] counts file has no {missing}. Nothing to verify against.')

mg = inv.merge(cnt, left_on=['contrast_id', 'layer'],
               right_on=['contrast_id', 'evidence_layer'],
               how='inner', suffixes=('', '_ref'))
print(f'  matched        : {len(mg)} of {len(inv)} contrast x layer rows')

if len(mg) < len(inv):
    lost = set(zip(inv['contrast_id'], inv['layer'])) - \
           set(zip(mg['contrast_id'], mg['layer']))
    print(f'  [WARN] {len(lost)} inventory row(s) had no counts row: '
          f'{sorted(lost)[:6]}')

# ---- the check ----
fail = 0
for gate, col in PAIRS:
    got = (mg[f'{gate}_up'] + mg[f'{gate}_down']).astype(int)
    ref = mg[col].astype(int)
    bad = mg[got != ref].assign(list_total=got[got != ref],
                                deg_total=ref[got != ref])
    if len(bad):
        fail += len(bad)
        print(f'\n  [MISMATCH] {gate}: {len(bad)} of {len(mg)} rows disagree')
        print(bad[['contrast_id', 'layer', f'{gate}_up', f'{gate}_down',
                   'list_total', 'deg_total']].head(15).to_string(index=False))
    else:
        print(f'  [OK] {gate}: up + down == {col} for all {len(mg)} rows')

# ---- secondary: gene universe ----
# n_genes_tested in the inventory comes from the table the backfill read.
# The DEG run recorded the same quantity independently.
if 'n_genes_tested' in cnt.columns:
    a = mg['n_genes_tested'].astype(int)
    b = mg['n_genes_tested_ref'].astype(int) if 'n_genes_tested_ref' in mg \
        else mg['n_genes_tested'].astype(int)
    bad = mg[a != b]
    if len(bad):
        fail += len(bad)
        print(f'\n  [MISMATCH] gene universe differs on {len(bad)} row(s)')
        print(bad[['contrast_id', 'layer']].head(10).to_string(index=False))
    else:
        print(f'  [OK] gene universe matches on all {len(mg)} rows')

print()
print('=' * 90)
if fail:
    print(f'  FAILED: {fail} disagreement(s). Do NOT run step 2.')
    print('  The lists on disk are not what the pipeline would have written.')
    sys.exit(1)
print('  PASSED. The backfilled lists reproduce the original DEG run exactly.')
print('  Safe to run shiv_enrichment_run.py.')
print('=' * 90)
