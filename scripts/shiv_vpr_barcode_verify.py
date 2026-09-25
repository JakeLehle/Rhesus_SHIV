#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Rhesus SHIV — vpr Barcode Recovery VERIFICATION (fuzzy re-scan + locus census)
=============================================================================
Purpose: before closing the single-cell barcode analysis, prove that EXACT flank
matching was not silently under-detecting. This re-scans the raw pool with a
mismatch-tolerant, seed-anchored flank match and reports how the recovered-read
count moves as we relax from 0 to MAX_MM mismatches, split infected vs negative.

The logic of the test
----------------------
  - If exact matching was the bottleneck, the count JUMPS as we allow mismatches
    -> the reads were there and we were missing them -> DO NOT close the analysis.
  - If the count stays flat in single digits and the negative floor stays clean,
    then exact matching was NOT the bottleneck -> the reads genuinely are not in
    the data -> the no-go is confirmed. Relaxing also inflates noise in the
    negatives IF we were matching noise, so this closes the noise concern too.

Also runs a cassette-locus coverage census off the mapped BAMs: the hard ceiling
on how many reads physically reach position ~5900, which bounds EVERY method
(alignment, surgery, demux) at once.

Seed-anchored fuzzy match keeps it fast: a true cassette read carries up to four
exact 6-mer seeds (two MluI scars + two inner flanks); we anchor the 34bp window
on any intact seed and score mismatches only at the fixed flank positions, so a
read is caught unless ALL its seeds are corrupted (negligible and untrustworthy).

Read-only. Dumps every hit's FULL read sequence so hits can be eyeballed/BLASTed.

Author: Jake Lehle / Kaushal Lab
Date: July 2026
"""

# %% ============================================================
# CELL 1 — CONFIG + HELPERS
# ============================================================
import os
import re
import glob
import subprocess
import shutil
from collections import Counter, defaultdict

import numpy as np
import pandas as pd

WORKING_DIR  = '/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV'
STARSOLO_DIR = os.path.join(WORKING_DIR, 'shiv_starsolo')
SOLO_DIR     = os.path.join(STARSOLO_DIR, 'solo')
RAW_BASE     = '/master/zwallis/WORKING/SC/10X_RAW'
BINHUA_XLSX  = os.path.join(WORKING_DIR, 'Binhua_FINAL_Barcode_4726.xlsx')
OUT_DIR      = os.path.join(WORKING_DIR, 'annotation_output', 'vpr_barcode_verify')

LIB_NAME = {
    'A': 'SHIV_Pre_PBMC',             'B': 'SHIV_Pre_LN',
    'C': 'SHIV_Pre_PBMC_IndianMixed', 'D': 'SHIV_21DPI_Pre_LN_Mixed',
    'E': 'SHIV_21DPI_PBMC',           'F': 'SHIV_21DPI_NonInf_LN_Mixed',
    'G': 'SHIV_Necropsy_PBMC',        'H': 'SHIV_Necropsy_LN',
}
LIBS = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
INFECTED = {'D', 'E', 'F', 'H'}

MAX_MM = 2                       # tolerate up to this many mismatches in the flank
MAX_READS_PER_LIB = None         # None = full; int to cap for a quick test
PROGRESS_EVERY = 25_000_000
CB_LEN, UMI_LEN = 16, 12

BC_LEN = 34
INNER_ORIENT1, INNER_ORIENT2 = 'AGCGTA', 'TAGCTG'

# cassette-locus coverage census (mapped BAMs); primer bracket 5790..6118
CENSUS_START1, CENSUS_END1 = 5791, 6118
CENSUS_MID0 = 5955               # midpoint reads must span to overlap the cassette
SAMTOOLS = os.environ.get('SAMTOOLS') or shutil.which('samtools')

os.makedirs(OUT_DIR, exist_ok=True)
pd.set_option('display.width', 200)
pd.set_option('display.max_columns', 60)

_COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
CIGAR_RE = re.compile(r'(\d+)([MIDNSHP=X])')


def revcomp(s):
    return s.translate(_COMP)[::-1]


def hr(t):
    print('\n' + '=' * 84); print(t); print('=' * 84)


def save_csv(df, name, index=False):
    p = os.path.join(OUT_DIR, name)
    df.to_csv(p, index=index)
    print(f'    saved: {name}  ({len(df)} rows)')


print(f'VERIFY: MAX_MM={MAX_MM}  samtools={SAMTOOLS}  out={OUT_DIR}')


# %% ============================================================
# CELL 2 — ground truth + flank templates + seed map
# ============================================================
hr('CELL 2 — Binhua cores + templates + seed anchors')

gt_core = defaultdict(set)
prefix_counts, suffix_counts = Counter(), Counter()


def canon_core_from_anchor(bc34):
    if len(bc34) != BC_LEN:
        return None
    inner, core = bc34[6:12], bc34[12:22]
    if inner == INNER_ORIENT1:
        return core
    if inner == INNER_ORIENT2:
        return revcomp(core)
    return None


if os.path.exists(BINHUA_XLSX):
    from openpyxl import load_workbook
    wb = load_workbook(BINHUA_XLSX, read_only=True, data_only=True)
    for sh in wb.sheetnames:
        animal = sh.split()[0]
        for row in wb[sh].iter_rows(values_only=True):
            if row and row[0] is not None and str(row[0]).strip().isdigit():
                bc = row[2]
                if isinstance(bc, str):
                    bc = bc.strip().upper()
                    if len(bc) == BC_LEN and set(bc) <= set('ACGTN'):
                        c = canon_core_from_anchor(bc) or canon_core_from_anchor(revcomp(bc))
                        if c:
                            gt_core[animal].add(c)
                        prefix_counts[bc[:12]] += 1
                        suffix_counts[bc[-12:]] += 1
    wb.close()
gt_core_all = set().union(*gt_core.values()) if gt_core else set()
print(f'  pooled cores: {len(gt_core_all)}')

# top constant flanks (data-driven, same as the demux)
top5 = [s for s, n in prefix_counts.most_common() if s.startswith('ACGCGC') and n >= 20][:2]
top3 = [s for s, n in suffix_counts.most_common() if s.endswith('GCGCGT') and n >= 20][:2]
print(f'  5-prime flanks: {top5}')
print(f'  3-prime flanks: {top3}')

# pair 5' and 3' flanks into the two orientation templates by inner flank identity
TEMPLATES = []   # (name, fixed{pos:base}, rule)
for f5 in top5:
    inner = f5[6:12]
    rule = 'asis' if inner == INNER_ORIENT1 else ('rc' if inner == INNER_ORIENT2 else None)
    if rule is None:
        continue
    # 3' partner shares the same orientation: inner2 = revcomp(inner) at its start
    want = revcomp(inner)
    f3 = next((s for s in top3 if s[:6] == want), None)
    if f3 is None:
        continue
    fixed = {i: ch for i, ch in enumerate(f5)}
    fixed.update({22 + j: ch for j, ch in enumerate(f3)})
    TEMPLATES.append((f'orient_{inner}', fixed, rule))

# seed anchors: exact 6-mers at fixed offsets 0,6,22,28 of each template
SEEDS = []   # (seed, template_index, offset)
for ti, (name, fixed, rule) in enumerate(TEMPLATES):
    tstr = ['?'] * BC_LEN
    for pos, ch in fixed.items():
        tstr[pos] = ch
    tstr = ''.join(tstr)
    for off in (0, 6, 22, 28):
        seed = tstr[off:off + 6]
        if '?' not in seed:
            SEEDS.append((seed, ti, off))
SEED_SET = sorted(set(s for s, _, _ in SEEDS))
SEED_SET_B = [s.encode() for s in SEED_SET]
print(f'  templates: {[t[0] for t in TEMPLATES]}')
print(f'  seed anchors ({len(SEEDS)}): {SEED_SET}')


def match_gt(cand10):
    if cand10 in gt_core_all:
        return cand10
    rc = revcomp(cand10)
    if rc in gt_core_all:
        return rc
    return None


# %% ============================================================
# CELL 3 — seed-anchored fuzzy extractor
# ============================================================
def fuzzy_extract(seq):
    """Return (mm, core, gt_bool, template, strand) for the best (fewest-mm) hit
    with mm<=MAX_MM at the fixed flank positions, else None."""
    best = None
    for s, strand in ((seq, '+'), (revcomp(seq), '-')):
        for seed, ti, off in SEEDS:
            start = 0
            while True:
                p = s.find(seed, start)
                if p < 0:
                    break
                start = p + 1
                w0 = p - off
                if w0 < 0 or w0 + BC_LEN > len(s):
                    continue
                window = s[w0:w0 + BC_LEN]
                name, fixed, rule = TEMPLATES[ti]
                mm = 0
                ok = True
                for pos, ch in fixed.items():
                    if window[pos] != ch:
                        mm += 1
                        if mm > MAX_MM:
                            ok = False
                            break
                if not ok:
                    continue
                core = window[12:22]
                if not (len(core) == 10 and set(core) <= set('ACGT')):
                    continue
                cand = core if rule == 'asis' else revcomp(core)
                g = match_gt(cand)
                hit = (mm, g or cand, g is not None, name, strand)
                if best is None or mm < best[0]:
                    best = hit
    return best


# %% ============================================================
# CELL 4 — raw FASTQ manifest
# ============================================================
hr('CELL 4 — raw FASTQ manifest')
DECOMP = 'pigz -dc -p 4' if shutil.which('pigz') else 'zcat'
manifest = []
for L in LIBS:
    d = os.path.join(RAW_BASE, f'GEX_LIBRARY_{L}')
    r1 = sorted(glob.glob(os.path.join(d, '*_R1_*.fastq.gz')))
    r2 = sorted(glob.glob(os.path.join(d, '*_R2_*.fastq.gz')))
    if r1 and r2:
        manifest.append({'lib': L, 'group': 'INF' if L in INFECTED else 'neg', 'r1': r1, 'r2': r2})
    else:
        print(f'  {L}: MISSING; skipped')
print(f'  {len(manifest)} libraries, decompressor {DECOMP}')


# %% ============================================================
# CELL 5 — fuzzy re-scan
# ============================================================
hr('CELL 5 — fuzzy re-scan (the decisive count-vs-mismatch test)')


def open_stream(files):
    return subprocess.Popen(DECOMP.split() + files, stdout=subprocess.PIPE, bufsize=1 << 22)


def scan(row):
    L = row['lib']
    p1, p2 = open_stream(row['r1']), open_stream(row['r2'])
    f1, f2 = p1.stdout, p2.stdout
    hits, n = [], 0
    try:
        while True:
            h1 = f1.readline()
            if not h1:
                break
            s1 = f1.readline(); f1.readline(); f1.readline()
            f2.readline(); s2 = f2.readline(); f2.readline(); f2.readline()
            n += 1
            if MAX_READS_PER_LIB and n > MAX_READS_PER_LIB:
                break
            if n % PROGRESS_EVERY == 0:
                print(f'    {L}: {n:,} reads, {len(hits)} hits')
            if not any(seed in s2 for seed in SEED_SET_B):
                continue
            res = fuzzy_extract(s2.decode(errors='ignore').strip())
            if res is None:
                continue
            mm, core, gt, tmpl, strand = res
            hits.append({'lib': L, 'group': row['group'], 'mm': mm, 'core': core,
                         'gt': int(gt), 'template': tmpl, 'strand': strand,
                         'cb_raw': s1[:CB_LEN].decode(errors='ignore'),
                         'umi': s1[CB_LEN:CB_LEN + UMI_LEN].decode(errors='ignore'),
                         'read_seq': s2.decode(errors='ignore').strip()})
    finally:
        for pp in (p1, p2):
            try:
                pp.stdout.close(); pp.wait()
            except Exception:
                pass
    return hits, n


all_hits = []
for row in manifest:
    hits, n = scan(row)
    all_hits.extend(hits)
    by_mm = Counter(h['mm'] for h in hits)
    gt_n = sum(h['gt'] for h in hits)
    print(f'  {row["lib"]} [{row["group"]}]: {n:,} reads | {len(hits)} hits '
          f'(mm0={by_mm.get(0,0)} mm1={by_mm.get(1,0)} mm2={by_mm.get(2,0)}) | GT {gt_n}')

hits_df = pd.DataFrame(all_hits)
if len(hits_df):
    save_csv(hits_df, 'verify_fuzzy_hits.csv')


# %% ============================================================
# CELL 6 — count-vs-mismatch table (the answer)
# ============================================================
hr('CELL 6 — recovered reads vs mismatch tolerance')

if len(hits_df):
    # cumulative: reads recoverable at each mm ceiling, split group x gt
    rows = []
    for k in range(MAX_MM + 1):
        sub = hits_df[hits_df['mm'] <= k]
        for grp in ('INF', 'neg'):
            g = sub[sub['group'] == grp]
            rows.append({'mm_ceiling': k, 'group': grp, 'total_hits': len(g),
                         'gt_matched': int(g['gt'].sum()),
                         'non_gt': int((g['gt'] == 0).sum())})
    ctab = pd.DataFrame(rows)
    print(ctab.to_string(index=False))
    save_csv(ctab, 'verify_count_vs_mismatch.csv')

    print('\n  READ THIS:')
    inf0 = ctab[(ctab.mm_ceiling == 0) & (ctab.group == 'INF')]['gt_matched'].iloc[0]
    inf2 = ctab[(ctab.mm_ceiling == MAX_MM) & (ctab.group == 'INF')]['gt_matched'].iloc[0]
    neg2 = ctab[(ctab.mm_ceiling == MAX_MM) & (ctab.group == 'neg')]['total_hits'].iloc[0]
    print(f'   - infected GT reads: {inf0} at mm=0 -> {inf2} at mm={MAX_MM}')
    print(f'   - negative hits at mm={MAX_MM} (any): {neg2}')
    if inf2 <= inf0 + 2 and neg2 == 0:
        print('   - VERDICT: relaxing the match did NOT unlock a hidden pool, and the')
        print('     negative floor stayed clean. Exact matching was not the bottleneck.')
        print('     The reads are not in the data. No-go CONFIRMED; safe to close/send.')
    elif inf2 > inf0 + 5:
        print('   - VERDICT: the count JUMPED with relaxed matching. Exact matching WAS')
        print('     under-detecting. DO NOT close yet — investigate the extra reads.')
    else:
        print('   - VERDICT: modest movement. Inspect verify_fuzzy_hits.csv read_seq')
        print('     column before deciding; confirm hits are cassette-in-SHIV-context.')
else:
    print('  No hits at all, even fuzzy. Strongest possible confirmation the reads')
    print('  are not present. No-go confirmed.')


# %% ============================================================
# CELL 7 — cassette-locus coverage census (the hard ceiling for ALL methods)
# ============================================================
hr('CELL 7 — mapped-BAM coverage over the cassette locus')

if SAMTOOLS:
    census = []
    for L in LIBS:
        bam = os.path.join(SOLO_DIR, f'{L}_{LIB_NAME[L]}_Aligned.sortedByCoord.out.bam')
        if not os.path.exists(bam):
            continue
        # max aligned depth across the primer bracket
        region = f'SHIVAD8EO:{CENSUS_START1}-{CENSUS_END1}'
        maxd = 0
        try:
            dp = subprocess.run([SAMTOOLS, 'depth', '-a', '-r', region, bam],
                                capture_output=True, text=True).stdout
            for line in dp.splitlines():
                f = line.split('\t')
                if len(f) >= 3:
                    maxd = max(maxd, int(f[2]))
        except Exception as e:
            print(f'    depth failed {L} ({e})')
        # reads whose alignment spans the cassette midpoint
        spanning = 0
        try:
            p = subprocess.Popen([SAMTOOLS, 'view', bam, region],
                                 stdout=subprocess.PIPE, text=True)
            for line in p.stdout:
                f = line.split('\t')
                if len(f) < 6:
                    continue
                pos1 = int(f[3])
                span = sum(int(n) for n, op in CIGAR_RE.findall(f[5]) if op in 'MDN=X')
                if pos1 <= CENSUS_MID0 + 1 <= pos1 + span - 1:
                    spanning += 1
            p.wait()
        except Exception as e:
            print(f'    view failed {L} ({e})')
        census.append({'lib': L, 'group': 'INF' if L in INFECTED else 'neg',
                       'max_depth_window': maxd, 'reads_spanning_cassette_mid': spanning})
    cdf = pd.DataFrame(census)
    print(cdf.to_string(index=False))
    save_csv(cdf, 'verify_locus_coverage.csv')
    ceil_inf = cdf[cdf.group == 'INF']['reads_spanning_cassette_mid'].sum() if len(cdf) else 0
    print(f'\n  Hard ceiling: {ceil_inf} infected reads physically span the cassette')
    print('  locus across all libraries. No method (alignment, surgery, demux) can')
    print('  recover a barcode from more reads than this, and only those carrying the')
    print('  full 34bp cassette are readable at all.')
else:
    print('  samtools not found; skipping locus census (activate sc_pre).')

hr('VERIFICATION COMPLETE — read-only. Nothing modified.')
