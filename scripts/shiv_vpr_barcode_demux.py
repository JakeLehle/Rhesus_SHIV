#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Rhesus SHIV — vpr Barcode Feature-Demux + Reservoir Map (Read-Only, Ceiling Test)
================================================================================
Stops treating the viral barcode as an alignment problem and treats it as a
feature-barcode problem. The 34bp cassette is engineered sequence: an MluI scar
plus a constant inner flank on each side, bracketing the 10bp core. That flank is
invariant across every quasispecies and every A3-hypermutated provirus, so we
anchor on the flank and let the divergent viral bases soft-clip / fall outside.
No genome alignment of the viral part happens anywhere in here.

This is the CEILING MEASUREMENT, not the production build. It answers one number:
how many clean, object-linked, plausible-target-cell barcodes can we recover from
the complete raw read pool, above the negative floor. That number decides whether
to productionize (cassette-only STARsolo) or to stop and recommend targeted VpxF1
amplicon enrichment.

Why raw FASTQ, not the unmapped pool
------------------------------------
The STARsolo run kept no Unmapped.out.mate (not launched with --outReadsUnmapped),
so unmapped reads were never written as FASTQ. Raw R2 is the complete superset: it
holds the unmapped reads AND the partially-mapped reads that clipped the cassette.

Evidence funnel (the confidence gate)
-------------------------------------
  TIER 1  full anchor, both invariant flanks bracket the core  -> TRUST the core
          even if it is not in Binhua (recovers edited / novel cores).
  TIER 2  one flank only                                       -> REQUIRE an exact
          Binhua core match (the single flank is weaker; GT supplies the 2nd gate).
  Negatives A/B/C/G run at full depth as the false-positive floor.

Output resolution: per animal x timepoint x tissue x FINAL cell type, plus a
per-animal-per-timepoint barcode-diversity table (the expand/contract signal) and
a Binhua coverage stat (of the known barcodes per animal, how many we recovered).

Read-only. Nothing on disk is modified. Spyder cells (# %%); runs under SLURM.

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
RAW_BASE     = '/master/zwallis/WORKING/SC/10X_RAW'
WHITELIST    = os.path.join(STARSOLO_DIR, 'ref', '3M-5pgex-jan-2023.txt')
BINHUA_XLSX  = os.path.join(WORKING_DIR, 'Binhua_FINAL_Barcode_4726.xlsx')
OBJECT       = os.path.join(STARSOLO_DIR, 'merged', 'shiv_host_final_embedded.h5ad')
OUT_DIR      = os.path.join(WORKING_DIR, 'annotation_output', 'vpr_barcode_demux')

LIB_NAME = {
    'A': 'SHIV_Pre_PBMC',             'B': 'SHIV_Pre_LN',
    'C': 'SHIV_Pre_PBMC_IndianMixed', 'D': 'SHIV_21DPI_Pre_LN_Mixed',
    'E': 'SHIV_21DPI_PBMC',           'F': 'SHIV_21DPI_NonInf_LN_Mixed',
    'G': 'SHIV_Necropsy_PBMC',        'H': 'SHIV_Necropsy_LN',
}
# Full run: all 8 (negatives are the floor). Smoke test: set LIBS=['E','H'] and
# MAX_READS_PER_LIB=10_000_000 first to sanity-check before the full pass.
LIBS = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
INFECTED = {'D', 'E', 'F', 'H'}

MAX_READS_PER_LIB = None        # None = full scan; int = cap for a quick test
PROGRESS_EVERY    = 25_000_000

# 10x 5' GEM-X read structure
CB_LEN, UMI_LEN = 16, 12

# barcode geometry
BC_LEN        = 34
ANCHOR_RE     = re.compile('ACGCGC.{22}GCGCGT')
INNER_ORIENT1 = 'AGCGTA'
INNER_ORIENT2 = 'TAGCTG'
PREFILTER = (b'ACGCGC', b'GCGCGT')   # either MluI scar, either strand

# obs columns (verified against shiv_host_final_embedded)
COL_LIBRARY, COL_ANIMAL, COL_TISSUE = 'library', 'animal_id', 'tissue'
COL_TIMEPOINT, COL_CELLTYPE, COL_TSUB = 'timepoint', 'tier1_celltype', 't_subtype'
COL_SHIVPOS, COL_SHIVUMI, COL_SUBSP = 'shiv_pos', 'shiv_umi', 'subspecies'

# cell types that are NOT plausible SHIV targets (ambient/enucleate/unresolved)
NON_TARGET = {'Erythrocytes', 'Platelets', 'Unassigned', 'Mast cells', '', 'nan'}

os.makedirs(OUT_DIR, exist_ok=True)
pd.set_option('display.width', 200)
pd.set_option('display.max_columns', 60)
pd.set_option('display.max_rows', 300)

_COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


def revcomp(s):
    return s.translate(_COMP)[::-1]


def canon_core_from_anchor(bc34):
    if len(bc34) != BC_LEN:
        return None
    inner, core = bc34[6:12], bc34[12:22]
    if inner == INNER_ORIENT1:
        return core
    if inner == INNER_ORIENT2:
        return revcomp(core)
    return None


def norm_cb(cb):
    if not cb:
        return ''
    m = re.search(r'[ACGT]{16}', cb)
    return m.group(0) if m else ''


def raw_cb_from_name(name):
    m = re.search(r'[ACGT]{16}', name)
    return m.group(0) if m else None


def hr(t):
    print('\n' + '=' * 84); print(t); print('=' * 84)


def save_csv(df, name, index=False):
    p = os.path.join(OUT_DIR, name)
    df.to_csv(p, index=index)
    print(f'    saved: {name}  ({len(df)} rows)')


# %% ============================================================
# CELL 2 — Binhua ground truth + flank probes
# ============================================================
hr('CELL 2 — Binhua ground truth cores + cassette flank probes')

gt_core = defaultdict(set)
prefix_counts, suffix_counts = Counter(), Counter()
have_gt = False

if os.path.exists(BINHUA_XLSX):
    try:
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
            print(f'  {animal}: {len(gt_core[animal])} canonical cores')
        wb.close()
        have_gt = True
    except Exception as e:
        print(f'  WARN: Binhua parse failed ({e})')
else:
    print(f'  WARN: Binhua not found ({BINHUA_XLSX}); Tier 2 disabled, Tier 1 novel cores only.')

gt_core_all = set().union(*gt_core.values()) if gt_core else set()
core_to_animals = defaultdict(set)
for a, cores in gt_core.items():
    for c in cores:
        core_to_animals[c].add(a)
print(f'  pooled unique cores: {len(gt_core_all)}')

FLANK5, FLANK3 = {}, {}
for seq, n in prefix_counts.most_common():
    if seq.startswith('ACGCGC') and n >= 20:
        FLANK5[seq] = 'asis' if seq[6:12] == INNER_ORIENT1 else ('rc' if seq[6:12] == INNER_ORIENT2 else None)
for seq, n in suffix_counts.most_common():
    if seq.endswith('GCGCGT') and n >= 20:
        eq = revcomp(seq[0:6])
        FLANK3[seq] = 'asis' if eq == INNER_ORIENT1 else ('rc' if eq == INNER_ORIENT2 else None)
FLANK5 = {k: v for k, v in FLANK5.items() if v}
FLANK3 = {k: v for k, v in FLANK3.items() if v}
print(f'  Tier-2 5-prime flanks: {FLANK5}')
print(f'  Tier-2 3-prime flanks: {FLANK3}')


def match_gt(cand10):
    if cand10 in gt_core_all:
        return cand10
    rc = revcomp(cand10)
    if rc in gt_core_all:
        return rc
    return None


def extract_core(seq):
    """(core, tier, gt) with the evidence funnel:
       tier 1 full anchor -> trusted (GT or novel); tier 2 one flank -> GT-gated."""
    for s in (seq, revcomp(seq)):
        m = ANCHOR_RE.search(s)
        if m:
            c = canon_core_from_anchor(m.group(0))
            if c:
                g = match_gt(c)
                return (g or c, 1, g is not None)
    for s in (seq, revcomp(seq)):
        for fl, rule in FLANK5.items():
            i = s.find(fl)
            if i != -1:
                raw = s[i + 12:i + 22]
                if len(raw) == 10 and set(raw) <= set('ACGT'):
                    g = match_gt(raw if rule == 'asis' else revcomp(raw))
                    if g:
                        return (g, 2, True)
        for fl, rule in FLANK3.items():
            i = s.find(fl)
            if i != -1:
                raw = s[i - 10:i]
                if len(raw) == 10 and set(raw) <= set('ACGT'):
                    g = match_gt(raw if rule == 'asis' else revcomp(raw))
                    if g:
                        return (g, 2, True)
    return (None, 0, False)


# %% ============================================================
# CELL 3 — whitelist + object lookup (final cell type)
# ============================================================
hr('CELL 3 — whitelist + embedded object lookup')

whitelist = set()
if os.path.exists(WHITELIST):
    with open(WHITELIST) as fh:
        for line in fh:
            b = line.strip()
            if b:
                whitelist.add(b)
    print(f'  whitelist: {len(whitelist):,} barcodes')
else:
    print(f'  *** whitelist not found ({WHITELIST}); CB correction disabled.')


def correct_cb(cb):
    """Exact, then 1-mismatch, against the whitelist. Run only on extraction hits."""
    if not whitelist:
        return cb
    if cb in whitelist:
        return cb
    for i in range(len(cb)):
        for b in 'ACGT':
            if b != cb[i]:
                alt = cb[:i] + b + cb[i + 1:]
                if alt in whitelist:
                    return alt
    return None


lookup, letter_to_libval, extra_present = {}, {}, {}
if os.path.exists(OBJECT):
    import scanpy as sc
    adata = sc.read_h5ad(OBJECT, backed='r')
    obs = adata.obs
    print(f'  n_obs: {adata.n_obs:,}')
    for c in [COL_LIBRARY, COL_ANIMAL, COL_TISSUE, COL_TIMEPOINT, COL_CELLTYPE,
              COL_TSUB, COL_SHIVPOS, COL_SHIVUMI, COL_SUBSP]:
        extra_present[c] = c in obs.columns
        if not extra_present[c]:
            print(f'  WARN: obs column missing: {c}')

    libvals = set(str(v) for v in obs[COL_LIBRARY].astype(str).unique()) if extra_present[COL_LIBRARY] else set()
    for L in LIBS:
        letter_to_libval[L] = L if L in libvals else (LIB_NAME[L] if LIB_NAME[L] in libvals else None)

    def final_ct(ct, ts):
        ct = str(ct)
        ts = str(ts)
        if ct == 'T cells' and ts not in ('nan', 'not T', 'Unassigned', '', 'None'):
            return ts
        return ct

    libcol = COL_LIBRARY
    for name, r in obs.iterrows():
        cb = raw_cb_from_name(name)
        if cb is None:
            continue
        ct = r[COL_CELLTYPE] if extra_present[COL_CELLTYPE] else None
        ts = r[COL_TSUB] if extra_present[COL_TSUB] else None
        lookup[(str(r[libcol]), cb)] = {
            'obs_name': name,
            'animal': r[COL_ANIMAL] if extra_present[COL_ANIMAL] else None,
            'tissue': r[COL_TISSUE] if extra_present[COL_TISSUE] else None,
            'timepoint': r[COL_TIMEPOINT] if extra_present[COL_TIMEPOINT] else None,
            'celltype': ct,
            'final_celltype': final_ct(ct, ts),
            't_subtype': ts,
            'subspecies': r[COL_SUBSP] if extra_present[COL_SUBSP] else None,
            'shiv_pos': r[COL_SHIVPOS] if extra_present[COL_SHIVPOS] else None,
            'shiv_umi': r[COL_SHIVUMI] if extra_present[COL_SHIVUMI] else None,
        }
    print(f'  lookup built: {len(lookup):,} keys')
    print(f'  letter->libval: {letter_to_libval}')
else:
    print(f'  *** object not found ({OBJECT}); linkage skipped, read/cell tables still emitted.')


# %% ============================================================
# CELL 4 — raw FASTQ manifest
# ============================================================
hr('CELL 4 — raw GEX FASTQ manifest')

DECOMP = 'pigz -dc -p 4' if shutil.which('pigz') else 'zcat'
print(f'  decompressor: {DECOMP}')
manifest = []
for L in LIBS:
    d = os.path.join(RAW_BASE, f'GEX_LIBRARY_{L}')
    r1 = sorted(glob.glob(os.path.join(d, '*_R1_*.fastq.gz')))
    r2 = sorted(glob.glob(os.path.join(d, '*_R2_*.fastq.gz')))
    if not r1 or not r2:
        print(f'  {L}: MISSING R1/R2 in {d}; skipped')
        continue
    manifest.append({'lib': L, 'group': 'INF' if L in INFECTED else 'neg', 'r1': r1, 'r2': r2})
    print(f'  {L} [{manifest[-1]["group"]}] R1={len(r1)} R2={len(r2)}  {os.path.basename(r1[0])}')


# %% ============================================================
# CELL 5 — scan raw reads: prefilter -> extract -> CB correct
# ============================================================
hr('CELL 5 — feature-demux scan (this is the heavy step)')


def open_stream(files):
    return subprocess.Popen(DECOMP.split() + files, stdout=subprocess.PIPE, bufsize=1 << 22)


def scan_library(row):
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
            if PREFILTER[0] not in s2 and PREFILTER[1] not in s2:
                continue
            core, tier, gt = extract_core(s2.decode(errors='ignore').strip())
            if tier == 0:
                continue
            cb_raw = s1[:CB_LEN].decode(errors='ignore')
            umi = s1[CB_LEN:CB_LEN + UMI_LEN].decode(errors='ignore')
            cb = correct_cb(cb_raw)
            hits.append({'lib': L, 'group': row['group'], 'tier': tier, 'core': core,
                         'gt': int(gt), 'cb_raw': cb_raw, 'cb': cb or '', 'umi': umi})
    finally:
        for pp in (p1, p2):
            try:
                pp.stdout.close(); pp.wait()
            except Exception:
                pass
    return hits, n


all_hits = []
for row in manifest:
    hits, n = scan_library(row)
    all_hits.extend(hits)
    t1 = sum(h['tier'] == 1 for h in hits)
    t2 = sum(h['tier'] == 2 for h in hits)
    with_cb = sum(bool(h['cb']) for h in hits)
    print(f'  {row["lib"]} [{row["group"]}]: {n:,} reads | {len(hits)} barcoded '
          f'(T1 {t1}, T2 {t2}) | {with_cb} CB-valid')

reads_df = pd.DataFrame(all_hits)
if len(reads_df):
    save_csv(reads_df, 'demux_barcoded_reads.csv')
    neg_gt = reads_df[(reads_df['group'] == 'neg') & (reads_df['gt'] == 1)]
    print(f"\n  specificity floor: {len(neg_gt)} GT core(s) in negatives "
          f"({'CLEAN' if len(neg_gt) == 0 else 'INSPECT — possible contamination'})")


# %% ============================================================
# CELL 6 — dedup molecules/cells + link to object
# ============================================================
hr('CELL 6 — molecules, cells, linkage')

cells_df = pd.DataFrame()
if len(reads_df):
    valid = reads_df[reads_df['cb'] != ''].copy()
    valid['umi_key'] = np.where(valid['umi'].str.len() > 0, valid['umi'],
                                'noUMI_' + valid.index.astype(str))
    mol = valid.drop_duplicates(['lib', 'cb', 'umi_key', 'core'])
    print(f'  barcoded reads (CB-valid): {len(valid)} | unique molecules: {len(mol)}')

    rows, matched = [], 0
    for (lib, cb), g in mol.groupby(['lib', 'cb']):
        cores = sorted(set(g['core']))
        meta = lookup.get((letter_to_libval.get(lib), cb)) if letter_to_libval.get(lib) else None
        if meta:
            matched += 1
        ct = meta['final_celltype'] if meta else ''
        rows.append({
            'lib': lib, 'group': 'INF' if lib in INFECTED else 'neg', 'sample': LIB_NAME[lib],
            'cb': cb, 'cores': ';'.join(cores), 'n_cores': len(cores),
            'n_umis': g['umi_key'].nunique(), 'best_tier': int(g['tier'].min()),
            'in_object': meta is not None,
            'animal': meta['animal'] if meta else '',
            'tissue': meta['tissue'] if meta else '',
            'timepoint': meta['timepoint'] if meta else '',
            'final_celltype': ct,
            't_subtype': meta['t_subtype'] if meta else '',
            'subspecies': meta['subspecies'] if meta else '',
            'shiv_pos': meta['shiv_pos'] if meta else '',
            'shiv_umi': meta['shiv_umi'] if meta else '',
            'is_target': bool(meta) and str(ct) not in NON_TARGET,
            'barcode_animals_binhua': ';'.join(sorted(set().union(
                *[core_to_animals.get(c, set()) for c in cores]))),
        })
    cells_df = pd.DataFrame(rows)
    n = len(cells_df)
    print(f'  barcode-positive cells: {n} | linked to object: {matched} '
          f'({100 * matched / n:.0f}%)' if n else '  no cells')
    if n and matched == 0:
        print('  *** 0% linkage — check letter->libval / obs_names before trusting the map.')
    if n:
        save_csv(cells_df, 'demux_barcode_by_cell.csv')


# %% ============================================================
# CELL 7 — high-resolution reservoir map + diversity + coverage
# ============================================================
hr('CELL 7 — reservoir map (animal x timepoint x tissue x final cell type)')

if len(cells_df):
    inf = cells_df[cells_df['group'] == 'INF']
    target = inf[(inf['in_object']) & (inf['is_target'])]
    print(f'  infected barcode+ cells: {len(inf)} | object-linked target cells: {len(target)}')

    if len(target):
        res = (target.groupby(['animal', 'timepoint', 'tissue', 'final_celltype'])
               .agg(cells=('cb', 'nunique'),
                    distinct_barcodes=('cores', lambda s: len(set(';'.join(s).split(';')))),
                    umis=('n_umis', 'sum')).reset_index()
               .sort_values(['animal', 'timepoint', 'cells'], ascending=[True, True, False]))
        print('\n  RESERVOIR MAP:')
        print(res.to_string(index=False))
        save_csv(res, 'demux_reservoir_map.csv')

        div = (target.groupby(['animal', 'timepoint'])
               .agg(cells=('cb', 'nunique'),
                    distinct_barcodes=('cores', lambda s: len(set(';'.join(s).split(';'))))).reset_index())
        print('\n  BARCODE DIVERSITY over time (the expand/contract signal):')
        print(div.to_string(index=False))
        save_csv(div, 'demux_barcode_diversity_by_time.csv')
    else:
        print('  No object-linked target cells to map (see per-library tally below).')

    # per-library specificity/tally regardless of linkage
    tally = (cells_df.groupby(['lib', 'sample', 'group'])
             .agg(cells=('cb', 'nunique'), linked=('in_object', 'sum'),
                  target=('is_target', 'sum'),
                  distinct_barcodes=('cores', lambda s: len(set(';'.join(s).split(';'))))).reset_index())
    print('\n  per-library tally:')
    print(tally.to_string(index=False))
    save_csv(tally, 'demux_cells_by_library.csv')

    # Binhua coverage: of known barcodes per animal, how many recovered in single cell
    if have_gt and len(target):
        rec = defaultdict(set)
        for _, r in target.iterrows():
            for c in r['cores'].split(';'):
                if c:
                    for a in core_to_animals.get(c, {r['animal']}):
                        rec[a].add(c)
        cov = pd.DataFrame([{'animal': a, 'known_barcodes': len(gt_core[a]),
                             'recovered_in_sc': len(rec.get(a, set())),
                             'pct': round(100 * len(rec.get(a, set())) / max(1, len(gt_core[a])), 2)}
                            for a in sorted(gt_core)])
        print('\n  Binhua coverage (single-cell recovery vs known repertoire):')
        print(cov.to_string(index=False))
        save_csv(cov, 'demux_binhua_coverage.csv')


# %% ============================================================
# CELL 8 — SUMMARY + GO/NO-GO
# ============================================================
hr('CELL 8 — SUMMARY + GO/NO-GO')

n_reads = len(reads_df) if len(reads_df) else 0
n_gt = int(reads_df['gt'].sum()) if len(reads_df) else 0
n_cells = len(cells_df) if len(cells_df) else 0
n_target = int(((cells_df['group'] == 'INF') & cells_df['in_object'] & cells_df['is_target']).sum()) if len(cells_df) else 0
n_neg = int((cells_df['group'] == 'neg').sum()) if len(cells_df) else 0
n_bc = 0
if len(cells_df):
    t = cells_df[(cells_df['group'] == 'INF') & cells_df['in_object'] & cells_df['is_target']]
    n_bc = len(set(';'.join(t['cores']).split(';'))) if len(t) else 0

print(f"""
  Barcoded reads (raw pool)         : {n_reads}   (Tier-1+2 GT: {n_gt})
  Barcode-positive cells (all)      : {n_cells}
  Negative-library cells (floor)    : {n_neg}
  Object-linked TARGET cells (INF)  : {n_target}
  Distinct barcodes in those cells  : {n_bc}
""")

verdict = os.path.join(OUT_DIR, 'GO_NOGO.txt')
with open(verdict, 'w') as fh:
    if n_neg > 0:
        msg = (f"CAUTION: {n_neg} barcode+ cell(s) in negative libraries. The floor is "
               "not clean; inspect before trusting any infected count.")
    elif n_target >= 20:
        msg = (f"GO: {n_target} clean object-linked target cells carrying {n_bc} distinct "
               "barcodes, negative floor clean. Productionize with the cassette-only "
               "STARsolo build and this becomes the reservoir dataset.")
    elif n_target >= 5:
        msg = (f"MARGINAL: {n_target} clean target cells, {n_bc} barcodes, floor clean. "
               "Enough to report as a proof-of-principle reservoir observation but thin "
               "for population dynamics. Weigh productionize vs targeted VpxF1 enrichment.")
    else:
        msg = (f"NO-GO for single-cell reservoir dynamics: only {n_target} clean target "
               "cell(s). The signal is structurally capped by 5-prime capture geometry, not "
               "by our extraction. Honest next step is targeted VpxF1 amplicon enrichment off "
               "the archived cDNA (new sequencing). Bulk barcode dynamics remain available in "
               "the Binhua file.")
    print('  ' + msg)
    fh.write(msg + '\n')

print(f'\n  wrote {verdict}')
print('=' * 84)
print('vpr FEATURE-DEMUX COMPLETE — read-only. Nothing modified.')
print('=' * 84)
