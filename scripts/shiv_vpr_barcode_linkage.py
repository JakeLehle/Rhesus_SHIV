#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Rhesus SHIV — vpr Barcode Extraction + Cell Linkage (Read-Only)
===============================================================
Takes the vpr-window reads we can already see in the mapped BAMs, extracts the
viral barcode core two ways, then links every barcode-positive cell to its
animal / tissue / timepoint / cell type off the embedded object. The output is
the barcode-by-cell table that every downstream viral-population plot is built
on. It is deliberately shaped so that when we later recover more reads from the
unmapped pool, the exact same script scales up without changes.

What it fixes from shiv_vpr_read_diagnostic.py
----------------------------------------------
  - Drops the fake "insertion locus" (no MluI in the window) and the borderline
    verdict it drove. This script does extraction + linkage, not surgery triage.
  - Strand-robust core matching: a read on either strand now canonicalizes to
    the same ground-truth core (the old first-hit logic could hand back the
    reverse-complement core and silently miss a real match).

What it adds
------------
  - TIER 2 one-sided flank extraction: reads carrying only ONE constant 12bp
    cassette flank (the other end clipped off) get their 10bp core read from the
    adjacent bases. GT-gated, so specificity is preserved. Unvalidated Tier 2
    candidates are logged separately, never counted as recovered.
  - Molecule and cell dedup across STARsolo + CellRanger (two alignments of the
    same FASTQs), so counts are true molecules/cells, not double-counted reads.
  - A SELF-VERIFYING join: Cell 3 prints the object's obs_names format and the
    per-library suffix scheme, then Cell 6 reports the match rate of barcoded CBs
    against the object under the chosen key. A zero match rate warns loudly
    instead of masquerading as "no virus-positive cells."

Read-only. No BAM, reference, or object is modified.
Conventions: Spyder cells (# %%), config in Cell 1, warn-and-continue.

Author: Jake Lehle / Kaushal Lab
Date: July 2026
"""

# %% ============================================================
# CELL 1 — CONFIG AND HELPERS
# ============================================================
import os
import re
import glob
import subprocess
import shutil
from collections import Counter, defaultdict

import numpy as np
import pandas as pd

WORKING_DIR    = '/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV'
STARSOLO_DIR   = os.path.join(WORKING_DIR, 'shiv_starsolo')
SOLO_DIR       = os.path.join(STARSOLO_DIR, 'solo')
CELLRANGER_DIR = os.path.join(WORKING_DIR, 'Analysis')
BINHUA_XLSX    = os.path.join(WORKING_DIR, 'Binhua_FINAL_Barcode_4726.xlsx')
OBJECT         = os.path.join(STARSOLO_DIR, 'merged', 'shiv_host_final_embedded.h5ad')
OUT_DIR        = os.path.join(WORKING_DIR, 'annotation_output', 'vpr_barcode_linkage')

LIB_NAME = {
    'A': 'SHIV_Pre_PBMC',             'B': 'SHIV_Pre_LN',
    'C': 'SHIV_Pre_PBMC_IndianMixed', 'D': 'SHIV_21DPI_Pre_LN_Mixed',
    'E': 'SHIV_21DPI_PBMC',           'F': 'SHIV_21DPI_NonInf_LN_Mixed',
    'G': 'SHIV_Necropsy_PBMC',        'H': 'SHIV_Necropsy_LN',
}
LIBS     = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
INFECTED = {'D', 'E', 'F', 'H'}
NEGATIVE = {'A', 'B', 'C', 'G'}

# vpr window (0-based) and SHIV contig sanity length
WIN_START1, WIN_END1 = 5701, 6210
SHIV_LEN_EXPECTED = 10262

# barcode geometry
BC_LEN        = 34
ANCHOR_RE     = re.compile('ACGCGC.{22}GCGCGT')
INNER_ORIENT1 = 'AGCGTA'
INNER_ORIENT2 = 'TAGCTG'

# CellRanger sources
INCLUDE_CELLRANGER_SAMPLE = True
INCLUDE_CELLRANGER_UNASSIGNED = True

# obs columns — set explicitly now that the object has been inspected.
# NOTE: shiv_host_final_embedded has NO tier2_subtype. Available annotation is
# tier1_celltype (broad lineage) + t_subtype (T cells only). Set overrides here.
COL_LIBRARY   = 'library'
COL_ANIMAL    = 'animal_id'
COL_TISSUE    = 'tissue'
COL_TIMEPOINT = 'timepoint'
COL_CELLTYPE  = 'tier1_celltype'
# extra obs columns carried into the table: finer T type, subspecies, and the
# independent SHIV-detection flag from the earlier recovery.
COL_EXTRA = ['t_subtype', 'subspecies', 'condition', 'sample_name', 'shiv_pos', 'shiv_umi']

# library-derived fallbacks when a barcoded CB is not a called cell in the object.
# Tissue is unambiguous per library; timepoint is None for the mixed libraries.
LIB_TISSUE    = {'A': 'PBMC', 'B': 'LN', 'C': 'PBMC', 'D': 'LN',
                 'E': 'PBMC', 'F': 'LN', 'G': 'PBMC', 'H': 'LN'}
LIB_TIMEPOINT = {'A': 'Pre', 'B': 'Pre', 'C': 'Pre', 'D': None,
                 'E': '21 DPI', 'F': None, 'G': 'Necropsy', 'H': 'Necropsy'}

SAMTOOLS = os.environ.get('SAMTOOLS') or shutil.which('samtools')

os.makedirs(OUT_DIR, exist_ok=True)
pd.set_option('display.width', 190)
pd.set_option('display.max_columns', 60)
pd.set_option('display.max_rows', 200)

_COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
CIGAR_RE = re.compile(r'(\d+)([MIDNSHP=X])')


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


def parse_cigar(cig):
    if cig == '*' or not cig:
        return []
    return [(int(n), op) for n, op in CIGAR_RE.findall(cig)]


def lead_sc(ops):
    return ops[0][0] if ops and ops[0][1] == 'S' else 0


def trail_sc(ops):
    return ops[-1][0] if ops and ops[-1][1] == 'S' else 0


def raw_cb_from_name(name):
    m = re.search(r'[ACGT]{16}', name)
    return m.group(0) if m else None


def norm_cb(cb):
    """Reduce a CB tag to its raw 16bp so STARsolo (no suffix) and CellRanger
    (…-1 gem suffix) collapse to the same cell, and so it matches obs_names."""
    if not cb:
        return ''
    m = re.search(r'[ACGT]{16}', cb)
    return m.group(0) if m else ''


def hr(t):
    print('\n' + '=' * 84); print(t); print('=' * 84)


def save_csv(df, name, index=True):
    p = os.path.join(OUT_DIR, name)
    df.to_csv(p, index=index)
    print(f'    saved: {name}  ({len(df)} rows)')


print('vpr barcode linkage — config')
print(f'  samtools : {SAMTOOLS}')
print(f'  object   : {OBJECT}')
print(f'  out      : {OUT_DIR}')


# %% ============================================================
# CELL 2 — Binhua ground truth + data-driven flank probes (for Tier 1 + Tier 2)
# ============================================================
hr('CELL 2 — Binhua ground truth cores + cassette flank probes')

gt_core = defaultdict(set)
prefix_counts, suffix_counts = Counter(), Counter()
have_gt = False

LIB_ANIMALS = {
    'A': ['39272', '40702', '40707', '41861', '41862'],
    'B': ['39272', '40702', '40707', '41861', '41862'],
    'C': ['34315', '41903'],
    'D': ['39272', '40702', '40707', '41861', '41862'],
    'E': ['39272', '40702', '40707', '41861', '41862'],
    'F': ['34315', '41903'],
    'G': ['40702', '40707', '41861', '41862'],
    'H': ['40702', '40707', '41861', '41862'],
}

if not os.path.exists(BINHUA_XLSX):
    print(f'  WARN: Binhua not found ({BINHUA_XLSX}); extraction cannot be GT-gated.')
else:
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

gt_core_all = set().union(*gt_core.values()) if gt_core else set()
# animal lookup for a core (for cross-check; a core is animal-private in Binhua)
core_to_animals = defaultdict(set)
for a, cores in gt_core.items():
    for c in cores:
        core_to_animals[c].add(a)
print(f'  pooled unique cores: {len(gt_core_all)}')

# Build Tier-2 flank probes data-driven: valid 5' flanks start ACGCGC, valid 3'
# flanks end GCGCGT. rule = how to canonicalize the adjacent 10bp core.
FLANK5, FLANK3 = {}, {}   # seq -> 'asis' | 'rc'
for seq, n in prefix_counts.most_common():
    if seq.startswith('ACGCGC') and n >= 20:
        inner = seq[6:12]
        if inner == INNER_ORIENT1:
            FLANK5[seq] = 'asis'
        elif inner == INNER_ORIENT2:
            FLANK5[seq] = 'rc'
for seq, n in suffix_counts.most_common():
    if seq.endswith('GCGCGT') and n >= 20:
        inner_equiv = revcomp(seq[0:6])
        if inner_equiv == INNER_ORIENT1:
            FLANK3[seq] = 'asis'
        elif inner_equiv == INNER_ORIENT2:
            FLANK3[seq] = 'rc'
print(f'  Tier-2 5-prime flank probes: {FLANK5}')
print(f'  Tier-2 3-prime flank probes: {FLANK3}')


def match_gt(cand10):
    """Return the GT-canonical core if cand10 or its revcomp is in ground truth."""
    if cand10 in gt_core_all:
        return cand10
    rc = revcomp(cand10)
    if rc in gt_core_all:
        return rc
    return None


def extract_core(seq):
    """
    Return (core, tier, gt_bool). tier 1 = full anchor, tier 2 = one-sided flank.
    GT-matched cores are returned with gt=True; a valid-but-unlisted Tier-2
    candidate returns gt=False and is NOT trusted downstream.
    """
    # Tier 1: full anchor, either strand
    for s in (seq, revcomp(seq)):
        m = ANCHOR_RE.search(s)
        if m:
            c = canon_core_from_anchor(m.group(0))
            if c:
                g = match_gt(c)
                return (g or c, 1, g is not None)
    # Tier 2: one-sided constant flank, either strand
    for s in (seq, revcomp(seq)):
        for fl, rule in FLANK5.items():
            i = s.find(fl)
            if i != -1:
                raw = s[i + 12:i + 22]
                if len(raw) == 10 and set(raw) <= set('ACGT'):
                    cand = raw if rule == 'asis' else revcomp(raw)
                    g = match_gt(cand)
                    if g:
                        return (g, 2, True)
        for fl, rule in FLANK3.items():
            i = s.find(fl)
            if i != -1:
                raw = s[i - 10:i]
                if len(raw) == 10 and set(raw) <= set('ACGT'):
                    cand = raw if rule == 'asis' else revcomp(raw)
                    g = match_gt(cand)
                    if g:
                        return (g, 2, True)
    return (None, 0, False)


# %% ============================================================
# CELL 3 — OBJECT INTROSPECTION + JOIN-KEY VERIFICATION (print before trusting)
# ============================================================
hr('CELL 3 — embedded object obs_names, columns, and per-library suffix scheme')

obs = None
lookup = {}                 # (libval, raw_CB) -> metadata dict
letter_to_libval = {}
col = {}

if not os.path.exists(OBJECT):
    print(f'  WARN: object not found ({OBJECT}); cell linkage will be skipped.')
else:
    import scanpy as sc
    adata = sc.read_h5ad(OBJECT, backed='r')
    obs = adata.obs
    names = list(obs.index[:12])
    print(f'  n_obs: {adata.n_obs:,}')
    print(f'  obs.columns ({len(obs.columns)}): {list(obs.columns)}')
    print('\n  sample obs_names:')
    for n in names:
        print(f'    {n}')

    # detect the CB + suffix scheme
    cb0 = raw_cb_from_name(names[0]) if names else None
    if cb0:
        suffixes = Counter(nm.replace(raw_cb_from_name(nm) or '', '|', 1).split('|')[-1]
                           for nm in obs.index[:20000] if raw_cb_from_name(nm))
        print('\n  detected 16bp CB inside obs_names; top suffixes (post-CB):')
        for suf, n in suffixes.most_common(12):
            print(f'    suffix [{suf}]  x{n}')
    else:
        print('\n  *** could not find a 16bp ACGT run in obs_names; '
              'set the join manually before trusting linkage.')

    # column auto-detect (override in Cell 1 if any is wrong)
    def pick(overr, keys, prefer=None):
        if overr:
            return overr
        cands = [c for c in obs.columns if any(k in c.lower() for k in keys)]
        if prefer:
            for p in prefer:
                for c in cands:
                    if p in c.lower():
                        return c
        return cands[0] if cands else None

    col['library']   = pick(COL_LIBRARY,   ['library', 'sample', 'gem', 'batch', 'channel'],
                            prefer=['library', 'gem_well', 'sample'])
    col['animal']    = pick(COL_ANIMAL,    ['animal', 'monkey', 'subject'])
    col['tissue']    = pick(COL_TISSUE,    ['tissue', 'compartment', 'organ'])
    col['timepoint'] = pick(COL_TIMEPOINT, ['timepoint', 'time', 'dpi', 'visit', 'day', 'week'])
    col['celltype']  = pick(COL_CELLTYPE,  ['tier2', 'subtype', 'celltype', 'cell_type', 'annotation'],
                            prefer=['tier2_subtype', 'tier2', 'subtype'])
    print('\n  chosen obs columns (override in Cell 1 if wrong):')
    for k, v in col.items():
        print(f'    {k:<10} -> {v}')
        if v is not None:
            vc = obs[v].astype(str).value_counts().head(8)
            for val, n in vc.items():
                print(f'         {val:<32} {n}')

    extra_cols = [c for c in COL_EXTRA if c in obs.columns]
    print(f'\n  extra columns carried into table: {extra_cols}')

    # map library letter -> obs library value. Prefer EXACT match; the obs values
    # here are the letters themselves. Never substring-match a single letter
    # (e.g. "B" is inside "PBMC"), which silently mismapped C/E/F/G before.
    if col['library']:
        libvals = [str(v) for v in obs[col['library']].astype(str).unique()]
        libval_set = set(libvals)
        for L in LIBS:
            if L in libval_set:
                letter_to_libval[L] = L
            elif LIB_NAME[L] in libval_set:
                letter_to_libval[L] = LIB_NAME[L]
            else:
                hit = None
                for v in libvals:
                    if (v.endswith('_' + L) or v.endswith('-' + L)
                            or v.startswith(L + '_') or v.startswith(L + '-')):
                        hit = v
                        break
                letter_to_libval[L] = hit
        print('\n  library letter -> obs library value:')
        for L in LIBS:
            flag = '' if letter_to_libval[L] else '   *** NO MATCH — set COL_LIBRARY / check mapping'
            print(f'    {L}  {LIB_NAME[L]:<28} -> {letter_to_libval[L]}{flag}')

        # build the lookup keyed by (libval, raw_CB)
        libcol = col['library']
        for name, r in obs.iterrows():
            cb = raw_cb_from_name(name)
            if cb is None:
                continue
            meta = {
                'obs_name': name,
                'animal': r[col['animal']] if col['animal'] else None,
                'tissue': r[col['tissue']] if col['tissue'] else None,
                'timepoint': r[col['timepoint']] if col['timepoint'] else None,
                'celltype': r[col['celltype']] if col['celltype'] else None,
            }
            for c in extra_cols:
                meta[c] = r[c]
            lookup[(str(r[libcol]), cb)] = meta
        print(f'\n  lookup built: {len(lookup):,} (libval, CB) keys')
    else:
        print('  *** no library column detected; linkage will fall back to CB-only '
              '(risk of cross-library CB collision). Set COL_LIBRARY in Cell 1.')


# %% ============================================================
# CELL 4 — BAM manifest (coord-valid SHIV BAMs only)
# ============================================================
hr('CELL 4 — BAM manifest')


def detect_shiv_contig(bam):
    try:
        hd = subprocess.run([SAMTOOLS, 'view', '-H', bam], capture_output=True, text=True).stdout
    except Exception as e:
        print(f'    header read failed {os.path.basename(bam)} ({e})')
        return None, None
    best = None
    for line in hd.splitlines():
        if not line.startswith('@SQ'):
            continue
        nm = ln = None
        for tok in line.split('\t'):
            if tok.startswith('SN:'):
                nm = tok[3:]
            elif tok.startswith('LN:'):
                ln = int(tok[3:])
        if nm and (('SHIV' in nm.upper()) or ('AD8' in nm.upper())):
            return nm, ln
        if ln and 9000 <= ln <= 12000:
            best = (nm, ln)
    return best if best else (None, None)


manifest = []
for L in LIBS:
    bam = os.path.join(SOLO_DIR, f'{L}_{LIB_NAME[L]}_Aligned.sortedByCoord.out.bam')
    if os.path.exists(bam):
        nm, ln = detect_shiv_contig(bam)
        if ln == SHIV_LEN_EXPECTED:
            manifest.append({'lib': L, 'group': 'INF' if L in INFECTED else 'neg',
                             'source': 'starsolo', 'bam': bam, 'contig': nm})

for L in LIBS:
    sample = LIB_NAME[L]
    base = os.path.join(CELLRANGER_DIR, sample, 'outs')
    paths = []
    if INCLUDE_CELLRANGER_SAMPLE:
        for p in glob.glob(os.path.join(base, 'per_sample_outs', sample, '**', '*.bam'),
                           recursive=True):
            if 'vdj' not in p.lower():
                paths.append(('cr_sample', p))
    if INCLUDE_CELLRANGER_UNASSIGNED:
        u = os.path.join(base, 'unassigned_alignments.bam')
        if os.path.exists(u):
            paths.append(('cr_unassigned', u))
    for src, p in paths:
        nm, ln = detect_shiv_contig(p)
        if ln == SHIV_LEN_EXPECTED:
            manifest.append({'lib': L, 'group': 'INF' if L in INFECTED else 'neg',
                             'source': src, 'bam': p, 'contig': nm})

man_df = pd.DataFrame(manifest)
if len(man_df):
    print(man_df[['lib', 'group', 'source', 'contig']].to_string(index=False))
else:
    print('  *** no coord-valid SHIV BAMs found.')


# %% ============================================================
# CELL 5 — pull window reads and extract cores (Tier 1 + Tier 2)
# ============================================================
hr('CELL 5 — extract barcode cores from vpr-window reads')


def scan(row):
    reads = []
    region = f'{row["contig"]}:{WIN_START1}-{WIN_END1}'
    try:
        p = subprocess.Popen([SAMTOOLS, 'view', row['bam'], region],
                             stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
        for line in p.stdout:
            f = line.rstrip('\n').split('\t')
            if len(f) < 11:
                continue
            seq = f[9]
            core, tier, gt = extract_core(seq)
            if tier == 0:
                continue
            cb = ub = None
            for tag in f[11:]:
                if tag.startswith('CB:Z:'):
                    cb = tag[5:]
                elif tag.startswith('UB:Z:'):
                    ub = tag[5:]
            reads.append({'lib': row['lib'], 'group': row['group'], 'source': row['source'],
                          'core': core, 'tier': tier, 'gt': int(gt),
                          'cb': norm_cb(cb) if (cb and cb != '-') else '',
                          'ub': (ub if ub and ub != '-' else '')})
        p.wait()
    except Exception as e:
        print(f'    view failed {os.path.basename(row["bam"])} ({e})')
    return reads


barcoded = []
for _, row in (man_df.iterrows() if len(man_df) else []):
    r = scan(row)
    barcoded.extend(r)
    if r:
        t1 = sum(x['tier'] == 1 for x in r)
        t2 = sum(x['tier'] == 2 for x in r)
        gtn = sum(x['gt'] for x in r)
        print(f'  {row["lib"]} [{row["group"]}] {row["source"]:<14}: '
              f'{len(r):3d} barcoded reads (T1 {t1}, T2 {t2}), {gtn} GT, '
              f'{sum(bool(x["cb"]) for x in r)} with CB')

reads_df = pd.DataFrame(barcoded)
if len(reads_df):
    save_csv(reads_df, 'vpr_barcoded_reads.csv', index=False)
    # specificity guard
    neg_gt = reads_df[(reads_df['group'] == 'neg') & (reads_df['gt'] == 1)]
    if len(neg_gt):
        print(f'\n  *** SPECIFICITY WARNING: {len(neg_gt)} GT-matched core(s) in NEGATIVE '
              'libraries. Possible index hopping / contamination — inspect before trusting.')
    else:
        print('\n  specificity OK: zero GT-matched cores in any negative library.')
else:
    print('  No barcoded reads extracted.')


# %% ============================================================
# CELL 6 — dedup to molecules and cells; link to the object (with match rate)
# ============================================================
hr('CELL 6 — molecules, cells, and cell linkage')

cell_rows = []
if len(reads_df):
    gt_reads = reads_df[reads_df['gt'] == 1].copy()
    print(f'  GT-matched barcoded reads: {len(gt_reads)}')

    # molecules = unique (lib, cb, ub, core); reads without a UB stay distinct
    mol = gt_reads[gt_reads['cb'] != ''].copy()
    mol['ub_key'] = np.where(mol['ub'] == '', 'noUMI_' + mol.index.astype(str), mol['ub'])
    mol_unique = mol.drop_duplicates(['lib', 'cb', 'ub_key', 'core'])
    print(f'  unique barcoded molecules (dedup across sources): {len(mol_unique)}')

    # cells = unique (lib, cb); aggregate cores + umi support + sources + tier
    matched = miss = 0
    for (lib, cb), grp in mol_unique.groupby(['lib', 'cb']):
        cores = sorted(set(grp['core']))
        n_umi = grp['ub_key'].nunique()
        sources = ','.join(sorted(set(grp['source'])))
        tier = int(grp['tier'].min())  # 1 if any full-anchor read supports the cell
        gt_animals = sorted(set().union(*[core_to_animals.get(c, set()) for c in cores]))

        meta = None
        if col.get('library') and letter_to_libval.get(lib):
            meta = lookup.get((letter_to_libval[lib], cb))
        if meta:
            matched += 1
        else:
            miss += 1

        # tissue/timepoint fall back to library-level truth when the CB is not a
        # called cell in the object. animal cannot fall back here: these cores are
        # multi-animal in Binhua, so only the object join resolves the animal.
        tissue = meta['tissue'] if meta else LIB_TISSUE.get(lib, '')
        timepoint = meta['timepoint'] if meta else (LIB_TIMEPOINT.get(lib) or '')
        meta_source = 'object' if meta else 'library_fallback'

        cell_rows.append({
            'lib': lib, 'group': 'INF' if lib in INFECTED else 'neg',
            'sample': LIB_NAME[lib], 'cb': cb,
            'cores': ';'.join(cores), 'n_cores': len(cores),
            'n_umis': n_umi, 'best_tier': tier, 'sources': sources,
            'barcode_animals_binhua': ';'.join(gt_animals),
            'in_object': meta is not None, 'meta_source': meta_source,
            'obs_name': meta['obs_name'] if meta else '',
            'animal': meta['animal'] if meta else '',
            'tissue': tissue, 'timepoint': timepoint,
            'celltype': meta['celltype'] if meta else '',
            't_subtype': meta.get('t_subtype') if meta else '',
            'subspecies': meta.get('subspecies') if meta else '',
            'shiv_pos': meta.get('shiv_pos') if meta else '',
            'shiv_umi': meta.get('shiv_umi') if meta else '',
        })

    cells_df = pd.DataFrame(cell_rows)
    n_cells = len(cells_df)
    print(f'\n  unique barcode-positive cells: {n_cells}')
    if n_cells:
        rate = 100 * matched / n_cells
        print(f'  linked to embedded object: {matched}/{n_cells} ({rate:.0f}%)   unlinked: {miss}')
        if matched == 0:
            print('  *** JOIN MATCH RATE 0% — the (library, CB) key is not hitting the object.')
            print('      Check Cell 3: obs_names suffix scheme, COL_LIBRARY, and the')
            print('      letter->libval map. Do NOT read the empty linkage as "no virus cells".')
        elif rate < 60:
            print('  *** LOW match rate — verify the join before trusting the table.')

        save_csv(cells_df, 'vpr_barcode_by_cell.csv', index=False)

        # independent cross-check: was a barcode-positive cell also flagged
        # SHIV-positive by the earlier detection stored in the object?
        if 'shiv_pos' in cells_df.columns:
            linked = cells_df[cells_df['in_object']]
            if len(linked):
                sp = linked['shiv_pos'].astype(str).isin(['True', '1', '1.0', 'True '])
                print(f'\n  independent SHIV flag: {int(sp.sum())}/{len(linked)} linked '
                      'barcode+ cells are also shiv_pos in the object')

        # animal check, honestly caveated: these cores are multi-animal in Binhua,
        # so this is a plausibility check, not a unique-animal validation.
        chk = cells_df[(cells_df['in_object']) & (cells_df['barcode_animals_binhua'] != '')]
        if len(chk):
            agree = chk.apply(
                lambda r: str(r['animal']) in r['barcode_animals_binhua'].split(';'), axis=1).sum()
            print(f'  animal plausibility (object animal within Binhua candidates): '
                  f'{agree}/{len(chk)}')
            print('    (cores here are shared across animals in Binhua, so this only '
                  'checks the object animal is among the candidates, not uniqueness)')
else:
    cells_df = pd.DataFrame()


# %% ============================================================
# CELL 7 — the deliverable break-out tables
# ============================================================
hr('CELL 7 — barcode-by-cell table + source break-outs')

if len(cells_df):
    show_cols = ['lib', 'sample', 'group', 'cb', 'cores', 'n_umis', 'best_tier',
                 'sources', 'in_object', 'meta_source', 'animal', 'tissue', 'timepoint',
                 'celltype', 't_subtype', 'subspecies', 'shiv_pos', 'shiv_umi',
                 'barcode_animals_binhua']
    show_cols = [c for c in show_cols if c in cells_df.columns]
    print('\n  BARCODE-BY-CELL (the headline table):')
    print(cells_df[show_cols].to_string(index=False))

    # break-out: virus+ cells per animal x tissue x timepoint x celltype (linked only)
    linked = cells_df[cells_df['in_object']]
    if len(linked):
        brk = (linked.groupby(['animal', 'tissue', 'timepoint', 'celltype'])
               .agg(cells=('cb', 'nunique'), umis=('n_umis', 'sum'),
                    cores=('cores', lambda s: len(set(';'.join(s).split(';')))))
               .reset_index().sort_values('cells', ascending=False))
        print('\n  VIRUS+ CELLS by animal x tissue x timepoint x cell type:')
        print(brk.to_string(index=False))
        save_csv(brk, 'vpr_cells_by_source.csv', index=False)
    else:
        print('\n  No cells linked to the object; source break-out skipped '
              '(fix the join in Cell 3 first).')

    # per-library tally regardless of linkage
    tally = (cells_df.groupby(['lib', 'sample', 'group'])
             .agg(cells=('cb', 'nunique'), umis=('n_umis', 'sum'),
                  linked=('in_object', 'sum'))
             .reset_index())
    print('\n  per-library virus+ cell tally:')
    print(tally.to_string(index=False))
    save_csv(tally, 'vpr_cells_by_library.csv', index=False)
else:
    print('  No barcode-positive cells to tabulate.')


# %% ============================================================
# CELL 8 — SUMMARY
# ============================================================
hr('CELL 8 — SUMMARY')
n_reads = len(reads_df) if len(reads_df) else 0
n_gt = int(reads_df['gt'].sum()) if len(reads_df) else 0
n_mol = len(mol_unique) if len(reads_df) else 0
n_cells = len(cells_df) if len(cells_df) else 0
n_linked = int(cells_df['in_object'].sum()) if len(cells_df) else 0
print(f"""
  Barcoded reads extracted     : {n_reads}   (GT-matched: {n_gt})
  Unique barcoded molecules    : {n_mol}
  Unique barcode-positive cells: {n_cells}   (linked to object: {n_linked})

  This is the current mapped-read ceiling. The table columns are already the
  shape we want at scale, so when the unmapped-pool recovery adds reads later
  the same script re-runs and the table just fills in. The two levers to raise
  these counts, both for next session:
    1. Reference surgery at the vpr cassette locus (recover clipped reads).
    2. A3-hypermutation / pol-divergence tolerant remapping (recover SHIV reads
       that never mapped at all), upstream of the barcode locus.

  Read vpr_barcode_by_cell.csv for the per-cell detail and vpr_cells_by_source.csv
  for the animal x tissue x timepoint x cell-type break-out.
""")
print('=' * 84)
print('vpr BARCODE LINKAGE COMPLETE — read-only. No BAM/reference/object modified.')
print('=' * 84)
