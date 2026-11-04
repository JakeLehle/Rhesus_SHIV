#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
00b — SHIV reference comparison and coordinate reconciliation (READ-ONLY)
=========================================================================

Purpose
-------
Settle which SHIV reference the rebuilt GEX reference should be built on, and
prove that its coordinate frame matches the one every existing barcode result
was measured in, BEFORE the GTF is written around it.

This exists because SHIV_build_BC_reference.py and the scripts downstream of it
disagree about which variant is production:

  SHIV_build_BC_reference.py   PRODUCTION_VARIANT = 'V4fullN'
  shiv_starsolo_r1_recovery.sh INDEX=idx_SHIVAD8EO_V2allN, CONTIG=..._V2allN
  shiv_starsolo_bc_recovery.sh REF_ORDER="V0_nocassette V2allN V3r1 V3r2"
  shiv_bc_core_extract.py      contig = 'SHIVAD8EO_V2allN'
  SESSION_CONTEXT.md           "V2allN is the production reference"

V4fullN is designated production in the builder and appears nowhere downstream.
Every measured barcode result in this project (13 cassette reads, 9 deduplicated
molecules, 2 to 3 barcoded cells, the 609 of 711 viral cell barcodes, the 0.131%
control background) was produced against V2allN. Cell 7 of this script settles
that empirically by reading the contig names out of the BAMs on disk rather than
trusting any of the five statements above.

What it reports
---------------
  Cell 2  every SHIV-bearing FASTA on disk: length, N content, md5
  Cell 3  insertion locus derived independently per reference, vs EXPECTED_LOCUS
  Cell 4  exact base-by-base differences between every pair of references
  Cell 5  one coordinate table: primers, Keele flanks, cassette, core, in every
          frame, so any hardcoded coordinate can be checked or translated
  Cell 6  alignment cost per reference against the real Binhua catalog: how many
          anchor bases each gives STAR, and what each catalog cassette pays
  Cell 7  which contig name each BAM on disk actually carries
  Cell 8  the collapsed GTF matching each reference, and a recommendation table

Nothing is written outside OUT_DIR. No reference is built or modified here;
building is 01_mkref_collapsed.sh, and it should not run until the table in
Cell 8 has been read.

Why the cassette goes in at all
-------------------------------
A read spanning the cassette cannot align contiguously to a reference that
lacks it, so the viral read pool is specifically depleted at the one locus the
barcode arm needs. Leaving it out means the only route to those reads is
string-scanning the raw FASTQs, which also fails on APOBEC3-edited provirus
because an exact flank probe does not survive G-to-A editing while an aligner
simply counts the mismatches. That argument is set out in full in the header of
SHIV_build_BC_reference.py and is not relitigated here.

What this script establishes, that was not previously written down
------------------------------------------------------------------
* The 44bp and 42bp Keele flanks are ALREADY present in MN816822.1 and abut
  each other there, which is how the builder derives the locus. So inserting
  the cassette is a +34 change, not +120. Only 34 bases go in.
* V2allN and V4fullN are coordinate-identical. Both are len(parental) + 34 and
  differ only at the 12 MluI half-site positions inside the cassette window.
  Choosing between them therefore costs nothing downstream: 5862, and the core
  at 5874..5883, hold in both.

Environment
-----------
  conda activate sc_pre
  export PATH="${CONDA_PREFIX}/bin:${PATH}"      # system samtools is broken

Author: Jake Lehle
Date:   October 2026
"""

# %% ============================================================
# CELL 1 — CONFIGURATION (every tunable lives here)
# ============================================================
import os
import re
import glob
import shutil
import hashlib
import subprocess
from collections import Counter

import pandas as pd

try:
    from shiv_config import WORKDIR, REFDIR, LIBS, LIB_NAME            # noqa: F401
except Exception:                                                       # noqa: BLE001
    WORKDIR = '/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV'
    REFDIR = '/master/jlehle/WORKING/SC/REF'
    LIBS = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
    LIB_NAME = {
        'A': 'Pre_PBMC', 'B': 'Pre_LN', 'C': 'Pre_PBMC_IndianMixed',
        'D': '21DPI_Pre_LN_Mixed', 'E': '21DPI_PBMC',
        'F': '21DPI_NonInf_LN_Mixed', 'G': 'Necropsy_PBMC', 'H': 'Necropsy_LN',
    }

BC_REF_DIR = os.path.join(WORKDIR, 'shiv_bc_reference')
STARSOLO = os.path.join(WORKDIR, 'shiv_starsolo')

# Every place a SHIV sequence might live. Missing ones are reported, not fatal.
REFERENCE_CANDIDATES = [
    ('parental_REF',      os.path.join(REFDIR, 'SHIVAD8EO.fna')),
    ('parental_starsolo', os.path.join(STARSOLO, 'ref', 'SHIVAD8EO.fa')),
    ('keele_asis',        os.path.join(WORKDIR, 'SHIVAD8EOM_reference.txt')),
]
# Plus every *.fa the barcode reference builder produced
REFERENCE_GLOBS = [os.path.join(BC_REF_DIR, '*.fa')]

# The SHIV contig inside the combined Cell Ranger reference, pulled out by name
COMBINED_FASTA = os.path.join(REFDIR, 'Mmul10_SHIVAD8EO_v10.fa')
COMBINED_CONTIG = 'SHIVAD8EO'

# --- ground truth geometry (from Keele via Zoey, Oct 2026) -----------------
KEELE_FLANK_5 = 'GGAGACCAGGACCTCCTCCTCCTCCCCCTCCAGGACTAGCATAA'   # 44 bp
KEELE_FLANK_3 = 'ATGGAAGAAAGACCTCCAGAAAATGAAGGACCACAAAGGGAA'     # 42 bp
VPXF1 = 'CTAGGGGAAGGACATGGGGCAGG'
VPR_RT = 'CAGGTTGGCCGATTCTGGAGT'

EXPECTED_LOCUS = 5862          # 0-based, parental MN816822.1 frame
BC_LEN, CFLANK, CORE_LEN = 34, 12, 10
assert BC_LEN == 2 * CFLANK + CORE_LEN
REGISTERS = {
    'r1': ('ACGCGCAGCGTA', 'CAGCTAGCGCGT'),
    'r2': ('ACGCGCTAGCTG', 'TACGCTGCGCGT'),
}
MLUI_5, MLUI_3 = 'ACGCGC', 'GCGCGT'

# Candidate cassette windows to score, including ones not built, so the choice
# is made against numbers rather than against the builder's stale designation.
CANDIDATE_WINDOWS = {
    'V2allN':  MLUI_5 + 'N' * 22 + MLUI_3,
    'V4fullN': 'N' * BC_LEN,
    'V3r1':    REGISTERS['r1'][0] + 'N' * CORE_LEN + REGISTERS['r1'][1],
    'V3r2':    REGISTERS['r2'][0] + 'N' * CORE_LEN + REGISTERS['r2'][1],
}

# --- barcode catalog --------------------------------------------------------
BINHUA_CANDIDATES = [
    os.path.join(WORKDIR, 'Binhua_FINAL_Barcode_4-7-26_Rhesus_SHIV_Barcodes.xlsx'),
    os.path.join(WORKDIR, 'Binhua_FINAL_Barcode_4726.xlsx'),
]
BC_SEQ_COL = 2        # 0-based column index of bc_seq in every animal sheet

# Which @SQ names count as the viral contig when fingerprinting BAMs.
# A substring test on "SHIV" is wrong: cellranger multi names every VDJ
# clonotype consensus contig after the sample id, so SHIV_Necropsy_LN_
# clonotype4119_consensus_1 and ~144,000 siblings all match it. Anchor on the
# reference contig stem instead.
VIRAL_CONTIG_RE = re.compile(r'^SHIVAD8EO', re.IGNORECASE)

# --- BAMs to fingerprint for contig name -----------------------------------
BAM_GLOBS = {
    'cellranger': os.path.join(WORKDIR, 'Analysis', '*', 'outs', '**', '*.bam'),
    'starsolo_r1': os.path.join(WORKDIR, 'shiv_starsolo_r1', 'solo', '*', '*.bam'),
    'starsolo_r2': os.path.join(WORKDIR, 'shiv_starsolo_bc', 'solo', '*', '*.bam'),
    'starsolo_gex': os.path.join(WORKDIR, 'shiv_starsolo', 'solo', '*', '*.bam'),
}

OUT_DIR = os.path.join(WORKDIR, 'pipeline_results', '00b_reference_comparison')
os.makedirs(OUT_DIR, exist_ok=True)
pd.set_option('display.width', 220)
pd.set_option('display.max_columns', 80)

_COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


def revcomp(s):
    return s.translate(_COMP)[::-1]


def hr(t):
    print('\n' + '=' * 94)
    print(t)
    print('=' * 94)


def save_csv(df, name, index=False):
    p = os.path.join(OUT_DIR, name)
    df.to_csv(p, index=index, lineterminator='\n')
    print(f'  wrote {name}  ({df.shape[0]} rows)')


def read_fasta_single(path):
    """Return (name, sequence).

    Tolerates Keele's file, which carries the header and the sequence on ONE
    line with no newline between them. A naive parser silently eats the first
    few bases of his sequence as part of the header, which would shift every
    coordinate derived from it. Lifted unchanged from
    SHIV_build_BC_reference.py so the two cannot drift apart.
    """
    raw = open(path).read().strip()
    if not raw.startswith('>'):
        raise ValueError(f'{path}: not FASTA (no leading ">")')
    body = raw[1:]
    if '\n' in body:
        head, rest = body.split('\n', 1)
        cleaned = ''.join(rest.split()).upper()
        if cleaned and re.fullmatch(r'[ACGTN]+', cleaned):
            return head.strip(), cleaned
    cut = 0
    for i, ch in enumerate(body[:100]):
        if ch.upper() not in 'ACGTN':
            cut = i + 1
    name = body[:cut].strip()
    seq = ''.join(body[cut:].split()).upper()
    if not re.fullmatch(r'[ACGTN]+', seq):
        raise ValueError(f'{path}: could not separate header from sequence')
    return name, seq


def read_fasta_contig(path, contig):
    """Pull one named contig out of a multi-contig FASTA."""
    seq, grabbing = [], False
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                if grabbing:
                    break
                grabbing = line[1:].split()[0] == contig
                continue
            if grabbing:
                seq.append(line.strip())
    return ''.join(seq).upper() if seq else None


SAMTOOLS = shutil.which('samtools')
if SAMTOOLS and SAMTOOLS.startswith('/usr/local/bin'):
    print('WARN  samtools resolves to /usr/local/bin, the broken system build.\n'
          '      conda activate sc_pre and '
          'export PATH="${CONDA_PREFIX}/bin:${PATH}"')


# %% ============================================================
# CELL 2 — INVENTORY EVERY SHIV SEQUENCE ON DISK
# ============================================================
hr('CELL 2 — inventory')

refs = {}      # label -> sequence
rows = []

paths = list(REFERENCE_CANDIDATES)
for g in REFERENCE_GLOBS:
    for p in sorted(glob.glob(g)):
        paths.append((os.path.splitext(os.path.basename(p))[0], p))

for label, path in paths:
    if not os.path.exists(path):
        print(f'  MISSING  {label:<22} {path}')
        rows.append({'label': label, 'path': path, 'found': False})
        continue
    try:
        name, seq = read_fasta_single(path)
    except Exception as e:                                          # noqa: BLE001
        print(f'  WARN     {label:<22} unreadable: {type(e).__name__}: {e}')
        rows.append({'label': label, 'path': path, 'found': False})
        continue
    refs[label] = seq
    nruns = [(m.start(), m.end()) for m in re.finditer(r'N+', seq)]
    rows.append({
        'label': label, 'path': path, 'found': True, 'contig_name': name,
        'length': len(seq), 'n_bases': seq.count('N'),
        'n_runs': len(nruns),
        'longest_n_run': max((e - s for s, e in nruns), default=0),
        'n_run_positions': ';'.join(f'{s}..{e}' for s, e in nruns[:5]),
        'md5': hashlib.md5(seq.encode()).hexdigest(),
    })
    print(f'  FOUND    {label:<22} {len(seq):>7,} bp  '
          f'N={seq.count("N"):>3}  {name}')

# the SHIV contig living inside the combined Cell Ranger reference
if os.path.exists(COMBINED_FASTA):
    s = read_fasta_contig(COMBINED_FASTA, COMBINED_CONTIG)
    if s:
        refs['combined_gex_contig'] = s
        rows.append({
            'label': 'combined_gex_contig', 'path': COMBINED_FASTA,
            'found': True, 'contig_name': COMBINED_CONTIG, 'length': len(s),
            'n_bases': s.count('N'), 'n_runs': len(re.findall(r'N+', s)),
            'longest_n_run': max((len(m) for m in re.findall(r'N+', s)),
                                 default=0),
            'n_run_positions': '', 'md5': hashlib.md5(s.encode()).hexdigest(),
        })
        print(f'  FOUND    {"combined_gex_contig":<22} {len(s):>7,} bp  '
              f'N={s.count("N"):>3}  (inside {os.path.basename(COMBINED_FASTA)})')
    else:
        print(f'  WARN     contig "{COMBINED_CONTIG}" not found in '
              f'{COMBINED_FASTA}')
else:
    print(f'  MISSING  combined reference FASTA  {COMBINED_FASTA}')

inv = pd.DataFrame(rows)
save_csv(inv, 'reference_inventory.csv')

if not refs:
    raise SystemExit('No SHIV references found. Check the paths in Cell 1.')

# Identical sequences under different names are worth knowing about up front
by_md5 = Counter(r['md5'] for r in rows if r.get('found') and r.get('md5'))
dups = [m for m, c in by_md5.items() if c > 1]
if dups:
    print('\n  identical sequences under different labels:')
    for m in dups:
        same = [r['label'] for r in rows if r.get('md5') == m]
        print(f'    {", ".join(same)}')


# %% ============================================================
# CELL 3 — DERIVE THE INSERTION LOCUS INDEPENDENTLY IN EACH REFERENCE
# ============================================================
hr('CELL 3 — locus derivation per reference')

loc_rows = []
for label, seq in refs.items():
    i5 = seq.find(KEELE_FLANK_5)
    i3 = seq.find(KEELE_FLANK_3)
    rec = {'label': label, 'flank5_at': i5, 'flank3_at': i3}
    if i5 < 0 or i3 < 0:
        rec['verdict'] = 'flank not found'
        loc_rows.append(rec)
        print(f'  {label:<22} flank missing (5\' {i5}, 3\' {i3})')
        continue
    p5_end = i5 + len(KEELE_FLANK_5)
    gap = i3 - p5_end
    rec['gap_between_flanks'] = gap
    rec['locus'] = p5_end
    if gap == 0:
        rec['state'] = 'unbarcoded (flanks abut)'
        rec['verdict'] = ('locus OK' if p5_end == EXPECTED_LOCUS
                          else f'locus {p5_end} != {EXPECTED_LOCUS}')
    elif gap == BC_LEN:
        rec['state'] = 'barcoded (34bp cassette present)'
        rec['cassette_window'] = seq[p5_end:i3]
        rec['verdict'] = ('locus OK' if p5_end == EXPECTED_LOCUS
                          else f'locus {p5_end} != {EXPECTED_LOCUS}')
    else:
        rec['state'] = f'unexpected gap of {gap} bp'
        rec['verdict'] = 'INVESTIGATE'
    # primer anchors, for the independent second route
    rec['vpxf1_at'] = seq.find(VPXF1)
    rec['vprrt_rc_at'] = seq.find(revcomp(VPR_RT))
    loc_rows.append(rec)
    print(f'  {label:<22} 5\'flank {i5:>6}  gap {gap:>3}  locus {p5_end:>6}  '
          f'{rec["state"]}  [{rec["verdict"]}]')

locs = pd.DataFrame(loc_rows)
save_csv(locs, 'locus_derivation.csv')

print('\n  The gap column is the whole argument for a +34 insertion rather')
print('  than +120: in the unbarcoded reference the two Keele flanks ABUT,')
print('  so both flanks are already present in MN816822.1 and only the 34')
print('  cassette bases are missing.')


# %% ============================================================
# CELL 4 — EXACT PAIRWISE DIFFERENCES
# ============================================================
hr('CELL 4 — base-by-base differences between references')

labels = list(refs)
diff_rows = []
for i, a in enumerate(labels):
    for b in labels[i + 1:]:
        sa, sb = refs[a], refs[b]
        rec = {'ref_a': a, 'ref_b': b, 'len_a': len(sa), 'len_b': len(sb),
               'len_diff': len(sa) - len(sb)}
        if len(sa) == len(sb):
            d = [k for k, (x, y) in enumerate(zip(sa, sb)) if x != y]
            rec['n_differing_positions'] = len(d)
            rec['first_positions'] = ';'.join(str(x) for x in d[:12])
            rec['identical'] = len(d) == 0
            if d:
                rec['example'] = (f'{d[0]}: {a}={sa[d[0]]} {b}={sb[d[0]]}')
        else:
            rec['n_differing_positions'] = None
            rec['identical'] = False
            rec['first_positions'] = '(different lengths, not compared)'
        diff_rows.append(rec)

diffs = pd.DataFrame(diff_rows)
save_csv(diffs, 'pairwise_differences.csv')

same_len = diffs[diffs['len_diff'] == 0]
if not same_len.empty:
    print('  equal-length pairs:')
    print(same_len[['ref_a', 'ref_b', 'n_differing_positions',
                    'first_positions']].to_string(index=False))
    print('\n  A pair differing at exactly 12 positions, all inside the')
    print('  cassette window, is the V2allN vs V4fullN comparison: the MluI')
    print('  half-sites. Those two are coordinate-identical, so every')
    print('  hardcoded viral coordinate in the barcode arm holds in both.')
print('\n  unequal-length pairs:')
print(diffs[diffs['len_diff'] != 0][['ref_a', 'ref_b', 'len_a', 'len_b',
                                     'len_diff']].to_string(index=False))


# %% ============================================================
# CELL 5 — ONE COORDINATE TABLE ACROSS EVERY FRAME
# ============================================================
hr('CELL 5 — coordinate map (use this to check any hardcoded coordinate)')

feat_rows = []
for label, seq in refs.items():
    i5 = seq.find(KEELE_FLANK_5)
    i3 = seq.find(KEELE_FLANK_3)
    locus = i5 + len(KEELE_FLANK_5) if i5 >= 0 else None
    barcoded = (i5 >= 0 and i3 >= 0 and i3 - locus == BC_LEN)
    feat_rows.append({
        'label': label,
        'length': len(seq),
        'VpxF1_fwd': seq.find(VPXF1),
        'vprRT_revcomp': seq.find(revcomp(VPR_RT)),
        'keele_flank5_start': i5,
        'cassette_start': locus if barcoded else None,
        'cassette_end': i3 if barcoded else None,
        'core_start': locus + CFLANK if barcoded else None,
        'core_end': locus + CFLANK + CORE_LEN if barcoded else None,
        'keele_flank3_start': i3,
        'insertion_point_if_unbarcoded': locus if not barcoded else None,
    })

coords = pd.DataFrame(feat_rows)
save_csv(coords, 'coordinate_map.csv')
print(coords.to_string(index=False))

print(f'\n  shiv_bc_core_extract.py hardcodes LOCUS={EXPECTED_LOCUS}, '
      f'CFLANK={CFLANK}, CORE_LEN={CORE_LEN}, so it expects')
print(f'  the core at bases {EXPECTED_LOCUS + CFLANK}..'
      f'{EXPECTED_LOCUS + CFLANK + CORE_LEN - 1} inclusive, which is the '
      f'half-open window')
print(f'  {EXPECTED_LOCUS + CFLANK}..{EXPECTED_LOCUS + CFLANK + CORE_LEN} '
      f'reported in core_start/core_end above. Any barcoded')
print('  reference whose core_start is not that number will silently extract')
print('  the wrong 10 bases, and the cores will still look like plausible')
print('  barcodes, so this is a failure that does not announce itself.')


# %% ============================================================
# CELL 6 — ALIGNMENT COST AGAINST THE REAL BINHUA CATALOG
# ============================================================
hr('CELL 6 — what each reference costs a genuine catalog barcode')

cat_path = next((p for p in BINHUA_CANDIDATES if os.path.exists(p)), None)
catalog = []
if cat_path is None:
    print('  WARN  Binhua catalog not found at either path; Cell 6 is skipped.')
    for p in BINHUA_CANDIDATES:
        print(f'          {p}')
else:
    try:
        from openpyxl import load_workbook
        wb = load_workbook(cat_path, read_only=True, data_only=True)
        for sh in wb.sheetnames:
            for row in wb[sh].iter_rows(values_only=True):
                v = row[BC_SEQ_COL] if len(row) > BC_SEQ_COL else None
                if not isinstance(v, str):
                    continue
                s = v.strip().upper()
                if re.fullmatch(r'[ACGT]+', s) and 20 <= len(s) <= 40:
                    catalog.append((sh.split()[0], s))
        print(f'  catalog: {os.path.basename(cat_path)}')
        print(f'  entries parsed                : {len(catalog):,}')
    except Exception as e:                                          # noqa: BLE001
        print(f'  WARN  workbook parse failed: {type(e).__name__}: {e}')

if catalog:
    seqs = [s for _, s in catalog]
    full = [s for s in seqs if len(s) == BC_LEN]
    lens = Counter(len(s) for s in seqs)
    print(f'  unique sequences              : {len(set(seqs)):,}')
    print(f'  length distribution           : {dict(sorted(lens.items()))}')
    print(f'  unique {BC_LEN}-mers                : {len(set(full)):,}')

    reg_counts = Counter()
    for s in full:
        hit = [k for k, (l, r) in REGISTERS.items()
               if s.startswith(l) and s.endswith(r)]
        reg_counts[hit[0] if len(hit) == 1 else 'neither'] += 1
    print(f'  {BC_LEN}-mer rows by register       : {dict(reg_counts)}')
    print(f'  rows starting {MLUI_5}          : '
          f'{sum(s.startswith(MLUI_5) for s in full):,} / {len(full):,}')
    print(f'  rows ending   {MLUI_3}          : '
          f'{sum(s.endswith(MLUI_3) for s in full):,} / {len(full):,}')
    print('  Entries that are not 34bp are deletion variants or spreadsheet')
    print('  off-by-one, not variable-length barcodes.')

    # Cost model: STAR does not penalise reference N (measured 2026-10-27 by
    # shiv_check_bc_reference_alignment.py), so a position costs a mismatch
    # only where the reference carries a real base that the read disagrees with.
    cost_rows = []
    for nm, win in CANDIDATE_WINDOWS.items():
        assert len(win) == BC_LEN, f'{nm}: window is {len(win)}bp'
        anchor = sum(c != 'N' for c in win)
        costs = [sum(1 for a, b in zip(win, s) if a != 'N' and a != b)
                 for s in full]
        cost_rows.append({
            'variant': nm,
            'cassette_window': win,
            'anchor_bases_for_STAR': anchor,
            'mean_mismatch': round(sum(costs) / len(costs), 3),
            'max_mismatch': max(costs),
            'pct_catalog_at_zero': round(100 * sum(c == 0 for c in costs)
                                         / len(costs), 1),
        })
    cost = pd.DataFrame(cost_rows)
    save_csv(cost, 'alignment_cost_by_variant.csv')
    print('\n' + cost[['variant', 'anchor_bases_for_STAR', 'mean_mismatch',
                       'max_mismatch', 'pct_catalog_at_zero']]
          .to_string(index=False))
    print('\n  Reading this: anchor_bases_for_STAR is real sequence the aligner')
    print('  can seed on inside the cassette window. pct_catalog_at_zero is the')
    print('  share of genuine catalog barcodes that pay nothing to align there.')
    print('  A variant that is high on both is strictly better than one that is')
    print('  high on only one. The register-specific pair scores badly because')
    print('  it charges every opposite-register read for the privilege.')


# %% ============================================================
# CELL 7 — WHICH REFERENCE DID THE EXISTING BAMs ACTUALLY USE
# ============================================================
hr('CELL 7 — contig names in the BAMs on disk (settles the V2 vs V4 question)')

bam_rows = []
if not SAMTOOLS:
    print('  samtools not on PATH; Cell 7 skipped.')
else:
    for src, pat in BAM_GLOBS.items():
        hits = sorted(glob.glob(pat, recursive=True))
        if not hits:
            print(f'  {src:<14} no BAMs found')
            continue
        seen = Counter()
        for bam in hits:
            try:
                p = subprocess.run([SAMTOOLS, 'view', '-H', bam],
                                   capture_output=True, text=True, check=True)
            except Exception as e:                                  # noqa: BLE001
                print(f'    WARN {os.path.basename(bam)}: {e}')
                continue
            for line in p.stdout.splitlines():
                if not line.startswith('@SQ'):
                    continue
                m = re.search(r'SN:(\S+)', line)
                if m and VIRAL_CONTIG_RE.match(m.group(1)):
                    ln = re.search(r'LN:(\d+)', line)
                    seen[(m.group(1), int(ln.group(1)) if ln else -1)] += 1
        for (contig, length), n in seen.items():
            bam_rows.append({'source': src, 'contig': contig,
                             'contig_length': length, 'n_bams': n})
            print(f'  {src:<14} {contig:<28} {length:>7,} bp   '
                  f'in {n} BAM(s)')

    if bam_rows:
        bams = pd.DataFrame(bam_rows)
        save_csv(bams, 'bam_contig_fingerprint.csv')
        print('\n  This is the empirical answer. Whatever contig name and')
        print('  length appear here are what every existing barcode result was')
        print('  measured against, regardless of what PRODUCTION_VARIANT says')
        print('  in SHIV_build_BC_reference.py.')
    else:
        print('\n  No SHIV contigs found in any BAM header.')


# %% ============================================================
# CELL 8 — COLLAPSED GTF PER REFERENCE, AND THE RECOMMENDATION
# ============================================================
hr('CELL 8 — matching collapsed GTF and recommendation')


def collapsed_gtf_lines(contig, length):
    """ONE feature spanning the whole contig, so every sense viral read counts
    regardless of which overlapping ORF it hits. Identical in design to
    write_collapsed_gtf() in SHIV_build_BC_reference.py and to the collapsed
    GTF in shiv_starsolo_recovery.sh."""
    attrs = (f'gene_id "{contig}"; transcript_id "{contig}_tx"; '
             f'gene_name "{contig}"; gene_biotype "protein_coding";')
    return [f'{contig}\tcustom\t{f}\t1\t{length}\t.\t+\t.\t{attrs}'
            for f in ('gene', 'transcript', 'exon')]


rec_rows = []
for label, seq in refs.items():
    i5 = seq.find(KEELE_FLANK_5)
    i3 = seq.find(KEELE_FLANK_3)
    locus = i5 + len(KEELE_FLANK_5) if i5 >= 0 else None
    barcoded = (i5 >= 0 and i3 >= 0 and i3 - locus == BC_LEN)
    window = seq[locus:i3] if barcoded else ''
    rec_rows.append({
        'label': label,
        'length': len(seq),
        'cassette_present': barcoded,
        'cassette_window': window,
        'anchor_bases_in_window': sum(c != 'N' for c in window),
        # True only when a cassette is actually present AND it sits where the
        # extractor's hardcoded LOCUS expects it. An unbarcoded reference whose
        # insertion POINT is 5862 is not the same thing and must not score True.
        'core_at_expected_offset': bool(barcoded and locus == EXPECTED_LOCUS),
        'full_length_LTRs': len(seq) >= 10000,
        'collapsed_gtf_span': f'1..{len(seq)}',
    })

rec = pd.DataFrame(rec_rows)
save_csv(rec, 'recommendation_table.csv')
print(rec.to_string(index=False))

print('\n  A reference is suitable as the rebuilt GEX reference when all of:')
print('    cassette_present          True   reads can align through the locus')
print('    core_at_expected_offset   True   the barcode arm keeps working')
print('                                     without re-deriving coordinates')
print('    full_length_LTRs          True   Keele\'s file is LTR-trimmed by 598')
print('                                     and 135 bp, and unannotated LTRs')
print('                                     were already a loss channel here')
print('  Among references meeting all three, prefer the higher')
print('  anchor_bases_in_window, provided Cell 6 shows it still costs the')
print('  catalog nothing.')

print('\n  Collapsed GTF for each candidate, to be written by '
      '01_mkref_collapsed.sh:')
for _, r in rec.iterrows():
    if not r['cassette_present']:
        continue
    print(f'\n  --- {r["label"]} ---')
    for line in collapsed_gtf_lines(r['label'], r['length']):
        print(f'    {line}')

print('\n  Note the GEX contig name must match whatever 01_mkref writes into')
print('  the combined FASTA, and the host Mmul10 GTF is concatenated after')
print('  this, unchanged.')

hr('DONE — 00b_reference_comparison.py')
print(f'  outputs: {OUT_DIR}')
print('  nothing outside that directory was written, and no reference was')
print('  built or modified.')
print('\n  Read recommendation_table.csv and bam_contig_fingerprint.csv '
      'together\n  before running 01_mkref_collapsed.sh.')
