#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Rhesus SHIV — VPR/VPX Barcode Recovery Feasibility (Read Level, Read-Only)
==========================================================================
Can we pull the 34bp Keele-lab viral barcode (the 10-base clonotype tag between
vpx and vpr) out of the 5' single-cell reads, link it to a cell, and track the
SHIV population across timepoints? This script answers the FEASIBILITY question
only. It counts, it does not build a barcode-by-cell matrix yet, and it writes
nothing except CSV summaries.

Barcode design (Khanal et al. 2020, J Virol; confirmed empirically against the
Binhua ground-truth file):
    34bp insert = [12bp flank][10bp random core][12bp flank]
    outer 6bp are the MluI ligation scar: ACGCGC ... GCGCGT (both orientations)
    -> anchor regex 'ACGCGC.{22}GCGCGT' is its own reverse complement at the
       outer level, so a forward-strand search catches either strand.
    -> the 10bp core is the clonotype ID; inner flank AGCGTA (orient 1) vs
       TAGCTG (orient 2) tells us whether to revcomp the core to canonicalize.

Modules
-------
  0  Ground truth + anchors : per-animal known-barcode set from the Binhua file;
                              lock the anchor / core geometry.
  1  Locate the vpx/vpr site : find the barcode anchor and the Khanal primers in
                              the SHIVAD8EO reference to bracket the insert.
  2  Three-tier anchor count : per library,
                                 T1 soft-clip hits in the STARsolo BAM (cell-linked, instant)
                                 T2 raw R2 fastq anchor hits (the true ceiling)
                                 T3 T2 hits whose paired R1 carries a valid whitelist CB
                                    (what is actually recoverable per cell)
  3  Feasibility readout      : recovered cores vs the Binhua ground truth
                                (exact-34 and canonical-10 match rates), infected
                                vs negative libraries as the specificity control.

Design
------
  - Read-only. No object touched, no barcode-by-cell matrix built here.
  - Warn-and-continue on any missing BAM / fastq / reference / whitelist.
  - Negatives (A/B/C/G) run alongside infected (D/E/F/H) as the background /
    specificity control; trim LIBS at the top if you want a faster pass.
  - Full fastq scan by default (MAX_READS_PER_LIB = None).
  - Spyder cells (# %%). Runs headless under SLURM (sc_pre env).

Author: Jake Lehle
Date: July 2026
"""

# %% ============================================================
# CELL 1 — CONFIG AND DIALS
# ============================================================
import os
import re
import gzip
import glob
import subprocess
import shutil
from collections import Counter, defaultdict
import numpy as np
import pandas as pd

WORKING_DIR  = '/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV'
STARSOLO_DIR = os.path.join(WORKING_DIR, 'shiv_starsolo')
REF_FA       = os.path.join(STARSOLO_DIR, 'ref', 'SHIVAD8EO.fa')
WHITELIST    = os.path.join(STARSOLO_DIR, 'ref', '3M-5pgex-jan-2023.txt')
BINHUA_XLSX  = os.path.join(WORKING_DIR, 'Binhua_FINAL_Barcode_4726.xlsx')
FASTQ_BASE   = '/master/zwallis/WORKING/SC/10X_RAW'
OUT_DIR      = os.path.join(WORKING_DIR, 'annotation_output', 'vpr_barcode_feasibility')

# Library -> CellRanger sample name (authoritative, from shiv_VDJ_integration.py)
LIB_NAME = {
    'A': 'SHIV_Pre_PBMC',             'B': 'SHIV_Pre_LN',
    'C': 'SHIV_Pre_PBMC_IndianMixed', 'D': 'SHIV_21DPI_Pre_LN_Mixed',
    'E': 'SHIV_21DPI_PBMC',           'F': 'SHIV_21DPI_NonInf_LN_Mixed',
    'G': 'SHIV_Necropsy_PBMC',        'H': 'SHIV_Necropsy_LN',
}
# All 8; negatives serve as the background/specificity control. Trim if desired.
LIBS = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
INFECTED = {'D', 'E', 'F', 'H'}   # Jake's classification from the STARsolo script
NEGATIVE = {'A', 'B', 'C', 'G'}

# Library -> animals contributing (from Group_H_CITEseq_Matrix_Info). Only the 5
# Chinese infected animals have ground-truth barcodes in the Binhua file.
LIB_ANIMALS = {
    'A': ['39272', '40702', '40707', '41861', '41862'],   # Pre PBMC (no virus yet)
    'B': ['39272', '40702', '40707', '41861', '41862'],   # Pre LN   (no virus yet)
    'C': ['34315', '41903'],                              # Pre Indian controls
    'D': ['39272', '40702', '40707', '41861', '41862'],   # 21 DPI LN
    'E': ['39272', '40702', '40707', '41861', '41862'],   # 21 DPI PBMC
    'F': ['34315', '41903'],                              # 21 DPI non-inf LN
    'G': ['40702', '40707', '41861', '41862'],            # Necropsy PBMC
    'H': ['40702', '40707', '41861', '41862'],            # Necropsy LN (reservoir)
}

# Barcode geometry
BC_LEN        = 34
ANCHOR_STR    = 'ACGCGC.{22}GCGCGT'
ANCHOR_RE     = re.compile(ANCHOR_STR)
ANCHOR_RE_B   = re.compile(ANCHOR_STR.encode())
PREFILTER_B   = b'ACGCGC'          # cheap substring gate before the regex
INNER_ORIENT1 = 'AGCGTA'           # core taken as-is
INNER_ORIENT2 = 'TAGCTG'           # core reverse-complemented to canonicalize

# 10x GEM-X 5' v3 read structure (from shiv_starsolo_recovery.sh)
CB_LEN  = 16
UMI_LEN = 12

# Khanal primers to bracket the vpx/vpr insert in the reference
VPXF1  = 'CTAGGGGAAGGACATGGGGCAGG'    # vpx forward
VPR_RT = 'CAGGTTGGCCGATTCTGGAGT'      # vpr cDNA (RT) primer

MAX_READS_PER_LIB = None    # None = full scan; set an int to cap for a quick test
PROGRESS_EVERY    = 25_000_000

os.makedirs(OUT_DIR, exist_ok=True)
pd.set_option('display.width', 170)
pd.set_option('display.max_columns', 50)

_COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


def revcomp(s):
    return s.translate(_COMP)[::-1]


def canon_core(bc34):
    """Given a forward-strand 34bp anchor match (ACGCGC..GCGCGT), return the
    orientation-canonical 10bp core, or None if the inner flank is corrupted."""
    if len(bc34) != BC_LEN:
        return None
    inner = bc34[6:12]
    core = bc34[12:22]
    if inner == INNER_ORIENT1:
        return core
    if inner == INNER_ORIENT2:
        return revcomp(core)
    return None


def hr(title):
    print("\n" + "=" * 80)
    print(title)
    print("=" * 80)


def save_csv(df, name, index=True):
    p = os.path.join(OUT_DIR, name)
    df.to_csv(p, index=index)
    print(f"    saved: {name}")


# %% ============================================================
# CELL 2 — MODULE 0 (ground-truth barcodes + anchor geometry)
# ============================================================
hr("CELL 2 — MODULE 0 (ground-truth barcode set + anchors)")
print(f"  anchor regex        : {ANCHOR_STR}  (self-revcomp at the outer 6bp)")
print(f"  barcode length      : {BC_LEN}  = 12 flank + 10 core + 12 flank")
print(f"  core canonicalized  : inner {INNER_ORIENT1} as-is / {INNER_ORIENT2} revcomp")

# per-animal known barcode sets: exact 34bp and canonical 10bp core
gt_exact = defaultdict(set)     # animal -> {34bp as stored, both orientations}
gt_core  = defaultdict(set)     # animal -> {canonical 10bp core}
have_gt = False

if not os.path.exists(BINHUA_XLSX):
    print(f"\n  WARN: Binhua file not found at {BINHUA_XLSX}")
    print("        Modules 0-2 still run; Module 3 cross-reference will be skipped.")
    print("        Copy Binhua_FINAL_Barcode_4726.xlsx into the working dir to enable it.")
else:
    try:
        from openpyxl import load_workbook
        wb = load_workbook(BINHUA_XLSX, read_only=True, data_only=True)
        for sh in wb.sheetnames:
            animal = sh.split()[0]
            ws = wb[sh]
            n = 0
            for row in ws.iter_rows(values_only=True):
                if row and row[0] is not None and str(row[0]).strip().isdigit():
                    bc = row[2]
                    if isinstance(bc, str):
                        bc = bc.strip().upper()
                        if len(bc) == BC_LEN and set(bc) <= set('ACGTN'):
                            # store both orientations for exact match robustness
                            gt_exact[animal].add(bc)
                            gt_exact[animal].add(revcomp(bc))
                            c = canon_core(bc) or canon_core(revcomp(bc))
                            if c:
                                gt_core[animal].add(c)
                            n += 1
            print(f"  {animal}: {n} ground-truth barcodes "
                  f"({len(gt_core[animal])} canonical cores)")
        wb.close()
        have_gt = True
    except Exception as e:
        print(f"  WARN: could not parse Binhua file ({e}); Module 3 cross-ref disabled.")
        print("        (needs openpyxl in sc_pre: conda install -n sc_pre openpyxl)")

# pooled sets for quick membership tests
gt_exact_all = set().union(*gt_exact.values()) if gt_exact else set()
gt_core_all  = set().union(*gt_core.values()) if gt_core else set()
if have_gt:
    print(f"\n  pooled ground truth: {len(gt_exact_all)//2} unique barcodes, "
          f"{len(gt_core_all)} unique canonical cores")


# %% ============================================================
# CELL 3 — MODULE 1 (locate the vpx/vpr barcode site in the reference)
# ============================================================
hr("CELL 3 — MODULE 1 (vpx/vpr insert location in SHIVAD8EO)")

if not os.path.exists(REF_FA):
    print(f"  WARN: reference not found at {REF_FA}; skipping Module 1.")
    ref_seq = ''
else:
    name, chunks = None, []
    with open(REF_FA) as fh:
        for line in fh:
            if line.startswith('>'):
                name = line[1:].strip().split()[0]
            else:
                chunks.append(line.strip().upper())
    ref_seq = ''.join(chunks)
    print(f"  contig {name}: {len(ref_seq):,} bp")

    def locate(label, probe):
        fwd = [m.start() for m in re.finditer(re.escape(probe), ref_seq)]
        rev = [m.start() for m in re.finditer(re.escape(revcomp(probe)), ref_seq)]
        print(f"  {label:<16} fwd {fwd if fwd else '-'}   rev {rev if rev else '-'}")
        return fwd, rev

    print("  probe positions (0-based):")
    anch_fwd = [m.start() for m in ANCHOR_RE.finditer(ref_seq)]
    print(f"  {'anchor 34bp':<16} at {anch_fwd if anch_fwd else '- (barcode region may be N-masked)'}")
    if anch_fwd:
        for p in anch_fwd:
            print(f"      -> reference barcode at {p}..{p+BC_LEN}: {ref_seq[p:p+BC_LEN]}")
    locate('VpxF1', VPXF1)
    locate('vpr RT primer', VPR_RT)
    print("  (barcode sits between the vpx forward primer and the vpr RT primer)")


# %% ============================================================
# CELL 4 — MODULE 2, TIER 1 (soft-clip anchor hits in STARsolo BAMs)
# ============================================================
hr("CELL 4 — MODULE 2 TIER 1 (BAM soft-clip anchor hits; cell-linked)")

samtools = shutil.which('samtools')
if samtools is None:
    print("  WARN: samtools not on PATH; Tier 1 skipped. "
          "(activate sc_pre and export CONDA_PREFIX/bin to PATH)")

tier1 = {}          # L -> dict(reads, cells, umis, cores Counter)
for L in LIBS:
    bam = os.path.join(STARSOLO_DIR, 'solo',
                       f'{L}_{LIB_NAME[L]}_Aligned.sortedByCoord.out.bam')
    rec = {'reads': 0, 'cells': set(), 'umis': set(), 'cores': Counter(),
           'exact34': Counter()}
    if samtools is None or not os.path.exists(bam):
        if samtools is not None:
            print(f"  {L}: BAM not found ({os.path.basename(bam)}); skipped")
        tier1[L] = rec
        continue
    try:
        p = subprocess.Popen([samtools, 'view', bam], stdout=subprocess.PIPE,
                             text=True, bufsize=1 << 20)
        for line in p.stdout:
            f = line.rstrip('\n').split('\t')
            if len(f) < 11:
                continue
            seq = f[9]
            m = ANCHOR_RE.search(seq)
            if not m:
                continue
            bc34 = m.group(0)
            cb = ub = None
            for tag in f[11:]:
                if tag.startswith('CB:Z:'):
                    cb = tag[5:]
                elif tag.startswith('UB:Z:'):
                    ub = tag[5:]
            rec['reads'] += 1
            rec['exact34'][bc34] += 1
            c = canon_core(bc34)
            if c:
                rec['cores'][c] += 1
            if cb and cb != '-':
                rec['cells'].add(cb)
                if ub and ub != '-':
                    rec['umis'].add((cb, ub))
        p.wait()
    except Exception as e:
        print(f"  {L}: WARN samtools/parse failed ({e})")
    tag = 'INF' if L in INFECTED else 'neg'
    print(f"  {L} [{tag}] {LIB_NAME[L]:32s}: {rec['reads']:5d} anchor reads, "
          f"{len(rec['cells']):4d} cells, {len(rec['umis']):4d} CB-UMI, "
          f"{len(rec['cores']):4d} distinct cores")
    tier1[L] = rec


# %% ============================================================
# CELL 5 — MODULE 2, TIER 2 + TIER 3 (raw fastq anchor scan + CB recovery)
# ============================================================
hr("CELL 5 — MODULE 2 TIER 2+3 (raw R2 anchor scan; R1 CB recovery)")

# whitelist (exact-match CB recovery; a floor, since 1MM would recover a few more)
whitelist = set()
if os.path.exists(WHITELIST):
    with open(WHITELIST) as fh:
        for line in fh:
            bc = line.strip()
            if bc:
                whitelist.add(bc)
    print(f"  whitelist loaded: {len(whitelist):,} barcodes  (exact-match CB, a floor)")
else:
    print(f"  WARN: whitelist not found at {WHITELIST}; Tier 3 CB check disabled.")


def open_zcat(files):
    """Stream decompressed bytes lines from a list of .fastq.gz via zcat (fast)."""
    p = subprocess.Popen(['zcat'] + files, stdout=subprocess.PIPE, bufsize=1 << 22)
    return p


def scan_library(L):
    fqdir = os.path.join(FASTQ_BASE, f'GEX_LIBRARY_{L}')
    r1s = sorted(glob.glob(os.path.join(fqdir, '*_R1_*.fastq.gz')))
    r2s = sorted(glob.glob(os.path.join(fqdir, '*_R2_*.fastq.gz')))
    out = {'total': 0, 't2': 0, 't3': 0, 'cores': Counter(), 'exact34': Counter(),
           'cells': set()}
    if not r1s or not r2s:
        print(f"  {L}: missing R1/R2 in {fqdir}; skipped")
        return out
    if len(r1s) != len(r2s):
        print(f"  {L}: WARN R1/R2 file count mismatch ({len(r1s)}/{len(r2s)})")

    p1, p2 = open_zcat(r1s), open_zcat(r2s)
    f1, f2 = p1.stdout, p2.stdout
    n = 0
    try:
        while True:
            # R1: name/seq/+/qual ; R2: name/seq/+/qual  (consume both to stay synced)
            h1 = f1.readline()
            if not h1:
                break
            s1 = f1.readline(); f1.readline(); f1.readline()
            f2.readline()
            s2 = f2.readline(); f2.readline(); f2.readline()
            n += 1
            if MAX_READS_PER_LIB and n > MAX_READS_PER_LIB:
                break
            if n % PROGRESS_EVERY == 0:
                print(f"    {L}: scanned {n:,} reads, T2={out['t2']:,} T3={out['t3']:,}")
            if PREFILTER_B not in s2:
                continue
            m = ANCHOR_RE_B.search(s2)
            if not m:
                continue
            out['t2'] += 1
            bc34 = m.group(0).decode()
            out['exact34'][bc34] += 1
            c = canon_core(bc34)
            if c:
                out['cores'][c] += 1
            if whitelist:
                cb = s1[:CB_LEN].decode(errors='ignore')
                if cb in whitelist:
                    out['t3'] += 1
                    out['cells'].add(cb)
    finally:
        out['total'] = n
        for pp in (p1, p2):
            try:
                pp.stdout.close(); pp.wait()
            except Exception:
                pass
    return out


tier23 = {}
for L in LIBS:
    tag = 'INF' if L in INFECTED else 'neg'
    print(f"\n  --- Library {L} [{tag}] {LIB_NAME[L]} ---")
    res = scan_library(L)
    tier23[L] = res
    rate = (1e6 * res['t2'] / res['total']) if res['total'] else 0
    print(f"  {L}: {res['total']:,} reads | T2 anchor {res['t2']:,} "
          f"({rate:.2f}/M) | T3 CB-valid {res['t3']:,} | "
          f"{len(res['cores']):,} distinct cores | {len(res['cells']):,} cells")


# %% ============================================================
# CELL 6 — MODULE 3 (feasibility readout + ground-truth cross-reference)
# ============================================================
hr("CELL 6 — MODULE 3 (feasibility + ground-truth cross-reference)")

rows = []
for L in LIBS:
    t1, t23 = tier1[L], tier23[L]
    # cross-reference the raw-recovered cores against this library's animals' truth
    lib_core_gt = set().union(*[gt_core.get(a, set()) for a in LIB_ANIMALS[L]]) \
        if have_gt else set()
    lib_exact_gt = set().union(*[gt_exact.get(a, set()) for a in LIB_ANIMALS[L]]) \
        if have_gt else set()

    core_hits = sum(v for c, v in t23['cores'].items() if c in lib_core_gt)
    core_tot  = sum(t23['cores'].values())
    exact_hits = sum(v for b, v in t23['exact34'].items()
                     if b in lib_exact_gt or revcomp(b) in lib_exact_gt)
    exact_tot  = sum(t23['exact34'].values())

    rows.append({
        'lib': L, 'group': 'INF' if L in INFECTED else 'neg',
        'sample': LIB_NAME[L],
        'T1_bam_reads': t1['reads'], 'T1_cells': len(t1['cells']),
        'T2_raw_anchor': t23['t2'], 'T3_cb_valid': t23['t3'],
        'T1_over_T2': round(t1['reads'] / t23['t2'], 3) if t23['t2'] else np.nan,
        'distinct_cores': len(t23['cores']),
        'core_match_gt': core_hits,
        'core_match_pct': round(100 * core_hits / core_tot, 1) if core_tot else np.nan,
        'exact34_match_gt': exact_hits,
        'exact34_match_pct': round(100 * exact_hits / exact_tot, 1) if exact_tot else np.nan,
        'gt_cores_available': len(lib_core_gt),
    })

feas = pd.DataFrame(rows).set_index('lib')
print("  per-library feasibility (T1=BAM soft-clip, T2=raw ceiling, T3=CB-recoverable):")
print(feas.to_string())
save_csv(feas, 'module3_feasibility_by_library.csv')

if have_gt:
    print("\n  READ:")
    print("   - T1/T2 near 1.0  -> most barcode reads soft-clip in the BAM already,")
    print("     so this is a clean cell-linked story off the STARsolo BAMs.")
    print("   - T1/T2 near 0    -> barcode reads mostly fully unmapped; recovery")
    print("     runs off the raw fastqs with the R1 CB rejoin (Tier 3 path).")
    print("   - core_match_pct on INFECTED high, on NEGATIVE ~0 -> anchor is")
    print("     specific and the recovered cores are real viral barcodes.")
    print("   - A/B carry infected animals but PRE-infection, so they should match")
    print("     ~0 despite a non-empty ground-truth pool (a second negative check).")

# infected-vs-negative background summary
if len(feas):
    inf = feas[feas['group'] == 'INF']
    neg = feas[feas['group'] == 'neg']
    print(f"\n  INFECTED total raw anchor (T2): {int(inf['T2_raw_anchor'].sum()):,}; "
          f"CB-valid (T3): {int(inf['T3_cb_valid'].sum()):,}")
    print(f"  NEGATIVE total raw anchor (T2): {int(neg['T2_raw_anchor'].sum()):,}; "
          f"CB-valid (T3): {int(neg['T3_cb_valid'].sum()):,}  (background floor)")

# dump the recovered core tables per library for later per-animal assignment
for L in LIBS:
    if tier23[L]['cores']:
        cdf = (pd.Series(tier23[L]['cores']).sort_values(ascending=False)
               .rename('reads').to_frame())
        cdf.index.name = 'canonical_core'
        save_csv(cdf, f'module3_cores_{L}.csv')


# %% ============================================================
# CELL 7 — SUMMARY
# ============================================================
hr("CELL 7 — SUMMARY (what to read, what it means for the VPR pipeline)")
print("""
  This was feasibility, not recovery. No barcode-by-cell matrix was built and the
  object was not touched. Read module3_feasibility_by_library.csv:

   1. Is there signal at all? INFECTED T2/T3 well above the NEGATIVE background
      floor means the barcode reads are recoverable from the 5' data.

   2. Which architecture? The T1/T2 ratio decides it:
        high -> harvest soft-clipped, already-cell-linked reads from the BAMs;
        low  -> scan raw fastqs and rejoin the R1 CB (Tier 3), heavier but works.

   3. Are the cores real? core_match_pct against the Binhua ground truth on the
      infected libraries, near-zero on negatives and on pre-infection A/B.

   4. Per-animal next step: the library pools are multi-animal, so assigning a
      recovered core to a specific animal needs the CB -> animal_id join off the
      embedded object. That is the follow-up once feasibility here is confirmed;
      module3_cores_*.csv holds the per-library core tables to feed it.
""")
print("=" * 80)
print("VPR BARCODE FEASIBILITY COMPLETE — read-only. Nothing written to the object.")
print("=" * 80)
