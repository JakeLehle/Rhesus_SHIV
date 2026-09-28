#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Rhesus SHIV — vpr Barcode Ground-Truth Search (Read-Only, Sensitivity Extension)
================================================================================
Successor to shiv_vpr_barcode_demux.py. Same read pool, same negative floor, same
linkage machinery. What is new is the matching layer.

WHY THIS EXISTS (read before running)
-------------------------------------
The prior demux anchored on the invariant cassette flank (ACGCGC.{22}GCGCGT) and
returned 3 barcoded reads across ~964M reads per library set. Zoey's action item
was to instead search for the literal barcode sequences from the Binhua workbook.

Searching the literal 34-mers CANNOT beat the flank anchor. The anchor fires on any
read carrying the cassette regardless of what the core says; a literal 34-mer only
fires when the core is present, unedited and error-free. Every read a 34-mer search
finds, the anchor already found. It is a strict subset. That tier (TIER 3) is
implemented here anyway, because it is the committed deliverable and it costs
nothing to carry, and because its agreement with TIER 1 is itself a useful check.

The tiers that can actually add reads are the new ones:

  TIER 4  PARTIAL  — the anchor needs all 34bp inside one read. A read whose end
                     lands mid-cassette carries real barcode information and the
                     prior code discarded it. Matches GT barcodes that overlap a
                     read end by >= MIN_PARTIAL_OVERLAP bp, requiring the overlap
                     to cover >= MIN_CORE_BASES of the 10bp core.
  TIER 5  MISMATCH — TIER 1-3 all require exact sequence. A3-edited or
                     mis-sequenced cassettes fail. Allows <= MAX_MISMATCH
                     mismatches against the KNOWN 1,245-member GT set, which is
                     safe precisely because the search space is closed.

Honest expectation: single-digit additional reads. The value is that this is a
strictly MORE sensitive search than the one that returned 3, so if it also returns
~3 you have two independent methods, one deliberately permissive, converging on the
same ceiling. That is the justification for VpxF1 amplicon enrichment.

PARSER FIX (this is the real change to the ground truth)
--------------------------------------------------------
The prior canon_core_from_anchor() hard-indexed bc34[6:12] and bc34[12:22], so any
off-register or off-length workbook entry returned None and vanished with NO
warning. Against Binhua_FINAL_Barcode_4-7-26 that silently dropped 21 of 1,264
entries. Cell 2 here rescues the five recoverable register slips, quarantines the
16 that are not safe to guess at, and PRINTS the full disposition of every row.

The quarantined set is ADE08M.991 / .4093 / .15266 / .1, which share a
CTCCAGGACTAGCATAA block, read in a register ending ACGCGT (the true MluI site
rather than the GCGCGT half-site), and recur across all or most animals. A barcode
present in five independently inoculated animals is not behaving like a clonotype.
Resolve with Binhua/Zoey before promoting these; set INCLUDE_QUARANTINED = True to
fold them in once resolved.

Read-only. Nothing on disk is modified. Spyder cells (# %%); runs under SLURM.

Author: Jake Lehle
Date: September 2026
"""

# %% ============================================================
# CELL 1 — CONFIG AND DIALS
# ============================================================
import os
import re
import glob
import gzip
import shutil
import subprocess
from collections import Counter, defaultdict

import numpy as np
import pandas as pd

import matplotlib
HEADLESS = True                      # False in Spyder for interactive plt.show()
if HEADLESS:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt

# --- paths -----------------------------------------------------------------
WORKING_DIR  = '/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV'
STARSOLO_DIR = os.path.join(WORKING_DIR, 'shiv_starsolo')
RAW_BASE     = '/master/zwallis/WORKING/SC/10X_RAW'
WHITELIST    = os.path.join(STARSOLO_DIR, 'ref', '3M-5pgex-jan-2023.txt')
OBJECT       = os.path.join(STARSOLO_DIR, 'merged', 'shiv_host_final_embedded.h5ad')

# The workbook Zoey sent 2026-09. Sequence-identical to Binhua_FINAL_Barcode_4726.xlsx
# but this is the canonical filename going forward. Five sheets, one per animal.
BINHUA_XLSX  = os.path.join(WORKING_DIR,
                            'Binhua_FINAL_Barcode_4-7-26_Rhesus_SHIV_Barcodes.xlsx')

OUT_DIR      = os.path.join(WORKING_DIR, 'annotation_output', 'vpr_barcode_gt_search')

# --- libraries -------------------------------------------------------------
LIB_NAME = {
    'A': 'SHIV_Pre_PBMC',             'B': 'SHIV_Pre_LN',
    'C': 'SHIV_Pre_PBMC_IndianMixed', 'D': 'SHIV_21DPI_Pre_LN_Mixed',
    'E': 'SHIV_21DPI_PBMC',           'F': 'SHIV_21DPI_NonInf_LN_Mixed',
    'G': 'SHIV_Necropsy_PBMC',        'H': 'SHIV_Necropsy_LN',
}
LIBS     = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']   # negatives ARE the floor; keep all 8
INFECTED = {'D', 'E', 'F', 'H'}

# --- SMOKE TEST FIRST ------------------------------------------------------
# Run once with MAX_READS_PER_LIB = 10_000_000 and LIBS = ['E','H'] to confirm the
# parse, the index build and the CB join all behave. Then set to None for the full
# pass. The full pass is the expensive one; do not discover a typo at hour six.
MAX_READS_PER_LIB = None             # None = full scan; int = cap per library
PROGRESS_EVERY    = 25_000_000

# --- matching tiers (turn on progressively; see runtime notes in Cell 6) ----
ENABLE_TIER3_EXACT_GT = True         # literal full 34-mer match (Zoey's ask)
ENABLE_TIER4_PARTIAL  = True         # GT barcode overlapping a read end
ENABLE_TIER5_MISMATCH = True         # GT barcode with <= MAX_MISMATCH mismatches

MIN_PARTIAL_OVERLAP = 16             # bp of cassette that must be visible for TIER 4
MIN_CORE_BASES      = 6              # of the 10bp core, how many must be covered
MAX_MISMATCH        = 2              # TIER 5 budget across the full 34bp

# --- ground-truth handling -------------------------------------------------
INCLUDE_QUARANTINED = False          # flip ONLY after Binhua/Zoey resolve the
                                     # ADE08M.991/.4093/.15266/.1 family

# Workbook 3' flank verification (Cell 2). These gate which malformed rows get
# rescued into the ground truth. Loosening them admits wrong cores; tightening
# them discards the recoverable AACGCG* family. Verified settings below.
FLANK3_MIN_COMPARE   = 6             # bases of 3' flank that must be present
FLANK3_MAX_MISMATCH  = 2             # mismatches allowed across those bases

# --- R1 geometry -----------------------------------------------------------
# 10x 5' GEM-X: R1 = 16bp CB + 12bp UMI (+ TSO). If R1 is longer than CB+UMI+TSO
# there is unscanned cDNA there. The prior demux never looked. Cell 5 probes the
# actual length and this dial decides what to do about it.
CB_LEN, UMI_LEN = 16, 12
R1_TSO_LEN      = 13                 # standard 5' TSO length past the UMI
SCAN_R1_TAIL    = True               # scan R1 past CB+UMI+TSO when long enough
R1_MIN_TAIL     = 20                 # only bother if at least this many bp remain

# --- cassette geometry -----------------------------------------------------
BC_LEN        = 34
CORE_SLICE    = (12, 22)             # core position within a well-formed 34-mer
ANCHOR_RE     = re.compile('ACGCGC.{22}GCGCGT')
INNER_ORIENT1 = 'AGCGTA'
INNER_ORIENT2 = 'TAGCTG'

# Byte prefilter. Any read containing one of these goes to the expensive stage.
# Widening this set increases sensitivity AND runtime roughly linearly; the two
# scar half-sites alone pass ~4% of reads, the full set ~12-15%.
PREFILTER_SEEDS = [b'ACGCGC', b'GCGCGT', b'AGCGTA', b'TAGCTG', b'CAGCTA', b'TACGCT']

# --- object obs columns (verified against shiv_host_final_embedded) --------
COL_LIBRARY, COL_ANIMAL, COL_TISSUE = 'library', 'animal_id', 'tissue'
COL_TIMEPOINT, COL_CELLTYPE, COL_TSUB = 'timepoint', 'tier1_celltype', 't_subtype'
COL_SHIVPOS, COL_SHIVUMI, COL_SUBSP = 'shiv_pos', 'shiv_umi', 'subspecies'

NON_TARGET = {'Erythrocytes', 'Platelets', 'Unassigned', 'Mast cells', '', 'nan'}

# --- figures ---------------------------------------------------------------
FIG_DPI   = 300
FONT_BASE = 28
HEX = {'t1': '#2C6E9B', 't2': '#5BA4CF', 't3': '#8FBF4D',
       't4': '#E8A33D', 't5': '#C4453C', 'neg': '#8A8A8A'}

os.makedirs(OUT_DIR, exist_ok=True)
pd.set_option('display.width', 220)
pd.set_option('display.max_columns', 60)
pd.set_option('display.max_rows', 400)

_COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


def revcomp(s):
    return s.translate(_COMP)[::-1]


def hr(t):
    print('\n' + '=' * 90)
    print(t)
    print('=' * 90)


def save_csv(df, name, index=False):
    p = os.path.join(OUT_DIR, name)
    df.to_csv(p, index=index)
    print(f'    saved: {name}  ({len(df)} rows)')


def save_fig(fig, stem):
    for ext in ('pdf', 'png'):
        p = os.path.join(OUT_DIR, f'{stem}.{ext}')
        fig.savefig(p, dpi=FIG_DPI, bbox_inches='tight')
    print(f'    saved: {stem}.pdf / .png')


# %% ============================================================
# CELL 2 — ROBUST WORKBOOK PARSE (the parser fix lives here)
# ============================================================
hr('CELL 2 — Binhua workbook parse with register rescue')

# Quarantine list: off-register family that recurs across animals. Not safe to
# promote without Binhua/Zoey. See module docstring.
QUARANTINE_NAMES = {'ADE08M.991', 'ADE08M.4093', 'ADE08M.15266', 'ADE08M.1'}


# Canonical cassette form. The flanks are synthetic and invariant, so a barcode is
# fully specified by its 10bp core and we can always rebuild the ideal 34-mer:
CANON_5 = 'ACGCGC' + INNER_ORIENT1          # ACGCGCAGCGTA
CANON_3 = 'CAGCTA' + 'GCGCGT'               # CAGCTAGCGCGT
ORIENT2_5 = 'ACGCGC' + INNER_ORIENT2        # ACGCGCTAGCTG
ORIENT2_3 = 'TACGCT' + 'GCGCGT'             # TACGCTGCGCGT


def build_canon34(core10):
    """Ideal canonical cassette for a given core. Both strands are indexed later."""
    return CANON_5 + core10 + CANON_3


def canon_core_from_anchor(bc34):
    """Canonical 10bp core from a WELL-FORMED 34-mer, or None.

    Kept for parity with the prior demux so Tier 1 behaviour is byte-identical.
    """
    if len(bc34) != BC_LEN:
        return None
    inner, core = bc34[6:12], bc34[CORE_SLICE[0]:CORE_SLICE[1]]
    if inner == INNER_ORIENT1:
        return core
    if inner == INNER_ORIENT2:
        return revcomp(core)
    return None


def core_from_sequence(seq):
    """Register-independent core extraction. Returns (canonical_core, note) or (None, reason).

    Locates the cassette by its invariant 12bp 5' block (scar + inner flank) rather
    than by hard-coded slice positions, so it is correct for entries that are
    shifted, padded or truncated in the workbook.

    This REPLACES the pad/trim approach, which was not just incomplete but wrong:
    padding a truncated entry and then reverse-complementing it shifts the core
    window by one base and silently returns a neighbouring 10-mer. Verified against
    ADE08M.11736, where pad-then-revcomp gave ACAGGGCGGC while the correct core is
    AACAGGGCGG.

    Canonical core is defined as the core read in orientation 1 (inner flank
    AGCGTA), matching the prior demux convention.
    """
    for s, flipped in ((seq, False), (revcomp(seq), True)):
        for five, three, orient in ((CANON_5, CANON_3, 'orient1'),
                                    (ORIENT2_5, ORIENT2_3, 'orient2')):
            i = s.find(five)
            if i == -1 or len(s) < i + 22:
                continue
            core = s[i + 12:i + 22]
            if len(core) != 10 or not set(core) <= set('ACGT'):
                continue
            # Verify the 3' flank, tolerantly. Without ANY check, a truncated entry
            # yields a core that has eaten part of the downstream flank and is
            # simply wrong (ADE08M.1409 returns TGCCCAGCTA, which is 5bp of core
            # plus 5bp of flank). With a STRICT prefix check, the large AACGCG*
            # family is wrongly rejected: those carry an intact 12bp 5' flank and
            # an in-frame core, with the damage confined to a one-base slip in the
            # 3' flank. So: allow the observed 3' block to align at offset 0 or 1
            # with a small mismatch budget, and require enough of it to be present
            # to be meaningful.
            obs3 = s[i + 22:i + 22 + len(three)]
            ok, tag3 = False, ''
            for off3 in (0, 1):
                ref3 = three[off3:]
                ncmp = min(len(obs3), len(ref3))
                if ncmp < FLANK3_MIN_COMPARE:
                    continue
                mm = sum(1 for a, b in zip(obs3[:ncmp], ref3[:ncmp]) if a != b)
                if mm <= FLANK3_MAX_MISMATCH:
                    ok = True
                    tag3 = f'_3p_off{off3}_mm{mm}_n{ncmp}'
                    break
            if not ok:
                continue
            tag = orient + ('_rc' if flipped else '') + tag3
            return (core if orient == 'orient1' else revcomp(core)), tag
    return None, 'no_invariant_flank'


gt_core       = defaultdict(set)     # animal -> {canonical 10bp cores}
gt_full       = {}                   # canonical 34-mer -> canonical core
gt_full_animal = defaultdict(set)    # canonical 34-mer -> {animals}
disposition   = []                   # per-row audit trail
have_gt       = False

if not os.path.exists(BINHUA_XLSX):
    print(f'  *** WORKBOOK NOT FOUND: {BINHUA_XLSX}')
    print('      Tier 2/3/4/5 all require ground truth. Fix the path in Cell 1.')
else:
    try:
        from openpyxl import load_workbook
        wb = load_workbook(BINHUA_XLSX, read_only=True, data_only=True)
        print(f'  workbook : {os.path.basename(BINHUA_XLSX)}')
        print(f'  sheets   : {wb.sheetnames}')
        if len(wb.sheetnames) != 5:
            print(f'  *** WARN: expected 5 animal sheets, found {len(wb.sheetnames)}. '
                  'Confirm the animal roster before trusting per-animal coverage.')

        for sh in wb.sheetnames:
            animal = sh.split()[0]
            n_ok = n_resc = n_quar = n_drop = 0
            for row in wb[sh].iter_rows(values_only=True):
                if len(row) < 3:
                    continue
                bc_name, bc = row[1], row[2]
                if not isinstance(bc, str):
                    continue
                bc = bc.strip().upper()
                if not bc or not set(bc) <= set('ACGTN'):
                    continue
                name = str(bc_name).strip() if bc_name is not None else ''

                if name in QUARANTINE_NAMES and not INCLUDE_QUARANTINED:
                    n_quar += 1
                    disposition.append({'animal': animal, 'bc_name': name, 'seq': bc,
                                        'len': len(bc), 'status': 'quarantined',
                                        'note': 'off-register family, unresolved'})
                    continue

                core, note = core_from_sequence(bc)
                if core is None:
                    n_drop += 1
                    disposition.append({'animal': animal, 'bc_name': name,
                                        'seq': bc, 'len': len(bc),
                                        'status': 'dropped', 'note': note})
                    continue

                # Rebuild the ideal cassette. For a well-formed entry this is the
                # entry itself (or its revcomp); for a malformed one this repairs
                # the flanks, which makes the T3/T4/T5 search stronger since we
                # then hunt the real cassette rather than a typo'd spreadsheet cell.
                canon34 = build_canon34(core)
                well_formed = (len(bc) == BC_LEN and
                               canon_core_from_anchor(bc) == core)
                if well_formed or (len(bc) == BC_LEN and
                                   canon_core_from_anchor(revcomp(bc)) == core):
                    status = 'parsed'
                    n_ok += 1
                else:
                    status = 'rescued'
                    note = f'{note}|rebuilt_from_core'
                    n_resc += 1

                gt_core[animal].add(core)
                gt_full[canon34] = core
                gt_full_animal[canon34].add(animal)
                disposition.append({'animal': animal, 'bc_name': name, 'seq': bc,
                                    'len': len(bc), 'status': status, 'note': note,
                                    'canon34': canon34, 'core': core})

            print(f'  {animal}: {n_ok} parsed, {n_resc} rescued, '
                  f'{n_quar} quarantined, {n_drop} dropped '
                  f'-> {len(gt_core[animal])} unique cores')
        wb.close()
        have_gt = len(gt_full) > 0
    except Exception as e:
        print(f'  *** WARN: workbook parse failed ({e}); continuing without ground truth.')

disp_df = pd.DataFrame(disposition)
if len(disp_df):
    print('\n  row disposition:')
    print(disp_df['status'].value_counts().to_string())
    save_csv(disp_df, 'gt_row_disposition.csv')
    drops = disp_df[disp_df['status'].isin(['dropped', 'quarantined'])]
    if len(drops):
        print('\n  NOT in ground truth (inspect these, do not ignore them):')
        print(drops[['animal', 'bc_name', 'seq', 'len', 'status', 'note']].to_string(index=False))

gt_core_all = set().union(*gt_core.values()) if gt_core else set()
core_to_animals = defaultdict(set)
for a, cores in gt_core.items():
    for c in cores:
        core_to_animals[c].add(a)

print(f'\n  ground truth: {len(gt_full)} canonical 34-mers, '
      f'{len(gt_core_all)} unique cores, {len(gt_core)} animals')

# Tier-2 flank probes, rebuilt from the PARSED set so a rescued entry contributes
FLANK5, FLANK3 = {}, {}
if have_gt:
    pre, suf = Counter(), Counter()
    for bc in gt_full:
        pre[bc[:12]] += 1
        suf[bc[-12:]] += 1
    for seq, n in pre.most_common():
        if seq.startswith('ACGCGC') and n >= 20:
            r = 'asis' if seq[6:12] == INNER_ORIENT1 else ('rc' if seq[6:12] == INNER_ORIENT2 else None)
            if r:
                FLANK5[seq] = r
    for seq, n in suf.most_common():
        if seq.endswith('GCGCGT') and n >= 20:
            eq = revcomp(seq[0:6])
            r = 'asis' if eq == INNER_ORIENT1 else ('rc' if eq == INNER_ORIENT2 else None)
            if r:
                FLANK3[seq] = r
    print(f'  tier-2 5-prime flanks: {FLANK5}')
    print(f'  tier-2 3-prime flanks: {FLANK3}')


# %% ============================================================
# CELL 3 — SEARCH INDEX (both strands)
# ============================================================
hr('CELL 3 — search index construction')

# Every GT barcode indexed on BOTH strands. A read can carry the cassette either
# way round, and we resolve back to the canonical core at match time.
gt_strands = {}                      # observed 34-mer (either strand) -> canonical core
for canon34, core in gt_full.items():
    gt_strands[canon34] = core
    gt_strands[revcomp(canon34)] = core
print(f'  full 34-mers indexed (both strands): {len(gt_strands)}')

# --- cassette anchoring ----------------------------------------------------
# RUNTIME NOTE. The naive implementations of T3/T4/T5 (scan every GT barcode per
# read; slide a seed over every read position; try every partial overlap length)
# cost hundreds of dict lookups or thousands of substring tests PER READ. Across
# ~120M prefiltered reads that is days of wall clock. Instead we anchor.
#
# Every cassette, in either written orientation, carries four invariant 6-mers at
# FIXED offsets within the 34bp window:
#
#   canonical     ACGCGC AGCGTA <core10> CAGCTA GCGCGT
#   orientation 2 ACGCGC TAGCTG <core10> TACGCT GCGCGT
#                 ^0     ^6               ^22    ^28
#
# So one occurrence of any of those 6-mers in a read implies a cassette start
# position. We find those occurrences once, derive a handful of candidate starts,
# and run T3/T4/T5 only at those starts. That turns per-read cost from hundreds of
# lookups into a few, with no loss of sensitivity that matters: a cassette so
# degraded that none of its 24bp of invariant flank survives intact for 6bp is not
# something we could confidently call anyway.
ANCHOR_OFFSETS = {
    'ACGCGC': 0, 'AGCGTA': 6, 'TAGCTG': 6,
    'CAGCTA': 22, 'TACGCT': 22, 'GCGCGT': 28,
}

# --- TIER 5 seed index -----------------------------------------------------
# Pigeonhole: with <= MAX_MISMATCH mismatches across 34bp, splitting into
# MAX_MISMATCH+1 disjoint chunks guarantees at least one chunk survives exactly.
SEED_N   = MAX_MISMATCH + 1
SEED_LEN = BC_LEN // SEED_N
seed_index = defaultdict(set)        # seed -> {(observed34, seed_offset)}
if ENABLE_TIER5_MISMATCH:
    for obs34 in gt_strands:
        for i in range(SEED_N):
            off = i * SEED_LEN
            seed_index[obs34[off:off + SEED_LEN]].add((obs34, off))
    print(f'  tier-5 seed index: {len(seed_index)} seeds of {SEED_LEN}bp '
          f'({SEED_N} per barcode, tolerates <= {MAX_MISMATCH} mismatches)')

# --- TIER 4 partial index --------------------------------------------------
# A read ending mid-cassette shows a PREFIX of the barcode at its 3' end.
# A read starting mid-cassette shows a SUFFIX of the barcode at its 5' end.
prefix_index = defaultdict(set)      # barcode prefix -> {observed34}
suffix_index = defaultdict(set)      # barcode suffix -> {observed34}
if ENABLE_TIER4_PARTIAL:
    for obs34 in gt_strands:
        for k in range(MIN_PARTIAL_OVERLAP, BC_LEN):
            prefix_index[obs34[:k]].add(obs34)
            suffix_index[obs34[-k:]].add(obs34)
    print(f'  tier-4 partial index: {len(prefix_index)} prefixes / '
          f'{len(suffix_index)} suffixes, overlap >= {MIN_PARTIAL_OVERLAP}bp')


def core_bases_covered(lo, hi):
    """How many of the 10 core bases fall inside cassette window [lo, hi)."""
    c_lo, c_hi = CORE_SLICE
    return max(0, min(hi, c_hi) - max(lo, c_lo))


def candidate_starts(s):
    """Cassette start positions implied by invariant 6-mer occurrences in s.

    Starts may be negative or run past the end of s; those are exactly the
    partial-overlap cases TIER 4 exists to catch, so they are kept.
    """
    starts = set()
    for six, off in ANCHOR_OFFSETS.items():
        p = s.find(six)
        while p != -1:
            starts.add(p - off)
            p = s.find(six, p + 1)
    return starts


def match_tiers(seq):
    """Return (core, tier, gt_flag, detail) for the best (lowest) tier that fires.

    Tier order is deliberate: 1 and 2 reproduce the prior demux exactly so the
    numbers are directly comparable; 3, 4, 5 are the new layers.
    """
    strands = (seq, revcomp(seq))

    # --- TIER 1: full flank anchor, core trusted even if novel ---------------
    for s in strands:
        m = ANCHOR_RE.search(s)
        if m:
            c = canon_core_from_anchor(m.group(0))
            if c:
                if c in gt_core_all:
                    return (c, 1, True, 'anchor')
                rcc = revcomp(c)
                if rcc in gt_core_all:
                    return (rcc, 1, True, 'anchor_rc')
                return (c, 1, False, 'anchor_novel')

    # --- TIER 2: single flank + core, ground-truth gated ---------------------
    for s in strands:
        for fl, rule in FLANK5.items():
            i = s.find(fl)
            if i != -1:
                raw = s[i + 12:i + 22]
                if len(raw) == 10 and set(raw) <= set('ACGT'):
                    cand = raw if rule == 'asis' else revcomp(raw)
                    if cand in gt_core_all:
                        return (cand, 2, True, 'flank5')
        for fl, rule in FLANK3.items():
            i = s.find(fl)
            if i >= 10:
                raw = s[i - 10:i]
                if len(raw) == 10 and set(raw) <= set('ACGT'):
                    cand = raw if rule == 'asis' else revcomp(raw)
                    if cand in gt_core_all:
                        return (cand, 2, True, 'flank3')

    # --- TIER 3/4/5 all run at anchored candidate starts ---------------------
    for s in strands:
        L = len(s)
        for start in candidate_starts(s):
            inside = (start >= 0 and start + BC_LEN <= L)

            # TIER 3: literal full 34-mer from the workbook. A strict subset of
            # TIER 1 by construction; carried as the committed deliverable and as
            # a cross-check that TIER 1 and the ground truth agree.
            if ENABLE_TIER3_EXACT_GT and inside:
                win = s[start:start + BC_LEN]
                hit = gt_strands.get(win)
                if hit is not None:
                    return (hit, 3, True, 'exact_gt34')

            # TIER 4: cassette runs off one end of the read.
            #
            # AMBIGUITY GATE. A partial view shows only some of the 10 core bases,
            # so two GT barcodes differing solely in the hidden bases are
            # indistinguishable. Returning an arbitrary one silently assigns the
            # wrong lineage, which is worse than recovering nothing: the entire
            # point of this is to say WHICH barcode is in a cell. Verified in
            # synthetic testing, where a 20bp overlap confused ACGGGCCTAC with
            # ACGGGCCTAG. So we accept only when the visible window resolves to a
            # single core, and count the rest as ambiguous.
            if ENABLE_TIER4_PARTIAL and not inside:
                if start < 0:
                    k = BC_LEN + start           # visible bases at the read head
                    if (MIN_PARTIAL_OVERLAP <= k <= L and
                            core_bases_covered(BC_LEN - k, BC_LEN) >= MIN_CORE_BASES):
                        cands = suffix_index.get(s[:k], ())
                        cores = {gt_strands[o] for o in cands}
                        if len(cores) == 1:
                            return (cores.pop(), 4, True, f'partial5p_{k}')
                        if len(cores) > 1:
                            return (None, -1, False, f'ambiguous5p_{k}_n{len(cores)}')
                else:
                    k = L - start                # visible bases at the read tail
                    if (k >= MIN_PARTIAL_OVERLAP and
                            core_bases_covered(0, k) >= MIN_CORE_BASES):
                        cands = prefix_index.get(s[start:], ())
                        cores = {gt_strands[o] for o in cands}
                        if len(cores) == 1:
                            return (cores.pop(), 4, True, f'partial3p_{k}')
                        if len(cores) > 1:
                            return (None, -1, False, f'ambiguous3p_{k}_n{len(cores)}')

            # TIER 5: full-length window, <= MAX_MISMATCH against the closed set.
            #
            # SCOPE NOTE (verified in synthetic testing). Tier 5 fires far less
            # often than it looks like it should, because a read with INTACT
            # flanks and a mutated core is already caught by Tier 1, which trusts
            # the anchor and reports the mutated core as novel. So Tier 5's real
            # job is the narrower case where the mutation lands in the FLANK
            # itself and breaks the Tier 1 anchor. That is a small population, and
            # it is the honest reason not to expect much from this tier.
            if ENABLE_TIER5_MISMATCH and inside:
                win = s[start:start + BC_LEN]
                seen = set()
                for i in range(SEED_N):
                    off = i * SEED_LEN
                    for obs34, soff in seed_index.get(win[off:off + SEED_LEN], ()):
                        if soff != off or obs34 in seen:
                            continue
                        seen.add(obs34)
                        mm = 0
                        for a, b in zip(win, obs34):
                            if a != b:
                                mm += 1
                                if mm > MAX_MISMATCH:
                                    break
                        if mm <= MAX_MISMATCH:
                            return (gt_strands[obs34], 5, True, f'mm{mm}')

    return (None, 0, False, '')


# %% ============================================================
# CELL 4 — whitelist + embedded object lookup
# ============================================================
hr('CELL 4 — whitelist + object lookup')

whitelist = set()
if os.path.exists(WHITELIST):
    with open(WHITELIST) as fh:
        for line in fh:
            b = line.strip()
            if b:
                whitelist.add(b)
    print(f'  whitelist: {len(whitelist):,} barcodes')
else:
    print(f'  *** WARN: whitelist not found ({WHITELIST}); CB correction disabled.')


def correct_cb(cb):
    """Exact, then 1-mismatch, against the whitelist. Only run on extraction hits."""
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


def raw_cb_from_name(name):
    m = re.search(r'[ACGT]{16}', name)
    return m.group(0) if m else None


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
        ct, ts = str(ct), str(ts)
        if ct == 'T cells' and ts not in ('nan', 'not T', 'Unassigned', '', 'None'):
            return ts
        return ct

    for name, r in obs.iterrows():
        cb = raw_cb_from_name(name)
        if cb is None:
            continue
        ct = r[COL_CELLTYPE] if extra_present[COL_CELLTYPE] else None
        ts = r[COL_TSUB] if extra_present[COL_TSUB] else None
        lookup[(str(r[COL_LIBRARY]), cb)] = {
            'obs_name': name,
            'animal':    r[COL_ANIMAL]    if extra_present[COL_ANIMAL]    else None,
            'tissue':    r[COL_TISSUE]    if extra_present[COL_TISSUE]    else None,
            'timepoint': r[COL_TIMEPOINT] if extra_present[COL_TIMEPOINT] else None,
            'celltype': ct,
            'final_celltype': final_ct(ct, ts),
            't_subtype': ts,
            'subspecies': r[COL_SUBSP]   if extra_present[COL_SUBSP]   else None,
            'shiv_pos':   r[COL_SHIVPOS] if extra_present[COL_SHIVPOS] else None,
            'shiv_umi':   r[COL_SHIVUMI] if extra_present[COL_SHIVUMI] else None,
        }
    print(f'  lookup built: {len(lookup):,} keys')
    print(f'  letter->libval: {letter_to_libval}')
else:
    print(f'  *** WARN: object not found ({OBJECT}); linkage skipped, read tables still emitted.')


# %% ============================================================
# CELL 5 — FASTQ manifest + R1 geometry probe
# ============================================================
hr('CELL 5 — FASTQ manifest + R1 geometry')

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
    manifest.append({'lib': L, 'group': 'INF' if L in INFECTED else 'neg',
                     'r1': r1, 'r2': r2})
    print(f'  {L} [{manifest[-1]["group"]}] R1={len(r1)} R2={len(r2)}  {os.path.basename(r1[0])}')


def probe_read_len(path, n=2000):
    """Modal read length from the first n records. Warn-and-continue on failure."""
    try:
        lens = Counter()
        with gzip.open(path, 'rt') as fh:
            for i, line in enumerate(fh):
                if i >= n * 4:
                    break
                if i % 4 == 1:
                    lens[len(line.strip())] += 1
        return lens.most_common(1)[0][0] if lens else None
    except Exception as e:
        print(f'      WARN: read-length probe failed on {os.path.basename(path)} ({e})')
        return None


R1_TAIL_START = CB_LEN + UMI_LEN + R1_TSO_LEN
r1_geom = []
for row in manifest:
    l1 = probe_read_len(row['r1'][0])
    l2 = probe_read_len(row['r2'][0])
    tail = (l1 - R1_TAIL_START) if l1 else 0
    usable = bool(SCAN_R1_TAIL and l1 and tail >= R1_MIN_TAIL)
    row['r1_len'], row['r2_len'], row['scan_r1'] = l1, l2, usable
    r1_geom.append({'lib': row['lib'], 'r1_len': l1, 'r2_len': l2,
                    'r1_cdna_tail_bp': max(0, tail), 'will_scan_r1': usable})
    print(f'  {row["lib"]}: R1={l1}bp R2={l2}bp | R1 cDNA tail past CB+UMI+TSO = '
          f'{max(0, tail)}bp | scan R1: {usable}')

geom_df = pd.DataFrame(r1_geom)
if len(geom_df):
    save_csv(geom_df, 'read_geometry.csv')
    if not geom_df['will_scan_r1'].any():
        print('\n  NOTE: no library has usable R1 cDNA tail. R1 is CB+UMI+TSO only, '
              'which confirms the prior demux lost nothing by reading R2 alone.')
    else:
        print('\n  NOTE: R1 carries unscanned cDNA in at least one library. This is '
              'read pool the prior demux never looked at.')


# %% ============================================================
# CELL 6 — SCAN
# ============================================================
hr('CELL 6 — ground-truth search scan (the heavy step)')
print(f'  tiers enabled: T1 anchor, T2 flank+GT'
      f'{", T3 exact-GT34" if ENABLE_TIER3_EXACT_GT else ""}'
      f'{", T4 partial" if ENABLE_TIER4_PARTIAL else ""}'
      f'{", T5 mismatch" if ENABLE_TIER5_MISMATCH else ""}')
print(f'  prefilter seeds: {[s.decode() for s in PREFILTER_SEEDS]}')
print('  RUNTIME: T4/T5 are the expensive tiers. If this is not progressing, set')
print('           ENABLE_TIER4/5 = False for a first pass, then re-enable.\n')


def open_stream(files):
    return subprocess.Popen(DECOMP.split() + files, stdout=subprocess.PIPE, bufsize=1 << 22)


def scan_library(row):
    L = row['lib']
    p1, p2 = open_stream(row['r1']), open_stream(row['r2'])
    f1, f2 = p1.stdout, p2.stdout
    hits, n, n_prefilter, n_ambig = [], 0, 0, 0
    scan_r1 = row.get('scan_r1', False)
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
                print(f'    {L}: {n:,} reads | {n_prefilter:,} prefiltered | {len(hits)} hits')

            s1s = s1.strip()
            targets = [('R2', s2.strip())]
            if scan_r1:
                targets.append(('R1', s1s[R1_TAIL_START:]))

            for src, raw in targets:
                if not any(sd in raw for sd in PREFILTER_SEEDS):
                    continue
                n_prefilter += 1
                core, tier, gt, detail = match_tiers(raw.decode(errors='ignore'))
                if tier == -1:
                    # partial match consistent with >1 barcode; counted, not used
                    n_ambig += 1
                    continue
                if tier == 0:
                    continue
                cb_raw = s1s[:CB_LEN].decode(errors='ignore')
                umi = s1s[CB_LEN:CB_LEN + UMI_LEN].decode(errors='ignore')
                cb = correct_cb(cb_raw)
                hits.append({'lib': L, 'group': row['group'], 'read_src': src,
                             'tier': tier, 'detail': detail, 'core': core,
                             'gt': int(gt), 'cb_raw': cb_raw, 'cb': cb or '', 'umi': umi})
    finally:
        for pp in (p1, p2):
            try:
                pp.stdout.close(); pp.wait()
            except Exception:
                pass
    return hits, n, n_prefilter, n_ambig


all_hits = []
for row in manifest:
    hits, n, npf, namb = scan_library(row)
    all_hits.extend(hits)
    tc = Counter(h['tier'] for h in hits)
    with_cb = sum(bool(h['cb']) for h in hits)
    print(f'  {row["lib"]} [{row["group"]}]: {n:,} reads | {npf:,} prefiltered | '
          f'{len(hits)} barcoded (T1 {tc[1]}, T2 {tc[2]}, T3 {tc[3]}, T4 {tc[4]}, T5 {tc[5]}) '
          f'| {namb} ambiguous-partial (discarded) | {with_cb} CB-valid')

reads_df = pd.DataFrame(all_hits)
if len(reads_df):
    save_csv(reads_df, 'gt_search_barcoded_reads.csv')
else:
    print('\n  No barcoded reads recovered in any library.')


# %% ============================================================
# CELL 7 — SPECIFICITY FLOOR (per tier — this is the gate on T4/T5)
# ============================================================
hr('CELL 7 — negative-library specificity floor, per tier')

if len(reads_df):
    floor = (reads_df.groupby(['tier', 'group']).size().unstack(fill_value=0)
             .reindex(columns=['INF', 'neg'], fill_value=0).reset_index())
    floor['fp_rate_vs_inf'] = np.where(floor['INF'] > 0,
                                       (floor['neg'] / floor['INF']).round(3), np.nan)
    print(floor.to_string(index=False))
    save_csv(floor, 'gt_search_specificity_floor.csv')

    dirty = floor[(floor['neg'] > 0)]['tier'].tolist()
    if dirty:
        print(f'\n  *** Tiers with a non-zero negative floor: {dirty}')
        print('      T1/T2 dirty  -> contamination, stop and inspect.')
        print('      T4/T5 dirty  -> expected; these are the permissive tiers. Treat')
        print('                      their infected counts as floor-corrected only.')
    else:
        print('\n  Negative floor CLEAN across all tiers.')
else:
    floor = pd.DataFrame()


# %% ============================================================
# CELL 8 — molecules, cells, linkage
# ============================================================
hr('CELL 8 — molecules, cells, object linkage')

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
            'lib': lib, 'group': 'INF' if lib in INFECTED else 'neg',
            'sample': LIB_NAME[lib], 'cb': cb,
            'cores': ';'.join(cores), 'n_cores': len(cores),
            'n_umis': g['umi_key'].nunique(),
            'best_tier': int(g['tier'].min()),
            'tiers': ';'.join(str(t) for t in sorted(set(g['tier']))),
            'read_srcs': ';'.join(sorted(set(g['read_src']))),
            'in_object': meta is not None,
            'animal':    meta['animal']     if meta else '',
            'tissue':    meta['tissue']     if meta else '',
            'timepoint': meta['timepoint']  if meta else '',
            'final_celltype': ct,
            't_subtype':  meta['t_subtype']  if meta else '',
            'subspecies': meta['subspecies'] if meta else '',
            'shiv_pos':   meta['shiv_pos']   if meta else '',
            'shiv_umi':   meta['shiv_umi']   if meta else '',
            'is_target': bool(meta) and str(ct) not in NON_TARGET,
            'barcode_animals_binhua': ';'.join(sorted(set().union(
                *[core_to_animals.get(c, set()) for c in cores]))) if cores else '',
        })
    cells_df = pd.DataFrame(rows)
    n = len(cells_df)
    if n:
        print(f'  barcode-positive cells: {n} | linked to object: {matched} ({100 * matched / n:.0f}%)')
        if matched == 0:
            print('  *** 0% linkage — check letter->libval / obs_names before trusting the map.')
        save_csv(cells_df, 'gt_search_barcode_by_cell.csv')
    else:
        print('  no CB-valid barcode-positive cells')


# %% ============================================================
# CELL 9 — reservoir map, diversity, Binhua coverage
# ============================================================
hr('CELL 9 — reservoir map (animal x timepoint x tissue x final cell type)')

if len(cells_df):
    inf = cells_df[cells_df['group'] == 'INF']
    target = inf[inf['in_object'] & inf['is_target']]
    print(f'  infected barcode+ cells: {len(inf)} | object-linked target cells: {len(target)}')

    if len(target):
        res = (target.groupby(['animal', 'timepoint', 'tissue', 'final_celltype'])
               .agg(cells=('cb', 'nunique'),
                    distinct_barcodes=('cores', lambda s: len(set(';'.join(s).split(';')))),
                    umis=('n_umis', 'sum'),
                    best_tier=('best_tier', 'min')).reset_index()
               .sort_values(['animal', 'timepoint', 'cells'], ascending=[True, True, False]))
        print('\n  RESERVOIR MAP:')
        print(res.to_string(index=False))
        save_csv(res, 'gt_search_reservoir_map.csv')

        div = (target.groupby(['animal', 'timepoint'])
               .agg(cells=('cb', 'nunique'),
                    distinct_barcodes=('cores', lambda s: len(set(';'.join(s).split(';'))))).reset_index())
        print('\n  BARCODE DIVERSITY over time:')
        print(div.to_string(index=False))
        save_csv(div, 'gt_search_diversity_by_time.csv')

        if have_gt:
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
            save_csv(cov, 'gt_search_binhua_coverage.csv')
    else:
        print('  No object-linked target cells to map.')

    tally = (cells_df.groupby(['lib', 'sample', 'group'])
             .agg(cells=('cb', 'nunique'), linked=('in_object', 'sum'),
                  target=('is_target', 'sum'),
                  distinct_barcodes=('cores', lambda s: len(set(';'.join(s).split(';'))))).reset_index())
    print('\n  per-library tally:')
    print(tally.to_string(index=False))
    save_csv(tally, 'gt_search_cells_by_library.csv')


# %% ============================================================
# CELL 10 — tier-yield figure
# ============================================================
hr('CELL 10 — tier yield figure')

if len(reads_df):
    try:
        piv = (reads_df.groupby(['lib', 'tier']).size().unstack(fill_value=0)
               .reindex(index=[l for l in LIBS if l in set(reads_df['lib'])], fill_value=0))
        fig, ax = plt.subplots(figsize=(16, 10))
        bottom = np.zeros(len(piv))
        for t in sorted(piv.columns):
            vals = piv[t].values
            ax.bar(piv.index, vals, bottom=bottom, label=f'Tier {t}',
                   color=HEX.get(f't{t}', '#999999'), edgecolor='white', linewidth=2)
            bottom += vals
        ax.set_xlabel('Library', fontsize=FONT_BASE + 4)
        ax.set_ylabel('Barcoded reads', fontsize=FONT_BASE + 4)
        ax.set_title('vpr barcode reads recovered by tier', fontsize=FONT_BASE + 6)
        ax.tick_params(labelsize=FONT_BASE)
        ax.legend(fontsize=FONT_BASE - 2, frameon=False)
        for sp in ('top', 'right'):
            ax.spines[sp].set_visible(False)
        save_fig(fig, 'gt_search_tier_yield')
        if not HEADLESS:
            plt.show()
        plt.close(fig)
    except Exception as e:
        print(f'  WARN: figure failed ({e})')


# %% ============================================================
# CELL 11 — SUMMARY + VERDICT
# ============================================================
hr('CELL 11 — SUMMARY + VERDICT')

n_reads = len(reads_df)
tier_counts = Counter(reads_df['tier']) if n_reads else Counter()
n_cells = len(cells_df)
n_target = int((cells_df['group'].eq('INF') & cells_df['in_object'] & cells_df['is_target']).sum()) if n_cells else 0
n_neg = int(cells_df['group'].eq('neg').sum()) if n_cells else 0
n_bc = 0
if n_cells:
    t = cells_df[cells_df['group'].eq('INF') & cells_df['in_object'] & cells_df['is_target']]
    n_bc = len(set(';'.join(t['cores']).split(';'))) if len(t) else 0

# Reads attributable ONLY to the new tiers. This is the whole question.
n_legacy = tier_counts[1] + tier_counts[2]
n_new = tier_counts[3] + tier_counts[4] + tier_counts[5]

print(f"""
  Ground truth used                 : {len(gt_full)} barcodes / {len(gt_core_all)} cores
  Barcoded reads (all tiers)        : {n_reads}
    Tier 1 anchor                   : {tier_counts[1]}
    Tier 2 flank+GT                 : {tier_counts[2]}
    Tier 3 exact GT 34-mer          : {tier_counts[3]}   (subset of T1 by construction)
    Tier 4 partial at read end      : {tier_counts[4]}   <- NEW
    Tier 5 mismatch <= {MAX_MISMATCH}          : {tier_counts[5]}   <- NEW
  Legacy-equivalent reads (T1+T2)   : {n_legacy}
  Reads only the new tiers found    : {tier_counts[4] + tier_counts[5]}
  Barcode-positive cells (all)      : {n_cells}
  Negative-library cells (floor)    : {n_neg}
  Object-linked TARGET cells (INF)  : {n_target}
  Distinct barcodes in those cells  : {n_bc}
""")

verdict_path = os.path.join(OUT_DIR, 'VERDICT.txt')
gained = tier_counts[4] + tier_counts[5]
with open(verdict_path, 'w') as fh:
    if n_neg > 0 and any(floor[floor['tier'].isin([1, 2])]['neg'] > 0) if len(floor) else False:
        msg = ('CAUTION: barcode+ cells in negative libraries at Tier 1/2. The floor is '
               'not clean at the strict tiers; inspect for contamination before trusting '
               'any infected count.')
    elif n_target >= 20:
        msg = (f'GO: {n_target} object-linked target cells carrying {n_bc} distinct barcodes. '
               'Productionize the reservoir map.')
    elif n_target >= 5:
        msg = (f'MARGINAL: {n_target} target cells, {n_bc} barcodes. Reportable as proof of '
               'principle, thin for population dynamics.')
    elif gained == 0:
        msg = (f'CEILING CONFIRMED: the permissive tiers (partial-end overlap >= '
               f'{MIN_PARTIAL_OVERLAP}bp, <= {MAX_MISMATCH} mismatches against a closed '
               f'{len(gt_full)}-member ground truth) added ZERO reads over the flank anchor. '
               'Two independent methods, one deliberately more sensitive than the other, '
               'converge on the same count. The limit is 5-prime capture geometry, not '
               'extraction. Targeted VpxF1 amplicon enrichment off the archived cDNA is the '
               'only route that raises the read count. Bulk barcode dynamics remain fully '
               'available from the Binhua SGA/NGS file.')
    else:
        msg = (f'MARGINAL GAIN: the permissive tiers added {gained} read(s) over the flank '
               f'anchor, for {n_target} object-linked target cell(s). Check the per-tier '
               'negative floor before claiming these: if T4/T5 also fire in the negatives, '
               'the gain is false-positive and the ceiling stands.')
    print('  ' + msg)
    fh.write(msg + '\n')

print(f'\n  wrote {verdict_path}')
print(f'  all outputs: {OUT_DIR}')
print('=' * 90)
print('vpr GROUND-TRUTH SEARCH COMPLETE — read-only. Nothing modified.')
print('=' * 90)
