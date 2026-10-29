#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Rhesus SHIV — R1 Structure Check (READ-ONLY)
============================================
Establishes, empirically, what is actually in R1 past the cell barcode, and
whether it is worth searching for the vpr cassette.

WHY THIS EXISTS
---------------
shiv_read_length_check.py found R1 = 150bp in every GEX library. For 10x 5'
chemistry R1 is not a barcode-only read; it is

    [CB 16][UMI 12][TSO 13]  then cDNA, sense strand, from the molecule's 5' end

which at 150bp leaves ~109bp of real transcript sequence per read. Every script
in this repo discards it:

    shiv_vpr_barcode_feasibility.py   cb = s1[:CB_LEN]           (rest dropped)
    SHIV_vpr_barcode_gt_search.py     same
    shiv_starsolo_recovery.sh         --soloCBlen 16 --soloUMIlen 12, and
                                      CB_UMI_Simple reads nothing past base 28

If that cDNA is real, it roughly doubles the searchable sequence per read pair,
AND a cassette hit in R1 puts the cell barcode and the viral barcode in the SAME
molecule 41 bases apart, which needs no pairing assumption and no barcode-
orientation decision. That is the strongest barcode-to-cell evidence this
dataset can produce.

WHAT THIS DOES NOT ASSUME
-------------------------
The 16/12/13 layout is the expected GEM-X 5' v3 geometry, not a measurement.
check_sample_whitelist.py already established CB=16 and UMI=12 empirically from
the BAMs; this script establishes the TSO offset and the cDNA start the same way,
by finding the TSO motif in the reads rather than trusting the product spec. If
the modal offset is not 28 the script says so and reports what it actually is.

THE REAL RISK IS POLY-G, NOT GEOMETRY
-------------------------------------
These are NovaSeq X reads (read IDs start LH00973; see
shiv_merge_barcode_diagnostic.py). Two-colour chemistry calls G on empty or
dark cycles, so a short insert reads as a poly-G tail. If the post-TSO portion
of R1 is mostly poly-G then the 109bp is an artifact and there is nothing to
search. CELL 4 measures that directly: poly-G fraction, base composition, and
a 4-mer entropy per position block. A verdict that says "poly-G dominated" kills
the R1 arm cleanly and cheaply, which is the point of running this first.

Read-only. Spyder cells (# %%), config in Cell 1, warn-and-continue.

Author: Jake Lehle
Date: October 2026
"""

# %% ============================================================
# CELL 1 — CONFIG AND DIALS
# ============================================================
import os
import glob
import gzip
import math
from collections import Counter, defaultdict

import pandas as pd

# --- paths -----------------------------------------------------------------
RAW_BASE  = '/master/zwallis/WORKING/SC/10X_RAW'
WORKING   = '/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV'
WHITELIST = os.path.join(WORKING, 'shiv_starsolo', 'ref',
                         '3M-5pgex-jan-2023.txt')
OUT_DIR   = os.path.join(WORKING, 'annotation_output', 'r1_structure_check')

# --- libraries -------------------------------------------------------------
LIBS     = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
INFECTED = {'D', 'E', 'F', 'H'}
LIB_DIR  = 'GEX_LIBRARY_{L}'

# --- expected geometry (GEM-X 5' v3; CB/UMI confirmed empirically by
#     check_sample_whitelist.py, TSO offset is what THIS script measures) -----
CB_LEN      = 16
UMI_LEN     = 12
TSO         = 'TTTCTTATATGGG'     # 10x 5' template-switch oligo
TSO_EXPECT_OFFSET = CB_LEN + UMI_LEN      # 28, the offset to confirm or refute
TSO_SEARCH_MAX    = 60            # how far into R1 to look for the TSO
TSO_MAX_MISMATCH  = 2             # tolerate sequencing error in the TSO

# --- sampling --------------------------------------------------------------
SAMPLE_READS = 200_000            # per library; seconds each, raise if you like

# --- poly-G / complexity dials ---------------------------------------------
POLYG_RUN       = 15              # a G-run this long counts as a poly-G tail
LOWCOMP_ENTROPY = 1.2             # bits per base below this = low complexity
BLOCK           = 20              # cDNA is profiled in blocks of this many bp

# --- whitelist validation --------------------------------------------------
WHITELIST_CHECK_N = 20_000        # barcodes to test against the whitelist
LOAD_WHITELIST    = True          # False to skip (big file, slow-ish read)

os.makedirs(OUT_DIR, exist_ok=True)
pd.set_option('display.width', 170)
pd.set_option('display.max_columns', 60)

_COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


def revcomp(s):
    return s.translate(_COMP)[::-1]


def hr(title):
    print("\n" + "=" * 84)
    print(title)
    print("=" * 84)


def find_motif(seq, motif, max_off, max_mm):
    """Leftmost offset where motif matches seq with <= max_mm mismatches.
    Returns (offset, n_mismatch) or (None, None). Plain loop on purpose: the
    motif is 13bp and the window is 60, so this is cheap and obvious."""
    m = len(motif)
    best = (None, None)
    for off in range(0, min(max_off, len(seq) - m) + 1):
        mm = 0
        for a, b in zip(seq[off:off + m], motif):
            if a != b:
                mm += 1
                if mm > max_mm:
                    break
        if mm <= max_mm:
            if best[0] is None or mm < best[1]:
                best = (off, mm)
            if mm == 0:
                return best
    return best


def shannon(seq):
    """Bits per base over ACGT; N ignored. 2.0 = uniform, 0.0 = homopolymer."""
    c = Counter(b for b in seq if b in 'ACGT')
    n = sum(c.values())
    if n == 0:
        return 0.0
    return -sum((v / n) * math.log2(v / n) for v in c.values())


def max_homopolymer(seq, base):
    best = run = 0
    for b in seq:
        run = run + 1 if b == base else 0
        best = max(best, run)
    return best


def iter_fastq_seqs(path, n_max):
    """Yield up to n_max sequence lines. Warn-and-continue on a bad stream."""
    try:
        with gzip.open(path, 'rt') as fh:
            n = 0
            for i, line in enumerate(fh):
                if i % 4 == 1:
                    yield line.rstrip('\n')
                    n += 1
                    if n >= n_max:
                        return
    except Exception as e:                       # noqa: BLE001
        print(f"    WARN  {os.path.basename(path)}: {type(e).__name__}: {e}")


# %% ============================================================
# CELL 2 — DISCOVER R1 FILES
# ============================================================
hr("CELL 2 — R1 discovery")

r1_files = {}
for L in LIBS:
    d = os.path.join(RAW_BASE, LIB_DIR.format(L=L))
    if not os.path.isdir(d):
        print(f"    WARN  {L}: directory not found ({d})")
        continue
    fs = sorted(glob.glob(os.path.join(d, '*_R1_*.fastq.gz')))
    if not fs:
        print(f"    WARN  {L}: no R1 files in {d}")
        continue
    r1_files[L] = fs
    print(f"  {L}: {len(fs)} file(s), "
          f"{sum(os.path.getsize(p) for p in fs) / 1e9:.1f} GB")

if not r1_files:
    print("\n  Nothing to measure. Check RAW_BASE / LIB_DIR.")


# %% ============================================================
# CELL 3 — WHERE DOES THE TSO SIT, AND THEREFORE WHERE DOES cDNA START
# ============================================================
hr(f"CELL 3 — TSO offset (expect {TSO_EXPECT_OFFSET}; measuring, not assuming)")
print(f"  TSO motif : {TSO}  ({len(TSO)}bp, <= {TSO_MAX_MISMATCH} mismatches)")

tso_rows = []
tso_offset_by_lib = {}

for L, fs in r1_files.items():
    off_hist = Counter()
    mm_hist = Counter()
    n = 0
    per_file = max(1, SAMPLE_READS // len(fs))
    for p in fs:
        for seq in iter_fastq_seqs(p, per_file):
            n += 1
            off, mm = find_motif(seq, TSO, TSO_SEARCH_MAX, TSO_MAX_MISMATCH)
            if off is None:
                off_hist['none'] += 1
            else:
                off_hist[off] += 1
                mm_hist[mm] += 1
    hits = n - off_hist.get('none', 0)
    numeric = {k: v for k, v in off_hist.items() if k != 'none'}
    if numeric:
        mode_off = max(numeric, key=numeric.get)
        mode_pct = 100.0 * numeric[mode_off] / n
    else:
        mode_off, mode_pct = None, 0.0
    tso_offset_by_lib[L] = mode_off
    tso_rows.append({
        'lib': L, 'group': 'INF' if L in INFECTED else 'neg',
        'n_reads': n, 'tso_found': hits,
        'tso_found_pct': round(100.0 * hits / n, 2) if n else 0.0,
        'modal_offset': mode_off,
        'modal_pct': round(mode_pct, 2),
        'cdna_start': (mode_off + len(TSO)) if mode_off is not None else None,
        'exact_match_pct': (round(100.0 * mm_hist.get(0, 0) / hits, 2)
                            if hits else 0.0),
    })
    top = ', '.join(f"{k}:{v}" for k, v in
                    sorted(numeric.items(), key=lambda kv: -kv[1])[:4])
    print(f"  {L}: {hits:,}/{n:,} reads carry the TSO "
          f"({100.0 * hits / n if n else 0:.1f}%)  top offsets -> {top}")

tso = pd.DataFrame(tso_rows)
if not tso.empty:
    tso.to_csv(os.path.join(OUT_DIR, 'r1_tso_offsets.csv'), index=False)
    print("\n" + tso.to_string(index=False))
    print("\n  saved: r1_tso_offsets.csv")

    offs = set(o for o in tso['modal_offset'] if o is not None)
    if offs == {TSO_EXPECT_OFFSET}:
        CDNA_START = TSO_EXPECT_OFFSET + len(TSO)
        print(f"\n  CONFIRMED: TSO at offset {TSO_EXPECT_OFFSET} in every "
              f"library. cDNA starts at base {CDNA_START} (0-based).")
    elif len(offs) == 1:
        o = offs.pop()
        CDNA_START = o + len(TSO)
        print(f"\n  *** TSO offset is {o}, NOT the expected "
              f"{TSO_EXPECT_OFFSET}. cDNA starts at base {CDNA_START}.")
        print("      Use the measured value everywhere downstream; do not "
              "assume CB+UMI=28.")
    else:
        CDNA_START = None
        print(f"\n  *** TSO offset is NOT uniform across libraries: "
              f"{sorted(offs)}")
        print("      Handle per library; a single global offset would be "
              "wrong for some.")
    if not tso.empty and tso['tso_found_pct'].min() < 50:
        print(f"      NOTE: lowest TSO detection is "
              f"{tso['tso_found_pct'].min():.1f}%. A low rate can mean "
              "trimmed reads or a different bead layout, not necessarily "
              "a problem; CELL 5 weighs it.")
else:
    CDNA_START = None


# %% ============================================================
# CELL 4 — IS THE POST-TSO SEQUENCE REAL, OR NOVASEQ-X POLY-G
# ============================================================
hr("CELL 4 — cDNA quality in R1 (poly-G and complexity)")

if CDNA_START is None:
    print("  cDNA start not established in CELL 3; skipping.")
    qual = pd.DataFrame()
else:
    print(f"  profiling R1[{CDNA_START}:] in {BLOCK}bp blocks\n")
    qual_rows = []
    block_rows = []
    for L, fs in r1_files.items():
        start = ((tso_offset_by_lib.get(L) or TSO_EXPECT_OFFSET) + len(TSO))
        n = 0
        n_polyg = 0
        n_lowcomp = 0
        ent_sum = 0.0
        comp = Counter()
        blocks = defaultdict(lambda: [0, 0.0])   # block idx -> [n, entropy sum]
        usable_len = Counter()
        per_file = max(1, SAMPLE_READS // len(fs))
        for p in fs:
            for seq in iter_fastq_seqs(p, per_file):
                cdna = seq[start:]
                if not cdna:
                    continue
                n += 1
                usable_len[len(cdna)] += 1
                comp.update(b for b in cdna if b in 'ACGTN')
                if max_homopolymer(cdna, 'G') >= POLYG_RUN:
                    n_polyg += 1
                e = shannon(cdna)
                ent_sum += e
                if e < LOWCOMP_ENTROPY:
                    n_lowcomp += 1
                for bi in range(0, len(cdna) // BLOCK):
                    seg = cdna[bi * BLOCK:(bi + 1) * BLOCK]
                    blocks[bi][0] += 1
                    blocks[bi][1] += shannon(seg)
        if not n:
            print(f"  {L}: no cDNA bases past offset {start}; skipped")
            continue
        tot_b = sum(v for k, v in comp.items() if k in 'ACGTN') or 1
        qual_rows.append({
            'lib': L, 'group': 'INF' if L in INFECTED else 'neg',
            'cdna_start': start,
            'cdna_len_mode': usable_len.most_common(1)[0][0],
            'n_reads': n,
            'polyg_pct': round(100.0 * n_polyg / n, 2),
            'lowcomplexity_pct': round(100.0 * n_lowcomp / n, 2),
            'mean_entropy': round(ent_sum / n, 3),
            'pct_G': round(100.0 * comp.get('G', 0) / tot_b, 2),
            'pct_N': round(100.0 * comp.get('N', 0) / tot_b, 2),
        })
        for bi, (bn, be) in sorted(blocks.items()):
            block_rows.append({'lib': L, 'block': bi,
                               'block_start': start + bi * BLOCK,
                               'n': bn, 'mean_entropy': round(be / bn, 3)})
        print(f"  {L}: {n:,} reads, cDNA {usable_len.most_common(1)[0][0]}bp, "
              f"poly-G {100.0 * n_polyg / n:5.2f}%, "
              f"low-complexity {100.0 * n_lowcomp / n:5.2f}%, "
              f"entropy {ent_sum / n:.3f}, G {100.0 * comp.get('G', 0) / tot_b:.1f}%")

    qual = pd.DataFrame(qual_rows)
    if not qual.empty:
        qual.to_csv(os.path.join(OUT_DIR, 'r1_cdna_quality.csv'), index=False)
        pd.DataFrame(block_rows).to_csv(
            os.path.join(OUT_DIR, 'r1_cdna_entropy_by_block.csv'), index=False)
        print("\n  saved: r1_cdna_quality.csv, r1_cdna_entropy_by_block.csv")
        print("\n  entropy by block (2.0 = uniform, <1.2 = low complexity):")
        piv = (pd.DataFrame(block_rows)
               .pivot(index='lib', columns='block', values='mean_entropy'))
        print(piv.to_string())
        print("\n  A clean left-to-right decay toward 0 means the tail is "
              "poly-G and\n  only the early blocks are searchable. A flat "
              "profile means the whole\n  cDNA stretch is real sequence.")


# %% ============================================================
# CELL 5 — CB SANITY AT OFFSET 0, AND THE VERDICT
# ============================================================
hr("CELL 5 — CB whitelist sanity + verdict")

wl = None
if LOAD_WHITELIST and os.path.exists(WHITELIST):
    try:
        with open(WHITELIST) as fh:
            wl = set(line.strip() for line in fh if line.strip())
        print(f"  whitelist loaded: {len(wl):,} barcodes")
    except Exception as e:                       # noqa: BLE001
        print(f"  WARN  could not read whitelist: {e}")
elif LOAD_WHITELIST:
    print(f"  WARN  whitelist not found at {WHITELIST}; CB check skipped")

if wl:
    print("\n  CB match rate at R1[0:16], direct and reverse-complement.")
    print("  (shiv_merge_barcode_diagnostic.py found CellRanger stores the "
          "revcomp\n   of what plain-whitelist STARsolo stores on NovaSeq X; "
          "this says which\n   orientation the RAW read is in, which the "
          "flank scan will need.)\n")
    cb_rows = []
    for L, fs in r1_files.items():
        direct = rc = n = 0
        for seq in iter_fastq_seqs(fs[0], WHITELIST_CHECK_N):
            cb = seq[:CB_LEN]
            if 'N' in cb:
                continue
            n += 1
            if cb in wl:
                direct += 1
            if revcomp(cb) in wl:
                rc += 1
        if n:
            cb_rows.append({'lib': L, 'n': n,
                            'direct_pct': round(100.0 * direct / n, 2),
                            'revcomp_pct': round(100.0 * rc / n, 2)})
            print(f"  {L}: direct {100.0 * direct / n:6.2f}%   "
                  f"revcomp {100.0 * rc / n:6.2f}%")
    if cb_rows:
        cb = pd.DataFrame(cb_rows)
        cb.to_csv(os.path.join(OUT_DIR, 'r1_cb_orientation.csv'), index=False)
        print("\n  saved: r1_cb_orientation.csv")
        d, r = cb['direct_pct'].mean(), cb['revcomp_pct'].mean()
        # Judge the RATIO, not an absolute rate. These are RAW barcodes: they
        # carry sequencing error and come from empty droplets as well as cells,
        # so 60-80% exact whitelist match is normal and healthy. An absolute
        # threshold here was a bug in the first version of this script; it
        # called a 2000:1 result "ambiguous". check_sample_whitelist.py sees
        # ~100% only because it reads CellRanger-CORRECTED barcodes from a BAM.
        FLOOR = 20.0        # below this, neither orientation is really matching
        RATIO = 10.0        # one orientation must beat the other this many times
        hi, lo = max(d, r), min(d, r)
        if hi < FLOOR:
            print(f"\n  *** Neither orientation matches the whitelist "
                  f"(direct {d:.2f}%, revcomp {r:.2f}%).\n      Wrong "
                  "whitelist file, wrong CB_LEN, or wrong read. Resolve before "
                  "scanning.")
        elif lo > 0 and hi / lo < RATIO:
            print(f"\n  *** Genuinely ambiguous: direct {d:.2f}%, revcomp "
                  f"{r:.2f}% (ratio {hi / max(lo, 1e-9):.1f}x).\n      Test "
                  "both at join time and report which wins, as\n      "
                  "shiv_merge_barcode_diagnostic.py did.")
        elif d > r:
            print(f"\n  Raw R1 CB is DIRECT, in whitelist orientation "
                  f"(direct {d:.2f}% vs revcomp {r:.2f}%).")
            print("  So there is NO FASTQ-level orientation flip. Note that "
                  "this contradicts\n  the NovaSeq X revcomp hypothesis in "
                  "shiv_merge_barcode_diagnostic.py's\n  header, so whatever "
                  "caused that merge to put 458 recovered reads on\n  only 2 "
                  "called cells, it was not this. Join direct, but still test "
                  "both\n  at join time until that history is explained.")
        else:
            print(f"\n  Raw R1 CB is REVERSE-COMPLEMENT (revcomp {r:.2f}% vs "
                  f"direct {d:.2f}%).\n  Reverse-complement R1[0:16] before "
                  "joining to the object's barcodes.")

print("\n  VERDICT on the R1 arm")
if CDNA_START is None:
    print("    CANNOT DECIDE. TSO offset not established. Resolve CELL 3 "
          "first.")
elif qual.empty:
    print("    CANNOT DECIDE. No cDNA profile produced.")
else:
    # The decision is driven by the INFECTED libraries, and within those by H,
    # which carries 864 of the ~947 SHIV reads in the whole dataset. A poly-G
    # problem confined to a negative library costs the scan nothing, so it must
    # not veto the arm. Per-library disposition is printed either way.
    inf = qual[qual['group'] == 'INF']
    decide = inf if not inf.empty else qual
    worst_pg = decide['polyg_pct'].max()
    mean_ent = decide['mean_entropy'].mean()
    cdna_len = int(qual['cdna_len_mode'].mode().iloc[0])
    print(f"    cDNA available per R1 : {cdna_len} bp (from base {CDNA_START})")
    print(f"    decision basis        : "
          f"{'infected libraries' if not inf.empty else 'all libraries'} "
          f"({', '.join(decide['lib'])})")
    print(f"    worst poly-G rate     : {worst_pg:.2f}%")
    print(f"    mean entropy          : {mean_ent:.3f} bits/base")
    h = qual[qual['lib'] == 'H']
    if not h.empty:
        print(f"    library H specifically : poly-G "
              f"{h['polyg_pct'].iloc[0]:.2f}%, entropy "
              f"{h['mean_entropy'].iloc[0]:.3f}  "
              f"(H holds 864 of ~947 SHIV reads, so H decides in practice)")
    print("\n    per-library disposition:")
    for _, r in qual.iterrows():
        if r['polyg_pct'] > 50 or r['mean_entropy'] < LOWCOMP_ENTROPY:
            d = 'unusable (poly-G / low complexity)'
        elif r['polyg_pct'] > 10:
            d = 'usable with poly-G trim'
        else:
            d = 'usable'
        print(f"      {r['lib']} [{r['group']}] poly-G {r['polyg_pct']:6.2f}%  "
              f"entropy {r['mean_entropy']:.3f}  -> {d}")
    if worst_pg > 50 or mean_ent < LOWCOMP_ENTROPY:
        print("\n    DROP THE R1 ARM. The post-TSO sequence is dominated by "
              "poly-G or\n    low-complexity runs, so it is instrument "
              "artifact, not transcript.\n    Scan R2 only. Record this in the "
              "scan script header so it is not\n    relitigated.")
    elif worst_pg > 10:
        print("\n    SCAN R1 WITH A POLY-G GUARD. Most reads carry real "
              "sequence but a\n    meaningful minority are poly-G tailed. "
              "Trim at the first G-run of\n    >= POLYG_RUN before matching, "
              "and report trimmed length per hit.")
    else:
        print("\n    SCAN R1. The post-TSO sequence is real transcript. This "
              "adds\n    ~{0}bp of searchable sequence per read pair on top of "
              "R2, and any\n    cassette hit here carries the cell barcode in "
              "the SAME molecule\n    {1} bases upstream, which needs no "
              "pairing assumption and no\n    barcode-orientation "
              "decision.".format(cdna_len, CDNA_START))

hr("R1 STRUCTURE CHECK COMPLETE — read-only. Nothing modified.")
