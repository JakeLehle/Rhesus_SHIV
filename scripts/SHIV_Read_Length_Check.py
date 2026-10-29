#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Rhesus SHIV — Read Length Check (READ-ONLY)
===========================================
Measures the actual R1 / R2 read lengths in the raw GEX FASTQs, per library,
before the flank-anchored barcode scan is written. Nothing is modified.

WHY THIS EXISTS
---------------
The expected yield of ANY barcode-recovery scan is set by read length, because a
read can only be called if it covers enough of the cassette. At 90bp the
flank-anchored scan is predicted to return ~8 reads across the whole dataset; at
150bp it is ~14. That difference decides whether the scan is worth the compute,
so it gets measured rather than assumed.

It also catches two things that would silently corrupt the scan:
  1. R1 shorter than CB+UMI (16+12=28). Then the cell barcode cannot be read at
     full length and every join degrades.
  2. Mixed read lengths within a library, which happens after trimming or after
     the repair_corrupt_fastqs.sh truncation pass. A single "first read length"
     report (bd_abseq_scan.sh Stage 0) would miss this, so this script samples a
     window of records and reports the full distribution.

WHAT IT DOES NOT DO
-------------------
Does not decompress whole files. It reads the first SAMPLE_RECORDS records of
each file and stops, so it is seconds per library, not hours. Read length is set
by the sequencing run and is constant within a file except for trimming, so a
head sample is sufficient. If you want certainty on a specific file, raise
SAMPLE_RECORDS for that one; the cost is linear.

Conventions: Spyder cells (# %%), config in Cell 1, warn-and-continue, read-only.

Author: Jake Lehle
Date: October 2026
"""

# %% ============================================================
# CELL 1 — CONFIG AND DIALS
# ============================================================
import os
import glob
import gzip
from collections import Counter, defaultdict

import pandas as pd

# --- paths -----------------------------------------------------------------
RAW_BASE   = '/master/zwallis/WORKING/SC/10X_RAW'
OUT_DIR    = ('/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV/'
              'annotation_output/read_length_check')

# --- libraries -------------------------------------------------------------
LIBS = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
INFECTED = {'D', 'E', 'F', 'H'}

# Library subdirectory pattern, matching shiv_vpr_barcode_feasibility.py CELL 5
# and SHIV_vpr_barcode_gt_search.py. PROTEIN / VDJ_RH are reported too because a
# length mismatch there explains ADT and TCR behaviour, but only GEX feeds the
# barcode scan.
LIB_DIR_PATTERN = {
    'GEX':     'GEX_LIBRARY_{L}',
    'PROTEIN': 'PROTEIN_LIBRARY_{L}',
    'VDJ_RH':  'VDJ_RH_LIBRARY_{L}',
}
SCAN_LIB_TYPES = ['GEX', 'PROTEIN', 'VDJ_RH']

# Which reads to measure. I1/I2 are included when present; they are indices and
# are not used by the scan, but a missing index file is worth seeing.
READ_TAGS = ['R1', 'R2', 'I1', 'I2']

# --- sampling --------------------------------------------------------------
SAMPLE_RECORDS = 20_000      # records per file; raise for certainty on one file
MAX_FILES_PER_TAG = None     # None = every lane; int = cap for a fast look

# --- expected geometry (GEM-X 5' v3, confirmed by check_sample_whitelist.py) -
CB_LEN  = 16
UMI_LEN = 12

# --- barcode geometry, for the yield prediction in CELL 4 -------------------
GENOME_LEN     = 10262       # SHIVAD8EO contig, MN816822.1
CASSETTE_LEN   = 34          # 12bp const flank + 10bp core + 12bp const flank
CASSETTE_FLANK = 12
MIN_VIRUS_FLANK = 20         # bp of virus sequence required to anchor
MIN_CORE_BASES  = 6          # of the 10bp core, how many must be visible

# SHIV-mapped read pool per library, from the baselines recorded in
# shiv_starsolo_recovery.sh (MM2_BASELINE) and shiv_denovo_assembly.sh.
# F is not recorded anywhere in the repo; left out, so the prediction is a FLOOR.
SHIV_READ_POOL = {'H': 864, 'E': 48, 'D': 35}

os.makedirs(OUT_DIR, exist_ok=True)
pd.set_option('display.width', 170)
pd.set_option('display.max_columns', 50)


def hr(title):
    print("\n" + "=" * 84)
    print(title)
    print("=" * 84)


def sample_lengths(path, n_records):
    """Return (Counter of seq lengths, n_sampled, error_or_None).

    Reads at most n_records FASTQ records. Warn-and-continue: a corrupt or
    truncated stream returns whatever was read plus the error string, rather
    than raising, because several of these files have been repaired and a hard
    stop here would hide the rest of the table.
    """
    lens = Counter()
    n = 0
    err = None
    try:
        with gzip.open(path, 'rt') as fh:
            for i, line in enumerate(fh):
                if i % 4 == 1:                 # sequence line
                    lens[len(line.rstrip('\n'))] += 1
                    n += 1
                    if n >= n_records:
                        break
    except Exception as e:                      # noqa: BLE001 - warn and continue
        err = f"{type(e).__name__}: {e}"
    return lens, n, err


# %% ============================================================
# CELL 2 — DISCOVER FILES
# ============================================================
hr("CELL 2 — FASTQ discovery")

found = []          # one row per file
missing = []

for L in LIBS:
    for lib_type in SCAN_LIB_TYPES:
        d = os.path.join(RAW_BASE, LIB_DIR_PATTERN[lib_type].format(L=L))
        if not os.path.isdir(d):
            missing.append((L, lib_type, d, 'directory not found'))
            continue
        for tag in READ_TAGS:
            fs = sorted(glob.glob(os.path.join(d, f'*_{tag}_*.fastq.gz')))
            if not fs:
                fs = sorted(glob.glob(os.path.join(d, f'*_{tag}_*.fq.gz')))
            if not fs:
                # R1/R2 absent is a real problem; I1/I2 absent is routine
                if tag in ('R1', 'R2'):
                    missing.append((L, lib_type, d, f'no {tag} files'))
                continue
            if MAX_FILES_PER_TAG:
                fs = fs[:MAX_FILES_PER_TAG]
            for p in fs:
                found.append({'lib': L, 'lib_type': lib_type, 'tag': tag,
                              'path': p,
                              'size_gb': round(os.path.getsize(p) / 1e9, 2)})

print(f"  found   : {len(found)} FASTQ files across {len(LIBS)} libraries")
print(f"  missing : {len(missing)}")
for L, lt, d, why in missing:
    print(f"    WARN  {L} {lt:8s} {why}  ({d})")

if not found:
    print("\n  Nothing to measure. Check RAW_BASE and LIB_DIR_PATTERN.")


# %% ============================================================
# CELL 3 — MEASURE
# ============================================================
hr(f"CELL 3 — length sampling ({SAMPLE_RECORDS:,} records per file)")

rows = []
for i, rec in enumerate(found, 1):
    lens, n, err = sample_lengths(rec['path'], SAMPLE_RECORDS)
    if err:
        print(f"    WARN  {os.path.basename(rec['path'])}: {err}")
    if not lens:
        rows.append({**rec, 'n_sampled': n, 'len_min': None, 'len_max': None,
                     'len_mode': None, 'mode_pct': None, 'n_distinct': 0,
                     'error': err})
        continue
    mode_len, mode_n = lens.most_common(1)[0]
    rows.append({**rec,
                 'n_sampled':  n,
                 'len_min':    min(lens),
                 'len_max':    max(lens),
                 'len_mode':   mode_len,
                 'mode_pct':   round(100.0 * mode_n / n, 2),
                 'n_distinct': len(lens),
                 'error':      err})
    if i % 20 == 0:
        print(f"    ...{i}/{len(found)} files")

df = pd.DataFrame(rows)
if df.empty:
    print("  no measurements; stopping before CELL 4")
else:
    df.to_csv(os.path.join(OUT_DIR, 'read_lengths_per_file.csv'), index=False)
    print(f"\n  saved: read_lengths_per_file.csv  ({len(df)} rows)")

    # --- per library x lib_type x tag rollup ------------------------------
    agg = (df.groupby(['lib_type', 'lib', 'tag'])
             .agg(n_files=('path', 'size'),
                  n_sampled=('n_sampled', 'sum'),
                  len_min=('len_min', 'min'),
                  len_max=('len_max', 'max'),
                  len_mode=('len_mode', lambda s: s.mode().iloc[0]
                            if not s.mode().empty else None),
                  n_distinct=('n_distinct', 'max'))
             .reset_index())
    agg.to_csv(os.path.join(OUT_DIR, 'read_lengths_summary.csv'), index=False)
    print("  saved: read_lengths_summary.csv")

    for lib_type in SCAN_LIB_TYPES:
        sub = agg[agg['lib_type'] == lib_type]
        if sub.empty:
            continue
        print(f"\n  --- {lib_type} ---")
        print(sub.drop(columns=['lib_type']).to_string(index=False))

    # --- flags ------------------------------------------------------------
    print("\n  FLAGS")
    flagged = False

    r1 = df[(df['lib_type'] == 'GEX') & (df['tag'] == 'R1')]
    short_r1 = r1[r1['len_min'].fillna(0) < CB_LEN + UMI_LEN]
    if not short_r1.empty:
        flagged = True
        print(f"    R1 shorter than CB+UMI ({CB_LEN}+{UMI_LEN}="
              f"{CB_LEN + UMI_LEN}) in {len(short_r1)} GEX file(s):")
        for _, r in short_r1.iterrows():
            print(f"      {r['lib']}  min={r['len_min']}  "
                  f"{os.path.basename(r['path'])}")

    mixed = df[(df['n_distinct'] > 1) & (df['mode_pct'].fillna(100) < 99.0)]
    if not mixed.empty:
        flagged = True
        print(f"    mixed read lengths (mode <99% of sampled) in "
              f"{len(mixed)} file(s):")
        for _, r in mixed.iterrows():
            print(f"      {r['lib']} {r['lib_type']:8s} {r['tag']}  "
                  f"{r['len_min']}..{r['len_max']} "
                  f"mode={r['len_mode']} ({r['mode_pct']}%)  "
                  f"{os.path.basename(r['path'])}")

    gex_r2 = df[(df['lib_type'] == 'GEX') & (df['tag'] == 'R2')]
    if not gex_r2.empty and gex_r2['len_mode'].nunique() > 1:
        flagged = True
        print("    GEX R2 length is NOT uniform across libraries:")
        print(gex_r2.groupby('lib')['len_mode']
              .agg(lambda s: sorted(set(s))).to_string())

    if not flagged:
        print("    none. Lengths uniform, R1 long enough for CB+UMI.")


# %% ============================================================
# CELL 4 — WHAT THAT LENGTH IMPLIES FOR THE BARCODE SCAN
# ============================================================
hr("CELL 4 — predicted barcode-scan yield at the measured R2 length")

if df.empty or df[(df['lib_type'] == 'GEX') & (df['tag'] == 'R2')].empty:
    print("  no GEX R2 measurement; cannot predict yield")
else:
    gex_r2 = df[(df['lib_type'] == 'GEX') & (df['tag'] == 'R2')]
    rl = int(gex_r2['len_mode'].mode().iloc[0])
    print(f"  GEX R2 modal length : {rl} bp")

    pool_total = sum(SHIV_READ_POOL.values())
    print(f"  SHIV read pool      : {pool_total} "
          f"({', '.join(f'{k}={v}' for k, v in SHIV_READ_POOL.items())}; "
          f"F not recorded, so this is a FLOOR)")
    print(f"  genome length       : {GENOME_LEN} bp")

    def starts_full(rl_):
        """Read must contain the entire 34bp cassette (the old anchor)."""
        return max(0, rl_ - CASSETTE_LEN + 1)

    def starts_flank(rl_, vflank, mincore):
        """Keele flank anchor: vflank bp of virus + 12bp cassette flank +
        mincore core bases, from either side, minus the double-counted reads
        that span both virus flanks."""
        need = vflank + CASSETTE_FLANK + mincore
        one  = max(0, rl_ - need + 1)
        both = max(0, rl_ - (vflank + CASSETTE_LEN + vflank) + 1)
        return 2 * one - both

    sf = starts_full(rl)
    sk = starts_flank(rl, MIN_VIRUS_FLANK, MIN_CORE_BASES)

    print(f"\n  acceptance windows at {rl} bp, over a {GENOME_LEN} bp genome:")
    print(f"    old cassette anchor (all {CASSETTE_LEN}bp in one read)      : "
          f"{sf:4d} start positions -> {pool_total * sf / GENOME_LEN:5.1f} reads")
    print(f"    flank anchor (>={MIN_VIRUS_FLANK}bp virus, >={MIN_CORE_BASES} "
          f"core bases) : {sk:4d} start positions -> "
          f"{pool_total * sk / GENOME_LEN:5.1f} reads")

    print("\n  per library (flank anchor):")
    for L, pool in sorted(SHIV_READ_POOL.items()):
        print(f"    {L}: pool {pool:4d} -> "
              f"{pool * sk / GENOME_LEN:5.2f} predicted reads")

    print("\n  MIN_CORE_BASES sensitivity:")
    for mc in (3, 4, 6, 8, 10):
        s = starts_flank(rl, MIN_VIRUS_FLANK, mc)
        print(f"    >={mc:2d} core bases : {s:4d} starts -> "
              f"{pool_total * s / GENOME_LEN:5.1f} reads")

    print("\n  depth required for a given yield (flank anchor, "
          f">={MIN_CORE_BASES} core bases):")
    for target in (25, 100, 1000):
        need = target * GENOME_LEN / sk if sk else float('inf')
        print(f"    {target:5d} barcode reads -> ~{need:11,.0f} SHIV reads "
              f"({need / pool_total:6.1f}x current pool)")

    print("\n  NOTE: these are predictions for reads ALIGNED to SHIV. The raw-"
          "FASTQ\n  scan is a superset, because it does not depend on any "
          "alignment decision,\n  so treat the numbers above as the floor "
          "rather than the expectation.")

hr("READ LENGTH CHECK COMPLETE — read-only. Nothing modified.")
