#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Rhesus SHIV — Are the R1 Alignments Virus or Noise (READ-ONLY)
==============================================================
Audits the reads that shiv_starsolo_r1_recovery.sh aligned, to decide whether a
35x apparent gain is signal or an artifact of loose filters.

WHY THIS EXISTS
---------------
The 2026-10-28 smoke test aligned 3,623 reads from 20M R1 reads of library H,
which is 181 per million against 5.2 per million for the R2 runs. Taken at face
value that is a 35x gain and ~61 cassette-spanning reads at full depth.

Three things in the same output argue it may not be real:

  STRAND TIE       Solo-counted Forward 834, Reverse 804 = 1.04x. On R2 the
                   same comparison was 337 vs 7 = 48x. The collapsed GTF puts
                   one gene on the plus strand and viral mRNA is plus strand,
                   so R1 reading sense out of the TSO should win decisively on
                   Forward. A tie is what alignments with NO strand preference
                   look like.
  MAPQ             47% of aligned reads sit below MAPQ30.
  BREADTH          97.9% of a 10kb contig covered by 3,623 reads, very even.

And the filters are loose for a 109bp read: --outFilterMatchNmin 25 with
--outFilterMismatchNoverLmax 0.2 means 25 matched bases and up to 21 mismatches
against a SHIV-ONLY index, where no host sequence competes for the read.

WHAT SEPARATES THE TWO
----------------------
                        real virus              spurious
  nM distribution       piled at 0-3            piled near the 21 ceiling
  matched length        most of the 109bp       near the 25bp floor
  strand                biased                  ~50/50
  MAPQ                  mostly high             much of it low
  position              uneven, structured      flat across the contig
  negatives             near-empty              as busy as the infected

The negative library is the decisive one and needs its own alignment run; this
script reads whatever BAMs exist and compares them if more than one is present.
Everything else here is computable from the BAM already on disk.

Read-only. Spyder cells (# %%), config in Cell 1, warn-and-continue.

Author: Jake Lehle
Date: October 2026
"""

# %% ============================================================
# CELL 1 — CONFIG AND DIALS
# ============================================================
import os
import re
import shutil
import subprocess
from collections import Counter

import pandas as pd

WORKING = '/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV'
SOLO    = os.path.join(WORKING, 'shiv_starsolo_r1', 'solo')
OUT_DIR = os.path.join(WORKING, 'annotation_output', 'r1_alignment_audit')

STRAND  = 'Forward'            # which run's BAMs to audit
LIBS    = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
INFECTED = {'D', 'E', 'F', 'H'}
LIB_NAME = {
    'A': 'SHIV_Pre_PBMC',             'B': 'SHIV_Pre_LN',
    'C': 'SHIV_Pre_PBMC_IndianMixed', 'D': 'SHIV_21DPI_Pre_LN_Mixed',
    'E': 'SHIV_21DPI_PBMC',           'F': 'SHIV_21DPI_NonInf_LN_Mixed',
    'G': 'SHIV_Necropsy_PBMC',        'H': 'SHIV_Necropsy_LN',
}

# Alignment settings these BAMs were produced under, for scoring against.
READ_LEN_CDNA   = 109          # 150 - 41 (CB+UMI+TSO)
MATCH_MIN_ABS   = 25           # the floor the AUDITED BAM was produced under
MISMATCH_MAX    = int(READ_LEN_CDNA * 0.2)      # 21
# The floor the runner now uses. Applied here in silico so the effect can be
# previewed on an existing BAM without waiting for a realignment, and so the
# threshold itself can be sanity-checked rather than assumed.
MIN_MATCHED_CLEAN = 70
LOCUS, BC_LEN   = 5862, 34
CONTIG_LEN      = 10296
POS_BINS        = 40

# Cassette / flank motifs, to check a read's own sequence independent of where
# it aligned.
# 12-mers ONLY. The first version also searched the 6bp MluI half-sites
# ACGCGC and GCGCGT, which occur by chance in about 5% of 109bp reads; the
# 2026-10-28 run returned 2.7%, i.e. BELOW chance, so that column carried no
# information at all. A specific 12-mer is expected in ~6e-6 of reads, so any
# hit here is meaningful.
MOTIFS = {
    'flank5': 'GGACTAGCATAA',     # cassette-proximal 12bp of Keele's 5' flank
    'flank3': 'ATGGAAGAAAGA',     # cassette-proximal 12bp of his 3' flank
    'cass_L': 'ACGCGCAGCGTA',     # register 1 left constant flank
    'cass_R': 'CAGCTAGCGCGT',     # register 1 right constant flank
    'cass_L2': 'ACGCGCTAGCTG',    # register 2 left
    'cass_R2': 'TACGCTGCGCGT',    # register 2 right
}

os.makedirs(OUT_DIR, exist_ok=True)
pd.set_option('display.width', 200)
pd.set_option('display.max_columns', 60)

_COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
CIGAR_RE = re.compile(r'(\d+)([MIDNSHP=X])')
revcomp = lambda s: s.translate(_COMP)[::-1]        # noqa: E731


def hr(t):
    print("\n" + "=" * 92)
    print(t)
    print("=" * 92)


SAMTOOLS = shutil.which('samtools')
if not SAMTOOLS:
    print("ERROR: samtools not on PATH. conda activate sc_pre and")
    print("       export PATH=\"${CONDA_PREFIX}/bin:${PATH}\"")
    raise SystemExit(1)


def matched_bases(cigar):
    return sum(int(n) for n, op in CIGAR_RE.findall(cigar) if op in 'M=X')


def spans_cassette(pos0, cigar):
    r = pos0
    for n, op in CIGAR_RE.findall(cigar):
        n = int(n)
        if op in 'M=X':
            if r <= LOCUS and r + n >= LOCUS + BC_LEN:
                return True
            r += n
        elif op in 'DN':
            r += n
    return False


# %% ============================================================
# CELL 2 — LOAD
# ============================================================
hr(f"CELL 2 — load R1 alignments (soloStrand {STRAND})")

rows = []
for lib in LIBS:
    bam = os.path.join(SOLO, STRAND,
                       f'{lib}_{LIB_NAME[lib]}_Aligned.sortedByCoord.out.bam')
    if not os.path.exists(bam):
        print(f"  (no BAM for {lib})")
        continue
    try:
        p = subprocess.run([SAMTOOLS, 'view', '-F', '0x904', bam],
                           capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:      # noqa: BLE001
        print(f"  WARN  {lib}: samtools failed: {e}")
        continue
    n = 0
    for line in p.stdout.splitlines():
        f = line.split('\t')
        if len(f) < 11:
            continue
        n += 1
        flag, pos0, mapq, cigar, seq = int(f[1]), int(f[3]) - 1, int(f[4]), \
            f[5], f[9]
        nm = None
        for t in f[11:]:
            if t.startswith(('nM:i:', 'NM:i:')):
                nm = int(t.split(':')[2])
                break
        hits = [k for k, m in MOTIFS.items()
                if m in seq or revcomp(m) in seq]
        rows.append({
            'lib': lib, 'group': 'INF' if lib in INFECTED else 'neg',
            'pos0': pos0, 'mapq': mapq, 'cigar': cigar, 'nM': nm,
            'matched': matched_bases(cigar),
            'reverse': bool(flag & 16),
            'spans': spans_cassette(pos0, cigar),
            'motifs': ','.join(hits),
        })
    print(f"  {lib} [{'INF' if lib in INFECTED else 'neg'}]: {n} aligned")

df = pd.DataFrame(rows)
if df.empty:
    print("\n  No alignments loaded. Check SOLO and STRAND.")
    raise SystemExit(1)
df.to_csv(os.path.join(OUT_DIR, 'r1_alignments.csv'), index=False)
print(f"\n  saved: r1_alignments.csv ({len(df)} rows)")


# %% ============================================================
# CELL 3 — THE FIVE SIGNATURES
# ============================================================
hr("CELL 3 — signal or noise, per library")

summ = []
for lib, g in df.groupby('lib'):
    nm = g['nM'].dropna()
    summ.append({
        'lib': lib, 'group': g['group'].iloc[0], 'n': len(g),
        'nM_median': nm.median() if len(nm) else None,
        'pct_nM_le3': round(100.0 * (nm <= 3).mean(), 1) if len(nm) else None,
        'pct_nM_ge15': round(100.0 * (nm >= 15).mean(), 1) if len(nm) else None,
        'matched_median': g['matched'].median(),
        'pct_matched_near_floor': round(
            100.0 * (g['matched'] <= MATCH_MIN_ABS + 5).mean(), 1),
        'pct_reverse': round(100.0 * g['reverse'].mean(), 1),
        'pct_mapq_ge30': round(100.0 * (g['mapq'] >= 30).mean(), 1),
        'spanning': int(g['spans'].sum()),
        'with_motif': int((g['motifs'] != '').sum()),
    })
s = pd.DataFrame(summ)
print(s.to_string(index=False))
s.to_csv(os.path.join(OUT_DIR, 'r1_audit_by_library.csv'), index=False)

print(f"\n  reference points: mismatch ceiling {MISMATCH_MAX}, matched floor "
      f"{MATCH_MIN_ABS}, full cDNA {READ_LEN_CDNA}bp")
print("  real virus  -> nM median low, matched near 109, strand biased, MAPQ high")
print(f"  spurious    -> nM piled near {MISMATCH_MAX}, matched near "
      f"{MATCH_MIN_ABS}, strand ~50/50, MAPQ low")


# %% ============================================================
# CELL 3b — THE CLEAN SUBSET, AND THE N-SINK TEST
# ============================================================
hr(f"CELL 3b — in-silico floor at matched >= {MIN_MATCHED_CLEAN}")

clean = df[df['matched'] >= MIN_MATCHED_CLEAN]
print(f"  {len(clean)} of {len(df)} reads survive "
      f"({100.0 * len(clean) / len(df):.0f}%)")

if clean.empty:
    print("  Nothing survives; the whole population is the noise mode.")
else:
    for label, sub in (('ALL', df), (f'matched>={MIN_MATCHED_CLEAN}', clean)):
        nm = sub['nM'].dropna()
        print(f"\n  {label}:")
        print(f"    n                {len(sub)}")
        print(f"    nM median        "
              f"{nm.median() if len(nm) else 'NA'}")
        print(f"    matched median   {sub['matched'].median():.0f}")
        print(f"    pct reverse      {100.0 * sub['reverse'].mean():.1f}%"
              f"   (50% = no strand preference)")
        print(f"    pct MAPQ>=30     {100.0 * (sub['mapq'] >= 30).mean():.1f}%")
        print(f"    cassette-spanning{int(sub['spans'].sum()):5d}")
        print(f"    carrying a 12-mer{int((sub['motifs'] != '').sum()):5d}")

    print("\n  If the surviving reads are real viral cDNA their reverse")
    print("  fraction should swing well away from 50%. If it stays at a coin")
    print("  flip after a length floor, the problem is not the filters.")

# --- do the spanning reads actually contain a cassette? ---------------------
print("\n  SPANNING vs MOTIF (the N-sink test)")
sp = df[df['spans']]
print(f"    reads spanning the cassette window : {len(sp)}")
if len(sp):
    print(f"    of those, carrying a 12-mer        : "
          f"{int((sp['motifs'] != '').sum())}")
    print(sp[['lib', 'pos0', 'cigar', 'nM', 'matched', 'reverse',
              'mapq', 'motifs']].to_string(index=False))
print("\n  A read that truly spans 5862..5896 covers BOTH 12bp constant")
print("  cassette flanks and should carry them. If it spans but carries")
print("  nothing, the 22 free N bases let a weak alignment sit there: the")
print("  reference is attracting junk rather than recovering barcodes.")
print("  (Hypermutation can break a 12-mer, so treat this as strong")
print("   evidence rather than proof.)")

# --- fine-grained position around the locus --------------------------------
lo, hi, binw = LOCUS - 300, LOCUS + 300, 20
print(f"\n  fine position {lo}..{hi}, {binw}bp bins "
      f"(cassette {LOCUS}..{LOCUS + BC_LEN}):")
for label, sub in (('ALL', df), (f'>={MIN_MATCHED_CLEAN}', clean)):
    near = sub[(sub['pos0'] >= lo) & (sub['pos0'] < hi)]
    hcnt = Counter(((near['pos0'] - lo) // binw).astype(int))
    mx = max(hcnt.values()) if hcnt else 1
    print(f"    --- {label} ({len(near)} reads in window) ---")
    for b in range((hi - lo) // binw):
        start = lo + b * binw
        mark = '  <-- CASSETTE' if start <= LOCUS < start + binw else ''
        print(f"      {start:6d} | "
              f"{'#' * min(40, hcnt.get(b, 0) * 40 // mx)} "
              f"{hcnt.get(b, 0)}{mark}")


# %% ============================================================
# CELL 4 — DISTRIBUTIONS
# ============================================================
hr("CELL 4 — nM, matched length, position")

nm = df['nM'].dropna()
if len(nm):
    print("  nM histogram (all libraries):")
    h = Counter(nm.astype(int))
    for v in sorted(h):
        print(f"    nM {v:3d} | {'#' * min(60, h[v] * 60 // max(h.values()))} "
              f"{h[v]}")

print("\n  matched-bases histogram, 10bp bins:")
hb = Counter((df['matched'] // 10) * 10)
for v in sorted(hb):
    print(f"    {v:3d}-{v+9:3d} | "
          f"{'#' * min(60, hb[v] * 60 // max(hb.values()))} {hb[v]}")

print(f"\n  position across the {CONTIG_LEN}bp contig, {POS_BINS} bins:")
w = CONTIG_LEN // POS_BINS
hp = Counter((df['pos0'] // w).clip(0, POS_BINS - 1))
mx = max(hp.values()) if hp else 1
for b in range(POS_BINS):
    star = '  <-- cassette' if b * w <= LOCUS < (b + 1) * w else ''
    print(f"    {b*w:6d} | {'#' * min(50, hp.get(b, 0) * 50 // mx)} "
          f"{hp.get(b, 0)}{star}")
print("\n  Flat across every bin is what a SHIV-only index does to random reads.")
print("  Real viral coverage is uneven and follows transcript structure.")


# %% ============================================================
# CELL 5 — THE FLOOR
# ============================================================
hr("CELL 5 — infected against negative")

if set(df['group']) < {'INF', 'neg'}:
    print("  Only one group present. The decisive comparison needs a negative")
    print("  library aligned the same way:")
    print("    LIBS='A H' MODE=smoke sbatch shiv_starsolo_r1_recovery.sh")
    print("  Until then nothing here distinguishes 35x signal from 35x noise.")
else:
    g = (df.groupby('group')
         .agg(reads=('pos0', 'size'), spanning=('spans', 'sum'),
              with_motif=('motifs', lambda s: (s != '').sum()),
              nM_median=('nM', 'median'),
              pct_mapq30=('mapq', lambda s: round(100.0 * (s >= 30).mean(), 1)),
              pct_reverse=('reverse', lambda s: round(100.0 * s.mean(), 1)))
         .reset_index())
    print(g.to_string(index=False))
    g.to_csv(os.path.join(OUT_DIR, 'r1_audit_by_group.csv'), index=False)
    inf = df[df['group'] == 'INF']
    neg = df[df['group'] == 'neg']
    n_inf = df[df['group'] == 'INF']['lib'].nunique() or 1
    n_neg = df[df['group'] == 'neg']['lib'].nunique() or 1
    ri, rn = len(inf) / n_inf, len(neg) / n_neg
    print(f"\n  per library: infected {ri:.0f} reads, negative {rn:.0f} reads")
    if rn == 0:
        print("  Negatives are empty. The infected signal is real.")
    else:
        print(f"  ratio {ri / rn:.2f}x")
        if ri / rn < 3:
            print("\n  *** The negatives are as busy as the infected libraries.")
            print("  *** These alignments are NOT virus. Raise "
                  "--outFilterMatchNmin")
            print("  *** toward 60-80 and tighten --outFilterMismatchNoverLmax")
            print("  *** to ~0.08, then re-run the smoke test.")
        else:
            print("\n  Infected exceeds the floor. Compare the cassette-spanning")
            print("  and motif-carrying counts too; those are the ones that")
            print("  matter for barcode recovery.")

hr("R1 ALIGNMENT AUDIT COMPLETE — read-only.")
