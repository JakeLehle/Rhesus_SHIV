#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Rhesus SHIV — Which Reads Did the New Reference Actually Change (READ-ONLY)
===========================================================================
Compares the READ POOLS, not the read counts, across the four references
aligned by shiv_starsolo_bc_recovery.sh.

WHY COUNTS CANNOT ANSWER THIS
-----------------------------
The 2026-10-27 run returned 2,188 / 2,192 / 2,190 / 2,192 aligned reads for
V0_nocassette / V2allN / V3r1 / V3r2. A net difference of 4 is equally
consistent with two very different worlds:

    identical pools, nothing changed
    forty reads lost and forty different ones gained

Those have opposite implications. The second would mean the cassette reference
is placing reads the old one could not, and those reads are the ones that
needed a cassette to be placed at all. Only the read IDs separate the two.

WHAT THIS LOOKS FOR
-------------------
  GAINED  aligned to a cassette reference, NOT to V0. This is the prize. A read
          here could not be placed against a reference lacking the insert, which
          is the signature of a read carrying one. CELL 4 checks each one for
          cassette structure and reports where it landed.
  LOST    aligned to V0, not to the cassette reference. Expected near the
          insertion point: a read of plain vpr sequence that spans 5862 in the
          UNBARCODED frame has to cross 34 inserted bases in the barcoded one.
  SHARED  the same read placed by both. If this is ~100%, the reference surgery
          changed nothing about which reads are recoverable, and the negative
          result from the count comparison is confirmed at read resolution.

Uses samtools through subprocess rather than pysam, because pysam lives in
sc_pre_R on this system and not in sc_pre. The BAMs hold hundreds of reads, so
parsing them whole in Python costs nothing.

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
from collections import defaultdict, Counter

import pandas as pd

WORKING = '/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV'
SOLO    = os.path.join(WORKING, 'shiv_starsolo_bc', 'solo')
OUT_DIR = os.path.join(WORKING, 'annotation_output', 'bc_reference_pool_diff')

BASELINE = 'V0_nocassette'
REFS     = ['V0_nocassette', 'V2allN', 'V3r1', 'V3r2']

LIB_NAME = {
    'A': 'SHIV_Pre_PBMC',             'B': 'SHIV_Pre_LN',
    'C': 'SHIV_Pre_PBMC_IndianMixed', 'D': 'SHIV_21DPI_Pre_LN_Mixed',
    'E': 'SHIV_21DPI_PBMC',           'F': 'SHIV_21DPI_NonInf_LN_Mixed',
    'G': 'SHIV_Necropsy_PBMC',        'H': 'SHIV_Necropsy_LN',
}
LIBS     = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
INFECTED = {'D', 'E', 'F', 'H'}

# Cassette window, 0-based half-open. In a cassette reference these 34 bases
# ARE the insert; in V0 the same coordinates are real vpr sequence, which is
# why the two are never compared on this window directly.
LOCUS       = 5862
BC_LEN      = 34
NEAR_LOCUS  = 200        # bp either side counted as "at the barcode locus"

# Cassette and flank motifs, for characterising a gained read's own sequence.
MLU_LEFT   = 'ACGCGC'
MLU_RIGHT  = 'GCGCGT'
FLANK_5    = 'GGACTAGCATAA'      # cassette-proximal 12bp of Keele's 5' flank
FLANK_3    = 'ATGGAAGAAAGA'      # cassette-proximal 12bp of his 3' flank

os.makedirs(OUT_DIR, exist_ok=True)
pd.set_option('display.width', 200)
pd.set_option('display.max_columns', 60)

_COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
CIGAR_RE = re.compile(r'(\d+)([MIDNSHP=X])')


def revcomp(s):
    return s.translate(_COMP)[::-1]


def hr(title):
    print("\n" + "=" * 92)
    print(title)
    print("=" * 92)


SAMTOOLS = shutil.which('samtools')
if not SAMTOOLS:
    print("ERROR: samtools not on PATH. Activate sc_pre and")
    print("       export PATH=\"${CONDA_PREFIX}/bin:${PATH}\" "
          "(the system samtools is broken).")
    raise SystemExit(1)


def bam_path(ref, lib):
    return os.path.join(SOLO, ref,
                        f'{lib}_{LIB_NAME[lib]}_Aligned.sortedByCoord.out.bam')


def read_bam(path):
    """Return {read_id: (pos0, cigar, seq, nM, flag)} for mapped records.
    Keeps the FIRST record per read id; a read with a secondary alignment is
    counted once, which is what a pool comparison wants."""
    out = {}
    if not os.path.exists(path):
        return out
    try:
        p = subprocess.run([SAMTOOLS, 'view', '-F', '0x904', path],
                           capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:      # noqa: BLE001
        print(f"    WARN  samtools failed on {os.path.basename(path)}: {e}")
        return out
    for line in p.stdout.splitlines():
        f = line.split('\t')
        if len(f) < 11:
            continue
        qname = f[0]
        if qname in out:
            continue
        nm = None
        for t in f[11:]:
            if t.startswith('nM:i:') or t.startswith('NM:i:'):
                nm = int(t.split(':')[2])
                break
        out[qname] = (int(f[3]) - 1, f[5], f[9], nm, int(f[1]))
    return out


def ref_end(pos0, cigar):
    r = pos0
    for n, op in CIGAR_RE.findall(cigar):
        if op in 'MDN=X':
            r += int(n)
    return r


def spans_window(pos0, cigar, s, e):
    """One aligned block containing [s,e). A gap over the window does not count."""
    r = pos0
    for n, op in CIGAR_RE.findall(cigar):
        n = int(n)
        if op in 'M=X':
            if r <= s and r + n >= e:
                return True
            r += n
        elif op in 'DN':
            r += n
    return False


def cassette_evidence(seq):
    """What does the READ's own sequence say, independent of where it aligned."""
    hits = []
    for name, motif in (('MluI_L', MLU_LEFT), ('MluI_R', MLU_RIGHT),
                        ('flank5', FLANK_5), ('flank3', FLANK_3)):
        if motif in seq or revcomp(motif) in seq:
            hits.append(name)
    return ','.join(hits) if hits else ''


# %% ============================================================
# CELL 2 — LOAD THE POOLS
# ============================================================
hr("CELL 2 — load read pools")

pools = defaultdict(dict)          # pools[ref][lib] = {qname: record}
for ref in REFS:
    for lib in LIBS:
        p = bam_path(ref, lib)
        if not os.path.exists(p):
            print(f"  WARN  missing {ref}/{lib}: {p}")
            continue
        pools[ref][lib] = read_bam(p)
    got = sum(len(v) for v in pools[ref].values())
    print(f"  {ref:16s} {got:6,} reads across {len(pools[ref])} libraries")

if BASELINE not in pools or not pools[BASELINE]:
    print(f"\n  ERROR: baseline {BASELINE} not loaded; nothing to compare against")
    raise SystemExit(1)


# %% ============================================================
# CELL 3 — SHARED, GAINED, LOST
# ============================================================
hr(f"CELL 3 — pool overlap against {BASELINE}")

rows = []
gained_records = []
for ref in REFS:
    if ref == BASELINE:
        continue
    for lib in LIBS:
        base = set(pools[BASELINE].get(lib, {}))
        new = set(pools[ref].get(lib, {}))
        if not base and not new:
            continue
        shared = base & new
        gained = new - base
        lost = base - new
        rows.append({
            'reference': ref, 'lib': lib,
            'group': 'INF' if lib in INFECTED else 'neg',
            'baseline_n': len(base), 'new_n': len(new),
            'shared': len(shared), 'gained': len(gained), 'lost': len(lost),
            'pct_shared': round(100.0 * len(shared) / len(new), 2) if new else 0,
        })
        for q in gained:
            pos0, cigar, seq, nm, flag = pools[ref][lib][q]
            gained_records.append({
                'reference': ref, 'lib': lib,
                'group': 'INF' if lib in INFECTED else 'neg',
                'read_id': q, 'pos0': pos0, 'cigar': cigar, 'nM': nm,
                'ref_end': ref_end(pos0, cigar),
                'at_locus': abs(pos0 - LOCUS) <= NEAR_LOCUS,
                'spans_cassette': spans_window(pos0, cigar, LOCUS,
                                               LOCUS + BC_LEN),
                'read_motifs': cassette_evidence(seq),
                'seq': seq,
            })

diff = pd.DataFrame(rows)
if diff.empty:
    print("  nothing to compare")
    raise SystemExit(0)
diff.to_csv(os.path.join(OUT_DIR, 'pool_diff_per_library.csv'), index=False)
print(diff.to_string(index=False))
print("\n  saved: pool_diff_per_library.csv")

tot = (diff.groupby('reference')
       .agg(baseline_n=('baseline_n', 'sum'), new_n=('new_n', 'sum'),
            shared=('shared', 'sum'), gained=('gained', 'sum'),
            lost=('lost', 'sum')).reset_index())
tot['pct_shared'] = (100.0 * tot['shared'] / tot['new_n']).round(2)
print("\n  totals:")
print(tot.to_string(index=False))
tot.to_csv(os.path.join(OUT_DIR, 'pool_diff_totals.csv'), index=False)

print("\n  pct_shared near 100 means the reference surgery changed WHICH reads")
print("  are recoverable hardly at all, and the flat count comparison is")
print("  confirmed at read resolution rather than only in aggregate.")


# %% ============================================================
# CELL 4 — WHAT ARE THE GAINED READS
# ============================================================
hr("CELL 4 — reads the old reference could not place")

if not gained_records:
    print("  No gained reads at all. Every read the cassette references placed")
    print("  was already placed by the cassette-less one. The insert is not")
    print("  what was keeping any read out, which closes the reference question")
    print("  at read resolution.")
else:
    g = pd.DataFrame(gained_records)
    g.drop(columns=['seq']).to_csv(
        os.path.join(OUT_DIR, 'pool_diff_gained_reads.csv'), index=False)
    print(f"  {len(g)} gained read(s) across all references and libraries")
    print(f"    at the barcode locus (+/-{NEAR_LOCUS}bp): "
          f"{int(g['at_locus'].sum())}")
    print(f"    spanning the cassette window               : "
          f"{int(g['spans_cassette'].sum())}")
    print(f"    carrying a cassette or flank motif in the read: "
          f"{int((g['read_motifs'] != '').sum())}")
    print(f"    in INFECTED libraries                      : "
          f"{int((g['group'] == 'INF').sum())}")
    print(f"    in NEGATIVE libraries                      : "
          f"{int((g['group'] == 'neg').sum())}")

    interesting = g[g['at_locus'] | g['spans_cassette'] |
                    (g['read_motifs'] != '')]
    if len(interesting):
        print("\n  gained reads at or near the locus, or carrying a motif:")
        print(interesting.drop(columns=['seq'])
              .sort_values(['reference', 'lib']).to_string(index=False))
    else:
        print("\n  None of the gained reads sit near the barcode locus or carry")
        print("  a cassette motif, so they are incidental alignment differences")
        print("  elsewhere on the contig, not barcode recovery.")

    print("\n  saved: pool_diff_gained_reads.csv")


# %% ============================================================
# CELL 5 — VERDICT
# ============================================================
hr("CELL 5 — verdict")

best = tot.sort_values('gained', ascending=False).iloc[0]
print(f"  most reads gained over {BASELINE}: {best['reference']} "
      f"({int(best['gained'])})")
if gained_records:
    gi = pd.DataFrame(gained_records)
    useful = gi[(gi['group'] == 'INF') &
                (gi['at_locus'] | gi['spans_cassette'] |
                 (gi['read_motifs'] != ''))]
    print(f"  of those, relevant to the barcode in infected libraries: "
          f"{len(useful)}")
    if len(useful) == 0:
        print("\n  The cassette reference is not recovering barcode-bearing")
        print("  reads that the old one missed. Combined with the flat count")
        print("  comparison, the reference was not the limiting factor.")
    else:
        print("\n  These are reads that could not be placed without the")
        print("  cassette. Pull their sequences from pool_diff_gained_reads.csv")
        print("  and run them through the same cassette parse the flank scan")
        print("  uses, then union on read id with flank_scan_hits.csv.")
else:
    print("  no gained reads; see CELL 4")

print("\n  REMEMBER what this comparison can and cannot see. STARsolo aligns R2")
print("  only; R1 is consumed for CB and UMI and never treated as cDNA. The")
print("  flank scan found 11 of its 16 hits in R1. So every number above is")
print("  computed on roughly half the available sequence, and no choice of")
print("  reference changes that. Aligning R1's cDNA is a separate experiment.")

hr("POOL DIFF COMPLETE — read-only.")
