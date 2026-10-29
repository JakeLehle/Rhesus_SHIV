#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Rhesus SHIV — Barcoded Reference Construction (Read-Only on Inputs)
===================================================================
Builds the SHIV reference variants needed to align reads THROUGH the vpr barcode
cassette, verifies Brandon Keele's reference against ours base by base, and
writes the STAR commands to index them. Reads nothing destructively; every
output is a new file.

WHY THE REFERENCE IS THE PRIORITY
---------------------------------
Everything downstream of alignment in this project inherits a reference that is
missing 34 bases at exactly the locus of interest. A read spanning the cassette
cannot align contiguously to a reference that lacks it, so the ~947-read SHIV
pool is specifically depleted at the barcode window, and every count built on it
inherits that hole. Improving the raw-read string scan (shiv_vpr_flank_scan.py)
works around the problem; this fixes it.

The decisive advantage is hypermutation. shiv_vpr_flank_scan.py requires an
EXACT 20bp flank match, so an APOBEC3-edited provirus (G->A on the plus strand)
breaks the match outright: the 3' probe ATGGAAGAAAGACCTCCAGA carries six G
positions. An aligner just counts mismatches. A read with 20% of its G's edited
carries ~8 extra mismatches in 150bp, which no exact probe survives and any
permissive aligner absorbs. Necropsy lymph node, which produced 14 of the 16
calls, is exactly where edited and defective proviruses concentrate, and the
Binhua catalog is plasma-derived from intact virus, so an edited sequence would
have failed BOTH the flank match and the catalog check. That is the blind spot
this reference closes.

THE THREE DESIGN DECISIONS, AND WHY
-----------------------------------
1. OURS WITH HIS CASSETTE, NOT HIS FILE. Keele's SHIVAD8EOM is LTR-trimmed:
   598bp off the 5' end and 135bp off the 3' relative to MN816822.1. Adopting it
   wholesale would lose LTR-mapping reads, and this project already established
   that unannotated LTRs were a major loss channel (shiv_mapping_diagnostic.py).
   The hybrid keeps our full-length contig and gains the cassette. His file is
   still built as variant V1, as a control that was actually run.

2. THE CORE IS N, THE CONSTANT FLANKS ARE NOT. This decides whether any of it
   works. 34 N's cost a spanning read 34 mismatches (0.227 of a 150bp read),
   past any sane filter. Writing in the 24 known constant bases and leaving only
   the 10bp core as N drops that to 10 (0.067).

3. TWO REGISTERS, TWO REFERENCES, UNION THE ALIGNMENTS. The cassette occurs in
   two registers in the Binhua data (ACGCGCAGCGTA..CAGCTAGCGCGT x755;
   ACGCGCTAGCTG..TACGCTGCGCGT x477), which differ at 10 of 12 flank positions.
   A read in the wrong register pays 10 N + 10 register = 20 mismatches. Rather
   than relax the budget to cover both with one reference, build both and union:
   a matched read then costs only 10, leaving headroom for ~30% hypermutation
   inside a 0.20 budget. The index is 10kb, so the second pass is free.

WHAT THIS SCRIPT DOES NOT DO
----------------------------
No alignment. STAR lives in the separate `star` env (CellRanger's bundled STAR
is too old for STARsolo flags), so this writes the exact commands to a shell
snippet rather than shelling out from sc_pre and failing on PATH. Run that
snippet, then point shiv_starsolo_recovery.sh at the new index.

Conventions: Spyder cells (# %%), config in Cell 1, warn-and-continue.

Author: Jake Lehle
Date: October 2026
"""

# %% ============================================================
# CELL 1 — CONFIG AND DIALS
# ============================================================
import os
import re
import math
from collections import Counter

import pandas as pd

# --- paths -----------------------------------------------------------------
WORKING   = '/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV'
REFDIR    = '/master/jlehle/WORKING/SC/REF'
STARSOLO  = os.path.join(WORKING, 'shiv_starsolo')

# Ours (MN816822.1). MKREF writes SHIVAD8EO.fna to REFDIR; the STARsolo recovery
# keeps its own copy. Try both, warn and continue.
OURS_CANDIDATES = [
    os.path.join(REFDIR, 'SHIVAD8EO.fna'),
    os.path.join(STARSOLO, 'ref', 'SHIVAD8EO.fa'),
]
# Brandon's barcoded reference, as sent Oct 2026.
KEELE_FA  = os.path.join(WORKING, 'SHIVAD8EOM_reference.txt')

OUT_DIR   = os.path.join(WORKING, 'shiv_bc_reference')

# --- Keele's flanking virus sequence (from his message, Oct 2026) ----------
KEELE_FLANK_5 = 'GGAGACCAGGACCTCCTCCTCCTCCCCCTCCAGGACTAGCATAA'   # 44bp
KEELE_FLANK_3 = 'ATGGAAGAAAGACCTCCAGAAAATGAAGGACCACAAAGGGAA'     # 42bp

# The locus is DERIVED in CELL 3, not taken from here. This is the assertion
# target: two independent routes (the abutting flanks, and the primer offsets)
# both give 5862, so a mismatch means an input changed and the build stops.
EXPECTED_LOCUS = 5862

# --- Khanal primers, for the comparison artifact ---------------------------
VPXF1  = 'CTAGGGGAAGGACATGGGGCAGG'
VPR_RT = 'CAGGTTGGCCGATTCTGGAGT'

# --- cassette geometry -----------------------------------------------------
# 34bp = 12bp constant flank + 10bp core + 12bp constant flank. CFLANK was
# missing here in the first version and only defined in
# shiv_check_bc_reference_alignment.py, which crashed CELL 6 after every output
# file had already been written. The assertion keeps the three constants
# consistent if any one of them is ever edited.
BC_LEN   = 34
CFLANK   = 12
CORE_LEN = 10
assert BC_LEN == 2 * CFLANK + CORE_LEN, (
    f"cassette geometry inconsistent: {CFLANK} + {CORE_LEN} + {CFLANK} "
    f"!= {BC_LEN}")
REGISTERS = {
    'r1': ('ACGCGCAGCGTA', 'CAGCTAGCGCGT'),   # 755 in the Binhua workbook
    'r2': ('ACGCGCTAGCTG', 'TACGCTGCGCGT'),   # 477
}

# --- variants to build -----------------------------------------------------
# MEASURED 2026-10-27 by shiv_check_bc_reference_alignment.py: STAR does NOT
# count reference N as a mismatch. A matched-register read spans V3 at nM=0 and
# the opposite register at exactly 10, which is the register flank difference
# and nothing for the 10 N core.
#
# That inverts the original design. More N is MORE permissive, not less, so the
# register-specific pair is now the only variant that penalises anything: it
# charges an opposite-register read 10 mismatches for no benefit. V4 is the
# production reference: our full-length contig with the whole cassette as N,
# which is Keele's treatment on our backbone without his LTR trimming. Both
# registers align at zero cost in ONE pass, and specificity is recovered
# downstream by checking the extracted 34 bases for cassette structure.
BUILD_V1_KEELE_ASIS = True    # control: his file as-is, LTR-trimmed
BUILD_V2_ALLN       = True    # conservative: 22 N, MluI half-sites kept
BUILD_V3_REGISTERS  = False   # superseded; set True only to reproduce the test
BUILD_V4_FULLN      = True    # PRODUCTION: full-length contig, 34 N cassette
PRODUCTION_VARIANT  = 'V4fullN'

# --- STARsolo settings to emit ---------------------------------------------
READ_LEN        = 150
MISMATCH_BUDGET = 0.20        # NoverLmax; 30 mismatches on 150bp
CONTIG_BASE     = 'SHIVAD8EO'

# --- synthetic spanning reads for the N-handling test ----------------------
MAKE_TEST_FASTQ = True
TEST_OFFSETS    = [10, 30, 58, 85, 105]   # bp of 5' flank before the cassette
BINHUA = os.path.join(WORKING,
                      'Binhua_FINAL_Barcode_4-7-26_Rhesus_SHIV_Barcodes.xlsx')
BINHUA_ALT = os.path.join(WORKING, 'Binhua_FINAL_Barcode_4726.xlsx')

os.makedirs(OUT_DIR, exist_ok=True)
pd.set_option('display.width', 170)

_COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


def revcomp(s):
    return s.translate(_COMP)[::-1]


def hr(title):
    print("\n" + "=" * 84)
    print(title)
    print("=" * 84)


def read_fasta_single(path):
    """Return (name, sequence). Tolerates Keele's file, which has the header and
    the sequence on ONE line with no newline between them, as well as ordinary
    multi-line FASTA. A naive parser silently eats the first ~9 bases of his
    sequence as part of the header, which would shift every coordinate."""
    raw = open(path).read().strip()
    if not raw.startswith('>'):
        raise ValueError(f"{path}: not FASTA (no leading '>')")
    body = raw[1:]
    # Standard case: a newline ends the header and the rest is sequence.
    if '\n' in body:
        head, rest = body.split('\n', 1)
        cleaned = ''.join(rest.split()).upper()
        if cleaned and re.fullmatch(r'[ACGTN]+', cleaned):
            return head.strip(), cleaned
    # Single-line case: the header name runs straight into the sequence. Find
    # the last character that cannot be a base within the first 100, and start
    # the sequence after it.
    cut = 0
    for i, ch in enumerate(body[:100]):
        if ch.upper() not in 'ACGTN':
            cut = i + 1
    name = body[:cut].strip()
    seq = ''.join(body[cut:].split()).upper()
    if not re.fullmatch(r'[ACGTN]+', seq):
        raise ValueError(f"{path}: could not separate header from sequence")
    return name, seq


def write_fasta(path, name, seq, width=60):
    with open(path, 'w') as fh:
        fh.write(f">{name}\n")
        for i in range(0, len(seq), width):
            fh.write(seq[i:i + width] + "\n")


def write_collapsed_gtf(path, contig, length):
    """ONE feature spanning the whole contig, so every sense viral read counts
    regardless of which overlapping ORF it hits. Same design as
    shiv_starsolo_recovery.sh, which is what took library H from 7 counted
    reads to ~337."""
    attrs = (f'gene_id "{contig}"; transcript_id "{contig}_tx"; '
             f'gene_name "{contig}"; gene_biotype "protein_coding";')
    with open(path, 'w') as fh:
        for feat in ('gene', 'transcript', 'exon'):
            fh.write(f"{contig}\tcustom\t{feat}\t1\t{length}\t.\t+\t.\t"
                     f"{attrs}\n")


# %% ============================================================
# CELL 2 — LOAD BOTH REFERENCES
# ============================================================
hr("CELL 2 — load references")

ours_path = next((p for p in OURS_CANDIDATES if os.path.exists(p)), None)
if ours_path is None:
    print("  ERROR: our SHIVAD8EO FASTA not found at any of:")
    for p in OURS_CANDIDATES:
        print(f"    {p}")
    raise SystemExit(1)
ours_name, OURS = read_fasta_single(ours_path)
print(f"  ours   : {ours_path}")
print(f"           name '{ours_name}', {len(OURS):,} bp")

if not os.path.exists(KEELE_FA):
    print(f"  ERROR: Keele reference not found at {KEELE_FA}")
    raise SystemExit(1)
keele_name, KEELE = read_fasta_single(KEELE_FA)
print(f"  Keele  : {KEELE_FA}")
print(f"           name '{keele_name}', {len(KEELE):,} bp")

n_run = re.search(r'N+', KEELE)
if not n_run:
    print("  ERROR: no N-run in the Keele reference; expected the 34bp cassette")
    raise SystemExit(1)
print(f"  Keele N-run: {n_run.start()}..{n_run.end()} "
      f"({n_run.end() - n_run.start()} bp)")
if OURS.count('N'):
    print(f"  NOTE: our reference contains {OURS.count('N')} N bases already")


# %% ============================================================
# CELL 3 — VERIFY, AND DERIVE THE INSERTION LOCUS INDEPENDENTLY
# ============================================================
hr("CELL 3 — reference comparison and locus derivation")

# Route 1: in OUR reference the two flanks abut, because the cassette is absent.
# Where they meet IS the insertion point. This uses no primer and no prior
# coordinate.
i5 = OURS.find(KEELE_FLANK_5)
i3 = OURS.find(KEELE_FLANK_3)
if i5 < 0 or i3 < 0:
    print(f"  ERROR: Keele flank not found in our reference "
          f"(5' at {i5}, 3' at {i3}).")
    print("         Our reference is not the virus his flanks came from, or an "
          "input changed.")
    raise SystemExit(1)
p5_end = i5 + len(KEELE_FLANK_5)
print(f"  our 5' flank : {i5}..{p5_end}")
print(f"  our 3' flank : {i3}..{i3 + len(KEELE_FLANK_3)}")
if p5_end != i3:
    print(f"  ERROR: flanks do not abut in our reference "
          f"(5' ends {p5_end}, 3' starts {i3}, gap {i3 - p5_end}).")
    print("         Our reference may already carry an insert at this locus.")
    raise SystemExit(1)
LOCUS = p5_end
print(f"\n  insertion locus, derived from abutting flanks : {LOCUS}")

# Route 2: the primer offsets, entirely independent of route 1.
k_vpxf1, o_vpxf1 = KEELE.find(VPXF1), OURS.find(VPXF1)
k_vprrt, o_vprrt = KEELE.find(revcomp(VPR_RT)), OURS.find(revcomp(VPR_RT))
up_off = o_vpxf1 - k_vpxf1 if (k_vpxf1 >= 0 and o_vpxf1 >= 0) else None
dn_off = o_vprrt - k_vprrt if (k_vprrt >= 0 and o_vprrt >= 0) else None
if up_off is not None:
    locus_from_primers = n_run.start() + up_off
    print(f"  insertion locus, derived from primer offsets  : "
          f"{locus_from_primers}")
    if locus_from_primers != LOCUS:
        print(f"  *** WARN the two routes disagree "
              f"({LOCUS} vs {locus_from_primers}). Resolve before building.")

if LOCUS != EXPECTED_LOCUS:
    print(f"\n  *** WARN derived locus {LOCUS} != EXPECTED_LOCUS "
          f"{EXPECTED_LOCUS}.")
    print("      An input changed. Check the FASTAs before trusting the build.")
else:
    print(f"  agrees with EXPECTED_LOCUS ({EXPECTED_LOCUS})")

cmp_rows = [
    {'feature': 'total length', 'keele': len(KEELE), 'ours': len(OURS),
     'offset': len(OURS) - len(KEELE)},
    {'feature': 'VpxF1 (fwd)', 'keele': k_vpxf1, 'ours': o_vpxf1,
     'offset': up_off},
    {'feature': 'cassette start', 'keele': n_run.start(), 'ours': LOCUS,
     'offset': LOCUS - n_run.start()},
    {'feature': 'cassette end', 'keele': n_run.end(), 'ours': LOCUS,
     'offset': LOCUS - n_run.end()},
    {'feature': 'vpr RT primer (rev)', 'keele': k_vprrt, 'ours': o_vprrt,
     'offset': dn_off},
]
cmp_df = pd.DataFrame(cmp_rows)
cmp_df.to_csv(os.path.join(OUT_DIR, 'reference_comparison.csv'), index=False)
print("\n" + cmp_df.to_string(index=False))
print("\n  saved: reference_comparison.csv")

if up_off is not None and dn_off is not None:
    print(f"\n  upstream offset {up_off}, downstream offset {dn_off}, "
          f"difference {up_off - dn_off}")
    if up_off - dn_off == BC_LEN:
        print(f"  The difference is exactly the {BC_LEN}bp cassette, so the "
              "two references are\n  the same virus differing only by LTR "
              "trimming plus the cassette, with no\n  detectable indel between "
              "the primer anchors.")
    else:
        print(f"  *** WARN expected the difference to be {BC_LEN}. An indel "
              "may sit between\n      the two primer anchors.")

# Does our sequence actually equal his over the shared span?
shared_lo = max(0, up_off) if up_off else 0
k_span = KEELE[:n_run.start()]
o_span = OURS[shared_lo:shared_lo + len(k_span)]
if len(o_span) == len(k_span):
    diffs = [i for i, (a, b) in enumerate(zip(k_span, o_span)) if a != b]
    print(f"\n  5' shared span ({len(k_span)} bp up to the cassette): "
          f"{len(diffs)} mismatch(es)")
    if diffs[:10]:
        print(f"    first positions (Keele coords): {diffs[:10]}")


# %% ============================================================
# CELL 4 — BUILD THE REFERENCE VARIANTS
# ============================================================
hr("CELL 4 — build variants")

variants = {}

if BUILD_V1_KEELE_ASIS:
    nm = f'{CONTIG_BASE}_V1keele'
    variants[nm] = KEELE
    print(f"  V1 {nm:28s} {len(KEELE):,} bp  (his file as-is, 34 N, "
          f"LTR-trimmed)")

if BUILD_V2_ALLN:
    cassette = 'ACGCGC' + 'N' * 6 + 'N' * CORE_LEN + 'N' * 6 + 'GCGCGT'
    assert len(cassette) == BC_LEN
    nm = f'{CONTIG_BASE}_V2allN'
    variants[nm] = OURS[:LOCUS] + cassette + OURS[LOCUS:]
    print(f"  V2 {nm:28s} {len(variants[nm]):,} bp  (22 N, MluI half-sites "
          f"kept as anchors)")

if BUILD_V3_REGISTERS:
    for reg, (left, right) in REGISTERS.items():
        cassette = left + 'N' * CORE_LEN + right
        assert len(cassette) == BC_LEN, f"{reg}: cassette is {len(cassette)}bp"
        nm = f'{CONTIG_BASE}_V3{reg}'
        variants[nm] = OURS[:LOCUS] + cassette + OURS[LOCUS:]
        print(f"  V3 {nm:28s} {len(variants[nm]):,} bp  "
              f"({left}..{right}, SUPERSEDED: charges the opposite register 10)")

if BUILD_V4_FULLN:
    cassette = 'N' * BC_LEN
    nm = f'{CONTIG_BASE}_{PRODUCTION_VARIANT}'
    variants[nm] = OURS[:LOCUS] + cassette + OURS[LOCUS:]
    print(f"  V4 {nm:28s} {len(variants[nm]):,} bp  "
          f"(34 N, both registers at nM=0, full LTRs)  <-- PRODUCTION")

print()
for nm, seq in variants.items():
    fa = os.path.join(OUT_DIR, f'{nm}.fa')
    gtf = os.path.join(OUT_DIR, f'{nm}.gtf')
    write_fasta(fa, nm, seq)
    write_collapsed_gtf(gtf, nm, len(seq))
    # self-check: the flanks must still abut the cassette in the built file
    if not nm.endswith('V1keele'):
        a = seq.find(KEELE_FLANK_5) + len(KEELE_FLANK_5)
        b = seq.find(KEELE_FLANK_3)
        ok = (b - a == BC_LEN)
        print(f"    {nm}.fa + .gtf   cassette at {a}..{b}  "
              f"{'OK' if ok else '*** WRONG, gap is ' + str(b - a)}")
    else:
        print(f"    {nm}.fa + .gtf")


# %% ============================================================
# CELL 5 — SYNTHETIC SPANNING READS (the N-handling test)
# ============================================================
hr("CELL 5 — synthetic reads that span the cassette")

if not MAKE_TEST_FASTQ:
    print("  skipped (MAKE_TEST_FASTQ = False)")
else:
    gt = {'r1': [], 'r2': []}
    gt_path = BINHUA if os.path.exists(BINHUA) else BINHUA_ALT
    if os.path.exists(gt_path):
        try:
            from openpyxl import load_workbook
            wb = load_workbook(gt_path, read_only=True, data_only=True)
            for sh in wb.sheetnames:
                for row in wb[sh].iter_rows(values_only=True):
                    v = row[2] if len(row) > 2 else None
                    if not isinstance(v, str):
                        continue
                    s = v.strip().upper()
                    if len(s) != BC_LEN or not re.fullmatch(r'[ACGT]+', s):
                        continue
                    for reg, (l, r) in REGISTERS.items():
                        if s.startswith(l) and s.endswith(r):
                            gt[reg].append(s)
            print(f"  real cassettes from {os.path.basename(gt_path)}: "
                  f"r1 {len(gt['r1'])}, r2 {len(gt['r2'])}")
        except Exception as e:                      # noqa: BLE001
            print(f"  WARN  workbook parse failed: {type(e).__name__}: {e}")
    if not (gt['r1'] or gt['r2']):
        print("  WARN  no real cassettes available; synthesising cores instead")
        gt = {reg: [l + 'ACGTACGTAC' + r] for reg, (l, r) in REGISTERS.items()}

    recs = []
    for reg in ('r1', 'r2'):
        for bc in gt[reg][:3]:
            for off in TEST_OFFSETS:
                left = OURS[LOCUS - off:LOCUS]
                need = READ_LEN - len(left) - BC_LEN
                if need < 0:
                    continue
                right = OURS[LOCUS:LOCUS + need]
                read = left + bc + right
                if len(read) != READ_LEN:
                    continue
                recs.append((f'{reg}_off{off}_{bc[12:22]}', read))
    fq = os.path.join(OUT_DIR, 'spanning_test_R2.fastq')
    with open(fq, 'w') as fh:
        for name, seq in recs:
            fh.write(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n")
    print(f"  wrote {len(recs)} synthetic 150bp reads spanning the cassette")
    print(f"    {fq}")
    print("\n  Expected nM, now that reference N is known to be FREE "
          "(measured 2026-10-27):")
    print(f"    V4 full-N, either register   0   <-- production")
    print(f"    V2 all-N, either register    0   (unless the MluI half-sites "
          f"are edited)")
    print(f"    V1 Keele as-is               0   (but LTR-trimmed)")
    print(f"    V3 matching register         0")
    print(f"    V3 opposite register        10   (the register flank "
          f"difference)")
    print(f"    V0 current reference         cannot align contiguously")
    print("\n  Score the result with shiv_check_bc_reference_alignment.py, "
          "never with a\n  record count: the relaxed filters emit a record for "
          "any read with 30\n  matchable bases, so soft-clipping at the "
          "cassette scores as a pass.")


# %% ============================================================
# CELL 6 — EMIT THE STAR COMMANDS
# ============================================================
hr("CELL 6 — STAR index + test alignment commands")

nmm = int(READ_LEN * MISMATCH_BUDGET)
lines = [
    '#!/bin/bash',
    '# Generated by shiv_build_bc_reference.py. Run in the `star` env:',
    '#   STAR 2.7.11b. CellRanger\'s bundled STAR is too old for STARsolo flags.',
    'set -euo pipefail',
    'source ~/anaconda3/bin/activate star',
    'export PATH="${CONDA_PREFIX}/bin:${PATH}"',
    f'cd {OUT_DIR}',
    '',
]
for nm, seq in variants.items():
    # STAR: min(14, log2(GenomeLength)/2 - 1). Hardcoding 6 drew
    # "--genomeSAindexNbases 6 is too large for the genome size, which may
    # cause seg-fault at the mapping step". It survived 30 synthetic reads.
    # Do not find out on 1e8.
    sa_n = min(14, int(math.log2(len(seq)) / 2 - 1))
    lines += [
        f'# --- {nm} ({len(seq):,} bp) ---',
        f'mkdir -p idx_{nm}',
        f'STAR --runMode genomeGenerate --runThreadN 8 \\',
        f'     --genomeDir idx_{nm} \\',
        f'     --genomeFastaFiles {nm}.fa \\',
        f'     --sjdbGTFfile {nm}.gtf \\',
        f'     --genomeSAindexNbases {sa_n}',
        '',
        f'# smoke test: do synthetic spanning reads survive the N core?',
        f'STAR --runThreadN 4 --genomeDir idx_{nm} \\',
        f'     --readFilesIn spanning_test_R2.fastq \\',
        f'     --outFileNamePrefix test_{nm}_ \\',
        f'     --outSAMtype SAM \\',
        f'     --outFilterMismatchNmax {nmm} \\',
        f'     --outFilterMismatchNoverLmax {MISMATCH_BUDGET} \\',
        f'     --outFilterMismatchNoverReadLmax {MISMATCH_BUDGET} \\',
        f'     --outFilterScoreMinOverLread 0 \\',
        f'     --outFilterMatchNminOverLread 0 \\',
        f'     --outFilterMatchNmin 30',
        f'grep -c -v "^@" test_{nm}_Aligned.out.sam '
        f'|| echo "  {nm}: NOTHING ALIGNED"',
        '',
    ]
sh = os.path.join(OUT_DIR, 'build_and_test_indices.sh')
with open(sh, 'w') as fh:
    fh.write("\n".join(lines) + "\n")
os.chmod(sh, 0o755)
print(f"  wrote {sh}")
print(f"\n  mismatch budget: NoverLmax {MISMATCH_BUDGET} = {nmm} mismatches on "
      f"a {READ_LEN}bp read")
print("  Safe to relax here because the index is SHIV-only: there is no host "
      "sequence\n  to cross-map to, the asymmetry shiv_starsolo_recovery.sh "
      "already relies on.")

prod = f'{CONTIG_BASE}_{PRODUCTION_VARIANT}'
print(f"\n  THEN, for the real pass: ONE index, idx_{prod}.")
print("  Reference N is free, so both registers align at nM=0 and the two-pass")
print("  union is unnecessary. Keep soloStrand REVERSE and the relaxed length")
print("  filters from shiv_starsolo_recovery.sh.")
print("\n  Off each aligned read, by COORDINATE and never by string match:")
print(f"    cassette    {LOCUS}..{LOCUS + BC_LEN}   all {BC_LEN} bases")
print(f"    inner flank {LOCUS + 6}..{LOCUS + CFLANK}   gives the register")
print(f"    core        {LOCUS + CFLANK}..{LOCUS + CFLANK + CORE_LEN}")
print("  The register comes off the READ, not the reference, which is why the")
print("  reference can be all N. Check cassette structure AFTER extraction as")
print("  an evidence tier, not before it as an alignment constraint.")

hr("REFERENCE BUILD COMPLETE — inputs untouched.")
