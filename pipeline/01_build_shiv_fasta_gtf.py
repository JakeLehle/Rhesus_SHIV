#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
01 — Build the rebuilt SHIV FASTA, verify it, then derive its GTF
==================================================================

Four stages, strictly in order, because the GTF must describe the FASTA that
actually landed on disk rather than the one we intended to write:

    CELL 4   build the FASTA
    CELL 5   verify the FASTA, re-read from disk, and STOP on any failure
    CELL 6   build the GTF FROM the verified on-disk FASTA
    CELL 7   verify the GTF against that same FASTA

Outputs:

    SHIVAD8EO_V4fullN.fa      10,296 bp, cassette inserted as 34 N
    SHIVAD8EO_V4fullN.gtf     ONE feature spanning the whole contig

02_mkref concatenates these with the unchanged Mmul10 FASTA and GTF and runs
`cellranger mkref`. Keeping that separate means the viral sequence and
annotation can be checked in seconds before committing hours to mkref and days
to `cellranger multi`.


WHY THE CASSETTE IS 34 N, NOT ACGCGC + 22 N + GCGCGT
-----------------------------------------------------
This reverses an earlier recommendation of V2allN, which was based on a wrong
model. Recording the reasoning so it is not relitigated.

The V2allN case was that `ACGCGC` and `GCGCGT`, the MluI ligation scar, are
conserved across the catalog (1,232 of 1,260 entries start with the first,
1,246 of 1,260 end with the second) and so cost a genuine barcode nothing
while giving STAR 12 real bases to seed on. Both halves of that are wrong.

1. The cost model was ungapped. It compared each catalog 34-mer against the
   34bp reference window position by position. That is the right model for the
   1,232 canonical entries and meaningless for the 32 that are not, because
   those molecules are SHORTER at this locus and align with a DELETION, not
   with substitutions. Measured substitution cost is not alignment cost.

2. Once a gap is involved, a real base in the window stops being an anchor and
   becomes a constraint on where the gap may be placed. Three catalog entries
   carry essentially no insert at all:

       ADE08M.15266   ACCTCCTCCTCCTCCCCCTCCAGGACTAGCATAA   34bp of pure viral
                                                            flank, zero insert
       ADE08M.4093    TCCTCCTCCCCCTCCAGGACTAGCATAAACGCGT   28bp flank + a
                                                            religated MluI site
       ADE08M.1       CTCCTCCCCCTCCAGGACTAGCATAAACGCGCGT   26bp flank + 8bp

   A read from ADE08M.15266 must delete through `ACGCGC` and `GCGCGT` under
   V2allN, bases that do not exist in that molecule. Under 34 N the window is
   free and the gap lands where the flanking virus supports it.

3. It demonstrably mis-registered a real barcode. ADE08M.11736 exists at 33bp,
   so a read from it carries one deletion. V2allN pins 5862..5867 as real, so
   the gap cannot open before 5868, which is the first free position and sits
   inside the register tag. The observed read was `41S38M1D71M` with `nM:0`
   and a deletion at exactly reference 5868, which made `canon_core` return
   ambiguous. The reference caused that.

4. The 12 bases were never needed for seeding. A 150bp read spanning a 34bp
   cassette carries roughly 58bp of real virus on each side, which is what
   STAR actually seeds on.

N is the conservative choice: it costs nothing and constrains nothing. Real
bases assert those bases ARE present, and for 2.5% of the known constructs
they are not. So the window is masked exactly as Brandon specified it:

    5'-GGAGACCAGGACCTCCTCCTCCTCCCCCTCCAGGACTAGCATAA
       NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
       ATGGAAGAAAGACCTCCAGAAAATGAAGGACCACAAAGGGAA-3'

Both 44bp and 42bp flanks are ALREADY present in MN816822.1 and abut each
other there, which is how CELL 3 derives the locus. Only the 34 N go in, so
this is a +34 change and not a +120 one.

CONSEQUENCE FOR THE BARCODE ARM, not fixed here
------------------------------------------------
shiv_starsolo_r1_recovery.sh and shiv_starsolo_bc_recovery.sh both ran against
V2allN, so every existing barcode alignment carries the same gap-placement
constraint described in point 3. Re-running those passes against this
reference is a separate job and it bears directly on whether the confirmed
barcode count is two or three.


WHY ONE COLLAPSED GTF FEATURE
------------------------------
Measured by 00_gtf_multiplicity_diagnostic.py on 2026-10-30. The current GTF
leaves 33.4% of the contig under two or more features and 11.6% under none,
because the converter in MKREF_Mmul10_SHIVAD8EO.sh collapses each GenBank
join() to outermost coordinates, inflating *tat* to 2,708bp and *rev* to
2,666bp so they swallow *vpu* and most of *env*. Those three counted zero
reads across all three alignment sources. Collapsing to one feature makes
every aligned viral read countable, a measured 2.08x / 2.76x / 1.69x in READS
across Cell Ranger, the R1 pass and the R2 pass.

NO NESTED CASSETTE FEATURE, DELIBERATELY. Reads are 150bp and the cassette is
34bp, so every cassette-covering read also covers flank. A nested feature
would make those reads overlap two features, and Cell Ranger discards anything
compatible with more than one gene, so it would make cassette reads LESS
countable. Barcode extraction stays a BAM operation.


NAMING
------
    contig     SHIVAD8EO_V4fullN
    gene_id    SHIVAD8EO_provirus
    gene_name  SHIVAD8EO

The bare name SHIVAD8EO already refers to three different sequences in this
project: the 10,262bp parental in 16 Cell Ranger BAMs, the same parental as V0
in 8 R2 BAMs, and Keele's 9,563bp file whose own header is ">SHIVAD8EO".
Reusing it for a fourth would mean a BAM header no longer identifies which
sequence produced it. gene_name stays the clean string SHIVAD8EO because it
becomes the var_name in the matrix, and no host gene collides with it.


DOWNSTREAM EDIT THIS FORCES (do not skip)
------------------------------------------
scanpy_analysis_script.py section 3 tags the virus with

    SHIV_GENES = ['gag','pol','vif','vpx','vpr','tat','rev','vpu','env','nef']
    adata.var['shiv'] = adata.var_names.isin(SHIV_GENES)

After the collapse none of those names exist, n_shiv becomes 0, and the
case-insensitive fallback below it matches host TAT, tyrosine aminotransferase,
already flagged by check_all_libraries_shiv.py as colliding with viral *tat*.
The run would NOT error. It would report per-cell SHIV UMI that are actually
tyrosine aminotransferase. Replace both with:

    adata.var['shiv'] = adata.var_names == 'SHIVAD8EO'


TO VERIFY IN THE H GATE, NOT ASSUMED HERE
------------------------------------------
The collapsed feature is plus-strand only, matching write_collapsed_gtf in
SHIV_build_BC_reference.py and shiv_starsolo_recovery.sh. That STARsolo run
needed `soloStrand REVERSE`. If Cell Ranger's 5' handling lands the other way,
a plus-only feature counts nothing and the gate reads as a GTF failure when it
is a strand setting. Check the antisense fraction on the viral contig first.

Environment
-----------
  conda activate sc_pre

Author: Jake Lehle
Date:   October 2026
"""

# %% ============================================================
# CELL 1 — CONFIGURATION (every tunable lives here)
# ============================================================
import os
import re
import hashlib
from collections import Counter

import pandas as pd

try:
    from shiv_config import WORKDIR, REFDIR                         # noqa: F401
except Exception:                                                   # noqa: BLE001
    WORKDIR = '/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV'
    REFDIR = '/master/jlehle/WORKING/SC/REF'

# --- input: the parental, unbarcoded MN816822.1 contig ---------------------
PARENTAL_CANDIDATES = [
    os.path.join(REFDIR, 'SHIVAD8EO.fna'),
    os.path.join(WORKDIR, 'shiv_starsolo', 'ref', 'SHIVAD8EO.fa'),
]

# --- cross-check: the already-built variant of the same design -------------
EXISTING_SAME_DESIGN = os.path.join(WORKDIR, 'shiv_bc_reference',
                                    'SHIVAD8EO_V4fullN.fa')

# --- barcode catalog, used only to justify the all-N window in the report --
BINHUA_CANDIDATES = [
    os.path.join(WORKDIR, 'Binhua_FINAL_Barcode_4-7-26_Rhesus_SHIV_Barcodes.xlsx'),
    os.path.join(WORKDIR, 'Binhua_FINAL_Barcode_4726.xlsx'),
]
BC_SEQ_COL, BC_NAME_COL = 2, 1

# --- naming ----------------------------------------------------------------
CONTIG_NAME = 'SHIVAD8EO_V4fullN'
GENE_ID = 'SHIVAD8EO_provirus'
GENE_NAME = 'SHIVAD8EO'
TRANSCRIPT_ID = f'{GENE_ID}_tx'
GTF_SOURCE = 'custom'           # must contain no spaces; mkref misparses them
GENE_BIOTYPE = 'protein_coding'
STRAND = '+'

# --- ground truth geometry (Keele via Zoey, Oct 2026) ----------------------
KEELE_FLANK_5 = 'GGAGACCAGGACCTCCTCCTCCTCCCCCTCCAGGACTAGCATAA'   # 44 bp
KEELE_FLANK_3 = 'ATGGAAGAAAGACCTCCAGAAAATGAAGGACCACAAAGGGAA'     # 42 bp
VPXF1 = 'CTAGGGGAAGGACATGGGGCAGG'
VPR_RT = 'CAGGTTGGCCGATTCTGGAGT'

BC_LEN, CFLANK, CORE_LEN = 34, 12, 10
assert BC_LEN == 2 * CFLANK + CORE_LEN, 'cassette geometry inconsistent'

# THE CASSETTE. Masked in full, exactly as Brandon specified it. See the
# header for why the MluI half-sites are NOT written in.
CASSETTE = 'N' * BC_LEN
assert len(CASSETTE) == BC_LEN

# Written out only so the verification can assert the window is fully masked.
MLUI_5, MLUI_3 = 'ACGCGC', 'GCGCGT'
REGISTERS = {
    'r1': ('ACGCGCAGCGTA', 'CAGCTAGCGCGT'),
    'r2': ('ACGCGCTAGCTG', 'TACGCTGCGCGT'),
}

# Assertion targets, not inputs. The locus is DERIVED; these are what it must
# come out as. A mismatch means an input changed and the build stops.
EXPECTED_LOCUS = 5862
EXPECTED_PARENTAL_LEN = 10262
EXPECTED_BUILT_LEN = EXPECTED_PARENTAL_LEN + BC_LEN     # 10,296
EXPECTED_CORE = (EXPECTED_LOCUS + CFLANK,
                 EXPECTED_LOCUS + CFLANK + CORE_LEN)    # 5874, 5884

FASTA_WIDTH = 60
OUT_DIR = os.path.join(WORKDIR, 'pipeline_results', '01_shiv_reference')

os.makedirs(OUT_DIR, exist_ok=True)
pd.set_option('display.width', 220)

_COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


def revcomp(s):
    return s.translate(_COMP)[::-1]


def hr(t):
    print('\n' + '=' * 94)
    print(t)
    print('=' * 94)


def read_fasta_single(path):
    """Return (name, sequence). Tolerates a header running straight into the
    sequence on one line, which is how Keele's file is formatted; a naive
    parser eats the first few bases as header and shifts every coordinate."""
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


def write_fasta(path, name, seq, width=FASTA_WIDTH):
    with open(path, 'w') as fh:
        fh.write(f'>{name}\n')
        for i in range(0, len(seq), width):
            fh.write(seq[i:i + width] + '\n')


MANIFEST = []
CHECKS = []


def note(key, value, source=''):
    MANIFEST.append({'key': key, 'value': value, 'source': source})


def check(stage, label, ok, detail=''):
    CHECKS.append({'stage': stage, 'check': label, 'pass': bool(ok),
                   'detail': detail})
    print(f'  [{"OK  " if ok else "FAIL"}] {label}'
          + (f'   {detail}' if detail else ''))
    return bool(ok)


def stop_if_failed(stage):
    bad = [c for c in CHECKS if c['stage'] == stage and not c['pass']]
    if bad:
        pd.DataFrame(CHECKS).to_csv(
            os.path.join(OUT_DIR, 'verification.csv'), index=False,
            lineterminator='\n')
        raise SystemExit(
            f'\n  STOP: {len(bad)} {stage} check(s) failed. '
            f'See verification.csv.\n'
            '        Nothing downstream of this point was written.')


# %% ============================================================
# CELL 2 — LOAD THE PARENTAL CONTIG
# ============================================================
hr('CELL 2 — parental contig')

par_path = next((p for p in PARENTAL_CANDIDATES if os.path.exists(p)), None)
if par_path is None:
    print('  ERROR: parental SHIVAD8EO FASTA not found at any of:')
    for p in PARENTAL_CANDIDATES:
        print(f'    {p}')
    raise SystemExit(1)

par_name, PARENTAL = read_fasta_single(par_path)
print(f'  source        : {par_path}')
print(f'  contig name   : {par_name}')
print(f'  length        : {len(PARENTAL):,} bp')
print(f'  N bases       : {PARENTAL.count("N")}')
print(f'  md5           : {hashlib.md5(PARENTAL.encode()).hexdigest()}')

note('parental_path', par_path)
note('parental_contig_name', par_name)
note('parental_length', len(PARENTAL))
note('parental_md5', hashlib.md5(PARENTAL.encode()).hexdigest())

if len(PARENTAL) != EXPECTED_PARENTAL_LEN:
    raise SystemExit(
        f'  STOP: parental is {len(PARENTAL):,} bp, expected '
        f'{EXPECTED_PARENTAL_LEN:,}.\n'
        '        Either the input changed, or this is already a barcoded\n'
        '        reference and building would insert a second cassette.')
if PARENTAL.count('N'):
    print(f'  WARN  parental already contains {PARENTAL.count("N")} N bases; '
          'expected 0')


# %% ============================================================
# CELL 3 — DERIVE THE INSERTION LOCUS (two independent routes)
# ============================================================
hr('CELL 3 — locus derivation, then assertion')

# Route 1: the two Keele flanks ABUT in the parental, because the cassette is
# absent. Where they meet IS the insertion point. Uses no primer, no prior.
i5 = PARENTAL.find(KEELE_FLANK_5)
i3 = PARENTAL.find(KEELE_FLANK_3)
if i5 < 0 or i3 < 0:
    raise SystemExit(
        f'  STOP: Keele flank not found in the parental (5\' {i5}, 3\' {i3}).\n'
        '        The parental is not the virus his flanks came from.')

p5_end = i5 + len(KEELE_FLANK_5)
gap = i3 - p5_end
print(f'  5\' flank  {i5:>6} .. {p5_end:<6} ({len(KEELE_FLANK_5)} bp)')
print(f'  3\' flank  {i3:>6} .. {i3 + len(KEELE_FLANK_3):<6} '
      f'({len(KEELE_FLANK_3)} bp)')
print(f'  gap between flanks : {gap}')

if gap != 0:
    raise SystemExit(
        f'  STOP: the flanks do not abut (gap {gap}). The parental may already\n'
        '        carry an insert at this locus. Building would double it.')

LOCUS = p5_end
print(f'\n  route 1, abutting flanks  : locus {LOCUS}')

o_vpxf1 = PARENTAL.find(VPXF1)
o_vprrt = PARENTAL.find(revcomp(VPR_RT))
print(f'  VpxF1 forward             : {o_vpxf1}')
print(f'  vpr RT primer (revcomp)   : {o_vprrt}')
if o_vpxf1 < 0 or o_vprrt < 0:
    print('  WARN  a Khanal primer is missing; route 2 cannot cross-check')
else:
    inside = o_vpxf1 < LOCUS < o_vprrt
    print(f'  route 2, primer anchors   : locus between them: {inside}')
    if not inside:
        print('  WARN  the derived locus is NOT between the two primers.')

if LOCUS != EXPECTED_LOCUS:
    raise SystemExit(
        f'  STOP: derived locus {LOCUS} != EXPECTED_LOCUS {EXPECTED_LOCUS}.\n'
        '        Every hardcoded coordinate in the barcode arm assumes 5862.')
print(f'  asserted against EXPECTED_LOCUS ({EXPECTED_LOCUS}): OK')

note('derived_locus', LOCUS, 'abutting Keele flanks')
note('vpxf1_parental', o_vpxf1)
note('vprrt_revcomp_parental', o_vprrt)


# %% ============================================================
# CELL 4 — BUILD THE FASTA
# ============================================================
hr('CELL 4 — build the FASTA')

print(f'  cassette      : {CASSETTE}')
print(f'  length        : {len(CASSETTE)} bp '
      f'({sum(c != "N" for c in CASSETTE)} real bases, '
      f'{CASSETTE.count("N")} N)')
print('  design        : fully masked, exactly as Brandon specified. No MluI')
print('                  half-sites, no register tags, no core. See header.')

BUILT = PARENTAL[:LOCUS] + CASSETTE + PARENTAL[LOCUS:]
cas_start, cas_end = LOCUS, LOCUS + BC_LEN
core_start, core_end = LOCUS + CFLANK, LOCUS + CFLANK + CORE_LEN

print(f'\n  built length  : {len(BUILT):,} bp  (parental + {BC_LEN})')
print(f'  cassette at   : {cas_start} .. {cas_end}  (half-open)')
print(f'  core at       : {core_start} .. {core_end}  '
      f'(bases {core_start}..{core_end - 1} inclusive)')

fa_path = os.path.join(OUT_DIR, f'{CONTIG_NAME}.fa')
write_fasta(fa_path, CONTIG_NAME, BUILT)
print(f'\n  wrote {os.path.basename(fa_path)}  '
      f'(contig {CONTIG_NAME}, {len(BUILT):,} bp)')

note('cassette_window', CASSETTE)
note('cassette_start', cas_start)
note('cassette_end', cas_end)
note('core_start', core_start)
note('core_end', core_end)
note('built_length', len(BUILT))
note('built_md5', hashlib.md5(BUILT.encode()).hexdigest())
note('fasta_path', fa_path)


# %% ============================================================
# CELL 5 — VERIFY THE FASTA (re-read from disk; stops on failure)
# ============================================================
hr('CELL 5 — verify the FASTA')

chk_name, FA = read_fasta_single(fa_path)
S = 'fasta'

check(S, 'contig name as configured', chk_name == CONTIG_NAME, chk_name)
check(S, 'round-trips byte for byte from disk', FA == BUILT)
check(S, f'length is {EXPECTED_BUILT_LEN:,}', len(FA) == EXPECTED_BUILT_LEN,
      f'{len(FA):,}')
check(S, 'only A/C/G/T/N present', bool(re.fullmatch(r'[ACGTN]+', FA)))

# --- the cassette window ---------------------------------------------------
win = FA[cas_start:cas_end]
check(S, f'cassette window is {BC_LEN} bp', len(win) == BC_LEN, str(len(win)))
check(S, 'cassette window is FULLY masked (Brandon spec)',
      set(win) == {'N'}, win)
check(S, 'no MluI half-site written into the window',
      MLUI_5 not in win and MLUI_3 not in win)
check(S, 'no register tag written into the window',
      not any(l in win or r in win for l, r in REGISTERS.values()))
check(S, 'total N in contig equals the cassette length',
      FA.count('N') == BC_LEN, str(FA.count('N')))
check(S, 'exactly one N run', len(re.findall(r'N+', FA)) == 1,
      str(len(re.findall(r'N+', FA))))

# --- the flanks and the rest of the contig ---------------------------------
check(S, "5' Keele flank immediately precedes the cassette",
      FA[cas_start - len(KEELE_FLANK_5):cas_start] == KEELE_FLANK_5)
check(S, "3' Keele flank immediately follows the cassette",
      FA[cas_end:cas_end + len(KEELE_FLANK_3)] == KEELE_FLANK_3)
check(S, 'sequence outside the cassette is unchanged from parental',
      FA[:cas_start] == PARENTAL[:LOCUS] and FA[cas_end:] == PARENTAL[LOCUS:])
check(S, 'VpxF1 unmoved (upstream of the insert)',
      FA.find(VPXF1) == o_vpxf1, str(FA.find(VPXF1)))
check(S, f'vpr RT primer shifted by exactly {BC_LEN}',
      FA.find(revcomp(VPR_RT)) == o_vprrt + BC_LEN,
      f'{o_vprrt} -> {FA.find(revcomp(VPR_RT))}')
check(S, 'core lands where shiv_bc_core_extract.py expects',
      (core_start, core_end) == EXPECTED_CORE,
      f'{core_start}..{core_end} vs {EXPECTED_CORE[0]}..{EXPECTED_CORE[1]}')

print(f'\n  cassette in context:')
print(f'    ...{FA[cas_start - 30:cas_start]}[{win}]{FA[cas_end:cas_end + 30]}...')

ctx_path = os.path.join(OUT_DIR, 'cassette_context.txt')
pad = 60
with open(ctx_path, 'w') as fh:
    fh.write(f'{CONTIG_NAME}  cassette locus in context\n')
    fh.write(f'cassette {cas_start}..{cas_end}, core {core_start}..{core_end}\n\n')
    fh.write(f'parental {LOCUS - pad}..{LOCUS + pad} (flanks abut)\n')
    fh.write(f'  {PARENTAL[LOCUS - pad:LOCUS]}|{PARENTAL[LOCUS:LOCUS + pad]}\n\n')
    fh.write(f'built {cas_start - pad}..{cas_end + pad}\n')
    fh.write(f'  {FA[cas_start - pad:cas_start]}[{win}]'
             f'{FA[cas_end:cas_end + pad]}\n')
print(f'  wrote {os.path.basename(ctx_path)}')

stop_if_failed(S)
print(f'\n  FASTA verified. The GTF in CELL 6 is derived from THIS file.')


# %% ============================================================
# CELL 5b — WHY FULLY MASKED, SHOWN AGAINST THE CATALOG
# ============================================================
hr('CELL 5b — the catalog entries the mask protects (reporting only)')

cat_path = next((p for p in BINHUA_CANDIDATES if os.path.exists(p)), None)
if cat_path is None:
    print('  catalog not found; this report is skipped. It changes nothing.')
else:
    try:
        from openpyxl import load_workbook
        wb = load_workbook(cat_path, read_only=True, data_only=True)
        rows = []
        for sh in wb.sheetnames:
            for r in wb[sh].iter_rows(values_only=True):
                nm = r[BC_NAME_COL] if len(r) > BC_NAME_COL else None
                v = r[BC_SEQ_COL] if len(r) > BC_SEQ_COL else None
                if isinstance(v, str):
                    s = v.strip().upper()
                    if re.fullmatch(r'[ACGT]+', s) and 20 <= len(s) <= 40:
                        rows.append((sh.split()[0], str(nm), s))

        canon = [r for r in rows
                 if len(r[2]) == BC_LEN and r[2].startswith(MLUI_5)
                 and r[2].endswith(MLUI_3)
                 and any(r[2].startswith(a) and r[2].endswith(b)
                         for a, b in REGISTERS.values())]
        noncanon = [r for r in rows if r not in canon]
        print(f'  catalog          : {os.path.basename(cat_path)}')
        print(f'  entries          : {len(rows):,}')
        print(f'  fully canonical  : {len(canon):,}')
        print(f'  NON-canonical    : {len(noncanon):,}  '
              f'({100 * len(noncanon) / max(1, len(rows)):.1f}%)')
        print()
        print('  These are the entries whose molecules are SHORTER than 34bp')
        print('  at this locus, so a read from them aligns with a DELETION.')
        print('  Any real base in the reference window would constrain where')
        print('  that gap may open. A fully masked window does not.')

        slid = [r for r in rows
                if any(r[2].startswith(KEELE_FLANK_5[-k:])
                       for k in range(8, BC_LEN))]
        print(f'\n  entries starting with >=8bp of the viral 5\' flank, i.e. '
              f'constructs\n  where little or no insert went in: {len(slid)}')
        for an, nm, s in sorted(slid, key=lambda x: -len(x[2]))[:6]:
            shown = s
            for k in range(BC_LEN - 1, 7, -1):
                if s.startswith(KEELE_FLANK_5[-k:]):
                    shown = f'[{s[:k]}]{s[k:]}'
                    break
            print(f'    {an}  {nm:<16} {shown}')

        lens = Counter(len(s) for _, _, s in rows)
        print(f'\n  length distribution : {dict(sorted(lens.items()))}')
        note('catalog_entries', len(rows))
        note('catalog_noncanonical', len(noncanon))
        pd.DataFrame(noncanon, columns=['animal', 'bc_name', 'bc_seq']).to_csv(
            os.path.join(OUT_DIR, 'catalog_noncanonical.csv'), index=False,
            lineterminator='\n')
        print('  wrote catalog_noncanonical.csv')
    except Exception as e:                                          # noqa: BLE001
        print(f'  WARN  catalog report failed: {type(e).__name__}: {e}')


# %% ============================================================
# CELL 6 — BUILD THE GTF FROM THE VERIFIED FASTA
# ============================================================
hr('CELL 6 — build the GTF from the verified on-disk FASTA')

# Re-read rather than reuse BUILT, so the GTF can only ever describe the file
# that is actually on disk. The span end comes from that length, never from a
# constant, so the two cannot drift apart.
gtf_name, GTF_SEQ = read_fasta_single(fa_path)
gtf_len = len(GTF_SEQ)
print(f'  read back     : {os.path.basename(fa_path)}')
print(f'  contig        : {gtf_name}')
print(f'  span will be  : 1..{gtf_len:,}  (from the FASTA, not hardcoded)')

attrs_gene = (f'gene_id "{GENE_ID}"; gene_name "{GENE_NAME}"; '
              f'gene_biotype "{GENE_BIOTYPE}";')
attrs_tx = (f'gene_id "{GENE_ID}"; transcript_id "{TRANSCRIPT_ID}"; '
            f'gene_name "{GENE_NAME}"; gene_biotype "{GENE_BIOTYPE}";')
attrs_exon = (f'gene_id "{GENE_ID}"; transcript_id "{TRANSCRIPT_ID}"; '
              f'gene_name "{GENE_NAME}"; gene_biotype "{GENE_BIOTYPE}"; '
              f'exon_number "1";')

gtf_path = os.path.join(OUT_DIR, f'{CONTIG_NAME}.gtf')
with open(gtf_path, 'w') as fh:
    for feat, attrs in (('gene', attrs_gene),
                        ('transcript', attrs_tx),
                        ('exon', attrs_exon)):
        fh.write(f'{gtf_name}\t{GTF_SOURCE}\t{feat}\t1\t{gtf_len}'
                 f'\t.\t{STRAND}\t.\t{attrs}\n')
print(f'\n  wrote {os.path.basename(gtf_path)}')
with open(gtf_path) as fh:
    for line in fh:
        print(f'    {line.rstrip()}')
note('gtf_path', gtf_path)
note('gtf_span_end', gtf_len, 'read from the built FASTA')


# %% ============================================================
# CELL 7 — VERIFY THE GTF AGAINST THAT FASTA
# ============================================================
hr('CELL 7 — verify the GTF')

S = 'gtf'
gtf_rows = [l.rstrip('\n').split('\t') for l in open(gtf_path) if l.strip()]

check(S, 'exactly 3 lines', len(gtf_rows) == 3, str(len(gtf_rows)))
check(S, 'all lines have 9 tab-separated fields',
      all(len(r) == 9 for r in gtf_rows))
check(S, 'feature types are gene / transcript / exon',
      [r[2] for r in gtf_rows] == ['gene', 'transcript', 'exon'])
check(S, 'seqname matches the FASTA contig exactly',
      all(r[0] == gtf_name for r in gtf_rows), gtf_name)
check(S, 'span starts at 1', all(int(r[3]) == 1 for r in gtf_rows))
check(S, 'span ends at the FASTA length',
      all(int(r[4]) == gtf_len for r in gtf_rows),
      f'{gtf_rows[0][4]} vs {gtf_len}')
check(S, 'no spaces in the source column (mkref misparses them)',
      all(' ' not in r[1] for r in gtf_rows))
check(S, 'strand is consistent', len({r[6] for r in gtf_rows}) == 1,
      gtf_rows[0][6])
check(S, 'every transcript has at least one exon',
      sum(1 for r in gtf_rows if r[2] == 'exon') >= 1)
check(S, 'gene_id present on every line',
      all(f'gene_id "{GENE_ID}"' in r[8] for r in gtf_rows))
check(S, 'transcript_id present on transcript and exon',
      all(f'transcript_id "{TRANSCRIPT_ID}"' in r[8]
          for r in gtf_rows if r[2] != 'gene'))
check(S, 'gene_name is the clean matrix name',
      all(f'gene_name "{GENE_NAME}"' in r[8] for r in gtf_rows), GENE_NAME)
check(S, 'gene_name does not collide with the contig name',
      GENE_NAME != gtf_name)
check(S, 'the single feature covers the cassette',
      1 <= cas_start + 1 and cas_end <= gtf_len)

stop_if_failed(S)

pd.DataFrame(CHECKS).to_csv(os.path.join(OUT_DIR, 'verification.csv'),
                            index=False, lineterminator='\n')
n_ok = sum(1 for c in CHECKS if c['pass'])
print(f'\n  {n_ok}/{len(CHECKS)} checks passed across both stages')


# %% ============================================================
# CELL 8 — CROSS-CHECK, MANIFEST, AND WHAT 02 DOES NEXT
# ============================================================
hr('CELL 8 — cross-check against the already-built variant')

if not os.path.exists(EXISTING_SAME_DESIGN):
    print(f'  not found: {EXISTING_SAME_DESIGN}')
    print('  Cross-check skipped. Not fatal, the FASTA is verified on its own.')
    note('crosscheck', 'SKIPPED, file not found')
else:
    ex_name, EXISTING = read_fasta_single(EXISTING_SAME_DESIGN)
    print(f'  existing : {EXISTING_SAME_DESIGN}')
    print(f'  contig   : {ex_name}   length {len(EXISTING):,} bp')
    if len(EXISTING) != len(FA):
        print(f'  *** length mismatch: {len(EXISTING):,} vs {len(FA):,}')
        note('crosscheck', 'LENGTH MISMATCH')
    else:
        d = [i for i, (a, b) in enumerate(zip(EXISTING, FA)) if a != b]
        if not d:
            print('\n  IDENTICAL, base for base.')
            print('  SHIV_build_BC_reference.py built this same sequence as')
            print('  V4fullN, so the design is reproduced independently.')
            note('crosscheck', 'IDENTICAL')
        else:
            print(f'\n  *** {len(d)} differing position(s): {d[:20]}')
            for i in d[:10]:
                print(f'      {i}: existing={EXISTING[i]} built={FA[i]}')
            note('crosscheck', f'{len(d)} differences')

hr('CELL 8 — manifest')
man = pd.DataFrame(MANIFEST)
man.to_csv(os.path.join(OUT_DIR, 'build_manifest.csv'), index=False,
           lineterminator='\n')
print(man.to_string(index=False))

hr('DONE — 01_build_shiv_fasta_gtf.py')
print(f'  outputs: {OUT_DIR}')
print(f'    {CONTIG_NAME}.fa     {len(FA):,} bp, cassette fully masked')
print(f'    {CONTIG_NAME}.gtf    1 feature spanning 1..{gtf_len:,}')
print('    verification.csv, build_manifest.csv, cassette_context.txt,')
print('    catalog_noncanonical.csv')
print('\n  Nothing in REF/ or Analysis/ was touched.')
print('\n  Next, 02_mkref:')
print('    cat GCF_003339765.1_Mmul_10_genomic.fna  this .fa   > combined.fa')
print('    cat GCF_003339765.1_Mmul_10_genomic.gtf  this .gtf  > combined.gtf')
print('    (host GTF stays full and unfiltered; the existing sed fixes for')
print('     spaces in the source column still apply to the host side)')
print('    cellranger mkref --genome=Mmul10_SHIVAD8EO_v11')
print('\n  And before stage 04: scanpy_analysis_script.py section 3 must become')
print("    adata.var['shiv'] = adata.var_names == 'SHIVAD8EO'")
print('  with the case-insensitive fallback deleted, or host TAT is counted')
print('  as virus without any error being raised.')
