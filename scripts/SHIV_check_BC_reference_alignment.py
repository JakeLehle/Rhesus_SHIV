#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Rhesus SHIV — Barcoded Reference Smoke Test, Properly Scored (READ-ONLY)
========================================================================
Reads the test SAM files written by build_and_test_indices.sh and answers the
question that script's `grep -c` could not.

WHY THIS EXISTS
---------------
The smoke test in build_and_test_indices.sh counted SAM records:

    V1keele 30 / V2allN 29 / V3r1 30 / V3r2 30   (of 30 synthetic reads)

That count is uninformative, and reading it as success would be a mistake. The
alignment ran with

    --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0
    --outFilterMatchNmin 30

so ANY read with 30 matchable bases on either side of the cassette emits a
record while soft-clipping the rest. Soft-clipping at the cassette is exactly
the failure we are trying to detect, and a record count scores it as a pass.
V1keele makes that obvious: 30 of 30 "aligned" against a 34 N cassette under a
30-mismatch ceiling, which cannot be a through-alignment.

WHAT ACTUALLY DECIDES IT
------------------------
Three things per read, all in the SAM already:

  1. Does the alignment SPAN the cassette, POS at or before the locus and the
     reference end at or past locus+34, with the cassette consumed by M rather
     than clipped or skipped.
  2. What is nM. A matched-register read against V3 should carry 10, one per N
     in the core. More means STAR is penalising reference N differently than the
     build assumed; fewer means it is not counting them at all.
  3. Does coordinate-based core extraction return the right 10 bases. Every
     synthetic read carries its true core in its QNAME, so this is checkable
     against a known answer. This is the exact operation the production
     extractor will perform, so testing it here on 30 reads with known truth is
     worth more than the alignment count.

A variant only passes if reads span the cassette AND the extracted core matches
the planted one. Spanning with a wrong core is worse than not spanning, because
it fails silently downstream.

Read-only. Spyder cells (# %%), config in Cell 1, warn-and-continue.

Author: Jake Lehle
Date: October 2026
"""

# %% ============================================================
# CELL 1 — CONFIG AND DIALS
# ============================================================
import os
import re
from collections import Counter, defaultdict

import pandas as pd

REF_DIR = ('/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV/shiv_bc_reference')
OUT_DIR = os.path.join(REF_DIR, 'smoke_test_report')

# Per-variant cassette start, 0-based. V1 is Keele's own frame; the rest are
# ours. Both were derived and cross-checked in shiv_build_bc_reference.py.
LOCUS = {
    'SHIVAD8EO_V1keele': 5264,
    'SHIVAD8EO_V2allN':  5862,
    'SHIVAD8EO_V3r1':    5862,
    'SHIVAD8EO_V3r2':    5862,
    'SHIVAD8EO_V4fullN': 5862,
}
# Updated 2026-10-27 after the first run of this script measured it: STAR does
# NOT charge a mismatch for reference N. Reference N is free; only real base
# differences are counted. So the only variant that costs anything is V3 read
# in the opposite register, which pays for the 10 register flank positions.
EXPECTED_NM = {            # for a MATCHED-register read
    'SHIVAD8EO_V1keele': 0,
    'SHIVAD8EO_V2allN':  0,
    'SHIVAD8EO_V3r1':    0,
    'SHIVAD8EO_V3r2':    0,
    'SHIVAD8EO_V4fullN': 0,
}
VARIANT_REGISTER = {       # which register each variant's cassette is written in
    'SHIVAD8EO_V3r1': 'r1',
    'SHIVAD8EO_V3r2': 'r2',
}

BC_LEN     = 34
CFLANK     = 12
CORE_LEN   = 10

os.makedirs(OUT_DIR, exist_ok=True)
pd.set_option('display.width', 200)
pd.set_option('display.max_columns', 60)

CIGAR_RE = re.compile(r'(\d+)([MIDNSHP=X])')


def hr(title):
    print("\n" + "=" * 90)
    print(title)
    print("=" * 90)


def ref_span(pos0, cigar):
    """Reference bases consumed, and the clipping at each end."""
    ops = CIGAR_RE.findall(cigar)
    ref = 0
    lead_clip = int(ops[0][0]) if ops and ops[0][1] in 'SH' else 0
    tail_clip = int(ops[-1][0]) if ops and ops[-1][1] in 'SH' else 0
    inner_skip = 0
    for n, op in ops:
        n = int(n)
        if op in 'MDN=X':
            ref += n
        if op in 'DN':
            inner_skip += n
    return pos0 + ref, lead_clip, tail_clip, inner_skip


def read_base_at(pos0, cigar, seq, target0):
    """Return the read base aligned to reference position target0, or None if
    that position is deleted, skipped or clipped. Walks the CIGAR rather than
    assuming an ungapped alignment, because a read that STAR chose to gap
    through the cassette must not silently yield shifted bases."""
    ops = CIGAR_RE.findall(cigar)
    r, q = pos0, 0
    for n, op in ops:
        n = int(n)
        if op in 'M=X':
            if r <= target0 < r + n:
                return seq[q + (target0 - r)]
            r += n
            q += n
        elif op in 'IS':
            q += n
        elif op in 'DN':
            if r <= target0 < r + n:
                return None
            r += n
        elif op == 'H':
            pass
    return None


def extract_core(pos0, cigar, seq, locus):
    """Pull the 10bp core BY COORDINATE. This is the production operation: it
    does not string-match the cassette, so an edited or divergent cassette still
    yields its core as long as the read aligned."""
    out = []
    for i in range(locus + CFLANK, locus + CFLANK + CORE_LEN):
        b = read_base_at(pos0, cigar, seq, i)
        if b is None:
            return None
        out.append(b)
    return ''.join(out)


# %% ============================================================
# CELL 2 — PARSE THE TEST SAMs
# ============================================================
hr("CELL 2 — parse test SAM files")

rows = []
for variant, locus in LOCUS.items():
    sam = os.path.join(REF_DIR, f'test_{variant}_Aligned.out.sam')
    if not os.path.exists(sam):
        print(f"  WARN  missing {os.path.basename(sam)}; skipped")
        continue
    n_rec = 0
    with open(sam) as fh:
        for line in fh:
            if line.startswith('@'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 11:
                continue
            n_rec += 1
            qname, flag, rname, pos, mapq, cigar = f[0], int(f[1]), f[2], \
                int(f[3]), f[4], f[5]
            seq = f[9]
            tags = {}
            for t in f[11:]:
                parts = t.split(':', 2)
                if len(parts) == 3:
                    tags[parts[0]] = parts[2]
            nm = tags.get('nM', tags.get('NM'))
            nm = int(nm) if nm is not None and nm.lstrip('-').isdigit() else None

            # QNAME is reg_offN_CORE, planted by shiv_build_bc_reference.py
            m = re.match(r'^(r[12])_off(\d+)_([ACGT]+)$', qname)
            reg, off, true_core = (m.group(1), int(m.group(2)), m.group(3)) \
                if m else (None, None, None)

            unmapped = bool(flag & 4)
            pos0 = pos - 1
            end0, lead, tail, skip = ref_span(pos0, cigar) if not unmapped \
                else (None, 0, 0, 0)
            spans = (not unmapped and pos0 <= locus
                     and end0 >= locus + BC_LEN and skip == 0)
            got_core = extract_core(pos0, cigar, seq, locus) \
                if (spans and not unmapped) else None

            rows.append({
                'variant': variant, 'qname': qname, 'register': reg,
                'offset': off, 'true_core': true_core,
                'unmapped': unmapped, 'pos0': pos0, 'cigar': cigar,
                'nM': nm, 'lead_clip': lead, 'tail_clip': tail,
                'ref_skip': skip, 'spans_cassette': spans,
                'core_extracted': got_core,
                'core_correct': (got_core == true_core) if got_core else False,
                'matched_register': (reg == VARIANT_REGISTER.get(variant))
                                    if reg and variant in VARIANT_REGISTER
                                    else None,
            })
    print(f"  {variant:22s} {n_rec} record(s)")

df = pd.DataFrame(rows)
if df.empty:
    print("\n  No SAM records parsed. Check REF_DIR.")
    raise SystemExit(1)
df.to_csv(os.path.join(OUT_DIR, 'smoke_test_per_read.csv'), index=False)
print(f"\n  saved: smoke_test_per_read.csv ({len(df)} rows)")


# %% ============================================================
# CELL 2b — WHERE DID THE MISSING READS GO (STAR's own accounting)
# ============================================================
hr("CELL 2b — mapping fate from Log.final.out")

# The 2026-10-27 run returned 30 / 29 / 27 records for V1keele / V2allN /
# V4fullN. V4 is the most permissive reference and produced the FEWEST records,
# which the design does not predict. V1keele and V4fullN carry an IDENTICAL 34 N
# cassette and differ only in backbone, so the gap cannot be about N handling.
# STAR already decomposed it; this reads its accounting rather than guessing.
LOG_FIELDS = [
    ('Number of input reads', 'input'),
    ('Uniquely mapped reads number', 'unique'),
    ('Number of reads mapped to multiple loci', 'multi'),
    ('Number of reads mapped to too many loci', 'too_many_loci'),
    ('Number of reads unmapped: too many mismatches', 'unmapped_mismatch'),
    ('Number of reads unmapped: too short', 'unmapped_short'),
    ('Number of reads unmapped: other', 'unmapped_other'),
]

log_rows = []
for variant in LOCUS:
    log = os.path.join(REF_DIR, f'test_{variant}_Log.final.out')
    if not os.path.exists(log):
        print(f"  WARN  missing {os.path.basename(log)}; skipped")
        continue
    vals = {'variant': variant}
    try:
        text = open(log).read()
    except Exception as e:                          # noqa: BLE001
        print(f"  WARN  {os.path.basename(log)}: {e}")
        continue
    for label, key in LOG_FIELDS:
        m = re.search(re.escape(label) + r'\s*\|\s*(\S+)', text)
        vals[key] = int(m.group(1)) if m and m.group(1).isdigit() else None
    log_rows.append(vals)

if not log_rows:
    print("  No Log.final.out files found; skipping the fate breakdown.")
else:
    lg = pd.DataFrame(log_rows)
    lg['accounted'] = lg[['unique', 'multi', 'too_many_loci',
                          'unmapped_mismatch', 'unmapped_short',
                          'unmapped_other']].fillna(0).sum(axis=1).astype(int)
    print(lg.to_string(index=False))
    lg.to_csv(os.path.join(OUT_DIR, 'smoke_test_mapping_fate.csv'),
              index=False)
    print("\n  saved: smoke_test_mapping_fate.csv")

    for _, r in lg.iterrows():
        lost = (r['input'] or 0) - (r['unique'] or 0)
        if lost:
            why = []
            if r['multi']:
                why.append(f"{r['multi']} multimapping")
            if r['too_many_loci']:
                why.append(f"{r['too_many_loci']} too many loci")
            if r['unmapped_short']:
                why.append(f"{r['unmapped_short']} too short")
            if r['unmapped_mismatch']:
                why.append(f"{r['unmapped_mismatch']} too many mismatches")
            if r['unmapped_other']:
                why.append(f"{r['unmapped_other']} other")
            print(f"\n  {r['variant']}: {lost} of {r['input']} not uniquely "
                  f"mapped -> {', '.join(why) if why else 'unexplained'}")

    print("\n  MULTIMAPPING is the one to watch for V4fullN. A 34bp run of N")
    print("  matches anything at no cost, so a read with little real flank on")
    print("  one side can place equally well at more than one offset inside the")
    print("  N block. That is a cost of maximum permissiveness, and it is why")
    print("  V2allN keeps the MluI half-sites: six fixed bases at each end pin")
    print("  the register without charging either register a mismatch.")
    print("  If V4 is losing reads to multimapping, V2allN is the better")
    print("  production reference, and CELL 6 ranks on that basis.")


# %% ============================================================
# CELL 3 — WHAT THE RECORD COUNT HID
# ============================================================
hr("CELL 3 — records vs through-alignments")

summ = (df.groupby('variant')
        .agg(records=('qname', 'size'),
             unmapped=('unmapped', 'sum'),
             spans_cassette=('spans_cassette', 'sum'),
             core_correct=('core_correct', 'sum'),
             mean_lead_clip=('lead_clip', 'mean'),
             mean_tail_clip=('tail_clip', 'mean'))
        .reset_index())
summ['clipped_not_spanning'] = summ['records'] - summ['spans_cassette']
print(summ.to_string(index=False))
summ.to_csv(os.path.join(OUT_DIR, 'smoke_test_by_variant.csv'), index=False)

print("\n  'records' is what build_and_test_indices.sh counted.")
print("  'spans_cassette' is what actually matters: POS at or before the "
      "locus, reference\n  end at or past locus+34, cassette consumed by M "
      "with no D or N skip.")
print("  'core_correct' is the production operation tested against the core "
      "planted in\n  each read name. A variant that spans but extracts the "
      "wrong core is WORSE than\n  one that does not span, because it fails "
      "silently downstream.")


# %% ============================================================
# CELL 4 — MISMATCH ACCOUNTING: HOW DOES STAR TREAT REFERENCE N
# ============================================================
hr("CELL 4 — nM against prediction")

sp = df[df['spans_cassette'] & df['nM'].notna()]
if sp.empty:
    print("  No spanning alignment carried an nM tag, so STAR's handling of")
    print("  reference N cannot be read off this run.")
else:
    g = (sp.groupby(['variant', 'matched_register'])['nM']
         .agg(['size', 'min', 'median', 'max']).reset_index())
    print(g.to_string(index=False))
    print("\n  predicted nM for a MATCHED-register spanning read:")
    for v, e in EXPECTED_NM.items():
        mr = sp['matched_register']
        keep = mr.isna() | (mr == True)       # noqa: E712 - object dtype
        obs = sp[(sp['variant'] == v) & keep]['nM']
        if len(obs):
            med = obs.median()
            verdict = ('as predicted' if abs(med - e) < 0.5 else
                       'STAR is NOT counting reference N as a mismatch'
                       if med < e - 0.5 else
                       'STAR penalises reference N MORE than assumed')
            print(f"    {v:22s} predicted {e:2d}  observed median {med:5.1f}  "
                  f"-> {verdict}")
        else:
            print(f"    {v:22s} predicted {e:2d}  no spanning reads")

print("\n  Why this matters for the real pass: the mismatch budget has to cover")
print("  the N core PLUS whatever divergence a real read carries. If STAR "
      "ignores\n  reference N, the budget is free for hypermutation. If it "
      "charges for them,\n  an APOBEC3-edited read spends 10 on the core "
      "before it spends anything on\n  its own G-to-A changes, and the budget "
      "has to be raised accordingly.")


# %% ============================================================
# CELL 5 — DOES IT DEPEND ON WHERE THE READ SITS
# ============================================================
hr("CELL 5 — spanning by offset (how much flank a read needs)")

piv = (df.pivot_table(index='offset', columns='variant',
                      values='spans_cassette', aggfunc='sum')
       .fillna(0).astype(int))
print(piv.to_string())
print("\n  offset = bp of 5' flank before the cassette in the synthetic read;")
print("  the remainder sits 3' of it. STAR seeds on the longer exact arm and")
print("  extends back through the cassette only while the extension score "
      "stays\n  positive, so a read with little flank on one side can fail to "
      "extend even\n  when the budget allows the mismatches. A floor here "
      "tells you the minimum\n  flank a real read needs, which is the same "
      "quantity MIN_VIRUS_FLANK was\n  guessing at in the string scan.")

if not df[df['spans_cassette']].empty:
    bad = df[df['spans_cassette'] & ~df['core_correct']]
    if len(bad):
        print(f"\n  *** {len(bad)} spanning read(s) extracted the WRONG core:")
        print(bad[['variant', 'qname', 'pos0', 'cigar', 'nM',
                   'true_core', 'core_extracted']].to_string(index=False))


# %% ============================================================
# CELL 6 — VERDICT
# ============================================================
hr("CELL 6 — verdict")

# Preference order, not a tie-break on core count. All variants that work will
# tie at 100%, so "most correct cores" picks whichever the dict yields first,
# which in the first version of this script recommended the LTR-trimmed
# control. Rank by what the reference is FOR: full-length backbone first, then
# permissiveness at the cassette.
PREFERENCE = ['SHIVAD8EO_V4fullN', 'SHIVAD8EO_V2allN',
              'SHIVAD8EO_V3r1', 'SHIVAD8EO_V3r2', 'SHIVAD8EO_V1keele']
NOT_FOR_PRODUCTION = {
    'SHIVAD8EO_V1keele': 'LTR-trimmed, 598bp short at the 5p end; control only',
    'SHIVAD8EO_V3r1': 'charges the opposite register 10 mismatches for nothing',
    'SHIVAD8EO_V3r2': 'charges the opposite register 10 mismatches for nothing',
}

working = []
for v in LOCUS:
    s = summ[summ['variant'] == v]
    if s.empty:
        continue
    n_sp = int(s['spans_cassette'].iloc[0])
    n_ok = int(s['core_correct'].iloc[0])
    n_rec = int(s['records'].iloc[0])
    note = NOT_FOR_PRODUCTION.get(v, '')
    print(f"  {v:22s} {n_rec:3d} records, {n_sp:3d} span, {n_ok:3d} "
          f"correct core{('   [' + note + ']') if note else ''}")
    if n_ok > 0:
        working.append(v)

print()
if working:
    # Rank on evidence first, preference only as a tiebreak. Uniquely-mapped
    # count is the measurement that separates the variants: V4fullN's 34bp N
    # run is maximally permissive but can let a read place at more than one
    # offset, and a multimapped read is a lost read. A variant that is excluded
    # from production for a structural reason (LTR-trimmed, or charges a
    # register) only wins if nothing else works.
    uniq = {}
    if log_rows:
        uniq = {r['variant']: (r.get('unique') or 0) for r in log_rows}

    def rank_key(v):
        return (
            v in NOT_FOR_PRODUCTION,                      # False sorts first
            -uniq.get(v, 0),                              # more unique = better
            PREFERENCE.index(v) if v in PREFERENCE else len(PREFERENCE),
        )

    ranked = sorted(working, key=rank_key)
    pick = ranked[0]
    print(f"  USE {pick}."
          + (f"  ({uniq[pick]} uniquely mapped)" if pick in uniq else ""))
    if len(working) > 1:
        others = ', '.join(f"{r}"
                           + (f" ({uniq[r]}u)" if r in uniq else "")
                           for r in ranked[1:])
        print(f"  Ranked on uniquely-mapped count, then on fitness for "
              f"production.\n  Also extracting correct cores: {others}")
    if uniq and pick in uniq and uniq[pick] < max(uniq.values()):
        best_u = max(uniq, key=uniq.get)
        print(f"\n  NOTE {best_u} mapped more reads uniquely "
              f"({uniq[best_u]} vs {uniq[pick]}) but is not suitable for\n  "
              f"production: {NOT_FOR_PRODUCTION.get(best_u, 'see above')}.")
    print("\n  Next: point shiv_starsolo_recovery.sh at the single index")
    print(f"  idx_{pick}, keep soloStrand REVERSE and the relaxed length")
    print("  filters, and pull the cassette by COORDINATE off each aligned")
    print("  read rather than string-matching it.")
else:
    print("  NO VARIANT PRODUCED A CORRECT CORE. Do not start the real pass.")
    print("  Read smoke_test_per_read.csv: if reads are clipping rather than")
    print("  spanning, STAR is refusing to extend through the N core and the")
    print("  answer is a local aligner (minimap2, already used by")
    print("  shiv_targeted_remap.sh) for the barcode locus rather than STAR.")

print("\n  CHECK BEFORE THE REAL PASS: the first build hardcoded")
print("  --genomeSAindexNbases 6 and STAR warned that 5 is correct for a ~10kb")
print("  genome, 'which may cause seg-fault at the mapping step'. It survived")
print("  30 synthetic reads. Do not find out on 1e8. shiv_build_bc_reference.py")
print("  now computes it from the contig length, so regenerate")
print("  build_and_test_indices.sh and rebuild the index you intend to use.")

hr("SMOKE TEST SCORED — read-only.")
