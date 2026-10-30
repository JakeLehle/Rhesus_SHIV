#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Rhesus SHIV — Barcode Core Extraction and Cell Linkage (READ-ONLY)
===================================================================
The capstone. Pulls the viral barcode by COORDINATE off every cassette-spanning
alignment, deduplicates it properly, collapses truncated recoveries, links each
read to a cell, and unions all three instruments into one table.

WHY BY COORDINATE AND NOT BY STRING MATCH
-----------------------------------------
shiv_vpr_flank_scan.py finds the cassette by exact 20bp flank match, so an
APOBEC3-edited read breaks the match and vanishes. Here the aligner has already
placed the read, so the cassette is wherever contig positions 5862..5896 fall in
that read, whatever its sequence says. An edited cassette still yields its core.
That is the entire reason the alignment route was built.

WHAT IT CORRECTS
----------------
  TRUNCATED CORES. The flank scan accepts partial cores down to 6 bases, so one
  barcode recovered at three completeness levels looks like three barcodes. The
  2026-10-29 run returned ACCATT / ACCATTTCT / ACCATTTCTT, which is one barcode,
  not three. CELL 5 collapses prefix- and suffix-compatible cores (suffix too,
  because canonicalisation reverse-complements one register).

  PCR DUPLICATES. 14 reads collapsing onto 4 barcodes in library H means repeat
  sampling, but repeats are only evidence of clonal expansion if they are
  INDEPENDENT MOLECULES. Same CB plus same UMI is one molecule amplified.
  CELL 4 separates the two, and nothing downstream should be called clonality
  before it does.

  DOUBLE COUNTING ACROSS INSTRUMENTS. R1 alignment, R2 alignment and the flank
  scan overlap; a read found by all three is one read. CELL 7 unions on read id.

THE CB ORIENTATION IS TESTED, NOT ASSUMED
-----------------------------------------
shiv_r1_structure_check.py measured raw R1 CB as direct in whitelist
orientation (~69% vs 0.03%), but shiv_merge_barcode_diagnostic.py once put 458
recovered reads onto 2 called cells and that history is still unexplained. The
BAM's CB tag is CellRanger-corrected, which is a third case again. CELL 6 tries
both and reports which wins rather than picking one.

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
OUT_DIR = os.path.join(WORKING, 'annotation_output', 'bc_core_extract')

# Alignment runs to harvest. Both carry CB/UB tags from --outSAMattributes.
SOURCES = {
    'R1': os.path.join(WORKING, 'shiv_starsolo_r1', 'solo', 'Forward'),
    'R2': os.path.join(WORKING, 'shiv_starsolo_bc', 'solo', 'V2allN'),
}
# Object CASCADE, most-annotated first. The 2026-10-29 run joined only 1 of 9
# molecules against shiv_host_final_embedded.h5ad alone, which holds 89,844
# QC-PASSING cells out of 126,252 loaded. A droplet carrying a single viral
# read plausibly failed cell calling, so its barcode is absent from the final
# object but still present upstream. Each molecule is reported with WHICH
# object matched it, so a cell recovered only from the pre-QC object is never
# silently presented as a called cell.
OBJECTS = [
    ('final_embedded', os.path.join(WORKING, 'shiv_starsolo', 'merged',
                                    'shiv_host_final_embedded.h5ad')),
    ('tier2_annotated', os.path.join(WORKING, 'shiv_starsolo', 'merged',
                                     'shiv_host_tier2_annotated.h5ad')),
    ('merged_preQC', os.path.join(WORKING, 'shiv_starsolo', 'merged',
                                  'shiv_host_merged.h5ad')),
    ('adata_qc', os.path.join(WORKING, 'Analysis', 'adata_qc.h5ad')),
]
FLANK_HITS  = os.path.join(WORKING, 'annotation_output', 'vpr_flank_scan',
                           'flank_scan_hits.csv')
BINHUA      = os.path.join(WORKING,
                           'Binhua_FINAL_Barcode_4-7-26_Rhesus_SHIV_Barcodes.xlsx')
BINHUA_ALT  = os.path.join(WORKING, 'Binhua_FINAL_Barcode_4726.xlsx')

LIBS     = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
INFECTED = {'D', 'E', 'F', 'H'}
LIB_NAME = {
    'A': 'SHIV_Pre_PBMC',             'B': 'SHIV_Pre_LN',
    'C': 'SHIV_Pre_PBMC_IndianMixed', 'D': 'SHIV_21DPI_Pre_LN_Mixed',
    'E': 'SHIV_21DPI_PBMC',           'F': 'SHIV_21DPI_NonInf_LN_Mixed',
    'G': 'SHIV_Necropsy_PBMC',        'H': 'SHIV_Necropsy_LN',
}

# Cassette geometry in V2allN coordinates, 0-based.
LOCUS, BC_LEN, CFLANK, CORE_LEN = 5862, 34, 12, 10
MIN_CORE_BASES = 6
INNER_ORIENT1, INNER_ORIENT2 = 'AGCGTA', 'TAGCTG'

os.makedirs(OUT_DIR, exist_ok=True)
pd.set_option('display.width', 220)
pd.set_option('display.max_columns', 80)

_COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
CIGAR_RE = re.compile(r'(\d+)([MIDNSHP=X])')
revcomp = lambda s: s.translate(_COMP)[::-1]        # noqa: E731


def hr(t):
    print("\n" + "=" * 94)
    print(t)
    print("=" * 94)


SAMTOOLS = shutil.which('samtools')
if not SAMTOOLS:
    print("ERROR: samtools not on PATH. conda activate sc_pre and")
    print("       export PATH=\"${CONDA_PREFIX}/bin:${PATH}\"")
    raise SystemExit(1)


def read_base_at(pos0, cigar, seq, target0):
    """Read base aligned to reference position target0, or None if that
    position is deleted, skipped or clipped. Walks the CIGAR rather than
    assuming an ungapped alignment: one of the 2026-10-29 spanning reads is
    41S12M3D87M10S, and a naive offset would read the wrong bases from it."""
    r, q = pos0, 0
    for n, op in CIGAR_RE.findall(cigar):
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
    return None


def pull_window(pos0, cigar, seq, start, end):
    """Read bases aligned to [start,end). '-' marks a deleted position so the
    register stays in frame rather than silently shifting."""
    out = []
    for i in range(start, end):
        b = read_base_at(pos0, cigar, seq, i)
        out.append(b if b else '-')
    return ''.join(out)


def canon_core(cassette):
    """Register off the READ's inner flank, not the reference, which is why the
    reference can carry N there. Returns (core, register, ok)."""
    if len(cassette) < BC_LEN:
        return None, None, False
    inner = cassette[CFLANK - 6:CFLANK]
    core = cassette[CFLANK:CFLANK + CORE_LEN]
    if inner == INNER_ORIENT1:
        return core, 'r1', '-' not in core
    if inner == INNER_ORIENT2:
        return revcomp(core), 'r2', '-' not in core
    return core, 'ambiguous', False


# %% ============================================================
# CELL 2 — HARVEST CASSETTE-SPANNING ALIGNMENTS
# ============================================================
hr("CELL 2 — pull spanning reads and extract the cassette by coordinate")

rows = []
for src, d in SOURCES.items():
    if not os.path.isdir(d):
        print(f"  WARN  {src}: directory not found ({d})")
        continue
    for lib in LIBS:
        bam = os.path.join(d, f'{lib}_{LIB_NAME[lib]}_Aligned.sortedByCoord.out.bam')
        if not os.path.exists(bam):
            continue
        contig = 'SHIVAD8EO_V2allN'
        try:
            p = subprocess.run(
                [SAMTOOLS, 'view', '-F', '0x904', bam,
                 f'{contig}:{LOCUS - 200}-{LOCUS + BC_LEN + 200}'],
                capture_output=True, text=True, check=True)
        except subprocess.CalledProcessError as e:      # noqa: BLE001
            print(f"  WARN  {src}/{lib}: {e}")
            continue
        n = 0
        for line in p.stdout.splitlines():
            f = line.split('\t')
            if len(f) < 11:
                continue
            pos0, cigar, seq = int(f[3]) - 1, f[5], f[9]
            cassette = pull_window(pos0, cigar, seq, LOCUS, LOCUS + BC_LEN)
            if cassette.count('-') > BC_LEN - CFLANK - MIN_CORE_BASES:
                continue                       # not enough of it visible
            tags = {}
            for t in f[11:]:
                parts = t.split(':', 2)
                if len(parts) == 3:
                    tags[parts[0]] = parts[2]
            core, reg, full = canon_core(cassette)
            if core is None or core.count('-') > CORE_LEN - MIN_CORE_BASES:
                continue
            n += 1
            rows.append({
                'source': src, 'lib': lib,
                'group': 'INF' if lib in INFECTED else 'neg',
                'read_id': f[0], 'pos0': pos0, 'cigar': cigar,
                'nM': int(tags.get('nM', tags.get('NM', -1))),
                'CB': tags.get('CB', '-'), 'UB': tags.get('UB', '-'),
                'cassette': cassette, 'core': core.replace('-', ''),
                'register': reg, 'core_full': full,
            })
        if n:
            print(f"  {src}/{lib}: {n} cassette read(s)")

df = pd.DataFrame(rows)
if df.empty:
    print("\n  No cassette-bearing alignments found. Check SOURCES and LOCUS.")
    raise SystemExit(0)
df.to_csv(os.path.join(OUT_DIR, 'cassette_reads_raw.csv'), index=False)
print(f"\n  {len(df)} cassette read(s) total")
print(df[['source', 'lib', 'group', 'pos0', 'cigar', 'nM', 'CB', 'UB',
          'core', 'register', 'core_full']].to_string(index=False))


# %% ============================================================
# CELL 3 — GROUND TRUTH, FOR VALIDATION ONLY
# ============================================================
hr("CELL 3 — Binhua catalog (validation, never a gate)")

gt = defaultdict(set)
gtp = BINHUA if os.path.exists(BINHUA) else BINHUA_ALT
if os.path.exists(gtp):
    try:
        from openpyxl import load_workbook
        wb = load_workbook(gtp, read_only=True, data_only=True)
        for sh in wb.sheetnames:
            animal = sh.split()[0]
            for row in wb[sh].iter_rows(values_only=True):
                v = row[2] if len(row) > 2 else None
                if not isinstance(v, str):
                    continue
                s = v.strip().upper()
                if len(s) != BC_LEN or not re.fullmatch(r'[ACGT]+', s):
                    continue
                c, _, ok = canon_core(s)
                if ok and c:
                    gt[c].add(animal)
        print(f"  {len(gt)} unique catalog cores from "
              f"{os.path.basename(gtp)}")
    except Exception as e:                              # noqa: BLE001
        print(f"  WARN  workbook parse failed: {e}")
else:
    print("  WARN  no Binhua workbook; GT columns stay empty")
GT_PRIVATE = {c: next(iter(a)) for c, a in gt.items() if len(a) == 1}


# %% ============================================================
# CELL 4 — DEDUPLICATION: READS, MOLECULES, INSTRUMENTS
# ============================================================
hr("CELL 4 — reads vs molecules")

df['molecule'] = df['CB'] + '|' + df['UB']
by_read = df.drop_duplicates(subset=['read_id'])
has_umi = by_read[(by_read['CB'] != '-') & (by_read['UB'] != '-')]
by_mol = has_umi.drop_duplicates(subset=['molecule'])

print(f"  alignment records            : {len(df)}")
print(f"  distinct reads               : {len(by_read)}")
print(f"  with both CB and UB          : {len(has_umi)}")
print(f"  distinct MOLECULES (CB+UMI)  : {len(by_mol)}")
print("\n  A read repeated at the same CB+UMI is one molecule amplified, not")
print("  independent sampling. Only the molecule count can speak to clonality.")

for g, sub in by_mol.groupby('group'):
    print(f"    {g}: {len(sub)} molecule(s) across "
          f"{sub['lib'].nunique()} librar(ies)")


# %% ============================================================
# CELL 5 — COLLAPSE TRUNCATED CORES
# ============================================================
hr("CELL 5 — barcode diversity, after collapsing truncations")

def collapse(cores):
    """Map each core onto the longest core it is a prefix OR suffix of. Suffix
    matters because canonicalisation reverse-complements register 2, so a
    truncated recovery can end up right-aligned rather than left-aligned."""
    uniq = sorted(set(cores), key=len, reverse=True)
    rep = {}
    for c in uniq:
        for longer in uniq:
            if len(longer) > len(c) and (longer.startswith(c)
                                         or longer.endswith(c)):
                rep[c] = rep.get(longer, longer)
                break
        else:
            rep[c] = rep.get(c, c)
    return rep

all_cores = list(by_mol['core']) if len(by_mol) else list(by_read['core'])
rep = collapse(all_cores)
by_mol = by_mol.copy()
by_mol['core_collapsed'] = by_mol['core'].map(rep)

merged = {k: v for k, v in rep.items() if k != v}
if merged:
    print("  truncated recoveries folded into their full-length parent:")
    for k, v in sorted(merged.items()):
        print(f"    {k:12s} -> {v}")
else:
    print("  no truncated cores to collapse")

print()
for g, sub in by_mol.groupby('group'):
    print(f"  {g}: {len(sub)} molecules, {sub['core'].nunique()} raw cores, "
          f"{sub['core_collapsed'].nunique()} distinct barcodes")
print("\n  per library:")
for (g, lib), sub in by_mol.groupby(['group', 'lib']):
    reps = Counter(sub['core_collapsed'])
    top = ', '.join(f'{c}x{n}' for c, n in reps.most_common())
    print(f"    {lib} [{g}]: {len(sub)} molecules -> "
          f"{len(reps)} barcode(s)   {top}")
print("\n  A barcode seen in >1 MOLECULE is independent re-sampling of the same")
print("  viral lineage. A barcode seen once is a single observation.")


# %% ============================================================
# CELL 6 — CELL LINKAGE
# ============================================================
hr("CELL 6 — link to cells")

FIELDS = ['animal_id', 'tissue', 'timepoint', 'tier1_celltype',
          'tier2_celltype', 'shiv_pos', 'has_tcr', 'leiden']

try:
    import anndata as ad
except ImportError:                                     # noqa: BLE001
    ad = None
    print("  WARN  anndata unavailable; cell linkage skipped")

# Load each object once and index it by (library, bare barcode).
loaded = []
for tag, path in OBJECTS:
    if ad is None:
        break
    if not os.path.exists(path):
        print(f"  (absent) {tag}: {path}")
        continue
    try:
        a = ad.read_h5ad(path)
    except Exception as e:                              # noqa: BLE001
        print(f"  WARN  {tag}: {type(e).__name__}: {e}")
        continue
    idx = defaultdict(set)
    for nm in a.obs_names:
        parts = nm.split('-')
        if len(parts) >= 3:
            idx[parts[-1]].add(parts[0])
    loaded.append((tag, a, idx))
    print(f"  {tag:16s} {a.shape[0]:>8,} cells   libraries "
          f"{''.join(sorted(idx))}")

linked = pd.DataFrame()
if not loaded:
    print("\n  No object loaded; cannot link.")
else:
    # Orientation, tested against the FIRST object that has any hits at all.
    res = {'direct': 0, 'revcomp': 0}
    for mode in res:
        for _, r in by_mol.iterrows():
            cb = r['CB'] if mode == 'direct' else revcomp(r['CB'])
            if any(cb in idx.get(r['lib'], ()) for _, _, idx in loaded):
                res[mode] += 1
    print(f"\n  CB join across the whole cascade: direct {res['direct']}"
          f"/{len(by_mol)}, revcomp {res['revcomp']}/{len(by_mol)}")
    MODE = max(res, key=res.get) if max(res.values()) else 'direct'
    print(f"  using {MODE}")

    recs = []
    for _, r in by_mol.iterrows():
        cb = r['CB'] if MODE == 'direct' else revcomp(r['CB'])
        rec = {**r.to_dict(), 'cb_joined': cb, 'found_in': None,
               'cell': None, 'is_called_cell': False}
        for tag, a, idx in loaded:
            if cb in idx.get(r['lib'], ()):
                key = f"{cb}-1-{r['lib']}"
                rec['found_in'] = tag
                rec['cell'] = key
                rec['is_called_cell'] = (tag == OBJECTS[0][0])
                if key in a.obs.index:
                    o = a.obs.loc[key]
                    for f in FIELDS:
                        if f in a.obs.columns and f not in rec:
                            rec[f] = o.get(f)
                break
        rec['gt_animal'] = GT_PRIVATE.get(r['core_collapsed'])
        recs.append(rec)

    linked = pd.DataFrame(recs)
    linked.to_csv(os.path.join(OUT_DIR, 'barcoded_cells.csv'), index=False)

    n_any = int(linked['found_in'].notna().sum())
    n_called = int(linked['is_called_cell'].sum())
    print(f"\n  molecules matched to a barcode in ANY object : "
          f"{n_any}/{len(linked)}")
    print(f"  of those, in the final QC-passing object      : {n_called}")
    print("  A molecule found only upstream sits in a droplet that FAILED cell")
    print("  calling. It is real sequence in a real barcode, but it is not a")
    print("  called cell and must not be reported as one.")

    show = [c for c in ['source', 'lib', 'group', 'CB', 'core_collapsed',
                        'register', 'core_full', 'found_in', 'animal_id',
                        'tissue', 'timepoint', 'tier1_celltype',
                        'tier2_celltype', 'shiv_pos', 'has_tcr', 'gt_animal']
            if c in linked.columns]
    print("\n  per molecule:")
    print(linked[show].to_string(index=False))


# %% ============================================================
# CELL 6b — CELL TYPE AND ANIMAL, BY LIBRARY
# ============================================================
hr("CELL 6b — cell type and animal per library")

if linked.empty:
    print("  nothing linked")
else:
    print("  molecules by library and animal:")
    if 'animal_id' in linked.columns:
        t = (linked.assign(animal_id=linked['animal_id'].fillna('unassigned'))
             .groupby(['group', 'lib', 'animal_id'])
             .agg(molecules=('CB', 'size'),
                  barcodes=('core_collapsed', 'nunique'),
                  cells=('cell', lambda s: s.notna().sum()))
             .reset_index())
        print(t.to_string(index=False))
        t.to_csv(os.path.join(OUT_DIR, 'by_library_animal.csv'), index=False)
    else:
        print("  (no animal_id column present in any loaded object)")

    if 'tier1_celltype' in linked.columns:
        print("\n  molecules by cell type:")
        ct = (linked.assign(
                  tier1_celltype=linked['tier1_celltype'].fillna('not a called cell'))
              .groupby(['group', 'lib', 'tier1_celltype'])
              .agg(molecules=('CB', 'size'),
                   barcodes=('core_collapsed', 'nunique'))
              .reset_index())
        print(ct.to_string(index=False))
        ct.to_csv(os.path.join(OUT_DIR, 'by_celltype.csv'), index=False)

    # Animal concordance, where a catalog-private core allows it.
    ok = linked[linked['gt_animal'].notna() & linked.get(
        'animal_id', pd.Series(dtype=object)).notna()] \
        if 'animal_id' in linked.columns else pd.DataFrame()
    if len(ok):
        conc = (ok['gt_animal'].astype(str) ==
                ok['animal_id'].astype(str)).mean()
        print(f"\n  animal concordance on catalog-private cores: "
              f"{conc:.0%} of {len(ok)}")
    else:
        print("\n  No linked core is catalog-private, so barcode-implied animal")
        print("  cannot be cross-checked against the genetic demultiplexing.")

    # A barcode in BOTH an infected and a negative library is contamination.
    shared = (linked.groupby('core_collapsed')['group']
              .agg(lambda s: set(s)))
    both = [c for c, g in shared.items() if {'INF', 'neg'} <= g]
    if both:
        print(f"\n  *** {len(both)} barcode(s) appear in BOTH an infected and a")
        print("  *** negative library, which is cross-library contamination:")
        for c in both:
            sub = linked[linked['core_collapsed'] == c]
            print(f"      {c}: " + ', '.join(
                f"{r['lib']}[{r['group']}]" for _, r in sub.iterrows()))
        print("  One flow cell, so index hopping is the simplest explanation.")
        print("  Use it to bound the contamination rate, and treat any single")
        print("  observation of such a barcode as unreliable.")


# %% ============================================================
# CELL 6c — BARCODE PROVENANCE: CELL, QC-DROPPED CELL, OR AMBIENT
# ============================================================
hr("CELL 6c — was each barcode ever a cell at all")

# The 2026-10-29 manual check found that 'not in the annotated object' covers
# two very different situations, and conflating them is the difference between
# a finding and a dead end:
#
#   CellRanger called it AND it is in the object  -> a cell, usable now
#   CellRanger called it, NOT in the object       -> YOUR QC dropped it,
#                                                    recoverable, and worth
#                                                    knowing which filter did it
#   only in raw_feature_bc_matrix                 -> never a cell; ambient RNA,
#                                                    almost certainly shed from
#                                                    a nearby infected cell
#
# Of 7 infected molecules, 1 was in the object, 1 was a CellRanger cell the QC
# removed, and 5 were ambient. That last group is not noise in the contamination
# sense: ACCATTTCTT appears in one real cell plus three ambient droplets, which
# is what one infected cell shedding RNA during dissociation looks like.
CELLRANGER_DIR = os.path.join(WORKING, 'Analysis')
SAMPLE_TABLE   = os.path.join(WORKING, 'sample_table.csv')

CR_SAMPLE = dict(LIB_NAME)      # library letter -> CellRanger sample name


def cr_barcodes(sample, kind):
    """Barcode set from a CellRanger matrix. kind is 'filtered' or 'raw'.
    Uses the top-level outs/ copies; the SC_MULTI_CS internals duplicate them."""
    p = os.path.join(CELLRANGER_DIR, sample, 'outs',
                     f'{kind}_feature_bc_matrix', 'barcodes.tsv.gz')
    if not os.path.exists(p):
        return None
    try:
        import gzip
        with gzip.open(p, 'rt') as fh:
            return {ln.strip().split('-')[0] for ln in fh if ln.strip()}
    except Exception as e:                              # noqa: BLE001
        print(f"    WARN  {sample}/{kind}: {e}")
        return None


cr_cache = {}
for lib in sorted(set(by_mol['lib'])):
    s = CR_SAMPLE[lib]
    cr_cache[lib] = {k: cr_barcodes(s, k) for k in ('filtered', 'raw')}
    f, r = cr_cache[lib]['filtered'], cr_cache[lib]['raw']
    print(f"  {lib} ({s}): filtered "
          f"{len(f) if f is not None else 'MISSING':>9}  raw "
          f"{len(r) if r is not None else 'MISSING':>9}")

# Provider genetic demultiplexing, keyed (library, bare barcode) -> animal.
# This is how scanpy_analysis_script.py assigns animal_id, and it works for a
# barcode whether or not that barcode survived QC, which is exactly what a
# QC-dropped cell needs.
bc2animal = {}
if os.path.exists(SAMPLE_TABLE):
    try:
        import scanpy as sc
        st = pd.read_csv(SAMPLE_TABLE)
        want = set(by_mol['lib'])
        for _, row in st.iterrows():
            lib = str(row.get('run_letter', '')).strip()
            if lib not in want:
                continue
            mp = row.get('matrix_path')
            if pd.isna(mp) or not os.path.exists(str(mp)):
                continue
            try:
                pa = sc.read_10x_h5(str(mp), gex_only=True)
                for bc in pa.obs_names:
                    bc2animal[(lib, bc.split('-')[0])] = str(row['animal_id'])
                del pa
            except Exception as e:                      # noqa: BLE001
                print(f"    WARN  provider {row.get('animal_id')}_{lib}: {e}")
        print(f"  provider demux entries for these libraries: "
              f"{len(bc2animal):,}")
    except Exception as e:                              # noqa: BLE001
        print(f"  WARN  provider demux unavailable: {type(e).__name__}: {e}")
else:
    print(f"  WARN  {SAMPLE_TABLE} not found; animal stays unassigned for")
    print("        any barcode that is not in the annotated object")

prov = []
for _, r in (linked if not linked.empty else by_mol).iterrows():
    cb = r.get('cb_joined', r['CB'])
    lib = r['lib']
    f = cr_cache.get(lib, {}).get('filtered')
    raw = cr_cache.get(lib, {}).get('raw')
    in_f = (cb in f) if f is not None else None
    in_r = (cb in raw) if raw is not None else None
    in_obj = bool(r.get('found_in'))
    if in_f and in_obj:
        status = 'cell (in object)'
    elif in_f and not in_obj:
        status = 'CELL, QC-DROPPED'
    elif in_r:
        status = 'ambient (raw only)'
    elif in_r is False:
        status = 'not in raw at all'
    else:
        status = 'unknown (matrix missing)'
    prov.append({
        'source': r['source'], 'lib': lib, 'group': r['group'], 'CB': cb,
        'core': r['core_collapsed'], 'register': r['register'],
        'cr_called_cell': in_f, 'in_raw': in_r, 'in_object': in_obj,
        'status': status,
        'animal_object': r.get('animal_id'),
        'animal_provider': bc2animal.get((lib, cb)),
        'celltype_object': r.get('tier1_celltype'),
    })

pv = pd.DataFrame(prov)
pv.to_csv(os.path.join(OUT_DIR, 'barcode_provenance.csv'), index=False)
print("\n  per molecule:")
print(pv.to_string(index=False))

print("\n  summary:")
print(pv.groupby(['group', 'status']).size().to_string())

qcd = pv[pv['status'] == 'CELL, QC-DROPPED']
if len(qcd):
    print(f"\n  *** {len(qcd)} barcode(s) were CellRanger-called cells that your")
    print("  *** QC removed. These are RECOVERABLE and should not be written")
    print("  *** off as ambient:")
    print(qcd[['lib', 'CB', 'core', 'animal_provider']].to_string(index=False))
    print("\n  Per ways-of-working the QC is per-sample one-sided MAD on counts")
    print("  and genes plus Scrublet, with MT% filtering OFF by default. So the")
    print("  likely cause is a low-count MAD cut or a doublet call. Pull these")
    print("  barcodes out of adata_qc.h5ad (or re-run QC with them whitelisted)")
    print("  to see which filter fired and whether the call was right.")

amb = pv[pv['status'] == 'ambient (raw only)']
if len(amb):
    print(f"\n  {len(amb)} ambient droplet(s). Ambient is not the same as")
    print("  contamination: a barcode shared between an ambient droplet and a")
    print("  real cell IN THE SAME LIBRARY is RNA shed by that cell, which")
    print("  supports the cell rather than undermining it. The same barcode in")
    print("  a NEGATIVE library is the contamination case.")
    for core, sub in pv.groupby('core'):
        cells = sub[sub['status'].str.startswith('cell') |
                    (sub['status'] == 'CELL, QC-DROPPED')]
        if len(cells) and len(sub) > len(cells):
            print(f"    {core}: {len(cells)} cell(s) + "
                  f"{len(sub) - len(cells)} ambient, libraries "
                  f"{sorted(set(sub['lib']))}")


# %% ============================================================
# CELL 7 — UNION WITH THE FLANK SCAN
# ============================================================
hr("CELL 7 — union across all three instruments")

if not os.path.exists(FLANK_HITS):
    print(f"  WARN  {FLANK_HITS} not found; alignment-only totals stand")
else:
    fs = pd.read_csv(FLANK_HITS)
    fs_called = fs[fs['tier'] != 'FLANK_ONLY'] if 'tier' in fs else fs
    aln_ids = set(by_read['read_id'])
    fs_ids = set(fs_called['read_id']) if 'read_id' in fs_called else set()
    both, aln_only, fs_only = (aln_ids & fs_ids, aln_ids - fs_ids,
                               fs_ids - aln_ids)
    print(f"  flank scan called reads : {len(fs_ids)}")
    print(f"  alignment cassette reads: {len(aln_ids)}")
    print(f"  found by BOTH           : {len(both)}")
    print(f"  alignment only          : {len(aln_only)}")
    print(f"  flank scan only         : {len(fs_only)}")
    print(f"  UNION                   : {len(aln_ids | fs_ids)}")
    print("\n  Alignment-only reads are the interesting ones: the exact 20bp")
    print("  flank match could not see them, which is what an APOBEC3-edited")
    print("  cassette looks like. Check their nM against the rest.")
    if aln_only:
        ao = by_read[by_read['read_id'].isin(aln_only)]
        print(f"\n  alignment-only, nM median {ao['nM'].median()}, "
              f"vs {by_read['nM'].median()} overall")

hr("EXTRACTION COMPLETE — read-only. Nothing modified.")
