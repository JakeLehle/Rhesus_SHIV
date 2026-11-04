#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
06 — Viral barcode detection from aligned reads
================================================

    in   pipeline_results/03_cellranger_multi/*/outs/**/*.bam
         pipeline_results/05_embed_annotate/adata_annotated.h5ad
         Binhua_FINAL_Barcode_*.xlsx                     (validation only)
    out  pipeline_results/06_viral_barcode/

THE METHOD, AND WHY IT IS NOW THIS SIMPLE
------------------------------------------
Every read that can carry a barcode is already in the BAM. The cassette is 34bp
inside a 10,296bp contig and reads are 150bp, so any read touching the cassette
carries at least 116bp of real virus sequence to anchor on. It aligns.

The old work had to scan raw FASTQ at full depth across eight libraries, and
that was never an alignment problem. It was the GTF: reads landed in the BAM
fine and were then discarded at COUNTING because they overlapped two or more
ORF features. Collapsing the contig to one gene feature in 01 fixed the
counting; the reads were always there. So the FASTQ scan was solving a problem
that no longer exists, and with it go MAX_READS_PER_LIB, PROGRESS_EVERY, the
smoke-test mode, the read caps and the hit cache. This stage reads a few
thousand SAM records and finishes in seconds.

Measured on library H (the richest, 34 of the 59 viral cells):

    sample_alignments.bam        862 records = 431 PAIRS, 13 near the cassette
    unassigned_alignments.bam     90 records =  45 pairs,  0 near the cassette

BOTH MATES ARE PRESENT AND THEY ARE ONE MOLECULE
-------------------------------------------------
431 R1 and 431 R2, exactly paired, because 03 ran SC5P-PE-v3 and Cell Ranger
aligns both mates in paired-end mode. The legacy scripts claim R1 is never
aligned, but that was true of their STARsolo run (CB_UMI_Simple reads nothing
past base 28), never of Cell Ranger multi in PE mode. Independent support: the
00 GTF diagnostic measured starsolo_r1 at 2.76x fold recovery, i.e. R1 reads
DID align to SHIV and were being lost at the GTF step.

R1 is a second look at the same molecule from the other end of the fragment,
which genuinely raises the chance of seeing the cassette per molecule. It does
NOT double the molecule count. Mates collapse on query name before anything is
counted.

CASSETTE GEOMETRY, VERIFIED NOT ASSUMED
----------------------------------------
01's cassette_context.txt shows Brandon's parental reference has NO cassette:
the flanks abut directly. Our 34 N is a landing pad we created, and Keele's
44bp and 42bp flanks are exactly the sequence either side of that junction.

    cassette  [5862, 5896)   34bp, all N in the reference
      MluI L  [5862, 5868)    6bp  ACGCGC    cloning scar, NOT used
      reg L   [5868, 5874)    6bp  AGCGTA or TAGCTG
      CORE    [5874, 5884)   10bp  <- the variable barcode
      reg R   [5884, 5890)    6bp  CAGCTA or TACGCT
      MluI R  [5890, 5896)    6bp  GCGCGT    cloning scar, NOT used

The MluI half-sites are out of the detection path on purpose. They are 6-mers
expected roughly every 4kb by chance, which is why the old
'ACGCGC.{22}GCGCGT' regex returned three reads in total. With the alignment
coordinate already localising the window they add nothing.

The register markers stay in, because they are what tells us which of the two
registers the core sits in: 707 catalog entries are register 1 and 435 are
register 2. They are read off the coordinate rather than searched for
genome-wide, so the specificity objection does not apply. They are masked in
the reference (everything is) but they are in the virus, present in 98.5% of
the catalog and documented in Khanal et al. 2019.

NO CANONICALISATION
-------------------
Register 1 and register 2 are the same cassette on opposite strands
(revcomp('AGCGTA') == 'TACGCT'), but the catalog holds 1,158 DISTINCT cores
across both, so they are different barcodes rather than one barcode written two
ways, and the catalog stores each 34mer AS-WRITTEN. So the extracted window is
matched directly. Both orientations are tried and which one hit is REPORTED, so
the first run tells us the convention instead of us assuming it.

MIN_CORE_BASES IS 10, NOT 6
----------------------------
The legacy scripts accepted partial cores down to 6 bases. Measured against the
catalog, a random fragment matches a real core:

      6 bases visible  41.7%          9 bases   0.88%
      7 bases          12.9%         10 bases   0.11%
      8 bases           3.4%

At six bases that is a coin flip, and it is how one barcode became three in the
2026-10-29 run (ACCATT / ACCATTTCT / ACCATTTCTT). Partials are reported and
never counted unless prefix- or suffix-compatible with a full-length core
already confirmed in the same animal.

THE CATALOG IS VALIDATION, NEVER A GATE, AND MATCHING IS EXACT
---------------------------------------------------------------
Brandon amplifies and mixes the whole stock before shipping, so the starting
pool is ~10,000 barcodes at unknown concentrations and nobody has the input
list. The Binhua workbook is a PLASMA readout at specific draw dates, 113 to
419 cores per animal. Most of what any animal carries was never in it.

Two consequences. A low catalog match rate is NOT evidence of false positives;
it is expected. And error-tolerant matching is off entirely: with 10,000
barcodes in a 10-mer space each has ~0.29 neighbours one base away, so roughly
a quarter of the real stock has a real neighbour at d=1 and a single sequencing
error converts one real barcode into another. Exact match only.

WHAT THIS STAGE IS FOR
-----------------------
The science question is how barcodes evolve over time, independent of the
starting mix. Viral reads exist only in D/E/F (21 DPI) and H (Necropsy), since
A/B/C are pre-infection and G came back with zero, so the longitudinal question
reduces to one thing: does any core appear at both timepoints in the same
animal?

At this depth the answer is probably no, and that is a finding rather than a
failure. Per-cell attribution is already known not to be publishable (05b: the
strongest B cell carried less virus than the strongest erythrocyte). What
survives is the per-animal barcode inventory, because the ambient RNA still
came out of genuinely infected cells in that animal. This stage ends with an
explicit verdict on whether targeted VpxF1 amplicon enrichment is required.

Environment:  conda activate sc_pre
Submit with:  sbatch run_06_viral_barcode.sh
Author: Jake Lehle
Date:   November 2026
"""

# %% ============================================================
# CELL 1 — CONFIGURATION (every tunable lives here)
# ============================================================
import os
os.environ.setdefault('MPLBACKEND', 'Agg')

import glob
import re
import subprocess
from collections import Counter, defaultdict

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import warnings
warnings.filterwarnings('ignore')

WORKDIR = '/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV'
MULTI_DIR = os.path.join(WORKDIR, 'pipeline_results', '03_cellranger_multi')
ANNOT_H5AD = os.path.join(WORKDIR, 'pipeline_results', '05_embed_annotate',
                          'adata_annotated.h5ad')
OUT_DIR = os.path.join(WORKDIR, 'pipeline_results', '06_viral_barcode')
FIG_DIR = os.path.join(OUT_DIR, 'figures')

# The catalog. Both filenames have been used; whichever exists wins.
CATALOG_CANDIDATES = [
    os.path.join(WORKDIR, 'Binhua_FINAL_Barcode_4-7-26_Rhesus_SHIV_Barcodes.xlsx'),
    os.path.join(WORKDIR, 'Binhua_FINAL_Barcode_4726.xlsx'),
]

# --- the contig ------------------------------------------------------------
# Read from the BAM header and CHECKED against this, never assumed. A silent
# contig rename is exactly the kind of thing that produces a believable zero.
CONTIG = 'SHIVAD8EO_V4fullN'

# --- cassette geometry, 0-based half-open, verified against 01 -------------
CASSETTE = (5862, 5896)        # 34bp
CORE = (5874, 5884)            # 10bp, the variable barcode
MLUI_L = (5862, 5868)          # ACGCGC, cloning scar, NOT used for detection
REG_L = (5868, 5874)           # AGCGTA or TAGCTG
REG_R = (5884, 5890)           # CAGCTA or TACGCT
MLUI_R = (5890, 5896)          # GCGCGT, cloning scar, NOT used
WINDOW_SLACK = 8               # bases either side of the cassette to extract

REGISTER_MARKERS = {'AGCGTA': 1, 'TAGCTG': 2}
USE_MLUI_AS_ANCHOR = False     # Jake's call 2026-11-03; see the header

# --- calling thresholds ----------------------------------------------------
MIN_CORE_BASES = 10            # was 6 in the legacy scripts; see the header
MIN_FLANK_ALIGNED = 15         # aligned (M) bases of virus either side
EXACT_MATCH_ONLY = True        # no d=1; 10K stock makes it unsafe
PROJECT_SOFTCLIPS = True       # cassette bases are usually soft-clipped

# --- libraries -------------------------------------------------------------
LIB_NAME = {
    'A': 'SHIV_Pre_PBMC',             'B': 'SHIV_Pre_LN',
    'C': 'SHIV_Pre_PBMC_IndianMixed', 'D': 'SHIV_21DPI_Pre_LN_Mixed',
    'E': 'SHIV_21DPI_PBMC',           'F': 'SHIV_21DPI_NonInf_LN_Mixed',
    'G': 'SHIV_Necropsy_PBMC',        'H': 'SHIV_Necropsy_LN',
}
# Measured in 05: viral cells appear only in D, E, F, H. A/B/C are
# pre-infection and G returned zero, so A/B/C/G are the negative floor and
# they run at full depth rather than being skipped.
INFECTED_LIBS = {'D', 'E', 'F', 'H'}
NEGATIVE_LIBS = {'A', 'B', 'C', 'G'}
CONTROL_ANIMALS = ['34315', '41903']   # Indian, uninfected, absent from catalog
TIME_ORDER = ['Pre', '21 DPI', 'Necropsy', 'Non-Infected', 'Unknown']

# --- SAM flags -------------------------------------------------------------
F_REVERSE = 0x10
F_READ1 = 0x40
F_READ2 = 0x80
F_SECONDARY = 0x100
F_SUPPLEMENTARY = 0x800
EXCLUDE_FLAGS = F_SECONDARY | F_SUPPLEMENTARY   # -F 0x900

# --- figures ---------------------------------------------------------------
FIG_DPI = 300
plt.rcParams.update({
    'font.size': 28, 'axes.titlesize': 34, 'axes.labelsize': 30,
    'xtick.labelsize': 28, 'ytick.labelsize': 28, 'legend.fontsize': 20,
    'figure.dpi': 100, 'pdf.fonttype': 42, 'ps.fonttype': 42,
})
HEX_ALIGNED, HEX_CLIP, HEX_GRID = '#1976D2', '#F9A825', '#E0E0E0'
HEX_INF, HEX_NEG = '#C62828', '#2E7D32'

os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)
pd.set_option('display.width', 250)
pd.set_option('display.max_columns', 60)
pd.set_option('display.max_rows', 300)

_COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


def revcomp(s):
    return s.translate(_COMP)[::-1]


def hr(t):
    print('\n' + '=' * 94, flush=True)
    print(t, flush=True)
    print('=' * 94, flush=True)


def save_csv(df, name, index=False):
    df.to_csv(os.path.join(OUT_DIR, name), index=index, lineterminator='\n')
    print(f'  wrote {name}  ({df.shape[0]:,} rows)', flush=True)


def savefig(fig, stem):
    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(FIG_DIR, f'{stem}.{ext}'), dpi=FIG_DPI,
                    bbox_inches='tight')
    plt.close(fig)
    print(f'  wrote {stem}.pdf / .png', flush=True)


def sh(cmd):
    """Run a command, return stdout. Raises with the stderr on failure."""
    p = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(f'command failed: {cmd}\n{p.stderr.strip()}')
    return p.stdout


_CIG = re.compile(r'(\d+)([MIDNSHP=X])')


def parse_cigar(c):
    return [(int(n), op) for n, op in _CIG.findall(c)] if c != '*' else []


def ref_to_query(pos1, cigar, project=True):
    """Map reference position (0-based) -> index into SEQ.

    Soft-clipped bases are PROJECTED onto the reference positions they would
    have occupied, which is the whole point: the cassette is N in the
    reference, so the aligner usually ends the alignment at the N-run boundary
    and puts the cassette bases into a soft clip. They are still in SEQ.

    The projection assumes no indels inside the clipped stretch, which is the
    only assumption available and is recorded per base as 'projected' so
    nothing downstream mistakes it for an aligned call.
    """
    m, src = {}, {}
    q = 0
    r = pos1 - 1                      # SAM POS is 1-based
    cig = parse_cigar(cigar)
    seen_aligned = False
    for i, (ln, op) in enumerate(cig):
        if op in 'M=X':
            for k in range(ln):
                m[r + k] = q + k
                src[r + k] = 'aligned'
            q += ln
            r += ln
            seen_aligned = True
        elif op == 'I':
            q += ln
        elif op in 'DN':
            r += ln
        elif op == 'S':
            if project:
                base = (r - ln) if not seen_aligned else r
                for k in range(ln):
                    if (base + k) not in m:
                        m[base + k] = q + k
                        src[base + k] = 'projected'
            q += ln
        elif op == 'H':
            pass                      # sequence absent from SEQ
        elif op == 'P':
            pass
    return m, src


def extract_window(seq, pos1, cigar, lo, hi, project=True):
    """Return (string, source_string, n_aligned) for reference [lo, hi)."""
    m, src = ref_to_query(pos1, cigar, project)
    out, osrc, n_al = [], [], 0
    for p in range(lo, hi):
        qi = m.get(p)
        if qi is None or qi >= len(seq):
            out.append('.')
            osrc.append('.')
        else:
            out.append(seq[qi])
            s = src.get(p, '.')
            osrc.append('A' if s == 'aligned' else 'p')
            if s == 'aligned':
                n_al += 1
    return ''.join(out), ''.join(osrc), n_al


def count_aligned_in(pos1, cigar, lo, hi):
    """Aligned (M/=/X) bases falling in reference [lo, hi). No projection."""
    m, src = ref_to_query(pos1, cigar, project=False)
    return sum(1 for p in range(lo, hi) if src.get(p) == 'aligned')


def register_anchored_core(seq, pos1, cigar, lo, hi, core_len, project=True):
    """Find a register marker in the RAW read sequence and take the core after
    it. Returns (register_number, core, query_index) or (0, '', -1).

    This reads the core out of SEQ directly rather than out of the
    reference-projected window, and the difference matters. The projected
    window has a '.' wherever a reference position has no query base, and
    stripping those dots splices sequence across a deletion, so a core read
    from the stripped window can be 10 bases that were never adjacent in the
    molecule. SEQ is contiguous by construction.

    The hit is constrained to the query indices the window actually spans, so
    a chance 6-mer elsewhere in a 150bp read cannot be picked up. That
    constraint is what makes a 6bp marker safe here; searching a whole read
    for one would not be.
    """
    m, _ = ref_to_query(pos1, cigar, project)
    qs = [m[p] for p in range(lo, hi) if p in m]
    if not qs:
        return 0, '', -1
    qlo, qhi = min(qs) - 6, max(qs) + 6
    best = (0, '', -1)
    for mk, rn in REGISTER_MARKERS.items():
        start = 0
        while True:
            j = seq.find(mk, start)
            if j < 0:
                break
            if qlo <= j <= qhi:
                cand = seq[j + len(mk):j + len(mk) + core_len]
                # a complete core wins outright; otherwise keep the longest
                if len(cand) == core_len:
                    return rn, cand, j
                if len(cand) > len(best[1]):
                    best = (rn, cand, j)
            start = j + 1
    return best


# %% ============================================================
# CELL 2 — DISCOVER THE BAMs
# ============================================================
hr('CELL 2 — BAM discovery')

# Globbed, never hardcoded. A guessed path already cost one round trip, and a
# Cell Ranger layout change should show up as a logged line rather than a
# traceback.
print('  searching for BAMs under', MULTI_DIR)
bams = []
for lib, sample in sorted(LIB_NAME.items()):
    sdir = os.path.join(MULTI_DIR, sample, 'outs')
    if not os.path.isdir(sdir):
        print(f'  library {lib} ({sample}): no outs/ directory, skipped')
        continue
    found = sorted(glob.glob(os.path.join(sdir, '**', '*.bam'), recursive=True))
    n_vdj = 0
    for b in found:
        # VDJ assembly BAMs (concat_ref.bam, consensus.bam under vdj_t/) have
        # one contig PER CLONOTYPE and no SHIV at all. The header check below
        # would exclude them anyway, but printing their near-miss list dumped
        # thousands of clonotype names into the 2026-11-04 log and buried
        # everything else. Skip them by path, and say how many were skipped.
        base = os.path.basename(b)
        if ('vdj' in b.split(os.sep)[-3:][0].lower()
                or base in ('concat_ref.bam', 'consensus.bam')
                or 'all_contig' in base or 'filtered_contig' in base):
            n_vdj += 1
            continue
        kind = ('sample' if b.endswith('sample_alignments.bam')
                else 'unassigned' if b.endswith('unassigned_alignments.bam')
                else 'other')
        bams.append({'library': lib, 'sample': sample, 'bam': b, 'kind': kind})
    if n_vdj:
        print(f'  library {lib}: skipped {n_vdj} VDJ assembly BAM(s) '
              '(per-clonotype contigs, no host or viral reference)')
bamdf = pd.DataFrame(bams)
if bamdf.empty:
    raise SystemExit(f'  STOP: no BAMs found under {MULTI_DIR}')

print(f'  {len(bamdf)} BAM(s) found across '
      f'{bamdf.library.nunique()} librar(ies)\n')

rows = []
for r in bamdf.itertuples():
    try:
        idx = sh(f'samtools idxstats {r.bam!r} 2>/dev/null')
    except RuntimeError:
        print(f'  {r.library} {os.path.basename(r.bam):<30} '
              'idxstats failed (unindexed?), skipped')
        continue
    contigs = [ln.split('\t')[0] for ln in idx.strip().split('\n') if ln]
    if CONTIG not in contigs:
        # Report a near-miss loudly, because a renamed contig is how a
        # believable zero gets produced. But TRUNCATE the list: a VDJ-style
        # header has one contig per clonotype and printing it in full buries
        # the log, which is what happened on 2026-11-04.
        near = [c for c in contigs if 'SHIV' in c.upper()
                and 'clonotype' not in c.lower()]
        if near:
            shown = ', '.join(near[:4]) + (' ...' if len(near) > 4 else '')
            print(f'  {r.library} {os.path.basename(r.bam):<30} '
                  f'WARN expected {CONTIG}, header has {len(near)} '
                  f'SHIV-like contig(s): {shown}')
        continue
    rows.append({**r._asdict(), 'has_contig': True})
usable = pd.DataFrame(rows).drop(columns=['Index'], errors='ignore')
if usable.empty:
    raise SystemExit(f'  STOP: no BAM contains contig {CONTIG}. '
                     'Check the reference used by 03.')
print(f'  {len(usable)} BAM(s) carry contig {CONTIG}')
print(usable[['library', 'kind', 'bam']].assign(
    bam=lambda d: d.bam.str.replace(MULTI_DIR, '...', regex=False)
).to_string(index=False))
save_csv(usable[['library', 'sample', 'kind', 'bam']], 'bams_used.csv')


# %% ============================================================
# CELL 3 — THE CATALOG (validation only, never a gate)
# ============================================================
hr('CELL 3 — Binhua catalog')

cat_path = next((p for p in CATALOG_CANDIDATES if os.path.exists(p)), None)
CATALOG = {}          # animal -> set of 34mers
CAT_CORE = {}         # animal -> set of 10bp cores
if cat_path is None:
    print('  WARN  no catalog file found. Looked for:')
    for p in CATALOG_CANDIDATES:
        print(f'          {p}')
    print('  Detection and counting are UNAFFECTED, because the catalog is')
    print('  validation only. The animal cross-check will be skipped.')
else:
    print(f'  reading {os.path.basename(cat_path)}')
    xl = pd.ExcelFile(cat_path)
    for sheet in xl.sheet_names:
        animal = sheet.split()[0]
        d = xl.parse(sheet)
        if 'bc_seq' not in d.columns:
            print(f'  sheet {sheet!r}: no bc_seq column, skipped')
            continue
        s = d['bc_seq'].dropna().astype(str).str.strip().str.upper()
        s = [x for x in s if set(x) <= set('ACGT')]
        CATALOG[animal] = set(s)
        CAT_CORE[animal] = {x[12:22] for x in s if len(x) == 34}
    allcore = set().union(*CAT_CORE.values()) if CAT_CORE else set()
    print(f'  animals in catalog : {sorted(CATALOG)}')
    print(f'  NOT in catalog     : {CONTROL_ANIMALS}  '
          '<- the uninfected Indian controls, as expected')
    rows = []
    cnt = Counter(c for v in CAT_CORE.values() for c in v)
    for a in sorted(CAT_CORE):
        priv = sum(1 for c in CAT_CORE[a] if cnt[c] == 1)
        rows.append({'animal': a, 'entries': len(CATALOG[a]),
                     'cores_10bp': len(CAT_CORE[a]), 'animal_private': priv,
                     'pct_private': round(100 * priv / max(len(CAT_CORE[a]), 1), 1)})
    cat_summary = pd.DataFrame(rows)
    print()
    print(cat_summary.to_string(index=False))
    save_csv(cat_summary, 'catalog_summary.csv')
    print(f'\n  {len(allcore):,} distinct cores across all animals, '
          f'{sum(1 for c in cnt if cnt[c] == 1):,} animal-private')
    print('  The stock Brandon ships is ~10,000 barcodes at unknown')
    print('  concentrations, so this is a PLASMA readout of a small fraction')
    print('  of what each animal carries. A core that is absent from it is')
    print('  not thereby suspect.')

    # register composition, which is why REGISTER_MARKERS stay in
    f34 = [x for v in CATALOG.values() for x in v if len(x) == 34]
    reg = Counter(x[6:12] for x in f34)
    print('\n  register composition of the catalog (inner 6-mer):')
    for k, n in reg.most_common(4):
        tag = f'  <- register {REGISTER_MARKERS[k]}' if k in REGISTER_MARKERS else ''
        print(f'    {k}  {n:>5}{tag}')


# %% ============================================================
# CELL 4 — HARVEST EVERY READ ALIGNED TO THE SHIV CONTIG
# ============================================================
hr('CELL 4 — harvest')

TAG_RE = re.compile(r'\b(CB|UB|NM):[AZi]:([^\t\s]+)')
recs = []
for r in usable.itertuples():
    sam = sh(f'samtools view -F 0x{EXCLUDE_FLAGS:X} {r.bam!r} {CONTIG}')
    n = 0
    for line in sam.split('\n'):
        if not line:
            continue
        f = line.split('\t')
        if len(f) < 11:
            continue
        tags = dict(TAG_RE.findall(line))
        flag = int(f[1])
        recs.append({
            'library': r.library, 'bam_kind': r.kind,
            'qname': f[0], 'flag': flag, 'pos': int(f[3]),
            'mapq': int(f[4]), 'cigar': f[5], 'seq': f[9],
            'mate': 1 if flag & F_READ1 else (2 if flag & F_READ2 else 0),
            'reverse': bool(flag & F_REVERSE),
            'CB': tags.get('CB', ''), 'UB': tags.get('UB', ''),
            'NM': int(tags.get('NM', -1)),
        })
        n += 1
    print(f'  {r.library} {r.kind:<11} {n:>6} records')

reads = pd.DataFrame(recs)
if reads.empty:
    print('\n  NO reads aligned to the contig anywhere. Nothing to do.')
    raise SystemExit(0)

reads['ref_end'] = [
    (p - 1) + sum(ln for ln, op in parse_cigar(c) if op in 'M=XDN')
    for p, c in zip(reads['pos'], reads['cigar'])]
reads['ref_start'] = reads['pos'] - 1

print(f'\n  {len(reads):,} records total')
mate = reads.groupby(['library', 'mate']).size().unstack(fill_value=0)
print('\n  records by library and mate (1 = R1, 2 = R2):')
print(mate.to_string())
print(f'\n  distinct query names (= PAIRS, i.e. molecules before UMI dedup): '
      f'{reads.qname.nunique():,}')
print('  R1 and R2 of a pair are the SAME molecule read from opposite ends.')
print('  R1 is a second look, not a second molecule.')

# cassette-proximal subset
prox = reads[(reads.ref_end >= CASSETTE[0] - 60)
             & (reads.ref_start <= CASSETTE[1] + 60)].copy()
print(f'\n  records within 60bp of the cassette: {len(prox):,} '
      f'({prox.qname.nunique():,} pairs)')
save_csv(reads.drop(columns=['seq']), 'reads_on_contig.csv')


# %% ============================================================
# CELL 5 — GEOMETRY: CAN THE CORE BE READ AT ALL?
# ============================================================
hr('CELL 5 — cassette geometry')

print('  The cassette is 34 N in the reference, so the aligner usually ends')
print('  the alignment at the N-run boundary and puts the cassette bases into')
print('  a soft clip. They are still in SEQ and still recoverable, but if')
print('  EVERY read stops dead at 5861 then the core comes entirely out of')
print('  projected clips and that is the finding, not a detail.')
print()

if prox.empty:
    print('  no cassette-proximal reads; nothing to characterise.')
    geom = pd.DataFrame()
else:
    g = []
    for r in prox.itertuples():
        core_s, core_src, core_al = extract_window(
            r.seq, r.pos, r.cigar, CORE[0], CORE[1], PROJECT_SOFTCLIPS)
        g.append({
            'library': r.library, 'qname': r.qname, 'mate': r.mate,
            'ref_start': r.ref_start, 'ref_end': r.ref_end,
            'cigar': r.cigar, 'reverse': r.reverse,
            'ends_at_cassette_L': abs(r.ref_end - CASSETTE[0]) <= 2,
            'starts_at_cassette_R': abs(r.ref_start - CASSETTE[1]) <= 2,
            'has_softclip': 'S' in r.cigar,
            'core_bases': sum(c != '.' for c in core_s),
            'core_aligned': core_al,
            'core_projected': sum(c == 'p' for c in core_src),
        })
    geom = pd.DataFrame(g)
    print(f'  {len(geom)} cassette-proximal records\n')
    print(f'  alignment ends within 2bp of the cassette 5\' edge (5862): '
          f'{int(geom.ends_at_cassette_L.sum())}')
    print(f'  alignment starts within 2bp of the 3\' edge (5896)       : '
          f'{int(geom.starts_at_cassette_R.sum())}')
    print(f'  records carrying any soft clip                          : '
          f'{int(geom.has_softclip.sum())}')
    print(f'\n  core bases recovered per record:')
    print(geom.core_bases.value_counts().sort_index().to_string())
    print(f'\n  of the recovered core bases: '
          f'{int(geom.core_aligned.sum())} aligned, '
          f'{int(geom.core_projected.sum())} projected from soft clips')
    if geom.core_projected.sum() > geom.core_aligned.sum():
        print('  MOSTLY PROJECTED. The CIGAR walk is doing all the work, which')
        print('  is expected with an N-run reference but means every core')
        print('  rests on the no-indel-in-the-clip assumption.')
    print(f'\n  records with a COMPLETE {CORE[1] - CORE[0]}bp core: '
          f'{int((geom.core_bases == CORE[1] - CORE[0]).sum())}')
    save_csv(geom, 'cassette_geometry.csv')

    if len(geom):
        fig, axes = plt.subplots(1, 2, figsize=(26, 10))
        axes[0].hist(geom.ref_end, bins=40, color=HEX_ALIGNED,
                     edgecolor='white')
        for x, lab in ((CASSETTE[0], "cassette 5'"), (CASSETTE[1], "cassette 3'")):
            axes[0].axvline(x, color=HEX_CLIP, lw=3, ls='--')
        axes[0].set_xlabel('alignment end (reference)')
        axes[0].set_ylabel('records')
        axes[0].set_title('Where alignments stop', fontsize=34)
        vc = geom.core_bases.value_counts().sort_index()
        axes[1].bar(vc.index, vc.values, color=HEX_ALIGNED, edgecolor='white')
        axes[1].axvline(MIN_CORE_BASES - 0.5, color=HEX_INF, lw=3, ls='--')
        axes[1].set_xlabel('core bases recovered')
        axes[1].set_ylabel('records')
        axes[1].set_title(f'Core completeness (gate at {MIN_CORE_BASES})',
                          fontsize=34)
        for ax in axes:
            ax.grid(color=HEX_GRID, lw=1)
            ax.set_axisbelow(True)
        savefig(fig, 'F06_cassette_geometry')


# %% ============================================================
# CELL 6 — EXTRACT CORES
# ============================================================
hr('CELL 6 — core extraction')

lo = CASSETTE[0] - WINDOW_SLACK
hi = CASSETTE[1] + WINDOW_SLACK
print(f'  extraction window: reference [{lo}, {hi})  '
      f'(cassette +/- {WINDOW_SLACK}bp)')
print('  the register marker is located IN THE READ\'S OWN SEQUENCE, with the')
print('  hit constrained to the query indices this window spans. A marker is')
print('  only 6bp, so searching a whole 150bp read for one would be unsafe;')
print('  the coordinate does the localising and the marker does the framing.')
print()
print('  NOTE ON STRAND. SAM stores SEQ in REFERENCE orientation, including')
print('  for reverse-strand alignments, so the CIGAR walk is identical for')
print('  both strands and nothing is reverse-complemented by FLAG. The design')
print('  notes said to revcomp minus-strand reads; that would have been a bug')
print('  and it is deliberately not done. The `reverse` column is recorded')
print('  for the log only. What DOES read as a reverse complement is a')
print('  cassette inserted in the opposite register in the virus itself, and')
print('  that is handled by trying both markers.')
print()

ext = []
for r in prox.itertuples():
    win, wsrc, _ = extract_window(r.seq, r.pos, r.cigar, lo, hi,
                                  PROJECT_SOFTCLIPS)
    core_s, core_src, core_al = extract_window(
        r.seq, r.pos, r.cigar, CORE[0], CORE[1], PROJECT_SOFTCLIPS)

    reg, core_from_reg, reg_at = register_anchored_core(
        r.seq, r.pos, r.cigar, lo, hi, CORE[1] - CORE[0], PROJECT_SOFTCLIPS)

    # prefer the register-anchored core when it is complete, because it is
    # read from contiguous SEQ and so survives an indel that shifts the
    # coordinate projection
    core = (core_from_reg if len(core_from_reg) == CORE[1] - CORE[0]
            else core_s.replace('.', ''))
    anchored = len(core_from_reg) == CORE[1] - CORE[0]

    ext.append({
        'library': r.library, 'bam_kind': r.bam_kind, 'qname': r.qname,
        'mate': r.mate, 'CB': r.CB, 'UB': r.UB, 'NM': r.NM, 'mapq': r.mapq,
        'reverse': r.reverse, 'cigar': r.cigar,
        'window': win, 'window_src': wsrc,
        'core': core.upper(), 'core_bases': len(core.replace('.', '')),
        'core_aligned': core_al, 'register': reg, 'register_anchored': anchored,
        'marker_qpos': reg_at,   # where in SEQ the marker sat; -1 = not found
        'flank5_aligned': count_aligned_in(r.pos, r.cigar,
                                           CASSETTE[0] - MIN_FLANK_ALIGNED,
                                           CASSETTE[0]),
        'flank3_aligned': count_aligned_in(r.pos, r.cigar, CASSETTE[1],
                                           CASSETTE[1] + MIN_FLANK_ALIGNED),
    })
E = pd.DataFrame(ext)
if E.empty:
    print('  no cassette-proximal records to extract from.')
else:
    print(f'  {len(E)} records processed')
    print(f'  register marker found : {int((E.register > 0).sum())}'
          f'   (reg1 {int((E.register == 1).sum())}, '
          f'reg2 {int((E.register == 2).sum())})')
    print(f'  register-anchored core: {int(E.register_anchored.sum())}')
    print(f'  complete {CORE[1] - CORE[0]}bp core : '
          f'{int((E.core_bases == CORE[1] - CORE[0]).sum())}')
    print(f'\n  flank support (aligned virus bases either side, '
          f'need >={MIN_FLANK_ALIGNED}):')
    print(f'    both flanks: {int(((E.flank5_aligned >= MIN_FLANK_ALIGNED) & (E.flank3_aligned >= MIN_FLANK_ALIGNED)).sum())}')
    print(f'    one flank  : {int(((E.flank5_aligned >= MIN_FLANK_ALIGNED) ^ (E.flank3_aligned >= MIN_FLANK_ALIGNED)).sum())}')
    print(f'    neither    : {int(((E.flank5_aligned < MIN_FLANK_ALIGNED) & (E.flank3_aligned < MIN_FLANK_ALIGNED)).sum())}')


# %% ============================================================
# CELL 7 — CONFIDENCE TIERS
# ============================================================
hr('CELL 7 — tiers')

CORE_LEN = CORE[1] - CORE[0]

# MIN_CORE_BASES has to be the thing the gate actually reads, or it is a knob
# that silently does nothing and the header's argument about 6 vs 10 bases is
# decoration. Validate it here rather than trusting it.
if MIN_CORE_BASES > CORE_LEN:
    raise SystemExit(f'  STOP: MIN_CORE_BASES={MIN_CORE_BASES} exceeds the '
                     f'core length {CORE_LEN}. Nothing could ever pass.')
if MIN_CORE_BASES < CORE_LEN:
    print(f'  WARN  MIN_CORE_BASES={MIN_CORE_BASES} is below the full core '
          f'length {CORE_LEN}, so PARTIAL cores will be counted.')
    print(f'        Measured random-match rate against the catalog at '
          f'{MIN_CORE_BASES} bases is roughly '
          f'{ {6: "41.7%", 7: "12.9%", 8: "3.4%", 9: "0.88%"}.get(MIN_CORE_BASES, "high") }.')
    print('        This is how one barcode became three in the 2026-10-29 run.')
print(f'  gate: core >= {MIN_CORE_BASES} of {CORE_LEN} bases, '
      f'register marker present, flanks >= {MIN_FLANK_ALIGNED}bp\n')

if E.empty:
    E['tier'] = []
else:
    both = ((E.flank5_aligned >= MIN_FLANK_ALIGNED)
            & (E.flank3_aligned >= MIN_FLANK_ALIGNED))
    one = ((E.flank5_aligned >= MIN_FLANK_ALIGNED)
           | (E.flank3_aligned >= MIN_FLANK_ALIGNED))
    full = (E.core_bases >= MIN_CORE_BASES) & (E.register > 0)
    E['tier'] = np.where(full & both, 1,
                         np.where(full & one, 2, 3))
    print(f'  TIER 1  core >={MIN_CORE_BASES} + register + BOTH flanks '
          f'>={MIN_FLANK_ALIGNED}bp : {int((E.tier == 1).sum())}')
    print(f'  TIER 2  core >={MIN_CORE_BASES} + register + ONE flank        : '
          f'{int((E.tier == 2).sum())}')
    print(f'  TIER 3  partial core or no register (REPORTED, NOT COUNTED): '
          f'{int((E.tier == 3).sum())}')
    n1, n2 = int((E.tier == 1).sum()), int((E.tier == 2).sum())
    if n1 == 0 and n2 > 0:
        print()
        print('  Tier 1 is empty while tier 2 is not. Both flanks can only be')
        print('  aligned if a read spans the 34-N cassette, so an empty tier 1')
        print('  means the aligner stopped at the N-run boundary on every')
        print('  read. Tier 2 is then the top tier available and this zero is')
        print('  not a quality problem.')
        print('  NOTE: the 2026-11-04 run on the real data produced 3 tier-1')
        print('  records, so STAR DOES align through the N run. An earlier')
        print('  version of this comment claimed tier 1 was structurally')
        print('  impossible. That was wrong.')
    print()
    print('  Tier 1 is trusted even when the core is absent from the catalog,')
    print('  because the stock is ~10,000 barcodes and the catalog holds a')
    print('  plasma-derived fraction of what each animal carries.')
    print('  Tier 3 exists so partials are visible, never to be counted: at')
    print('  6 visible bases a random fragment matches a real catalog core')
    print('  41.7% of the time.')

    # WHICH MECHANISM ACTUALLY RECOVERED THE COUNTABLE CORES?
    # This decides a sentence in the methods. Soft-clip projection was designed
    # as the crux of this stage, on the assumption STAR would never align
    # through 34 N. If the countable cores all come out fully aligned instead,
    # the projection is insurance rather than the mechanism, and the methods
    # should say so. The aggregate 'aligned vs projected' count in CELL 5
    # cannot answer this, because it pools partials with complete cores.
    cc = E[E.tier <= 2]
    if len(cc):
        full_al = int((cc.core_aligned == CORE_LEN).sum())
        part_al = int(((cc.core_aligned > 0) & (cc.core_aligned < CORE_LEN)).sum())
        no_al = int((cc.core_aligned == 0).sum())
        print(f'\n  HOW THE {len(cc)} COUNTABLE CORE(S) WERE RECOVERED:')
        print(f'    all {CORE_LEN} bases ALIGNED through the N run : {full_al}')
        print(f'    part aligned, part projected from a clip  : {part_al}')
        print(f'    entirely PROJECTED from a soft clip       : {no_al}')
        if no_al == 0 and full_al == len(cc):
            print('    Every countable core came from a read the aligner took')
            print('    straight through the 34-N cassette. The soft-clip')
            print('    projection recovered no complete core, so it is')
            print('    insurance here, not the mechanism. Say that in methods')
            print('    rather than describing projection as load-bearing.')
        elif no_al:
            print(f'    {no_al} core(s) rest ENTIRELY on the projection, which')
            print('    assumes no indel inside the clipped stretch. Those are')
            print('    the ones to check by hand before publishing.')
    save_csv(E, 'core_extractions.csv')


# %% ============================================================
# CELL 8 — MATES TO MOLECULES, THEN UMI DEDUP
# ============================================================
hr('CELL 8 — molecules')

countable = E[(E.tier <= 2)].copy() if not E.empty else E.copy()
print(f'  countable records (tier 1 or 2): {len(countable)}')

if countable.empty:
    MOL = pd.DataFrame(columns=['library', 'CB', 'UB', 'core', 'tier',
                                'n_records', 'n_pairs'])
    print('  nothing countable; no molecules to form.')
else:
    # 1. mates of one pair are one observation of one molecule
    pair = (countable.sort_values('tier')
            .groupby(['library', 'qname', 'core'], as_index=False)
            .agg(tier=('tier', 'min'), CB=('CB', 'first'), UB=('UB', 'first'),
                 register=('register', 'first'), n_mates=('mate', 'nunique')))
    print(f'  after collapsing mates on query name: {len(pair)} pair-level '
          f'observations')
    print(f'    pairs contributing BOTH mates: '
          f'{int((pair.n_mates == 2).sum())}')

    # DID R1 EARN ITS PLACE? The argument for reading R1 is that it is a
    # second look at the same molecule from the other end, so it can see a
    # cassette R2 missed. That is a claim about this dataset, and it is only
    # supported if some molecule was seen by R1 ALONE. Measure it rather than
    # asserting it: a pair where only mate 1 is countable is a molecule that
    # an R2-only analysis would have lost outright.
    bym = (countable.groupby(['library', 'qname'])['mate']
           .agg(lambda s: set(s.dropna())))
    r1_only = int(sum(1 for s in bym if s == {1}))
    r2_only = int(sum(1 for s in bym if s == {2}))
    print(f'    seen by R1 only: {r1_only}   by R2 only: {r2_only}   '
          f'by both: {int(sum(1 for s in bym if s == {1, 2}))}')
    if r1_only:
        print(f'    {r1_only} molecule(s) would have been LOST by an R2-only')
        print('    analysis. That is the measured justification for reading')
        print('    both mates, and it belongs in the methods.')
    else:
        print('    No molecule was seen by R1 alone, so on this data reading')
        print('    R1 cost nothing and added nothing. It stays in because it')
        print('    is free and the next run may differ, but do not claim it')
        print('    recovered anything here.')

    # 2. same CB + same UMI is one molecule amplified
    has_umi = (pair.CB != '') & (pair.UB != '')
    print(f'    with both CB and UB tags: {int(has_umi.sum())} of {len(pair)}')
    if (~has_umi).any():
        print(f'    WARN {int((~has_umi).sum())} observation(s) lack CB or UB '
              'and cannot be deduped; kept and flagged.')
    MOL = (pair.groupby(['library', 'CB', 'UB', 'core'], as_index=False)
           .agg(tier=('tier', 'min'), register=('register', 'first'),
                n_pairs=('qname', 'nunique')))
    print(f'  after (CB, UMI) dedup: {len(MOL)} MOLECULES')
    print('  Repeats are only evidence of clonal expansion if the molecules')
    print('  are independent. Same CB plus same UMI is one molecule amplified.')

    # 3. truncation collapse, prefix AND suffix
    # Suffix as well as prefix because the two registers are reverse
    # complements of each other, so a truncated core can be compatible from
    # either end depending on which register it came from.
    cores = sorted(MOL.core.unique(), key=len, reverse=True)
    canon = {}
    for c in cores:
        hit = next((f for f in canon.values()
                    if len(f) > len(c) and (f.startswith(c) or f.endswith(c))),
                   None)
        canon[c] = hit if hit else c
    n_coll = sum(1 for k, v in canon.items() if k != v)
    if n_coll:
        print(f'\n  truncation collapse merged {n_coll} core(s) into longer '
              'compatible ones:')
        for k, v in canon.items():
            if k != v:
                print(f'    {k} -> {v}')
        MOL['core'] = MOL.core.map(canon)
        MOL = (MOL.groupby(['library', 'CB', 'UB', 'core'], as_index=False)
               .agg(tier=('tier', 'min'), register=('register', 'first'),
                    n_pairs=('n_pairs', 'sum')))
        print(f'  molecules after collapse: {len(MOL)}')
    else:
        print('\n  no truncated cores to collapse (expected, since only '
              'complete cores are countable)')


# %% ============================================================
# CELL 9 — LINK TO CELLS AND ANIMALS
# ============================================================
hr('CELL 9 — cell and animal linkage')

OBS = None
if not os.path.exists(ANNOT_H5AD):
    print(f'  WARN  {ANNOT_H5AD} not found; cells cannot be linked.')
else:
    import anndata as ad
    a = ad.read_h5ad(ANNOT_H5AD, backed='r')
    keep = [c for c in ('library', 'animal_id', 'timepoint', 'tissue',
                        'tier1_celltype', 'tier2_subtype', 'qc_exclude',
                        'shiv_umi_cellranger', 'shiv_call', 'has_tcr')
            if c in a.obs.columns]
    OBS = a.obs[keep].copy()
    a.file.close()
    print(f'  loaded obs for {len(OBS):,} cells, columns {keep}')
    # 04 wrote GEX barcodes as '<bare>-1-<lib>'; the BAM CB tag is '<bare>-1'
    print("  join key: CB + '-' + library, matching 04's suffix convention")

if OBS is not None and not MOL.empty:
    MOL['cell'] = MOL.CB + '-' + MOL.library
    hit = MOL.cell.isin(OBS.index)
    print(f'  molecules whose cell is in the annotated object: '
          f'{int(hit.sum())} of {len(MOL)}')
    if (~hit).any():
        print(f'  {int((~hit).sum())} did not match. Those are barcodes that')
        print('  never became a called cell, or were removed by QC. They still')
        print('  tell us the core exists; they just cannot be attributed.')
        # the reverse CB orientation has an unexplained history; test it
        rc_cell = MOL.CB.map(lambda s: revcomp(s.split('-')[0]) + '-1') \
            + '-' + MOL.library
        n_rc = int(rc_cell.isin(OBS.index).sum())
        print(f'  sanity check, reverse-complemented CB matches: {n_rc}')
        if n_rc > int(hit.sum()):
            print('  WARN  the reverse-complemented CB matches MORE cells than')
            print('        the direct one. That is the 458-reads-to-2-cells')
            print('        history resurfacing. Stop and investigate.')
    # EVERY linkage column gets filled, not just animal_id. pandas groupby
    # drops rows with a NaN key by default, so leaving timepoint or tissue as
    # NaN silently deleted unlinked molecules from the inventory table in the
    # 2026-11-04 run: 3 barcodes were found and 1 was printed. The molecules
    # are real and the core is the finding; only the attribution is missing.
    LINK_COLS = ('animal_id', 'timepoint', 'tissue', 'tier1_celltype',
                 'tier2_subtype')
    for c in LINK_COLS:
        if c in OBS.columns:
            MOL[c] = MOL.cell.map(OBS[c]).astype(object).fillna('unlinked')
    if 'qc_exclude' in OBS.columns:
        MOL['qc_exclude'] = MOL.cell.map(OBS['qc_exclude']).astype(object)

    # ---- ATTRIBUTION, IN THREE TIERS OF EVIDENCE --------------------------
    # A molecule whose cell barcode never became a called cell is not therefore
    # unattributable. Three different things were conflated under 'unlinked' in
    # the first run, and they are not equally strong, so each gets its own
    # provenance and they are never mixed inside one column.
    #
    #   cell     the cell barcode is in the annotated object. Strongest.
    #   library  the library itself contains exactly one animal, so any
    #            molecule from it belongs to that animal regardless of whether
    #            its cell survived QC. As strong as the demultiplexing.
    #   catalog  the core is in the plasma catalog, which is animal-private
    #            77-90% of the time, so the catalog names an animal. WEAKEST,
    #            and it is kept out of the cross-check (see below).
    #
    # timepoint and tissue are NOT in this hierarchy at all, because they are
    # properties of the LIBRARY and are therefore always known. Printing them
    # as 'unlinked' threw away information LIB_NAME already holds.
    print('\n  ATTRIBUTION')
    # PLACEHOLDERS ARE NOT VALUES. 05 writes the string 'Unknown' for cells
    # that failed demultiplexing, and dropna() does not remove a string. In
    # the 2026-11-04 run that made library H's timepoints {Necropsy, Unknown},
    # length 2, so the fill declined and two barcodes from a library plainly
    # labelled Necropsy came out with timepoint 'unknown'. Strip the
    # placeholders first, then test uniqueness.
    PLACEHOLDERS = {'unknown', 'unlinked', 'none', 'nan', 'na', 'n/a', ''}

    def real_values(s):
        v = s.dropna().astype(str)
        return sorted(set(v[~v.str.lower().isin(PLACEHOLDERS)]))

    lib_tp, lib_tis, lib_animal = {}, {}, {}
    if {'library', 'timepoint', 'tissue'} <= set(OBS.columns):
        for lib, g in OBS.groupby('library', observed=True):
            tp, tis = real_values(g.timepoint), real_values(g.tissue)
            if len(tp) == 1:
                lib_tp[lib] = tp[0]
            if len(tis) == 1:
                lib_tis[lib] = tis[0]
    if {'library', 'animal_id'} <= set(OBS.columns):
        print('  distinct animals per library, read from the annotated object.')
        print('  Placeholder values are excluded from the count, since a cell')
        print('  that failed demultiplexing is not an extra animal.')
        for lib in sorted(OBS.library.unique()):
            g = OBS.loc[OBS.library == lib]
            animals = real_values(g.animal_id)
            n_ph = int(len(g.animal_id.dropna()) - g.animal_id.dropna().astype(str)
                       .str.lower().isin(PLACEHOLDERS).eq(False).sum())
            tag = ('  <- single animal, so attributable by library'
                   if len(animals) == 1 else '')
            extra = f'  (+{n_ph:,} placeholder cells)' if n_ph else ''
            print(f'    {lib}  {len(animals)} animal(s): '
                  f'{";".join(animals)}{extra}{tag}')
            if len(animals) == 1:
                lib_animal[lib] = animals[0]
        if not lib_animal:
            print('  No library is single-animal, so a molecule without a')
            print('  called cell cannot be attributed by library.')

    # timepoint and tissue: always from the library
    n_fix = 0
    for col, mapping in (('timepoint', lib_tp), ('tissue', lib_tis)):
        if col in MOL.columns and mapping:
            miss = MOL[col].isin(['unlinked']) | MOL[col].isna()
            n_fix += int((miss & MOL.library.isin(mapping)).sum())
            MOL.loc[miss, col] = MOL.loc[miss, 'library'].map(mapping)
            MOL[col] = MOL[col].fillna('unknown')
    if n_fix:
        print(f'\n  filled timepoint/tissue from the library on {n_fix} '
              'field(s).')
        print('  These are library properties, never cell properties, so they')
        print('  are known whether or not a cell survived QC.')
    # A library that contributed molecules but has NO cells in the annotated
    # object cannot be resolved by any of the three tiers, because all three
    # read from obs. That is the only way timepoint or tissue can stay
    # 'unknown', so name the library rather than leaving a bare placeholder.
    if 'library' in OBS.columns:
        absent = sorted(set(MOL.library) - set(OBS.library.unique()))
        if absent:
            print(f'\n  WARN  librar(ies) {absent} contributed molecules but '
                  'have no cells')
            print('        in the annotated object, so their timepoint, tissue')
            print('        and animal cannot be filled from it. Check that 04')
            print('        and 05 actually ingested them.')

    # animal: cell, then library, then catalog (catalog applied in CELL 10,
    # after the matching runs, so the column is created there)
    MOL['animal_source'] = np.where(
        MOL.get('animal_id', 'unlinked') != 'unlinked', 'cell', 'unattributed')
    if lib_animal:
        pick = (MOL.animal_source == 'unattributed') & MOL.library.isin(lib_animal)
        if pick.any():
            MOL.loc[pick, 'animal_id'] = MOL.loc[pick, 'library'].map(lib_animal)
            MOL.loc[pick, 'animal_source'] = 'library'
            print(f'\n  attributed {int(pick.sum())} molecule(s) by library.')

    print('\n  molecules by animal and evidence:')
    print(MOL.groupby(['animal_id', 'animal_source'], dropna=False)
          .size().rename('molecules').to_string())
    n_un = int((MOL.animal_source == 'unattributed').sum())
    if n_un:
        print(f'\n  {n_un} molecule(s) still unattributed after cell and')
        print('  library. CELL 10 tries the catalog on these, as the weakest')
        print('  tier. They are KEPT in every table either way, because the')
        print('  barcode existing is the finding and only the animal is')
        print('  missing.')


# %% ============================================================
# CELL 10 — CATALOG MATCH (exact, validation only)
# ============================================================
hr('CELL 10 — catalog validation')

if MOL.empty or not CAT_CORE:
    print('  nothing to validate.' if MOL.empty else
          '  no catalog loaded; skipping.')
else:
    allcore = set().union(*CAT_CORE.values())
    def match_animals(core):
        fwd = {a for a, v in CAT_CORE.items() if core in v}
        rev = {a for a, v in CAT_CORE.items() if revcomp(core) in v}
        return fwd, rev
    rows = []
    for r in MOL.itertuples():
        fwd, rev = match_animals(r.core)
        rows.append({'catalog_fwd': ';'.join(sorted(fwd)),
                     'catalog_rev': ';'.join(sorted(rev)),
                     'in_catalog': bool(fwd or rev)})
    MOL = pd.concat([MOL.reset_index(drop=True), pd.DataFrame(rows)], axis=1)
    n_in = int(MOL.in_catalog.sum())
    print(f'  molecules whose core is in the catalog (EXACT, either '
          f'orientation): {n_in} of {len(MOL)}')
    print(f'    matched as written      : {int((MOL.catalog_fwd != "").sum())}')
    print(f'    matched reverse-complemented: '
          f'{int((MOL.catalog_rev != "").sum())}')
    print('  Which orientation wins tells us the catalog convention; it is')
    print('  reported rather than assumed, and it gates nothing.')
    print()
    print('  A LOW MATCH RATE IS NOT EVIDENCE OF FALSE POSITIVES. The stock is')
    print('  ~10,000 barcodes and the catalog is a plasma readout of 113 to')
    print('  419 cores per animal, so most of what an animal carries was')
    print('  never in it.')

    # ---- THE CATALOG AS A LAST-RESORT ANIMAL, AND THE CIRCULARITY RULE ----
    # For a molecule with no called cell and no single-animal library, the
    # catalog is the only remaining route to an animal, and it is a real one:
    # 77-90% of catalog cores are animal-private.
    #
    # But it MUST NOT feed the cross-check below. Using the catalog to assign
    # the animal and then testing whether the catalog agrees with the animal
    # is circular and would manufacture a perfect score. So the catalog animal
    # goes in its own column, animal_source records it as 'catalog', and the
    # cross-check runs ONLY on animal_source == 'cell'.
    MOL['catalog_animal'] = [
        (r.catalog_fwd or r.catalog_rev) if (r.catalog_fwd or r.catalog_rev)
        else '' for r in MOL.itertuples()]
    if 'animal_source' in MOL.columns:
        cand = ((MOL.animal_source == 'unattributed')
                & (MOL.catalog_animal != '')
                & (~MOL.catalog_animal.str.contains(';')))   # unambiguous only
        if cand.any():
            MOL.loc[cand, 'animal_id'] = MOL.loc[cand, 'catalog_animal']
            MOL.loc[cand, 'animal_source'] = 'catalog'
            print(f'\n  attributed {int(cand.sum())} molecule(s) by CATALOG '
                  'alone, where the core maps to exactly one animal.')
            print('  This is the WEAKEST tier. It says the core is on that')
            print('  animal\'s plasma list, not that this molecule came from')
            print('  that animal, and it is excluded from the cross-check')
            print('  below to keep that test from being circular.')
        amb = ((MOL.animal_source == 'unattributed')
               & MOL.catalog_animal.str.contains(';'))
        if amb.any():
            print(f'\n  {int(amb.sum())} molecule(s) match a core shared by '
                  'several animals; left unattributed rather than guessed.')

    if 'animal_id' in MOL.columns:
        # 'cell' only. See the circularity rule above.
        chk = MOL[MOL.in_catalog & (MOL.animal_source == 'cell')].copy()
        if len(chk):
            chk['catalog_agrees'] = [
                r.animal_id in (set(r.catalog_fwd.split(';'))
                                | set(r.catalog_rev.split(';')))
                for r in chk.itertuples()]
            n_ok = int(chk.catalog_agrees.sum())
            print(f'\n  ANIMAL CROSS-CHECK on {len(chk)} molecule(s)')
            print('    (cell-attributed only, so the test is not circular):')
            print(f'    catalog animal agrees with demultiplexed animal: '
                  f'{n_ok} of {len(chk)}')
            print(f'    null expectation if matching were random: '
                  f'~{len(chk) / max(len(CAT_CORE), 1):.1f}')
            print('    This is a SANITY CHECK, not a statistical test. With a')
            print('    handful of molecules there is no power. A core matching')
            print('    the WRONG animal is a genuine red flag; agreement is')
            print('    reassurance and nothing more.')
            bad = chk[~chk.catalog_agrees]
            if len(bad):
                print(f'\n    WRONG-ANIMAL MATCHES ({len(bad)}). Read these;')
                print('    agreement is weak evidence but disagreement is')
                print('    strong. Either the cell was demultiplexed into the')
                print('    wrong animal, or the core is a false positive that')
                print('    happens to be a real barcode from another animal.')
                cols = [c for c in ('library', 'core', 'tier', 'animal_id',
                                    'catalog_fwd', 'catalog_rev')
                        if c in bad.columns]
                print(bad[cols].to_string(index=False))
            save_csv(chk, 'animal_cross_check.csv')
        else:
            print('\n  no catalog-matching molecules with a linked animal; '
                  'cross-check not possible.')


# %% ============================================================
# CELL 11 — THE NEGATIVE FLOOR
# ============================================================
hr('CELL 11 — negative floor')

print(f'  negative libraries (pre-infection or zero viral reads): '
      f'{sorted(NEGATIVE_LIBS)}')
print(f'  control animals (uninfected, absent from the catalog): '
      f'{CONTROL_ANIMALS}')
print('  Any core called here is a false positive by construction, so this is')
print('  the empirical error rate the tiers get judged against.')
print()

if MOL.empty:
    print('  no molecules; the floor is trivially zero.')
else:
    MOL['is_negative_lib'] = MOL.library.isin(NEGATIVE_LIBS)
    MOL['is_control_animal'] = MOL.get(
        'animal_id', pd.Series('unlinked', index=MOL.index)
    ).isin(CONTROL_ANIMALS)
    MOL['is_false_positive'] = MOL.is_negative_lib | MOL.is_control_animal
    rows = []
    for lab, m in (('negative library', MOL.is_negative_lib),
                   ('control animal', MOL.is_control_animal),
                   ('EITHER (false positive)', MOL.is_false_positive),
                   ('infected lib + infected animal', ~MOL.is_false_positive)):
        rows.append({'class': lab, 'molecules': int(m.sum()),
                     'distinct_cores': MOL.loc[m, 'core'].nunique()})
    floor = pd.DataFrame(rows)
    print(floor.to_string(index=False))
    save_csv(floor, 'negative_floor.csv')
    nfp = int(MOL.is_false_positive.sum())
    if nfp:
        print(f'\n  {nfp} molecule(s) in a negative. Those cores and their')
        print('  tiers are listed below; treat the tier that produced them as')
        print('  unreliable for the real calls too.')
        print(MOL.loc[MOL.is_false_positive,
                      ['library', 'core', 'tier', 'animal_id']
                      if 'animal_id' in MOL.columns else
                      ['library', 'core', 'tier']].to_string(index=False))
    else:
        print('\n  ZERO molecules in any negative. That is the good outcome:')
        print('  no measurable false-positive floor at these tiers.')


# %% ============================================================
# CELL 12 — INVENTORY AND THE LONGITUDINAL QUESTION
# ============================================================
hr('CELL 12 — per-animal inventory')

real = MOL[~MOL.is_false_positive].copy() if not MOL.empty else MOL.copy()
print(f'  molecules passing the negative filter: {len(real)}')

if real.empty:
    print('\n  No barcode molecules survived. That is the result, and it is')
    print('  consistent with the depth: 05 found 59 cells carrying any viral')
    print('  UMI at all, median 1 molecule, and the cassette is 34bp of a')
    print('  10,296bp contig.')
    inv = pd.DataFrame()
else:
    save_csv(real, 'barcode_molecules.csv')
    keys = [k for k in ('animal_id', 'timepoint', 'tissue') if k in real.columns]
    # dropna=False as belt and braces. Cell 9 now fills these, but a groupby
    # that silently deletes rows is not something to rely on a fill for.
    def joinset(s):
        return ';'.join(sorted({str(x) for x in s.dropna()})) or 'unlinked'

    aggs = dict(molecules=('UB', 'nunique'), tier=('tier', 'min'),
                libraries=('library', lambda s: ';'.join(sorted(set(s)))),
                in_catalog=('in_catalog', 'max'))
    if 'animal_source' in real.columns:
        aggs['animal_evidence'] = ('animal_source', joinset)
    # The cell type call. Available ONLY for cell-attributed molecules, since
    # a molecule with no called cell has no cell type by definition. It is
    # descriptive, not evidence: see the caveat printed under the table.
    for c in ('tier1_celltype', 'tier2_subtype'):
        if c in real.columns:
            aggs[c] = (c, joinset)
    inv = (real.groupby(keys + ['core'], as_index=False, dropna=False)
           .agg(**aggs))
    print('\n  barcode inventory:')
    print(inv.to_string(index=False))
    if 'tier1_celltype' in inv.columns:
        n_ct = int((inv.tier1_celltype != 'unlinked').sum())
        print(f'\n  CELL TYPE IS DESCRIPTIVE, NOT EVIDENCE OF INFECTION.')
        print(f'  {n_ct} of {len(inv)} inventory row(s) carry a cell type, and')
        print('  only cell-attributed molecules can have one. 05b established')
        print('  that per-cell viral attribution below shiv_confident is not')
        print('  publishable: the strongest B cell carried less virus than the')
        print('  strongest erythrocyte. So a barcode appearing in a T cell')
        print('  says the barcode is in that ANIMAL and that this ambient')
        print('  molecule happened to land in a droplet that also held a T')
        print('  cell. Do not build a tropism claim on this column.')
    # The inventory must account for every molecule that got here. If it does
    # not, a grouping key ate some, and the table is wrong rather than short.
    if int(inv.molecules.sum()) != real.UB.nunique():
        print(f'  WARN  inventory holds {int(inv.molecules.sum())} molecules '
              f'but {real.UB.nunique()} arrived. A grouping key dropped rows.')
    if inv.core.nunique() != real.core.nunique():
        print(f'  WARN  inventory shows {inv.core.nunique()} distinct cores '
              f'but {real.core.nunique()} arrived.')
    save_csv(inv, 'barcode_inventory.csv')

    if 'animal_id' in real.columns:
        dkeys = [k for k in keys if k != 'tissue']
        dagg = dict(distinct_cores=('core', 'nunique'),
                    molecules=('UB', 'nunique'))
        if 'animal_source' in real.columns:
            # Split the diversity count by evidence tier. A per-animal
            # diversity number that silently pools cell-attributed with
            # catalog-attributed molecules is not one number, it is two.
            dagg['by_cell'] = ('animal_source', lambda s: int((s == 'cell').sum()))
            dagg['by_library'] = ('animal_source', lambda s: int((s == 'library').sum()))
            dagg['by_catalog'] = ('animal_source', lambda s: int((s == 'catalog').sum()))
        div = (real.groupby(dkeys, dropna=False).agg(**dagg).reset_index())
        print('\n  diversity per animal and timepoint:')
        print(div.to_string(index=False))
        if 'by_catalog' in div.columns and int(div.by_catalog.sum()):
            print('  by_catalog counts rest on the weakest attribution tier.')
            print('  Quote a diversity figure from by_cell plus by_library, and')
            print('  report the catalog-only ones separately.')
        save_csv(div, 'barcode_diversity.csv')

    # --- the longitudinal question -------------------------------------
    print('\n  THE LONGITUDINAL QUESTION')
    print('  Viral reads exist only at 21 DPI (D, E, F) and Necropsy (H),')
    print('  since A/B/C are pre-infection and G returned zero. So the only')
    print('  longitudinal statement available is whether a core appears at')
    print('  BOTH timepoints in the SAME animal.')
    if {'animal_id', 'timepoint'} <= set(real.columns):
        # TWO GUARDS, because this is the headline longitudinal claim and the
        # weakest evidence tier must not be able to create it.
        #
        # 1. 'unknown' and 'unlinked' are NOT timepoints. A core seen at
        #    21 DPI and again at an unknown timepoint has been seen at one
        #    known timepoint, not two, and counting the placeholder would
        #    invent persistence out of a missing value.
        # 2. A catalog-attributed molecule says only that the core is on some
        #    animal's plasma list. Letting that tier supply the second
        #    timepoint would mean the catalog, not the data, produced the
        #    persistence result. Those are reported separately as candidates.
        BAD_TP = {'unknown', 'unlinked', 'nan', 'None'}
        r = real.copy()
        r['tp_known'] = ~r.timepoint.astype(str).isin(BAD_TP)
        strong = (r.animal_source.isin(['cell', 'library'])
                  if 'animal_source' in r.columns
                  else pd.Series(True, index=r.index))

        def tps(df):
            out = (df[df.tp_known].groupby(['animal_id', 'core'])['timepoint']
                   .agg(lambda s: sorted({str(x) for x in s})).reset_index())
            out['n_timepoints'] = out.timepoint.map(len)
            return out

        pers = tps(r[strong])
        shared = pers[pers.n_timepoints > 1]
        n_drop = int((~r.tp_known).sum())
        if n_drop:
            print(f'\n  {n_drop} molecule(s) have no known timepoint and are')
            print('  excluded from this test. A placeholder is not a timepoint.')
        print(f'\n  cores at more than one KNOWN timepoint in one animal, on '
              f'cell or library evidence: {len(shared)}')
        if len(shared):
            print(shared.to_string(index=False))
            print('  Persistence across timepoints at this depth would be')
            print('  notable. Check the molecule counts before believing it,')
            print('  and confirm the two observations are separate molecules')
            print('  rather than one molecule counted twice.')
        else:
            print('  None. At this depth that is the expected answer and it is')
            print('  a finding, not a failure: there is not enough coverage to')
            print('  observe persistence even if it is occurring.')

        # the weaker set, reported and never merged into the line above
        if 'animal_source' in r.columns and (~strong).any():
            pers_all = tps(r)
            weak = pers_all[(pers_all.n_timepoints > 1)].merge(
                shared[['animal_id', 'core']], on=['animal_id', 'core'],
                how='left', indicator=True)
            weak = weak[weak._merge == 'left_only'].drop(columns=['_merge'])
            if len(weak):
                print(f'\n  {len(weak)} further candidate(s) appear at two')
                print('  timepoints only once CATALOG-attributed molecules are')
                print('  included. These are NOT persistence results: the')
                print('  catalog supplied the animal, so the catalog supplied')
                print('  the claim. Listed for follow-up, not for reporting.')
                print(weak.to_string(index=False))
                save_csv(weak, 'barcode_persistence_catalog_only.csv')
            pers = pers_all
        save_csv(pers, 'barcode_persistence.csv')

    # tissue observation that is free to make
    print('\n  FREE OBSERVATION, already in the 05 data: at necropsy the virus')
    print('  is in LN (library H, 34 viral cells) and absent from PBMC')
    print('  (library G, zero). Confounded by depth, but worth recording.')


# %% ============================================================
# CELL 13 — VERDICT
# ============================================================
hr('CELL 13 — verdict')

n_prox = len(prox) if 'prox' in dir() and not prox.empty else 0
n_full = int((E.core_bases == CORE_LEN).sum()) if not E.empty else 0
n_count = int((E.tier <= 2).sum()) if not E.empty else 0
n_mol = len(MOL) if not MOL.empty else 0
n_real = len(real) if not real.empty else 0
n_cores = real.core.nunique() if not real.empty else 0

funnel = pd.DataFrame([
    {'stage': 'records aligned to contig', 'n': len(reads)},
    {'stage': 'pairs (molecules pre-UMI)', 'n': reads.qname.nunique()},
    {'stage': 'records within 60bp of cassette', 'n': n_prox},
    {'stage': f'records with complete {CORE_LEN}bp core', 'n': n_full},
    {'stage': 'records in tier 1 or 2', 'n': n_count},
    {'stage': 'molecules after (CB, UMI) dedup', 'n': n_mol},
    {'stage': 'molecules passing negative filter', 'n': n_real},
    {'stage': 'DISTINCT BARCODES', 'n': n_cores},
])
print(funnel.to_string(index=False))
save_csv(funnel, 'detection_funnel.csv')

print()
if n_cores == 0:
    print('  VERDICT: no barcodes recoverable from this data.')
elif n_cores <= 5:
    print(f'  VERDICT: {n_cores} barcode(s) recovered. Enough to report as')
    print('  individual observations, NOT enough for a diversity or')
    print('  longitudinal analysis.')
else:
    print(f'  VERDICT: {n_cores} barcodes recovered, which is more than prior')
    print('  passes found. Re-read the negative floor and the tier breakdown')
    print('  before treating it as a diversity measurement.')

print()
print('  ON TARGETED ENRICHMENT, which is the decision this stage exists to')
print('  inform. The cassette is 34bp of a 10,296bp contig, so even with the')
print('  GTF fixed and both mates aligned, whole-transcriptome 5\' capture')
print('  puts a vanishing fraction of reads across the barcode. Nothing in')
print('  the analysis can create coverage that was never sequenced.')
if n_cores <= 5:
    print()
    print('  RECOMMENDATION: targeted VpxF1 amplicon enrichment on the same')
    print('  libraries is the only route to a barcode diversity or')
    print('  longitudinal result. This pipeline has now extracted everything')
    print('  the existing data contains, and the limit is the library, not')
    print('  the method.')

print()
print('  AND FOR THE RECORD, carried forward from 05b: per-cell viral')
print('  attribution below shiv_confident is not publishable. The ambient')
print('  rate puts virus in erythrocytes. A barcode recovered from a cell is')
print('  evidence that the barcode exists in that ANIMAL, not that the cell')
print('  was infected.')

hr('DONE — 06_viral_barcode.py')
print(f'  outputs: {OUT_DIR}')
print('    detection_funnel.csv        read it top to bottom first')
print('    cassette_geometry.csv       whether the core is readable at all')
print('    core_extractions.csv        every record, every tier')
print('    barcode_molecules.csv       after mate collapse and UMI dedup')
print('    barcode_inventory.csv       per animal, timepoint, tissue')
print('    barcode_persistence.csv     the longitudinal question')
print('    negative_floor.csv          the false-positive rate')
print('    animal_cross_check.csv      catalog animal vs demultiplexed animal')
print('\n  Next: 07_summary.py')
