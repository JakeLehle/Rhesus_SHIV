#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
00 — SHIV GTF feature-multiplicity diagnostic (READ-ONLY)
=========================================================

Purpose
-------
Size the GTF fix BEFORE committing a five-day Cell Ranger rerun.

`shiv_mapping_diagnostic.py` established in September that the SHIV undercount
is an annotation failure and not an alignment failure: reads align to the virus
correctly and then fail to land in a counted feature. This script measures how
much of that loss the GTF actually causes, and how much a single collapsed
provirus feature would recover, using the GTF that built the reference now in
use plus the viral reads that already exist on disk.

Cell Ranger's counting rule is the whole basis of this script:

  * a read overlapping exactly ONE gene feature  -> counted
  * a read overlapping TWO OR MORE gene features -> ambiguous, DISCARDED
  * a read overlapping NO gene feature           -> intergenic, DISCARDED

So every aligned viral read is classified into one of those three buckets
against the current GTF, and then re-classified against a hypothetical single
collapsed feature spanning the whole contig. The ratio between the two
"counted" numbers is the measured version of the projected ~4x.

Mechanism this is expected to confirm
-------------------------------------
The GTF converter in MKREF_Mmul10_SHIVAD8EO.sh derives gene spans with
`parse_location`, which takes min() and max() of every number inside a
GenBank join(). *tat* and *rev* are multiply spliced in the SIVmac239
backbone, so collapsing their join() to outermost coordinates produces spans
of roughly 3 kb that swallow *vpu*, *env* and *nef* whole. `flush_feature`
then widens each gene name again across all its CDS entries. The prediction
is therefore that loss is concentrated in the 3' half of the genome (high
multiplicity) plus the LTRs (zero multiplicity), and NOT spread evenly.

If the measured recoverable fraction comes back near 1.0, the GTF was never
the bottleneck and the rerun is not justified on these grounds. That is a real
possible outcome of this script and the reason it runs first.

Dropped / rejected approaches, recorded so they are not relitigated
-------------------------------------------------------------------
* Per-base depth (`samtools depth`) was rejected as the primary metric.
  Cell Ranger's ambiguity rule applies per READ, not per base, so a per-base
  multiplicity track cannot tell you how many reads are discarded. The
  per-base track is still computed, but only for the figure.
* Nesting the 34 bp cassette as its own GTF feature inside a collapsed
  provirus feature was rejected on first principles. Reads are 150 bp and the
  cassette is 34 bp, so every cassette-covering read also covers flanking
  sequence; if the flank belongs to another feature the read overlaps two
  features and is discarded. That would make cassette reads LESS countable,
  not more. See the collapsed-feature recommendation printed by Cell 6.
* `include-introns` does not rescue the LTRs. With this GTF the gene,
  transcript and exon entries share identical spans, so there are no introns,
  and the LTRs fall outside every gene span regardless of the setting.

Inputs (all read-only)
----------------------
  SHIV GTF          ${REFDIR}/SHIVAD8EO.gtf      (the SHIV-only GTF)
  SHIV FASTA        ${REFDIR}/SHIVAD8EO.fna      (for contig length only)
  Viral alignments  whichever of these exist, each reported separately:
                      Cell Ranger per-library BAMs under Analysis/*/outs/
                      STARsolo R1 BAMs  shiv_starsolo_r1/solo/Forward/
                      STARsolo R2 BAMs  shiv_starsolo_bc/solo/V2allN/

Outputs (all under OUT_DIR, nothing else is written)
----------------------------------------------------
  gtf_gene_spans.csv              every gene feature with its span and width
  gtf_pairwise_overlap.csv        gene x gene overlap in bp
  gtf_base_multiplicity.csv       per-base feature count along the contig
  read_classification.csv         per source x library read bucket counts
  read_classification_by_gene.csv counted reads attributed per gene
  recovery_projection.csv         counted now vs counted under collapse
  F00_gtf_multiplicity.pdf/.png   multiplicity track + read depth

Environment
-----------
  conda activate sc_pre
  export PATH="${CONDA_PREFIX}/bin:${PATH}"      # system samtools is broken

Author: Jake Lehle
Date:   October 2026
"""

# %% ============================================================
# CELL 1 — CONFIGURATION (every tunable lives here)
# ============================================================
import os
import re
import glob
import shutil
import subprocess
from collections import defaultdict

import numpy as np
import pandas as pd

os.environ.setdefault('MPLBACKEND', 'Agg')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# --- shared project config, with standalone fallback so this still opens
#     and runs in Spyder without the pipeline/config directory on sys.path ---
try:
    from shiv_config import WORKDIR, REFDIR, LIBS, LIB_NAME        # noqa: F401
except Exception:                                                   # noqa: BLE001
    WORKDIR = '/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV'
    REFDIR = '/master/jlehle/WORKING/SC/REF'
    LIBS = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
    LIB_NAME = {
        'A': 'Pre_PBMC',
        'B': 'Pre_LN',
        'C': 'Pre_PBMC_IndianMixed',
        'D': '21DPI_Pre_LN_Mixed',
        'E': '21DPI_PBMC',
        'F': '21DPI_NonInf_LN_Mixed',
        'G': 'Necropsy_PBMC',
        'H': 'Necropsy_LN',
    }

# --- inputs ---
SHIV_GTF = os.path.join(REFDIR, 'SHIVAD8EO.gtf')
SHIV_FASTA = os.path.join(REFDIR, 'SHIVAD8EO.fna')

# Contig name per alignment source. The STARsolo barcode references use a
# different contig name than the Cell Ranger reference, so each source
# declares its own rather than assuming one global name.
# The globs are deliberately permissive and the contig is confirmed from each
# BAM header, because the 2026-10-30 run found no STARsolo BAMs at all with
# tighter patterns while 00b's solo/*/*.bam found 10 R1 and 32 R2 BAMs.
# Mate/reference subdirectories vary, so do not pin them.
SOURCES = {
    'cellranger': {
        'glob': os.path.join(WORKDIR, 'Analysis', 'SHIV_{libname}', 'outs', '**', '*.bam'),
        'contig': 'SHIVAD8EO',
    },
    'starsolo_r1': {
        'glob': os.path.join(WORKDIR, 'shiv_starsolo_r1', 'solo', '**', '*.bam'),
        'contig': 'SHIVAD8EO_V2allN',
    },
    'starsolo_r2': {
        'glob': os.path.join(WORKDIR, 'shiv_starsolo_bc', 'solo', '**', '*.bam'),
        'contig': 'SHIVAD8EO_V2allN',
    },
}

# Classify EVERY candidate BAM per library rather than only the largest.
# The 2026-10-30 run kept sample_alignments.bam and discarded
# unassigned_alignments.bam, which holds reads from droplets not called as
# cells. That is exactly where ambient viral material and the QC-dropped
# positives sit, so discarding it undercounts the available viral reads.
# Each BAM is reported on its own row, keyed by source and basename.
CLASSIFY_ALL_CANDIDATE_BAMS = True

# --- behaviour ---
# 0x904 = unmapped + secondary + supplementary. Primary alignments only, which
# is what Cell Ranger would be counting.
SAM_EXCLUDE_FLAG = '0x904'
MIN_ALIGNED_BASES = 1       # a read must have at least this many M/=/X bases
REQUIRE_CB = False          # True restricts to reads carrying a CB tag

# --- output ---
OUT_DIR = os.path.join(WORKDIR, 'pipeline_results', '00_gtf_multiplicity')

# --- figure conventions ---
FIG_DPI = 300
FS_TITLE, FS_LABEL, FS_TICK, FS_LEGEND = 34, 30, 28, 28
HEX_COUNTED = '#2E7D32'     # multiplicity 1, countable
HEX_AMBIG = '#C62828'       # multiplicity >= 2, discarded as ambiguous
HEX_INTERGENIC = '#9E9E9E'  # multiplicity 0, discarded as intergenic
HEX_DEPTH = '#1565C0'       # observed read depth
HEX_GRID = '#E0E0E0'

os.makedirs(OUT_DIR, exist_ok=True)
pd.set_option('display.width', 220)
pd.set_option('display.max_columns', 80)

CIGAR_RE = re.compile(r'(\d+)([MIDNSHP=X])')


def hr(t):
    print('\n' + '=' * 94)
    print(t)
    print('=' * 94)


def save_csv(df, name, index=True):
    p = os.path.join(OUT_DIR, name)
    df.to_csv(p, index=index, lineterminator='\n')
    print(f'  wrote {name}  ({df.shape[0]} rows)')


def savefig(fig, stem):
    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(OUT_DIR, f'{stem}.{ext}'),
                    dpi=FIG_DPI, bbox_inches='tight')
    plt.close(fig)
    print(f'  wrote {stem}.pdf / {stem}.png')


SAMTOOLS = shutil.which('samtools')
if SAMTOOLS and SAMTOOLS.startswith('/usr/local/bin'):
    print('WARN  samtools resolves to /usr/local/bin, which is the broken '
          'system build.\n      conda activate sc_pre and '
          'export PATH="${CONDA_PREFIX}/bin:${PATH}"')


# %% ============================================================
# CELL 2 — PARSE THE SHIV GTF INTO GENE SPANS
# ============================================================
hr('CELL 2 — gene spans from the GTF that built the reference in use')

if not os.path.exists(SHIV_GTF):
    raise SystemExit(f'GTF not found: {SHIV_GTF}\n'
                     'Point SHIV_GTF at the SHIV-only GTF used by mkref.')

# Contig length from the FASTA. Falls back to the GTF maximum so the script
# still runs if the FASTA is not beside it.
contig_len = None
if os.path.exists(SHIV_FASTA):
    n = 0
    with open(SHIV_FASTA) as fh:
        for line in fh:
            if not line.startswith('>'):
                n += len(line.strip())
    contig_len = n
    print(f'  contig length from FASTA : {contig_len:,} bp')
else:
    print(f'  WARN  FASTA not found ({SHIV_FASTA}); '
          'falling back to GTF max coordinate')

gene_rows = []
with open(SHIV_GTF) as fh:
    for line in fh:
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'gene':
            continue
        m = re.search(r'gene_name "([^"]+)"', f[8])
        gname = m.group(1) if m else None
        m = re.search(r'gene_id "([^"]+)"', f[8])
        gid = m.group(1) if m else gname
        gene_rows.append({
            'gene_name': gname or gid,
            'gene_id': gid,
            'contig': f[0],
            # GTF is 1-based inclusive; convert to 0-based half-open
            'start0': int(f[3]) - 1,
            'end0': int(f[4]),
            'strand': f[6],
        })

if not gene_rows:
    raise SystemExit('No "gene" features found in the GTF. Nothing to measure.')

genes = pd.DataFrame(gene_rows)
genes['width'] = genes['end0'] - genes['start0']
genes = genes.sort_values('start0').reset_index(drop=True)

if contig_len is None:
    contig_len = int(genes['end0'].max())
    print(f'  contig length from GTF    : {contig_len:,} bp')

print(f'  gene features in GTF      : {len(genes)}')
print(f'  contigs named in GTF      : {sorted(genes["contig"].unique())}')
print()
print(genes[['gene_name', 'start0', 'end0', 'width', 'strand']].to_string(index=False))
save_csv(genes, 'gtf_gene_spans.csv', index=False)

# Widest features first. If tat/rev top this list with multi-kb spans, the
# join()-collapse mechanism in the header is confirmed from the GTF alone.
print('\n  widest features (the join()-collapse signature):')
print(genes.nlargest(min(5, len(genes)), 'width')
      [['gene_name', 'start0', 'end0', 'width']].to_string(index=False))


# %% ============================================================
# CELL 3 — PER-BASE FEATURE MULTIPLICITY AND PAIRWISE OVERLAP
# ============================================================
hr('CELL 3 — how much of the contig is countable, ambiguous or invisible')

mult = np.zeros(contig_len, dtype=np.int16)
for _, g in genes.iterrows():
    s, e = max(0, g['start0']), min(contig_len, g['end0'])
    if e > s:
        mult[s:e] += 1

n0 = int((mult == 0).sum())
n1 = int((mult == 1).sum())
n2 = int((mult >= 2).sum())
print(f'  bases in no feature   (intergenic, never counted) : {n0:6,}  '
      f'({100 * n0 / contig_len:5.1f}%)')
print(f'  bases in one feature  (countable)                 : {n1:6,}  '
      f'({100 * n1 / contig_len:5.1f}%)')
print(f'  bases in 2+ features  (ambiguous, discarded)      : {n2:6,}  '
      f'({100 * n2 / contig_len:5.1f}%)')
print(f'  max features over any single base                 : {int(mult.max())}')

base_mult = pd.DataFrame({'pos0': np.arange(contig_len), 'n_features': mult})
save_csv(base_mult, 'gtf_base_multiplicity.csv', index=False)

# Pairwise overlap, so it is visible which genes are swallowing which
names = genes['gene_name'].tolist()
ov = pd.DataFrame(0, index=names, columns=names, dtype=int)
for i, a in genes.iterrows():
    for j, b in genes.iterrows():
        if i >= j:
            continue
        o = max(0, min(a['end0'], b['end0']) - max(a['start0'], b['start0']))
        if o:
            ov.iat[i, j] = o
            ov.iat[j, i] = o
save_csv(ov, 'gtf_pairwise_overlap.csv')

pairs = []
for i, a in genes.iterrows():
    for j, b in genes.iterrows():
        if i >= j:
            continue
        o = ov.iat[i, j]
        if o:
            pairs.append((a['gene_name'], b['gene_name'], int(o)))
pairs.sort(key=lambda x: -x[2])
if pairs:
    print('\n  largest pairwise overlaps (bp):')
    for a, b, o in pairs[:12]:
        print(f'    {a:>8} x {b:<8} {o:6,}')
else:
    print('\n  no pairwise overlaps at all. The GTF is already disjoint, and '
          'ambiguity\n  is not the loss mechanism. Read Cell 5 before '
          'committing to the rerun.')


# %% ============================================================
# CELL 4 — DISCOVER ALIGNMENT SOURCES
# ============================================================
hr('CELL 4 — locate viral alignments on disk (nothing is modified)')


def contig_in_header(bam, contig):
    """True if `contig` appears as an @SQ SN in this BAM header."""
    try:
        p = subprocess.run([SAMTOOLS, 'view', '-H', bam],
                           capture_output=True, text=True, check=True)
    except Exception as e:                                          # noqa: BLE001
        print(f'    WARN  header read failed for {os.path.basename(bam)}: {e}')
        return False
    return any(line.startswith('@SQ') and f'SN:{contig}\t' in line + '\t'
               for line in p.stdout.splitlines())


LIB_FROM_NAME = re.compile(r'^([A-H])_')


def lib_from_path(path):
    """Library letter from a STARsolo BAM basename, which is
    {LIB}_{LIB_NAME}_Aligned.sortedByCoord.out.bam."""
    m = LIB_FROM_NAME.match(os.path.basename(path))
    return m.group(1) if m else '?'


found = []   # (source, lib, bam_path, contig)
if not SAMTOOLS:
    print('  ERROR samtools not on PATH. Cells 4-6 are skipped; Cells 2-3 '
          'results above\n        are still valid and already written.')
else:
    for src, spec in SOURCES.items():
        print(f'  {src}:')
        tmpl = spec['glob']
        # A pattern carrying {lib}/{libname} is expanded once per library. A
        # pattern without them is library-agnostic and must be globbed ONCE,
        # with the library read off each filename. Expanding such a pattern in
        # the per-library loop matches every BAM eight times over and inflates
        # every absolute count by 8x, which is what the 2026-10-30 run did.
        per_lib = ('{lib}' in tmpl) or ('{libname}' in tmpl)
        if per_lib:
            cand = []
            for lib in LIBS:
                pat = tmpl.format(lib=lib, libname=LIB_NAME[lib])
                for h in sorted(glob.glob(pat, recursive=True)):
                    cand.append((lib, h))
        else:
            cand = [(lib_from_path(h), h)
                    for h in sorted(glob.glob(tmpl, recursive=True))]

        # Cell Ranger multi writes several BAMs under outs/; keep only those
        # whose header actually carries the viral contig.
        cand = [(lib, h) for lib, h in cand
                if contig_in_header(h, spec['contig'])]
        # Belt and braces: never classify the same file twice in one source.
        seen_paths = set()
        cand = [(lib, h) for lib, h in cand
                if not (h in seen_paths or seen_paths.add(h))]

        if not CLASSIFY_ALL_CANDIDATE_BAMS and cand:
            by_lib = {}
            for lib, h in cand:
                if lib not in by_lib or os.path.getsize(h) > os.path.getsize(by_lib[lib]):
                    by_lib[lib] = h
            cand = sorted(by_lib.items())

        for lib, h in cand:
            found.append((src, lib, h, spec['contig']))
            print(f'    {lib}: {os.path.relpath(h, WORKDIR)}')
        if not cand:
            print('    none found (skipping this source)')

    if not found:
        print('\n  WARN  no BAMs carrying a SHIV contig were found. Check the '
              'SOURCES\n        globs in Cell 1 against your Analysis/ layout. '
              'Cells 5-6 will be empty.')


# %% ============================================================
# CELL 5 — CLASSIFY EVERY ALIGNED VIRAL READ BY CELL RANGER'S RULE
# ============================================================
hr("CELL 5 — per-read buckets: counted, ambiguous, intergenic")


def ref_blocks(pos0, cigar):
    """Reference intervals this alignment actually covers.

    Only M/=/X consume both read and reference, so only those are real
    coverage. D and N advance the reference without read bases and are NOT
    counted as coverage, which is what lets a deletion variant inside the
    cassette still be attributed correctly.
    """
    blocks, r = [], pos0
    for n, op in CIGAR_RE.findall(cigar):
        n = int(n)
        if op in 'M=X':
            blocks.append((r, r + n))
            r += n
        elif op in 'DN':
            r += n
        # I, S, H, P consume no reference
    return blocks


def genes_hit(blocks, spans):
    """Set of gene names any aligned block overlaps."""
    out = set()
    for bs, be in blocks:
        for gname, gs, ge in spans:
            if bs < ge and gs < be:
                out.add(gname)
    return out


spans = [(r['gene_name'], r['start0'], r['end0']) for _, r in genes.iterrows()]

cls_rows = []
gene_rows_counted = []
depth = np.zeros(contig_len, dtype=np.int64)

for src, lib, bam, contig in found:
    try:
        p = subprocess.run(
            [SAMTOOLS, 'view', '-F', SAM_EXCLUDE_FLAG, bam, contig],
            capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:                      # noqa: BLE001
        print(f'  WARN  {src}/{lib}: samtools view failed ({e}); skipping')
        continue

    n_tot = n_cnt = n_amb = n_int = n_skip = 0
    per_gene = defaultdict(int)

    for line in p.stdout.splitlines():
        f = line.split('\t')
        if len(f) < 11:
            continue
        if REQUIRE_CB and 'CB:Z:' not in line:
            continue
        pos0 = int(f[3]) - 1
        blocks = ref_blocks(pos0, f[5])
        aligned = sum(e - s for s, e in blocks)
        if aligned < MIN_ALIGNED_BASES:
            n_skip += 1
            continue
        n_tot += 1
        for s, e in blocks:
            s, e = max(0, s), min(contig_len, e)
            if e > s:
                depth[s:e] += 1
        hits = genes_hit(blocks, spans)
        if len(hits) == 1:
            n_cnt += 1
            per_gene[next(iter(hits))] += 1
        elif len(hits) == 0:
            n_int += 1
        else:
            n_amb += 1

    cls_rows.append({
        'source': src, 'lib': lib, 'library_name': LIB_NAME[lib],
        # Path, not basename: the R1 pass writes the same filename under both
        # solo/Forward and solo/Reverse, and keying on basename silently
        # merges the two strands into one row.
        'bam': os.path.relpath(bam, WORKDIR),
        'reads_total': n_tot,
        'counted_current_gtf': n_cnt,
        'ambiguous_multi_feature': n_amb,
        'intergenic_no_feature': n_int,
        'skipped_no_aligned_bases': n_skip,
        # Under one collapsed feature spanning the contig, every read with any
        # aligned base overlaps exactly one feature, so all of them count.
        'counted_collapsed_feature': n_tot,
    })
    for g, c in per_gene.items():
        gene_rows_counted.append({'source': src, 'lib': lib,
                                  'gene_name': g, 'reads_counted': c})
    if n_tot:
        tag = os.path.relpath(bam, WORKDIR)
        if len(tag) > 56:
            tag = '...' + tag[-53:]
        print(f'  {src}/{lib} {tag:<56} '
              f'{n_tot:7,} reads | counted {n_cnt:7,} '
              f'| ambiguous {n_amb:7,} | intergenic {n_int:7,}')

cls = pd.DataFrame(cls_rows)
if cls.empty:
    print('  no reads classified. Cell 6 is skipped.')
else:
    save_csv(cls, 'read_classification.csv', index=False)
    if gene_rows_counted:
        gdf = (pd.DataFrame(gene_rows_counted)
               .pivot_table(index='gene_name', columns='source',
                            values='reads_counted', aggfunc='sum')
               .fillna(0).astype(int))
        save_csv(gdf, 'read_classification_by_gene.csv')
        print('\n  counted reads attributed per gene:')
        print(gdf.to_string())


# %% ============================================================
# CELL 6 — RECOVERY PROJECTION AND THE COLLAPSED-FEATURE RECOMMENDATION
# ============================================================
hr('CELL 6 — what a single collapsed provirus feature would recover')

if cls.empty:
    print('  nothing to project.')
else:
    proj = (cls.groupby('source')[
        ['reads_total', 'counted_current_gtf', 'ambiguous_multi_feature',
         'intergenic_no_feature', 'counted_collapsed_feature']]
        .sum())
    proj['fold_recovery'] = np.where(
        proj['counted_current_gtf'] > 0,
        proj['counted_collapsed_feature'] / proj['counted_current_gtf'],
        np.nan)
    proj['pct_lost_ambiguous'] = (100 * proj['ambiguous_multi_feature']
                                  / proj['reads_total'].replace(0, np.nan))
    proj['pct_lost_intergenic'] = (100 * proj['intergenic_no_feature']
                                   / proj['reads_total'].replace(0, np.nan))
    print(proj.round(2).to_string())
    save_csv(proj.round(4), 'recovery_projection.csv')

    # Broken out per BAM as well, because within cellranger the
    # sample_alignments and unassigned_alignments BAMs mean different things:
    # the first is cell-associated reads, the second is droplets never called
    # as cells, which is where ambient virus and QC-dropped positives sit.
    per_bam = (cls.groupby(['source', 'bam'])[
        ['reads_total', 'counted_current_gtf', 'ambiguous_multi_feature',
         'intergenic_no_feature']].sum())
    per_bam = per_bam[per_bam['reads_total'] > 0]
    if not per_bam.empty:
        print('\n  broken out per BAM:')
        print(per_bam.to_string())
        save_csv(per_bam, 'recovery_projection_by_bam.csv')

    # ---- background inflation, the part that bears on the biology ----------
    # The rebuild plan's anchor claim is that the negative controls carry zero
    # viral signal. That claim is measured against the CURRENT GTF, which
    # discards the tat/rev/env overlap. The negatives' reads sit almost
    # entirely in exactly that overlap, so collapsing to one feature converts
    # invisible background into visible background. This table makes the size
    # of that effect explicit BEFORE the rerun, so the claim is re-measured
    # rather than inherited.
    #
    # Library-level grouping is approximate on purpose and flagged as such:
    # C carries 1,636 Pre cells from 41862 alongside control cells, and F
    # carries 8,626 genuinely infected 21 DPI cells alongside 931 control
    # cells. Any cell-level statement must use animal_id and timepoint.
    CONTROL_LIBS = ['A', 'B', 'G']        # unambiguously virus-free libraries
    INFECTED_LIBS_APPROX = ['D', 'E', 'H']  # unambiguously infected libraries
    cls2 = cls.copy()
    cls2['group'] = cls2['lib'].map(
        lambda x: 'control_ABG' if x in CONTROL_LIBS
        else ('infected_DEH' if x in INFECTED_LIBS_APPROX else 'mixed_CF'))
    bg = (cls2.groupby(['source', 'group'])[
        ['reads_total', 'counted_current_gtf']].sum())
    bg = bg.rename(columns={'counted_current_gtf': 'counted_now'})
    bg['counted_after_collapse'] = bg['reads_total']
    bg = bg[bg['reads_total'] > 0]
    if not bg.empty:
        print('\n  background inflation under the collapsed GTF:')
        print(bg.to_string())
        save_csv(bg, 'background_inflation.csv')
        for src_name in bg.index.get_level_values(0).unique():
            sub = bg.loc[src_name]
            if 'control_ABG' not in sub.index or 'infected_DEH' not in sub.index:
                continue
            c_now, c_aft = sub.loc['control_ABG', ['counted_now',
                                                   'counted_after_collapse']]
            i_now, i_aft = sub.loc['infected_DEH', ['counted_now',
                                                    'counted_after_collapse']]
            r_now = (i_now / c_now) if c_now else float('inf')
            r_aft = (i_aft / c_aft) if c_aft else float('inf')
            print(f'    {src_name:<12} infected:control counted reads  '
                  f'now {r_now:>8.1f}x   after collapse {r_aft:>8.1f}x')
        print('\n    If the ratio FALLS, the collapse recovers proportionally')
        print('    more background than signal, and the separation that the')
        print('    "clean negatives" claim rests on gets worse, not better.')
        print('    These are READS. UMI and cells may behave differently')
        print('    after CB+UMI dedup, which is what the H gate must measure.')
        print()
        print('    Read the sources separately, they do not mean the same')
        print('    thing. cellranger is the instrument the rebuild changes:')
        print('    if its control libraries carry zero viral reads at all,')
        print('    not merely zero COUNTED reads, then collapsing its GTF')
        print('    cannot manufacture background out of nothing and the')
        print('    "clean negatives" claim survives the GTF change intact.')
        print('    The STARsolo passes used a permissive viral-only index')
        print('    with relaxed length filters, deliberately, so reads in')
        print('    THEIR negatives reflect that stringency choice and not')
        print('    the GTF. Their background does not transfer to the rerun.')
        print()
        print('    The residual risk is the chemistry switch, not the GTF.')
        print('    SC5P-PE-v3 aligns R1\'s ~109 cDNA bases, which Cell Ranger')
        print('    has never aligned. If the R1 pass recruits reads in the')
        print('    negatives that the R2-only pass does not, paired-end may')
        print('    import some of that. Compare the control rows between')
        print('    cellranger and starsolo_r1 above to size it, and run a')
        print('    control library beside H in the gate rather than H alone.')

    print('\n  Reading this table:')
    print('    fold_recovery is the measured version of the projected ~4x.')
    print('    Near 1.0 means the GTF was NOT the bottleneck and the rerun is')
    print('    not justified on annotation grounds. Well above 1.0 means it was.')
    print('    pct_lost_ambiguous vs pct_lost_intergenic says whether the fix')
    print('    is ORF collapse (ambiguous) or LTR annotation (intergenic).')
    print('    Note fold_recovery is in READS. SHIV-positive CELLS scale more')
    print('    slowly, since recovered reads concentrate in cells already')
    print('    positive, and CB+UMI dedup collapses duplicates.')

# Where in the genome the loss sits, as a prediction check. This needs no
# reads, so it is reported whether or not any BAM was found.
thirds = {'5prime_0-33pct': (0, contig_len // 3),
          'middle_33-66pct': (contig_len // 3, 2 * contig_len // 3),
          '3prime_66-100pct': (2 * contig_len // 3, contig_len)}
print('\n  multiplicity by genome third (expected: 3prime is worst):')
for lab, (s, e) in thirds.items():
    seg = mult[s:e]
    print(f'    {lab:<18} mean features/base {seg.mean():4.2f} '
          f'| ambiguous {100 * (seg >= 2).mean():5.1f}% '
          f'| invisible {100 * (seg == 0).mean():5.1f}%')

print('\n  RECOMMENDED collapsed SHIV GTF for 01_mkref_collapsed.sh')
print('  (one gene == one transcript == one exon, spanning the whole contig,')
print('   LTRs included, no nested cassette feature):')
print()
_gid = 'SHIVAD8EO_provirus'
for feat in ('gene', 'transcript', 'exon'):
    attrs = f'gene_id "{_gid}"; '
    if feat != 'gene':
        attrs += f'transcript_id "{_gid}_t"; '
    attrs += f'gene_name "{_gid}";'
    if feat == 'exon':
        attrs += ' exon_number "1";'
    print(f'    SHIVAD8EO\tSHIVAD8EO\t{feat}\t1\t{contig_len}\t.\t+\t.\t{attrs}')
print()
print(f'  Span is 1..{contig_len} on the CURRENT FASTA. If the 34 bp cassette')
print('  is inserted, the contig lengthens and this end coordinate must move')
print('  with it. 01_mkref_collapsed.sh derives it from the FASTA rather than')
print('  hardcoding it, so the two cannot drift apart.')


# %% ============================================================
# CELL 7 — FIGURE: MULTIPLICITY TRACK WITH OBSERVED READ DEPTH
# ============================================================
hr('CELL 7 — figure')

fig, axes = plt.subplots(2, 1, figsize=(26, 14), sharex=True,
                         gridspec_kw={'height_ratios': [1, 1.3]})

ax = axes[0]
x = np.arange(contig_len)
ax.fill_between(x, 0, mult, where=(mult == 0), step='mid',
                color=HEX_INTERGENIC, label='0 features (invisible)')
ax.fill_between(x, 0, np.maximum(mult, 0), where=(mult == 1), step='mid',
                color=HEX_COUNTED, label='1 feature (countable)')
ax.fill_between(x, 0, mult, where=(mult >= 2), step='mid',
                color=HEX_AMBIG, label='2+ features (ambiguous)')
ax.set_ylabel('features\nper base', fontsize=FS_LABEL)
ax.set_title('SHIV GTF feature multiplicity and observed viral read depth',
             fontsize=FS_TITLE, pad=18)
ax.tick_params(labelsize=FS_TICK)
ax.legend(fontsize=FS_LEGEND, loc='upper left', frameon=False)
ax.grid(axis='y', color=HEX_GRID, lw=1)
ax.set_axisbelow(True)

# gene spans as stacked bars under the track
y0 = -0.6
step = 0.55
for k, (_, g) in enumerate(genes.iterrows()):
    yy = y0 - step * (k % 4)
    ax.plot([g['start0'], g['end0']], [yy, yy], lw=8, solid_capstyle='butt',
            color=HEX_AMBIG if g['width'] > 1500 else HEX_COUNTED)
    ax.text((g['start0'] + g['end0']) / 2, yy - 0.18, g['gene_name'],
            ha='center', va='top', fontsize=FS_TICK - 6)
ax.set_ylim(y0 - step * 4 - 0.6, max(2, int(mult.max())) + 0.6)

ax = axes[1]
if depth.max() > 0:
    ax.fill_between(x, 0, depth, step='mid', color=HEX_DEPTH, alpha=0.85)
    ax.set_yscale('log')
    ax.set_ylabel('aligned read depth\n(all sources, log)', fontsize=FS_LABEL)
else:
    ax.text(0.5, 0.5, 'no aligned viral reads found', ha='center',
            va='center', transform=ax.transAxes, fontsize=FS_LABEL)
    ax.set_ylabel('aligned read depth', fontsize=FS_LABEL)
ax.set_xlabel('position on SHIVAD8EO contig (bp)', fontsize=FS_LABEL)
ax.tick_params(labelsize=FS_TICK)
ax.grid(color=HEX_GRID, lw=1)
ax.set_axisbelow(True)
ax.set_xlim(0, contig_len)

savefig(fig, 'F00_gtf_multiplicity')

hr('DONE — 00_gtf_multiplicity_diagnostic.py')
print(f'  outputs: {OUT_DIR}')
print('  nothing outside that directory was written.')
print('\n  Next: read recovery_projection.csv before running '
      '01_mkref_collapsed.sh.')
