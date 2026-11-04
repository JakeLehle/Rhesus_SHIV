#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
03b — Gate evaluation: does the v11 rebuild justify the full array? (READ-ONLY)
===============================================================================

Answers the gate criteria with computed numbers rather than remembered ones,
measuring v11 and the v10 baseline side by side with the same code.

RUNTIME: seconds. The first version of this script took hours and that was a
bug, recorded here so it is not reintroduced.

WHAT WAS SLOW, AND WHY IT IS NOT ANY MORE
------------------------------------------
The original asked whether the reads were paired like this:

    samtools view -c -F 0x900 bam            # NO region argument
    samtools view -c -f 0x1 -F 0x900 bam     # NO region argument

With no region, `samtools view -c` decompresses and scans the ENTIRE BAM. A
20,000-cell 5' library is tens to low hundreds of GB, and that was two full
passes per BAM across H and G in both versions: eight complete scans to answer
one yes/no question. Everything else in the script was already seconds.

Three replacements, each turning a full scan into an index or header lookup:

  1. CHEMISTRY comes from web_summary.html, which is Cell Ranger's own record
     of the chemistry it resolved. Reading it is instant. A bounded sample of
     the first SAMPLE_READS alignments confirms it independently, because
     pairing is a property of every read and the first few thousand settle it
     as well as all of them.
  2. VIRAL READ TOTALS come from `samtools idxstats`, which reads the .bai and
     never touches the BAM body.
  3. FILE DISCOVERY tries the known `multi` output paths before falling back
     to a glob. A recursive `**` walk over a multi outs tree crosses hundreds
     of thousands of VDJ consensus files.

Set EXACT_READ_COUNTS = True only if you want the slow exhaustive counts back.

HOW THE CHEMISTRY IS VERIFIED
------------------------------
Not from metrics_summary.csv: `cellranger multi` emits no chemistry row there,
so its absence means nothing. Two independent sources are used instead. The
web summary states what Cell Ranger resolved. The BAM flags show what the
aligner did, since SC5P-PE-v3 aligns both mates and SC5P-R2-v3 aligns R2 only,
so v11 must carry paired alignments and v10 must not.

FEATURE NAMES DIFFER BETWEEN THE REFERENCES, deliberately
----------------------------------------------------------
  v10   ten per-ORF features: gag pol vif vpx vpr tat rev vpu env nef
  v11   one collapsed feature: SHIVAD8EO
The viral total is summed across the ten for v10 and read off the one for v11.
Using a single name for both silently reports zero for whichever reference
lacks it. Matching is CASE-SENSITIVE: host TAT, tyrosine aminotransferase,
differs from viral *tat* only by case, and is reported separately.

Outputs (all under OUT_DIR, nothing is modified)
  chemistry_check.csv, viral_strand.csv, gate_summary.csv,
  viral_umi_distribution.csv, host_metrics.csv, gate_verdict.txt

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
from collections import OrderedDict

import numpy as np
import pandas as pd

try:
    from shiv_config import WORKDIR                                 # noqa: F401
except Exception:                                                   # noqa: BLE001
    WORKDIR = '/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV'

V11_ROOT = os.path.join(WORKDIR, 'pipeline_results', '03_cellranger_multi')
V10_ROOT = os.path.join(WORKDIR, 'Analysis')

LIBS = OrderedDict([
    ('H', 'Necropsy_LN'),      # the signal library
    ('G', 'Necropsy_PBMC'),    # the control library
])
SIGNAL_LIB, CONTROL_LIB = 'H', 'G'

# --- speed dials -----------------------------------------------------------
SAMPLE_READS = 5000          # reads sampled from the head of each BAM
EXACT_READ_COUNTS = False    # True restores the full-BAM scans. Hours. Avoid.

# --- viral features, per reference version ---------------------------------
V11_VIRAL_FEATURES = ['SHIVAD8EO']
V10_VIRAL_FEATURES = ['gag', 'pol', 'vif', 'vpx', 'vpr',
                      'tat', 'rev', 'vpu', 'env', 'nef']
HOST_COLLISION = ['TAT']     # host tyrosine aminotransferase, NOT viral

V11_CONTIG = 'SHIVAD8EO_V4fullN'
V10_CONTIG = 'SHIVAD8EO'

MIN_SHIV_UMI = 2             # a cell is SHIV+ at >= this many viral UMI

# --- gate thresholds, from the rebuild plan --------------------------------
GATE_MIN_POS_CELLS_H = 20
GATE_BASELINE_POS_H = 3
GATE_BASELINE_UMI_H = 29
ANTISENSE_ALARM = 0.80

OUT_DIR = os.path.join(WORKDIR, 'pipeline_results', '03b_gate_evaluation')
os.makedirs(OUT_DIR, exist_ok=True)
pd.set_option('display.width', 220)
pd.set_option('display.max_columns', 60)


def hr(t):
    print('\n' + '=' * 94, flush=True)
    print(t, flush=True)
    print('=' * 94, flush=True)


def save_csv(df, name, index=False):
    df.to_csv(os.path.join(OUT_DIR, name), index=index, lineterminator='\n')
    print(f'  wrote {name}  ({df.shape[0]} rows)', flush=True)


# samtools: the system build at /usr/local/bin is broken and is skipped
SAMTOOLS = None
for cand in [os.path.join(os.environ.get('CONDA_PREFIX', ''), 'bin', 'samtools'),
             shutil.which('samtools')]:
    if cand and os.path.exists(cand) and not cand.startswith('/usr/local/bin'):
        SAMTOOLS = cand
        break
if SAMTOOLS is None:
    print('WARN  no usable samtools. The system build at /usr/local/bin is\n'
          '      broken and is deliberately skipped. Run:\n'
          '        conda activate sc_pre\n'
          '        export PATH="${CONDA_PREFIX}/bin:${PATH}"\n'
          '      BAM checks will be skipped; matrix checks still run.')
else:
    print(f'  samtools: {SAMTOOLS}')


# --- fast file discovery ----------------------------------------------------
# Known cellranger multi layout first. The glob fallback is bounded to two
# directory levels rather than a recursive walk, because a multi outs tree
# holds a VDJ consensus file per clonotype and `**` crosses all of them.
KNOWN = 'outs/per_sample_outs/{sid}/{name}'


def find_one(root, sid, name):
    direct = os.path.join(root, sid, KNOWN.format(sid=sid, name=name))
    if os.path.exists(direct):
        return direct
    for pat in (os.path.join(root, sid, 'outs', name),
                os.path.join(root, sid, 'outs', '*', name),
                os.path.join(root, sid, 'outs', '*', '*', name)):
        hits = sorted(glob.glob(pat))
        if hits:
            return hits[0]
    return None


def sh(args):
    """Run a command, return stdout or None. Never raises."""
    try:
        return subprocess.run(args, capture_output=True, text=True,
                              check=True).stdout
    except Exception:                                               # noqa: BLE001
        return None


def head_sample(bam, n=SAMPLE_READS):
    """First n alignment records, without reading the whole BAM.

    Streams `samtools view` and stops after n lines, so only the leading BGZF
    blocks are decompressed. Pairing is a property of every read, so n in the
    thousands settles it exactly as well as a full scan would.
    """
    try:
        p = subprocess.Popen([SAMTOOLS, 'view', bam], stdout=subprocess.PIPE,
                             stderr=subprocess.DEVNULL, text=True)
    except Exception:                                               # noqa: BLE001
        return []
    out = []
    try:
        for i, line in enumerate(p.stdout):
            if i >= n:
                break
            out.append(line)
    finally:
        p.kill()
        p.wait()
    return out


# %% ============================================================
# CELL 2 — LOCATE THE OUTPUTS
# ============================================================
hr('CELL 2 — locate outputs')

FOUND = {}
for lib, cond in LIBS.items():
    sid = f'SHIV_{cond}'
    for ver, root in (('v11', V11_ROOT), ('v10', V10_ROOT)):
        FOUND[(lib, ver)] = {
            'sid': sid,
            'bam': find_one(root, sid, 'sample_alignments.bam'),
            'h5': find_one(root, sid, 'sample_filtered_feature_bc_matrix.h5'),
            'met': find_one(root, sid, 'metrics_summary.csv'),
            'web': find_one(root, sid, 'web_summary.html'),
        }
        have = [k for k in ('bam', 'h5', 'met', 'web') if FOUND[(lib, ver)][k]]
        miss = [k for k in ('bam', 'h5', 'met', 'web') if not FOUND[(lib, ver)][k]]
        print(f'  {lib} {ver}: found {",".join(have) or "nothing"}'
              + (f'   missing {",".join(miss)}' if miss else ''), flush=True)


# %% ============================================================
# CELL 3 — CHEMISTRY, FROM TWO INDEPENDENT SOURCES (instant)
# ============================================================
hr('CELL 3 — chemistry')

CHEM_RE = re.compile(r'SC5P-[A-Za-z0-9\-]+|Single Cell 5. ?(?:PE|R2)[^<",]{0,12}')

chem_rows = []
for (lib, ver), f in FOUND.items():
    rec = {'lib': lib, 'version': ver}

    # source 1: Cell Ranger's own web summary
    if f['web']:
        try:
            with open(f['web'], errors='ignore') as fh:
                txt = fh.read()
            hits = sorted(set(m.group(0).strip() for m in CHEM_RE.finditer(txt)))
            rec['web_summary_says'] = '; '.join(hits[:4]) if hits else 'no match'
        except Exception as e:                                      # noqa: BLE001
            rec['web_summary_says'] = f'unreadable ({type(e).__name__})'
    else:
        rec['web_summary_says'] = 'web_summary.html not found'

    # source 2: what the aligner actually did, from a bounded read sample
    if f['bam'] and SAMTOOLS:
        rows = head_sample(f['bam'])
        if rows:
            flags = [int(r.split('\t')[1]) for r in rows if '\t' in r]
            paired = sum(1 for fl in flags if fl & 0x1)
            rec['reads_sampled'] = len(flags)
            rec['sampled_paired'] = paired
            rec['pct_paired_sampled'] = round(100 * paired / len(flags), 1)
            rec['inferred'] = ('paired-end' if paired > 0.5 * len(flags)
                               else 'single-end')
        else:
            rec['inferred'] = 'no reads sampled'
    else:
        rec['inferred'] = 'bam or samtools missing'

    if EXACT_READ_COUNTS and f['bam'] and SAMTOOLS:
        # Opt-in only. This is the full-BAM scan that made the first version
        # of this script take hours.
        tot = sh([SAMTOOLS, 'view', '-c', '-F', '0x900', f['bam']])
        pr = sh([SAMTOOLS, 'view', '-c', '-f', '0x1', '-F', '0x900', f['bam']])
        if tot and pr:
            rec['reads_primary_exact'] = int(tot.strip())
            rec['reads_paired_exact'] = int(pr.strip())

    chem_rows.append(rec)
    print(f'  {lib} {ver}: {rec.get("inferred")} '
          f'({rec.get("pct_paired_sampled", "NA")}% of '
          f'{rec.get("reads_sampled", 0)} sampled)   '
          f'web: {rec["web_summary_says"]}', flush=True)

chem = pd.DataFrame(chem_rows)
save_csv(chem, 'chemistry_check.csv')

v11_inf = chem[chem.version == 'v11']['inferred'].astype(str)
v10_inf = chem[chem.version == 'v10']['inferred'].astype(str)
chem_ok = len(v11_inf) > 0 and v11_inf.str.startswith('paired-end').all()
if chem_ok:
    print('\n  VERDICT: v11 is paired-end. The chemistry change took effect.')
    if len(v10_inf) and v10_inf.str.startswith('single-end').all():
        print('           v10 is single-end, as expected for SC5P-R2-v3.')
else:
    print('\n  *** VERDICT: v11 does NOT look paired-end. The config line was')
    print('      not honoured, and any viral comparison confounds the GTF')
    print('      change with a chemistry that never applied.')


# %% ============================================================
# CELL 4 — VIRAL READS AND STRAND (index lookups, instant)
# ============================================================
hr('CELL 4 — viral reads on the contig, and which strand')

strand_rows = []
for (lib, ver), f in FOUND.items():
    contig = V11_CONTIG if ver == 'v11' else V10_CONTIG
    rec = {'lib': lib, 'version': ver, 'contig': contig}
    if f['bam'] and SAMTOOLS:
        # idxstats reads the .bai only; it never touches the BAM body
        idx = sh([SAMTOOLS, 'idxstats', f['bam']])
        if idx:
            for line in idx.splitlines():
                p = line.split('\t')
                if len(p) >= 4 and p[0] == contig:
                    rec['contig_len'] = int(p[1])
                    rec['mapped_idxstats'] = int(p[2])
                    break
            else:
                rec['note'] = f'contig {contig} not in this BAM'
        # region queries use the index and touch only the viral slice
        fwd = sh([SAMTOOLS, 'view', '-c', '-F', '0x914', f['bam'], contig])
        rev = sh([SAMTOOLS, 'view', '-c', '-F', '0x904', '-f', '0x10',
                  f['bam'], contig])
        if fwd is not None and rev is not None:
            fwd, rev = int(fwd.strip()), int(rev.strip())
            rec.update(viral_fwd=fwd, viral_rev=rev, viral_total=fwd + rev)
            rec['pct_antisense'] = (round(100 * rev / (fwd + rev), 1)
                                    if fwd + rev else np.nan)
    else:
        rec['note'] = 'bam or samtools missing'
    strand_rows.append(rec)

strand = pd.DataFrame(strand_rows)
print(strand.to_string(index=False), flush=True)
save_csv(strand, 'viral_strand.csv')

bad = strand[strand.get('pct_antisense', pd.Series(dtype=float))
             > 100 * ANTISENSE_ALARM]
if len(bad):
    print('\n  *** Over 80% antisense against a plus-strand-only GTF feature.')
    print('      A low viral count would then be a STRAND problem, not a GTF')
    print('      problem. Resolve before reading the next criteria.')
else:
    print('\n  No strand alarm.')


# %% ============================================================
# CELL 5 — VIRAL UMI AND SHIV+ CELLS FROM THE MATRICES
# ============================================================
hr('CELL 5 — viral UMI and SHIV+ cells')

try:
    import scanpy as sc
except ImportError:
    raise SystemExit('scanpy not available. conda activate sc_pre')


def viral_from_matrix(h5, features):
    """Viral UMI and SHIV+ cells from a filtered matrix.

    Case-SENSITIVE name matching. Host TAT differs from viral *tat* only by
    case, so a case-insensitive match would fold tyrosine aminotransferase
    into the viral total and inflate the v10 baseline.
    """
    a = sc.read_10x_h5(h5)
    a.var_names_make_unique()
    present = [f for f in features if f in a.var_names]
    missing = [f for f in features if f not in a.var_names]
    collision = [f for f in HOST_COLLISION if f in a.var_names]
    if not present:
        return dict(n_cells=int(a.n_obs), features_found=0,
                    features_missing=';'.join(missing), viral_umi=0.0,
                    cells_any=0, cells_pos=0, max_umi_in_a_cell=0.0,
                    host_TAT_in_matrix=bool(collision),
                    note='no viral feature present', _per_cell=np.array([]))
    per_cell = np.asarray(a[:, present].X.sum(axis=1)).ravel()
    return dict(
        n_cells=int(a.n_obs), features_found=len(present),
        features_present=';'.join(present), features_missing=';'.join(missing),
        host_TAT_in_matrix=bool(collision),
        viral_umi=float(per_cell.sum()),
        cells_any=int((per_cell > 0).sum()),
        cells_pos=int((per_cell >= MIN_SHIV_UMI).sum()),
        max_umi_in_a_cell=float(per_cell.max()) if per_cell.size else 0.0,
        _per_cell=per_cell)


rows, dist_rows = [], []
for (lib, ver), f in FOUND.items():
    feats = V11_VIRAL_FEATURES if ver == 'v11' else V10_VIRAL_FEATURES
    if not f['h5']:
        print(f'  WARN  {lib} {ver}: no filtered matrix found', flush=True)
        rows.append({'lib': lib, 'version': ver, 'note': 'matrix not found'})
        continue
    r = viral_from_matrix(f['h5'], feats)
    pc = r.pop('_per_cell')
    r.update(lib=lib, version=ver)
    rows.append(r)
    for v in np.unique(pc[pc > 0]):
        dist_rows.append({'lib': lib, 'version': ver,
                          'viral_umi_in_cell': float(v),
                          'n_cells': int((pc == v).sum())})
    print(f"  {lib} {ver}: {r['n_cells']:,} cells | viral UMI "
          f"{r['viral_umi']:,.0f} | cells >=1 {r['cells_any']} | "
          f"SHIV+ (>={MIN_SHIV_UMI}) {r['cells_pos']}", flush=True)

res = pd.DataFrame(rows)
cols = [c for c in ['lib', 'version', 'n_cells', 'viral_umi', 'cells_any',
                    'cells_pos', 'max_umi_in_a_cell', 'features_found',
                    'features_present', 'host_TAT_in_matrix', 'note']
        if c in res.columns]
print()
print(res[cols].to_string(index=False), flush=True)
save_csv(res, 'gate_summary.csv')
if dist_rows:
    save_csv(pd.DataFrame(dist_rows), 'viral_umi_distribution.csv')


# %% ============================================================
# CELL 6 — HOST METRICS (criterion 5)
# ============================================================
hr('CELL 6 — host metrics')


def gex_metrics(path):
    """Gene Expression rows from a multi metrics_summary.csv.

    Six columns: Category, Library Type, Grouped By, Group Name, Metric Name,
    Metric Value. Filtering on Library Type is the step 03 skipped, which is
    why it reported Antibody Capture numbers as if they were GEX.
    """
    if path is None:
        return {}
    df = pd.read_csv(path)
    if 'Library Type' not in df.columns:
        return {'_note': 'unexpected columns: ' + ','.join(df.columns)}
    gex = df[df['Library Type'].astype(str).str.contains('Gene Expression',
                                                         na=False)]
    out = {}
    for want in ['Cells', 'Median UMI counts per cell', 'Median genes per cell',
                 'Mean reads per cell', 'Sequencing saturation',
                 'Reads mapped confidently to transcriptome',
                 'Reads mapped antisense to gene']:
        hit = gex[gex['Metric Name'].astype(str)
                  .str.fullmatch(want, case=False, na=False)]
        if len(hit):
            out[want] = str(hit.iloc[0]['Metric Value'])
    return out


host_rows = []
for (lib, ver), f in FOUND.items():
    m = gex_metrics(f['met'])
    m.update(lib=lib, version=ver)
    host_rows.append(m)
host = pd.DataFrame(host_rows)
front = ['lib', 'version']
host = host[front + [c for c in host.columns
                     if c not in front and not c.startswith('_')]]
print(host.to_string(index=False), flush=True)
save_csv(host, 'host_metrics.csv')


# %% ============================================================
# CELL 7 — VERDICT
# ============================================================
hr('CELL 7 — gate verdict')


def get(lib, ver, col):
    s = res[(res.lib == lib) & (res.version == ver)]
    return (s.iloc[0][col] if len(s) and col in s.columns
            and pd.notna(s.iloc[0][col]) else None)


lines = []


def verdict(n, text, state, detail=''):
    line = f'  [{ {"PASS":"PASS","FAIL":"FAIL","CHECK":"????"}[state] }] {n}. {text}'
    if detail:
        line += f'\n         {detail}'
    lines.append(line)
    print(line, flush=True)


h11p, h10p = get(SIGNAL_LIB, 'v11', 'cells_pos'), get(SIGNAL_LIB, 'v10', 'cells_pos')
h11u, h10u = get(SIGNAL_LIB, 'v11', 'viral_umi'), get(SIGNAL_LIB, 'v10', 'viral_umi')
g11u, g10u = get(CONTROL_LIB, 'v11', 'viral_umi'), get(CONTROL_LIB, 'v10', 'viral_umi')

verdict(0, 'chemistry: v11 aligns paired-end',
        'PASS' if chem_ok else 'FAIL',
        '' if chem_ok else 'the viral comparison below is not interpretable')

if len(bad):
    verdict(1, 'strand sanity on the viral contig', 'FAIL',
            'over 80% antisense against a plus-strand-only feature')
else:
    verdict(1, 'strand sanity on the viral contig', 'PASS')

if h11p is None:
    verdict(2, f'{SIGNAL_LIB} SHIV+ cells', 'CHECK', 'v11 matrix not readable')
else:
    base = h10p if h10p is not None else GATE_BASELINE_POS_H
    verdict(2, f'{SIGNAL_LIB} SHIV+ cells {base:.0f} -> {h11p:.0f} '
               f'(threshold {GATE_MIN_POS_CELLS_H})',
            'PASS' if h11p >= GATE_MIN_POS_CELLS_H else 'FAIL',
            '' if h11p >= GATE_MIN_POS_CELLS_H else
            'below threshold: the GTF was not the bottleneck and the full '
            'array is not justified on these grounds')

if h11u is None:
    verdict(3, f'{SIGNAL_LIB} viral UMI', 'CHECK', 'v11 matrix not readable')
else:
    base = h10u if h10u is not None else GATE_BASELINE_UMI_H
    verdict(3, f'{SIGNAL_LIB} viral UMI {base:,.0f} -> {h11u:,.0f}',
            'PASS' if h11u > base else 'FAIL')

if g11u is None or h11u is None:
    verdict(4, f'{CONTROL_LIB} stays near zero', 'CHECK', 'matrices not readable')
else:
    r11 = (h11u / g11u) if g11u else float('inf')
    r10 = (h10u / g10u) if (h10u and g10u) else float('inf')
    verdict(4, f'{CONTROL_LIB} viral UMI {(g10u or 0):,.0f} -> {g11u:,.0f}; '
               f'{SIGNAL_LIB}:{CONTROL_LIB} {r10:.1f}x -> {r11:.1f}x',
            'PASS' if (g11u == 0 or r11 >= r10) else 'FAIL',
            '' if (g11u == 0 or r11 >= r10) else
            'separation got WORSE: the collapse recovers background at least '
            'as fast as signal, and the clean-negatives claim needs restating')

verdict(5, 'host cell count and median genes per cell stay sane', 'CHECK',
        'read the Cell 6 table; a large shift means QC thresholds moved and '
        'annotation needs rework')

print('\n  Criterion 6, cassette-spanning reads, is not computed here. It needs')
print('  the flank-anchored extraction rather than the matrix, so run')
print('  shiv_bc_core_extract.py against the v11 BAMs for that one.')

with open(os.path.join(OUT_DIR, 'gate_verdict.txt'), 'w') as fh:
    fh.write('Gate verdict, v11 vs v10\n' + '=' * 60 + '\n')
    fh.write('\n'.join(lines) + '\n')
print('\n  wrote gate_verdict.txt')

hr('DONE — 03b_gate_evaluation.py')
print(f'  outputs: {OUT_DIR}')
print('  nothing was modified.')
print('\n  If 0 through 4 pass, run:  sbatch --array=0-5 03_cellranger_multi.sh')
