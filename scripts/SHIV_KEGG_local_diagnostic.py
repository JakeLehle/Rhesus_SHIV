#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Rhesus SHIV - why did the local KEGG engine return zero?
=========================================================
The run finished with `local=0` on every one of the 35 gene lists, while
g:Profiler's KEGG returned 51 significant terms over the same lists with the
same backgrounds. That gap is the whole reason the second engine exists, so it
gets explained rather than shrugged off.

The local engine is not obviously broken. Its self-test passed, reproducing a
hand-computed hypergeometric to 1e-12, and the background-sensitivity check
confirmed the background is actually applied. So the question is not whether
the arithmetic is right, it is where the signal disappears between the gene
list and a corrected p value.

There are five candidates, and they are distinguishable by counting:

  A. SYMBOL -> ENTREZ COLLAPSE. Only 75.3% of symbols map, and the query is hit
     twice, once directly and once through the background.
  B. KEGG COVERAGE. Only 45.4% of mapped genes are in any mcc pathway, so the
     effective query is smaller again.
  C. MIN_OVERLAP = 3. With a small query and a median pathway of 86 genes, the
     expected overlap can sit below 3, and every term is dropped before a p
     value is ever computed.
  D. BH OVER 368 PATHWAYS. Raw p values could be fine while nothing survives
     correction.
  E. AN ACTUAL BUG, e.g. an identifier namespace mismatch that makes every
     intersection empty.

E is the one that matters, and it has a signature the others do not: if the
namespaces disagree, the observed overlap is exactly zero for EVERY pathway in
EVERY contrast. Any non-zero overlap anywhere rules it out. That is the first
thing this script checks.

It then walks one contrast through the funnel and prints the surviving count at
each stage, so the real cause is read off a table instead of argued about.

Read only. Seconds to run. No network.

Author: Jake Lehle
Date: September 2026
"""

import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import hypergeom

WORKING_DIR = '/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV'
ANN = os.path.join(WORKING_DIR, 'annotation_output')
ENR = os.path.join(ANN, 'enrichment')
MAP_CSV = os.path.join(ENR, '00_gene_id_mapping', 'gene_id_map.csv')
GMT = os.path.join(ENR, '01_gene_sets', 'KEGG_mcc.gmt')
SIG_CSV = os.path.join(ENR, '03_summary', 'enrichment_significant.csv')
ORA_DIR = os.path.join(ANN, 'deg_comparisons', 'ora_gene_lists')

GATE, LAYER = 'FDR', 'pseudobulk'
ALPHA, MIN_SET, MAX_SET, MIN_OVERLAP = 0.05, 5, 500, 3
N_CONTRASTS = 6          # walk this many of the biggest lists through the funnel


def norm_id(x):
    """Same normalisation the fixed run script applies."""
    s = str(x).strip()
    if s.endswith('.0') and s[:-2].isdigit():
        s = s[:-2]
    for pre in ('mcc:', 'ncbi-geneid:', 'hsa:', 'ncbi-proteinid:'):
        if s.startswith(pre):
            s = s[len(pre):]
    return s


def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.clip(q, 0, 1)


def hr(t):
    print('\n' + '=' * 90)
    print(t)
    print('=' * 90)


hr('CELL 1 - load')
for p in (MAP_CSV, GMT):
    if not os.path.exists(p):
        sys.exit(f'[FATAL] missing {p}')

id_map = pd.read_csv(MAP_CSV)
raw2ez = {str(r.symbol).strip(): str(r.entrez) for r in id_map.itertuples()
          if pd.notna(getattr(r, 'entrez', None))}
sym2ez = {k: norm_id(v) for k, v in raw2ez.items()}
print(f'  symbols in map      : {len(id_map):,}')
print(f'  with an entrez id   : {len(sym2ez):,}')

sets, names, raw_gmt = {}, {}, set()
with open(GMT) as fh:
    for line in fh:
        parts = line.rstrip('\n').split('\t')
        if len(parts) > 2:
            sets[parts[0]] = {norm_id(g) for g in parts[2:]}
            names[parts[0]] = parts[1]
            raw_gmt.update(parts[2:])
print(f'  KEGG mcc pathways   : {len(sets):,}')

all_set_genes = set().union(*sets.values()) if sets else set()
print(f'  distinct genes in the GMT : {len(all_set_genes):,}')

hr('CELL 2 - THE ONE THAT MATTERS: do the two namespaces even intersect?')
# If the GMT stores ids one way (e.g. "mcc:694372") and the map another
# ("694372"), every intersection is empty and the engine silently returns
# nothing forever. This is candidate E.
# Check RAW first, then normalised. If raw fails and normalised succeeds, the
# cause is a string-formatting mismatch, which is exactly what happened here:
# pandas types the entrez column float64 because it has blanks, so str() gives
# '694372.0' while the GMT holds '694372'.
raw_shared = set(raw2ez.values()) & raw_gmt
ez_vals = set(sym2ez.values())
shared = ez_vals & all_set_genes
print(f'  RAW intersection        : {len(raw_shared):,}')
print(f'    raw map sample  : {sorted(set(raw2ez.values()))[:3]}')
print(f'    raw GMT sample  : {sorted(raw_gmt)[:3]}')
if not raw_shared and shared:
    print('\n  >>> CAUSE FOUND: candidate E, an identifier FORMATTING mismatch.')
    print('      Raw ids share nothing; normalised ids share '
          f'{len(shared):,}. So the')
    print('      engine compared \'694372.0\'-style strings against the GMT\'s')
    print('      \'694372\', every overlap was zero, and it reported nothing.')
    print('      This is a BUG and it is fixed in the updated run script.')
    print('      g:Profiler KEGG and all GO results are UNAFFECTED: that')
    print('      engine never touches these identifiers.')
    print('\n  Re-run SHIV_enrichment_run.py to get the KEGG cross-check.')
    sys.exit(0)
print(f'  entrez ids from gene_id_map : {len(ez_vals):,}')
print(f'  gene ids inside the GMT     : {len(all_set_genes):,}')
print(f'  INTERSECTION                : {len(shared):,}')
print(f'  sample from map : {sorted(list(ez_vals))[:5]}')
print(f'  sample from GMT : {sorted(list(all_set_genes))[:5]}')

if len(shared) == 0:
    print('\n  >>> CAUSE FOUND: candidate E, namespace mismatch.')
    print('      The two id spaces do not overlap at all, so every pathway')
    print('      overlap is zero by construction and no p value is ever')
    print('      computed. This is a BUG, not biology. The g:Profiler KEGG')
    print('      results are unaffected, because that engine never touches')
    print('      these identifiers.')
    sys.exit(0)

print(f'\n  Namespaces agree ({len(shared):,} ids in common), so candidate E is')
print('  ruled out. The loss is somewhere in the funnel below.')

hr('CELL 3 - walk the biggest lists through the funnel')
lists = sorted([f for f in os.listdir(ORA_DIR)
                if f.endswith(('_up.csv', '_down.csv'))
                and f'_{LAYER}_{GATE}_' in f])
if not lists:
    sys.exit(f'[FATAL] no {GATE}/{LAYER} lists in {ORA_DIR}')

sized = []
for fn in lists:
    n = len(pd.read_csv(os.path.join(ORA_DIR, fn)))
    sized.append((n, fn))
sized.sort(reverse=True)

rows = []
for n_in, fn in sized[:N_CONTRASTS]:
    cid = fn.split('_')[0]
    stem = fn[:-4]
    bgf = os.path.join(ORA_DIR, f'{stem.rsplit("_" + GATE, 1)[0]}_background.csv')
    if not os.path.exists(bgf):
        cand = [f for f in os.listdir(ORA_DIR)
                if f.startswith(cid) and f.endswith('_background.csv')
                and f'_{LAYER}_' in f]
        if not cand:
            continue
        bgf = os.path.join(ORA_DIR, cand[0])

    genes = pd.read_csv(os.path.join(ORA_DIR, fn))['gene'].astype(str).tolist()
    bg = pd.read_csv(bgf)['gene'].astype(str).tolist()

    q_e = {sym2ez[g] for g in genes if g in sym2ez}
    b_e = {sym2ez[g] for g in bg if g in sym2ez}
    q_e &= b_e
    M, N = len(b_e), len(q_e)

    n_in_any = len(q_e & all_set_genes)

    tested = kept_overlap = 0
    best_p = np.nan
    ps, ks = [], []
    for sid, gs in sets.items():
        gset = gs & b_e
        if len(gset) < MIN_SET or len(gset) > MAX_SET:
            continue
        tested += 1
        k = len(q_e & gset)
        ks.append(k)
        if k < MIN_OVERLAP:
            continue
        kept_overlap += 1
        ps.append(float(hypergeom.sf(k - 1, M, len(gset), N)))

    n_raw = n_fdr = 0
    if ps:
        ps = np.array(ps)
        best_p = float(ps.min())
        n_raw = int((ps < ALPHA).sum())
        n_fdr = int((bh(ps) < ALPHA).sum())

    rows.append({
        'contrast': cid,
        'direction': 'up' if fn.endswith('_up.csv') else 'down',
        'genes_in': n_in,
        'A_mapped': N,
        'B_in_any_pathway': n_in_any,
        'bg_mapped': M,
        'pathways_tested': tested,
        'max_overlap': int(max(ks)) if ks else 0,
        'C_pass_min_overlap': kept_overlap,
        'best_raw_p': best_p,
        'D_raw_sig': n_raw,
        'D_fdr_sig': n_fdr,
    })

res = pd.DataFrame(rows)
pd.set_option('display.width', 220)
print(res.to_string(index=False))

hr('CELL 4 - verdict')
if res.empty:
    print('  no lists could be evaluated')
    sys.exit(0)

print('  Reading the funnel left to right, the count that collapses names the')
print('  cause.\n')

if res['max_overlap'].max() == 0:
    print('  >>> Every pathway overlap is ZERO despite the namespaces sharing')
    print('      ids. Something is wrong in how the query is built. BUG.')
elif res['C_pass_min_overlap'].max() == 0:
    print(f'  >>> CAUSE: candidate C. Overlaps are non-zero (max '
          f'{int(res["max_overlap"].max())}) but none reaches MIN_OVERLAP = '
          f'{MIN_OVERLAP},')
    print('      so no term is ever tested. Not a bug: the query is simply too')
    print('      small once it is restricted to genes in KEGG mcc pathways.')
    print('      Lowering MIN_OVERLAP to 2 would produce terms, but 2-gene')
    print('      overlaps are not worth reporting, so I would leave it.')
elif res['D_raw_sig'].max() == 0:
    print('  >>> CAUSE: candidate A+B. Terms are tested, but the effective')
    print('      query is so eroded by symbol mapping and KEGG coverage that')
    print('      nothing is even nominally significant. Real, not a bug.')
elif res['D_fdr_sig'].max() == 0:
    print('  >>> CAUSE: candidate D, the correction. Raw p values are')
    print(f'      significant (best {res["best_raw_p"].min():.2e}) but nothing')
    print(f'      survives BH over {int(res["pathways_tested"].median())} '
          f'pathways.')
    print('      This is the honest answer and it means the local engine is')
    print('      simply stricter than g:Profiler\'s g_SCS on this data. It does')
    print('      NOT invalidate the g:Profiler KEGG terms, but it does mean')
    print('      they are not independently corroborated, and the KEGG figure')
    print('      should be presented with that stated.')
else:
    print('  >>> The engine DOES produce significant terms here. The zero in')
    print('      the main run then came from a filter applied afterwards, not')
    print('      from the statistics. Send me this table.')

print('\n  Either way, the GO results and the g:Profiler KEGG p values are')
print('  untouched by this. Only the independent KEGG cross-check is affected.')
