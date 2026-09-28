#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Rhesus SHIV - Over-Representation Analysis: KEGG and GO, per contrast, up and down
===================================================================================
Step 2 of 2. Runs KEGG and GO over-representation on the DEG gene lists agreed
with Zoey (2026-09-28): |log2FC| > 1.5 and FDR < 0.05, per-animal layer.

This is OVER-REPRESENTATION (a hypergeometric test on a gene LIST against a
background universe), not GSEA in the Subramanian ranked-list sense. Up and down
genes are tested SEPARATELY, because a pathway enriched among up genes and one
enriched among down genes are different results and pooling them hides direction.

TWO INDEPENDENT ENGINES, deliberately
--------------------------------------
  A. g:Profiler, organism 'mmulatta'. One call returns GO:BP, GO:MF, GO:CC,
     KEGG, Reactome and WikiPathways. Server-side statistics with the g:SCS
     correction, which is g:Profiler's own multiple-testing method and is less
     conservative than Bonferroni while still accounting for term overlap.
  B. A local hypergeometric against the KEGG mcc GMT built in step 1. Implemented
     here rather than delegated, for one reason: the background universe has to
     be THIS contrast's tested gene set, and a wrong background silently inflates
     every p value. Thirty lines of scipy we can verify beats a library call we
     cannot inspect. CELL 3 validates the implementation against a worked example
     with a known answer before any real data touches it.

Where the two disagree on KEGG, that is worth knowing: g:Profiler is Ensembl-
backed and KEGG REST is KEGG-backed, so the gene-to-pathway assignments differ
slightly. Agreement is reassurance; disagreement is a flag, not an error.

Nothing here converts rhesus genes to human orthologs. Both engines run on
Macaca mulatta annotations natively.

Read-only w.r.t. everything upstream. Spyder cells (# %%).

Author: Jake Lehle
Date: September 2026
"""

# %% ============================================================
# CELL 1 - CONFIG AND DIALS
# ============================================================
import os
import re
import json
import time
import urllib.request
from collections import defaultdict, Counter

import numpy as np
import pandas as pd
from scipy.stats import hypergeom

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

WORKING_DIR = '/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV'
ENRICH_DIR  = os.path.join(WORKING_DIR, 'annotation_output', 'enrichment')
MAP_DIR     = os.path.join(ENRICH_DIR, '00_gene_id_mapping')
SETS_DIR    = os.path.join(ENRICH_DIR, '01_gene_sets')
RES_DIR     = os.path.join(ENRICH_DIR, '02_enrichment')
SUM_DIR     = os.path.join(ENRICH_DIR, '03_summary')
FIG_DIR     = os.path.join(ENRICH_DIR, '04_figures')
LOG_DIR     = os.path.join(ENRICH_DIR, 'logs')
ORA_IN_DIR  = os.path.join(WORKING_DIR, 'annotation_output', 'deg_comparisons',
                           'ora_gene_lists')

GATE, LAYER = 'FDR', 'pseudobulk'
GPROFILER_ORG, KEGG_ORG = 'mmulatta', 'mcc'

GP_GOST = 'https://biit.cs.ut.ee/gprofiler/api/gost/profile/'
GP_SOURCES = ['GO:BP', 'GO:MF', 'GO:CC', 'KEGG', 'REAC', 'WP']

# ---- statistics ----
ENRICH_ALPHA      = 0.05      # on the corrected p value
GP_CORRECTION     = 'g_SCS'   # 'g_SCS' | 'fdr' | 'bonferroni'
MIN_GENES_IN_LIST = 10        # per direction, below this an ORA is not run
MIN_SET_SIZE      = 5         # ignore tiny pathways
MAX_SET_SIZE      = 500       # ignore giant unspecific terms
MIN_OVERLAP       = 3         # a term needs this many of our genes to be reported

RUN_GPROFILER = True
RUN_LOCAL_KEGG = True
USE_CUSTOM_BACKGROUND = True  # the whole point; False falls back to whole-genome

HTTP_TIMEOUT, HTTP_RETRIES, HTTP_BACKOFF = 180, 3, 4.0
SLEEP_BETWEEN_CALLS = 0.4     # be polite to the public endpoint

N_TOP_REPORT = 10
FIG_DPI = 300
HEX = {'up': '#C4453C', 'down': '#2C6E9B', 'grid': '#D8D8D8'}
plt.rcParams.update({
    'font.size': 28, 'axes.titlesize': 34, 'axes.labelsize': 30,
    'xtick.labelsize': 24, 'ytick.labelsize': 24, 'legend.fontsize': 20,
    'pdf.fonttype': 42, 'ps.fonttype': 42,
})

for d in (RES_DIR, SUM_DIR, FIG_DIR, LOG_DIR):
    os.makedirs(d, exist_ok=True)
for lib in ('KEGG', 'GO_BP', 'GO_MF', 'GO_CC', 'REAC', 'WP', 'KEGG_local'):
    os.makedirs(os.path.join(RES_DIR, lib), exist_ok=True)
pd.set_option('display.width', 220)
pd.set_option('display.max_columns', 40)

_report = []


def log(msg=''):
    print(msg, flush=True)
    _report.append(str(msg))


def hr(t):
    log('')
    log('=' * 90)
    log(t)
    log('=' * 90)


def bh_fdr(p):
    p = np.asarray(p, dtype=float)
    out = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if not ok.any():
        return out
    pv = p[ok]
    n = pv.size
    o = np.argsort(pv)
    q = pv[o] * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    res = np.empty(n)
    res[o] = np.clip(q, 0, 1)
    out[ok] = res
    return out


def http_post_json(url, payload):
    data = json.dumps(payload).encode()
    last = ''
    for a in range(HTTP_RETRIES):
        try:
            req = urllib.request.Request(
                url, data=data,
                headers={'Content-Type': 'application/json',
                         'User-Agent': 'rhesus-shiv-enrichment/1.0'})
            with urllib.request.urlopen(req, timeout=HTTP_TIMEOUT) as r:
                return json.loads(r.read().decode()), ''
        except Exception as e:
            last = f'{type(e).__name__}: {e}'
            if a < HTTP_RETRIES - 1:
                time.sleep(HTTP_BACKOFF * (a + 1))
    return None, last


def slug(s):
    out = str(s)
    for a, b in [(' ', '_'), ('/', '-'), ('+', 'pos'), (',', ''), ('(', ''),
                 (')', ''), (':', ''), ('.', '')]:
        out = out.replace(a, b)
    return out


# %% ============================================================
# CELL 2 - LOAD INPUTS FROM STEP 1
# ============================================================
hr('CELL 2 - load gene lists, id map and gene sets')

map_path = os.path.join(MAP_DIR, 'gene_id_map.csv')
if not os.path.exists(map_path):
    raise FileNotFoundError(f'{map_path} not found. Run shiv_enrichment_prep.py first.')
id_map = pd.read_csv(map_path)
# --- identifier normalisation ---------------------------------------------
# BUG FIXED 2026-09-28, and it is the reason the local KEGG engine returned
# local=0 on all 35 lists while g:Profiler's KEGG found 51 terms.
#
# gene_id_map.csv has blanks in its entrez column, so pandas types that column
# float64. str(694372.0) is '694372.0'. The GMT holds '694372'. The two never
# match, every pathway overlap is zero, no p value is ever computed, and the
# engine reports nothing at all rather than erroring. The self-test could not
# catch it because it builds its own ids and never touches either file.
#
# Normalise both sides, and then ASSERT that they intersect, so a future
# namespace change fails loudly instead of silently returning zero.
def norm_id(x):
    s = str(x).strip()
    if s.endswith('.0') and s[:-2].isdigit():
        s = s[:-2]
    for pre in ('mcc:', 'ncbi-geneid:', 'hsa:', 'ncbi-proteinid:'):
        if s.startswith(pre):
            s = s[len(pre):]
    return s


sym2entrez = {str(r.symbol).strip(): norm_id(r.entrez)
              for r in id_map.itertuples()
              if pd.notna(getattr(r, 'entrez', None))}
log(f'  id map: {len(id_map):,} symbols, {len(sym2entrez):,} with an NCBI GeneID')

pat_list = re.compile(rf'^(C\d+)_(.+)_{LAYER}_{GATE}_(up|down)\.csv$')
pat_bg = re.compile(rf'^(C\d+)_(.+)_{LAYER}_background\.csv$')
lists, backgrounds = [], {}
for fn in sorted(os.listdir(ORA_IN_DIR)):
    m = pat_list.match(fn)
    if m:
        cid, stem, direction = m.groups()
        g = pd.read_csv(os.path.join(ORA_IN_DIR, fn))
        lists.append({'contrast_id': cid, 'stem': stem, 'direction': direction,
                      'genes': g['gene'].astype(str).tolist()})
        continue
    m = pat_bg.match(fn)
    if m:
        backgrounds[m.group(1)] = pd.read_csv(
            os.path.join(ORA_IN_DIR, fn))['gene'].astype(str).tolist()
log(f'  gene lists: {len(lists)}   backgrounds: {len(backgrounds)}')

kegg_sets, kegg_names = {}, {}
gmt = os.path.join(SETS_DIR, f'KEGG_{KEGG_ORG}.gmt')
if RUN_LOCAL_KEGG and os.path.exists(gmt):
    with open(gmt) as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 3:
                kegg_sets[parts[0]] = {norm_id(g) for g in parts[2:]}
                kegg_names[parts[0]] = parts[1]
    log(f'  local KEGG GMT: {len(kegg_sets):,} pathways')
    _gmt_ids = set().union(*kegg_sets.values()) if kegg_sets else set()
    _shared = set(sym2entrez.values()) & _gmt_ids
    log(f'  id namespaces : {len(_shared):,} ids shared between the gene map '
        f'and the GMT')
    log(f'    map sample  : {sorted(set(sym2entrez.values()))[:3]}')
    log(f'    GMT sample  : {sorted(_gmt_ids)[:3]}')
    if not _shared:
        raise SystemExit(
            '[FATAL] the gene map and the KEGG GMT share NO identifiers, so '
            'every\n  pathway overlap would be zero and the local engine '
            'would silently\n  report nothing. Fix norm_id() for the '
            'namespaces printed above.')
    if len(_shared) < 0.25 * len(_gmt_ids):
        log(f'  [WARN] only {100*len(_shared)/max(1,len(_gmt_ids)):.1f}% of GMT '
            f'ids are reachable from the gene map.')
else:
    if RUN_LOCAL_KEGG:
        log(f'  WARN: {gmt} missing; local KEGG cross-check disabled')
    RUN_LOCAL_KEGG = False


# %% ============================================================
# CELL 3 - LOCAL HYPERGEOMETRIC, AND ITS SELF-TEST
# ============================================================
hr('CELL 3 - local over-representation engine + self-test')


def ora_hypergeometric(query, background, gene_sets, set_names=None):
    """One-sided hypergeometric over-representation.

    M = background size, n = genes in the set AND in the background,
    N = query size (restricted to the background), k = overlap.
    p = P(X >= k) = hypergeom.sf(k-1, M, n, N).

    The set is intersected with the background FIRST. Skipping that is the
    classic error: it counts pathway genes that were never testable in this
    experiment and makes every p value too small.
    """
    bg = set(background)
    q = set(query) & bg
    M, N = len(bg), len(q)
    if M == 0 or N == 0:
        return pd.DataFrame()
    rows = []
    for sid, genes in gene_sets.items():
        gset = set(genes) & bg
        n = len(gset)
        if n < MIN_SET_SIZE or n > MAX_SET_SIZE:
            continue
        hits = q & gset
        k = len(hits)
        if k < MIN_OVERLAP:
            continue
        p = float(hypergeom.sf(k - 1, M, n, N))
        rows.append({
            'term_id': sid,
            'term_name': (set_names or {}).get(sid, sid),
            'p_value': p,
            'overlap': k, 'term_size': n, 'query_size': N, 'background_size': M,
            'fold_enrichment': (k / N) / (n / M) if n and N else np.nan,
            'genes': ';'.join(sorted(hits)),
        })
    if not rows:
        return pd.DataFrame()
    out = pd.DataFrame(rows)
    out['p_adj'] = bh_fdr(out['p_value'].values)
    return out.sort_values('p_value').reset_index(drop=True)


# --- self-test with a hand-checkable answer -------------------------------
_bg = [f'g{i}' for i in range(1000)]
_set = {f'g{i}' for i in range(50)}                       # 50-gene pathway
_query = [f'g{i}' for i in range(20)] + [f'g{i}' for i in range(900, 980)]
_res = ora_hypergeometric(_query, _bg, {'P1': _set}, {'P1': 'test pathway'})
_expect_p = float(hypergeom.sf(20 - 1, 1000, 50, 100))
_got_p = float(_res.loc[0, 'p_value'])
log(f'  self-test: 1000-gene background, 50-gene pathway, 100-gene query,')
log(f'             20 of the query inside the pathway')
log(f'    expected p = {_expect_p:.3e}')
log(f'    computed p = {_got_p:.3e}')
assert abs(_got_p - _expect_p) < 1e-12, 'hypergeometric self-test FAILED'
assert int(_res.loc[0, 'overlap']) == 20
assert int(_res.loc[0, 'term_size']) == 50
assert int(_res.loc[0, 'query_size']) == 100
log('    PASS (overlap, term size, query size and p all correct)')

# background must actually shrink the set: verify the intersection happens
_res2 = ora_hypergeometric(_query, [f'g{i}' for i in range(500)],
                           {'P1': _set}, {'P1': 'test'})
log(f'  background sensitivity: halving the background changes p from '
    f'{_got_p:.2e} to {float(_res2.loc[0, "p_value"]):.2e}')
assert float(_res2.loc[0, 'p_value']) != _got_p, 'background is being ignored'
log('    PASS (background is respected, not silently dropped)')


# %% ============================================================
# CELL 4 - g:PROFILER ENGINE
# ============================================================
hr('CELL 4 - g:Profiler over-representation')


# --- how g:Profiler actually reports the hit genes -------------------------
# BUG FIXED 2026-09-28. The first version did
#       ';'.join(t.get('intersections', []) or [])
# on the assumption that `intersections` is a list of gene names. It is not.
# With no_evidences = False, g:Profiler returns `intersections` as a list of
# per-gene EVIDENCE CODE LISTS, e.g.
#       [[], ['IEA'], [], ['IDA','IEA'], ...]
# positionally aligned to the RESOLVED query genes, which the response carries
# separately at meta.genes_metadata.query.<query_name>.ensgs. Position i is
# non-empty if and only if resolved query gene i belongs to that term. Joining
# it directly raised
#       TypeError: sequence item 0: expected str instance, list found
# which is the error you hit, and it fired on the very first term of the very
# first contrast.
#
# The reason this needs care rather than a one-line patch: the positional
# alignment is an ASSUMPTION about an API I cannot reach from here to confirm.
# If it is wrong, the failure mode is silent and nasty, because every pathway
# would be annotated with a plausible but WRONG set of driver genes, and that
# is exactly the column someone reads to decide what the pathway means. So the
# extractor below validates the alignment against intersection_size, which
# g:Profiler computes independently, and returns nothing rather than guessing
# when the two disagree. An empty genes column is recoverable. A confidently
# wrong one is not.

_GP_SHAPE_LOGGED = {'done': False}
_GP_GENE_WARN = {'n': 0, 'reasons': Counter()}


def _resolved_query_genes(resp):
    """The query genes g:Profiler actually resolved, in the order the
    `intersections` positions refer to. Returns (ids, symbol_map) where
    symbol_map turns a resolved id back into the symbol we submitted, when the
    response gives us enough to do that."""
    meta = resp.get('meta', {}) or {}
    gm = (meta.get('genes_metadata', {}) or {}).get('query', {}) or {}
    if not isinstance(gm, dict) or not gm:
        return None, {}
    first = next(iter(gm.values()))
    if not isinstance(first, dict):
        return None, {}
    ids = first.get('ensgs')

    # `mapping` is {submitted_symbol: [resolved_id, ...]}. Invert it so the
    # genes column comes back in OUR namespace (rhesus symbols) rather than as
    # ENSMMUG accessions, which are not what Zoey wants to read.
    sym = {}
    mapping = first.get('mapping')
    if isinstance(mapping, dict):
        for submitted, resolved in mapping.items():
            if isinstance(resolved, str):
                resolved = [resolved]
            for rid in (resolved or []):
                sym.setdefault(rid, submitted)
    return ids, sym


def _hit_genes(term, qids, symmap):
    """Hit gene symbols for one term, or '' when the shape cannot be trusted."""
    inter = term.get('intersections')
    if not inter:
        return ''

    # Shape B, defensive: a plain list of gene-name strings. Older and some
    # alternate g:Profiler responses do this. Join it directly.
    if all(isinstance(x, str) for x in inter):
        return ';'.join(inter)

    if not all(isinstance(x, list) for x in inter):
        _GP_GENE_WARN['reasons']['mixed_types'] += 1
        return ''

    # Shape A: per-gene evidence code lists, positionally aligned.
    if qids is None:
        _GP_GENE_WARN['reasons']['no_resolved_gene_list'] += 1
        return ''
    if len(inter) != len(qids):
        _GP_GENE_WARN['reasons']['length_mismatch'] += 1
        return ''

    hits = [qids[i] for i, ev in enumerate(inter) if ev]

    # Independent confirmation. g:Profiler computed intersection_size itself,
    # without reference to the positional layout, so if our reading of the
    # layout is right these two MUST agree.
    isize = term.get('intersection_size')
    if isize is not None and len(hits) != isize:
        _GP_GENE_WARN['reasons']['overlap_disagrees'] += 1
        return ''

    return ';'.join(symmap.get(g, g) for g in hits)


def _log_shape_once(resp):
    """Write the real structure of the first response to disk. I could not
    reach the API from my side, so this is how the actual shape gets recorded
    rather than inferred."""
    if _GP_SHAPE_LOGGED['done']:
        return
    _GP_SHAPE_LOGGED['done'] = True
    try:
        res = (resp.get('result') or [])[:1]
        qids, symmap = _resolved_query_genes(resp)
        lines = ['g:Profiler response shape, first successful call', '']
        lines.append(f'top-level keys      : {sorted(resp.keys())}')
        lines.append(f'resolved query ids  : '
                     f'{len(qids) if qids is not None else "NOT FOUND"}')
        lines.append(f'  first 5           : {(qids or [])[:5]}')
        lines.append(f'symbol back-map     : {len(symmap)} entries')
        lines.append(f'  sample            : {list(symmap.items())[:5]}')
        if res:
            t = res[0]
            lines.append('')
            lines.append(f'first term keys     : {sorted(t.keys())}')
            lines.append(f'  native            : {t.get("native")}')
            lines.append(f'  intersection_size : {t.get("intersection_size")}')
            inter = t.get('intersections')
            lines.append(f'  intersections len : '
                         f'{len(inter) if inter is not None else None}')
            if inter:
                lines.append(f'  element types     : '
                             f'{sorted({type(x).__name__ for x in inter})}')
                lines.append(f'  first 8 elements  : {inter[:8]}')
                nz = sum(1 for x in inter if x)
                lines.append(f'  non-empty count   : {nz}  '
                             f'(should equal intersection_size)')
        p = os.path.join(LOG_DIR, 'gprofiler_response_shape.txt')
        with open(p, 'w') as fh:
            fh.write('\n'.join(str(x) for x in lines))
        log(f'  [shape probe] wrote {p}')
    except Exception as e:                      # never let the probe kill a run
        log(f'  [shape probe] failed harmlessly: {e}')


def run_gprofiler(query, background):
    payload = {
        'organism': GPROFILER_ORG,
        'query': list(query),
        'sources': GP_SOURCES,
        'user_threshold': ENRICH_ALPHA,
        'significance_threshold_method': GP_CORRECTION,
        'no_evidences': False,
        'all_results': False,
    }
    if USE_CUSTOM_BACKGROUND and background:
        payload['domain_scope'] = 'custom'
        payload['background'] = list(background)
    r, err = http_post_json(GP_GOST, payload)
    if r is None:
        return None, err
    _log_shape_once(r)
    res = r.get('result', [])
    if not res:
        return pd.DataFrame(), ''

    qids, symmap = _resolved_query_genes(r)

    rows = []
    for t in res:
        g = _hit_genes(t, qids, symmap)
        if not g and t.get('intersection_size'):
            _GP_GENE_WARN['n'] += 1
        rows.append({
            'source': t.get('source'),
            'term_id': t.get('native'),
            'term_name': t.get('name'),
            'p_adj': t.get('p_value'),      # already corrected by g:Profiler
            'overlap': t.get('intersection_size'),
            'term_size': t.get('term_size'),
            'query_size': t.get('query_size'),
            'background_size': t.get('effective_domain_size'),
            'precision': t.get('precision'),
            'recall': t.get('recall'),
            'genes': g,
        })
    df = pd.DataFrame(rows)
    df = df[(df['term_size'] >= MIN_SET_SIZE) &
            (df['term_size'] <= MAX_SET_SIZE) &
            (df['overlap'] >= MIN_OVERLAP)]
    return df.sort_values('p_adj').reset_index(drop=True), ''


# %% ============================================================
# CELL 5 - RUN EVERY LIST
# ============================================================
hr('CELL 5 - run enrichment for every contrast and direction')

SRC_DIR = {'GO:BP': 'GO_BP', 'GO:MF': 'GO_MF', 'GO:CC': 'GO_CC',
           'KEGG': 'KEGG', 'REAC': 'REAC', 'WP': 'WP'}
all_rows, run_log = [], []

for i, item in enumerate(lists, 1):
    cid, stem, direction = item['contrast_id'], item['stem'], item['direction']
    genes = item['genes']
    bg = backgrounds.get(cid, [])
    tag = f'{cid}_{slug(stem)}_{LAYER}_{GATE}_{direction}'

    rec = {'contrast_id': cid, 'stem': stem, 'direction': direction,
           'n_genes': len(genes), 'n_background': len(bg)}

    if len(genes) < MIN_GENES_IN_LIST:
        rec.update(status='skipped_small_list', n_terms_gp=0, n_terms_local=0)
        run_log.append(rec)
        continue
    if not bg:
        rec.update(status='skipped_no_background', n_terms_gp=0, n_terms_local=0)
        run_log.append(rec)
        continue

    n_gp = n_loc = 0

    # --- engine A: g:Profiler (symbols, Ensembl-backed) ---
    if RUN_GPROFILER:
        gp_df, err = run_gprofiler(genes, bg)
        time.sleep(SLEEP_BETWEEN_CALLS)
        if gp_df is None:
            rec['gprofiler_error'] = err
            log(f'  [{i}/{len(lists)}] {tag}: g:Profiler FAILED ({err})')
        elif len(gp_df):
            for src, sub in gp_df.groupby('source'):
                d = SRC_DIR.get(src)
                if not d:
                    continue
                sub.to_csv(os.path.join(RES_DIR, d, f'{tag}.tsv'),
                           sep='\t', index=False)
            gp_df['engine'] = 'gprofiler'
            gp_df['contrast_id'] = cid
            gp_df['stem'] = stem
            gp_df['direction'] = direction
            all_rows.append(gp_df)
            n_gp = len(gp_df)

    # --- engine B: local hypergeometric on KEGG mcc (Entrez) ---
    if RUN_LOCAL_KEGG:
        q_e = [sym2entrez[g] for g in genes if g in sym2entrez]
        b_e = [sym2entrez[g] for g in bg if g in sym2entrez]
        if len(q_e) >= MIN_GENES_IN_LIST and b_e:
            loc = ora_hypergeometric(q_e, b_e, kegg_sets, kegg_names)
            if len(loc):
                loc = loc[loc['p_adj'] < ENRICH_ALPHA]
            if len(loc):
                loc.to_csv(os.path.join(RES_DIR, 'KEGG_local', f'{tag}.tsv'),
                           sep='\t', index=False)
                loc['engine'] = 'kegg_local'
                loc['source'] = 'KEGG_local'
                loc['contrast_id'] = cid
                loc['stem'] = stem
                loc['direction'] = direction
                all_rows.append(loc)
                n_loc = len(loc)
        rec['n_genes_entrez'] = len(q_e)

    rec.update(status='ok', n_terms_gp=n_gp, n_terms_local=n_loc)
    run_log.append(rec)
    if i % 10 == 0 or i == len(lists):
        log(f'  [{i}/{len(lists)}] {tag}: gp={n_gp} local={n_loc}')

runs = pd.DataFrame(run_log)
runs.to_csv(os.path.join(SUM_DIR, 'enrichment_run_log.csv'), index=False)
log('\n  run status:')
log(runs['status'].value_counts().to_string())

# --- did the hit-gene extraction hold up? ---
# The p values, term IDs and overlap counts come straight from g:Profiler and
# are unaffected by any of this. Only the `genes` column depends on our reading
# of the positional layout, so report it separately and honestly.
log('\n  g:Profiler hit-gene extraction:')
if _GP_GENE_WARN['n'] == 0:
    log('    [OK] every term with a non-zero overlap got its gene list, and')
    log('         each one was confirmed against g:Profiler\'s own')
    log('         intersection_size.')
else:
    log(f'    [WARN] {_GP_GENE_WARN["n"]} term(s) have an overlap but an EMPTY')
    log('           genes column. Deliberate: the layout did not validate, so')
    log('           the driver genes were withheld rather than guessed.')
    for k, v in _GP_GENE_WARN['reasons'].most_common():
        log(f'             {k}: {v}')
    log('           Enrichment p values and overlaps are UNAFFECTED.')
    log(f'           See {os.path.join(LOG_DIR, "gprofiler_response_shape.txt")}')
    log('           for the real response structure, and send it to me.')


# %% ============================================================
# CELL 6 - SUMMARY TABLES
# ============================================================
hr('CELL 6 - summary')

if all_rows:
    res = pd.concat(all_rows, ignore_index=True)
    keep = ['contrast_id', 'stem', 'direction', 'engine', 'source', 'term_id',
            'term_name', 'p_adj', 'overlap', 'term_size', 'query_size',
            'background_size', 'genes']
    res = res[[c for c in keep if c in res.columns]]
    res.to_csv(os.path.join(SUM_DIR, 'enrichment_all_results.csv'), index=False)
    sig = res[res['p_adj'] < ENRICH_ALPHA]
    sig.to_csv(os.path.join(SUM_DIR, 'enrichment_significant.csv'), index=False)
    log(f'  terms returned   : {len(res):,}')
    log(f'  significant terms: {len(sig):,}  (corrected p < {ENRICH_ALPHA})')

    log('\n  significant terms by source and direction:')
    log(pd.crosstab(sig['source'], sig['direction']).to_string())

    bypop = (sig.groupby(['stem', 'direction']).size()
             .rename('n_terms').reset_index()
             .pivot(index='stem', columns='direction', values='n_terms')
             .fillna(0).astype(int))
    log('\n  significant terms per contrast and direction (top 20):')
    log(bypop.assign(total=bypop.sum(axis=1)).nlargest(20, 'total').to_string())
    bypop.to_csv(os.path.join(SUM_DIR, 'terms_by_contrast_matrix.csv'))

    top = (sig.sort_values('p_adj')
           .groupby(['contrast_id', 'direction', 'source'])
           .head(N_TOP_REPORT))
    top.to_csv(os.path.join(SUM_DIR, 'top_terms_per_contrast.csv'), index=False)
    log(f'\n  wrote top_terms_per_contrast.csv ({len(top):,} rows)')

    # KEGG agreement between the two engines, where both ran
    if {'gprofiler', 'kegg_local'} <= set(res['engine']):
        a = set(zip(sig[sig['source'] == 'KEGG']['contrast_id'],
                    sig[sig['source'] == 'KEGG']['direction'],
                    sig[sig['source'] == 'KEGG']['term_name']))
        b = set(zip(sig[sig['source'] == 'KEGG_local']['contrast_id'],
                    sig[sig['source'] == 'KEGG_local']['direction'],
                    sig[sig['source'] == 'KEGG_local']['term_name']))
        if a or b:
            log(f'\n  KEGG engine agreement: {len(a & b)} shared, '
                f'{len(a - b)} g:Profiler only, {len(b - a)} local only')
            log('  Disagreement is expected in part: g:Profiler is Ensembl-backed')
            log('  and the local GMT is KEGG-backed, so gene-to-pathway membership')
            log('  differs. Treat shared terms as the confident set.')

    try:
        fig, ax = plt.subplots(figsize=(18, 10))
        cnt = sig.groupby(['source', 'direction']).size().unstack(fill_value=0)
        idx = np.arange(len(cnt))
        w = 0.38
        for j, (d, col) in enumerate([('up', HEX['up']), ('down', HEX['down'])]):
            if d in cnt.columns:
                ax.bar(idx + (j - 0.5) * w, cnt[d].values, w, label=d, color=col,
                       edgecolor='white', linewidth=2)
        ax.set_xticks(idx)
        ax.set_xticklabels(cnt.index, rotation=30, ha='right')
        ax.set_ylabel('significant terms', fontsize=30)
        ax.set_title(f'Enrichment yield, corrected p < {ENRICH_ALPHA}', fontsize=34)
        ax.legend(fontsize=22, frameon=False)
        for sp in ('top', 'right'):
            ax.spines[sp].set_visible(False)
        for ext in ('pdf', 'png'):
            fig.savefig(os.path.join(FIG_DIR, f'enrichment_yield.{ext}'),
                        dpi=FIG_DPI, bbox_inches='tight')
        plt.close(fig)
        log('  saved enrichment_yield.pdf / .png')
    except Exception as e:
        log(f'  WARN: figure failed ({e})')
else:
    res = sig = pd.DataFrame()
    log('  No enrichment results at all. Check the run log and the mapping rate')
    log('  from step 1 before concluding this is biological.')

log(f"""
  SUMMARY
  -------
  gene lists processed : {len(lists)}
  ran successfully     : {int((runs['status'] == 'ok').sum()) if len(runs) else 0}
  skipped, small list  : {int((runs['status'] == 'skipped_small_list').sum()) if len(runs) else 0}
  skipped, no background: {int((runs['status'] == 'skipped_no_background').sum()) if len(runs) else 0}
  significant terms    : {len(sig):,}

  Everything is in 03_summary/. enrichment_significant.csv is the one to read
  first; top_terms_per_contrast.csv is the one to put in front of Zoey.
""")
with open(os.path.join(LOG_DIR, 'enrichment_run_report.txt'), 'w') as fh:
    fh.write('\n'.join(_report))
log('=' * 90)
log('ENRICHMENT COMPLETE')
log('=' * 90)
