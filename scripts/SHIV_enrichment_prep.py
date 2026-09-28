#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Rhesus SHIV - Enrichment Prep: connectivity probe, gene ID mapping, KEGG gene sets
===================================================================================
Step 1 of 2. Does NO enrichment. It answers the question that decides whether the
enrichment is worth running at all: how many of our rhesus gene symbols can the
annotation databases actually resolve.

WHY NOT gseapy.enrichr()
------------------------
Enrichr accepts human and mouse gene lists only. Its companion modEnrichr adds
fly, worm, yeast and zebrafish, but no macaque. So the enrichr() path in the old
KEGG_On_DEGs.py template cannot be reused here. gseapy is still useful, but only
its OFFLINE enrich() against a local GMT, which is what CELL 5 builds.

NHP-NATIVE ROUTES, so no human-ortholog conversion is needed
-------------------------------------------------------------
  KEGG      organism code 'mcc' = Macaca mulatta. KEGG maintains rhesus pathway
            maps natively. CELL 5 pulls them over the KEGG REST API and writes a
            GMT keyed on NCBI GeneIDs.
  g:Profiler organism 'mmulatta' = Macaca mulatta, served from Ensembl. Covers
            GO:BP, GO:MF, GO:CC, KEGG, Reactome and WikiPathways in one call, and
            accepts a custom background, which we need.
  (For the record, R users have org.Mmu.eg.db on Bioconductor and
   clusterProfiler::enrichKEGG(organism = "mcc"). Same annotations, different
   language. We stay in Python.)

THE REAL RISK IS SYMBOL MAPPING, NOT TOOL SUPPORT
--------------------------------------------------
The Mmul10 annotation in this object carries several symbol classes that are not
standard HGNC-style gene names, and every one of them will silently drop out of
an enrichment unless we look:
    KEG06_*        mitochondrial locus tags (MT- equivalents), ~23 genes
    MAMU-*         MHC class I/II, the rhesus equivalent of HLA-
    LOC*           NCBI GeneIDs with no assigned symbol
    *_1            duplicate-resolved symbols, e.g. RPS27A_1, FCGRT_1, CD9_1
    C##H##orf##    macaque ortholog naming for human Cxorf genes
CELL 4 classifies every unmapped symbol into these buckets and reports the rate
per gene list, so a low-yield enrichment is explained rather than mysterious.

CONNECTIVITY
------------
CELL 2 probes every endpoint before anything else and prints exactly what came
back. If titan cannot reach biit.cs.ut.ee or rest.kegg.jp, you find out in the
first ten seconds rather than after the mapping work. Both hosts are plain HTTPS
with no key required.

Read-only w.r.t. the DEG output. Spyder cells (# %%).

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
import glob
import urllib.request
import urllib.error
from collections import Counter, defaultdict

import numpy as np
import pandas as pd

WORKING_DIR = '/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV'
DEG_DIR     = os.path.join(WORKING_DIR, 'annotation_output', 'deg_comparisons')
ORA_IN_DIR  = os.path.join(DEG_DIR, 'ora_gene_lists')

OUT_DIR   = os.path.join(WORKING_DIR, 'annotation_output', 'enrichment')
MAP_DIR   = os.path.join(OUT_DIR, '00_gene_id_mapping')
SETS_DIR  = os.path.join(OUT_DIR, '01_gene_sets')
LOG_DIR   = os.path.join(OUT_DIR, 'logs')

# ---- organism ----
GPROFILER_ORG = 'mmulatta'     # Macaca mulatta in g:Profiler (Ensembl-backed)
KEGG_ORG      = 'mcc'          # Macaca mulatta in KEGG

# ---- endpoints ----
GP_BASE   = 'https://biit.cs.ut.ee/gprofiler/api'
GP_GOST   = f'{GP_BASE}/gost/profile/'
GP_CONV   = f'{GP_BASE}/convert/convert/'
KEGG_BASE = 'https://rest.kegg.jp'

# ---- which gate's lists to prepare ----
# The DEG script writes both rawP_* and FDR_* lists. Zoey and Jake agreed on
# FDR < 0.05 with |log2FC| > 1.5 (2026-09-28), so FDR is the one we prepare.
GATE = 'FDR'                   # 'FDR' or 'rawP'
LAYER = 'pseudobulk'           # 'pseudobulk' (per-animal, credible) or 'wilcoxon'

# ---- network behaviour ----
HTTP_TIMEOUT   = 120
HTTP_RETRIES   = 3
HTTP_BACKOFF   = 3.0
CONVERT_CHUNK  = 2000          # g:Convert query size per request
CACHE_OK       = True          # reuse a previous pull instead of re-hitting the API

# ---- unmapped symbol classification ----
SYMBOL_CLASSES = [
    ('mitochondrial_locus_tag', re.compile(r'^KEG06_')),
    ('MHC_MAMU',                re.compile(r'^MAMU-')),
    ('LOC_no_symbol',           re.compile(r'^LOC\d+$')),
    ('duplicate_suffix',        re.compile(r'_\d+$')),
    ('macaque_orf_naming',      re.compile(r'^C\d+H\d+orf\d+$', re.I)),
    ('ensembl_id',              re.compile(r'^ENSMMUG\d+$')),
]

for d in (OUT_DIR, MAP_DIR, SETS_DIR, LOG_DIR):
    os.makedirs(d, exist_ok=True)
pd.set_option('display.width', 200)
pd.set_option('display.max_columns', 40)
pd.set_option('display.max_rows', 200)

_report = []


def log(msg=''):
    print(msg, flush=True)
    _report.append(str(msg))


def hr(t):
    log('')
    log('=' * 90)
    log(t)
    log('=' * 90)


def http_post_json(url, payload):
    """POST JSON with retries. Returns (parsed, error_string)."""
    data = json.dumps(payload).encode()
    last = ''
    for attempt in range(HTTP_RETRIES):
        try:
            req = urllib.request.Request(
                url, data=data,
                headers={'Content-Type': 'application/json',
                         'User-Agent': 'rhesus-shiv-enrichment/1.0'})
            with urllib.request.urlopen(req, timeout=HTTP_TIMEOUT) as r:
                return json.loads(r.read().decode()), ''
        except Exception as e:
            last = f'{type(e).__name__}: {e}'
            if attempt < HTTP_RETRIES - 1:
                time.sleep(HTTP_BACKOFF * (attempt + 1))
    return None, last


def http_get_text(url):
    last = ''
    for attempt in range(HTTP_RETRIES):
        try:
            req = urllib.request.Request(
                url, headers={'User-Agent': 'rhesus-shiv-enrichment/1.0'})
            with urllib.request.urlopen(req, timeout=HTTP_TIMEOUT) as r:
                return r.read().decode(), ''
        except Exception as e:
            last = f'{type(e).__name__}: {e}'
            if attempt < HTTP_RETRIES - 1:
                time.sleep(HTTP_BACKOFF * (attempt + 1))
    return None, last


# %% ============================================================
# CELL 2 - CONNECTIVITY PROBE (runs first, deliberately)
# ============================================================
hr('CELL 2 - connectivity probe')
log(f'  g:Profiler organism : {GPROFILER_ORG}')
log(f'  KEGG organism       : {KEGG_ORG}')

probe_rows = []

# --- g:Profiler gost, with a small known-good rhesus gene list ---
_probe_genes = ['LYZ', 'CD68', 'C1QA', 'C1QB', 'C1QC', 'APOE', 'CTSB', 'CTSD',
                'TYROBP', 'FCER1G', 'CSF1R', 'AIF1', 'PSAP', 'CST3', 'MRC1']
gost_ok = False
r, err = http_post_json(GP_GOST, {
    'organism': GPROFILER_ORG, 'query': _probe_genes,
    'sources': ['GO:BP', 'KEGG'], 'user_threshold': 0.05,
    'significance_threshold_method': 'g_SCS', 'no_evidences': True})
if r is None:
    log(f'  g:Profiler gost     : UNREACHABLE  ({err})')
    probe_rows.append({'endpoint': 'gprofiler_gost', 'ok': False, 'detail': err})
else:
    res = r.get('result', [])
    gost_ok = True
    srcs = sorted({t.get('source') for t in res})
    log(f'  g:Profiler gost     : OK, {len(res)} terms on a 15-gene probe, '
        f'sources {srcs}')
    for t in res[:3]:
        log(f'                        {t.get("source")} {t.get("native")} '
            f'p={t.get("p_value"):.2e} {str(t.get("name"))[:48]}')
    probe_rows.append({'endpoint': 'gprofiler_gost', 'ok': True,
                       'detail': f'{len(res)} terms'})

# --- g:Convert ---
conv_ok = False
r2, err2 = http_post_json(GP_CONV, {
    'organism': GPROFILER_ORG, 'query': ['LYZ', 'CD68'],
    'target': 'ENTREZGENE_ACC'})
if r2 is None:
    log(f'  g:Profiler convert  : UNREACHABLE  ({err2})')
    probe_rows.append({'endpoint': 'gprofiler_convert', 'ok': False, 'detail': err2})
else:
    conv_ok = True
    log(f'  g:Profiler convert  : OK, {len(r2.get("result", []))} rows returned')
    probe_rows.append({'endpoint': 'gprofiler_convert', 'ok': True, 'detail': 'ok'})

# --- KEGG REST ---
kegg_ok = True
for label, path in [('pathway list', f'/list/pathway/{KEGG_ORG}'),
                    ('gene-pathway link', f'/link/pathway/{KEGG_ORG}'),
                    ('geneid conv', f'/conv/ncbi-geneid/{KEGG_ORG}')]:
    txt, e = http_get_text(KEGG_BASE + path)
    if txt is None:
        kegg_ok = False
        log(f'  KEGG {label:18s}: UNREACHABLE  ({e})')
        probe_rows.append({'endpoint': f'kegg_{label}', 'ok': False, 'detail': e})
    else:
        n = len([l for l in txt.strip().split('\n') if l])
        log(f'  KEGG {label:18s}: OK, {n:,} lines')
        probe_rows.append({'endpoint': f'kegg_{label}', 'ok': True, 'detail': f'{n} lines'})

pd.DataFrame(probe_rows).to_csv(os.path.join(LOG_DIR, 'connectivity_probe.csv'),
                                index=False)

if not (gost_ok or kegg_ok):
    log('\n  *** Neither g:Profiler nor KEGG is reachable from this host.')
    log('      Both are plain HTTPS on port 443 with no API key. If titan is')
    log('      behind a proxy, set https_proxy before running. Nothing below')
    log('      can work without at least one of them, so stopping here.')
    raise SystemExit(1)
if not gost_ok:
    log('\n  WARN: g:Profiler unreachable. KEGG-only enrichment is still possible')
    log('        from the local GMT, but there will be no GO results.')
if not kegg_ok:
    log('\n  WARN: KEGG REST unreachable. g:Profiler still serves KEGG terms for')
    log('        mmulatta, so enrichment can proceed without the local GMT, but')
    log('        the independent offline cross-check will be missing.')


# %% ============================================================
# CELL 3 - COLLECT THE GENE LISTS AND BACKGROUNDS
# ============================================================
hr('CELL 3 - collect DEG gene lists')

if not os.path.isdir(ORA_IN_DIR):
    raise FileNotFoundError(
        f'ORA gene lists not found at {ORA_IN_DIR}\n'
        'Run shiv_deg_comparisons.py with MANIFEST_ONLY = False first.')

pat_list = re.compile(rf'^(C\d+)_(.+)_{LAYER}_{GATE}_(up|down)\.csv$')
pat_bg   = re.compile(rf'^(C\d+)_(.+)_{LAYER}_background\.csv$')

lists, backgrounds = [], {}
for fn in sorted(os.listdir(ORA_IN_DIR)):
    m = pat_list.match(fn)
    if m:
        cid, stem, direction = m.groups()
        g = pd.read_csv(os.path.join(ORA_IN_DIR, fn))
        lists.append({'contrast_id': cid, 'stem': stem, 'direction': direction,
                      'file': fn, 'n_genes': len(g),
                      'genes': g['gene'].astype(str).tolist()})
        continue
    m = pat_bg.match(fn)
    if m:
        cid, stem = m.groups()
        g = pd.read_csv(os.path.join(ORA_IN_DIR, fn))
        backgrounds[cid] = g['gene'].astype(str).tolist()

log(f'  gate  : {GATE}   layer : {LAYER}')
log(f'  gene lists found     : {len(lists)}')
log(f'  backgrounds found    : {len(backgrounds)}')
if not lists:
    # Say WHICH of the three possible causes it is, rather than listing them.
    all_csv = [f for f in os.listdir(ORA_IN_DIR) if f.endswith('.csv')]
    tab_dir = os.path.join(DEG_DIR, 'tables')
    n_tab = len([f for f in os.listdir(tab_dir)
                 if f.endswith('.csv')]) if os.path.isdir(tab_dir) else 0

    log(f'  files in ora_gene_lists : {len(all_csv)}')
    log(f'  per-contrast DEG tables : {n_tab}')

    if not all_csv and n_tab:
        raise RuntimeError(
            f'ora_gene_lists is EMPTY but tables/ holds {n_tab} DEG tables.\n'
            '  The ORA-list writer never ran. It sits inside the per-contrast\n'
            '  test loop in the DEG script, and a MANIFEST_ONLY = True pass\n'
            '  creates the output directories and then exits before that loop.\n'
            '  Fix without recomputing anything:\n'
            '      python SHIV_ORA_list_backfill.py\n'
            '  It rebuilds the lists from tables/ in seconds.')
    if not all_csv and not n_tab:
        raise RuntimeError(
            'Neither ora_gene_lists nor tables/ has content. The DEG script has\n'
            '  not completed a real run. Set MANIFEST_ONLY = False and run it.')

    # There ARE files, just none matching this gate/layer. Show what is there.
    other = sorted({re.sub(r'^C\d+_.+?_(pseudobulk|wilcoxon)_',
                           r'\1_', f) for f in all_csv})[:12]
    raise RuntimeError(
        f'{len(all_csv)} file(s) in ora_gene_lists, but none match '
        f'GATE={GATE} LAYER={LAYER}.\n'
        f'  Suffixes present: {other}\n'
        '  Either set GATE/LAYER to one of those, or, if the threshold you want\n'
        '  produced no list, no contrast cleared it in both directions.')

ld = pd.DataFrame([{k: v for k, v in x.items() if k != 'genes'} for x in lists])
paired = (ld.groupby('contrast_id')['direction'].nunique() == 2)
log(f'  contrasts with BOTH directions: {int(paired.sum())} of '
    f'{ld["contrast_id"].nunique()}')
log('\n  gene list sizes:')
log(ld.groupby('direction')['n_genes'].describe()[['count', 'min', '50%', 'max']]
    .round(1).to_string())
ld.to_csv(os.path.join(MAP_DIR, 'gene_lists_inventory.csv'), index=False)

missing_bg = sorted({x['contrast_id'] for x in lists} - set(backgrounds))
if missing_bg:
    log(f'\n  *** {len(missing_bg)} contrast(s) have a gene list but NO background '
        f'file: {missing_bg[:6]}')
    log('      Over-representation without the right background inflates every')
    log('      p value. These will be skipped in step 2 unless the background is')
    log('      regenerated.')


# %% ============================================================
# CELL 4 - GENE ID MAPPING (the diagnostic that decides everything)
# ============================================================
hr('CELL 4 - map rhesus symbols to database identifiers')

all_symbols = sorted({g for x in lists for g in x['genes']} |
                     {g for bg in backgrounds.values() for g in bg})
log(f'  unique symbols across lists + backgrounds: {len(all_symbols):,}')

cache_path = os.path.join(MAP_DIR, 'gene_id_map.csv')
id_map = None
if CACHE_OK and os.path.exists(cache_path):
    cached = pd.read_csv(cache_path)
    if set(all_symbols) <= set(cached['symbol'].astype(str)):
        id_map = cached
        log(f'  reusing cached map ({len(cached):,} rows). Set CACHE_OK=False to refresh.')

if id_map is None and conv_ok:
    rows = []
    for i in range(0, len(all_symbols), CONVERT_CHUNK):
        chunk = all_symbols[i:i + CONVERT_CHUNK]
        r, err = http_post_json(GP_CONV, {
            'organism': GPROFILER_ORG, 'query': chunk,
            'target': 'ENTREZGENE_ACC', 'numeric_ns': 'ENTREZGENE_ACC'})
        if r is None:
            log(f'  WARN: convert chunk {i // CONVERT_CHUNK} failed ({err}); '
                'those symbols are marked unmapped')
            rows += [{'symbol': s, 'entrez': None, 'ensg': None, 'name': None}
                     for s in chunk]
            continue
        for x in r.get('result', []):
            conv = x.get('converted')
            rows.append({'symbol': x.get('incoming'),
                         'entrez': None if conv in (None, 'None', 'N/A') else str(conv),
                         'ensg': x.get('converted_alternative') or None,
                         'name': x.get('name')})
        log(f'  converted {min(i + CONVERT_CHUNK, len(all_symbols)):,} '
            f'of {len(all_symbols):,}')
    id_map = pd.DataFrame(rows).drop_duplicates('symbol')
    id_map.to_csv(cache_path, index=False)
    log(f'  wrote {os.path.basename(cache_path)}')
elif id_map is None:
    log('  g:Convert unreachable; building a symbol-only map (KEGG path will')
    log('  match on symbol via the GMT built in CELL 5).')
    id_map = pd.DataFrame({'symbol': all_symbols, 'entrez': None,
                           'ensg': None, 'name': None})

id_map['mapped'] = id_map['entrez'].notna()
n_map = int(id_map['mapped'].sum())
log(f'\n  MAPPING RATE: {n_map:,} of {len(id_map):,} symbols '
    f'({100 * n_map / max(len(id_map), 1):.1f}%) resolved to an NCBI GeneID')


def classify(sym):
    for name, rx in SYMBOL_CLASSES:
        if rx.search(sym):
            return name
    return 'standard_symbol_unmatched'


unmapped = id_map[~id_map['mapped']].copy()
if len(unmapped):
    unmapped['class'] = unmapped['symbol'].map(classify)
    log('\n  unmapped symbols by class:')
    cls = unmapped['class'].value_counts()
    log(cls.to_string())
    unmapped.sort_values(['class', 'symbol']).to_csv(
        os.path.join(MAP_DIR, 'unmapped_genes.csv'), index=False)
    log(f'  wrote unmapped_genes.csv ({len(unmapped):,} rows)')
    n_std = int(cls.get('standard_symbol_unmatched', 0))
    if n_std > 0.25 * len(id_map):
        log(f'\n  *** {n_std:,} ordinary-looking symbols failed to map. That is high')
        log('      enough to suspect a namespace mismatch rather than annotation')
        log('      gaps. Check a few by hand in g:Profiler before trusting the')
        log('      enrichment.')

# per gene list mapping rate, which is what actually limits each enrichment
mp = dict(zip(id_map['symbol'], id_map['mapped']))
rate_rows = []
for x in lists:
    n_ok = sum(1 for g in x['genes'] if mp.get(g, False))
    rate_rows.append({'contrast_id': x['contrast_id'], 'stem': x['stem'],
                      'direction': x['direction'], 'n_genes': x['n_genes'],
                      'n_mapped': n_ok,
                      'pct_mapped': round(100 * n_ok / max(x['n_genes'], 1), 1),
                      'n_background': len(backgrounds.get(x['contrast_id'], [])),
                      'n_background_mapped': sum(
                          1 for g in backgrounds.get(x['contrast_id'], [])
                          if mp.get(g, False))})
rate = pd.DataFrame(rate_rows)
rate['ora_viable'] = rate['n_mapped'] >= 10
rate.to_csv(os.path.join(MAP_DIR, 'mapping_rate_by_contrast.csv'), index=False)
log('\n  mapping rate per gene list (worst 12):')
log(rate.nsmallest(12, 'pct_mapped')[
    ['contrast_id', 'direction', 'n_genes', 'n_mapped', 'pct_mapped',
     'ora_viable']].to_string(index=False))
log(f'\n  gene lists still viable after mapping (>=10 genes): '
    f'{int(rate["ora_viable"].sum())} of {len(rate)}')
viable_both = (rate[rate['ora_viable']].groupby('contrast_id')['direction']
               .nunique() == 2)
log(f'  contrasts viable in BOTH directions: {int(viable_both.sum())}')


# %% ============================================================
# CELL 5 - BUILD THE KEGG mcc GENE SETS (local GMT)
# ============================================================
hr('CELL 5 - build KEGG gene sets for Macaca mulatta')

gmt_path = os.path.join(SETS_DIR, f'KEGG_{KEGG_ORG}.gmt')
if kegg_ok:
    txt_paths, _ = http_get_text(f'{KEGG_BASE}/list/pathway/{KEGG_ORG}')
    txt_links, _ = http_get_text(f'{KEGG_BASE}/link/pathway/{KEGG_ORG}')
    txt_conv, _ = http_get_text(f'{KEGG_BASE}/conv/ncbi-geneid/{KEGG_ORG}')

    pname = {}
    for line in (txt_paths or '').strip().split('\n'):
        if '\t' in line:
            pid, nm = line.split('\t', 1)
            pname[pid.replace('path:', '')] = re.sub(
                r'\s*-\s*Macaca mulatta.*$', '', nm).strip()

    # kegg gene id -> ncbi geneid
    kegg2ncbi = {}
    for line in (txt_conv or '').strip().split('\n'):
        if '\t' in line:
            k, n = line.split('\t', 1)
            kegg2ncbi[k.strip()] = n.strip().replace('ncbi-geneid:', '')

    sets = defaultdict(set)
    for line in (txt_links or '').strip().split('\n'):
        if '\t' not in line:
            continue
        gene, path = line.split('\t', 1)
        pid = path.replace('path:', '').strip()
        ncbi = kegg2ncbi.get(gene.strip())
        if ncbi:
            sets[pid].add(ncbi)

    with open(gmt_path, 'w') as fh:
        for pid, genes in sorted(sets.items()):
            fh.write('\t'.join([pid, pname.get(pid, pid)] + sorted(genes)) + '\n')
    sizes = [len(v) for v in sets.values()]
    log(f'  KEGG pathways for {KEGG_ORG}: {len(sets):,}')
    log(f'  genes per pathway: median {int(np.median(sizes))}, '
        f'range {min(sizes)}-{max(sizes)}')
    log(f'  distinct genes across all sets: {len(set().union(*sets.values())):,}')
    log(f'  wrote {os.path.basename(gmt_path)}')
    pd.DataFrame([{'pathway_id': k, 'name': pname.get(k, k), 'n_genes': len(v)}
                  for k, v in sorted(sets.items())]).to_csv(
        os.path.join(SETS_DIR, f'KEGG_{KEGG_ORG}_pathway_names.tsv'),
        sep='\t', index=False)

    # how much of OUR mapped background actually appears in KEGG
    kegg_genes = set().union(*sets.values()) if sets else set()
    our = set(id_map.loc[id_map['mapped'], 'entrez'].astype(str))
    log(f'\n  our mapped genes that appear in any KEGG {KEGG_ORG} pathway: '
        f'{len(our & kegg_genes):,} of {len(our):,} '
        f'({100 * len(our & kegg_genes) / max(len(our), 1):.1f}%)')
    log('  That fraction is the ceiling on what a KEGG over-representation can see.')
else:
    log('  KEGG REST unreachable; no local GMT written. g:Profiler still supplies')
    log('  KEGG terms in step 2.')


# %% ============================================================
# CELL 6 - SUMMARY
# ============================================================
hr('CELL 6 - summary')

log(f"""
  endpoints reachable : gost={gost_ok}  convert={conv_ok}  kegg_rest={kegg_ok}
  gene lists          : {len(lists)}  ({GATE}, {LAYER})
  symbols seen        : {len(all_symbols):,}
  symbols mapped      : {n_map:,} ({100 * n_map / max(len(all_symbols), 1):.1f}%)
  lists >=10 mapped   : {int(rate['ora_viable'].sum())} of {len(rate)}
  contrasts both dirs : {int(viable_both.sum())}

  NEXT: read mapping_rate_by_contrast.csv and unmapped_genes.csv before running
  step 2. If a whole symbol class is missing (MAMU-*, KEG06_*) that is expected
  and harmless for pathway analysis. If 'standard_symbol_unmatched' is large,
  stop and check the namespace, because the enrichment would be running on a
  silently truncated gene universe.

  Then: shiv_enrichment_run.py
""")
with open(os.path.join(LOG_DIR, 'enrichment_prep_report.txt'), 'w') as fh:
    fh.write('\n'.join(_report))
log('=' * 90)
log('ENRICHMENT PREP COMPLETE - no enrichment run.')
log('=' * 90)
