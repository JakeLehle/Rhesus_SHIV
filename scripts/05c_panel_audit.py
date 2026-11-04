#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
05c — Panel audit: is every marker gene actually in the matrix?
================================================================
READ-ONLY. Loads only `var` from the annotated object, so it runs in seconds
and uses almost no memory. Writes nothing into any path a pipeline stage owns.

WHY THIS EXISTS
---------------
A panel gene absent from Mmul10 under its human name silently shortens the
panel, and a panel reduced to a non-specific gene produces CONFIDENT WRONG
answers, which is worse than Unresolved.

The live example that prompted this. The Monocyte Intermediate panel is
CD74 + FCGR3 + HLA-DRA and the 05 log reported:

    Intermediate   2/3   MISSING ['HLA-DRA']

HLA-DRA does not exist in Mmul10. The rhesus MHC class II alpha gene is
MAMU-DRA, and it appears throughout the top-marker tables alongside MAMU-DRB1,
MAMU-DRB5, MAMU-DPA, MAMU-DQA1 and MAMU-DQB1. So the panel ran on CD74 plus
FCGR3, and CD74 is one of the most abundant transcripts in B cells. That is
why "Intermediate" scored 1.716 in a Monocyte subcluster whose top markers are
CD79B, MS4A1, BANK1 and CD79A: the panel had become an MHC-II score wearing a
monocyte-subset label.

Eight more genes were reported missing across the panels and none of them has
been chased down:

    Tier 1   HBA2 (Erythrocytes), PF4 (Platelets)
    Tier 2   IL17A (T Th17)
             IGHD (B Naive), CD27 (B Memory), CD24 (B Transitional)
             LILRB2 (Monocyte NonClassical), HLA-DRA (Monocyte Intermediate)

The three B-cell ones matter most: IGHD, CD27 and CD24 are the discriminating
markers for naive, memory and transitional B cells, so losing all three is why
the B lineage cannot separate at all.

WHAT THIS DOES, AND WHAT IT DELIBERATELY DOES NOT DO
-----------------------------------------------------
For every missing gene it searches var_names and prints candidates. It does
NOT pick one. Gene-name matching is a judgement call that depends on knowing
the biology, and a wrong substitution is exactly the failure this script
exists to catch. The output is a shortlist for you to confirm, after which the
confirmed pairs go into GENE_ALIASES in 05 Cell 1.

Search strategies, in the order reported:
  1. case-insensitive exact match
  2. the MAMU- substitution for any HLA- gene
  3. the gene name as a whole-word-ish substring of a longer symbol
  4. a curated alias list (gene symbols that were renamed, or whose rhesus
     ortholog carries a different approved symbol)
  5. LOC identifiers are reported separately, since Mmul10 uses them heavily
     for genes with no approved symbol and they cannot be matched by name

Inputs   pipeline_results/05_embed_annotate/adata_annotated.h5ad
         (falls back to 04_ingest_qc/adata_qc.h5ad)
Outputs  pipeline_results/05c_panel_audit/
         panel_gene_status.csv     every panel gene, present or not
         missing_gene_candidates.csv   the shortlist to confirm

Environment:  conda activate sc_pre
Author: Jake Lehle
Date:   November 2026
"""

# %% ============================================================
# CELL 1 — CONFIGURATION
# ============================================================
import os
import re

import pandas as pd
import anndata as ad

import warnings
warnings.filterwarnings('ignore')

WORKDIR = '/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV'
CANDIDATES = [
    os.path.join(WORKDIR, 'pipeline_results', '05_embed_annotate',
                 'adata_annotated.h5ad'),
    os.path.join(WORKDIR, 'pipeline_results', '04_ingest_qc',
                 'adata_qc.h5ad'),
]
OUT_DIR = os.path.join(WORKDIR, 'pipeline_results', '05c_panel_audit')

# Copied verbatim from 05 Cell 1 so the audit tests the panels that ran.
TIER1_MARKERS = {
    'T cells':      ['CD3D', 'CD3E', 'CD3G'],
    'B cells':      ['CD79A', 'MS4A1', 'CD19', 'BANK1', 'PAX5'],
    'NK cells':     ['KLRD1', 'KLRC3', 'GNLY', 'PRF1', 'GZMB', 'NKG7'],
    'Monocytes':    ['S100A8', 'S100A9', 'VCAN', 'CD14', 'LYZ'],
    'Macrophages':  ['CD68', 'CD163', 'MRC1', 'MERTK', 'APOC1'],
    'DCs':          ['FLT3', 'CLEC9A', 'CLEC10A', 'IL3RA', 'CLEC4C'],
    'Mast cells':   ['CPA3', 'KIT', 'HPGD', 'ENPP3'],
    'Erythrocytes': ['HBB', 'HBA1', 'HBA2'],
    'Platelets':    ['PPBP', 'PF4', 'GP9'],
}
PANELS = {
    'T cells': {
        'Naive':         ['CCR7', 'SELL', 'TCF7', 'LEF1', 'IL7R'],
        'CentralMem':    ['CCR7', 'SELL', 'IL7R', 'CD28', 'GPR183'],
        'EffectorMem':   ['GZMK', 'CCL5', 'CST7', 'KLRG1'],
        'Cytotoxic':     ['GZMB', 'GNLY', 'PRF1', 'NKG7', 'GZMH', 'FGFBP2'],
        'Th1':           ['TBX21', 'IFNG', 'CXCR3', 'STAT1'],
        'Th17':          ['RORC', 'IL17A', 'CCR6', 'IL23R', 'KLRB1'],
        'Tfh':           ['CXCR5', 'PDCD1', 'BCL6', 'ICOS', 'IL21'],
        'Treg':          ['FOXP3', 'IL2RA', 'CTLA4', 'IKZF2'],
        'Proliferating': ['MKI67', 'TOP2A'],
    },
    'B cells': {
        'Naive':          ['TCL1A', 'IGHD', 'FCER2', 'IL4R'],
        'Memory':         ['CD27', 'TNFRSF13B', 'AIM2', 'TXNIP'],
        'Transitional':   ['MME', 'VPREB3', 'SOX4', 'CD24'],
        'Plasma':         ['MZB1', 'XBP1', 'PRDM1', 'JCHAIN', 'SDC1'],
        'GerminalCenter': ['BCL6', 'AICDA', 'RGS13', 'MEF2B'],
        'Activated':      ['CD69', 'CD83', 'NR4A1'],
    },
    'Monocytes': {
        'Classical':    ['CD14', 'S100A8', 'S100A9', 'VCAN'],
        'NonClassical': ['FCGR3', 'CDKN1C', 'LILRB2', 'MS4A7'],
        'Intermediate': ['CD74', 'FCGR3', 'HLA-DRA'],
        'Macrophage':   ['CD68', 'CD163', 'MRC1', 'MERTK', 'APOC1', 'C1QA'],
        'DC':           ['FLT3', 'CLEC9A', 'CLEC10A', 'IL3RA'],
    },
}
OTHER_GENES = {
    'CD3 gate':      ['CD3D', 'CD3E', 'CD3G'],
    'NK exclusive':  ['KLRD1', 'KLRC3', 'NCAM1'],
    'NK shared':     ['GNLY', 'PRF1', 'GZMB', 'NKG7'],
    'CD4/CD8':       ['CD4', 'CD8A', 'CD8B'],
    'Proliferation': ['MKI67', 'TOP2A'],
    'Activation':    ['FOS', 'FOSB', 'JUN', 'JUNB', 'JUND',
                      'EGR1', 'NR4A1', 'DUSP1'],
    'Viral entry':   ['CD4', 'CCR5', 'CXCR4'],
}

# Curated candidate aliases. These are SUGGESTIONS to check against the search
# output, never applied automatically. A gene with no entry here is still
# searched by pattern.
ALIAS_HINTS = {
    'HLA-DRA':  ['MAMU-DRA'],
    'HLA-DRB1': ['MAMU-DRB1', 'MAMU-DRB5'],
    'HLA-DQA1': ['MAMU-DQA1'],
    'HLA-DPA1': ['MAMU-DPA'],
    'FCGR3':    ['FCGR3A', 'FCGR3B', 'CD16'],
    'HBA2':     ['HBA1', 'HBA'],
    'PF4':      ['CXCL4', 'PF4V1'],
    'CD27':     ['TNFRSF7'],
    'CD24':     ['CD24A'],
    'IGHD':     ['IGHD1', 'IGD'],
    'IL17A':    ['IL17', 'CTLA8'],
    'LILRB2':   ['LIR2', 'ILT4', 'MAMU-LILRB2'],
    'NCAM1':    ['CD56'],
    'CLEC4C':   ['BDCA2', 'CD303'],
    'IL3RA':    ['CD123'],
    'ENPP3':    ['CD203C'],
    'GP9':      ['CD42B', 'GP1BB'],
    'TNFRSF13B': ['TACI'],
    'FCER2':    ['CD23'],
    'MME':      ['CD10'],
    'SDC1':     ['CD138'],
}

MAX_CANDIDATES = 8
os.makedirs(OUT_DIR, exist_ok=True)
pd.set_option('display.width', 220)
pd.set_option('display.max_rows', 400)
pd.set_option('display.max_colwidth', 70)


def hr(t):
    print('\n' + '=' * 94, flush=True)
    print(t, flush=True)
    print('=' * 94, flush=True)


def save_csv(df, name, index=False):
    df.to_csv(os.path.join(OUT_DIR, name), index=index, lineterminator='\n')
    print(f'  wrote {name}  ({df.shape[0]:,} rows)', flush=True)


# %% ============================================================
# CELL 2 — LOAD var_names ONLY
# ============================================================
hr('CELL 2 — gene universe')

src = next((c for c in CANDIDATES if os.path.exists(c)), None)
if src is None:
    raise SystemExit('No object found. Looked for:\n  ' + '\n  '.join(CANDIDATES))
print(f'  reading var from {src}')

# backed='r' keeps the matrix on disk; we only need var_names.
a = ad.read_h5ad(src, backed='r')
VAR = list(a.var_names)
a.file.close()
VAR_SET = set(VAR)
UPPER = {v.upper(): v for v in VAR}
print(f'  {len(VAR):,} features in the matrix')
n_loc = sum(1 for v in VAR if v.startswith('LOC'))
n_mamu = sum(1 for v in VAR if v.upper().startswith('MAMU'))
print(f'  LOC identifiers : {n_loc:,} ({100 * n_loc / len(VAR):.1f}%)')
print(f'  MAMU- symbols   : {n_mamu}')
print('  A high LOC fraction is why a missing gene often is not missing, it')
print('  just has no approved symbol in this annotation.')
if n_mamu:
    print(f'\n  MAMU- symbols present: '
          f'{sorted(v for v in VAR if v.upper().startswith("MAMU"))}')


# %% ============================================================
# CELL 3 — EVERY PANEL GENE, PRESENT OR NOT
# ============================================================
hr('CELL 3 — panel gene status')

rows = []
for tier, group in (('Tier1', {'_': TIER1_MARKERS}),
                    ('Tier2', PANELS),
                    ('Other', {'_': OTHER_GENES})):
    for lineage, panels in group.items():
        for panel, genes in panels.items():
            for g in genes:
                hit = g if g in VAR_SET else UPPER.get(g.upper())
                rows.append({
                    'tier': tier,
                    'lineage': lineage if lineage != '_' else panel,
                    'panel': panel,
                    'gene': g,
                    'present': hit is not None,
                    'matched_as': hit or '',
                    'case_mismatch': bool(hit and hit != g),
                })
status = pd.DataFrame(rows)
save_csv(status, 'panel_gene_status.csv')

summ = (status.groupby(['tier', 'lineage', 'panel'])
        .agg(n_genes=('gene', 'size'), n_present=('present', 'sum'))
        .reset_index())
summ['n_missing'] = summ.n_genes - summ.n_present
summ['usable'] = summ.n_present / summ.n_genes
print(summ.sort_values(['tier', 'usable']).to_string(index=False))

cm = status[status.case_mismatch]
if len(cm):
    print('\n  CASE MISMATCHES, which score_genes treats as absent:')
    print(cm[['panel', 'gene', 'matched_as']].to_string(index=False))

thin = summ[summ.n_present <= 2]
if len(thin):
    print('\n  Panels down to 2 or fewer genes. Each of these can be carried')
    print('  by a single gene, and if that gene is not lineage-specific the')
    print('  subtype call will be confidently wrong:')
    print(thin[['tier', 'lineage', 'panel', 'n_present',
                'n_genes']].to_string(index=False))


# %% ============================================================
# CELL 4 — CANDIDATES FOR EVERY MISSING GENE
# ============================================================
hr('CELL 4 — candidates for the missing genes')

missing = sorted(set(status.loc[~status.present, 'gene']))
print(f'  {len(missing)} distinct gene(s) missing: {missing}\n')

if not missing:
    print('  Nothing missing. No substitutions needed.')
    cand = pd.DataFrame()
else:
    def search(g):
        """Candidate symbols for a missing gene, best strategy first."""
        out = []

        def add(sym, how):
            if sym and sym in VAR_SET and sym not in [o[0] for o in out]:
                out.append((sym, how))

        for h in ALIAS_HINTS.get(g, []):
            add(h, 'curated hint')
            add(UPPER.get(h.upper()), 'curated hint (case)')

        if g.upper().startswith('HLA-'):
            add(UPPER.get('MAMU-' + g.upper()[4:]), 'HLA-> MAMU-')

        # the symbol as a prefix or suffix of a longer one, e.g. FCGR3 -> FCGR3A
        pat = re.compile(r'(^|[^A-Z0-9])' + re.escape(g.upper())
                         + r'($|[^A-Z0-9])')
        for v in VAR:
            if len(out) >= MAX_CANDIDATES:
                break
            vu = v.upper()
            if vu.startswith(g.upper()) or vu.endswith(g.upper()):
                add(v, 'prefix/suffix')
            elif pat.search(vu):
                add(v, 'word match')

        # trailing-digit variants: HBA2 -> HBA1, HBA
        base = g.rstrip('0123456789')
        if base != g:
            for alt in (base, base + '1', base + '2', base + '3'):
                if alt != g:
                    add(UPPER.get(alt.upper()), f'digit variant of {base}')
        return out[:MAX_CANDIDATES]

    rows = []
    for g in missing:
        used_in = sorted(set(
            status.loc[status.gene == g, 'panel'] + ' ('
            + status.loc[status.gene == g, 'lineage'] + ')'))
        hits = search(g)
        print(f'  {g}')
        print(f'    used by: {", ".join(used_in)}')
        if hits:
            for sym, how in hits:
                print(f'    candidate: {sym:<22} [{how}]')
        else:
            print('    no candidate found by name. If this gene matters it is')
            print('    probably a LOC identifier, which cannot be matched by')
            print('    symbol. Find it by locus or by correlation with a')
            print('    known partner gene instead.')
        print()
        for sym, how in (hits or [('', 'none found')]):
            rows.append({'missing_gene': g, 'used_by': '; '.join(used_in),
                         'candidate': sym, 'strategy': how, 'confirmed': ''})
    cand = pd.DataFrame(rows)
    save_csv(cand, 'missing_gene_candidates.csv')


# %% ============================================================
# CELL 5 — WHAT TO DO WITH THIS
# ============================================================
hr('DONE — 05c_panel_audit.py')
print(f'  outputs: {OUT_DIR}')
print('  nothing was modified.')
print()
print('  Fill the `confirmed` column of missing_gene_candidates.csv for the')
print('  substitutions you accept, then add them to GENE_ALIASES in 05 Cell 1')
print('  and re-run 05 with RESUME_FROM_TIER1 = False.')
print()
print('  Nothing here is applied automatically. A wrong substitution is')
print('  exactly the failure this script exists to catch, and deciding that')
print('  MAMU-DRA is the right stand-in for HLA-DRA needs someone who knows')
print('  the biology, not a regex.')
print()
print('  Priority order, by what it unblocks:')
print('    HLA-DRA  the Monocyte Intermediate panel is currently an MHC-II')
print('             score, and it mislabelled a 582-cell B-cell subcluster')
print('    CD27     }  the three discriminating markers for naive, memory')
print('    CD24     }  and transitional B cells. Losing all three is why the')
print('    IGHD     }  B lineage will not separate at all')
print('    LILRB2   NonClassical monocytes are down to 3 of 4')
print('    IL17A    Th17 is down to 4 of 5; least consequential, since Th17')
print('             barely varies across subclusters anyway')
print('    PF4      Platelets run on PPBP + GP9')
print('    HBA2     Erythrocytes run on HBB + HBA1')
