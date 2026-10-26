#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Rhesus SHIV - Annotation Reference Collector (marker panels + cluster breakout)
===============================================================================
Zoey's asks 1 and 2 (2026-10-09): one complete reference table of the tiered
marker genes, and one table saying what each cluster was called.

Neither artifact exists in the repo today. The marker panels live as hardcoded
dicts in four different scripts and they are NOT all the same set, and the
cluster-to-label mapping lives only inside the two annotation_write scripts and
is never written out. This script collects both into CSVs, with provenance, so
what goes to Zoey is traceable back to the line it came from.

HOW THE PANELS ARE HARVESTED, and why it matters. This script does NOT keep its
own copy of any panel. It parses each source script with `ast` and reads the
assignment node, so the source scripts stay the single source of truth and the
table cannot silently drift from the code that actually produced the labels. No
source script is imported or executed; nothing in the pipeline runs.

  Two entries are not plain literals and are handled explicitly:
    - MYELOID_PANELS['cDC2'] is a conditional on CDC2_KEEP_CLASS_II. The flag is
      parsed from the same file and the expression evaluated against it, so the
      table reflects the dial as it currently sits. The note column says so.
    - QC_PANELS['QC_mitochondrial'] and ['QC_ribosomal'] are None in source and
      filled at runtime from the KEG06_ and RPS/RPL prefixes. They appear in the
      table as runtime-derived with no gene list, which is accurate.

GENE PRESENCE is only checked when an .h5ad is available. The resolver mirrors
SHIV_tier2_myeloid_diagnostic.py CELL 2: exact match, then RHESUS_ALIASES
substitution, then suggestions from the MHC / paralog / LOC-column heuristics.
Suggestions are reported, never silently applied, same as that script. Without
an object the table still builds and the presence columns read 'not_checked'.

THE TIER 3 QUESTION. There is no Tier 3 cell type annotation in this repo. The
only 'Tier 3' in the codebase is the exact 34-mer stratum in the VPR barcode
search, which is unrelated. What exists below Tier 2 is a set of per-cell
attributes deliberately kept OUT of the labels (CD4/CD8 surface call,
proliferation, activation, the two myeloid state scores). CELL 6 writes those to
their own table so Zoey can see what is actually available and decide whether
that is what she meant.

Outputs, all to OUT_DIR:
  tiered_marker_genes_master.csv   every panel, every gene, with provenance
  tiered_marker_panel_summary.csv  one row per panel, for a quick read
  cluster_label_breakout.csv       cluster -> Tier 2 label, both compartments
  per_cell_attributes.csv          the sub-Tier-2 scores, i.e. the not-Tier-3
  annotation_reference_README.txt  what each column means, for sending on

Read-only throughout. Spyder cells (# %%). Headless-safe (sc_pre env).

Author: Jake Lehle
Date: October 2026
"""

# %% ============================================================
# CELL 1 - CONFIG AND DIALS
# ============================================================
import os
import ast
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings("ignore")

WORKING_DIR = '/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV'
MERGED      = os.path.join(WORKING_DIR, 'shiv_starsolo', 'merged')

# Where the pipeline scripts live. Set SCRIPT_DIR explicitly to override; left
# as None it is auto-located from the candidates below, because the GitHub repo
# keeps these in a scripts/ subdirectory while a working copy on titan often has
# them flat in the project directory, and guessing wrong silently produced an
# empty table in the 2026-10-09 run.
SCRIPT_DIR = None
SCRIPT_DIR_CANDIDATES = [
    os.path.dirname(os.path.abspath(__file__)),   # next to this script
    os.getcwd(),                                  # where you launched it
    os.path.join(WORKING_DIR, 'scripts'),
    WORKING_DIR,
    os.path.join(os.path.dirname(os.path.abspath(__file__)), 'scripts'),
    os.path.join(os.getcwd(), 'scripts'),
]
# a directory only counts if it actually holds the pipeline. Requiring several
# stops a directory with one stray copy from being mistaken for the repo.
SCRIPT_DIR_SENTINELS = [
    'scanpy_analysis_script.py',
    'shiv_tier2_diagnostic_v2.py',
    'SHIV_tier2_myeloid_diagnostic.py',
    'shiv_tier2_annotation_write.py',
    'SHIV_tier2_myeloid_annotation_write.py',
    'shiv_comprehensive_diagnostic.py',
]
MIN_SENTINELS = 4       # below this the directory is not the pipeline
FAIL_ON_EMPTY_HARVEST = True   # refuse to write an empty reference table

OUT_DIR     = os.path.join(WORKING_DIR, 'annotation_output', 'annotation_reference')

# Object used ONLY to check gene presence and to count cells per cluster.
# Optional: without it the tables still build, presence reads 'not_checked'.
GLOBAL_CANDIDATES = [
    os.path.join(MERGED, 'shiv_host_tier2_myeloid_annotated.h5ad'),
    os.path.join(MERGED, 'shiv_host_tier2_annotated.h5ad'),
    os.path.join(MERGED, 'shiv_host_final_tcrTrecovered.h5ad'),
]
TCELL_SUBSET   = os.path.join(MERGED, 'shiv_Tcell_annotated.h5ad')
MYELOID_SUBSET = os.path.join(MERGED, 'shiv_Myeloid_annotated.h5ad')

CHECK_PRESENCE   = True    # False skips loading any object entirely
VERIFY_AGAINST_OBJECTS = True   # recount clusters from the objects, not the dicts

# ---- where each panel comes from. (file, dict name, tier, how to read it) ----
# 'scope' is the lineage the panel subdivides. 'stage' says which pass of the
# pipeline the panel belongs to, because the same lineage has different panels
# at the global pass and at the Tier 2 reclustering pass, and that difference is
# the single most confusing thing about this table.
PANEL_SOURCES = [
    # file,                              dict,              tier,  stage,              scope
    ('scanpy_analysis_script.py',        'TIER1_MARKERS',   'Tier 1', 'global pass',      None),
    ('scanpy_analysis_script.py',        'TIER2_MARKERS',   'Tier 2', 'global pass',      None),
    ('shiv_tier2_diagnostic_v2.py',      'PANELS',          'Tier 2', 'reclustering',     None),
    ('SHIV_tier2_myeloid_diagnostic.py', 'MYELOID_PANELS',  'Tier 2', 'reclustering',     'Myeloid'),
    ('SHIV_tier2_myeloid_diagnostic.py', 'EXTRA_PANELS',    'Tier 2', 'reclustering',     'Myeloid'),
    ('SHIV_tier2_myeloid_diagnostic.py', 'QC_PANELS',       'QC',     'reclustering',     'Myeloid'),
    ('SHIV_tier2_myeloid_diagnostic.py', 'CONTAM_PANELS',   'QC',     'reclustering',     'Myeloid'),
    ('shiv_comprehensive_diagnostic.py', 'RIRA_PANEL',      'overlay', 'RIRA cross-check', None),
]

# flat gene-list constants worth carrying, same (file, name, tier, stage, scope,
# panel_name) shape flattened into the table as single-panel entries
FLAT_SOURCES = [
    ('shiv_tier2_diagnostic_v2.py',      'CD3_GENES',  'Tier 1', 'reclustering', 'T cells', 'CD3 gate'),
    ('shiv_tier2_diagnostic_v2.py',      'CD4_GENES',  'Tier 2', 'reclustering', 'T cells', 'CD4 call'),
    ('shiv_tier2_diagnostic_v2.py',      'CD8_GENES',  'Tier 2', 'reclustering', 'T cells', 'CD8 call'),
    ('shiv_tier2_annotation_write.py',   'PROLIF_GENES',     'attribute', 'annotation', 'T cells', 'Proliferation gate'),
    ('shiv_tier2_annotation_write.py',   'ACTIVATION_GENES', 'attribute', 'annotation', 'T cells', 'CD4 activation (soft)'),
    ('SHIV_tier2_myeloid_diagnostic.py', 'MYELOID_ANCHORS',  'QC',        'reclustering', 'Myeloid', 'Myeloid purity anchors'),
]

# ---- cluster -> label maps ----
LABEL_SOURCES = [
    # file, dict, compartment, qc dict, sizes dict, subtype col, subset path
    ('shiv_tier2_annotation_write.py', 'SUBTYPE_MAP', 'T cells',
     None, None, 'tier2_subtype', TCELL_SUBSET),
    ('SHIV_tier2_myeloid_annotation_write.py', 'LABEL_MAP', 'Myeloid',
     'QC_FLAG', 'EXPECTED_SIZES', 'tier2_myeloid_subtype', MYELOID_SUBSET),
]

# ---- per-cell attributes: the things that are NOT a Tier 3 ----
ATTR_SPEC = [
    ('T cells', 'cd48_adt',            'CD4/CD8 surface call (CD8 CLR > CD4 CLR)',
     'kept as an attribute so the naive pool stays one label; split with Zoey later'),
    ('T cells', 'is_proliferating',    'proliferation gate (MKI67/TOP2A)',
     'the clean proliferating cluster dissolved on the recovered pool'),
    ('T cells', 'sig_Proliferating',   'continuous proliferation score', ''),
    ('T cells', 'sig_Activation',      'CD4 activation, immediate-early genes',
     'partly a dissociation artifact; soft score, never a label'),
    ('Myeloid', 'state_Antigen_Presenting_Myeloid', 'MHC-II-high state score',
     'a state any myeloid lineage can enter, not a lineage'),
    ('Myeloid', 'state_IFN_Stimulated_Myeloid',     'interferon-stimulated state score',
     'a state any myeloid lineage can enter, not a lineage'),
    ('Myeloid', 'qc_exclude_myeloid',  'QC exclusion flag', ''),
    ('Myeloid', 'myeloid_analysis_ready', 'myeloid AND not QC AND subtype resolved', ''),
]

os.makedirs(OUT_DIR, exist_ok=True)
pd.set_option('display.width', 240)
pd.set_option('display.max_columns', 60)
pd.set_option('display.max_rows', 500)


def hr(t):
    print('\n' + '=' * 90)
    print(t)
    print('=' * 90)


def save_csv(df, name, index=False):
    df.to_csv(os.path.join(OUT_DIR, name), index=index)
    print(f'    saved: {name}  ({len(df)} rows)')


# %% ============================================================
# CELL 2 - AST HARVEST (parse, never import or execute)
# ============================================================
hr('CELL 2 - harvest dicts from the source scripts by parsing, not importing')


def locate_script_dir():
    """Find the directory holding the pipeline scripts, or say clearly why not.

    A directory qualifies only if it contains at least MIN_SENTINELS of the
    known pipeline files. Scoring rather than taking the first hit means a
    folder with one stray copy of a script does not win over the real one.
    """
    if SCRIPT_DIR:
        found = [f for f in SCRIPT_DIR_SENTINELS
                 if os.path.exists(os.path.join(SCRIPT_DIR, f))]
        print(f'  SCRIPT_DIR set explicitly: {SCRIPT_DIR}')
        print(f'    pipeline scripts found there: {len(found)} of '
              f'{len(SCRIPT_DIR_SENTINELS)}')
        if len(found) < MIN_SENTINELS:
            raise FileNotFoundError(
                f'SCRIPT_DIR is set to {SCRIPT_DIR} but only {len(found)} of '
                f'{len(SCRIPT_DIR_SENTINELS)} pipeline scripts are there. '
                'Point it at your Rhesus_SHIV checkout.')
        return SCRIPT_DIR

    scored = []
    seen = set()
    for d in SCRIPT_DIR_CANDIDATES:
        d = os.path.abspath(d)
        if d in seen or not os.path.isdir(d):
            continue
        seen.add(d)
        n = sum(os.path.exists(os.path.join(d, f)) for f in SCRIPT_DIR_SENTINELS)
        scored.append((n, d))
    scored.sort(key=lambda t: -t[0])

    print('  locating the pipeline scripts:')
    for n, d in scored:
        mark = '  <-- using' if (n >= MIN_SENTINELS and d == scored[0][1]) else ''
        print(f'    {n}/{len(SCRIPT_DIR_SENTINELS)}  {d}{mark}')

    if not scored or scored[0][0] < MIN_SENTINELS:
        raise FileNotFoundError(
            'Could not find the Rhesus_SHIV pipeline scripts. None of the '
            f'candidate directories holds at least {MIN_SENTINELS} of '
            f'{SCRIPT_DIR_SENTINELS}.\n'
            'Set SCRIPT_DIR at the top of this script to your checkout, e.g.\n'
            "    SCRIPT_DIR = '/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV'\n"
            'Nothing was written. This script reads the panels OUT of those '
            'files, so without them there is no table to build.')
    return scored[0][1]


SCRIPT_DIR = locate_script_dir()

_tree_cache = {}


def parse_file(fname):
    """Parse a source script to an AST once and cache it."""
    if fname in _tree_cache:
        return _tree_cache[fname]
    path = os.path.join(SCRIPT_DIR, fname)
    if not os.path.exists(path):
        print(f'  WARN: source script not found, skipping: {fname}')
        _tree_cache[fname] = None
        return None
    with open(path, 'r', encoding='utf-8') as fh:
        tree = ast.parse(fh.read(), filename=fname)
    _tree_cache[fname] = tree
    return tree


def find_assign(tree, name):
    """Return the value node of the LAST top-level `name = ...` assignment."""
    if tree is None:
        return None
    found = None
    for node in tree.body:
        if isinstance(node, ast.Assign):
            for tgt in node.targets:
                if isinstance(tgt, ast.Name) and tgt.id == name:
                    found = node.value
    return found


def simple_flags(tree):
    """Top-level `NAME = <literal>` bindings, for resolving conditionals."""
    env = {}
    if tree is None:
        return env
    for node in tree.body:
        if isinstance(node, ast.Assign) and len(node.targets) == 1 \
           and isinstance(node.targets[0], ast.Name):
            try:
                env[node.targets[0].id] = ast.literal_eval(node.value)
            except Exception:
                pass
    return env


def eval_node(node, env, where):
    """literal_eval first; fall back to a restricted eval for conditionals.

    The fallback exists for exactly one case in this repo, MYELOID_PANELS['cDC2'],
    which is a list concatenation gated on CDC2_KEEP_CLASS_II. The namespace is
    the literal flags parsed from the same file and nothing else: no builtins, no
    imports. Every use of the fallback is printed and noted in the output table.
    """
    try:
        return ast.literal_eval(node), ''
    except Exception:
        pass
    try:
        code = compile(ast.Expression(body=node), '<panel>', 'eval')
        val = eval(code, {'__builtins__': {}}, dict(env))  # noqa: S307
        flags = sorted(set(env) & set(
            n.id for n in ast.walk(node) if isinstance(n, ast.Name)))
        note = f'conditional on {", ".join(flags)}' if flags else 'conditional'
        print(f'    {where}: not a literal -> evaluated against {flags or "{}"}')
        return val, note
    except Exception as e:
        print(f'    WARN {where}: could not evaluate ({e})')
        return None, 'UNPARSED'


def get_dict(fname, dname):
    """Harvest one dict, entry by entry, so one bad entry does not lose the rest."""
    tree = parse_file(fname)
    node = find_assign(tree, dname)
    if node is None:
        print(f'  WARN: {dname} not found in {fname}')
        return {}, {}
    env = simple_flags(tree)
    if not isinstance(node, ast.Dict):
        val, note = eval_node(node, env, f'{fname}:{dname}')
        return (val or {}), {}
    out, notes = {}, {}
    for k_node, v_node in zip(node.keys, node.values):
        try:
            key = ast.literal_eval(k_node)
        except Exception:
            continue
        val, note = eval_node(v_node, env, f'{fname}:{dname}[{key!r}]')
        out[key] = val
        if note:
            notes[key] = note
    return out, notes


def get_list(fname, name, optional=False):
    """Harvest a top-level list constant. optional=True for names that legitimately
    exist in only one of the two annotation scripts (QC_LABELS is T-side only,
    NON_MYELOID_CLUSTERS is myeloid-side only), so absence is not a warning."""
    tree = parse_file(fname)
    node = find_assign(tree, name)
    if node is None:
        if not optional:
            print(f'  WARN: {name} not found in {fname}')
        return None
    val, _ = eval_node(node, simple_flags(tree), f'{fname}:{name}')
    return val


harvested = {}
for fname, dname, tier, stage, scope in PANEL_SOURCES:
    d, notes = get_dict(fname, dname)
    harvested[(fname, dname)] = (d, notes)
    print(f'  {fname:36s} {dname:18s} -> {len(d)} panel(s)')

# the collinearity merge structure and the alias table, carried as context
TIER2_GROUPS, _ = get_dict('SHIV_tier2_myeloid_diagnostic.py', 'TIER2_GROUPS')
RHESUS_ALIASES, _ = get_dict('SHIV_tier2_myeloid_diagnostic.py', 'RHESUS_ALIASES')
STATE_PANELS = get_list('SHIV_tier2_myeloid_diagnostic.py', 'STATE_PANELS') or []
print(f'\n  context: {len(TIER2_GROUPS)} Tier-2 myeloid groups, '
      f'{len(RHESUS_ALIASES)} rhesus aliases, {len(STATE_PANELS)} state panels')

# ---- precondition: did anything actually come back? ----
# The 2026-10-09 run harvested zero panels from zero files and still wrote the
# CSVs and printed a success banner. A 0-row reference table looks like a
# deliverable and is the kind of thing that gets forwarded to a collaborator.
# Warn-and-continue is right for ONE missing panel; it is wrong when the whole
# input set is missing, because then there is no output, only an empty shell.
_total_panels = sum(len(d) for d, _ in harvested.values())
if _total_panels == 0:
    msg = (f'ZERO panels harvested from {len(PANEL_SOURCES)} source entries in '
           f'{SCRIPT_DIR}.\nEvery dict lookup failed, which means the files are '
           'not there or are not the expected scripts.\nNothing has been '
           'written. Fix SCRIPT_DIR and rerun.')
    print('\n  *** ' + msg.replace('\n', '\n      '))
    if FAIL_ON_EMPTY_HARVEST:
        raise RuntimeError(msg)
print(f'  harvested {_total_panels} panel(s) in total')

# which fine panel rolls up into which group, for the group column
PANEL_TO_GROUP = {}
for grp, members in (TIER2_GROUPS or {}).items():
    for m in (members or []):
        PANEL_TO_GROUP[m] = grp


# %% ============================================================
# CELL 3 - OPTIONAL OBJECT LOAD + SYMBOL RESOLVER
# ============================================================
hr('CELL 3 - gene presence resolver')

adata = None
var_names = []
var_set, var_upper = set(), {}

if CHECK_PRESENCE:
    try:
        import anndata as ad
        path = next((p for p in GLOBAL_CANDIDATES if os.path.exists(p)), None)
        if path is None:
            print('  no global object found; presence will read not_checked')
            print('  looked for:\n    ' + '\n    '.join(GLOBAL_CANDIDATES))
        else:
            adata = ad.read_h5ad(path, backed='r')
            var_names = list(adata.var_names)
            var_set = set(var_names)
            var_upper = {v.upper(): v for v in var_names}
            print(f'  loaded (backed, read-only): {os.path.basename(path)}')
            print(f'    {adata.n_obs:,} cells x {adata.n_vars:,} genes')
    except Exception as e:
        print(f'  WARN: could not load an object for presence checking ({e})')
else:
    print('  CHECK_PRESENCE is False; skipping object load')


def suggest_alternatives(g, limit=4):
    """Candidates for a symbol absent from Mmul10.

    Trimmed from SHIV_tier2_myeloid_diagnostic.py CELL 2. Suggestion-only by
    design: rhesus renames MHC to MAMU- and paralog suffixes drift, so this
    surfaces options for Zoey to confirm rather than silently picking one.
    """
    if not var_names:
        return []
    g_up = g.upper()
    cands = []
    if g_up.startswith('HLA-'):
        stem = g_up[4:]
        cands += [v for v in var_names if v.upper() == f'MAMU-{stem}']
        cands += [v for v in var_names
                  if v.upper().startswith('MAMU-') and stem[:2] in v.upper()]
    base = g_up.rstrip('ABCDEFG0123456789')
    if len(base) >= 3:
        cands += [v for v in var_names if v.upper().startswith(base)]
    if any(ch.isdigit() for ch in g_up):
        fam = g_up.rstrip('0123456789')
        if len(fam) >= 4:
            cands += [v for v in var_names if v.upper().startswith(fam)]
    if adata is not None:
        for col in adata.var.columns:
            try:
                ser = adata.var[col].astype(str)
            except Exception:
                continue
            hit = [v for v, s in zip(var_names, ser) if s.upper() == g_up]
            if hit:
                cands += [f'{h}[{col}]' for h in hit[:2]]
    seen, out = set(), []
    for c in cands:
        if c not in seen and c.upper() != g_up:
            seen.add(c)
            out.append(c)
    return out[:limit]


def resolve(gene):
    """(status, resolved_symbol, suggestions) for one requested symbol."""
    if not var_set:
        return 'not_checked', '', ''
    if gene in var_set:
        return 'present', gene, ''
    if gene.upper() in var_upper:
        return 'present_case', var_upper[gene.upper()], ''
    alias_hits = [a for a in (RHESUS_ALIASES or {}).get(gene, []) if a in var_set]
    if alias_hits:
        return 'aliased', '+'.join(alias_hits), ''
    return 'ABSENT', '', '; '.join(suggest_alternatives(gene))


# %% ============================================================
# CELL 4 - BUILD THE MASTER MARKER TABLE
# ============================================================
hr('CELL 4 - master marker gene table')

rows = []


def add_panel(tier, stage, scope, lineage, panel, genes, source_file,
              source_dict, note=''):
    group = PANEL_TO_GROUP.get(panel, '')
    is_state = 'state, not lineage' if panel in STATE_PANELS else ''
    if genes is None:
        rows.append({
            'tier': tier, 'stage': stage, 'lineage_scope': scope or '',
            'lineage': lineage or '', 'panel': panel,
            'tier2_group': group, 'gene_as_written': '',
            'status': 'runtime_derived', 'resolved_symbol': '',
            'suggestions': '', 'source_script': source_file,
            'source_dict': source_dict,
            'note': (note + '; ' if note else '') +
                    'gene list built at runtime from a prefix, not literal in source',
        })
        return
    seen = set()
    for g in genes:
        dup = g in seen
        seen.add(g)
        status, resolved, sugg = resolve(g)
        n = [x for x in (note, is_state) if x]
        if dup:
            n.append('listed twice in the source panel')
        rows.append({
            'tier': tier, 'stage': stage, 'lineage_scope': scope or '',
            'lineage': lineage or '', 'panel': panel,
            'tier2_group': group, 'gene_as_written': g,
            'status': status, 'resolved_symbol': resolved,
            'suggestions': sugg, 'source_script': source_file,
            'source_dict': source_dict, 'note': '; '.join(n),
        })


for fname, dname, tier, stage, scope in PANEL_SOURCES:
    d, notes = harvested.get((fname, dname), ({}, {}))
    for key, val in (d or {}).items():
        note = notes.get(key, '')
        # TIER1_MARKERS is {lineage: {'markers': [...]}}; everything else is
        # either {panel: [genes]} or {lineage: {panel: [genes]}}
        if isinstance(val, dict) and 'markers' in val:
            add_panel(tier, stage, scope, key, key, val['markers'],
                      fname, dname, note)
        elif isinstance(val, dict):
            for panel, genes in val.items():
                add_panel(tier, stage, scope, key, panel, genes,
                          fname, dname, note)
        else:
            add_panel(tier, stage, scope, scope, key, val, fname, dname, note)

for fname, name, tier, stage, scope, panel in FLAT_SOURCES:
    genes = get_list(fname, name)
    if genes is not None:
        add_panel(tier, stage, scope, scope, panel, genes, fname, name)

master = pd.DataFrame(rows)
if len(master):
    master = master.sort_values(
        ['tier', 'stage', 'lineage_scope', 'lineage', 'panel', 'gene_as_written'],
        kind='stable').reset_index(drop=True)
save_csv(master, 'tiered_marker_genes_master.csv')

print(f'\n  {len(master)} gene entries across '
      f'{master["panel"].nunique() if len(master) else 0} panels')
if len(master):
    print('\n  by tier and stage:')
    print(master.groupby(['tier', 'stage']).size().rename('genes').to_string())
    print('\n  resolution status:')
    print(master['status'].value_counts().to_string())

    # Sanity check of the same class as the empty-harvest one. A handful of
    # ABSENT symbols is the real finding (FCN1, SPP1, IFITM3 and friends). Most
    # of them absent is not a biology result, it means the resolver is looking
    # at the wrong thing: a var_names / gene_symbol mismatch, or an object whose
    # var is not the gene axis at all. Flagged rather than reported as 300
    # missing genes.
    _checked = master[master['status'].isin(
        ['present', 'present_case', 'aliased', 'ABSENT'])]
    if len(_checked):
        _frac = (_checked['status'] == 'ABSENT').mean()
        if _frac > 0.5:
            print(f'\n  *** {100 * _frac:.0f}% of symbols resolved as ABSENT '
                  f'({int((_checked["status"] == "ABSENT").sum())} of '
                  f'{len(_checked)}).')
            print('  That is too many to be real. Check that the object in '
                  'GLOBAL_CANDIDATES is the gene-level annotated object and '
                  'that var_names hold gene symbols.')
            print('  Treat the presence columns in the CSV as unreliable until '
                  'this is resolved.')

    absent = master[master['status'] == 'ABSENT']
    if len(absent):
        print(f'\n  *** {len(absent)} symbol(s) NOT in the reference. '
              'These scored on fewer genes than the panel claims:')
        for _, r in absent.iterrows():
            s = f'   suggestions: {r["suggestions"]}' if r['suggestions'] else ''
            print(f'    {r["panel"]:<30s} {r["gene_as_written"]:<10s}{s}')
        print('  Suggestions are NOT applied. Confirm with Zoey, then edit the')
        print('  panel in its source script, not here.')

# one row per panel, the version most people will actually read
if len(master):
    g = master[master['status'] != 'runtime_derived'].groupby(
        ['tier', 'stage', 'lineage_scope', 'lineage', 'panel', 'tier2_group'],
        dropna=False)
    summary = g.agg(
        n_genes=('gene_as_written', 'nunique'),
        n_present=('status', lambda s: int((s != 'ABSENT').sum())),
        n_absent=('status', lambda s: int((s == 'ABSENT').sum())),
        genes=('gene_as_written', lambda s: ', '.join(dict.fromkeys(s))),
        absent_genes=('gene_as_written', lambda s: ''),
        source_script=('source_script', 'first'),
        source_dict=('source_dict', 'first'),
    ).reset_index()
    ab = (master[master['status'] == 'ABSENT']
          .groupby('panel')['gene_as_written']
          .agg(lambda s: ', '.join(dict.fromkeys(s))))
    summary['absent_genes'] = summary['panel'].map(ab).fillna('')
    save_csv(summary, 'tiered_marker_panel_summary.csv')


# %% ============================================================
# CELL 5 - CLUSTER -> LABEL BREAKOUT
# ============================================================
hr('CELL 5 - cluster to Tier 2 label breakout')

lab_rows = []
for fname, dname, compartment, qc_name, size_name, subtype_col, subset_path \
        in LABEL_SOURCES:
    lmap, _ = get_dict(fname, dname)
    if not lmap:
        print(f'  WARN: no label map harvested from {fname}:{dname}')
        continue
    qc_map, _ = get_dict(fname, qc_name) if qc_name else ({}, {})
    sizes, _ = get_dict(fname, size_name) if size_name else ({}, {})
    nonmy = get_list(fname, 'NON_MYELOID_CLUSTERS', optional=True) or []
    qc_labels = get_list(fname, 'QC_LABELS', optional=True) or []

    # recount from the object where we can, so the table is verified not asserted
    obs_counts, obs_ready = {}, {}
    if VERIFY_AGAINST_OBJECTS and os.path.exists(subset_path):
        try:
            import anndata as ad
            sb = ad.read_h5ad(subset_path, backed='r')
            ck = 'sub_leiden'
            if ck in sb.obs.columns:
                obs_counts = (sb.obs[ck].astype(str).value_counts().to_dict())
            if 'myeloid_analysis_ready' in sb.obs.columns:
                obs_ready = (sb.obs.groupby(sb.obs[ck].astype(str))
                             ['myeloid_analysis_ready'].any().to_dict())
            print(f'  {compartment}: recounted from '
                  f'{os.path.basename(subset_path)} ({sb.n_obs:,} cells)')
        except Exception as e:
            print(f'  WARN: could not recount {compartment} from object ({e})')
    else:
        print(f'  {compartment}: object not available, using the map as written')

    for cl in sorted(lmap, key=lambda x: int(x) if str(x).isdigit() else 1e9):
        label = lmap[cl]
        n_map = sizes.get(cl, '')
        n_obs = obs_counts.get(str(cl), '')
        drift = ''
        if n_map != '' and n_obs != '' and int(n_map) != int(n_obs):
            drift = f'SIZE DRIFT: map expects {n_map}, object has {n_obs}'
        lab_rows.append({
            'compartment': compartment,
            'cluster': cl,
            'tier2_label': label,
            'n_cells_expected': n_map,
            'n_cells_observed': n_obs,
            'qc_flagged': bool(qc_map.get(cl)) if qc_map else (label in qc_labels),
            'qc_reason': qc_map.get(cl, '') if qc_map else (
                'QC label' if label in qc_labels else ''),
            'excluded_as_non_lineage': cl in nonmy or str(label).startswith('EXCLUDE'),
            'analysis_ready': obs_ready.get(str(cl), ''),
            'subtype_column': subtype_col,
            'source_script': fname,
            'source_dict': dname,
            'note': drift,
        })

breakout = pd.DataFrame(lab_rows)
save_csv(breakout, 'cluster_label_breakout.csv')
if len(breakout):
    print()
    print(breakout[['compartment', 'cluster', 'tier2_label', 'n_cells_expected',
                    'n_cells_observed', 'qc_flagged',
                    'excluded_as_non_lineage']].to_string(index=False))
    drifted = breakout[breakout['note'] != '']
    if len(drifted):
        print(f'\n  *** {len(drifted)} cluster(s) differ in size from the label map.')
        print('  The label maps key on cluster ID. If the clustering shifted, the')
        print('  labels are wrong. Resolve before this table goes to anyone.')
        for _, r in drifted.iterrows():
            print(f'    {r["compartment"]} cluster {r["cluster"]}: {r["note"]}')

    print('\n  collapsed to unique subtypes:')
    for comp, grp in breakout.groupby('compartment'):
        print(f'\n   {comp}:')
        for lab, sub in grp.groupby('tier2_label'):
            cls = ', '.join(str(c) for c in sub['cluster'])
            ns = [int(x) for x in sub['n_cells_observed'] if x != '']
            if not ns:
                ns = [int(x) for x in sub['n_cells_expected'] if x != '']
            tot = f'  n={sum(ns):,}' if ns else ''
            print(f'     {lab:<38s} clusters {cls}{tot}')


# %% ============================================================
# CELL 6 - PER-CELL ATTRIBUTES (what is NOT a Tier 3)
# ============================================================
hr('CELL 6 - per-cell attributes below Tier 2')

attr_rows = []
present_cols = {}
if VERIFY_AGAINST_OBJECTS:
    for comp, p in [('T cells', TCELL_SUBSET), ('Myeloid', MYELOID_SUBSET)]:
        if os.path.exists(p):
            try:
                import anndata as ad
                sb = ad.read_h5ad(p, backed='r')
                present_cols[comp] = set(sb.obs.columns)
            except Exception as e:
                print(f'  WARN: could not read {comp} subset ({e})')

for comp, col, desc, note in ATTR_SPEC:
    found = ''
    if comp in present_cols:
        found = 'yes' if col in present_cols[comp] else 'NOT FOUND'
    attr_rows.append({'compartment': comp, 'obs_column': col,
                      'description': desc, 'present_on_object': found,
                      'note': note})
attrs = pd.DataFrame(attr_rows)
save_csv(attrs, 'per_cell_attributes.csv')
print(attrs[['compartment', 'obs_column', 'present_on_object',
             'description']].to_string(index=False))
print('\n  These are continuous scores and flags carried per cell, deliberately')
print('  kept out of the subtype labels. They are NOT a third annotation tier.')
print('  No Tier 3 cell type annotation exists in this repo.')


# %% ============================================================
# CELL 7 - README FOR THE TABLES
# ============================================================
hr('CELL 7 - README')

readme = f"""Rhesus SHIV annotation reference
Generated {pd.Timestamp.today().date()} by SHIV_annotation_reference_collector.py
All content parsed directly from the pipeline scripts in {SCRIPT_DIR}.

FILES
-----
tiered_marker_genes_master.csv   one row per gene per panel
tiered_marker_panel_summary.csv  one row per panel
cluster_label_breakout.csv       cluster -> Tier 2 label
per_cell_attributes.csv          per-cell scores below Tier 2

READ THIS FIRST
---------------
The same lineage has DIFFERENT panels at different stages. 'stage' tells them
apart. 'global pass' is scanpy_analysis_script.py, which produced tier1_celltype
and tier2_celltype on the whole object. 'reclustering' is the per-lineage Tier 2
diagnostics, which produced the fine labels actually used for DEG. Where the two
disagree, the reclustering panel is the one behind the current labels.

COLUMNS, tiered_marker_genes_master.csv
---------------------------------------
tier             Tier 1 / Tier 2 / QC / overlay / attribute
stage            which pipeline pass the panel belongs to
lineage_scope    the Tier 1 lineage the panel subdivides, where applicable
lineage          the lineage the panel sits under in its source dict
panel            the subtype or signature name
tier2_group      the merged group this panel rolls into, for myeloid only.
                 Panels correlating above ~0.8 across subclusters were grouped
                 because no clustering can separate signatures measuring largely
                 the same genes. The call is made at group level.
gene_as_written  the symbol exactly as it appears in the source script
status           present / present_case / aliased / ABSENT / not_checked /
                 runtime_derived
resolved_symbol  what it resolved to in Mmul10
suggestions      candidate symbols for an ABSENT gene. NOT applied. Confirm,
                 then edit the panel in its source script.
note             state-vs-lineage flags, duplicates, conditional panels

CAVEATS
-------
- An ABSENT gene means the panel scored on fewer genes than it lists. Known
  open cases carried from the v2 myeloid run: FCN1, SPP1, IFITM3, all likely
  filed under LOC identifiers in Mmul10.
- Rhesus renames MHC to MAMU- and paralog suffixes drift (FCGR3A vs FCGR3).
  Substitutions are listed under 'aliased' and are reversible by setting
  USE_RHESUS_ALIASES = False in the myeloid diagnostic.
- Antigen Presenting Myeloid and IFN-Stimulated Myeloid describe a STATE, not a
  lineage, and are excluded from the lineage call for that reason.
- There is no Tier 3 cell type annotation. See per_cell_attributes.csv.
"""
with open(os.path.join(OUT_DIR, 'annotation_reference_README.txt'), 'w') as fh:
    fh.write(readme)
print('    saved: annotation_reference_README.txt')

print(f"""
  SUMMARY
  -------
  marker gene entries   : {len(master):,}
  panels                : {master['panel'].nunique() if len(master) else 0}
  symbols ABSENT        : {int((master['status'] == 'ABSENT').sum()) if len(master) else 0}
  clusters in breakout  : {len(breakout)}
  per-cell attributes   : {len(attrs)}
  output                : {OUT_DIR}
""")
print('=' * 90)
# Only claim success if there is something in the tables. An empty run that
# prints a success banner is how a 0-row CSV ends up being forwarded to someone.
_ok = len(master) > 0 and len(breakout) > 0
if _ok:
    print('ANNOTATION REFERENCE COLLECTED - nothing was imported, executed or '
          'modified.')
else:
    print('INCOMPLETE - the tables below are empty or partial, DO NOT SEND THEM:')
    if len(master) == 0:
        print('  tiered_marker_genes_master.csv   : 0 rows (no panels harvested)')
    if len(breakout) == 0:
        print('  cluster_label_breakout.csv       : 0 rows (no label maps found)')
    print(f'  Check SCRIPT_DIR. Currently: {SCRIPT_DIR}')
print('=' * 90)
