#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Rhesus SHIV - Tier 2 (v2) cluster adjudication
==============================================
Read-only. Reports per sub_leiden: doublet score (existing QC flag if present,
else a fresh Scrublet pass on the subset counts), CD3 / B-gene / MHC-II /
proliferation fractions. Purpose: decide whether the CD79/MHC-II-high subcluster
(cluster 9 in the v2 run) is residual B-doublet contamination (quarantine the
high-score cells) or a genuine MHC-II-high activated T subtype (label it).

Run WITH python (sc_pre env).
"""
import os
import numpy as np
import pandas as pd
import scanpy as sc
import warnings
warnings.filterwarnings("ignore")

SUBSET = ('/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV/shiv_starsolo/merged/'
          'shiv_Tcell_subset_reclustered_v2.h5ad')
OUT_DIR = ('/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV/annotation_output/'
           'tier2_Tcell_v2')
CLUSTER_KEY = 'sub_leiden'

sub = sc.read_h5ad(SUBSET)
print(f"loaded: {sub.n_obs:,} cells x {sub.n_vars:,} genes")
lay = 'counts' if 'counts' in sub.layers else None
print(f"counts layer for gene calls: {lay}")


def anypos(genes):
    """Boolean per cell: any of `genes` detected (>0) on the counts layer."""
    g = [x for x in genes if x in sub.var_names]
    if not g:
        return pd.Series(False, index=sub.obs_names)
    m = sub[:, g].layers[lay] if lay else sub[:, g].X
    m = m.toarray() if hasattr(m, 'toarray') else np.asarray(m)
    return pd.Series((m > 0).any(axis=1), index=sub.obs_names)


Bpos   = anypos(['CD79A', 'CD79B', 'MS4A1', 'BANK1'])
CD3pos = anypos(['CD3D', 'CD3E', 'CD3G'])
MHC2   = anypos(['CD74', 'MAMU-DRA', 'MAMU-DRB1'])
Prolif = anypos(['MKI67', 'STMN1', 'TOP2A'])

# --- doublet score: use existing QC flag, else run a fresh Scrublet pass ---
have_dblt = 'doublet_score' in sub.obs.columns
if not have_dblt:
    print("\n'doublet_score' NOT in subset obs (QC flag did not propagate).")
    print("Running a fresh Scrublet pass on the subset counts as a fallback...")
    try:
        import scrublet as scr
        counts = sub.layers['counts'] if lay else sub.X
        scrub = scr.Scrublet(counts, expected_doublet_rate=0.08)
        scores, preds = scrub.scrub_doublets(min_counts=2, min_cells=3,
                                             verbose=False)
        sub.obs['doublet_score'] = scores
        sub.obs['predicted_doublet'] = preds
        have_dblt = True
        print(f"  fresh Scrublet done; predicted doublets: {int(preds.sum()):,} "
              f"({100*preds.mean():.1f}%)")
    except Exception as e:
        print(f"  WARN fresh Scrublet failed: {e}")
else:
    print("\nUsing existing 'doublet_score' from QC.")

# --- per-cluster table ---
rows = {}
for cl, g in sub.obs.groupby(CLUSTER_KEY, observed=True):
    idx = g.index
    r = {'n': len(g)}
    if have_dblt:
        r['dblt_mean'] = round(float(sub.obs.loc[idx, 'doublet_score'].mean()), 3)
        r['dblt_p90'] = round(float(sub.obs.loc[idx, 'doublet_score'].quantile(0.9)), 3)
    if 'predicted_doublet' in sub.obs.columns:
        r['pct_pred_dblt'] = round(100 * sub.obs.loc[idx, 'predicted_doublet'].mean(), 1)
    r['pct_CD3']    = round(100 * CD3pos.loc[idx].mean(), 1)
    r['pct_Bgene']  = round(100 * Bpos.loc[idx].mean(), 1)
    r['pct_MHCII']  = round(100 * MHC2.loc[idx].mean(), 1)
    r['pct_prolif'] = round(100 * Prolif.loc[idx].mean(), 1)
    rows[cl] = r

tbl = pd.DataFrame(rows).T
tbl.index.name = CLUSTER_KEY
print("\nper-subcluster adjudication table:")
print(tbl.to_string())
tbl.to_csv(os.path.join(OUT_DIR, 'tier2_Tcell_v2_cluster_adjudication.csv'))
print(f"\nsaved: {os.path.join(OUT_DIR, 'tier2_Tcell_v2_cluster_adjudication.csv')}")
print("""
read:
  cluster 9 dblt_mean/dblt_p90 ELEVATED vs the others -> residual B doublets;
    quarantine the high-score cells, label the clean remainder.
  cluster 9 dblt FLAT but pct_Bgene high with high pct_MHCII/pct_CD3 -> the CD79
    signal is ambient; cluster 9 is a real MHC-II-high activated T subtype.
""")
