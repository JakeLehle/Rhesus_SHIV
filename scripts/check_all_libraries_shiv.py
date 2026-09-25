#!/usr/bin/env python3
"""
check_all_libraries_shiv.py

Scan all 8 library h5 matrices for non-zero SHIV gene counts.
Also checks for host gene name collisions (TAT vs tat).
"""

import scanpy as sc
import numpy as np
from scipy.sparse import issparse
import os

ANALYSIS_DIR = "/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV/Analysis"

LIBRARIES = [
    ("A", "SHIV_Pre_PBMC",              "Pre PBMC (Chinese)"),
    ("B", "SHIV_Pre_LN",                "Pre LN (Chinese)"),
    ("C", "SHIV_Pre_PBMC_IndianMixed",  "Pre PBMC (Indian + 41862)"),
    ("D", "SHIV_21DPI_Pre_LN_Mixed",    "21DPI LN (3 Chinese) + Pre LN (41862)"),
    ("E", "SHIV_21DPI_PBMC",            "21 DPI PBMC (Chinese)"),
    ("F", "SHIV_21DPI_NonInf_LN_Mixed", "21DPI LN (41861,41862) + NonInf LN (Indian)"),
    ("G", "SHIV_Necropsy_PBMC",         "Necropsy PBMC (Chinese)"),
    ("H", "SHIV_Necropsy_LN",           "Necropsy LN (Chinese)"),
]

# SHIV genes — lowercase (from SHIV GTF) and uppercase TAT (host gene)
SHIV_GENES = ['gag', 'pol', 'vif', 'vpx', 'vpr', 'tat', 'rev', 'vpu', 'env', 'nef']
# TAT (uppercase) is a host gene (tyrosine aminotransferase) — track it separately
HOST_COLLISION = ['TAT']

print("=" * 90)
print("SHIV Gene Counts Across All 8 Libraries")
print("=" * 90)

# Collect results
results = []

for lib_letter, lib_name, description in LIBRARIES:
    # Try sample-level h5 first, then library-level
    h5_sample = f"{ANALYSIS_DIR}/{lib_name}/outs/per_sample_outs/{lib_name}/sample_filtered_feature_bc_matrix.h5"
    h5_lib = f"{ANALYSIS_DIR}/{lib_name}/outs/filtered_feature_bc_matrix.h5"

    h5_path = h5_sample if os.path.exists(h5_sample) else h5_lib
    if not os.path.exists(h5_path):
        print(f"\n[{lib_letter}] {lib_name} — h5 NOT FOUND")
        continue

    adata = sc.read_10x_h5(h5_path, gex_only=False)
    n_cells = adata.shape[0]

    print(f"\n[{lib_letter}] {lib_name} — {description}")
    print(f"    {n_cells} cells × {adata.shape[1]} features")

    lib_result = {"lib": lib_letter, "name": lib_name, "n_cells": n_cells, "genes": {}}

    # Check each SHIV gene
    for gene in SHIV_GENES:
        # Exact match on var_names (case-sensitive)
        if gene in adata.var_names:
            idx = list(adata.var_names).index(gene)
            col = adata.X[:, idx]
            if issparse(col):
                col = col.toarray().flatten()
            n_pos = int((col > 0).sum())
            total = float(col.sum())
            max_c = float(col.max())
            lib_result["genes"][gene] = {"cells": n_pos, "total_umi": total, "max": max_c}
            if n_pos > 0:
                print(f"    {gene:6s}: {n_pos:4d} cells, {total:6.0f} UMIs, max={max_c:.0f}  *** DETECTED ***")
            # Don't print zeros to keep output clean
        else:
            lib_result["genes"][gene] = {"cells": 0, "total_umi": 0, "max": 0}

    # Check host TAT separately
    for gene in HOST_COLLISION:
        if gene in adata.var_names:
            idx = list(adata.var_names).index(gene)
            col = adata.X[:, idx]
            if issparse(col):
                col = col.toarray().flatten()
            n_pos = int((col > 0).sum())
            if n_pos > 0:
                print(f"    {gene:6s}: {n_pos:4d} cells (HOST tyrosine aminotransferase, NOT SHIV)")

    # Count total SHIV+ cells (any SHIV gene > 0)
    shiv_positive = 0
    shiv_idx = [list(adata.var_names).index(g) for g in SHIV_GENES if g in adata.var_names]
    if shiv_idx:
        shiv_mat = adata.X[:, shiv_idx]
        if issparse(shiv_mat):
            shiv_mat = shiv_mat.toarray()
        shiv_positive = int((shiv_mat.sum(axis=1) > 0).sum())

    lib_result["shiv_positive"] = shiv_positive
    if shiv_positive == 0:
        print(f"    → No SHIV+ cells in this library")
    else:
        print(f"    → TOTAL SHIV+ cells (any gene): {shiv_positive}")

    results.append(lib_result)

# ─── Summary table ──────────────────────────────────────────────────────────
print("\n" + "=" * 90)
print("SUMMARY TABLE")
print("=" * 90)
print(f"{'Lib':<4} {'Name':<35} {'Cells':>7} {'SHIV+':>7} {'Rate':>8}  Top genes")
print("-" * 90)

total_cells = 0
total_shiv = 0
for r in results:
    total_cells += r["n_cells"]
    total_shiv += r["shiv_positive"]
    rate = f"{100*r['shiv_positive']/r['n_cells']:.3f}%" if r["n_cells"] > 0 else "N/A"
    top = ", ".join(f"{g}({d['cells']})" for g, d in r["genes"].items() if d["cells"] > 0)
    if not top:
        top = "—"
    print(f"{r['lib']:<4} {r['name']:<35} {r['n_cells']:>7} {r['shiv_positive']:>7} {rate:>8}  {top}")

print("-" * 90)
rate_total = f"{100*total_shiv/total_cells:.4f}%" if total_cells > 0 else "N/A"
print(f"{'':4} {'TOTAL':<35} {total_cells:>7} {total_shiv:>7} {rate_total:>8}")

# ─── Key question: do ANY libraries have non-zero SHIV counts? ──────────────
print("\n" + "=" * 90)
if total_shiv == 0:
    print("RESULT: ZERO SHIV+ cells across ALL 8 libraries.")
    print("")
    print("This means the 46 SHIV+ cells reported by the Scanpy v2 script")
    print("did NOT come from the cellranger count matrices. Possible sources:")
    print("")
    print("  1. The Scanpy script loaded data from a DIFFERENT cellranger run")
    print("     (e.g., an earlier run with a different reference or settings)")
    print("")
    print("  2. The SHIV detection in Scanpy used a different method than")
    print("     checking gene expression counts (e.g., searching raw reads)")
    print("")
    print("  3. The h5 files we're reading here are from a re-run that")
    print("     overwrote the original outputs")
    print("")
    print("NEXT STEP: Check the Scanpy v2 script to see exactly how it")
    print("identified SHIV+ cells and which input files it used.")
else:
    print(f"RESULT: {total_shiv} SHIV+ cells found across {sum(1 for r in results if r['shiv_positive'] > 0)} libraries.")
    print("")
    print("The SHIV genes ARE being detected in the count matrix.")
    print("The zero-BAM-reads issue may be because cellranger strips")
    print("low-confidence alignments from the BAM but keeps the counts.")
    print("")
    print("NEXT STEP: Check raw_molecule_info.h5 for the actual read evidence")
    print("behind these counts, and examine the raw BAM (not per-sample) for")
    print("SHIV-mapped reads.")
