#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SHIV Recovery -> Metadata-Preserving Merge  (v2, corrected)
===========================================================
FIX (July 2026)
---------------
The prior merge re-read raw CellRanger cells with sc.read_10x_h5() and
concatenated them. That produced a clean host matrix but DROPPED every piece
of pipeline metadata: the recovered object carried only
    obs = [shiv_umi, shiv_pos, library, condition, lib_key]
on 126,252 raw (pre-QC) cells x 40,317 unfiltered genes, with no animal_id,
no timepoint/tissue, no QC metrics, no ADT, no uns. That is why the T/NK
diagnostic's per-animal (Cell 8), ADT (Cell 6), and cluster (Cell 7) panels
all went dark.

This version instead LOADS an existing pipeline object that already carries
the metadata and ATTACHES the recovered SHIV signal onto it as obs, touching
nothing else. The host transcriptome, ADT, QC, and every obs/var/uns/obsm
entry are preserved exactly.

Which object to point at (TARGET_IN)
------------------------------------
Point it at the furthest-downstream object that has what you need:
  adata_qc.h5ad        -> animal_id, tissue, timepoint, subspecies, QC, ADT (raw
                          protein_counts + uns['adt_names']). NO clusters/UMAP/labels.
  adata_pp.h5ad        -> the above + embedding + clusters + CLR ADT (protein_clr).
  adata_annotated.h5ad -> the above + tier1/tier2 labels. NOTE its SHIV columns
                          (shiv_positive / total_counts_shiv) are the CellRanger
                          UNDERCOUNT; the recovered shiv_umi/shiv_pos we attach
                          here supersede them, use these downstream.

Default is adata_qc.h5ad (that's the object you flagged as holding the metadata).
That unblocks the per-animal and ADT diagnostic panels immediately. If you want
the cluster cross-tab and the real saved labels in the diagnostic too, repoint
TARGET_IN at adata_pp.h5ad or adata_annotated.h5ad, it is a one-line change and
the join logic is identical.

Join key
--------
(library, bare_barcode). CellRanger and STARsolo used the same whitelist +
1MM_CR correction, so a recovered read's barcode equals the cell barcode
CellRanger assigned. Library comes from obs['library']; the bare barcode is the
obs_name up to the first '-' ('ACGT-1-A' -> 'ACGT'), which matches the bare
barcode STARsolo writes.

Because the default target is the QC'd object, any recovered SHIV+ cell that QC
removed will not appear here, by design. Retention is reported so nothing is
silently lost (NHP samples are irreplaceable; we want that number explicit).

Run WITH python (sc_pre env). Spyder cells (# %%).

Author: Jake Lehle / Kaushal Lab
Date: July 2026
"""

# %% Cell 1 — Configuration
import os
import gzip
import numpy as np
import pandas as pd
import scipy.io
import scanpy as sc

WORKING_DIR = "/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV"
SOLO_DIR    = os.path.join(WORKING_DIR, "shiv_starsolo", "solo")
OUT_DIR     = os.path.join(WORKING_DIR, "shiv_starsolo", "merged")

# --- TARGET: the object whose metadata we preserve (see header for options) ---
TARGET_IN   = os.path.join(WORKING_DIR, "adata_qc.h5ad")
OUT_H5AD    = os.path.join(OUT_DIR, "shiv_host_merged.h5ad")
SUMMARY_TSV = os.path.join(OUT_DIR, "shiv_merge_summary.tsv")

MIN_SHIV_UMI = 2          # cell is SHIV+ at >= this many recovered UMIs
LIBRARY_COL  = "library"  # obs column giving the library letter on the target

# letter -> CellRanger sample name (== STARsolo prefix stem)
LIB_NAME = {
    "A": "SHIV_Pre_PBMC",             "B": "SHIV_Pre_LN",
    "C": "SHIV_Pre_PBMC_IndianMixed", "D": "SHIV_21DPI_Pre_LN_Mixed",
    "E": "SHIV_21DPI_PBMC",           "F": "SHIV_21DPI_NonInf_LN_Mixed",
    "G": "SHIV_Necropsy_PBMC",        "H": "SHIV_Necropsy_LN",
}
LIBS = ["A", "B", "C", "D", "E", "F", "G", "H"]

os.makedirs(OUT_DIR, exist_ok=True)


# %% Cell 2 — Build the (library, bare_barcode) -> shiv_umi lookup
def load_shiv_counts(solo_gene_raw_dir):
    """Return {bare_barcode: shiv_umi} from a STARsolo raw Gene matrix.
    Memory-light: only nonzero barcodes are pulled from barcodes.tsv."""
    def pick(base):
        for ext in ("", ".gz"):
            p = os.path.join(solo_gene_raw_dir, base + ext)
            if os.path.exists(p):
                return p
        return None
    mtx = pick("matrix.mtx")
    bcs = pick("barcodes.tsv")
    if not mtx or not bcs:
        return None
    M = scipy.io.mmread(mtx).tocoo()   # (n_genes, n_barcodes); n_genes == 1
    counts_by_idx = {}
    for c, v in zip(M.col, M.data):
        counts_by_idx[c] = counts_by_idx.get(c, 0) + int(v)
    needed = set(counts_by_idx)
    idx2bc = {}
    opener = gzip.open if bcs.endswith(".gz") else open
    with opener(bcs, "rt") as fh:
        for i, line in enumerate(fh):
            if i in needed:
                idx2bc[i] = line.strip().split("\t")[0].split("-")[0]
    return {idx2bc[i]: counts_by_idx[i] for i in counts_by_idx}


print("=" * 72)
print("Building (library, barcode) -> SHIV UMI lookup from STARsolo output")
print("=" * 72)

shiv_lookup = {}                 # {(L, bare_bc): umi}
recovered_per_lib = {}           # {L: {bare_bc: umi}} for retention accounting
for L in LIBS:
    name = LIB_NAME[L]
    solo_raw = os.path.join(SOLO_DIR, f"{L}_{name}_Solo.out", "Gene", "raw")
    counts = load_shiv_counts(solo_raw)
    if counts is None:
        print(f"  [WARN] {L} {name}: no STARsolo matrix at {solo_raw}")
        recovered_per_lib[L] = {}
        continue
    recovered_per_lib[L] = counts
    for bc, u in counts.items():
        shiv_lookup[(L, bc)] = u
    n1 = sum(1 for u in counts.values() if u >= 1)
    n2 = sum(1 for u in counts.values() if u >= MIN_SHIV_UMI)
    print(f"  {L} {name:32s}: {len(counts):>6,} SHIV barcodes "
          f"(>=1umi={n1}, >={MIN_SHIV_UMI}umi={n2})")

print(f"\n  total recovered SHIV barcodes across libraries: {len(shiv_lookup):,}")


# %% Cell 3 — Load the metadata-bearing target and attach SHIV (host untouched)
print("\n" + "=" * 72)
print(f"Attaching SHIV onto: {TARGET_IN}")
print("=" * 72)

adata = sc.read_h5ad(TARGET_IN)
print(f"  target loaded: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
print(f"  obs preserved: {list(adata.obs.columns)}")
print(f"  obsm preserved: {list(adata.obsm.keys())}")
print(f"  uns preserved : {list(adata.uns.keys())}")

# Resolve per-cell (library, bare_barcode)
if LIBRARY_COL in adata.obs.columns:
    libs = adata.obs[LIBRARY_COL].astype(str).values
    print(f"  library from obs['{LIBRARY_COL}']")
else:
    # fall back: library is the last '-'-delimited token of the barcode
    libs = np.array([bc.rsplit("-", 1)[-1] for bc in adata.obs_names])
    print(f"  WARN no '{LIBRARY_COL}' column; parsing library from obs_names suffix")

bare = np.array([bc.split("-")[0] for bc in adata.obs_names])

umi = np.array([shiv_lookup.get((l, b), 0) for l, b in zip(libs, bare)], dtype=int)
adata.obs["shiv_umi"] = umi
adata.obs["shiv_pos"] = umi >= MIN_SHIV_UMI

n_pos = int(adata.obs["shiv_pos"].sum())
n_any = int((umi >= 1).sum())
print(f"\n  attached shiv_umi / shiv_pos (>= {MIN_SHIV_UMI} UMI)")
print(f"  SHIV+ cells retained on target: {n_pos:,} (>=1 umi: {n_any:,})")

# Retention accounting: of the recovered SHIV+ barcodes, how many survived QC
# and landed on a target cell? (Only meaningful if the target is QC-filtered.)
print("\n  retention of recovered SHIV+ barcodes onto the target object:")
target_keys = set(zip(libs, bare))
tot_rec2 = tot_match2 = 0
for L in LIBS:
    rec2 = {bc for bc, u in recovered_per_lib[L].items() if u >= MIN_SHIV_UMI}
    match2 = {bc for bc in rec2 if (L, bc) in target_keys}
    tot_rec2 += len(rec2); tot_match2 += len(match2)
    if rec2:
        dropped = len(rec2) - len(match2)
        flag = f"  <-- {dropped} dropped by QC" if dropped else ""
        print(f"    {L}: {len(match2)}/{len(rec2)} recovered SHIV+ retained{flag}")
print(f"    TOTAL: {tot_match2}/{tot_rec2} recovered SHIV+ (>= {MIN_SHIV_UMI} umi) "
      f"retained after QC")
if tot_match2 < tot_rec2:
    print("    NOTE: dropped cells failed QC upstream (counts/mito/doublet). If you")
    print("          need every SHIV+ cell, attach onto a pre-QC object instead.")


# %% Cell 4 — Write merged object + per-library summary
rows = []
for L in LIBS:
    m = (libs == L)
    n_cells = int(m.sum())
    if n_cells == 0:
        continue
    u = umi[m]
    cond = (adata.obs.loc[m, "condition"].iloc[0]
            if "condition" in adata.obs.columns else LIB_NAME[L])
    rows.append({
        "library": L, "condition": cond, "n_cells": n_cells,
        "shiv_cells_>=1umi": int((u >= 1).sum()),
        f"shiv_cells_>={MIN_SHIV_UMI}umi": int((u >= MIN_SHIV_UMI).sum()),
        "pct_shiv_pos": round(100 * (u >= MIN_SHIV_UMI).sum() / n_cells, 3),
        "max_shiv_umi": int(u.max()) if n_cells else 0,
        "total_shiv_umi": int(u.sum()),
    })

adata.write_h5ad(OUT_H5AD)
print(f"\n[SAVED] {OUT_H5AD}")
print(f"        {adata.n_obs:,} cells x {adata.n_vars:,} genes, "
      f"metadata preserved, SHIV in obs['shiv_umi']/['shiv_pos']")

summary = pd.DataFrame(rows).set_index("library")
summary.to_csv(SUMMARY_TSV, sep="\t")
print("\n" + "=" * 72)
print(f"SHIV MERGE SUMMARY  (SHIV+ defined as >= {MIN_SHIV_UMI} UMI)")
print("=" * 72)
print(summary.to_string())
print(f"\n[SAVED] {SUMMARY_TSV}")
print("\nHost transcriptome, ADT, and all metadata unchanged. SHIV lives only in")
print("obs['shiv_umi']/['shiv_pos']; re-threshold from shiv_umi without re-merging.")
