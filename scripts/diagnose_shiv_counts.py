#!/usr/bin/env python3
"""
diagnose_shiv_counts.py

Work backwards from the count matrix to understand where the 46 "SHIV+"
cells came from, given that samtools idxstats shows 0 reads on the
SHIVAD8EO contig.

Questions this script answers:
  1. Which features in the h5 matrix match SHIV gene names?
  2. Are those features on the SHIVAD8EO contig, or are they host genes
     with the same name?
  3. How many cells have non-zero counts for each?
  4. What do the feature IDs/metadata look like?
"""

import scanpy as sc
import numpy as np
import pandas as pd
import os

# ─── Configuration ──────────────────────────────────────────────────────────
ANALYSIS_DIR = "/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV/Analysis"

# All 8 libraries — check all of them
LIBRARIES = [
    "SHIV_Pre_PBMC",              # A
    "SHIV_Pre_LN",                # B
    "SHIV_Pre_PBMC_IndianMixed",  # C
    "SHIV_21DPI_Pre_LN_Mixed",    # D
    "SHIV_21DPI_PBMC",            # E
    "SHIV_21DPI_NonInf_LN_Mixed", # F
    "SHIV_Necropsy_PBMC",         # G
    "SHIV_Necropsy_LN",           # H
]

# Expected SHIV gene names from the GTF
SHIV_GENE_NAMES = ['gag', 'pol', 'vif', 'vpx', 'vpr', 'tat', 'rev', 'vpu', 'env', 'nef']

# ─── Step 1: Inspect the feature reference from one library ─────────────────
print("=" * 70)
print("STEP 1: Inspect feature space in the count matrix")
print("=" * 70)

# Load Library G (Necropsy PBMC) as the reference
h5_path = f"{ANALYSIS_DIR}/SHIV_Necropsy_PBMC/outs/per_sample_outs/SHIV_Necropsy_PBMC/sample_filtered_feature_bc_matrix.h5"

print(f"\nLoading: {h5_path}")
adata = sc.read_10x_h5(h5_path, gex_only=False)
print(f"Matrix shape: {adata.shape[0]} cells × {adata.shape[1]} features")
print(f"Feature types: {adata.var['feature_types'].value_counts().to_dict()}")

# ─── Step 2: Search for SHIV gene names in all features ────────────────────
print("\n" + "=" * 70)
print("STEP 2: Search for SHIV gene names across ALL features")
print("=" * 70)

print(f"\nAll columns in adata.var: {list(adata.var.columns)}")
print(f"\nFirst 5 features:")
print(adata.var.head())
print(f"\nLast 20 features (SHIV genes would be appended at end):")
print(adata.var.tail(20))

# Search by gene name (var_names / index)
print("\n--- Searching var_names (index) for SHIV gene names ---")
for gene in SHIV_GENE_NAMES:
    # Exact match
    exact = adata.var_names[adata.var_names == gene]
    # Case-insensitive partial match
    partial = adata.var_names[adata.var_names.str.contains(gene, case=False, na=False)]

    if len(exact) > 0:
        print(f"  EXACT match: '{gene}' found in var_names")
    elif len(partial) > 0:
        print(f"  PARTIAL match for '{gene}': {list(partial)}")
    else:
        print(f"  NO match for '{gene}'")

# Search in gene_ids column if it exists
if 'gene_ids' in adata.var.columns:
    print("\n--- Searching gene_ids column for SHIV-related IDs ---")
    for gene in SHIV_GENE_NAMES:
        matches = adata.var[adata.var['gene_ids'].str.contains(gene, case=False, na=False)]
        if len(matches) > 0:
            print(f"  '{gene}' found in gene_ids:")
            print(f"    {matches[['gene_ids']].to_dict()}")

# Also search for the contig name directly
print("\n--- Searching all var columns for 'SHIV' or 'MN816822' ---")
for col in adata.var.columns:
    if adata.var[col].dtype == object:
        shiv_hits = adata.var[adata.var[col].str.contains('SHIV|MN816822|shiv', case=False, na=False)]
        if len(shiv_hits) > 0:
            print(f"  Column '{col}': {len(shiv_hits)} features contain SHIV/MN816822")
            print(shiv_hits)

# ─── Step 3: Check for host gene name collisions ───────────────────────────
print("\n" + "=" * 70)
print("STEP 3: Check if SHIV gene names collide with host (Mmul10) genes")
print("=" * 70)

# The key question: does the feature matrix contain genes named 'gag', 'pol',
# etc. that are actually HOST genes, not SHIV genes?
gex_features = adata.var[adata.var['feature_types'] == 'Gene Expression']
adt_features = adata.var[adata.var['feature_types'] == 'Antibody Capture']

print(f"\nGene Expression features: {len(gex_features)}")
print(f"Antibody Capture features: {len(adt_features)}")

print("\n--- SHIV gene name search in GEX features only ---")
found_features = {}
for gene in SHIV_GENE_NAMES:
    # Check both index (gene name) and gene_ids
    in_names = gex_features.index[gex_features.index.str.lower() == gene.lower()]
    in_ids = pd.Index([])
    if 'gene_ids' in gex_features.columns:
        in_ids = gex_features.index[gex_features['gene_ids'].str.lower() == gene.lower()]

    all_matches = in_names.union(in_ids)
    if len(all_matches) > 0:
        for m in all_matches:
            row = gex_features.loc[m]
            gene_id = row.get('gene_ids', 'N/A')
            print(f"  {gene:6s} → var_name='{m}', gene_id='{gene_id}'")
            found_features[m] = row
    else:
        print(f"  {gene:6s} → NOT FOUND in GEX features")

# ─── Step 4: Count cells with expression for each found feature ────────────
print("\n" + "=" * 70)
print("STEP 4: Cell counts with non-zero expression for SHIV-named features")
print("=" * 70)

if len(found_features) > 0:
    from scipy.sparse import issparse

    for feat_name, feat_info in found_features.items():
        idx = list(adata.var_names).index(feat_name)
        col = adata.X[:, idx]
        if issparse(col):
            col = col.toarray().flatten()
        else:
            col = np.asarray(col).flatten()

        n_nonzero = (col > 0).sum()
        total_counts = col.sum()
        max_count = col.max()

        gene_id = feat_info.get('gene_ids', 'N/A')
        print(f"  {feat_name:15s} (ID: {gene_id})")
        print(f"    Cells with >0 counts: {n_nonzero}")
        print(f"    Total UMI counts:     {total_counts:.0f}")
        print(f"    Max count in a cell:  {max_count:.0f}")
        print()
else:
    print("  No SHIV-named features found in GEX — checking if the Scanpy")
    print("  pipeline used a different method to identify SHIV+ cells.")

# ─── Step 5: Check the raw h5 at the library level too ─────────────────────
print("\n" + "=" * 70)
print("STEP 5: Cross-check with library-level filtered matrix")
print("=" * 70)

lib_h5 = f"{ANALYSIS_DIR}/SHIV_Necropsy_PBMC/outs/filtered_feature_bc_matrix.h5"
if os.path.exists(lib_h5):
    adata_lib = sc.read_10x_h5(lib_h5, gex_only=False)
    print(f"Library-level matrix: {adata_lib.shape[0]} cells × {adata_lib.shape[1]} features")

    # Same search
    gex_lib = adata_lib.var[adata_lib.var['feature_types'] == 'Gene Expression']
    print(f"GEX features: {len(gex_lib)}")

    for gene in SHIV_GENE_NAMES:
        in_names = gex_lib.index[gex_lib.index.str.lower() == gene.lower()]
        if len(in_names) > 0:
            for m in in_names:
                idx = list(adata_lib.var_names).index(m)
                col = adata_lib.X[:, idx]
                if issparse(col):
                    col = col.toarray().flatten()
                n_pos = (col > 0).sum()
                gene_id = gex_lib.loc[m].get('gene_ids', 'N/A')
                print(f"  {m:15s} (ID: {gene_id}): {n_pos} cells with >0 counts")
else:
    print(f"  Library-level h5 not found: {lib_h5}")

# ─── Step 6: Check the GTF for gene ID context ─────────────────────────────
print("\n" + "=" * 70)
print("STEP 6: Look up SHIV gene names in the reference GTF")
print("=" * 70)

gtf_path = "/master/jlehle/cellranger-10.0.0/Mmul10_SHIVAD8EO_v10/genes/genes.gtf"
if os.path.exists(gtf_path):
    print(f"Scanning GTF: {gtf_path}")
    print("(Looking for lines with SHIV gene names...)\n")

    for gene in SHIV_GENE_NAMES:
        # grep-like search in GTF
        hits = []
        with open(gtf_path) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                # Check if this gene name appears in the attributes
                if f'gene_name "{gene}"' in line or f'gene_id "{gene}"' in line:
                    fields = line.strip().split('\t')
                    chrom = fields[0]
                    feature_type = fields[2]
                    start = fields[3]
                    end = fields[4]
                    attrs = fields[8]
                    hits.append((chrom, feature_type, start, end, attrs[:120]))

        if hits:
            print(f"  {gene}: {len(hits)} GTF entries")
            # Show unique chromosomes
            chroms = set(h[0] for h in hits)
            print(f"    Chromosomes: {chroms}")
            # Show first gene-level entry
            gene_entries = [h for h in hits if h[1] == 'gene']
            for h in gene_entries[:3]:
                print(f"    {h[0]}:{h[2]}-{h[3]} ({h[1]}) {h[4]}")
        else:
            print(f"  {gene}: NOT in GTF")
else:
    print(f"  GTF not found at: {gtf_path}")
    print("  Try: find /master/jlehle/cellranger-10.0.0/Mmul10_SHIVAD8EO_v10 -name '*.gtf'")

# ─── Summary ────────────────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
print("""
Key questions answered:
  1. Are SHIV gene names present in the feature matrix?
  2. If yes — are they on the SHIVAD8EO contig or a host chromosome?
  3. How many cells had non-zero counts for each?

If the gene names are on HOST chromosomes:
  → The 46 'SHIV+' cells were a false positive from name collision.
  → Host genes named gag/pol/env/etc exist in mammalian genomes
     (e.g., endogenous retroviruses, retrotransposons).
  → The Scanpy script needs to filter by chromosome/contig, not just gene name.

If the gene names are on the SHIVAD8EO contig:
  → The counts are real but reads were filtered from the BAM
     (cellranger multi may strip low-confidence alignments from the
      per-sample BAM while keeping counts in the matrix).
  → Need to check the library-level BAMs or the raw molecule_info.h5.
""")
