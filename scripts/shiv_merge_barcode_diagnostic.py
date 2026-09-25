#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SHIV Merge — Barcode Join Diagnostic (Library H)
================================================
The merge attached almost no SHIV to cells (H: 458 recovered reads, but only 2
landed on called cells, vs the 27 SHIV+ cells CellRanger itself found). That's a
barcode-space mismatch, not biology.

Prime suspect: NovaSeq X (read IDs start 'LH00973'). The X reads the barcode in
the opposite orientation; CellRanger auto-detects and reverse-complements the cell
barcode in its output, STARsolo (plain whitelist) does not. Both still match the
whitelist at high rate, so nothing looked wrong upstream, but a given cell is
stored as a sequence by CellRanger and its reverse-complement by STARsolo, so the
join key never matches.

This tests it directly: overlap of STARsolo SHIV barcodes with CellRanger cells,
DIRECT vs REVERSE-COMPLEMENT. Whichever wins names the fix.

Run WITH python (sc_pre). Light.
Author: Jake Lehle / Kaushal Lab — June 2026
"""

# %% Cell 1 — Config + loaders
import os, gzip, subprocess
import scipy.io
import scanpy as sc

SOLO_DIR     = "/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV/shiv_starsolo/solo"
ANALYSIS_DIR = "/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV/Analysis"
LIB, NAME = "H", "SHIV_Necropsy_LN"   # richest library = clearest test

_comp = str.maketrans("ACGTN", "TGCAN")
def rc(s): return s.translate(_comp)[::-1]

def shell(c): return subprocess.run(c, shell=True, capture_output=True, text=True).stdout

def load_shiv_barcodes(solo_gene_raw_dir):
    def pick(b):
        for e in ("", ".gz"):
            p = os.path.join(solo_gene_raw_dir, b + e)
            if os.path.exists(p): return p
        return None
    mtx, bcs = pick("matrix.mtx"), pick("barcodes.tsv")
    M = scipy.io.mmread(mtx).tocoo()
    cnt = {}
    for c, v in zip(M.col, M.data):
        cnt[c] = cnt.get(c, 0) + int(v)
    needed = set(cnt)
    out = {}
    opener = gzip.open if bcs.endswith(".gz") else open
    with opener(bcs, "rt") as fh:
        for i, line in enumerate(fh):
            if i in needed:
                out[line.strip().split("\t")[0].split("-")[0]] = cnt[i]
    return out


# %% Cell 2 — Load both barcode sets
solo_raw = os.path.join(SOLO_DIR, f"{LIB}_{NAME}_Solo.out", "Gene", "raw")
shiv = load_shiv_barcodes(solo_raw)                       # {bare_bc: shiv_umi}
print(f"STARsolo SHIV barcodes: {len(shiv):,}  (total umi {sum(shiv.values()):,})")

h5 = shell(f"find -L '{ANALYSIS_DIR}/{NAME}/outs' "
           f"-name 'sample_filtered_feature_bc_matrix.h5' 2>/dev/null").split()[0]
cells = sc.read_10x_h5(h5)
cell_bcs = set(bc.split("-")[0] for bc in cells.obs_names)
print(f"CellRanger called cells:  {len(cell_bcs):,}")


# %% Cell 3 — The decisive test: direct vs reverse-complement overlap
shiv_bcs = set(shiv)
direct = shiv_bcs & cell_bcs
rcomp  = {rc(b) for b in shiv_bcs} & cell_bcs

print("\n" + "=" * 60)
print("OVERLAP OF SHIV BARCODES WITH CELLRANGER CELLS")
print("=" * 60)
print(f"   direct match            : {len(direct):,} / {len(shiv_bcs):,}")
print(f"   reverse-complement match: {len(rcomp):,} / {len(shiv_bcs):,}")

print("\n--- example STARsolo SHIV barcodes (top by umi) ---")
for b, u in sorted(shiv.items(), key=lambda x: -x[1])[:8]:
    print(f"   {b}  umi={u}  rc={rc(b)}  "
          f"{'DIRECT-in-cells' if b in cell_bcs else ''}"
          f"{'  RC-in-cells' if rc(b) in cell_bcs else ''}")
print("\n--- example CellRanger cell barcodes ---")
for b in list(cell_bcs)[:5]:
    print(f"   {b}")

print("\nREAD:")
print("  reverse-complement >> direct  -> NovaSeq X orientation; RC the STARsolo")
print("    barcodes in the merge (or use a translated whitelist). Easy fix.")
print("  direct >> rc                  -> barcodes already aligned; the SHIV")
print("    barcodes are just mostly ambient (not cells) — different conversation.")
print("  both ~0                       -> deeper issue; inspect raw barcode lists.")
