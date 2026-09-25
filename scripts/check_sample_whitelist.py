#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SHIV Chemistry + Whitelist Confirmation
=======================================
Goal: confirm, from the CellRanger outputs themselves, the three things STARsolo
needs to get barcode handling exactly right:
    1. CB length      -> --soloCBlen   (expected 16)
    2. UMI length     -> --soloUMIlen  (10 == 5' v2 ; 12 == GEM-X 5' v3)
    3. the whitelist  -> --soloCBwhitelist  (the exact file CellRanger used)

Why empirical and not from the chemistry name: the CB tags in the BAM were
error-corrected by CellRanger against whatever whitelist it auto-detected, so
the barcodes ARE ground truth. We measure lengths directly and identify the
whitelist by matching observed barcodes against every whitelist file in the
install. The correct file matches ~100% of observed barcodes; wrong ones ~0%.

This resolves the "GEM-X vs v3" ambiguity with data, not product naming.

Run WITH python (not bash). Spyder cells (# %%). Memory-light: observed barcode
set is small and whitelist files are STREAMED (safe on the shared titan node).

Author: Jake Lehle / Kaushal Lab
Date: June 2026
"""

# %% Cell 1 — Configuration + helpers
import os
import re
import gzip
import glob
import shutil
import subprocess

ANALYSIS_DIR   = "/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV/Analysis"
CELLRANGER_DIR = "/master/jlehle/cellranger-10.0.0"

# Libraries to sample barcodes from (chemistry is library-independent, but
# sampling a couple guards against one odd library). Good-signal libs are fine.
SAMPLE_LIBS = ["SHIV_Necropsy_LN", "SHIV_21DPI_PBMC"]

SAMTOOLS = None                 # None => auto-detect (system one is broken)
N_SAMPLE_READS = 200000         # reads to scan per BAM for CB/UMI lengths
N_MATCH_BC     = 3000           # unique observed CBs to test against whitelists

_CB_RE = re.compile(r"CB:Z:(\S+)")
_UB_RE = re.compile(r"UB:Z:(\S+)")


def find_samtools():
    if SAMTOOLS:
        return SAMTOOLS
    cp = os.environ.get("CONDA_PREFIX", "")
    if cp and os.path.exists(os.path.join(cp, "bin", "samtools")):
        return os.path.join(cp, "bin", "samtools")
    f = shutil.which("samtools")
    return f if (f and "/usr/local/bin" not in f) else (f or "samtools")


SAMTOOLS_BIN = find_samtools()


def shell(cmd):
    return subprocess.run(cmd, shell=True, capture_output=True, text=True).stdout


def find_bam(lib):
    d = os.path.join(ANALYSIS_DIR, lib, "outs")
    hits = shell(f"find -L '{d}' -name 'sample_alignments.bam' 2>/dev/null").split()
    return hits[0] if hits else None


print(f"[CONFIG] samtools     : {SAMTOOLS_BIN}")
print(f"[CONFIG] sample libs  : {SAMPLE_LIBS}")


# %% Cell 2 — Detected chemistry string from CellRanger metrics
print("\n" + "=" * 78)
print("DETECTED CHEMISTRY (from CellRanger metrics)")
print("=" * 78)
for lib in SAMPLE_LIBS:
    d = os.path.join(ANALYSIS_DIR, lib, "outs")
    csvs = shell(f"find -L '{d}' -name 'metrics_summary.csv' "
                 f"-o -name 'qc_*metrics*.csv' 2>/dev/null").split()
    print(f"\n--- {lib} ---")
    found = False
    for c in csvs:
        out = shell(f"grep -i -h 'chemistry' '{c}' 2>/dev/null")
        if out.strip():
            found = True
            print(f"   {os.path.basename(c)}:")
            for line in out.strip().splitlines():
                print(f"     {line.strip()}")
    if not found:
        print("   [no chemistry line found in metrics CSVs; rely on Cell 3 lengths]")


# %% Cell 3 — Measure CB and UMI lengths empirically from BAM tags
print("\n" + "=" * 78)
print("EMPIRICAL CB / UMI LENGTHS (from BAM tags)")
print("=" * 78)

observed_cbs = set()      # 16bp barcodes, suffix stripped
cb_len_counts = {}
umi_len_counts = {}
cb_suffixes = {}

for lib in SAMPLE_LIBS:
    bam = find_bam(lib)
    if not bam:
        print(f"   [WARN] no BAM for {lib}")
        continue
    # subsample for diversity, cap with head
    out = shell(f"{SAMTOOLS_BIN} view -s 42.01 '{bam}' 2>/dev/null | "
                f"head -n {N_SAMPLE_READS}")
    n_cb = n_umi = 0
    for line in out.splitlines():
        mcb = _CB_RE.search(line)
        if mcb:
            cb = mcb.group(1)
            # split the GEM-well suffix (e.g. -1)
            core, _, suf = cb.partition("-")
            cb_len_counts[len(core)] = cb_len_counts.get(len(core), 0) + 1
            if suf:
                cb_suffixes[suf] = cb_suffixes.get(suf, 0) + 1
            if len(core) == 16 and len(observed_cbs) < N_MATCH_BC:
                observed_cbs.add(core)
            n_cb += 1
        mub = _UB_RE.search(line)
        if mub:
            ub = mub.group(1)
            umi_len_counts[len(ub)] = umi_len_counts.get(len(ub), 0) + 1
            n_umi += 1
    print(f"   {lib}: scanned, {n_cb:,} CB tags, {n_umi:,} UB tags")

print("\n--- CB core length distribution (expect 16) ---")
for L, c in sorted(cb_len_counts.items()):
    print(f"   len {L:>3} : {c:,}")
print("--- CB GEM-well suffix distribution ---")
for s, c in sorted(cb_suffixes.items()):
    print(f"   -{s:>3} : {c:,}")
print("--- UMI length distribution (10 == 5' v2 ; 12 == GEM-X 5' v3) ---")
for L, c in sorted(umi_len_counts.items()):
    print(f"   len {L:>3} : {c:,}")

CB_LEN = max(cb_len_counts, key=cb_len_counts.get) if cb_len_counts else None
UMI_LEN = max(umi_len_counts, key=umi_len_counts.get) if umi_len_counts else None
if UMI_LEN == 10:
    chem_label = "5' v2"
elif UMI_LEN == 12:
    chem_label = "GEM-X 5' v3"
else:
    chem_label = "??"
print(f"\n   => CB length  : {CB_LEN}")
print(f"   => UMI length : {UMI_LEN}  ({chem_label})")
print(f"   collected {len(observed_cbs):,} unique observed barcodes for matching")


# %% Cell 4 — Locate whitelist files in the CellRanger install
print("\n" + "=" * 78)
print("WHITELIST FILES IN INSTALL")
print("=" * 78)
wl_files = shell(
    f"find -L '{CELLRANGER_DIR}' -ipath '*barcode*' "
    f"\\( -name '*.txt' -o -name '*.txt.gz' \\) 2>/dev/null").split()
# de-dup and drop obviously-tiny files
wl_files = sorted(set(wl_files))
for w in wl_files:
    sz = os.path.getsize(w) if os.path.exists(w) else 0
    print(f"   {sz/1e6:8.1f} MB  {w}")
if not wl_files:
    print("   [WARN] no whitelist files found under a 'barcode' path; widening search")
    wl_files = sorted(set(shell(
        f"find -L '{CELLRANGER_DIR}' -name '3M-*.txt*' -o -name '737K-*.txt*' "
        f"2>/dev/null").split()))
    for w in wl_files:
        print(f"   {w}")


# %% Cell 5 — Identify the whitelist by matching observed barcodes (streamed)
print("\n" + "=" * 78)
print("WHITELIST MATCH (observed barcodes vs each whitelist)")
print("=" * 78)
if not observed_cbs:
    print("   [WARN] no observed barcodes collected; cannot match")
else:
    n_obs = len(observed_cbs)
    results = []
    for w in wl_files:
        try:
            opener = gzip.open if w.endswith(".gz") else open
            hit = 0
            with opener(w, "rt") as fh:
                for line in fh:
                    bc = line.strip().split("\t")[0].split("-")[0]
                    if bc in observed_cbs:
                        hit += 1
            frac = 100.0 * hit / n_obs
            results.append((frac, hit, w))
            print(f"   {frac:6.1f}%  ({hit:,}/{n_obs:,})  {os.path.basename(w)}")
        except Exception as e:
            print(f"   [skip] {os.path.basename(w)}: {e}")
    results.sort(reverse=True)
    if results and results[0][0] > 50:
        best = results[0]
        print(f"\n   => MATCHED WHITELIST: {best[2]}")
        print(f"      ({best[0]:.1f}% of observed barcodes present)")
        MATCHED_WL = best[2]
    else:
        print("\n   [WARN] no whitelist matched well; inspect manually")
        MATCHED_WL = None


# %% Cell 6 — STARsolo-ready barcode config summary
print("\n" + "=" * 78)
print("STARsolo BARCODE CONFIG (derived from the above)")
print("=" * 78)
print(f"   --soloType CB_UMI_Simple")
print(f"   --soloCBstart 1   --soloCBlen {CB_LEN or '??'}")
print(f"   --soloUMIstart {(CB_LEN + 1) if CB_LEN else '??'}   "
      f"--soloUMIlen {UMI_LEN or '??'}")
try:
    print(f"   --soloCBwhitelist {MATCHED_WL or '<<unresolved>>'}")
except NameError:
    print(f"   --soloCBwhitelist <<unresolved>>")
print(f"   --soloStrand Forward   # 5' GEX: R2 cDNA, confirm against recovery")
print(f"   --soloFeatures Gene")
print("\n   NOTE: CB+UMI live on R1, cDNA on R2. The relaxed-stringency flags")
print("   (outFilterScoreMinOverLread / outFilterMatchNminOverLread) get added")
print("   on top of this once the barcode config is confirmed correct.")
print("\nDone. Send Cell 3 (lengths), Cell 5 (matched whitelist), Cell 6.")
