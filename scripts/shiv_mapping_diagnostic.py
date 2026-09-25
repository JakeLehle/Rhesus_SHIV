#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SHIV Mapping Mechanism Diagnostic — v3 (all libraries)
======================================================
Rhesus SHIV scRNA-seq — full cross-library view of SHIV read availability and loss.

WHAT CHANGED FROM v2 (and why):
  v2 nailed the per-sample BAM for Library H and showed the loss is annotation
  (overlapping ORFs + unannotated LTRs), not STAR stringency. Two loose ends:
    (a) Jake wants the picture across ALL 8 libraries, not just Necropsy LN, to
        see if SHIV sparsity is uniform or concentrated (e.g. 21 DPI).
    (b) Library H showed a suspicious perfect 477/477 strand split AND counted
        reads on both strands, which is impossible if only sense reads count.
        Most likely the BAM holds BOTH mates of each pair. v3 decomposes the
        SAM FLAG to settle that directly before we trust any per-read tally.

This script:
  1. Discovers every SHIV_* multi output under the Analysis dir.
  2. For each, pulls SHIV-contig reads from sample_alignments.bam (counted pop)
     and unassigned_alignments.bam (ambient pop).
  3. Decomposes the FLAG (paired / read1 / read2 / proper-pair / reverse) so we
     know whether "reads" are fragments or mates.
  4. Classifies each read by position (LTR / single-gene / multi-gene-overlap)
     using gene intervals parsed live from the GTF.
  5. Emits a per-library summary table + pooled view, and saves the summary TSV.

Counting rule (10x): counted == xf:i:25 (CONF_MAPPED + UMI_COUNT + CONF_FEATURE),
i.e. uniquely one gene, sense, exon-consistent. Antisense and multi-gene reads
are correctly NOT counted.

Run WITH python (not bash). Spyder cells (# %%). Light job.

Author: Jake Lehle / Kaushal Lab
Date: June 2026
"""

# %% Cell 1 — Configuration + helpers
import os
import re
import glob
import shutil
import subprocess
import pandas as pd

ANALYSIS_DIR = "/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV/Analysis"
GEX_REF_DIR  = "/master/jlehle/WORKING/SC/REF/Mmul10_SHIVAD8EO_v10"
GEX_REF_ALT  = "/master/jlehle/cellranger-10.0.0/Mmul10_SHIVAD8EO_v10"
VIRAL_CONTIG = "SHIVAD8EO"
SUMMARY_TSV  = "/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV/shiv_mapping_summary_all_libs.tsv"

SAMTOOLS = None  # None => auto-detect (system one is broken)
WANTED_TAGS = ["NH", "nM", "xf", "GX", "GN", "RE", "CB", "UB"]
HIGH_NM = 9          # threshold to flag the divergent/spurious mismatch tail
_CIGAR_REF_OPS = set("MDN=X")
_CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


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


def cigar_ref_len(c):
    if c == "*":
        return 0
    return sum(int(n) for n, op in _CIGAR_RE.findall(c) if op in _CIGAR_REF_OPS)


print(f"[CONFIG] samtools     : {SAMTOOLS_BIN}")
print(f"[CONFIG] analysis dir : {ANALYSIS_DIR}")
print(f"[CONFIG] contig       : {VIRAL_CONTIG}")


# %% Cell 2 — Parse gene intervals from GTF (for position classification)
ref_dir = GEX_REF_DIR if os.path.isdir(GEX_REF_DIR) else GEX_REF_ALT
gtf_candidates = (glob.glob(os.path.join(ref_dir, "genes", "genes.gtf*")) +
                  glob.glob(os.path.join(ref_dir, "genes", "*.gtf*")))
genes = []
if gtf_candidates:
    gtf = gtf_candidates[0]
    opener = "zcat" if gtf.endswith(".gz") else "cat"
    out = shell(f"{opener} '{gtf}' | awk -F'\\t' '$1==\"{VIRAL_CONTIG}\" && $3==\"gene\"'")
    for line in out.splitlines():
        f = line.split("\t")
        if len(f) >= 9:
            m = re.search(r'gene_name "([^"]+)"', f[8])
            genes.append((m.group(1) if m else "?", int(f[3]), int(f[4]), f[6]))
    genes.sort(key=lambda g: g[1])

MIN_GENE_START = min((g[1] for g in genes), default=None)
MAX_GENE_END   = max((g[2] for g in genes), default=None)
print(f"[GTF] {len(genes)} genes on {VIRAL_CONTIG}; "
      f"span {MIN_GENE_START}-{MAX_GENE_END}")


def classify_position(start, end):
    hits = [g[0] for g in genes if (g[1] <= end and g[2] >= start)]
    if not hits:
        if MIN_GENE_START is not None and end < MIN_GENE_START:
            return "5prime_LTR_UTR"
        if MAX_GENE_END is not None and start > MAX_GENE_END:
            return "3prime_LTR_UTR"
        return "intergenic_gap"
    return f"single:{hits[0]}" if len(hits) == 1 else "multi_gene_overlap"


# %% Cell 3 — Discover libraries + BAMs, parse all SHIV-contig reads
lib_dirs = sorted(glob.glob(os.path.join(ANALYSIS_DIR, "SHIV_*", "outs")))
print(f"\n[DISCOVER] {len(lib_dirs)} SHIV_* output(s):")
for d in lib_dirs:
    print("   ", d)


def parse_sam_line(line):
    f = line.rstrip("\n").split("\t")
    flag = int(f[1])
    start = int(f[3])
    rlen = cigar_ref_len(f[5])
    rec = {
        "flag": flag, "pos": start, "mapq": int(f[4]),
        "end": start + rlen - 1 if rlen else start,
        "is_secondary": bool(flag & 0x100), "is_supplementary": bool(flag & 0x800),
        "is_paired": bool(flag & 0x1), "is_proper_pair": bool(flag & 0x2),
        "is_read1": bool(flag & 0x40), "is_read2": bool(flag & 0x80),
        "is_reverse": bool(flag & 0x10),
    }
    tv = {t: None for t in WANTED_TAGS}
    for fld in f[11:]:
        p = fld.split(":", 2)
        if len(p) == 3 and p[0] in tv:
            tv[p[0]] = p[2]
    rec.update(tv)
    return rec


all_rows = []
for d in lib_dirs:
    cond = os.path.basename(os.path.dirname(d)).replace("SHIV_", "")
    bams = shell(f"find -L '{d}' -name 'sample_alignments.bam' "
                 f"-o -name 'unassigned_alignments.bam' 2>/dev/null").split()
    for b in bams:
        src = "per_sample" if "sample_alignments" in b else "unassigned"
        out = shell(f"{SAMTOOLS_BIN} view '{b}' '{VIRAL_CONTIG}' 2>/dev/null")
        for line in out.splitlines():
            r = parse_sam_line(line)
            r["condition"] = cond
            r["bam_source"] = src
            all_rows.append(r)
    print(f"   [{cond}] parsed")

df = pd.DataFrame(all_rows)
if df.empty:
    raise RuntimeError("No SHIV reads found in any library.")

dp = df[~df["is_secondary"] & ~df["is_supplementary"]].copy()
for c in ["NH", "nM", "xf"]:
    dp[c] = pd.to_numeric(dp[c], errors="coerce")
dp["has_GX"]    = dp["GX"].notna() & (dp["GX"] != "-")
dp["multi_GX"]  = dp["GX"].fillna("").str.contains(";")
dp["counted"]   = (dp["xf"] == 25)
dp["pos_class"] = [classify_position(s, e) for s, e in zip(dp["pos"], dp["end"])]
dp["in_LTR"]    = dp["pos_class"].isin(["5prime_LTR_UTR", "3prime_LTR_UTR"])
dp["high_nM"]   = dp["nM"] >= HIGH_NM


# %% Cell 4 — Resolve the paired/mate question (this gates how we read counts)
print("\n" + "=" * 78)
print("FLAG DECOMPOSITION — are 'reads' fragments or mates?")
print("=" * 78)
n = len(dp)
print(f"   primary alignment records : {n:,}")
print(f"   paired (FLAG 0x1)         : {dp['is_paired'].sum():,}")
print(f"   proper pair (0x2)         : {dp['is_proper_pair'].sum():,}")
print(f"   read1 (0x40)              : {dp['is_read1'].sum():,}")
print(f"   read2 (0x80)              : {dp['is_read2'].sum():,}")
print(f"   reverse-strand (0x10)     : {dp['is_reverse'].sum():,}")
print("\n   read1 x reverse-strand:")
print(pd.crosstab(dp["is_read1"], dp["is_reverse"]))
print("\n   INTERPRETATION GUIDE:")
print("     - if ~half read1 and ~half read2 -> BAM holds BOTH mates;")
print("       a 50/50 strand split is just mate orientation, count FRAGMENTS.")
print("     - if nearly all read2 (or unpaired) -> single aligned read per")
print("       fragment; the strand split is real sense/antisense (assay artifact).")


# %% Cell 5 — Per-library summary table (+ save TSV)
print("\n" + "=" * 78)
print("PER-LIBRARY SUMMARY")
print("=" * 78)


def lib_summary(g):
    ps = g[g["bam_source"] == "per_sample"]
    return pd.Series({
        "total_reads":      len(g),
        "per_sample_reads": len(ps),
        "unassigned_reads": len(g) - len(ps),
        "counted_xf25":     int(g["counted"].sum()),
        "has_gene":         int(g["has_GX"].sum()),
        "multi_gene_GX":    int(g["multi_GX"].sum()),
        "multimapped":      int((g["NH"] > 1).sum()),
        "in_LTR_UTR":       int(g["in_LTR"].sum()),
        "multi_overlap":    int((g["pos_class"] == "multi_gene_overlap").sum()),
        f"high_nM(>={HIGH_NM})": int(g["high_nM"].sum()),
        "fwd_strand":       int((~g["is_reverse"]).sum()),
        "rev_strand":       int(g["is_reverse"].sum()),
    })


summary = dp.groupby("condition").apply(lib_summary).astype(int)
# recoverable-if-annotation-fixed = multi_overlap + LTR reads not already counted
summary["recoverable_est"] = summary["multi_overlap"] + summary["in_LTR_UTR"]
summary = summary.sort_values("counted_xf25", ascending=False)
print(summary.to_string())

summary.to_csv(SUMMARY_TSV, sep="\t")
print(f"\n[SAVED] {SUMMARY_TSV}")

print("\n--- POOLED TOTALS ---")
print(f"   libraries with SHIV reads : {dp['condition'].nunique()}")
print(f"   total SHIV primary reads  : {len(dp):,}")
print(f"   counted molecules (xf25)  : {int(dp['counted'].sum()):,}")
print(f"   recoverable est (overlap+LTR): {int(summary['recoverable_est'].sum()):,}")


# %% Cell 6 — Pooled position + assignment breakdown (where the loss is)
print("\n" + "=" * 78)
print("POOLED — WHERE READS LAND / WHY NOT COUNTED")
print("=" * 78)
print("\nposition class (all libraries):")
print(dp["pos_class"].value_counts(dropna=False).to_string())
print("\nGX gene assignment (top 15; ';' == multi-gene = uncountable):")
print(dp["GX"].value_counts(dropna=False).head(15).to_string())
print("\nxf (25 == counted):")
print(dp["xf"].value_counts(dropna=False).to_string())
print(f"\nmismatch nM >= {HIGH_NM} (divergent/spurious tail): "
      f"{int(dp['high_nM'].sum()):,} / {len(dp):,}")
print("\nDone. The per-library table is the full view Jake asked for:")
print("  - is signal concentrated (e.g. 21DPI / LN) or uniformly sparse?")
print("  - recoverable_est = reads an annotation fix (collapse ORFs + cover LTRs)")
print("    would make countable, separate from signal that simply isn't there.")
