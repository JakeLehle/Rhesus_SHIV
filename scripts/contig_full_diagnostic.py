#!/usr/bin/env python3
# %% Cell 1 — parameters (all adjustable knobs at the top)
import pandas as pd
import numpy as np
import re
import os
import time

DB_ROOT  = "/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV/viral_id_databases"
RESULTS  = f"{DB_ROOT}/search_results"
BLASTN   = f"{RESULTS}/blastn_contigs_vs_lentivirus.tsv"
DIAMOND  = f"{RESULTS}/diamond_blastx_contigs_vs_lentivirus.tsv"
OUTDIR   = f"{RESULTS}/summary"

BLASTN_REAL_MIN_LEN   = 80
BLASTN_DIVERGENT_HI   = 95.0
BLASTN_DIVERGENT_LO   = 70.0

DIAMOND_REAL_MIN_LEN  = 40       # aa
DIAMOND_DIVERGENT_HI  = 90.0
DIAMOND_DIVERGENT_LO  = 30.0
DIAMOND_IDENTITY_FLOOR = 30.0    # NEW: drop sub-30% protein noise before summarizing

os.makedirs(OUTDIR, exist_ok=True)

COLS_BLASTN = ["qseqid","sseqid","pident","length","mismatch","gapopen",
               "qstart","qend","sstart","send","evalue","bitscore",
               "qlen","slen","staxids","stitle"]
COLS_DIAMOND = ["qseqid","sseqid","pident","length","mismatch","gapopen",
                "qstart","qend","sstart","send","evalue","bitscore",
                "qlen","slen","qframe","stitle"]

# %% Cell 2 — lineage classifier
def classify_lineage(title):
    t = title.lower()
    if "simian-human" in t or re.search(r"\bshiv\b", t):
        return "SHIV"
    if "simian immunodeficiency" in t:
        return "SIV"
    if "immunodeficiency virus 2" in t or "hiv-2" in t:
        return "HIV-2"
    if "hiv-1" in t or "human immunodeficiency virus 1" in t:
        return "HIV-1"
    return "other"

# %% Cell 3 — protein-class classifier + NEW gene-region parser
VIRAL_PROT = re.compile(
    r"\b(gag|pol|env|nef|tat|rev|vif|vpr|vpu|vpx|"
    r"integrase|protease|reverse transcriptase|gp120|gp41|gp160|"
    r"capsid|matrix|polyprotein|p24|p17)\b", re.I)
HOST_CONTAM = re.compile(
    r"(fab|light chain|heavy chain|immunoglobulin|antibody|"
    r"\bchain [a-z]\b|neutralizing|adp-ribosylation|kinase|"
    r"\breceptor\b|\bcd4\b|\btcr\b|\bmhc\b|fragment antigen)", re.I)

def classify_protein(title):
    if HOST_CONTAM.search(title):
        return "host_contam"
    if VIRAL_PROT.search(title):
        return "viral_protein"
    return "unclassified"

# NEW: which viral gene a hit points to, parsed from the title.
# Env-family terms collapse to env; pol-family (integrase/protease/RT) to pol.
def classify_gene(title):
    t = title.lower()
    if re.search(r"\b(env|gp120|gp41|gp160)\b", t):                 return "env"
    if re.search(r"\b(pol|integrase|protease|reverse transcriptase)\b", t): return "pol"
    if re.search(r"\b(gag|capsid|matrix|p24|p17)\b", t):            return "gag"
    if re.search(r"\bnef\b", t):  return "nef"
    if re.search(r"\bvif\b", t):  return "vif"
    if re.search(r"\bvpr\b", t):  return "vpr"
    if re.search(r"\bvpx\b", t):  return "vpx"
    if re.search(r"\bvpu\b", t):  return "vpu"
    if re.search(r"\btat\b", t):  return "tat"
    if re.search(r"\brev\b", t):  return "rev"
    if re.search(r"\bpolyprotein\b", t): return "polyprotein"
    return "other"

def map_on_uniques(series, fn):
    uniq = pd.Series(series.unique())
    lut = dict(zip(uniq, uniq.map(fn)))
    return series.map(lut)

# %% Cell 4 — load blastn (UNCHANGED — this side reads correctly)
print("="*60)
t0 = time.time()
bn = pd.read_csv(BLASTN, sep="\t", names=COLS_BLASTN, dtype={"staxids":str})
print(f"blastn loaded in {time.time()-t0:.1f}s | rows: {len(bn):,} | "
      f"distinct contigs: {bn.qseqid.nunique():,}")
bn["lineage"]     = map_on_uniques(bn.stitle, classify_lineage)
bn["is_real_len"] = bn.length >= BLASTN_REAL_MIN_LEN
bn["is_divergent"] = (bn.pident >= BLASTN_DIVERGENT_LO) & (bn.pident < BLASTN_DIVERGENT_HI)

# %% Cell 5 — blastn row-level (UNCHANGED)
print("\n### BLASTN: all hit ROWS by lineage ###")
print(bn.lineage.value_counts().to_string())
print("\n### BLASTN: alignment-length distribution (all rows) ###")
print(bn.length.describe().to_string())
short = (bn.length < BLASTN_REAL_MIN_LEN).sum()
print(f"  rows < {BLASTN_REAL_MIN_LEN}bp (short-motif noise): {short:,} ({100*short/len(bn):.1f}%)")
print(f"\n### BLASTN: REAL-length rows (>={BLASTN_REAL_MIN_LEN}bp) by lineage ###")
print(bn[bn.is_real_len].lineage.value_counts().to_string())

# %% Cell 6 — blastn distinct contigs (UNCHANGED)
print("\n### BLASTN: DISTINCT CONTIGS with a real-length hit, by lineage ###")
real = bn[bn.is_real_len]
for lin in ["SHIV","SIV","HIV-1","HIV-2","other"]:
    sub = real[real.lineage == lin]
    n = sub.qseqid.nunique()
    if n == 0:
        print(f"  {lin:7s}: 0 contigs"); continue
    div = sub[sub.is_divergent].qseqid.nunique()
    print(f"  {lin:7s}: {n:4d} contigs | {div} in divergent band | "
          f"pident {sub.pident.min():.1f}-{sub.pident.max():.1f} | max aln {sub.length.max()}bp")

# %% Cell 7 — blastn best hit per contig (UNCHANGED)
print("\n### BLASTN: best hit per contig (by bitscore), lineage of best hit ###")
best_bn = bn.sort_values("bitscore", ascending=False).drop_duplicates("qseqid")
print(best_bn.lineage.value_counts().to_string())

# %% Cell 8 — load DIAMOND, classify, then REDUCE to one best hit per contig
print("\n" + "="*60)
t0 = time.time()
dm = pd.read_csv(DIAMOND, sep="\t", names=COLS_DIAMOND,
                 usecols=["qseqid","pident","length","evalue","bitscore",
                          "qlen","slen","stitle"])
print(f"DIAMOND loaded in {time.time()-t0:.1f}s | rows: {len(dm):,} | "
      f"distinct contigs: {dm.qseqid.nunique():,}")

# FIX 3: drop sub-floor identity noise BEFORE any summary
before = len(dm)
dm = dm[dm.pident >= DIAMOND_IDENTITY_FLOOR]
print(f"  dropped {before-len(dm):,} rows below {DIAMOND_IDENTITY_FLOOR}% identity "
      f"({len(dm):,} remain)")

t0 = time.time()
dm["protclass"] = map_on_uniques(dm.stitle, classify_protein)
dm["lineage"]   = map_on_uniques(dm.stitle, classify_lineage)
dm["gene"]      = map_on_uniques(dm.stitle, classify_gene)   # NEW
dm["is_real_len"] = dm.length >= DIAMOND_REAL_MIN_LEN
print(f"  classified in {time.time()-t0:.1f}s")

# FIX 1 + FIX 2: collapse to ONE row per contig = its best hit by bitscore.
# Every per-contig summary below reads off this, so counts are exclusive.
dmin = dm[dm.is_real_len].sort_values("bitscore", ascending=False).drop_duplicates("qseqid")
print(f"  reduced to best-hit-per-contig: {len(dmin)} contigs "
      f"(real-length >={DIAMOND_REAL_MIN_LEN}aa)")

# %% Cell 9 — DIAMOND protein class, now per-contig EXCLUSIVE (FIX 2)
print("\n### DIAMOND: DISTINCT CONTIGS by DOMINANT (best-hit) protein class ###")
print(dmin.protclass.value_counts().to_string())
print("  ^ each contig counted once, by its single best hit")

# %% Cell 10 — DIAMOND viral-protein contigs by EXCLUSIVE best-hit lineage (FIX 1)
print("\n### DIAMOND: VIRAL-PROTEIN contigs by best-hit lineage (exclusive) ###")
viral = dmin[dmin.protclass == "viral_protein"]
print(f"  total viral-protein contigs: {len(viral)}")
print(viral.lineage.value_counts().to_string())

print("\n### DIAMOND: viral-protein contigs by best-hit GENE (exclusive) ###")  # NEW
print(viral.gene.value_counts().to_string())

print("\n### DIAMOND: identity range within each best-hit lineage ###")
for lin, g in viral.groupby("lineage"):
    print(f"  {lin:7s}: {len(g):3d} contigs | pident {g.pident.min():.1f}-{g.pident.max():.1f} "
          f"| max aln {g.length.max()}aa")

# %% Cell 11 — cross-method (set logic, already sound — kept)
print("\n### CROSS-METHOD: viral-protein contigs NOT in blastn real-length set ###")
bn_real_contigs  = set(real.qseqid.unique())
dm_viral_contigs = set(viral.qseqid.unique())
diamond_only = dm_viral_contigs - bn_real_contigs
both = dm_viral_contigs & bn_real_contigs
print(f"  viral-protein contigs in DIAMOND: {len(dm_viral_contigs)}")
print(f"  also found by blastn (real-len):  {len(both)}")
print(f"  DIAMOND-ONLY (nt diverged, protein conserved): {len(diamond_only)}")

# NEW: gene breakdown of the DIAMOND-only divergence set specifically
print("\n### Gene breakdown of the DIAMOND-ONLY divergent contigs ###")
print(viral[viral.qseqid.isin(diamond_only)].gene.value_counts().to_string())

# %% Cell 12 — write candidate lists + gene-annotated best-hit tables
print("\n### Writing outputs ###")
candidates = sorted(dm_viral_contigs | bn_real_contigs)
with open(f"{OUTDIR}/candidate_viral_contigs.txt","w") as f:
    f.write("\n".join(candidates) + "\n")
print(f"  {len(candidates)} candidate contigs -> candidate_viral_contigs.txt")

best_bn.to_csv(f"{OUTDIR}/blastn_best_per_contig.tsv", sep="\t", index=False)
dmin.to_csv(f"{OUTDIR}/diamond_best_per_contig.tsv", sep="\t", index=False)
# the divergence set, gene-annotated, for the reference-build decision
viral[viral.qseqid.isin(diamond_only)].to_csv(
    f"{OUTDIR}/diamond_only_divergent_contigs.tsv", sep="\t", index=False)
print(f"  best-hit + divergence tables written to {OUTDIR}")

print("\n" + "="*60)
print("Diagnostic complete.")
print("="*60)
