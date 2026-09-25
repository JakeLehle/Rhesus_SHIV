#!/usr/bin/env python3
# %% Cell 1 — parameters
import pandas as pd, numpy as np, re, os
from Bio import SeqIO

BASE    = "/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV"
OUTDIR  = f"{BASE}/viral_id_databases/search_results/summary/ref_mapping"
REF_GB  = f"{OUTDIR}/SHIVAD8EO_MN816822.1.gb"
HITS    = f"{OUTDIR}/candidates_vs_SHIVAD8EO.tsv"
MIN_LEN = 80   # bp, same real-length floor as the nucleotide side

COLS = ["qseqid","sseqid","pident","length","mismatch","gapopen",
        "qstart","qend","sstart","send","evalue","bitscore","qlen","slen"]

# %% Cell 2 — build the gene map from the GenBank annotation (exact, traceable)
print("### Parsing SHIVAD8-EO gene map from GenBank ###")
rec = SeqIO.read(REF_GB, "genbank")
genes = []
for feat in rec.features:
    if feat.type in ("gene","CDS","LTR","mat_peptide","misc_feature"):
        name = (feat.qualifiers.get("gene") or feat.qualifiers.get("product")
                or feat.qualifiers.get("note") or [feat.type])[0]
        genes.append((int(feat.location.start)+1, int(feat.location.end),
                      feat.type, name))
gene_df = pd.DataFrame(genes, columns=["start","end","ftype","name"]).sort_values("start")
print(f"  reference length: {len(rec.seq)} bp")
print(f"  annotated features: {len(gene_df)}")
print(gene_df.to_string(index=False))

# Prefer 'gene' features for assignment; fall back to CDS. Build a clean set.
assign_feats = gene_df[gene_df.ftype.isin(["gene","CDS"])].copy()

def assign_gene(sstart, send):
    mid = (sstart + send) / 2.0
    hits = assign_feats[(assign_feats.start <= mid) & (assign_feats.end >= mid)]
    if len(hits) == 0:
        return "intergenic/LTR"
    # smallest containing feature = most specific gene
    return hits.assign(span=hits.end-hits.start).sort_values("span").iloc[0]["name"]

# %% Cell 3 — load hits, normalize subject coords (blast can report send<sstart)
df = pd.read_csv(HITS, sep="\t", names=COLS)
print(f"\n### Hits loaded: {len(df):,} rows | "
      f"{df.qseqid.nunique()} contigs hit SHIVAD8-EO of 706 candidates ###")
df["s_lo"] = df[["sstart","send"]].min(axis=1)
df["s_hi"] = df[["sstart","send"]].max(axis=1)
df = df[df.length >= MIN_LEN].copy()
print(f"  real-length (>={MIN_LEN}bp) hit rows: {len(df):,}")

df["gene"] = [assign_gene(lo, hi) for lo, hi in zip(df.s_lo, df.s_hi)]

# %% Cell 4 — per-contig best hit, then gene distribution (exclusive)
best = df.sort_values("bitscore", ascending=False).drop_duplicates("qseqid")
print("\n### DISTINCT CONTIGS by best-hit gene on SHIVAD8-EO ###")
print(best.gene.value_counts().to_string())

print("\n### Identity distribution per gene (best hit per contig) ###")
for g, sub in best.groupby("gene"):
    print(f"  {g:18s}: {len(sub):3d} contigs | "
          f"pident {sub.pident.min():.1f}-{sub.pident.max():.1f} | "
          f"max aln {sub.length.max()}bp")

# %% Cell 5 — construct coverage: which positions are hit, which are dark
print("\n### SHIVAD8-EO coverage by candidate contigs ###")
covered = np.zeros(len(rec.seq)+1, dtype=bool)
for lo, hi in zip(df.s_lo, df.s_hi):
    covered[int(lo):int(hi)+1] = True
cov_bp = covered.sum()
print(f"  positions covered: {cov_bp:,} / {len(rec.seq):,} "
      f"({100*cov_bp/len(rec.seq):.1f}%)")

print("\n### Per-gene coverage (fraction of each gene hit by any contig) ###")
for _, row in assign_feats.iterrows():
    s, e = int(row.start), int(row.end)
    frac = covered[s:e+1].sum() / (e-s+1)
    print(f"  {row['name']:18s} [{s:5d}-{e:5d}] {e-s+1:5d}bp: {100*frac:5.1f}% covered")

# %% Cell 6 — write the gene-annotated per-contig table
best.to_csv(f"{OUTDIR}/candidates_genemap_best_per_contig.tsv", sep="\t", index=False)
df.to_csv(f"{OUTDIR}/candidates_genemap_all_hits.tsv", sep="\t", index=False)
print(f"\nWrote gene-annotated tables to {OUTDIR}")
