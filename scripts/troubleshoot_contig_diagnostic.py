# %% Cell 1 — probe: where is the time actually going?
import pandas as pd, numpy as np, re, time

DIAMOND = "/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV/viral_id_databases/search_results/diamond_blastx_contigs_vs_lentivirus.tsv"

COLS_DIAMOND = ["qseqid","sseqid","pident","length","mismatch","gapopen",
                "qstart","qend","sstart","send","evalue","bitscore",
                "qlen","slen","qframe","stitle"]

t0 = time.time()
dm = pd.read_csv(DIAMOND, sep="\t", names=COLS_DIAMOND,
                 usecols=["qseqid","pident","length","evalue","bitscore",
                          "qlen","slen","stitle"])
print(f"load time:            {time.time()-t0:6.1f}s")
print(f"rows:                 {len(dm):,}")
print(f"distinct contigs:     {dm.qseqid.nunique():,}")
print(f"distinct stitles:     {dm.stitle.nunique():,}")
print(f"  ^ this is the number that matters. If it's a few thousand,")
print(f"    classifying on uniques is the whole fix.")

# %% Cell 2 — time the SLOW way on a 200k-row sample (do NOT run on all 5.7M)
VIRAL_PROT = re.compile(
    r"\b(gag|pol|env|nef|tat|rev|vif|vpr|vpu|vpx|integrase|protease|"
    r"reverse transcriptase|gp120|gp41|gp160|capsid|matrix|polyprotein|p24|p17)\b", re.I)
HOST_CONTAM = re.compile(
    r"(fab|light chain|heavy chain|immunoglobulin|antibody|\bchain [a-z]\b|"
    r"neutralizing|adp-ribosylation|kinase|\breceptor\b|\bcd4\b|\btcr\b|\bmhc\b)", re.I)

def classify_protein(t):
    if HOST_CONTAM.search(t): return "host_contam"
    if VIRAL_PROT.search(t):  return "viral_protein"
    return "unclassified"

sample = dm.stitle.head(200_000)
t0 = time.time()
_ = sample.map(classify_protein)
per_row = (time.time()-t0)/len(sample)
print(f"\nslow per-row classify: {per_row*1e6:.1f} microsec/row")
print(f"projected for 5.7M rows x2 maps: {per_row*len(dm)*2/60:.1f} min")

# %% Cell 3 — time the FAST way: classify uniques, map back
t0 = time.time()
uniq = pd.Series(dm.stitle.unique())
pmap = dict(zip(uniq, uniq.map(classify_protein)))
res = dm.stitle.map(pmap)
print(f"\nfast classify-on-unique total: {time.time()-t0:.1f}s")
print(f"speedup factor: ~{(per_row*len(dm))/(time.time()-t0):.0f}x")
