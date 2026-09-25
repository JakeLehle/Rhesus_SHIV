# %% Cell 1 — Parameters (STARsolo / unmapped / MEGAHIT read-pool discovery)
# Read-only orientation pass. Nothing is written or modified.
import os, time, fnmatch

SEARCH_ROOTS = [
    "/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV/shiv_starsolo",
    "/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV",
]

MIN_BYTES   = 1_000_000                       # ignore anything under 1 MB (logs, indices, noise)
EXCLUDE_EXT = {".h5ad", ".h5mu", ".h5", ".png", ".pdf", ".bai", ".idx", ".log"}

# Filename patterns grouped by candidate read pool (case-insensitive).
POOL_PATTERNS = {
    "starsolo_bam":      ["Aligned.sortedByCoord.out.bam", "Aligned.out.bam"],
    "starsolo_solo_out": ["*solo.out*"],
    "recovered_shiv":    ["*shiv*read*", "*shiv*.fastq*", "*shiv*.fq*", "*shiv*.bam",
                          "*shiv*.sam", "*shiv*barcode*"],
    "unmapped_pool":     ["*unmapped.out.mate*", "*unmapped*.fastq*", "*unmapped*.fq*",
                          "*unmapped*.bam"],
    "megahit":           ["final.contigs.fa*", "*megahit*contigs*", "*contigs*.fa*"],
}

SKIP_DIRS = {"_STARtmp", "vdj_reference", "viral_id_databases", ".snakemake"}

# %% Cell 2 — Walk, classify, rank by write time
def human(n):
    for u in ["B", "KB", "MB", "GB", "TB"]:
        if n < 1024: return f"{n:6.1f}{u}"
        n /= 1024
    return f"{n:.1f}PB"

def classify(name):
    if os.path.splitext(name)[1].lower() in EXCLUDE_EXT:
        return None
    low = name.lower()
    for pool, pats in POOL_PATTERNS.items():
        if any(fnmatch.fnmatch(low, p) for p in pats):
            return pool
    return None

hits = {}
for root in SEARCH_ROOTS:
    if not os.path.isdir(root):
        print(f"[warn] root not found, skipping: {root}")
        continue
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = [d for d in dirnames if d not in SKIP_DIRS]
        for fn in filenames:
            pool = classify(fn)
            if pool is None:
                continue
            fp = os.path.join(dirpath, fn)
            try:
                st = os.stat(fp)
            except OSError as e:
                print(f"[warn] cannot stat {fp}: {e}")
                continue
            if st.st_size >= MIN_BYTES:
                hits[fp] = (st.st_mtime, st.st_size, pool, fp)

ranked = sorted(hits.values(), reverse=True)   # newest first

if not ranked:
    print("No candidate read-pool files found under the search roots.")
else:
    print(f"{'MTIME':<18}{'SIZE':>10}  {'POOL':<18}PATH")
    print("-" * 100)
    for m, s, pool, fp in ranked:
        print(f"{time.strftime('%Y-%m-%d %H:%M', time.localtime(m)):<18}"
              f"{human(s):>10}  {pool:<18}{fp}")

    print("\n--- Newest file per pool ---")
    seen = set()
    for m, s, pool, fp in ranked:
        if pool in seen:
            continue
        seen.add(pool)
        print(f"  {pool:<18}{time.strftime('%Y-%m-%d %H:%M', time.localtime(m))}   {fp}")

    m, s, pool, fp = ranked[0]
    print(f"\nMost recently written overall: {pool}  "
          f"({time.strftime('%Y-%m-%d %H:%M', time.localtime(m))})\n  {fp}")
