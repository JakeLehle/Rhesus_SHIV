#!/usr/bin/env python3
"""
bd_abseq_summarize.py
=====================
Step 3 of 3. Read-only. Turns the raw grep output from bd_abseq_scan.sh into
tables and a verdict.

Nothing here filters the inputs. Every hit line the scan produced is parsed,
every probe is located in every hit read, and the filtering happens on the
evidence side, at the end, where you can see it.

Outputs (into OUT_DIR):
  scan_totals.csv          file, lib_type, n_reads, n_hit_reads, hits_per_million
  bd_hits_detail.csv       one row per (read, probe, offset), with flanks
  bd_hits_summary.csv      file x probe: n_hits, n_confident, offset mode
  bd_offset_histogram.csv  file x probe x offset: n
  ctrl_summary.csv         file x control probe: n_hits, offset mode, pct at mode
  verdict.txt              plain-language readout of the three hypotheses
"""

# %% Cell 1 - Parameters
# =============================================================================
import os
import csv
import glob
from collections import Counter, defaultdict

OUT_DIR = "/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV/adt_probe_scan"
RAW_DIR = os.path.join(OUT_DIR, "raw")

# A read is CONFIDENT for a target+orientation if it carries the full 36-mer,
# or >= MIN_ANCHORS distinct tiled 20-mers of that target+orientation. Two
# overlapping anchors imply >= 28bp of contiguous match, which no chance hit
# in a mammalian read pool will produce.
MIN_ANCHORS = 2

# Flank length pulled either side of a confident match. The 5' flank of a
# forward hit (or the 3' flank of a reverse hit, revcomped) is where BD's PCR
# handle lives, if any BD read exists at all.
FLANK = 24

# Expected chance hits are computed against this, purely for the report.
ANCHOR_LEN = 20

MANIFEST = os.path.join(OUT_DIR, "scan_manifest.tsv")
PROBE_MAP = os.path.join(OUT_DIR, "probe_map.csv")

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s):
    return s.translate(_COMP)[::-1]


# %% Cell 2 - Load probes and manifest
# =============================================================================
with open(PROBE_MAP) as fh:
    probes = list(csv.DictReader(fh))
for p in probes:
    p["length"] = int(p["length"])

bd_probes = [p for p in probes if p["kind"].startswith("bd_")]
ctrl_probes = [p for p in probes if p["kind"] == "ctrl"]

with open(MANIFEST) as fh:
    manifest = list(csv.DictReader(fh, delimiter="\t"))
meta = {m["tag"]: m for m in manifest}

n_anchor_probes = sum(1 for p in bd_probes if p["kind"] == "bd_anchor")
print(f"BD probes   : {len(bd_probes)}  ({n_anchor_probes} anchors)")
print(f"Ctrl probes : {len(ctrl_probes)}")
print(f"Files       : {len(manifest)}")


# %% Cell 3 - Totals
# =============================================================================
totals = {}
for tag in meta:
    n_path = os.path.join(RAW_DIR, f"{tag}.nreads")
    n_reads = None
    if os.path.exists(n_path):
        txt = open(n_path).read().strip()
        n_reads = int(txt) if txt.isdigit() else None
    hits_path = os.path.join(RAW_DIR, f"{tag}.bdhits.txt")
    n_hit_reads = sum(1 for _ in open(hits_path)) if os.path.exists(hits_path) else 0
    totals[tag] = {"n_reads": n_reads, "n_hit_reads": n_hit_reads}

with open(os.path.join(OUT_DIR, "scan_totals.csv"), "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["tag", "lib_type", "n_reads", "n_hit_reads", "hits_per_million"])
    for tag, t in totals.items():
        hpm = ""
        if t["n_reads"]:
            hpm = round(1e6 * t["n_hit_reads"] / t["n_reads"], 4)
        w.writerow([tag, meta[tag]["lib_type"], t["n_reads"], t["n_hit_reads"], hpm])

print("\n--- Totals ---")
for tag, t in totals.items():
    nr = f"{t['n_reads']:,}" if t["n_reads"] else "NA"
    print(f"  {tag:<14s} {meta[tag]['lib_type']:<13s} {nr:>15s} reads   {t['n_hit_reads']:>7,} hit reads")


# %% Cell 4 - Locate every probe in every BD hit read
# =============================================================================
def find_all(hay, needle):
    out, i = [], hay.find(needle)
    while i != -1:
        out.append(i)
        i = hay.find(needle, i + 1)
    return out


detail = []
offset_hist = defaultdict(Counter)     # (tag, probe_name) -> Counter(offset)
per_probe = defaultdict(int)           # (tag, probe_name) -> n hits
confident_reads = defaultdict(int)     # (tag, target, orientation) -> n reads

for tag in meta:
    hits_path = os.path.join(RAW_DIR, f"{tag}.bdhits.txt")
    if not os.path.exists(hits_path):
        continue
    with open(hits_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if ":" not in line:
                continue
            idx, seq = line.split(":", 1)

            found = []                                  # (probe, offset)
            by_group = defaultdict(set)                 # (target, orient) -> {probe_name}
            has_full = set()                            # (target, orient)
            for p in bd_probes:
                for off in find_all(seq, p["sequence"]):
                    found.append((p, off))
                    per_probe[(tag, p["probe_name"])] += 1
                    offset_hist[(tag, p["probe_name"])][off] += 1
                    if p["kind"] == "bd_anchor":
                        by_group[(p["target"], p["orientation"])].add(p["probe_name"])
                    else:
                        has_full.add((p["target"], p["orientation"]))

            conf_groups = set(has_full)
            for grp, names in by_group.items():
                if len(names) >= MIN_ANCHORS:
                    conf_groups.add(grp)
            for grp in conf_groups:
                confident_reads[(tag, grp[0], grp[1])] += 1

            for p, off in found:
                grp = (p["target"], p["orientation"])
                lo = max(0, off - FLANK)
                flank5 = seq[lo:off]
                flank3 = seq[off + p["length"] : off + p["length"] + FLANK]
                detail.append(
                    {
                        "tag": tag,
                        "lib_type": meta[tag]["lib_type"],
                        "read_index": idx,
                        "probe_name": p["probe_name"],
                        "target": p["target"],
                        "orientation": p["orientation"],
                        "kind": p["kind"],
                        "offset": off,
                        "read_len": len(seq),
                        "confident": int(grp in conf_groups),
                        "flank5": flank5,
                        "flank3": flank3,
                        "read": seq,
                    }
                )

with open(os.path.join(OUT_DIR, "bd_hits_detail.csv"), "w", newline="") as fh:
    cols = ["tag", "lib_type", "read_index", "probe_name", "target", "orientation",
            "kind", "offset", "read_len", "confident", "flank5", "flank3", "read"]
    w = csv.DictWriter(fh, fieldnames=cols)
    w.writeheader()
    w.writerows(detail)

with open(os.path.join(OUT_DIR, "bd_offset_histogram.csv"), "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["tag", "probe_name", "offset", "n"])
    for (tag, pname), c in sorted(offset_hist.items()):
        for off, n in sorted(c.items()):
            w.writerow([tag, pname, off, n])

with open(os.path.join(OUT_DIR, "bd_hits_summary.csv"), "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["tag", "lib_type", "probe_name", "target", "orientation", "kind",
                "n_hits", "offset_mode", "n_distinct_offsets"])
    for p in bd_probes:
        for tag in meta:
            n = per_probe.get((tag, p["probe_name"]), 0)
            c = offset_hist.get((tag, p["probe_name"]), Counter())
            mode = c.most_common(1)[0][0] if c else ""
            w.writerow([tag, meta[tag]["lib_type"], p["probe_name"], p["target"],
                        p["orientation"], p["kind"], n, mode, len(c)])

print(f"\n--- BD probe hits ---")
print(f"  detail rows : {len(detail):,}")
if confident_reads:
    for (tag, target, orient), n in sorted(confident_reads.items()):
        print(f"  CONFIDENT  {tag:<14s} {target}  orientation {orient}  {n:,} reads")
else:
    print("  No confident BD reads in any file.")


# %% Cell 5 - Positive controls
# =============================================================================
ctrl_by_seq = {p["sequence"]: p for p in ctrl_probes}
ctrl_rows = []

for tag in [t for t in meta if meta[t]["lib_type"] == "PROTEIN"]:
    counts_path = os.path.join(RAW_DIR, f"{tag}.ctrlcounts.tsv")
    hits_path = os.path.join(RAW_DIR, f"{tag}.ctrlhits.txt")
    if not os.path.exists(counts_path):
        continue

    offs = defaultdict(Counter)
    if os.path.exists(hits_path):
        with open(hits_path) as fh:
            for line in fh:
                if ":" not in line:
                    continue
                _, seq = line.rstrip("\n").split(":", 1)
                for cseq, p in ctrl_by_seq.items():
                    for off in find_all(seq, cseq):
                        offs[p["probe_name"]][off] += 1

    with open(counts_path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 4:
                continue
            _, seq, n, n_sub = parts
            p = ctrl_by_seq.get(seq)
            if p is None:
                continue
            c = offs.get(p["probe_name"], Counter())
            if c:
                mode, mode_n = c.most_common(1)[0]
                pct_at_mode = round(100 * mode_n / sum(c.values()), 1)
            else:
                mode, pct_at_mode = "", ""
            ctrl_rows.append(
                {"tag": tag, "probe_name": p["probe_name"], "target": p["target"],
                 "orientation": p["orientation"], "n_hits": int(n),
                 "n_reads_sampled": int(n_sub), "offset_mode": mode,
                 "pct_at_mode": pct_at_mode}
            )

with open(os.path.join(OUT_DIR, "ctrl_summary.csv"), "w", newline="") as fh:
    cols = ["tag", "probe_name", "target", "orientation", "n_hits",
            "n_reads_sampled", "offset_mode", "pct_at_mode"]
    w = csv.DictWriter(fh, fieldnames=cols)
    w.writeheader()
    w.writerows(ctrl_rows)

print("\n--- Positive controls (PROTEIN R2 subsample) ---")
for r in ctrl_rows:
    if r["orientation"] != "F":
        continue
    print(f"  {r['tag']:<14s} {r['probe_name']:<14s} {r['n_hits']:>8,} hits / "
          f"{r['n_reads_sampled']:,} reads   offset mode {r['offset_mode']} "
          f"({r['pct_at_mode']}% of hits)")


# %% Cell 6 - Verdict
# =============================================================================
lines = []


def say(s=""):
    print(s)
    lines.append(s)


ctrl_fwd = [r for r in ctrl_rows if r["orientation"] == "F"]
ctrl_ok = sum(r["n_hits"] for r in ctrl_fwd) > 0
ctrl_modes = Counter(r["offset_mode"] for r in ctrl_fwd if r["n_hits"] > 0 and r["offset_mode"] != "")

prot_conf = sum(n for (tag, _, _), n in confident_reads.items() if meta[tag]["lib_type"] == "PROTEIN")
gex_conf = sum(n for (tag, _, _), n in confident_reads.items() if meta[tag]["lib_type"] == "GEX")
und_conf = sum(n for (tag, _, _), n in confident_reads.items() if meta[tag]["lib_type"] == "UNDETERMINED")

say("=" * 68)
say("VERDICT")
say("=" * 68)

if not ctrl_ok:
    say("CONTROLS DID NOT FIRE. Nothing below is interpretable. The TotalSeq-C")
    say("15-mers should be abundant in PROTEIN R2. Check the manifest paths and")
    say("that Stage 2 ran on the right files before reading anything else.")
else:
    obs_offset = ctrl_modes.most_common(1)[0][0] if ctrl_modes else "?"
    say(f"Controls fired. TotalSeq-C barcodes sit at offset {obs_offset} in PROTEIN R2.")
    say(f"  feature_reference.csv encodes this as 5P + {obs_offset} N's + (BC).")
    say(f"  The BD rows currently carry 36 N's, i.e. offset 36. That is wrong")
    say(f"  independent of everything below, and should be corrected to {obs_offset}.")
    say("")

    if prot_conf > 0:
        say(f"(a) PATTERN BUG supported. {prot_conf:,} confident BD reads in PROTEIN R2.")
        say("    The tags were sequenced. Fix the offset in feature_reference.csv,")
        say("    re-run cellranger multi counting, and CD163/CD169 come back.")
    else:
        say("(a) PATTERN BUG not supported. Zero confident BD reads in PROTEIN R2,")
        say("    in a file where the TotalSeq-C controls are abundant. No pattern")
        say("    change can recover reads that were never sequenced.")
    say("")

    if gex_conf > 0 or und_conf > 0:
        say(f"(c) GEX CARRYOVER supported. {gex_conf:,} confident BD reads in GEX R2, "
            f"{und_conf:,} in Undetermined.")
        say("    Check bd_offset_histogram.csv: a FLAT offset distribution means")
        say("    fragmented cDNA, i.e. real carryover, and these reads carry a valid")
        say("    CB + UMI in R1 at the same read index. Recovery path is a STARsolo")
        say("    run against a two-contig reference holding just the CD163 and CD169")
        say("    36-mers, with the 5' whitelist, giving per-cell counts directly.")
        say("    A TIGHT spike instead means something is anchoring them, which would")
        say("    be surprising in a GEX library and should be looked at before trusting.")
        say("    Also read flank5/flank3 in bd_hits_detail.csv: the bases adjacent to")
        say("    the ABC give us BD's PCR handle for free, which becomes a shared probe.")
    else:
        say("(c) GEX CARRYOVER not supported. Zero confident BD reads in GEX R2.")
    say("")

    if prot_conf == 0 and gex_conf == 0 and und_conf == 0:
        say("(b) HANDLE MISMATCH is the surviving explanation. The BD AbSeq oligo")
        say("    carries a poly(A) tail and BD's own 5' PCR handle. The 10x Cell")
        say("    Surface Protein library is amplified off the TotalSeq-C handle, so")
        say("    a BD molecule is never indexed and never sequenced. This is a")
        say("    reagent chemistry incompatibility, not a bioinformatics problem.")
        say("    Consequence: CD163 and CD169 are unrecoverable from this run.")
        say("    Myeloid Tier 2 proceeds RNA-driven on CD163, SIGLEC1, MRC1, MARCO,")
        say("    VSIG4, with CD11b/CD68/CD206/CD1c/CX3CR1 carrying the surface side.")

say("")
say("Weak evidence note: a single 20-mer anchor with no partner is at the chance")
say(f"level for a mammalian read pool and is NOT called confident (MIN_ANCHORS={MIN_ANCHORS}).")
say("Filter bd_hits_detail.csv on confident == 1 before believing anything.")
say("=" * 68)

with open(os.path.join(OUT_DIR, "verdict.txt"), "w") as fh:
    fh.write("\n".join(lines) + "\n")
print(f"\nWrote {os.path.join(OUT_DIR, 'verdict.txt')}")
