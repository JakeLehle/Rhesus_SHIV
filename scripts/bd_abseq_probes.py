#!/usr/bin/env python3
"""
bd_abseq_probes.py
==================
Step 1 of 3. Read-only. Builds the probe set and the FASTQ manifest for the
BD AbSeq (CD163 / CD169) recovery diagnostic.

Question being tested
---------------------
The two BD AbSeq tags return zero counts. Three non-exclusive explanations:

  (a) PATTERN BUG   - the tag is in the PROTEIN library but feature_reference.csv
                      encodes 36 N's as the offset, so CellRanger looks for the
                      36bp barcode starting at base 37 of R2.
  (b) HANDLE MISMATCH - BD's 5' PCR handle is not the TotalSeq-C handle the 10x
                      Cell Surface Protein library PCR primes off, so the molecule
                      is never amplified or indexed. Not in the PROTEIN FASTQ at all.
  (c) GEX CARRYOVER - a poly(A)-tailed BD oligo can be primed by the solution
                      poly(dT) RT primer and template-switch onto the bead TSO,
                      picking up a real 10x CB + UMI. That product is short and
                      should partition into the SPRI supernatant, but size
                      selection leaks. If so it sits in the GEX FASTQ, in
                      reverse-complement orientation, and is recoverable.

Readouts that separate them
---------------------------
  PROTEIN R2 hits, tight offset spike     -> (a). Fix the pattern, recount.
  PROTEIN R2 zero, controls fire          -> (b). Reagent incompatibility, banked.
  GEX R2 hits, flat offset distribution   -> (c). Recoverable via STARsolo.

Outputs (all into OUT_DIR):
  probes_bd.txt      one sequence per line, for grep -F -f
  probes_ctrl.txt    TotalSeq-C positive controls, for grep -F -f
  probe_map.csv      sequence -> name / target / orientation / kind / length
  scan_manifest.tsv  tag, library, lib_type, read, path

Run in Spyder cell by cell, or as a script. Nothing is modified on disk except
the four files above.
"""

# %% Cell 1 - Parameters
# =============================================================================
import os
import csv
from pathlib import Path

RAW_ROOT = "/master/zwallis/WORKING/SC/10X_RAW"
OUT_DIR = "/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV/adt_probe_scan"

# BD AbSeq antibody clone-specific barcodes (ABC), 36bp, from BD product pages.
# CD163 = cat. 940058 (AHS0062). CD169 / SIGLEC1 = cat. 940223.
BD_TAGS = {
    "CD163": "TATTATGTGCGAACTATGGTATCCGTATTGAGGGCT",
    "CD169": "CATTAAGCACGAAGGGTATAGGTAGGAACGGTTGGC",
}

# TotalSeq-C 15bp barcodes already in feature_reference.csv. These are the
# internal positive controls. Without them a zero on the BD probes is just an
# absence, not a result.
CTRL_TAGS = {
    "CD4": "GAGGTTAGTGATGGA",
    "CD8": "GCGCAACTTGATGAT",
    "CD68": "CGGTGTTTGTAGCAA",
    "CD28": "TGAGAACGACCCTAA",
}

# Sub-anchors tolerate sequencing error and 3' truncation of the 36-mer.
# Tiled at these starts, each ANCHOR_LEN long.
ANCHOR_LEN = 20
ANCHOR_STARTS = [0, 8, 16]

# R1 in 10x 5' v2 is 26bp (16bp CB + 10bp UMI) and structurally cannot hold a
# 36bp tag. Cell 4 of the scan script verifies this from the data rather than
# trusting it. Leave as R2 unless that check says otherwise.
SCAN_READ = "R2"

# Undetermined reads are index-unmatched but still clustered. Cheap to include.
INCLUDE_UNDETERMINED = True

# Directories inside PROTEIN_LIBRARY_B / VDJ_RH_LIBRARY_H that hold pre-repair
# duplicates. Must never enter the manifest.
EXCLUDE_DIRS = {"BACKUP", "REPAIR_BACKUP", "repair_tmp"}

os.makedirs(OUT_DIR, exist_ok=True)
print(f"RAW_ROOT : {RAW_ROOT}")
print(f"OUT_DIR  : {OUT_DIR}")


# %% Cell 2 - Probe construction
# =============================================================================
_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s):
    return s.translate(_COMP)[::-1]


probe_rows = []


def _add(seq, name, target, orientation, kind):
    probe_rows.append(
        {
            "sequence": seq,
            "probe_name": name,
            "target": target,
            "orientation": orientation,
            "kind": kind,
            "length": len(seq),
        }
    )


for target, tag in BD_TAGS.items():
    assert len(tag) == 36, f"{target} is {len(tag)}bp, expected 36"
    _add(tag, f"{target}_full_F", target, "F", "bd_full")
    _add(revcomp(tag), f"{target}_full_R", target, "R", "bd_full")
    for start in ANCHOR_STARTS:
        sub = tag[start : start + ANCHOR_LEN]
        assert len(sub) == ANCHOR_LEN
        _add(sub, f"{target}_a{start}_F", target, "F", "bd_anchor")
        _add(revcomp(sub), f"{target}_a{start}_R", target, "R", "bd_anchor")

for target, tag in CTRL_TAGS.items():
    _add(tag, f"{target}_ctrl_F", target, "F", "ctrl")
    _add(revcomp(tag), f"{target}_ctrl_R", target, "R", "ctrl")

# Sanity: no probe sequence should appear twice, and no BD anchor should collide
# with a control.
seqs = [r["sequence"] for r in probe_rows]
dupes = {s for s in seqs if seqs.count(s) > 1}
if dupes:
    raise SystemExit(f"Duplicate probe sequences, resolve before scanning: {dupes}")

probe_map = os.path.join(OUT_DIR, "probe_map.csv")
with open(probe_map, "w", newline="") as fh:
    w = csv.DictWriter(
        fh, fieldnames=["sequence", "probe_name", "target", "orientation", "kind", "length"],
        lineterminator="\n",
    )
    w.writeheader()
    w.writerows(probe_rows)

bd_seqs = [r["sequence"] for r in probe_rows if r["kind"].startswith("bd_")]
ctrl_seqs = [r["sequence"] for r in probe_rows if r["kind"] == "ctrl"]

with open(os.path.join(OUT_DIR, "probes_bd.txt"), "w") as fh:
    fh.write("\n".join(bd_seqs) + "\n")
with open(os.path.join(OUT_DIR, "probes_ctrl.txt"), "w") as fh:
    fh.write("\n".join(ctrl_seqs) + "\n")

print(f"\nBD probes    : {len(bd_seqs)}  (2 antibodies x [1 full + {len(ANCHOR_STARTS)} anchors] x 2 orientations)")
print(f"Ctrl probes  : {len(ctrl_seqs)}")
print(f"probe_map    : {probe_map}")
for r in probe_rows:
    print(f"  {r['probe_name']:>16s}  {r['kind']:<10s} {r['length']:>2d}bp  {r['sequence']}")


# %% Cell 3 - FASTQ manifest
# =============================================================================
root = Path(RAW_ROOT)
manifest = []


def _clean(paths):
    keep = []
    for p in paths:
        if EXCLUDE_DIRS & set(p.parts):
            continue
        keep.append(p)
    return sorted(keep)


# GEX libraries A-H
for lib_dir in sorted(root.glob("GEX_LIBRARY_?")):
    letter = lib_dir.name.split("_")[-1]
    hits = _clean(lib_dir.glob(f"GEX_LIBRARY_{letter}-*_{SCAN_READ}_001.fastq.gz"))
    for p in hits:
        manifest.append(
            {"tag": f"GEX_{letter}_{SCAN_READ}", "library": letter, "lib_type": "GEX",
             "read": SCAN_READ, "path": str(p)}
        )

# PROTEIN (Antibody Capture) libraries A-H
for lib_dir in sorted(root.glob("PROTEIN_LIBRARY_?")):
    letter = lib_dir.name.split("_")[-1]
    hits = _clean(lib_dir.glob(f"PROTEIN_LIBRARY_{letter}-*_{SCAN_READ}_001.fastq.gz"))
    for p in hits:
        manifest.append(
            {"tag": f"PROT_{letter}_{SCAN_READ}", "library": letter, "lib_type": "PROTEIN",
             "read": SCAN_READ, "path": str(p)}
        )

# Undetermined. Novogene naming: _1 = R1, _2 = R2.
if INCLUDE_UNDETERMINED:
    und = _clean((root / "Undetermined").glob("Undetermined_*_L8_2.fq.gz"))
    for p in und:
        manifest.append(
            {"tag": "UND_R2", "library": "-", "lib_type": "UNDETERMINED",
             "read": "R2", "path": str(p)}
        )

manifest_path = os.path.join(OUT_DIR, "scan_manifest.tsv")
with open(manifest_path, "w", newline="") as fh:
    w = csv.DictWriter(
        fh, fieldnames=["tag", "library", "lib_type", "read", "path"], delimiter="\t",
        lineterminator="\n",
    )
    w.writeheader()
    w.writerows(manifest)

n_gex = sum(1 for m in manifest if m["lib_type"] == "GEX")
n_prot = sum(1 for m in manifest if m["lib_type"] == "PROTEIN")
n_und = sum(1 for m in manifest if m["lib_type"] == "UNDETERMINED")

print(f"\nmanifest : {manifest_path}")
print(f"  GEX          : {n_gex}")
print(f"  PROTEIN      : {n_prot}")
print(f"  UNDETERMINED : {n_und}")
print(f"  TOTAL        : {len(manifest)}")

missing = [m for m in manifest if not os.path.exists(m["path"])]
if missing:
    print("\n*** WARNING: paths in manifest that do not exist:")
    for m in missing:
        print(f"    {m['path']}")

if n_gex != 8 or n_prot != 8:
    print("\n*** WARNING: expected 8 GEX and 8 PROTEIN R2 files. Check EXCLUDE_DIRS "
          "and the glob before submitting the scan.")

for m in manifest:
    print(f"  {m['tag']:<12s} {m['lib_type']:<13s} {m['path']}")

print("\nNext: sbatch bd_abseq_scan.sh")
