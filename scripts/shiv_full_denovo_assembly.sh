#!/bin/bash
#SBATCH --job-name=shiv_megahit
#SBATCH --output=/master/jlehle/WORKING/LOGS/shiv_megahit_%j.out
#SBATCH --error=/master/jlehle/WORKING/LOGS/shiv_megahit_%j.err
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=96
#SBATCH --mem=990G
#SBATCH --time=7-00:00:00
#SBATCH --partition=normal

###############################################################################
# shiv_megahit_assembly.sh
#
# Purpose: De novo assemble ALL ~280M unmapped reads from Library H with
#          MEGAHIT, reference-free, to recover SHIV contigs INCLUDING any
#          in vivo divergent quasispecies that escaped CellRanger/minimap2.
#
# Scope:   ASSEMBLY ONLY. No BLAST, no recruitment pre-filter, no subsampling.
#          BLAST identification against a comprehensive DB is a SEPARATE script
#          designed after inspecting the contig landscape produced here.
#
# Rationale for key choices:
#   - MEGAHIT: tolerates the extreme uneven coverage of scRNA-derived reads
#     (host transcripts at huge depth, SHIV at ~9x or lower) and scales to
#     280M reads at low memory. Single-end native (-r).
#   - --min-count 1: do NOT prune low-multiplicity k-mers. A divergent
#     low-frequency variant could sit at 2-5x coverage; pruning would erase
#     exactly the signal we are hunting. Over-assemble now, filter later.
#   - --presets meta-sensitive: sensitive k-mer schedule tuned for recovering
#     low-abundance / divergent strains in metagenomes. Overrides manual
#     k-sweep (do NOT also pass --k-min/--k-max/--k-step).
#   - Reference-free: assembly does not consult SHIVAD8EO, so divergent SHIV
#     assembles regardless of similarity to the canonical reference. This is
#     the whole point — it sidesteps the recruitment blind spot.
#
# Resources: 96 cores, 990GB SLURM alloc (MEGAHIT reads physical RAM for
#            --memory 0.9 -> ~900GB on this 1TB node).
#
# Author: Jake Lehle
# Date:   June 14, 2026
###############################################################################

set -euo pipefail
source ~/anaconda3/bin/activate
eval "$(conda shell.bash hook)"
conda activate sc_pre
export PATH="${CONDA_PREFIX}/bin:${PATH}"

# --- Configuration ----------------------------------------------------------
THREADS=${SLURM_CPUS_PER_TASK:-96}
WORKDIR="/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV"
OUTDIR="${WORKDIR}/shiv_megahit_assembly"

# Source: all unmapped reads from Library H (already extracted, single-end cDNA)
UNMAPPED_FQ="${WORKDIR}/shiv_denovo_analysis/reads/all_unmapped.fastq"

# MEGAHIT output dir MUST NOT pre-exist (MEGAHIT refuses to overwrite).
MEGAHIT_OUT="${OUTDIR}/megahit_out"

mkdir -p "${OUTDIR}/reports"

LOGFILE="${OUTDIR}/reports/megahit_assembly.log"
exec > >(tee -a "${LOGFILE}") 2>&1

echo "============================================================"
echo "SHIV MEGAHIT Assembly (assembly-only) - $(date)"
echo "Threads: ${THREADS}   Node RAM target: --memory 0.9"
echo "============================================================"
echo ""

# --- Dependency check -------------------------------------------------------
echo "[CHECK] Dependencies..."
if ! command -v megahit &>/dev/null; then
    echo "  ERROR: megahit not found on PATH inside conda env sc_pre."
    echo "  Expected at: ${CONDA_PREFIX}/bin/megahit"
    echo "  Install: conda install -c bioconda megahit"
    exit 1
fi
echo "  megahit:  $(megahit --version 2>&1 || true)"
echo "  PATH:     ${CONDA_PREFIX}/bin (front of PATH)"
echo ""

# --- Verify input -----------------------------------------------------------
if [[ ! -f "${UNMAPPED_FQ}" ]]; then
    echo "ERROR: Unmapped FASTQ not found: ${UNMAPPED_FQ}"
    exit 1
fi
if [[ ! -s "${UNMAPPED_FQ}" ]]; then
    echo "ERROR: Unmapped FASTQ is empty: ${UNMAPPED_FQ}"
    exit 1
fi

echo "Input FASTQ:  ${UNMAPPED_FQ}"
echo "Output dir:   ${OUTDIR}"
echo "MEGAHIT out:  ${MEGAHIT_OUT}"
echo ""

###############################################################################
# PHASE 0: Pre-flight diagnostics (WARN AND CONTINUE — never hard-stop)
###############################################################################
echo "============================================================"
echo "PHASE 0: Pre-flight diagnostics (informational; non-blocking)"
echo "============================================================"

# --- Total read count -------------------------------------------------------
echo "  Counting reads (line count / 4)..."
TOTAL_LINES=$(wc -l < "${UNMAPPED_FQ}")
TOTAL_READS=$(( TOTAL_LINES / 4 ))
echo "  Total reads: ${TOTAL_READS}"
if (( TOTAL_LINES % 4 != 0 )); then
    echo "  WARNING: line count (${TOTAL_LINES}) is not a multiple of 4."
    echo "           FASTQ may be truncated or malformed. Continuing anyway."
fi
echo ""

# --- Read-length histogram on first 100k reads ------------------------------
# Confirms clean cDNA R2 (expect ~122bp and ~150bp). Presence of a 26bp spike
# would indicate R1 barcode/UMI contamination -> would inflate assembly with
# barcode junk. WARN ONLY.
echo "  Read-length histogram (first 100,000 reads):"
awk 'NR%4==2 { print length($0) }' "${UNMAPPED_FQ}" \
    | head -100000 \
    | sort -n \
    | uniq -c \
    | sort -rn \
    | head -20 \
    | awk '{ printf "    len=%-5s  count=%s\n", $2, $1 }' \
    || echo "    (length histogram step returned nothing — continuing)"

# Explicit check for a 26bp barcode spike in the sampled reads.
BC_SPIKE=$(awk 'NR%4==2 { if (length($0)==26) c++ } NR>=400000 { print c+0; exit } END { print c+0 }' "${UNMAPPED_FQ}" 2>/dev/null | head -1 || echo "0")
BC_SPIKE=${BC_SPIKE:-0}
if [[ "${BC_SPIKE}" =~ ^[0-9]+$ ]] && (( BC_SPIKE > 0 )); then
    echo "  WARNING: found ${BC_SPIKE} reads of length 26 in the sample."
    echo "           Possible R1 (barcode+UMI) contamination. Assembly will"
    echo "           proceed, but inspect contigs for barcode artifacts."
else
    echo "  No 26bp barcode-length spike detected in sample (good — looks like cDNA)."
fi
echo ""

# --- Dedup status note ------------------------------------------------------
echo "  NOTE: This pipeline does NOT deduplicate reads. scRNA-seq carries"
echo "        heavy PCR amplification, so apparent contig coverage is NOT a"
echo "        clean proxy for true viral abundance. Interpret coverage of any"
echo "        divergent contig with that caveat (assembly itself is unaffected)."
echo ""

###############################################################################
# PHASE 1: De novo assembly — ALL reads, MEGAHIT
###############################################################################
echo "============================================================"
echo "PHASE 1: MEGAHIT de novo assembly of ${TOTAL_READS} reads"
echo "============================================================"
echo "  Start: $(date)"
echo "  Settings: --presets meta-sensitive  --min-count 1"
echo "            --min-contig-len 150  --memory 0.9  -t ${THREADS}"
echo ""

# MEGAHIT refuses to write into an existing output dir; clear it first.
rm -rf "${MEGAHIT_OUT}"

megahit \
    -r "${UNMAPPED_FQ}" \
    -o "${MEGAHIT_OUT}" \
    -t "${THREADS}" \
    --presets meta-sensitive \
    --min-count 1 \
    --min-contig-len 150 \
    --memory 0.9 \
    2>&1 | tee "${OUTDIR}/reports/megahit_run.log"

echo ""
echo "  Assembly finished: $(date)"

ASSEMBLY="${MEGAHIT_OUT}/final.contigs.fa"

if [[ ! -f "${ASSEMBLY}" || ! -s "${ASSEMBLY}" ]]; then
    echo "  ERROR: Assembly produced no output (missing/empty):"
    echo "    ${ASSEMBLY}"
    echo "  Check logs:"
    echo "    ${OUTDIR}/reports/megahit_run.log"
    echo "    ${MEGAHIT_OUT}/log"
    exit 1
fi

echo "  Contigs file: ${ASSEMBLY}"
echo ""

###############################################################################
# PHASE 2: Assembly statistics (no reference; pure descriptive)
###############################################################################
echo "============================================================"
echo "PHASE 2: Assembly statistics"
echo "============================================================"

# Guarded contig count (zero-safe under pipefail).
NCONTIGS=$(grep -c "^>" "${ASSEMBLY}" || true)
NCONTIGS=${NCONTIGS:-0}
echo "  Total contigs: ${NCONTIGS}"
echo ""

if (( NCONTIGS == 0 )); then
    echo "  WARNING: zero contigs parsed. Skipping length statistics."
else
    awk '
        /^>/ {
            if (len > 0) {
                n++; sum += len; lens[n] = len
                if (len > max) max = len
                if (min == 0 || len < min) min = len
                if (len >= 500)   n500++
                if (len >= 1000)  n1k++
                if (len >= 5000)  n5k++
                if (len >= 10000) n10k++
            }
            len = 0; next
        }
        { len += length($0) }
        END {
            if (len > 0) {
                n++; sum += len; lens[n] = len
                if (len > max) max = len
                if (min == 0 || len < min) min = len
                if (len >= 500)   n500++
                if (len >= 1000)  n1k++
                if (len >= 5000)  n5k++
                if (len >= 10000) n10k++
            }
            if (n == 0) { print "    (no contigs)"; exit }
            # N50
            asort(lens); cumul = 0; n50 = 0
            for (i = n; i >= 1; i--) {
                cumul += lens[i]
                if (cumul >= sum/2) { n50 = lens[i]; break }
            }
            printf "    Total contigs:     %d\n", n
            printf "    Total length:      %d bp\n", sum
            printf "    Longest contig:    %d bp\n", max
            printf "    Shortest contig:   %d bp\n", min
            printf "    Mean length:       %d bp\n", (n > 0 ? sum/n : 0)
            printf "    N50:               %d bp\n", n50
            printf "    Contigs >= 500bp:  %d\n", n500 + 0
            printf "    Contigs >= 1kb:    %d\n", n1k + 0
            printf "    Contigs >= 5kb:    %d\n", n5k + 0
            printf "    Contigs >= 10kb:   %d\n", n10k + 0
        }' "${ASSEMBLY}"
fi
echo ""

###############################################################################
# Summary
###############################################################################
echo "============================================================"
echo "FINAL SUMMARY"
echo "============================================================"
echo ""
echo "Input reads:    ${TOTAL_READS} (all unmapped, Library H, no subsampling)"
echo "Assembler:      MEGAHIT v1.2.9 (meta-sensitive, min-count 1)"
echo "Total contigs:  ${NCONTIGS}"
echo ""
echo "KEY OUTPUT:"
echo "  Contigs:  ${ASSEMBLY}"
echo "  Run log:  ${OUTDIR}/reports/megahit_run.log"
echo "  Full log: ${LOGFILE}"
echo ""
echo "NEXT STEP (separate script, after inspecting contigs):"
echo "  Design BLAST identification against a comprehensive database to"
echo "  classify contigs (SHIV / SIV-HIV lentivirus / host / contaminant)."
echo "  Do NOT pre-filter against the canonical SHIVAD8EO reference."
echo ""
echo "Done: $(date)"
