#!/bin/bash
#SBATCH -J shiv_starsolo_bc
#SBATCH -o /master/jlehle/WORKING/LOGS/shiv_starsolo_bc_recovery.o.log
#SBATCH -e /master/jlehle/WORKING/LOGS/shiv_starsolo_bc_recovery.e.log
#SBATCH --partition=normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=24:00:00
###############################################################################
# shiv_starsolo_bc_recovery.sh
#
# Barcode-aware SHIV recovery against the CASSETTE-CONTAINING references built
# by shiv_build_bc_reference.py. Adapted from shiv_starsolo_recovery.sh; the
# barcode/UMI logic, the strand call and the relaxed length filters are carried
# over unchanged because they are what took library H from 7 counted reads to
# ~337. Only three things differ.
#
# 1. THE REFERENCE. The old index was MN816822.1 with NO cassette, so a read
#    spanning the barcode could not align contiguously and the ~947-read SHIV
#    pool is specifically depleted at exactly the locus we care about.
#
# 2. THE MISMATCH BUDGET, 0.1 -> 0.2. Measured 2026-10-27: STAR does NOT charge
#    a mismatch for reference N, so the whole budget is available for real
#    divergence. At 150bp that is 30 mismatches, enough for roughly 30% G-to-A
#    hypermutation on top of the cassette. This is the point of the exercise:
#    an APOBEC3-edited provirus breaks an exact string match outright but an
#    aligner just counts mismatches.
#
# 3. SEVERAL REFERENCES, ONE PASS EACH. The smoke test could not separate the
#    variants: all extracted 100% correct cores, and they differed only by 27
#    to 30 uniquely-mapped out of 30 synthetic reads, in STAR's "other" bucket
#    which we do not understand. That is below the resolution to pick on
#    principle, so run them and let real yield decide.
#
# V0 (the OLD cassette-less index) is included deliberately, at the SAME 0.2
# budget as the rest. Without it, a gain could be the reference OR the relaxed
# filter and there would be no way to tell them apart. V0-at-0.2 isolates the
# filter; everything else against V0-at-0.2 isolates the reference.
#
# Negatives A/B/C/G run at full depth as the false-positive floor. Do not trim
# them to save time; a count with no floor is not a count.
#
#   sbatch shiv_starsolo_bc_recovery.sh
#
# Author: Jake Lehle
# Date: October 2026
###############################################################################

set -euo pipefail
source ~/anaconda3/bin/activate
eval "$(conda shell.bash hook)"
conda activate sc_pre
export PATH="${CONDA_PREFIX}/bin:${PATH}"   # conda samtools; the system one is broken

# ─── Configuration ───────────────────────────────────────────────────────────
THREADS=${SLURM_CPUS_PER_TASK:-16}
WORKDIR="/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV"
FASTQ_BASE="/master/zwallis/WORKING/SC/10X_RAW"
REFDIR="${WORKDIR}/shiv_bc_reference"
OLD_INDEX="${WORKDIR}/shiv_starsolo/index"      # V0: the cassette-less index
OUTDIR="${WORKDIR}/shiv_starsolo_bc"

STAR_ENV_BIN="${HOME}/anaconda3/envs/star/bin/STAR"
STAR_BIN="$(command -v STAR || true)"
[[ -z "${STAR_BIN}" || ! -x "${STAR_BIN}" ]] && STAR_BIN="${STAR_ENV_BIN}"

# --- Confirmed barcode config (check_sample_whitelist.py; do not change) -----
CB_LEN=16
UMI_LEN=12
WHITELIST="${WORKDIR}/shiv_starsolo/ref/3M-5pgex-jan-2023.txt"

# --- Strand: REVERSE for 10x 5'. R2 is ~89% minus-strand while the collapsed
#     SHIV gene is on +. Forward counted only the wrong-strand minority. ------
STRAND="Reverse"

# --- Relaxed filters, carried over, except the mismatch budget ---------------
SCORE_MIN_OVER_L=0
MATCH_MIN_OVER_L=0
MATCH_MIN_ABS=30
MISMATCH_OVER_L=0.2          # was 0.1; reference N is free, so spend it on divergence

# --- References to run. Comment out to narrow once real yield has spoken. ----
#     Name -> index dir. V0 must stay for the comparison to mean anything.
declare -A REFS=(
    [V0_nocassette]="${OLD_INDEX}"
    [V2allN]="${REFDIR}/idx_SHIVAD8EO_V2allN"
    [V3r1]="${REFDIR}/idx_SHIVAD8EO_V3r1"
    [V3r2]="${REFDIR}/idx_SHIVAD8EO_V3r2"
)
REF_ORDER="V0_nocassette V2allN V3r1 V3r2"

LIBS="A B C D E F G H"
declare -A LIB_NAME=(
    [A]="SHIV_Pre_PBMC"             [B]="SHIV_Pre_LN"
    [C]="SHIV_Pre_PBMC_IndianMixed" [D]="SHIV_21DPI_Pre_LN_Mixed"
    [E]="SHIV_21DPI_PBMC"           [F]="SHIV_21DPI_NonInf_LN_Mixed"
    [G]="SHIV_Necropsy_PBMC"        [H]="SHIV_Necropsy_LN"
)
# Historical baselines at the OLD reference and the OLD 0.1 budget.
declare -A MM2_BASELINE=( [H]=553 [D]=35 [E]=48 )

# SMOKE TEST FIRST: set LIBS="E H" and REF_ORDER="V0_nocassette V2allN" for one
# pass to confirm the plumbing, then restore. Do not discover a typo at hour six.

mkdir -p "${OUTDIR}"/{solo,reports}
cd "${OUTDIR}"
LOGFILE="${OUTDIR}/reports/starsolo_bc_recovery.log"
exec > >(tee -a "${LOGFILE}") 2>&1

echo "============================================================"
echo "SHIV STARsolo BC Recovery — $(date)"
echo "STAR            : ${STAR_BIN}"
if [[ ! -x "${STAR_BIN}" ]]; then
    echo "ERROR: no STAR at ${STAR_BIN}"
    echo "  install: conda create -y -n star -c bioconda -c conda-forge star=2.7.11b"
    exit 1
fi
echo "STAR version    : $(${STAR_BIN} --version)"
echo "mismatch budget : ${MISMATCH_OVER_L} (was 0.1 on the old reference)"
echo "strand          : ${STRAND}"
echo "references      : ${REF_ORDER}"
echo "libraries       : ${LIBS}"
echo "============================================================"

SUMMARY="${OUTDIR}/reports/bc_recovery_summary.tsv"
printf 'reference\tlib\tsample\tgroup\taligned\tmapq30\tbreadth_pct\tcassette_span\n' > "${SUMMARY}"

for REF in ${REF_ORDER}; do
    INDEX="${REFS[$REF]}"
    if [[ ! -d "${INDEX}" ]]; then
        echo "[WARN] index not found for ${REF}: ${INDEX} — skipping reference"
        continue
    fi
    # Cassette window, 0-based half-open, in THIS reference's coordinates.
    # V0 has no cassette; the window is where it WOULD be, so the count there is
    # by construction reads that merely reach the locus.
    CASSETTE_START=5862
    CASSETTE_END=5896
    CONTIG="SHIVAD8EO"
    [[ "${REF}" != "V0_nocassette" ]] && CONTIG="SHIVAD8EO_${REF}"

    echo ""
    echo "############################################################"
    echo "# REFERENCE ${REF}   index ${INDEX}"
    echo "############################################################"

    for L in ${LIBS}; do
        NAME="${LIB_NAME[$L]}"
        GROUP="neg"; [[ "DEFH" == *"${L}"* ]] && GROUP="INF"
        FQDIR="${FASTQ_BASE}/GEX_LIBRARY_${L}"
        PREFIX="${OUTDIR}/solo/${REF}/${L}_${NAME}_"
        mkdir -p "$(dirname "${PREFIX}")"

        echo ""
        echo "--- ${REF} / library ${L}: ${NAME} [${GROUP}] ---"
        if [[ ! -d "${FQDIR}" ]]; then
            echo "  [WARN] FASTQ dir not found, skipping: ${FQDIR}"
            continue
        fi
        R2S=$(ls "${FQDIR}"/*_R2_*.fastq.gz 2>/dev/null | sort | paste -sd, -)
        R1S=$(ls "${FQDIR}"/*_R1_*.fastq.gz 2>/dev/null | sort | paste -sd, -)
        if [[ -z "${R2S}" || -z "${R1S}" ]]; then
            echo "  [WARN] missing R1/R2 in ${FQDIR}, skipping"
            continue
        fi

        "${STAR_BIN}" \
            --runMode alignReads \
            --genomeDir "${INDEX}" \
            --readFilesIn "${R2S}" "${R1S}" \
            --readFilesCommand zcat \
            --runThreadN "${THREADS}" \
            --outFileNamePrefix "${PREFIX}" \
            --soloType CB_UMI_Simple \
            --soloCBwhitelist "${WHITELIST}" \
            --soloCBstart 1 --soloCBlen ${CB_LEN} \
            --soloUMIstart $((CB_LEN + 1)) --soloUMIlen ${UMI_LEN} \
            --soloBarcodeReadLength 0 \
            --soloStrand ${STRAND} \
            --soloFeatures Gene \
            --soloCellFilter None \
            --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
            --soloUMIdedup 1MM_CR \
            --soloUMIfiltering MultiGeneUMI_CR \
            --clipAdapterType CellRanger4 \
            --outFilterScoreMinOverLread ${SCORE_MIN_OVER_L} \
            --outFilterMatchNminOverLread ${MATCH_MIN_OVER_L} \
            --outFilterMatchNmin ${MATCH_MIN_ABS} \
            --outFilterMismatchNoverLmax ${MISMATCH_OVER_L} \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattributes NH HI AS nM CB UB \
            --limitBAMsortRAM 60000000000

        BAM="${PREFIX}Aligned.sortedByCoord.out.bam"
        samtools index "${BAM}"

        ALN=$(samtools view -c -F 0x4 "${BAM}")
        HQ=$(samtools view -c -q 30 "${BAM}")
        BREADTH=$(samtools depth -a "${BAM}" \
                  | awk '{n++; if($3>0)c++} END{printf "%.1f", (n>0?100*c/n:0)}')
        # Reads whose alignment covers the WHOLE cassette window. This is the
        # number the whole exercise is about, and it is 0 by construction on V0.
        # A read counts only if ONE aligned block covers the whole cassette.
        # Spanning end-to-end is not enough: a CIGAR like 58M34N58M reaches past
        # the window while skipping every cassette base, so it carries no
        # barcode. Requiring a single M/=/X block that contains [s,e) is the
        # same test shiv_check_bc_reference_alignment.py applies.
        SPAN=$(samtools view "${BAM}" "${CONTIG}:$((CASSETTE_START+1))-${CASSETTE_END}" 2>/dev/null \
               | awk -v s=${CASSETTE_START} -v e=${CASSETTE_END} '
                   { r=$4-1; c=$6; ok=0;
                     while (match(c, /^[0-9]+[MIDNSHP=X]/)) {
                       n=substr(c,1,RLENGTH-1)+0; op=substr(c,RLENGTH,1);
                       if (op ~ /[M=X]/) { if (r<=s && r+n>=e) ok=1; r+=n }
                       else if (op ~ /[DN]/) r+=n;
                       c=substr(c,RLENGTH+1) }
                     if (ok) k++ }
                   END{print k+0}')

        echo "  aligned ${ALN}, MAPQ30 ${HQ}, breadth ${BREADTH}%, cassette-spanning ${SPAN}"
        if [[ -n "${MM2_BASELINE[$L]:-}" ]]; then
            echo "  (minimap2 baseline on the old reference: ${MM2_BASELINE[$L]})"
        fi
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "${REF}" "${L}" "${NAME}" "${GROUP}" "${ALN}" "${HQ}" "${BREADTH}" "${SPAN}" \
            >> "${SUMMARY}"
    done
done

echo ""
echo "============================================================"
echo "DONE — $(date)"
echo "summary: ${SUMMARY}"
echo ""
echo "Read the cassette_span column FIRST, by reference:"
echo "  V0_nocassette should be ~0. If it is not, the window is being reached"
echo "  by reads that clip rather than span, and the comparison is not clean."
echo "  Any cassette-containing reference should exceed it. The size of that"
echo "  gap, in the INFECTED libraries and NOT in the negatives, is the"
echo "  measurement this whole reference exercise was for."
echo ""
echo "THEN run shiv_bc_core_extract.py against the winning reference's BAMs to"
echo "pull the 34bp cassette by coordinate, read the register off the read, and"
echo "join CB to the object."
echo "============================================================"
