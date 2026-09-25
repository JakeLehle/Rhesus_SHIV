#!/bin/bash
#SBATCH --job-name=shiv_remap
#SBATCH --output=/master/jlehle/WORKING/LOGS/shiv_remap_%j.out
#SBATCH --error=/master/jlehle/WORKING/LOGS/shiv_remap_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=6:00:00
#SBATCH --partition=normal

###############################################################################
# shiv_targeted_remap.sh
#
# Context: Phase A from the previous script showed the SHIVAD8EO reference
#          is NOT divergent — 864 mapped reads cover 97.4% of the genome
#          with nM:i:0 (zero mismatches). The low SHIV+ cell count (46/126K)
#          is likely biological (low viral transcription), not technical.
#
# But: We should still check whether there are additional SHIV reads hiding
#      in the 280M unmapped reads that cellranger's STAR aligner missed.
#
# Strategy (fast and targeted):
#   1. Re-use the unmapped FASTQs from the previous run (already extracted)
#   2. Remap against SHIV-only with minimap2 (permissive params)
#   3. Compare: how many additional reads map vs the 864 cellranger found?
#   4. If substantial: characterize them (which genes, what quality)
#   5. If negligible: confirm the detection rate is biological, not technical
#
# Also runs on Libraries D and E (21 DPI) for comparison.
#
# Author: Jake Lehle
# Date: June 2, 2026
###############################################################################

set -euo pipefail
source ~/anaconda3/bin/activate
eval "$(conda shell.bash hook)"
conda activate sc_pre
export PATH="${CONDA_PREFIX}/bin:${PATH}"

# ─── Configuration ──────────────────────────────────────────────────────────
THREADS=${SLURM_CPUS_PER_TASK:-16}
WORKDIR="/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV"
ANALYSIS_DIR="${WORKDIR}/Analysis"
OUTDIR="${WORKDIR}/shiv_denovo_analysis"
SHIV_REF="${OUTDIR}/SHIVAD8EO_reference.fasta"
SHIV_CONTIG="SHIVAD8EO"

# Libraries to process
declare -A LIBRARIES=(
    ["H"]="SHIV_Necropsy_LN"
    ["D"]="SHIV_21DPI_Pre_LN_Mixed"
    ["E"]="SHIV_21DPI_PBMC"
)

mkdir -p "${OUTDIR}"/{reads,remap,reports}
cd "${OUTDIR}"

LOGFILE="${OUTDIR}/reports/targeted_remap.log"
exec > >(tee -a "${LOGFILE}") 2>&1

echo "============================================================"
echo "SHIV Targeted Remap — $(date)"
echo "============================================================"

# ─── Verify SHIV reference exists ──────────────────────────────────────────
if [[ ! -f "${SHIV_REF}" ]]; then
    echo "Extracting SHIV reference..."
    CR_FASTA="/master/jlehle/cellranger-10.0.0/Mmul10_SHIVAD8EO_v10/fasta/genome.fa"
    samtools faidx "${CR_FASTA}" "${SHIV_CONTIG}" > "${SHIV_REF}"
fi
echo "SHIV ref: ${SHIV_REF} ($(grep -v '^>' "${SHIV_REF}" | tr -d '\n' | wc -c) bp)"

# ─── Build minimap2 index ──────────────────────────────────────────────────
MMI="${OUTDIR}/remap/shiv_only.mmi"
echo "Building minimap2 index..."
minimap2 -d "${MMI}" "${SHIV_REF}" 2>/dev/null
echo ""

# ─── Process each library ──────────────────────────────────────────────────
for LIB_KEY in H D E; do
    LIB_NAME="${LIBRARIES[${LIB_KEY}]}"
    LIB_DIR="${ANALYSIS_DIR}/${LIB_NAME}/outs"
    SAMPLE_BAM="${LIB_DIR}/per_sample_outs/${LIB_NAME}/sample_alignments.bam"
    UNASSIGNED_BAM="${LIB_DIR}/unassigned_alignments.bam"

    echo "============================================================"
    echo "Library ${LIB_KEY}: ${LIB_NAME}"
    echo "============================================================"

    if [[ ! -f "${SAMPLE_BAM}" ]]; then
        echo "  WARNING: BAM not found, skipping: ${SAMPLE_BAM}"
        continue
    fi

    # ── Step 1: How many reads did cellranger map to SHIV? ──────────────
    CR_MAPPED=$(samtools view -c "${SAMPLE_BAM}" "${SHIV_CONTIG}")
    CR_UNMAPPED=$(samtools view -c -f 4 "${SAMPLE_BAM}")
    CR_TOTAL=$(samtools view -c "${SAMPLE_BAM}")
    echo "  CellRanger alignment stats:"
    echo "    Total reads in BAM:  ${CR_TOTAL}"
    echo "    Mapped to SHIV:      ${CR_MAPPED}"
    echo "    Unmapped:            ${CR_UNMAPPED}"

    # ── Step 2: Extract unmapped reads ──────────────────────────────────
    UNMAPPED_FQ="${OUTDIR}/reads/unmapped_${LIB_KEY}.fastq"

    # Re-use from previous run if Library H already extracted
    if [[ "${LIB_KEY}" == "H" && -f "${OUTDIR}/reads/all_unmapped.fastq" ]]; then
        echo "  Re-using previously extracted unmapped reads for Library H"
        UNMAPPED_FQ="${OUTDIR}/reads/all_unmapped.fastq"
    else
        echo "  Extracting unmapped reads..."
        samtools view -@ "${THREADS}" -f 4 -F 0x900 "${SAMPLE_BAM}" | \
            awk '{print "@"$1"\n"$10"\n+\n"$11}' \
            > "${UNMAPPED_FQ}"

        # Also grab from unassigned BAM
        if [[ -f "${UNASSIGNED_BAM}" ]]; then
            samtools view -@ "${THREADS}" -f 4 -F 0x900 "${UNASSIGNED_BAM}" | \
                awk '{print "@"$1"\n"$10"\n+\n"$11}' \
                >> "${UNMAPPED_FQ}"
        fi
    fi

    N_UNMAPPED=$(( $(wc -l < "${UNMAPPED_FQ}") / 4 ))
    echo "  Unmapped reads: ${N_UNMAPPED}"

    # ── Step 3: Permissive remap against SHIV-only ──────────────────────
    echo ""
    echo "  Remapping unmapped reads against SHIV-only (minimap2, permissive)..."
    echo "  (This handles ${N_UNMAPPED} reads — may take 10-30 min)"

    REMAP_BAM="${OUTDIR}/remap/recovered_${LIB_KEY}.bam"
    REMAP_LOG="${OUTDIR}/remap/minimap2_${LIB_KEY}.log"
    REMAP_UNSORTED="${OUTDIR}/remap/recovered_${LIB_KEY}.unsorted.bam"

    # Step 3a: minimap2 → samtools view (filter mapped only) → unsorted BAM
    # Keep minimap2 stderr in a log file so we can see any errors
    set +eo pipefail  # Temporarily relax — pipe may exit non-zero if 0 reads map
    minimap2 -a -x sr \
        -k 10 -w 5 \
        -A 1 -B 2 -O 2,16 -E 1,0 \
        --score-N 0 \
        --secondary=no \
        -t "${THREADS}" \
        "${SHIV_REF}" \
        "${UNMAPPED_FQ}" \
        2>"${REMAP_LOG}" | \
        samtools view -@ 4 -b -F 4 -o "${REMAP_UNSORTED}" -
    PIPE_EXIT=$?
    set -eo pipefail  # Re-enable

    echo "  minimap2 log: ${REMAP_LOG}"
    tail -3 "${REMAP_LOG}" | sed 's/^/    /'

    # Step 3b: Sort and index if we got output
    if [[ -s "${REMAP_UNSORTED}" ]]; then
        samtools sort -@ "${THREADS}" -o "${REMAP_BAM}" "${REMAP_UNSORTED}"
        samtools index "${REMAP_BAM}"
        rm -f "${REMAP_UNSORTED}"
        RECOVERED=$(samtools view -c "${REMAP_BAM}")
    else
        RECOVERED=0
        echo "  WARNING: No output from minimap2 (pipe exit: ${PIPE_EXIT})"
        echo "  Check minimap2 log: ${REMAP_LOG}"
    fi
    echo "  Recovered SHIV-like reads: ${RECOVERED}"

    # ── Step 4: Quality assessment of recovered reads ───────────────────
    if [[ "${RECOVERED}" -gt 0 ]]; then
        echo ""
        echo "  MAPQ distribution of recovered reads:"
        samtools view "${REMAP_BAM}" | \
            awk '{print $5}' | sort -n | uniq -c | sort -rn | head -10 | \
            awk '{printf "    MAPQ %3d: %d reads\n", $2, $1}'

        # Alignment scores
        echo ""
        echo "  Alignment score (AS tag) distribution:"
        AS_DATA=$(samtools view "${REMAP_BAM}" | \
            grep -oP 'AS:i:\K-?[0-9]+' || true)
        if [[ -n "${AS_DATA}" ]]; then
            echo "${AS_DATA}" | sort -n | \
                awk 'BEGIN{n=0; sum=0}
                     {vals[n++]=$1; sum+=$1}
                     END{
                       printf "    n=%d, mean=%.1f, median=%d, min=%d, max=%d\n",
                         n, sum/n, vals[int(n/2)], vals[0], vals[n-1]
                     }'
        fi

        # Coverage on SHIV genome from recovered reads
        echo ""
        echo "  Coverage from recovered reads:"
        samtools depth -a "${REMAP_BAM}" \
            > "${OUTDIR}/reports/recovered_coverage_${LIB_KEY}.txt"

        awk '
        BEGIN { sum=0; n=0; zero=0 }
        {
            sum += $3; n++
            if ($3 == 0) zero++
        }
        END {
            printf "    Positions: %d\n", n
            printf "    Mean depth: %.2f\n", (n>0 ? sum/n : 0)
            printf "    Breadth:    %.1f%%\n", (n>0 ? 100.0*(n-zero)/n : 0)
        }' "${OUTDIR}/reports/recovered_coverage_${LIB_KEY}.txt"

        # How many are high-quality vs noise?
        HQ=$(samtools view -c -q 30 "${REMAP_BAM}")
        LQ=$(( RECOVERED - HQ ))
        echo ""
        echo "  Quality breakdown:"
        echo "    MAPQ >= 30 (high confidence): ${HQ}"
        echo "    MAPQ < 30 (low confidence):   ${LQ}"
    fi

    # ── Step 5: Compare ────────────────────────────────────────────────
    echo ""
    echo "  ── COMPARISON ──"
    echo "  CellRanger found:   ${CR_MAPPED} reads on SHIV"
    echo "  Minimap2 recovered: ${RECOVERED} additional reads from unmapped pool"
    if [[ "${RECOVERED}" -gt 0 ]]; then
        RATIO=$(echo "scale=1; ${RECOVERED} * 100 / ${CR_MAPPED}" | bc 2>/dev/null || echo "N/A")
        echo "  Recovery ratio:     ${RATIO}% additional vs cellranger"
    fi
    echo ""

done

###############################################################################
# Summary
###############################################################################
echo ""
echo "============================================================"
echo "OVERALL SUMMARY"
echo "============================================================"
echo ""
echo "Library | CR Mapped | Recovered | Total SHIV reads"
echo "--------|-----------|-----------|------------------"

for LIB_KEY in H D E; do
    LIB_NAME="${LIBRARIES[${LIB_KEY}]}"
    SAMPLE_BAM="${ANALYSIS_DIR}/${LIB_NAME}/outs/per_sample_outs/${LIB_NAME}/sample_alignments.bam"
    REMAP_BAM="${OUTDIR}/remap/recovered_${LIB_KEY}.bam"

    if [[ -f "${SAMPLE_BAM}" ]]; then
        CR=$(samtools view -c "${SAMPLE_BAM}" "${SHIV_CONTIG}" 2>/dev/null || echo "0")
    else
        CR="N/A"
    fi

    if [[ -f "${REMAP_BAM}" ]]; then
        REC=$(samtools view -c "${REMAP_BAM}" 2>/dev/null || echo "0")
    else
        REC="0"
    fi

    TOTAL="N/A"
    if [[ "${CR}" != "N/A" && "${REC}" != "0" ]]; then
        TOTAL=$(( CR + REC ))
    elif [[ "${CR}" != "N/A" ]]; then
        TOTAL="${CR}"
    fi

    printf "   %s     | %9s | %9s | %s\n" "${LIB_KEY}" "${CR}" "${REC}" "${TOTAL}"
done

echo ""
echo "INTERPRETATION:"
echo "  If recovery is <10% of CellRanger count:"
echo "    → Reference is adequate. Low SHIV detection is biological."
echo "    → Proceed with Scanpy analysis using current counts."
echo ""
echo "  If recovery is >50% of CellRanger count:"
echo "    → There ARE hidden SHIV reads that CellRanger missed."
echo "    → Consider rebuilding reference or post-hoc counting."
echo ""
echo "  If recovery is >200% of CellRanger count:"
echo "    → Likely many false positives from permissive mapping."
echo "    → Filter by MAPQ and alignment score before trusting."
echo ""
echo "KEY FILES:"
echo "  ${OUTDIR}/reports/targeted_remap.log"
echo "  ${OUTDIR}/remap/recovered_H.bam  (Library H recovered reads)"
echo "  ${OUTDIR}/remap/recovered_D.bam  (Library D recovered reads)"
echo "  ${OUTDIR}/remap/recovered_E.bam  (Library E recovered reads)"
echo ""
echo "Done: $(date)"
