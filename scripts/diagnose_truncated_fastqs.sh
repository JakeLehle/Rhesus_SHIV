#!/bin/bash
###############################################################################
# diagnose_fastqs_v2.sh
#
# Fast diagnostic for truncated 10x FASTQ files.
# Counts reads in all files for affected libraries and reports mismatches.
#
# NOTE ON 10X FASTQ ORDERING:
#   10x Chromium FASTQs (from bcl2fastq/mkfastq) are generated in LOCKSTEP.
#   Read #1 in R1 == Read #1 in R2 == Read #1 in I1 == Read #1 in I2.
#   They are NOT randomly shuffled. When a file is truncated, the missing
#   reads are ALL at the tail end. Truncating the longer files to match
#   the shortest file IS the correct paired-end repair — every kept read
#   has its mate in all four files.
#
# Usage: bash diagnose_fastqs_v2.sh
###############################################################################

# Don't exit on error — we need to handle corrupted files gracefully
set +e

FASTQ_BASE="/master/zwallis/WORKING/SC/10X_RAW"

echo "============================================================"
echo " FASTQ Truncation Diagnostic v2"
echo " $(date)"
echo "============================================================"
echo ""

count_reads() {
    # Count reads in a FASTQ.gz file
    # Uses unpigz if available (much faster), falls back to zcat
    local FQ="$1"

    if command -v unpigz &>/dev/null; then
        LINES=$(unpigz -c "${FQ}" 2>/dev/null | wc -l)
    else
        LINES=$(zcat "${FQ}" 2>/dev/null | wc -l)
    fi

    local EXIT_CODE=$?
    local READS=$((LINES / 4))
    local REMAINDER=$((LINES % 4))

    echo "${READS}|${REMAINDER}|${EXIT_CODE}|${LINES}"
}

check_library() {
    local LIB_DIR="$1"
    local LIB_NAME=$(basename "${LIB_DIR}")

    echo "=== ${LIB_NAME} ==="
    echo ""

    if [ ! -d "${LIB_DIR}" ]; then
        echo "  Directory not found: ${LIB_DIR}"
        echo ""
        return
    fi

    # List files with sizes first (instant)
    echo "  File sizes:"
    for FQ in "${LIB_DIR}"/*.fastq.gz; do
        [ -f "${FQ}" ] || continue
        FNAME=$(basename "${FQ}")
        FSIZE=$(stat -c%s "${FQ}" 2>/dev/null || echo "0")
        FSIZE_MB=$((FSIZE / 1024 / 1024))
        READ_TYPE=$(echo "${FNAME}" | grep -oP '(R[12]|I[12])(?=_001)' || echo "??")
        echo "    ${READ_TYPE}: ${FSIZE_MB} MB  (${FNAME})"
    done
    echo ""

    # Check gzip integrity (fast — doesn't decompress fully)
    echo "  Gzip integrity:"
    for FQ in "${LIB_DIR}"/*.fastq.gz; do
        [ -f "${FQ}" ] || continue
        FNAME=$(basename "${FQ}")
        READ_TYPE=$(echo "${FNAME}" | grep -oP '(R[12]|I[12])(?=_001)' || echo "??")
        if gzip -t "${FQ}" 2>/dev/null; then
            echo "    ${READ_TYPE}: OK"
        else
            echo "    ${READ_TYPE}: *** CORRUPT ***"
        fi
    done
    echo ""

    # Count reads (slow for large files — this is the critical step)
    echo "  Counting reads (this may take several minutes for large files)..."
    declare -A COUNTS
    MIN_READS=999999999999
    MAX_READS=0

    for FQ in "${LIB_DIR}"/*.fastq.gz; do
        [ -f "${FQ}" ] || continue
        FNAME=$(basename "${FQ}")
        READ_TYPE=$(echo "${FNAME}" | grep -oP '(R[12]|I[12])(?=_001)' || echo "??")

        echo -n "    ${READ_TYPE}: counting... "
        RESULT=$(count_reads "${FQ}")
        READS=$(echo "${RESULT}" | cut -d'|' -f1)
        REMAINDER=$(echo "${RESULT}" | cut -d'|' -f2)
        EXIT_CODE=$(echo "${RESULT}" | cut -d'|' -f3)

        COUNTS["${READ_TYPE}"]="${READS}"

        STATUS=""
        if [ "${EXIT_CODE}" -ne 0 ]; then
            STATUS=" (decompression error!)"
        elif [ "${REMAINDER}" -ne 0 ]; then
            STATUS=" (INCOMPLETE: ${REMAINDER} trailing lines)"
        fi

        echo "${READS} reads${STATUS}"

        if [ "${READS}" -lt "${MIN_READS}" ]; then
            MIN_READS="${READS}"
        fi
        if [ "${READS}" -gt "${MAX_READS}" ]; then
            MAX_READS="${READS}"
        fi
    done
    echo ""

    # Report
    if [ "${MIN_READS}" -eq "${MAX_READS}" ]; then
        echo "  RESULT: All files match — ${MIN_READS} reads each. No repair needed."
    else
        DIFF=$((MAX_READS - MIN_READS))
        PCT_LOST=$(echo "scale=4; ${DIFF} * 100 / ${MAX_READS}" | bc)
        echo "  RESULT: *** MISMATCH ***"
        echo "    Max reads: ${MAX_READS}"
        echo "    Min reads: ${MIN_READS}"
        echo "    Difference: ${DIFF} reads (${PCT_LOST}% of total)"
        echo ""
        echo "  TO REPAIR (truncate longer files to match shortest):"
        echo "    bash repair_truncated_fastqs.sh ${LIB_DIR} ${MIN_READS}"
        echo ""
        echo "  This is safe because 10x FASTQs are ordered — read N in R1"
        echo "  pairs with read N in R2/I1/I2. The missing reads are all at"
        echo "  the tail end. Truncating preserves all complete mate pairs."
    fi

    unset COUNTS
    echo ""
    echo ""
}

# MD5 check (fast, do this first for all libraries)
echo "--- Quick MD5 check on affected files ---"
echo ""

for LIB in PROTEIN_LIBRARY_B VDJ_RH_LIBRARY_H; do
    LIB_DIR="${FASTQ_BASE}/${LIB}"
    MD5_FILE="${LIB_DIR}/MD5.txt"

    echo "  ${LIB}:"
    if [ -f "${MD5_FILE}" ]; then
        while IFS= read -r line; do
            # Handle both "md5  filename" and "md5 filename" formats
            EXPECTED=$(echo "${line}" | awk '{print $1}')
            FNAME=$(echo "${line}" | awk '{print $NF}')
            FPATH="${LIB_DIR}/${FNAME}"

            if [ -f "${FPATH}" ]; then
                ACTUAL=$(md5sum "${FPATH}" | awk '{print $1}')
                if [ "${EXPECTED}" = "${ACTUAL}" ]; then
                    echo "    ${FNAME}: MATCH (original/intact)"
                else
                    echo "    ${FNAME}: *** MISMATCH *** (modified or corrupted)"
                fi
            else
                echo "    ${FNAME}: file not found"
            fi
        done < "${MD5_FILE}"
    else
        echo "    No MD5.txt found"
    fi

    # Check for backup directory
    if [ -d "${LIB_DIR}/BACKUP" ]; then
        echo "    BACKUP/ directory exists — originals may be saved there"
        ls -lh "${LIB_DIR}/BACKUP/" 2>/dev/null | grep -v "^total" | awk '{print "      " $NF ": " $5}'
    fi
    echo ""
done

echo ""
echo "--- Full read count analysis (may take 10-20 min per library) ---"
echo ""

# Check the two failed libraries
check_library "${FASTQ_BASE}/PROTEIN_LIBRARY_B"
check_library "${FASTQ_BASE}/VDJ_RH_LIBRARY_H"

echo "============================================================"
echo " SUMMARY"
echo "============================================================"
echo ""
echo " If BACKUP/ contains the original corrupted files and the"
echo " current files have been trimmed to end on complete records,"
echo " the repair script will make all files match the shortest."
echo ""
echo " After repair, rerun only the two failed libraries:"
echo "   sbatch --array=1,7 CELLRANGER_MULTI_RH_SHIV.sh"
echo ""
echo " (array=1 is Library B = SHIV_Pre_LN)"
echo " (array=7 is Library H = SHIV_Necropsy_LN)"
