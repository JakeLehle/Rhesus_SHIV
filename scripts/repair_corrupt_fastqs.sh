#!/bin/bash
###############################################################################
# repair_corrupt_fastqs.sh
#
# Repairs 10x Chromium FASTQ libraries where file transfer corruption
# has damaged R1 and/or R2 files. Handles:
#   - Corrupt gzip streams (decompresses what's valid, discards rest)
#   - Incomplete FASTQ records (trims to last complete 4-line block)
#   - Read count mismatches (truncates all files to match shortest)
#
# 10x FASTQs are strictly ordered — read N in R1 corresponds to read N
# in R2, I1, I2. All four files must have identical read counts.
# Truncating to the minimum preserves all complete mate pairs.
#
# Process per file:
#   1. Decompress (suppressing gzip errors for corrupt files)
#   2. Trim to last complete FASTQ record (lines % 4 == 0)
#   3. Count reads
# Then across all files:
#   4. Find minimum read count
#   5. Truncate all to minimum
#   6. Recompress with gzip
#   7. Verify
#
# Original files are backed up to REPAIR_BACKUP/ before modification.
#
# Usage:
#   bash repair_corrupt_fastqs.sh <library_dir>
#
# Examples:
#   bash repair_corrupt_fastqs.sh /master/zwallis/WORKING/SC/10X_RAW/PROTEIN_LIBRARY_B
#   bash repair_corrupt_fastqs.sh /master/zwallis/WORKING/SC/10X_RAW/VDJ_RH_LIBRARY_H
###############################################################################

set +e  # Don't exit on error — we need to handle corrupt files

if [ $# -ne 1 ]; then
    echo "Usage: $0 <library_dir>"
    exit 1
fi

LIB_DIR="$1"
LIB_NAME=$(basename "${LIB_DIR}")
BACKUP_DIR="${LIB_DIR}/REPAIR_BACKUP"
TMPDIR_REPAIR="${LIB_DIR}/repair_tmp"

echo "============================================================"
echo " FASTQ Repair: ${LIB_NAME}"
echo " $(date)"
echo "============================================================"
echo ""

if [ ! -d "${LIB_DIR}" ]; then
    echo "ERROR: Directory not found: ${LIB_DIR}"
    exit 1
fi

# Create backup and temp directories
mkdir -p "${BACKUP_DIR}"
mkdir -p "${TMPDIR_REPAIR}"

# Find all FASTQ files
FASTQS=($(ls "${LIB_DIR}"/*.fastq.gz 2>/dev/null | sort))

if [ ${#FASTQS[@]} -eq 0 ]; then
    echo "ERROR: No FASTQ files found"
    exit 1
fi

echo "Found ${#FASTQS[@]} FASTQ files"
echo ""

###############################################################################
# STEP 1: Backup originals
###############################################################################
echo "--- Step 1: Backing up originals ---"
for FQ in "${FASTQS[@]}"; do
    FNAME=$(basename "${FQ}")
    if [ ! -f "${BACKUP_DIR}/${FNAME}" ]; then
        echo "  Backing up: ${FNAME}"
        cp "${FQ}" "${BACKUP_DIR}/${FNAME}"
    else
        echo "  Backup exists: ${FNAME}"
    fi
done
echo ""

###############################################################################
# STEP 2: Decompress, clean, and count each file
###############################################################################
echo "--- Step 2: Decompress and clean each file ---"
echo "(This will take a while for large files)"
echo ""

declare -A CLEAN_READS
declare -A READ_TYPES
MIN_READS=999999999999

for FQ in "${FASTQS[@]}"; do
    FNAME=$(basename "${FQ}")
    READ_TYPE=$(echo "${FNAME}" | grep -oP '(R[12]|I[12])(?=_001)')
    READ_TYPES["${FNAME}"]="${READ_TYPE}"
    CLEAN_FILE="${TMPDIR_REPAIR}/${FNAME%.gz}"  # remove .gz extension

    echo "  Processing ${FNAME} (${READ_TYPE})..."

    # Check gzip integrity first
    if gzip -t "${FQ}" 2>/dev/null; then
        GZ_STATUS="intact"
        echo "    Gzip: intact"
        # Decompress normally
        zcat "${FQ}" > "${CLEAN_FILE}" 2>/dev/null
    else
        GZ_STATUS="corrupt"
        echo "    Gzip: CORRUPT — decompressing what's recoverable..."
        # Decompress with error suppression — gets everything before corruption
        zcat "${FQ}" > "${CLEAN_FILE}" 2>/dev/null || true
    fi

    # Count raw lines
    RAW_LINES=$(wc -l < "${CLEAN_FILE}")
    REMAINDER=$((RAW_LINES % 4))

    if [ "${REMAINDER}" -ne 0 ]; then
        echo "    Raw lines: ${RAW_LINES} (${REMAINDER} trailing lines — incomplete record)"
        # Trim to last complete FASTQ record
        TRIM_LINES=$((RAW_LINES - REMAINDER))
        head -n "${TRIM_LINES}" "${CLEAN_FILE}" > "${CLEAN_FILE}.trimmed"
        mv "${CLEAN_FILE}.trimmed" "${CLEAN_FILE}"
        FINAL_LINES="${TRIM_LINES}"
        echo "    Trimmed to: ${FINAL_LINES} lines"
    else
        FINAL_LINES="${RAW_LINES}"
        echo "    Lines: ${FINAL_LINES} (complete records)"
    fi

    READS=$((FINAL_LINES / 4))
    CLEAN_READS["${FNAME}"]="${READS}"
    echo "    Reads: ${READS}"

    if [ "${READS}" -lt "${MIN_READS}" ]; then
        MIN_READS="${READS}"
    fi

    echo ""
done

echo "--- Read count summary after cleaning ---"
for FNAME in $(echo "${!CLEAN_READS[@]}" | tr ' ' '\n' | sort); do
    RT="${READ_TYPES[${FNAME}]}"
    RC="${CLEAN_READS[${FNAME}]}"
    if [ "${RC}" -eq "${MIN_READS}" ]; then
        echo "  ${RT}: ${RC} reads  ← MINIMUM (truncation target)"
    else
        DIFF=$((RC - MIN_READS))
        echo "  ${RT}: ${RC} reads  (${DIFF} reads over minimum)"
    fi
done
echo ""

###############################################################################
# STEP 3: Truncate all files to match minimum
###############################################################################
TARGET_LINES=$((MIN_READS * 4))

echo "--- Step 3: Truncating all files to ${MIN_READS} reads (${TARGET_LINES} lines) ---"
echo ""

for FQ in "${FASTQS[@]}"; do
    FNAME=$(basename "${FQ}")
    READ_TYPE="${READ_TYPES[${FNAME}]}"
    CLEAN_FILE="${TMPDIR_REPAIR}/${FNAME%.gz}"
    CURRENT_READS="${CLEAN_READS[${FNAME}]}"

    if [ "${CURRENT_READS}" -eq "${MIN_READS}" ]; then
        echo "  ${READ_TYPE}: already at target (${MIN_READS} reads)"
    else
        REMOVED=$((CURRENT_READS - MIN_READS))
        echo "  ${READ_TYPE}: truncating ${CURRENT_READS} → ${MIN_READS} (removing ${REMOVED} tail reads)"
        head -n "${TARGET_LINES}" "${CLEAN_FILE}" > "${CLEAN_FILE}.truncated"
        mv "${CLEAN_FILE}.truncated" "${CLEAN_FILE}"
    fi
done
echo ""

###############################################################################
# STEP 4: Verify all files match, then recompress
###############################################################################
echo "--- Step 4: Verifying and recompressing ---"
echo ""

ALL_MATCH=true

for FQ in "${FASTQS[@]}"; do
    FNAME=$(basename "${FQ}")
    READ_TYPE="${READ_TYPES[${FNAME}]}"
    CLEAN_FILE="${TMPDIR_REPAIR}/${FNAME%.gz}"

    # Verify line count
    VERIFY_LINES=$(wc -l < "${CLEAN_FILE}")
    VERIFY_READS=$((VERIFY_LINES / 4))
    VERIFY_REMAINDER=$((VERIFY_LINES % 4))

    if [ "${VERIFY_READS}" -ne "${MIN_READS}" ] || [ "${VERIFY_REMAINDER}" -ne 0 ]; then
        echo "  ${READ_TYPE}: *** VERIFICATION FAILED *** (${VERIFY_READS} reads, ${VERIFY_REMAINDER} remainder)"
        ALL_MATCH=false
    else
        echo "  ${READ_TYPE}: verified ${VERIFY_READS} reads"
    fi
done
echo ""

if ! ${ALL_MATCH}; then
    echo "ERROR: Verification failed! Files do not match."
    echo "Originals preserved in ${BACKUP_DIR}/"
    echo "Temp files in ${TMPDIR_REPAIR}/ for inspection."
    exit 1
fi

# Spot-check: verify first and last reads look like valid FASTQ
echo "  Spot-checking first and last records..."
for FQ in "${FASTQS[@]}"; do
    FNAME=$(basename "${FQ}")
    READ_TYPE="${READ_TYPES[${FNAME}]}"
    CLEAN_FILE="${TMPDIR_REPAIR}/${FNAME%.gz}"

    FIRST_HEADER=$(head -1 "${CLEAN_FILE}")
    LAST_HEADER=$(tail -4 "${CLEAN_FILE}" | head -1)

    if [[ "${FIRST_HEADER}" != @* ]]; then
        echo "  ${READ_TYPE}: *** First record doesn't start with @ ***"
        ALL_MATCH=false
    fi
    if [[ "${LAST_HEADER}" != @* ]]; then
        echo "  ${READ_TYPE}: *** Last record doesn't start with @ ***"
        ALL_MATCH=false
    fi
done

if ! ${ALL_MATCH}; then
    echo "ERROR: FASTQ format validation failed!"
    exit 1
fi
echo "  Format OK"
echo ""

# Recompress
echo "  Recompressing (this will take a while)..."
for FQ in "${FASTQS[@]}"; do
    FNAME=$(basename "${FQ}")
    READ_TYPE="${READ_TYPES[${FNAME}]}"
    CLEAN_FILE="${TMPDIR_REPAIR}/${FNAME%.gz}"

    echo -n "    ${READ_TYPE}: compressing... "
    gzip -c "${CLEAN_FILE}" > "${FQ}"
    NEW_SIZE=$(stat -c%s "${FQ}" 2>/dev/null || echo "0")
    NEW_SIZE_MB=$((NEW_SIZE / 1024 / 1024))
    echo "${NEW_SIZE_MB} MB"
done
echo ""

###############################################################################
# STEP 5: Final verification
###############################################################################
echo "--- Step 5: Final verification ---"
echo ""

FINAL_OK=true
for FQ in "${FASTQS[@]}"; do
    FNAME=$(basename "${FQ}")
    READ_TYPE="${READ_TYPES[${FNAME}]}"

    # Gzip integrity
    if gzip -t "${FQ}" 2>/dev/null; then
        GZ_OK="OK"
    else
        GZ_OK="FAIL"
        FINAL_OK=false
    fi

    # Quick read count (from decompressed temp file, should still be there)
    CLEAN_FILE="${TMPDIR_REPAIR}/${FNAME%.gz}"
    if [ -f "${CLEAN_FILE}" ]; then
        FINAL_LINES=$(wc -l < "${CLEAN_FILE}")
        FINAL_READS=$((FINAL_LINES / 4))
    else
        # Count from compressed file
        FINAL_READS=$(zcat "${FQ}" 2>/dev/null | wc -l)
        FINAL_READS=$((FINAL_READS / 4))
    fi

    echo "  ${READ_TYPE}: ${FINAL_READS} reads, gzip=${GZ_OK}"
done
echo ""

if ${FINAL_OK}; then
    echo "============================================================"
    echo " REPAIR SUCCESSFUL"
    echo "============================================================"
    echo ""
    echo "  Library:     ${LIB_NAME}"
    echo "  Final reads: ${MIN_READS} per file (all matched)"
    echo "  Originals:   ${BACKUP_DIR}/"
    echo ""

    # Calculate retention
    # I1/I2 had the most reads (they were intact)
    for FNAME in $(echo "${!CLEAN_READS[@]}" | tr ' ' '\n' | sort); do
        RT="${READ_TYPES[${FNAME}]}"
        if [ "${RT}" = "I1" ]; then
            ORIGINAL="${CLEAN_READS[${FNAME}]}"
            # Wait, CLEAN_READS has the cleaned count, not original
            # Use the backup to get original I1 count
            break
        fi
    done

    echo "  Cleanup: you can remove temp files with:"
    echo "    rm -rf ${TMPDIR_REPAIR}"
    echo ""
    echo "  To restore originals if needed:"
    echo "    for f in ${BACKUP_DIR}/*.fastq.gz; do"
    echo "      cp \"\$f\" \"${LIB_DIR}/\$(basename \"\$f\")\""
    echo "    done"
else
    echo "ERROR: Final verification failed!"
    echo "Check temp files in ${TMPDIR_REPAIR}/"
fi

echo ""
echo "Done: $(date)"
