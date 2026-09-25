#!/usr/bin/env bash
# =============================================================================
# bd_abseq_scan.sh
# =============================================================================
# Step 2 of 3. Read-only against the raw FASTQs.
#
# Stage 0  read-length report for every file in the manifest (verify, don't assume)
# Stage 1  full-file scan of every R2 for the 16 BD probes  -> expected sparse
# Stage 2  subsampled scan of PROTEIN R2 for TotalSeq-C controls -> expected dense
#
# Stage 2 is what makes a zero in Stage 1 interpretable. Without a control that
# fires in the same file with the same machinery, "no hits" is an absence, not
# a finding.
#
# Single job, no arrays. Parallelism is internal via xargs -P so any failure is
# a single named output file that is missing or empty.
#
#   sbatch bd_abseq_scan.sh
# =============================================================================
#SBATCH --job-name=bd_abseq_scan
#SBATCH --partition=normal
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=/master/jlehle/WORKING/LOGS/bd_abseq_scan_%j.out

set -u

# ---------------------------------------------------------------------------
# Cell 1 - Parameters
# ---------------------------------------------------------------------------
OUT_DIR="/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV/adt_probe_scan"
RAW_DIR="${OUT_DIR}/raw"

# Source inputs from bd_abseq_probes.py. These may have CRLF line endings
# (Python's csv writer defaults to \r\n even on Linux). A trailing \r on a path
# makes the file "not exist"; a trailing \r on a probe makes grep -F search for
# SEQ\r and silently return zero. Both failure modes are invisible in the counts.
# So we strip \r from all three into cleaned working copies and use only those.
MANIFEST_SRC="${OUT_DIR}/scan_manifest.tsv"
PROBES_BD_SRC="${OUT_DIR}/probes_bd.txt"
PROBES_CTRL_SRC="${OUT_DIR}/probes_ctrl.txt"

MANIFEST="${RAW_DIR}/scan_manifest.clean.tsv"
PROBES_BD="${RAW_DIR}/probes_bd.clean.txt"
PROBES_CTRL="${RAW_DIR}/probes_ctrl.clean.txt"

NJOBS=8               # files scanned concurrently
DECOMP_THREADS=4      # pigz threads per file; NJOBS * DECOMP_THREADS <= cpus-per-task
CTRL_READS=2000000    # reads subsampled from each PROTEIN R2 for the control scan
CTRL_MAX_HITS=50000   # cap on control hit lines written (offsets only need a sample)
KEEP_CTRL_TEMP=0      # 1 to keep the subsampled sequence files

mkdir -p "${RAW_DIR}"

for f in "${MANIFEST_SRC}" "${PROBES_BD_SRC}" "${PROBES_CTRL_SRC}"; do
    [[ -s "${f}" ]] || { echo "MISSING or EMPTY: ${f}. Run bd_abseq_probes.py first."; exit 1; }
done

# Strip CRLF -> LF. Also drop any blank probe line, which would make grep -F
# match every read.
tr -d '\r' < "${MANIFEST_SRC}"   > "${MANIFEST}"
tr -d '\r' < "${PROBES_BD_SRC}"  | sed '/^$/d' > "${PROBES_BD}"
tr -d '\r' < "${PROBES_CTRL_SRC}" | sed '/^$/d' > "${PROBES_CTRL}"

echo "Sanitized inputs (CRLF stripped):"
echo "  manifest : ${MANIFEST}"
echo "  bd       : ${PROBES_BD}  ($(wc -l < "${PROBES_BD}") probes)"
echo "  ctrl     : ${PROBES_CTRL}  ($(wc -l < "${PROBES_CTRL}") probes)"

# Guard: confirm the first manifest path now opens, before launching 17 jobs.
first_path=$(awk -F'\t' 'NR==2 {print $5; exit}' "${MANIFEST}")
if [[ ! -f "${first_path}" ]]; then
    echo "*** First manifest path still not found after sanitization:"
    echo "    [${first_path}]"
    echo "    Check the manifest was written for this filesystem."
    exit 1
fi
echo "  first path opens OK: ${first_path}"

# Prefer pigz. Fall back to zcat, which is slower but correct.
if command -v pigz >/dev/null 2>&1; then
    DECOMP="pigz -dc -p ${DECOMP_THREADS}"
else
    DECOMP="zcat"
    echo "NOTE: pigz not found, falling back to zcat (single threaded)."
fi

echo "================================================================"
echo "BD AbSeq probe scan"
echo "  started   : $(date)"
echo "  host      : $(hostname)"
echo "  manifest  : ${MANIFEST}  ($(( $(wc -l < "${MANIFEST}") - 1 )) files)"
echo "  bd probes : $(wc -l < "${PROBES_BD}")"
echo "  ctrl probe: $(wc -l < "${PROBES_CTRL}")"
echo "  decomp    : ${DECOMP}"
echo "  out       : ${RAW_DIR}"
echo "================================================================"

export OUT_DIR RAW_DIR PROBES_BD PROBES_CTRL DECOMP CTRL_READS CTRL_MAX_HITS KEEP_CTRL_TEMP


# ---------------------------------------------------------------------------
# Cell 2 - Stage 0: read-length report
# ---------------------------------------------------------------------------
echo ""
echo "--- Stage 0: read lengths (first read of each file) ---"
RL="${RAW_DIR}/read_lengths.tsv"
printf 'tag\tlib_type\tfirst_read_len\tpath\n' > "${RL}"
tail -n +2 "${MANIFEST}" | while IFS=$'\t' read -r tag library lib_type rd path; do
    len=$(${DECOMP} "${path}" 2>/dev/null | head -n 4 | awk 'NR==2 {print length($0)}')
    printf '%s\t%s\t%s\t%s\n' "${tag}" "${lib_type}" "${len:-NA}" "${path}" >> "${RL}"
    printf '  %-14s %-13s %sbp\n' "${tag}" "${lib_type}" "${len:-NA}"
done
echo "  wrote ${RL}"
echo "  A PROTEIN R2 shorter than 46bp cannot hold a 36bp tag at offset 10."
echo "  A GEX R2 shorter than 36bp cannot hold the full-length probe at all."


# ---------------------------------------------------------------------------
# Cell 3 - Stage 1: full-file BD probe scan
# ---------------------------------------------------------------------------
# grep -n numbers lines of the sequence-only stream, so the line number IS the
# 1-based read index within the file. That lets us go back to the read later.
scan_one() {
    local tag="$1" path="$2"
    local hits="${RAW_DIR}/${tag}.bdhits.txt"
    local nreads="${RAW_DIR}/${tag}.nreads"

    ${DECOMP} "${path}" \
      | awk -v out="${nreads}" 'NR % 4 == 2 { n++; print } END { print (n + 0) > out }' \
      | LC_ALL=C grep -n -F -f "${PROBES_BD}" > "${hits}"

    local rc=$?
    # grep exit 1 = no match, which is a legitimate result here, not a failure.
    if [[ ${rc} -gt 1 ]]; then
        echo "  *** grep failed on ${tag} (rc=${rc})"
        return 1
    fi
    printf '  %-14s %12s reads   %8s probe hits\n' \
        "${tag}" "$(cat "${nreads}" 2>/dev/null || echo NA)" "$(wc -l < "${hits}")"
}
export -f scan_one

echo ""
echo "--- Stage 1: full-file BD probe scan (${NJOBS} files at a time) ---"
tail -n +2 "${MANIFEST}" | cut -f1,5 \
  | xargs -P "${NJOBS}" -n 2 bash -c 'scan_one "$0" "$1"'
echo "  Stage 1 complete."


# ---------------------------------------------------------------------------
# Cell 4 - Stage 2: subsampled positive-control scan (PROTEIN R2 only)
# ---------------------------------------------------------------------------
ctrl_one() {
    local tag="$1" path="$2"
    local sub="${RAW_DIR}/${tag}.ctrl_seqs.tmp"
    local hits="${RAW_DIR}/${tag}.ctrlhits.txt"
    local counts="${RAW_DIR}/${tag}.ctrlcounts.tsv"

    ${DECOMP} "${path}" 2>/dev/null \
      | head -n $(( CTRL_READS * 4 )) \
      | awk 'NR % 4 == 2' > "${sub}"

    local n_sub
    n_sub=$(wc -l < "${sub}")

    LC_ALL=C grep -n -F -f "${PROBES_CTRL}" "${sub}" | head -n "${CTRL_MAX_HITS}" > "${hits}"

    : > "${counts}"
    while read -r seq; do
        [[ -z "${seq}" ]] && continue
        c=$(LC_ALL=C grep -c -F "${seq}" "${sub}" || true)
        printf '%s\t%s\t%s\t%s\n' "${tag}" "${seq}" "${c}" "${n_sub}" >> "${counts}"
    done < "${PROBES_CTRL}"

    [[ "${KEEP_CTRL_TEMP}" -eq 1 ]] || rm -f "${sub}"
    printf '  %-14s %12s reads sampled   %8s control hits\n' \
        "${tag}" "${n_sub}" "$(wc -l < "${hits}")"
}
export -f ctrl_one

echo ""
echo "--- Stage 2: positive-control scan, first ${CTRL_READS} reads of each PROTEIN R2 ---"
awk -F'\t' 'NR > 1 && $3 == "PROTEIN" { print $1"\t"$5 }' "${MANIFEST}" \
  | xargs -P "${NJOBS}" -n 2 bash -c 'ctrl_one "$0" "$1"'
echo "  Stage 2 complete."


# ---------------------------------------------------------------------------
# Cell 5 - Wrap up
# ---------------------------------------------------------------------------
echo ""
echo "================================================================"
echo "  finished : $(date)"
echo "  raw outputs in ${RAW_DIR}"
echo "  Next: python bd_abseq_summarize.py"
echo "================================================================"
