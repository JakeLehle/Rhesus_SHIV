#!/usr/bin/env bash
# =============================================================================
# run_shiv_vpr_barcode_gt_search.sh
# =============================================================================
# One job, no arrays. Streams the raw GEX FASTQs, searches for the Binhua
# ground-truth barcodes across five matching tiers, links to the object, writes
# the reservoir map and a per-tier negative floor. Read-only.
#
#   sbatch run_shiv_vpr_barcode_gt_search.sh
#
# SMOKE TEST FIRST. Edit Cell 1 of the Python script:
#     LIBS = ['E', 'H']
#     MAX_READS_PER_LIB = 10_000_000
# Confirm the workbook parse, the index build and the CB join all behave, then
# revert to all 8 libraries and MAX_READS_PER_LIB = None for the full pass.
#
# Requires Binhua_FINAL_Barcode_4-7-26_Rhesus_SHIV_Barcodes.xlsx to be present in
# /master/jlehle/WORKING/SC/fastq/Rhesus_SHIV before launching.
# =============================================================================
#SBATCH --job-name=shiv_bc_gt
#SBATCH --partition=normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=300G
#SBATCH --time=4-00:00:00
#SBATCH --output=/master/jlehle/WORKING/LOGS/shiv_bc_gt_search_%j.out

set -u
source ~/anaconda3/bin/activate sc_pre
export PATH="${CONDA_PREFIX}/bin:${PATH}"

WORKDIR=/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV
XLSX="${WORKDIR}/Binhua_FINAL_Barcode_4-7-26_Rhesus_SHIV_Barcodes.xlsx"

echo "host    : $(hostname)"
echo "python  : $(which python)"
echo "pigz    : $(which pigz 2>/dev/null || echo 'not found (falls back to zcat)')"
echo "started : $(date)"

# Warn and continue, per convention. The script degrades to Tier 1 only without
# the workbook, which would quietly reproduce the old result and look like a
# negative finding rather than a missing input.
if [[ ! -f "${XLSX}" ]]; then
    echo "*** WARNING: ground-truth workbook not found at ${XLSX}"
    echo "*** Tiers 2-5 all require it. Without it this run is Tier 1 only and"
    echo "*** its output is NOT comparable to the prior demux. Copy the file in"
    echo "*** and relaunch unless you specifically intend a Tier 1 rerun."
else
    echo "workbook: ${XLSX}"
fi

cd "${WORKDIR}"
python SHIV_vpr_barcode_gt_search.py

echo "finished: $(date)"
