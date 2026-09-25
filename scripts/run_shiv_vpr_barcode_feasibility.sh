#!/bin/bash
#SBATCH --job-name=shiv_vpr_feas
#SBATCH --output=/master/jlehle/WORKING/LOGS/shiv_vpr_feas_%j.out
#SBATCH --error=/master/jlehle/WORKING/LOGS/shiv_vpr_feas_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=3-00:00:00
#SBATCH --partition=normal

###############################################################################
# run_shiv_vpr_barcode_feasibility.sh
#
# Submits the read-level VPR/VPX barcode recovery feasibility scan (Script B).
# Read-only: counts anchor hits across the STARsolo BAMs and the raw fastqs for
# all 8 libraries, cross-references recovered cores to the Binhua ground truth,
# and writes only CSV summaries. The full raw-fastq scan across all libraries
# (negatives included as the specificity control) is the long part; 24h is
# generous headroom.
#
# Usage:  sbatch run_shiv_vpr_barcode_feasibility.sh
#
# Author: Jake Lehle
# Date: July 2026
###############################################################################

set -euo pipefail

# ─── conda / env ─────────────────────────────────────────────────────────────
source ~/anaconda3/bin/activate
eval "$(conda shell.bash hook)"
conda activate sc_pre
export PATH="${CONDA_PREFIX}/bin:${PATH}"   # conda samtools + zcat ahead of the broken system ones

SCRIPT_DIR="/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV"
cd "${SCRIPT_DIR}"

echo "============================================================"
echo "SHIV VPR barcode feasibility — $(date)"
echo "  host     : $(hostname)"
echo "  python   : $(which python)"
echo "  samtools : $(which samtools || echo 'NOT FOUND')"
echo "  env      : ${CONDA_DEFAULT_ENV}"
echo "============================================================"

python shiv_vpr_barcode_feasibility.py

echo "============================================================"
echo "DONE — $(date)"
echo "  CSV summaries: ${SCRIPT_DIR}/annotation_output/vpr_barcode_feasibility/"
echo "============================================================"
