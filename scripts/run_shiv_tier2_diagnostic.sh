#!/bin/bash
#SBATCH --job-name=shiv_tier2_diag
#SBATCH --output=/master/jlehle/WORKING/LOGS/shiv_tier2_diag_%j.out
#SBATCH --error=/master/jlehle/WORKING/LOGS/shiv_tier2_diag_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=04:00:00
#SBATCH --partition=normal

###############################################################################
# run_shiv_tier2_diagnostic.sh
#
# Submits the Tier 2 reclustering + signature diagnostic for one lineage
# (set LINEAGE at the top of the .py; T cells first). Reclusters the lineage
# with BBKNN, scores functional signatures, reports subcluster profiles and
# animal/timepoint/tissue/infection enrichment. Writes a NEW subset object plus
# CSVs and figures; the canonical object is not modified.
#
# Usage:  sbatch run_shiv_tier2_diagnostic.sh
#
# Author: Jake Lehle
# Date: July 2026
###############################################################################

set -euo pipefail

source ~/anaconda3/bin/activate
eval "$(conda shell.bash hook)"
conda activate sc_pre
export PATH="${CONDA_PREFIX}/bin:${PATH}"

SCRIPT_DIR="/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV"
cd "${SCRIPT_DIR}"

echo "============================================================"
echo "SHIV Tier 2 diagnostic — $(date)"
echo "  host   : $(hostname)"
echo "  python : $(which python)"
echo "  env    : ${CONDA_DEFAULT_ENV}"
echo "============================================================"

python shiv_tier2_diagnostic.py

echo "============================================================"
echo "DONE — $(date)"
echo "  Outputs under: ${SCRIPT_DIR}/annotation_output/tier2_*/"
echo "============================================================"
