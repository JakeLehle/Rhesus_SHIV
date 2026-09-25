#!/bin/bash
#SBATCH --job-name=shiv_comp_diag
#SBATCH --output=/master/jlehle/WORKING/LOGS/shiv_comp_diag_%j.out
#SBATCH --error=/master/jlehle/WORKING/LOGS/shiv_comp_diag_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --partition=normal

###############################################################################
# run_shiv_comprehensive_diagnostic.sh
#
# Submits the read-only object-level diagnostic (Script A) on the canonical
# embedded object. Nothing is written back to the object; only CSV summaries
# land in annotation_output/comprehensive_diagnostic/.
#
# Usage:  sbatch run_shiv_comprehensive_diagnostic.sh
#
# Author: Jake Lehle
# Date: July 2026
###############################################################################

set -euo pipefail

# ─── conda / env ─────────────────────────────────────────────────────────────
source ~/anaconda3/bin/activate
eval "$(conda shell.bash hook)"
conda activate sc_pre
export PATH="${CONDA_PREFIX}/bin:${PATH}"   # conda tools ahead of the broken system ones

# ─── run ─────────────────────────────────────────────────────────────────────
SCRIPT_DIR="/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV"
cd "${SCRIPT_DIR}"

echo "============================================================"
echo "SHIV comprehensive diagnostic — $(date)"
echo "  host   : $(hostname)"
echo "  python : $(which python)"
echo "  env    : ${CONDA_DEFAULT_ENV}"
echo "============================================================"

python shiv_comprehensive_diagnostic.py

echo "============================================================"
echo "DONE — $(date)"
echo "  CSV summaries: ${SCRIPT_DIR}/annotation_output/comprehensive_diagnostic/"
echo "============================================================"
