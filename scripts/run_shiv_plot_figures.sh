#!/bin/bash
#SBATCH --job-name=shiv_plots
#SBATCH --output=/master/jlehle/WORKING/LOGS/shiv_plots_%j.out
#SBATCH --error=/master/jlehle/WORKING/LOGS/shiv_plots_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --partition=normal

###############################################################################
# run_shiv_plot_figures.sh
#
# Submits the read-only figure set (cell types, clusters, surface markers,
# metadata, composition) off shiv_host_final_embedded.h5ad. Figures land in
# annotation_output/figures/ as PDF + PNG at 300 DPI. The object is not modified.
#
# Usage:  sbatch run_shiv_plot_figures.sh
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
echo "SHIV plotting — $(date)"
echo "  host   : $(hostname)"
echo "  python : $(which python)"
echo "  env    : ${CONDA_DEFAULT_ENV}"
echo "============================================================"

python shiv_plot_figures.py

echo "============================================================"
echo "DONE — $(date)"
echo "  Figures: ${SCRIPT_DIR}/annotation_output/figures/"
echo "============================================================"
