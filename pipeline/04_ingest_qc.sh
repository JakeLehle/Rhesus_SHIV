#!/bin/bash
#SBATCH -J SHIV_04_ingest_qc
#SBATCH -o /master/jlehle/WORKING/LOGS/shiv_04_ingest_qc_%j.o.log
#SBATCH -e /master/jlehle/WORKING/LOGS/shiv_04_ingest_qc_%j.e.log
#SBATCH -t 08:00:00
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem 300GB

###############################################################################
# Submitter for 04_ingest_qc.py
#
# 04 runs fine interactively, so this exists for two cases:
#   1. REUSE_PREFILTER = False, where re-ingesting all 8 libraries is the
#      expensive part and the peak memory is highest
#   2. your session drops and you would rather collect a log
#
# A run with REUSE_PREFILTER = True takes minutes, since it loads
# adata_prefilter.h5ad and goes straight to QC.
#
# 300GB is for case 1: eight filtered matrices plus the concatenation exist at
# the same time before anything is filtered. Scrublet is also per-library and
# builds a simulated-doublet matrix.
#
# Usage:  sbatch run_04_ingest_qc.sh
#
# Author: Jake Lehle
# Date:   November 2026
###############################################################################

set -euo pipefail

PIPELINE_DIR="/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV/Rhesus_SHIV/pipeline"
RESULTS_DIR="/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV/pipeline_results"
SCRIPT="04_ingest_qc.py"
CONDA_ENV="sc_pre"

echo "=============================================================="
echo "SHIV 04 — ingestion and QC"
echo "  job     : ${SLURM_JOB_ID:-interactive}"
echo "  node    : $(hostname)"
echo "  cpus    : ${SLURM_CPUS_PER_TASK:-unset}"
echo "  started : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "=============================================================="

source ~/anaconda3/bin/activate "${CONDA_ENV}"
export PATH="${CONDA_PREFIX}/bin:${PATH}"

echo "  conda   : ${CONDA_PREFIX}"
echo "  python  : $(command -v python)"

cd "${PIPELINE_DIR}"
[ -f "${SCRIPT}" ] || { echo "ERROR ${SCRIPT} not found in ${PIPELINE_DIR}"; exit 1; }

# The two settings that decide what this run actually does. Echo them, because
# a run that silently reused the prefilter, or silently deleted doublets, is a
# confusing thing to reconstruct from a log three weeks later.
echo
echo "  reuse prefilter : $(grep -m1 '^REUSE_PREFILTER' "${SCRIPT}" || echo 'unset')"
echo "  remove doublets : $(grep -m1 '^REMOVE_DOUBLETS' "${SCRIPT}" || echo 'unset')"
echo "  mito threshold  : $(grep -m1 '^MT_THRESHOLD' "${SCRIPT}" || echo 'unset')"
echo

PRE="${RESULTS_DIR}/04_ingest_qc/adata_prefilter.h5ad"
if [ -f "${PRE}" ]; then
    echo "  prefilter exists: ${PRE}"
    echo "                    $(stat -c '%s bytes, modified %y' "${PRE}")"
else
    echo "  no prefilter object yet; this will be a full re-ingestion."
fi
echo

python -u "${SCRIPT}"
RC=$?

echo
echo "=============================================================="
echo "  finished: $(date -u +%Y-%m-%dT%H:%M:%SZ)   exit ${RC}"
echo "  outputs : ${RESULTS_DIR}/04_ingest_qc/"
echo "  next    : sbatch run_05_embed_annotate.sh"
echo "=============================================================="
exit ${RC}
