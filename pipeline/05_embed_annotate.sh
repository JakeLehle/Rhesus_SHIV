#!/bin/bash
#SBATCH -J SHIV_05_annotate
#SBATCH -o /master/jlehle/WORKING/LOGS/shiv_05_embed_annotate_%j.o.log
#SBATCH -e /master/jlehle/WORKING/LOGS/shiv_05_embed_annotate_%j.e.log
#SBATCH -t 24:00:00
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem 256GB

###############################################################################
# Submitter for 05_embed_annotate.py
#
# WALLTIME AND MEMORY
# -------------------
# The long poles on 115,402 cells x 18,582 genes are the two UMAPs (the
# pre-BBKNN one and the corrected one), Leiden, and then three within-lineage
# subclusterings each with their own PCA, BBKNN, UMAP, Leiden and a Wilcoxon
# rank_genes_groups. Expect roughly 4 to 8 hours. 24h is slack, not an
# estimate.
#
# 256GB is sized for the peak, which is the HVG step: it copies the matrix
# minus the excluded features, so two dense-ish copies coexist briefly. The
# steady state is far lower.
#
# A RESUME_FROM_TIER1 = True run skips the embedding entirely and finishes in
# minutes, so you can submit it with -t 02:00:00 --mem 128GB if you want to be
# polite to the queue.
#
# THREADS
# -------
# scanpy is largely single-threaded outside BLAS, so 16 cores is for the BLAS
# calls in PCA and the neighbor search, not for a parallel speedup across the
# whole script. The script reads SLURM_CPUS_PER_TASK for sc.settings.n_jobs,
# so this number and the one scanpy uses cannot drift apart.
#
# Usage:  sbatch run_05_embed_annotate.sh
#         then read the log at /master/jlehle/WORKING/LOGS/
#         tail -f it if you want to watch; python runs with -u so the log
#         fills as it goes rather than only at exit.
#
# Author: Jake Lehle
# Date:   November 2026
###############################################################################

set -euo pipefail

PIPELINE_DIR="/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV/Rhesus_SHIV/pipeline"
RESULTS_DIR="/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV/pipeline_results"
SCRIPT="05_embed_annotate.py"
CONDA_ENV="sc_pre"

echo "=============================================================="
echo "SHIV 05 — embedding, Tier 1, TCR recovery, Tier 2"
echo "  job     : ${SLURM_JOB_ID:-interactive}"
echo "  node    : $(hostname)"
echo "  cpus    : ${SLURM_CPUS_PER_TASK:-unset}"
echo "  mem     : ${SLURM_MEM_PER_NODE:-unset} MB"
echo "  started : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "=============================================================="

source ~/anaconda3/bin/activate "${CONDA_ENV}"
export PATH="${CONDA_PREFIX}/bin:${PATH}"

echo "  conda   : ${CONDA_PREFIX}"
echo "  python  : $(command -v python)"

# Fail fast on a missing dependency rather than 5 hours in, after the
# embedding is already done, at the import inside the Tier 2 loop.
echo
echo "  dependency check:"
python - <<'PY'
import importlib, sys
need = ['numpy', 'pandas', 'scanpy', 'anndata', 'scipy', 'matplotlib']
soft = {'bbknn': 'batch correction falls back to uncorrected neighbors',
        'scirpy': 'the whole TCR step is skipped; no CD3-dropout recovery',
        'igraph': 'Leiden falls back to the slower legacy flavor',
        'leidenalg': 'needed by the legacy Leiden flavor'}
bad = []
for m in need:
    try:
        mod = importlib.import_module(m)
        print(f"    {m:<12} {getattr(mod, '__version__', 'ok')}")
    except ImportError:
        bad.append(m)
        print(f"    {m:<12} MISSING  (required)")
for m, consequence in soft.items():
    try:
        mod = importlib.import_module(m)
        print(f"    {m:<12} {getattr(mod, '__version__', 'ok')}")
    except ImportError:
        print(f"    {m:<12} MISSING  ->  {consequence}")
if bad:
    sys.exit(f"  ERROR required packages missing: {bad}")
PY
echo

# The input has to exist, and it has to be the post-04 object rather than the
# prefilter one, which looks similar and is a different thing.
QC_H5AD="${RESULTS_DIR}/04_ingest_qc/adata_qc.h5ad"
if [ ! -f "${QC_H5AD}" ]; then
    echo "  ERROR ${QC_H5AD} not found. Run 04_ingest_qc.py first."
    exit 1
fi
echo "  input   : ${QC_H5AD}"
echo "            $(stat -c '%s bytes, modified %y' "${QC_H5AD}")"

cd "${PIPELINE_DIR}"
[ -f "${SCRIPT}" ] || { echo "ERROR ${SCRIPT} not found in ${PIPELINE_DIR}"; exit 1; }

# Echo the resume setting, because a run that silently resumed when you meant
# a full rebuild is a confusing thing to debug from a log.
echo "  resume  : $(grep -m1 '^RESUME_FROM_TIER1' "${SCRIPT}" || echo 'unset')"
echo "  batch   : $(grep -m1 '^BATCH_KEY' "${SCRIPT}" || echo 'unset')"
echo

# -u so the log fills as the script runs rather than only at exit.
python -u "${SCRIPT}"
RC=$?

echo
echo "=============================================================="
echo "  finished: $(date -u +%Y-%m-%dT%H:%M:%SZ)   exit ${RC}"
echo "  outputs : ${RESULTS_DIR}/05_embed_annotate/"
echo "  read    : tier2_*_subtype_table.csv  (the review tables)"
echo "            tier1_transition_matrix.csv"
echo "            figures/F05_umap_prebbknn_by_library vs F05_umap_by_library"
echo "=============================================================="
exit ${RC}
