#!/bin/bash
#SBATCH -J SHIV_07_report
#SBATCH -o /master/jlehle/WORKING/LOGS/shiv_07_object_report_%j.o.log
#SBATCH -e /master/jlehle/WORKING/LOGS/shiv_07_object_report_%j.e.log
#SBATCH -t 02:00:00
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem 128GB

###############################################################################
# Submitter for 07_object_report.py
#
# WALLTIME AND MEMORY
# -------------------
# This stage reads obs and var and never touches the matrix, so it is cheap:
# expect a couple of minutes. 2 hours is slack.
#
# 128GB is sized for the obs read on a 1.6GB h5ad with ~100 columns across
# 116,536 cells, plus the crosstabs. The object is opened backed='r' and X is
# never loaded, so the steady state is far below this.
#
# NO FIGURES
# ----------
# 07 produces no plots and does not import matplotlib. That is deliberate:
# the figure set is a separate decision, and automating a figure nobody has
# chosen yet produces work that gets thrown away. The report this writes is
# the input to that decision.
#
# Usage:  sbatch run_07_object_report.sh
#         then read pipeline_results/07_object_report/07_object_report.md
#
# Author: Jake Lehle
# Date:   November 2026
###############################################################################

set -euo pipefail

PIPELINE_DIR="/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV/Rhesus_SHIV/pipeline"
RESULTS_DIR="/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV/pipeline_results"
SCRIPT="07_object_report.py"
CONDA_ENV="sc_pre"

echo "=============================================================="
echo "SHIV 07 — object report (no figures)"
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

echo
echo "  dependency check:"
python - <<'PY'
import importlib, sys
need = ['numpy', 'pandas', 'anndata']
bad = []
for m in need:
    try:
        mod = importlib.import_module(m)
        print(f"    {m:<12} {getattr(mod, '__version__', 'ok')}")
    except ImportError:
        bad.append(m)
        print(f"    {m:<12} MISSING  (required)")
if bad:
    sys.exit(f"  ERROR required packages missing: {bad}")
# Stated rather than implied: this stage deliberately has no plotting
# dependency, so a missing or broken matplotlib cannot stop the report.
print("    matplotlib   not required (07 produces no figures by design)")
print("    scanpy       not required (obs/var only, no matrix operations)")
PY

echo
echo "  inputs:"
ANNOT="${RESULTS_DIR}/05_embed_annotate/adata_annotated.h5ad"
if [ ! -f "${ANNOT}" ]; then
    echo "    ERROR ${ANNOT} not found. Run 05_embed_annotate.py first."
    exit 1
fi
echo "    annot   : ${ANNOT}"
echo "              $(stat -c '%s bytes, modified %y' "${ANNOT}")"

# 06's inventory is optional: without it the report simply omits the barcode
# section and says so, rather than failing.
INV="${RESULTS_DIR}/06_viral_barcode/barcode_inventory.csv"
if [ -f "${INV}" ]; then
    echo "    barcodes: ${INV} ($(wc -l < "${INV}") lines)"
else
    echo "    WARN  ${INV} not found; the barcode section will be omitted."
    echo "          Run 06_viral_barcode.py for a complete report."
fi

cd "${PIPELINE_DIR}"
[ -f "${SCRIPT}" ] || { echo "ERROR ${SCRIPT} not found in ${PIPELINE_DIR}"; exit 1; }

echo
echo "  thresholds that decide the answerable-questions section:"
for k in MIN_CELLS_PER_GROUP MIN_CELLS_PER_CELLTYPE CONTROL_ANIMALS; do
    line=$(grep -m1 -E "^${k} *=" "${SCRIPT}" | sed 's/  *#.*//')
    echo "    ${line:-${k} unset}"
done
echo

python -u "${SCRIPT}"
RC=$?

echo
echo "=============================================================="
echo "  finished: $(date -u +%Y-%m-%dT%H:%M:%SZ)   exit ${RC}"
echo "  outputs : ${RESULTS_DIR}/07_object_report/"
echo
echo "  THE REPORT:  07_object_report.md   <- read this, top to bottom"
echo
echo "  Sections 5 (the design), 6 (composition) and 11 (what the object"
echo "  can answer) are the ones that decide what is worth plotting."
echo "  Every table is also a CSV in the same directory."
echo "=============================================================="
exit ${RC}
