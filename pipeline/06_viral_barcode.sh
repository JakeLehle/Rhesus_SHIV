#!/bin/bash
#SBATCH -J SHIV_06_barcode
#SBATCH -o /master/jlehle/WORKING/LOGS/shiv_06_viral_barcode_%j.o.log
#SBATCH -e /master/jlehle/WORKING/LOGS/shiv_06_viral_barcode_%j.e.log
#SBATCH -t 01:00:00
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem 64GB

###############################################################################
# Submitter for 06_viral_barcode.py
#
# WHY THIS ONE IS TINY COMPARED TO 03 AND 05
# -------------------------------------------
# The old barcode work scanned raw FASTQ at full depth across eight libraries
# and took hours. This stage reads only the records already aligned to the SHIV
# contig. Measured on library H, the richest of the eight:
#
#     sample_alignments.bam        862 records (431 pairs),  13 near the cassette
#     unassigned_alignments.bam     90 records ( 45 pairs),   0 near the cassette
#
# So the whole eight-library sweep is a few thousand SAM lines. 1 hour of
# walltime is slack, not an estimate, and it will almost certainly finish in
# under a minute.
#
# The 64GB is NOT for the SAM parsing. It is for Cell 9, which opens
# adata_annotated.h5ad to map cell barcodes onto animal, timepoint, tissue and
# cell type. It opens it backed='r' and pulls only the obs columns it needs, so
# the real peak is modest, but h5ad reads have surprised us before and this
# stage is cheap enough that there is no reason to run it tight.
#
# 4 cores is for samtools, which is the only thing here that threads.
#
# SAMTOOLS ON TITAN
# -----------------
# /usr/local/bin/samtools on titan is a broken build and will fail with a
# library error partway through. The conda one in sc_pre works. The PATH export
# below puts CONDA_PREFIX/bin first for exactly this reason, and the check
# further down prints which binary won so a wrong one shows up in the log
# rather than in a confusing traceback.
#
# Usage:  sbatch run_06_viral_barcode.sh
#         then read detection_funnel.csv top to bottom before anything else.
#
# Author: Jake Lehle
# Date:   November 2026
###############################################################################

set -euo pipefail

PIPELINE_DIR="/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV/Rhesus_SHIV/pipeline"
RESULTS_DIR="/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV/pipeline_results"
SCRIPT="06_viral_barcode.py"
CONDA_ENV="sc_pre"

echo "=============================================================="
echo "SHIV 06 — viral barcode detection from aligned reads"
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

# samtools is the hard dependency of this stage, and WHICH samtools matters.
echo
echo "  samtools check:"
SAMTOOLS="$(command -v samtools || true)"
if [ -z "${SAMTOOLS}" ]; then
    echo "    ERROR no samtools on PATH. conda install -c bioconda samtools"
    exit 1
fi
echo "    binary  : ${SAMTOOLS}"
if [ "${SAMTOOLS}" = "/usr/local/bin/samtools" ]; then
    echo "    ERROR this is the broken titan build. CONDA_PREFIX/bin did not"
    echo "          win the PATH race. Check that ${CONDA_ENV} actually has"
    echo "          samtools installed."
    exit 1
fi
if ! samtools --version >/dev/null 2>&1; then
    echo "    ERROR samtools is on PATH but will not run (library error?)."
    samtools --version || true
    exit 1
fi
echo "    version : $(samtools --version | head -1)"

echo
echo "  python dependency check:"
python - <<'PY'
import importlib, sys
need = ['numpy', 'pandas', 'matplotlib', 'anndata']
soft = {'openpyxl': 'the Binhua catalog cannot be read; detection still runs, '
                    'the animal cross-check is skipped'}
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

# Note what is deliberately NOT required: pysam. The script shells out to
# samtools through subprocess instead, because pysam is not confirmed present
# in sc_pre and there is no reason to add a build dependency to parse a few
# thousand SAM lines.
echo "    pysam        not required (samtools via subprocess by design)"

echo
echo "  inputs:"

# 03's BAMs. Globbed, because the script globs too and a layout change should
# surface here rather than as an empty result.
MULTI_DIR="${RESULTS_DIR}/03_cellranger_multi"
if [ ! -d "${MULTI_DIR}" ]; then
    echo "    ERROR ${MULTI_DIR} not found. Run 03_cellranger_multi.sh first."
    exit 1
fi
N_BAM=$(find "${MULTI_DIR}" -name '*_alignments.bam' 2>/dev/null | wc -l)
N_BAI=$(find "${MULTI_DIR}" -name '*_alignments.bam.bai' 2>/dev/null | wc -l)
echo "    BAMs    : ${N_BAM} found under 03_cellranger_multi/"
echo "    indexes : ${N_BAI}"
if [ "${N_BAM}" -eq 0 ]; then
    echo "    ERROR no *_alignments.bam anywhere under ${MULTI_DIR}."
    exit 1
fi
if [ "${N_BAI}" -lt "${N_BAM}" ]; then
    echo "    WARN  $(( N_BAM - N_BAI )) BAM(s) are unindexed. A region query"
    echo "          needs an index, so those will be skipped with a logged"
    echo "          message. Index them with: samtools index <bam>"
fi

# 05's annotated object. Optional: without it the cores are still extracted and
# counted, they just cannot be attributed to an animal, so the negative floor
# loses the control-animal half of its definition. Warn, do not stop.
ANNOT="${RESULTS_DIR}/05_embed_annotate/adata_annotated.h5ad"
if [ -f "${ANNOT}" ]; then
    echo "    annot   : ${ANNOT}"
    echo "              $(stat -c '%s bytes, modified %y' "${ANNOT}")"
else
    echo "    WARN  ${ANNOT} not found."
    echo "          Detection runs, but molecules cannot be linked to animals,"
    echo "          timepoints or tissues, and the control-animal half of the"
    echo "          negative floor is unavailable. Run 05 first for the full"
    echo "          result."
fi

# The catalog, which is validation only and never gates a call.
CAT=""
for c in "/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV/Binhua_FINAL_Barcode_4-7-26_Rhesus_SHIV_Barcodes.xlsx" \
         "/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV/Binhua_FINAL_Barcode_4726.xlsx"; do
    [ -f "${c}" ] && CAT="${c}" && break
done
if [ -n "${CAT}" ]; then
    echo "    catalog : $(basename "${CAT}")  (validation only)"
else
    echo "    WARN  no Binhua catalog found. Detection is unaffected; the"
    echo "          animal cross-check is skipped."
fi

cd "${PIPELINE_DIR}"
[ -f "${SCRIPT}" ] || { echo "ERROR ${SCRIPT} not found in ${PIPELINE_DIR}"; exit 1; }

# Echo the settings that change what counts, so a log alone tells you which
# thresholds produced the numbers in it.
echo
echo "  settings:"
# The '=' is required. Without it the pattern matched prose headings in the
# module docstring ("MIN_CORE_BASES IS 10, NOT 6", "CASSETTE GEOMETRY,
# VERIFIED NOT ASSUMED") and the log reported those instead of the values.
for k in CONTIG CASSETTE CORE MIN_CORE_BASES MIN_FLANK_ALIGNED \
         EXACT_MATCH_ONLY PROJECT_SOFTCLIPS USE_MLUI_AS_ANCHOR \
         WINDOW_SLACK REGISTER_MARKERS; do
    line=$(grep -m1 -E "^${k} *=" "${SCRIPT}" | sed 's/  *#.*//')
    echo "    ${line:-${k} unset}"
done
echo

python -u "${SCRIPT}"
RC=$?

echo
echo "=============================================================="
echo "  finished: $(date -u +%Y-%m-%dT%H:%M:%SZ)   exit ${RC}"
echo "  outputs : ${RESULTS_DIR}/06_viral_barcode/"
echo
echo "  READ IN THIS ORDER:"
echo "    detection_funnel.csv      where the reads go and where they stop"
echo "    cassette_geometry.csv     whether the core is readable at all"
echo "    negative_floor.csv        the false-positive rate; check it BEFORE"
echo "                              believing any positive call"
echo "    barcode_inventory.csv     per animal, timepoint, tissue"
echo "    barcode_persistence.csv   the longitudinal question"
echo
echo "  The verdict block at the end of the log says whether targeted VpxF1"
echo "  enrichment is needed. That is the decision this stage exists to make."
echo "=============================================================="
exit ${RC}
