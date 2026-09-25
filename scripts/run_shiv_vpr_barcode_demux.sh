#!/usr/bin/env bash
# =============================================================================
# run_shiv_vpr_barcode_demux.sh
# =============================================================================
# One job, no arrays. Streams the raw GEX FASTQs, feature-demuxes the vpr
# cassette, links to the object, writes the reservoir map. Read-only.
#
#   sbatch run_shiv_vpr_barcode_demux.sh
#
# Smoke test first: edit the script's Cell 1 to LIBS=['E','H'] and
# MAX_READS_PER_LIB=10_000_000, confirm the plumbing, then revert to all 8.
# =============================================================================
#SBATCH --job-name=shiv_bc_demux
#SBATCH --partition=normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=300G
#SBATCH --time=3-00:00:00
#SBATCH --output=/master/jlehle/WORKING/LOGS/shiv_bc_demux_%j.out

set -u
source ~/anaconda3/bin/activate sc_pre
export PATH="${CONDA_PREFIX}/bin:${PATH}"

echo "host   : $(hostname)"
echo "python : $(which python)"
echo "pigz   : $(which pigz 2>/dev/null || echo 'not found (falls back to zcat)')"
echo "started: $(date)"

cd /master/jlehle/WORKING/SC/fastq/Rhesus_SHIV
python shiv_vpr_barcode_demux.py

echo "finished: $(date)"
