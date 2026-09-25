#!/bin/bash
#SBATCH --job-name=shiv_diamond
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=160
#SBATCH --mem=900G
#SBATCH --time=7-00:00:00
#SBATCH --output=/master/jlehle/WORKING/LOGS/shiv_diamond_%j.out

# %% Cell 1 — parameters (all adjustable knobs live here)
set -o pipefail

CONTIGS="/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV/shiv_megahit_assembly/megahit_out/final.contigs.fa"
DB_ROOT="/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV/viral_id_databases"
PROT_DB="${DB_ROOT}/diamond_prot_db/lentivirus_prot.dmnd"
OUTDIR="${DB_ROOT}/search_results"
THREADS=160
EVALUE="1e-5"

# --very-sensitive for remote homology. No max-target-seqs cap (0 = unlimited).
# e-value is the selection, same philosophy as blastn.

mkdir -p "${OUTDIR}"

OUT_TSV="${OUTDIR}/diamond_blastx_contigs_vs_lentivirus.tsv"

# %% Cell 2 — environment
source ~/anaconda3/bin/activate
eval "$(conda shell.bash hook)"
conda activate viral_id
export PATH="${CONDA_PREFIX}/bin:${PATH}"

echo "=================================================="
echo "SHIV contig DIAMOND blastx search"
echo "Started:         $(date)"
echo "diamond:         $(which diamond)"
echo "diamond version: $(diamond version 2>&1 | head -1)"
echo "Query contigs:   ${CONTIGS}"
echo "DB:              ${PROT_DB}"
echo "Threads:         ${THREADS}"
echo "E-value:         ${EVALUE}"
echo "Output:          ${OUT_TSV}"
echo "=================================================="

# %% Cell 3 — sanity checks before the run
if [ ! -s "${CONTIGS}" ]; then
    echo "ERROR: contigs FASTA missing or empty: ${CONTIGS}"
    exit 1
fi
if [ ! -s "${PROT_DB}" ]; then
    echo "ERROR: DIAMOND protein DB missing or empty: ${PROT_DB}"
    exit 1
fi

NCONTIGS=$(grep -c '^>' "${CONTIGS}")
echo "Contig count in query: ${NCONTIGS}"
echo ""

# %% Cell 4 — DIAMOND blastx
# Column set mirrors the blastn table as closely as blastx allows, plus
# sframe so we can see which frame each hit landed in. stitle for region/parent.
echo "### Running DIAMOND blastx ... $(date) ###"

diamond blastx \
    --query "${CONTIGS}" \
    --db "${PROT_DB}" \
    --evalue "${EVALUE}" \
    --very-sensitive \
    --max-target-seqs 0 \
    --threads "${THREADS}" \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qframe stitle \
    --out "${OUT_TSV}"

RC=$?
echo ""
if [ ${RC} -ne 0 ]; then
    echo "ERROR: diamond exited with code ${RC} at $(date)"
    echo "Partial output (if any) is at: ${OUT_TSV}"
    exit ${RC}
fi

# %% Cell 5 — quick post-run summary
NLINES=$(wc -l < "${OUT_TSV}")
NQ=$(cut -f1 "${OUT_TSV}" | sort -u | wc -l)
echo "=================================================="
echo "DIAMOND blastx complete: $(date)"
echo "Total hit rows:        ${NLINES}"
echo "Distinct contigs hit:  ${NQ} of ${NCONTIGS}"
echo "Output: ${OUT_TSV}"
echo "=================================================="
