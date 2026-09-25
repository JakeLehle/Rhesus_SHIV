#!/bin/bash
#SBATCH --job-name=shiv_blastn
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=160
#SBATCH --mem=900G
#SBATCH --time=7-00:00:00
#SBATCH --output=/master/jlehle/WORKING/LOGS/shiv_blastn_%j.out

# %% Cell 1 — parameters (all adjustable knobs live here)
set -o pipefail

CONTIGS="/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV/shiv_megahit_assembly/megahit_out/final.contigs.fa"
DB_ROOT="/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV/viral_id_databases"
NUCL_DB="${DB_ROOT}/blast_nucl_db/lentivirus_nucl"
OUTDIR="${DB_ROOT}/search_results"
THREADS=160
EVALUE="1e-5"

# No max_target_seqs on purpose. e-value is the selection. Full table.

mkdir -p "${OUTDIR}"

OUT_TSV="${OUTDIR}/blastn_contigs_vs_lentivirus.tsv"

# %% Cell 2 — environment
source ~/anaconda3/bin/activate
eval "$(conda shell.bash hook)"
conda activate viral_id
export PATH="${CONDA_PREFIX}/bin:${PATH}"

echo "=================================================="
echo "SHIV contig blastn search"
echo "Started:        $(date)"
echo "blastn:         $(which blastn)"
echo "blastn version: $(blastn -version | head -1)"
echo "Query contigs:  ${CONTIGS}"
echo "DB:             ${NUCL_DB}"
echo "Threads:        ${THREADS}"
echo "E-value:        ${EVALUE}"
echo "Output:         ${OUT_TSV}"
echo "=================================================="

# %% Cell 3 — sanity checks before the long run
if [ ! -s "${CONTIGS}" ]; then
    echo "ERROR: contigs FASTA missing or empty: ${CONTIGS}"
    exit 1
fi
if [ ! -f "${NUCL_DB}.nin" ] && [ ! -f "${NUCL_DB}.nal" ]; then
    echo "ERROR: nucleotide DB not found at ${NUCL_DB} (no .nin/.nal)"
    exit 1
fi

NCONTIGS=$(grep -c '^>' "${CONTIGS}")
echo "Contig count in query: ${NCONTIGS}"
echo ""

# %% Cell 4 — blastn (the long pole)
# Column set: alignment geometry + qlen/slen for coverage + staxids/stitle
# so we can read SIV-half vs HIV-half and gene region at analysis time.
echo "### Running blastn ... $(date) ###"

blastn \
    -query "${CONTIGS}" \
    -db "${NUCL_DB}" \
    -evalue "${EVALUE}" \
    -num_threads "${THREADS}" \
    -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen staxids stitle' \
    -out "${OUT_TSV}"

RC=$?
echo ""
if [ ${RC} -ne 0 ]; then
    echo "ERROR: blastn exited with code ${RC} at $(date)"
    echo "Partial output (if any) is at: ${OUT_TSV}"
    exit ${RC}
fi

# %% Cell 5 — quick post-run summary (not analysis, just a pulse check)
NLINES=$(wc -l < "${OUT_TSV}")
NQ=$(cut -f1 "${OUT_TSV}" | sort -u | wc -l)
echo "=================================================="
echo "blastn complete: $(date)"
echo "Total hit rows:        ${NLINES}"
echo "Distinct contigs hit:  ${NQ} of ${NCONTIGS}"
echo "Output: ${OUT_TSV}"
echo "=================================================="
