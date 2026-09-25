#!/bin/bash
#SBATCH --job-name=shiv_candidates_vs_ref
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/master/jlehle/WORKING/LOGS/shiv_candidates_vs_ref_%j.out

# %% Cell 1 — parameters
set -o pipefail

BASE="/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV"
CONTIGS="${BASE}/shiv_megahit_assembly/megahit_out/final.contigs.fa"
RESULTS="${BASE}/viral_id_databases/search_results"
SUMMARY="${RESULTS}/summary"
CANDS="${SUMMARY}/candidate_viral_contigs.txt"
OUTDIR="${SUMMARY}/ref_mapping"
REF_ACC="MN816822.1"
THREADS=16

mkdir -p "${OUTDIR}"

REF_FASTA="${OUTDIR}/SHIVAD8EO_${REF_ACC}.fasta"
REF_GB="${OUTDIR}/SHIVAD8EO_${REF_ACC}.gb"
CAND_FASTA="${OUTDIR}/candidate_contigs.fasta"
OUT_TSV="${OUTDIR}/candidates_vs_SHIVAD8EO.tsv"

# %% Cell 2 — environment
source ~/anaconda3/bin/activate
eval "$(conda shell.bash hook)"
conda activate viral_id
export PATH="${CONDA_PREFIX}/bin:${PATH}"
if [ -f ~/.ncbi_api_key ]; then source ~/.ncbi_api_key; export NCBI_API_KEY; fi

echo "=================================================="
echo "Candidate contigs vs SHIVAD8-EO (${REF_ACC})"
echo "Started: $(date)"
echo "=================================================="

# %% Cell 3 — fetch the reference: FASTA for blast, GenBank for the gene map
echo "### Fetching ${REF_ACC} (FASTA + GenBank) ###"
efetch -db nuccore -id "${REF_ACC}" -format fasta    > "${REF_FASTA}" 2>"${OUTDIR}/efetch_fasta.err"
efetch -db nuccore -id "${REF_ACC}" -format gb       > "${REF_GB}"    2>"${OUTDIR}/efetch_gb.err"
if [ ! -s "${REF_FASTA}" ] || [ ! -s "${REF_GB}" ]; then
    echo "ERROR: reference fetch failed. Check efetch err files in ${OUTDIR}"
    exit 1
fi
REFLEN=$(grep -v '^>' "${REF_FASTA}" | tr -d '\n' | wc -c)
echo "  reference length: ${REFLEN} bp"

# %% Cell 4 — extract the 706 candidate contigs by name
echo "### Extracting candidate contigs ###"
NCAND=$(wc -l < "${CANDS}")
echo "  candidate IDs in list: ${NCAND}"
seqkit grep -f "${CANDS}" "${CONTIGS}" > "${CAND_FASTA}" 2>"${OUTDIR}/seqkit.err"
NGOT=$(grep -c '^>' "${CAND_FASTA}")
echo "  contigs extracted:     ${NGOT}"
if [ "${NGOT}" -ne "${NCAND}" ]; then
    echo "  WARNING: extracted count != list count. Check ${OUTDIR}/seqkit.err"
fi

# %% Cell 5 — make a tiny blast db of the single reference and blast candidates
echo "### Building single-reference blast db + blasting candidates ###"
makeblastdb -in "${REF_FASTA}" -dbtype nucl -out "${OUTDIR}/ref_db" \
    > "${OUTDIR}/makeblastdb.log" 2>&1

# permissive e-value; no max_target_seqs (only one subject anyway).
# coordinate columns are the whole point: sstart/send map onto the gene map.
blastn \
    -query "${CAND_FASTA}" \
    -db "${OUTDIR}/ref_db" \
    -evalue 1e-5 \
    -num_threads "${THREADS}" \
    -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' \
    -out "${OUT_TSV}"

RC=$?
NHITS=$(wc -l < "${OUT_TSV}")
NQ=$(cut -f1 "${OUT_TSV}" | sort -u | wc -l)
echo "=================================================="
echo "blastn rc=${RC} | hit rows: ${NHITS} | contigs with a hit to ref: ${NQ} of ${NGOT}"
echo "Outputs in ${OUTDIR}"
echo "Done: $(date)"
echo "=================================================="
