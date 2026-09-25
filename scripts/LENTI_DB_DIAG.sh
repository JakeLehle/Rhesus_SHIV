#!/bin/bash
#SBATCH --job-name=lenti_db_diag
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --time=12:00:00
#SBATCH --output=/master/jlehle/WORKING/LOGS/lenti_db_diag_%j.out

# %% Cell 1 — parameters (all adjustable knobs live here)
set -o pipefail

DB_ROOT="/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV/viral_id_databases"
NUCL_DB="${DB_ROOT}/blast_nucl_db/lentivirus_nucl"
OUTDIR="${DB_ROOT}/diagnostics"
THREADS=32

# References to probe redundancy depth with. Measurement cap is wide open
# on purpose so nothing gets clipped while we are measuring.
MEASURE_MAXTARGETS=5000
PIDENT_BINS="70 80 90 95 97 99"   # report hit counts at/above each identity

mkdir -p "${OUTDIR}"

# %% Cell 2 — environment
source ~/anaconda3/bin/activate
eval "$(conda shell.bash hook)"
conda activate viral_id
export PATH="${CONDA_PREFIX}/bin:${PATH}"

if [ -f ~/.ncbi_api_key ]; then
    source ~/.ncbi_api_key
    export NCBI_API_KEY
fi

echo "=================================================="
echo "Lentivirus DB diagnostic"
echo "Started: $(date)"
echo "blastn: $(which blastn)"
echo "blastn version: $(blastn -version | head -1)"
echo "DB: ${NUCL_DB}"
echo "=================================================="

# %% Cell 3 — Q1: lineage composition of the whole database
# Use the loose pattern from the notes. Literal SIVmac returns almost nothing
# because records title themselves "Simian immunodeficiency virus (Mac251)" etc.
echo ""
echo "### Q1: Database composition by lineage ###"

ALL_TITLES="${OUTDIR}/all_titles.txt"
blastdbcmd -db "${NUCL_DB}" -entry all -outfmt '%t' > "${ALL_TITLES}" 2>"${OUTDIR}/blastdbcmd.err"

TOTAL=$(wc -l < "${ALL_TITLES}")
SIV=$(grep -icE "simian immunodeficiency" "${ALL_TITLES}")
HIV1=$(grep -icE "HIV-1|human immunodeficiency virus 1|human immunodeficiency virus type 1" "${ALL_TITLES}")
SHIV=$(grep -icE "\bSHIV\b|simian.human immunodeficiency" "${ALL_TITLES}")
MAC=$(grep -icE "mac239|mac251|mac32|smm|sooty|\(mac\)" "${ALL_TITLES}")

echo "Total records:                 ${TOTAL}"
echo "Simian immunodeficiency (SIV): ${SIV}"
echo "HIV-1:                         ${HIV1}"
echo "SHIV (chimera deposits):       ${SHIV}"
echo "mac/smm/sooty lineage subset:  ${MAC}"

# %% Cell 4 — Q3: length distribution of database records
# (Doing length before the redundancy probe so the references are already
#  validated as present when we blast them.)
echo ""
echo "### Q3: Database record length distribution ###"

LENFILE="${OUTDIR}/seq_lengths.txt"
blastdbcmd -db "${NUCL_DB}" -entry all -outfmt '%l' > "${LENFILE}" 2>>"${OUTDIR}/blastdbcmd.err"

echo "Length summary (bp):"
sort -n "${LENFILE}" | awk '
{
    a[NR]=$1; sum+=$1
}
END {
    n=NR
    printf "  count:  %d\n", n
    printf "  min:    %d\n", a[1]
    printf "  median: %d\n", a[int(n/2)]
    printf "  mean:   %.1f\n", sum/n
    printf "  max:    %d\n", a[n]
    # how many full-genome-scale records exist
    full=0; for(i=1;i<=n;i++) if(a[i]>=8000) full++
    printf "  >=8kb (genome-scale): %d\n", full
    frag=0; for(i=1;i<=n;i++) if(a[i]<1000) frag++
    printf "  <1kb (fragments):     %d\n", frag
}'

# %% Cell 5 — Q2: redundancy depth a real contig will encounter
# Pull three references, blast each with the cap wide open, count how many
# DB entries each hits at/above a series of identity thresholds. This is the
# number that decides max_target_seqs for the real search.
echo ""
echo "### Q2: Redundancy depth per reference (cap wide open at ${MEASURE_MAXTARGETS}) ###"

declare -A REFS=(
    ["SHIVAD8EO"]="MN816822.1"
    ["SIVmac239"]="M33262.1"
    ["HIV1_HXB2"]="K03455.1"
)

REFFASTA="${OUTDIR}/refs.fasta"
> "${REFFASTA}"

for name in "${!REFS[@]}"; do
    acc="${REFS[$name]}"
    echo ""
    echo "--- fetching ${name} (${acc}) ---"
    efetch -db nuccore -id "${acc}" -format fasta > "${OUTDIR}/${name}.fasta" 2>"${OUTDIR}/${name}.efetch.err"
    if [ ! -s "${OUTDIR}/${name}.fasta" ]; then
        echo "WARNING: ${name} (${acc}) fetch returned empty. Check ${OUTDIR}/${name}.efetch.err"
        continue
    fi
    cat "${OUTDIR}/${name}.fasta" >> "${REFFASTA}"

    HITS="${OUTDIR}/${name}_vs_db.tsv"
    blastn -query "${OUTDIR}/${name}.fasta" -db "${NUCL_DB}" \
        -evalue 1e-5 -num_threads "${THREADS}" \
        -max_target_seqs "${MEASURE_MAXTARGETS}" \
        -outfmt '6 qseqid sseqid pident length evalue bitscore slen stitle' \
        > "${HITS}" 2>"${OUTDIR}/${name}.blastn.err"

    nhits=$(wc -l < "${HITS}")
    echo "${name}: ${nhits} total hits returned (cap ${MEASURE_MAXTARGETS})"
    if [ "${nhits}" -eq "${MEASURE_MAXTARGETS}" ]; then
        echo "  NOTE: hit the measurement cap exactly. True depth may exceed ${MEASURE_MAXTARGETS}."
    fi
    echo "  hits at/above each identity threshold:"
    for thr in ${PIDENT_BINS}; do
        c=$(awk -F'\t' -v t="${thr}" '$3>=t {n++} END{print n+0}' "${HITS}")
        printf "    >=%s%%: %d\n" "${thr}" "${c}"
    done
done

# %% Cell 6 — done
echo ""
echo "=================================================="
echo "Diagnostic complete: $(date)"
echo "Outputs in: ${OUTDIR}"
echo "Read Q2 hit-depth tables to set max_target_seqs for the real search."
echo "=================================================="
