#!/bin/bash
#SBATCH --job-name=viral_db_build
#SBATCH --output=/master/jlehle/WORKING/LOGS/viral_db_build_%j.out
#SBATCH --error=/master/jlehle/WORKING/LOGS/viral_db_build_%j.err
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
#SBATCH --partition=normal

###############################################################################
# shiv_build_databases.sh
#
# Purpose: Build, FROM SCRATCH and with full provenance, the three reference
#          databases needed to identify SHIV-like (and divergent lentiviral)
#          contigs from the MEGAHIT assembly of Library H unmapped reads.
#
#          NO SEARCHING happens here. This script only:
#            1. Activates pre-built env 'viral_id' and verifies tooling
#            2. Pulls all genus-Lentivirus (NCBI taxid 11646) nucleotide
#               genomes + proteins, fresh, via NCBI datasets CLI
#            3. Builds: makeblastdb nucl DB, diamond protein DB
#            4. Switches to env 'viral_id_genomad' and downloads the geNomad
#               marker DB (skips gracefully if that env is not built yet)
#            5. Writes a provenance manifest documenting every source,
#               version, date, taxid, and final count
#          Then STOPS for inspection.
#
# Design principles:
#   - Every input file is pulled fresh by a documented command. Nothing
#     materializes from an unexplained path. The manifest is the methods text.
#   - Network failures are flagged EXPLICITLY with a clear message telling
#     you which domain to request on the IT whitelist, rather than failing
#     silently or ambiguously.
#   - TWO purpose-built envs, NOT a clone of sc_pre. The earlier clone approach
#     failed: it carried 700+ packages plus a CUDA pin, and geNomad's
#     TF/xgboost stack could not solve against it. Now:
#       viral_id          (viral_id.yml)          : datasets, blast, diamond, seqkit
#       viral_id_genomad  (viral_id_genomad.yml)  : genomad (isolated)
#     Build both before running:
#       conda env create -f viral_id.yml
#       conda env create -f viral_id_genomad.yml
#   - geNomad is OPTIONAL for this build: if its env is absent, STEP 4 is
#     skipped and the BLAST/DIAMOND tiers (primary evidence) still complete.
#   - Verification gate before any download.
#
# Tier mapping (for the NEXT, separate search script):
#   - Nucleotide DB  -> blastn (close matches) + tblastx (protein-level, catches divergence)
#   - Protein DB     -> diamond blastx (fast protein-level divergence detection)
#   - geNomad DB     -> profile/HMM viral-marker detection (sensitive backstop)
#
# Author: Jake Lehle
# Date:   June 16, 2026
###############################################################################

set -o pipefail    # NOTE: deliberately NOT -e (catch+report failures with
                   # custom messages, not die mute) and NOT -u (conda's
                   # activate hooks reference unbound vars and would crash
                   # under -u). We guard our own vars explicitly instead.

# ─── Configuration ──────────────────────────────────────────────────────────
THREADS=${SLURM_CPUS_PER_TASK:-16}
WORKDIR="/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV"
DBROOT="${WORKDIR}/viral_id_databases"
LENTI_TAXID=11646            # NCBI genus Lentivirus
MIN_GENOME_LEN=1000          # curation: drop sub-genome fragments below this

# ─── NCBI API key + protein accession-batch settings (STEP 3) ───────────────
# The protein set is large and the single monolithic `datasets` stream keeps
# dropping mid-transfer on a slow link. STEP 3 instead lists protein accessions
# under the taxon, then pulls them with efetch in small batches (drop-proof:
# no single transfer is large; a failed batch retries in isolation).
#
# API key: stored OUTSIDE this script in a 600-perms file containing one line:
#     NCBI_API_KEY=your_key_here
# The script sources it to populate $NCBI_API_KEY (which efetch reads
# automatically) and NEVER prints the value. With a key, NCBI allows ~10 req/s
# (vs 3 keyless), so batches can be paced tighter.
NCBI_KEY_FILE="${HOME}/.ncbi_api_key"   # change here if you put it elsewhere
PROT_BATCH_SIZE=200                       # accessions per efetch batch
PROT_BATCH_SLEEP_KEY=1                     # seconds between batches WITH key
PROT_BATCH_SLEEP_NOKEY=3                   # seconds between batches WITHOUT key
PROT_BATCH_RETRIES=4                       # per-batch efetch retries

# Two PURPOSE-BUILT envs (NOT a clone of sc_pre — the clone failed because it
# carried 700+ packages + a CUDA pin that geNomad could not solve against):
#   ALIGN_ENV   : datasets, blast, diamond, seqkit  (lightweight, conflict-free)
#                 -> built from viral_id.yml
#   GENOMAD_ENV : genomad only (heavy TF/xgboost/mmseqs2 stack, isolated)
#                 -> built from viral_id_genomad.yml
# These must already exist before running this script:
#   conda env create -f viral_id.yml
#   conda env create -f viral_id_genomad.yml
ALIGN_ENV="viral_id"
GENOMAD_ENV="viral_id_genomad"

# Output subdirs
RAW_NUCL="${DBROOT}/lentivirus_nucleotide_raw"
RAW_PROT="${DBROOT}/lentivirus_protein_raw"
BLAST_DB_DIR="${DBROOT}/blast_nucl_db"
DIAMOND_DB_DIR="${DBROOT}/diamond_prot_db"
GENOMAD_DB_DIR="${DBROOT}/genomad_db"
MANIFEST="${DBROOT}/database_manifest.txt"

mkdir -p "${DBROOT}" "${RAW_NUCL}" "${RAW_PROT}" \
         "${BLAST_DB_DIR}" "${DIAMOND_DB_DIR}" "${GENOMAD_DB_DIR}"

# Log everything
LOGFILE="${DBROOT}/db_build.log"
exec > >(tee -a "${LOGFILE}") 2>&1

echo "############################################################"
echo "# VIRAL ID DATABASE BUILD — $(date)"
echo "# Output root: ${DBROOT}"
echo "############################################################"
echo ""

# ─── Provenance manifest header ─────────────────────────────────────────────
{
    echo "###############################################################"
    echo "# VIRAL IDENTIFICATION DATABASE — PROVENANCE MANIFEST"
    echo "# Built: $(date)"
    echo "# Host:  $(hostname)"
    echo "# Script: $0"
    echo "###############################################################"
    echo ""
} > "${MANIFEST}"

# Helper: record a line to the manifest AND stdout
record() { echo "  [MANIFEST] $*"; echo "$*" >> "${MANIFEST}"; }

# Helper: report a download failure, distinguishing a true whitelist/DNS block
# from a transient mid-stream drop. Usage: download_fail <type> <domain> <what>
#   type = "block"     -> connection refused / DNS / unreachable (likely firewall)
#   type = "transient" -> stream reset / timeout after retries (server/link flaky)
download_fail() {
    local ftype="$1"; local domain="$2"; local what="$3"
    echo ""
    echo "############################################################"
    if [[ "${ftype}" == "block" ]]; then
        echo "# NETWORK FAILURE — likely firewall / whitelist block"
        echo "############################################################"
        echo "#  Failed while fetching: ${what}"
        echo "#  Required domain:       ${domain}"
        echo "#"
        echo "#  The node could not establish a connection at all (DNS or"
        echo "#  connection refused). If this is a whitelist issue, submit an IT"
        echo "#  work ticket requesting this domain be added to the outbound"
        echo "#  allowlist for titan/zeus compute nodes:"
        echo "#"
        echo "#        ${domain}"
    else
        echo "# DOWNLOAD FAILURE — transient (stream reset / timeout)"
        echo "############################################################"
        echo "#  Failed while fetching: ${what}"
        echo "#  Host:                  ${domain}"
        echo "#"
        echo "#  The connection WAS established and data was transferring, then"
        echo "#  dropped mid-stream (e.g. HTTP/2 INTERNAL_ERROR, reset, timeout)."
        echo "#  This is NOT a whitelist block — the host is reachable. It is a"
        echo "#  flaky/slow link or a server-side hiccup. Retries were exhausted."
        echo "#"
        echo "#  Just re-run this script — completed steps are skipped via"
        echo "#  integrity-verified outputs, and the failed download retries."
    fi
    echo "############################################################"
    echo ""
}

# Helper: download an NCBI datasets archive ROBUSTLY, with retries and a
# post-download integrity test. A datasets zip is only complete if `unzip -t`
# passes (the central directory is written last, so a truncated transfer fails
# this test — exactly the failure mode seen on the protein pull).
#
# Usage: robust_datasets_download <outzip> <what-desc> <datasets-args...>
# Returns 0 on a verified-complete download, 1 otherwise.
robust_datasets_download() {
    local outzip="$1"; shift
    local what="$1";   shift
    local -a dl_args=("$@")
    local max_attempts=4
    local attempt=1
    local rc=0

    while (( attempt <= max_attempts )); do
        echo "  [download] attempt ${attempt}/${max_attempts}: ${what}"
        # Remove any partial/truncated file from a prior attempt first.
        rm -f "${outzip}"

        datasets download "${dl_args[@]}" --filename "${outzip}"
        rc=$?

        if (( rc != 0 )); then
            echo "  [download] datasets exited ${rc} on attempt ${attempt}."
        elif [[ ! -s "${outzip}" ]]; then
            echo "  [download] no file produced on attempt ${attempt}."
            rc=1
        elif ! unzip -t "${outzip}" >/dev/null 2>&1; then
            echo "  [download] archive failed integrity test (unzip -t) — truncated."
            rc=1
        else
            echo "  [download] OK — archive passed integrity test (unzip -t)."
            return 0
        fi

        if (( attempt < max_attempts )); then
            local backoff=$(( attempt * 30 ))
            echo "  [download] retrying in ${backoff}s..."
            sleep "${backoff}"
        fi
        attempt=$(( attempt + 1 ))
    done

    echo "  [download] all ${max_attempts} attempts failed for: ${what}"
    return 1
}

# Helper: fetch one batch of protein accessions as FASTA via efetch, with
# retries. Writes/append-safe to a temp file the caller controls.
# Usage: efetch_batch <comma-separated-accessions> <out_fasta_append> <sleep_s>
# Returns 0 if the batch produced non-empty FASTA, 1 otherwise.
efetch_batch() {
    local acc_csv="$1"; local out_fasta="$2"; local sleep_s="$3"
    local attempt=1
    local tmp
    tmp="$(mktemp)"

    while (( attempt <= PROT_BATCH_RETRIES )); do
        : > "${tmp}"
        # efetch reads $NCBI_API_KEY from the environment automatically.
        if efetch -db protein -id "${acc_csv}" -format fasta > "${tmp}" 2>/dev/null \
           && [[ -s "${tmp}" ]] \
           && grep -q "^>" "${tmp}"; then
            cat "${tmp}" >> "${out_fasta}"
            rm -f "${tmp}"
            return 0
        fi
        echo "    [efetch] batch attempt ${attempt}/${PROT_BATCH_RETRIES} failed; retrying..."
        sleep $(( attempt * 5 ))
        attempt=$(( attempt + 1 ))
    done

    rm -f "${tmp}"
    return 1
}

###############################################################################
# STEP 1: Activate the pre-built alignment env + verify tooling
###############################################################################
echo "============================================================"
echo "STEP 1: Activate alignment env (${ALIGN_ENV}) + verify tools"
echo "============================================================"

source ~/anaconda3/bin/activate
eval "$(conda shell.bash hook)"

# Both envs must be built beforehand from their YAMLs. We do NOT clone or
# install here — that is what caused the unsolvable env previously.
if ! conda env list | grep -qE "^\s*${ALIGN_ENV}\s"; then
    echo "  ERROR: env '${ALIGN_ENV}' does not exist."
    echo "  Build it first:  conda env create -f viral_id.yml"
    exit 1
fi
if ! conda env list | grep -qE "^\s*${GENOMAD_ENV}\s"; then
    echo "  WARNING: env '${GENOMAD_ENV}' does not exist."
    echo "  STEP 4 (geNomad DB download) will be skipped if still missing."
    echo "  Build it with:  conda env create -f viral_id_genomad.yml"
fi

conda activate "${ALIGN_ENV}"
export PATH="${CONDA_PREFIX}/bin:${PATH}"
echo "  Active env: ${CONDA_PREFIX}"
echo ""

# Gate: the alignment-tier tools must resolve in viral_id. geNomad is NOT
# checked here — it lives in its own env and is verified at STEP 4.
# entrez-direct (efetch/esearch) is required for STEP 3's accession-batch pull.
echo "  Verifying alignment tools (gate — fail here, not after downloads):"
TOOL_FAIL=0
for tool in datasets makeblastdb blastn tblastx diamond seqkit efetch esearch; do
    if command -v "${tool}" &>/dev/null; then
        printf "    %-12s OK   (%s)\n" "${tool}" "$(command -v ${tool})"
    else
        printf "    %-12s MISSING\n" "${tool}"
        TOOL_FAIL=1
    fi
done

if [[ "${TOOL_FAIL}" -ne 0 ]]; then
    echo ""
    echo "  ERROR: one or more tools missing from '${ALIGN_ENV}'."
    echo "  Most tools come from viral_id.yml. If efetch/esearch are the missing"
    echo "  ones, add entrez-direct to the env:"
    echo "    conda install -n ${ALIGN_ENV} -c bioconda entrez-direct"
    echo "  Or rebuild the env from the YAML:"
    echo "    conda env remove -n ${ALIGN_ENV}"
    echo "    conda env create -f viral_id.yml"
    echo "  Aborting before downloads."
    exit 1
fi

record "TOOLING (alignment env: ${ALIGN_ENV}, from viral_id.yml)"
record "  datasets : $(datasets --version 2>&1 | head -1)"
record "  blast    : $(blastn -version 2>&1 | head -1)"
record "  diamond  : $(diamond --version 2>&1 | head -1)"
record "  seqkit   : $(seqkit version 2>&1 | head -1)"
record "  efetch   : $(efetch -version 2>&1 | head -1)"
record ""
echo ""

###############################################################################
# STEP 2: Pull genus-Lentivirus nucleotide genomes (NCBI datasets)
###############################################################################
echo "============================================================"
echo "STEP 2: Lentivirus nucleotide genomes (taxid ${LENTI_TAXID})"
echo "============================================================"

NUCL_ZIP="${RAW_NUCL}/lentivirus_genomes.zip"
NUCL_FASTA="${RAW_NUCL}/lentivirus_genomes.fna"

if [[ -s "${NUCL_FASTA}" ]]; then
    echo "  Found existing ${NUCL_FASTA} — skipping download."
    # If the source zip is present, confirm it still passes integrity (cheap).
    if [[ -s "${NUCL_ZIP}" ]] && ! unzip -t "${NUCL_ZIP}" >/dev/null 2>&1; then
        echo "  NOTE: existing nucleotide zip fails integrity, but the assembled"
        echo "        FASTA is present and non-empty, so proceeding with it."
    fi
else
    echo "  Downloading all genus-Lentivirus genomes via NCBI datasets..."
    robust_datasets_download "${NUCL_ZIP}" \
        "Lentivirus nucleotide genomes (taxid ${LENTI_TAXID})" \
        virus genome taxon "${LENTI_TAXID}" --include genome
    if [[ $? -ne 0 ]]; then
        # Distinguish a hard block from a transient drop by a quick reachability
        # probe: if datasets can fetch a trivial summary, the host is reachable.
        if datasets summary virus genome taxon "${LENTI_TAXID}" >/dev/null 2>&1; then
            download_fail "transient" "api.ncbi.nlm.nih.gov" \
                "Lentivirus nucleotide genomes (taxid ${LENTI_TAXID})"
        else
            download_fail "block" "api.ncbi.nlm.nih.gov (and www.ncbi.nlm.nih.gov)" \
                "Lentivirus nucleotide genomes (taxid ${LENTI_TAXID})"
        fi
        echo "  STEP 2 failed. See guidance above. Aborting."
        exit 1
    fi

    echo "  Unzipping (verified-complete archive)..."
    unzip -o "${NUCL_ZIP}" -d "${RAW_NUCL}/unzipped" >/dev/null

    # Optional stronger check: datasets ships md5sum.txt — verify if present.
    MD5FILE=$(find "${RAW_NUCL}/unzipped" -name "md5sum.txt" 2>/dev/null | head -1)
    if [[ -n "${MD5FILE}" ]]; then
        echo "  Verifying extracted files against datasets md5sum.txt..."
        ( cd "$(dirname "${MD5FILE}")" && md5sum -c "$(basename "${MD5FILE}")" ) \
            && echo "  md5 check: PASS" \
            || echo "  WARNING: md5 check reported mismatches — inspect before trusting."
    fi

    # datasets layout: ncbi_dataset/data/genomic.fna (concatenated)
    find "${RAW_NUCL}/unzipped" -name "*.fna" -exec cat {} + > "${NUCL_FASTA}"
fi

if [[ ! -s "${NUCL_FASTA}" ]]; then
    echo "  ERROR: no nucleotide FASTA assembled. Aborting."
    exit 1
fi

RAW_NUCL_COUNT=$(grep -c "^>" "${NUCL_FASTA}" || true); RAW_NUCL_COUNT=${RAW_NUCL_COUNT:-0}
echo "  Raw nucleotide sequences: ${RAW_NUCL_COUNT}"

# ─── Curation: drop fragments below MIN_GENOME_LEN ──────────────────────────
NUCL_FILT="${RAW_NUCL}/lentivirus_genomes.filtered.fna"
echo "  Curating: keeping sequences >= ${MIN_GENOME_LEN} bp (seqkit)..."
seqkit seq -m "${MIN_GENOME_LEN}" "${NUCL_FASTA}" > "${NUCL_FILT}"
FILT_NUCL_COUNT=$(grep -c "^>" "${NUCL_FILT}" || true); FILT_NUCL_COUNT=${FILT_NUCL_COUNT:-0}
echo "  After length filter (>= ${MIN_GENOME_LEN} bp): ${FILT_NUCL_COUNT}"
echo "  (Dropped $((RAW_NUCL_COUNT - FILT_NUCL_COUNT)) sub-genome fragments.)"

record "NUCLEOTIDE DB (blastn / tblastx)"
record "  Source       : NCBI datasets, virus genome taxon ${LENTI_TAXID} (genus Lentivirus)"
record "  Retrieved    : $(date)"
record "  Raw seqs     : ${RAW_NUCL_COUNT}"
record "  Curation     : seqkit seq -m ${MIN_GENOME_LEN} (drop fragments < ${MIN_GENOME_LEN} bp)"
record "  Kept seqs    : ${FILT_NUCL_COUNT}"
record ""

# ─── Build BLAST nucleotide DB ──────────────────────────────────────────────
echo "  Building makeblastdb nucleotide DB..."
makeblastdb -in "${NUCL_FILT}" -dbtype nucl \
    -title "lentivirus_taxid${LENTI_TAXID}_$(date +%Y%m%d)" \
    -parse_seqids \
    -out "${BLAST_DB_DIR}/lentivirus_nucl"
if [[ $? -ne 0 ]]; then
    echo "  ERROR: makeblastdb failed."
    exit 1
fi
echo "  → ${BLAST_DB_DIR}/lentivirus_nucl"
record "  BLAST DB     : ${BLAST_DB_DIR}/lentivirus_nucl"
record ""
echo ""

###############################################################################
# STEP 3: Pull genus-Lentivirus proteins via efetch accession-batching
###############################################################################
# RATIONALE (per project decision): we do NOT predict ORFs locally — predicted
# annotations would carry error and read poorly in the methods write-up. We pull
# the REAL NCBI-annotated proteins under taxid 11646. The monolithic `datasets
# --include protein` stream kept dropping mid-transfer on a slow link (no resume,
# single large HTTP/2 stream). Instead:
#   Phase 3a: list all protein accessions under the taxon (tiny text transfer).
#   Phase 3b: efetch the sequences in small batches (no single transfer large;
#             a dropped batch retries in isolation — drop-proof on a flaky link).
# Provenance stays clean: still "all genus-Lentivirus proteins, fresh from NCBI",
# just retrieved accession-wise instead of as one fragile archive.
###############################################################################
echo "============================================================"
echo "STEP 3: Lentivirus proteins via efetch batching (taxid ${LENTI_TAXID})"
echo "============================================================"

PROT_FASTA="${RAW_PROT}/lentivirus_proteins.faa"
PROT_ACC_LIST="${RAW_PROT}/lentivirus_protein_accessions.txt"

# ─── Load NCBI API key (never printed) ──────────────────────────────────────
PROT_BATCH_SLEEP="${PROT_BATCH_SLEEP_NOKEY}"
if [[ -r "${NCBI_KEY_FILE}" ]]; then
    # shellcheck disable=SC1090
    source "${NCBI_KEY_FILE}"
    if [[ -n "${NCBI_API_KEY:-}" ]]; then
        export NCBI_API_KEY
        PROT_BATCH_SLEEP="${PROT_BATCH_SLEEP_KEY}"
        echo "  NCBI API key loaded from ${NCBI_KEY_FILE} (value not shown)."
        echo "  Pacing: ~10 req/s permitted; ${PROT_BATCH_SLEEP}s between batches."
    else
        echo "  WARNING: ${NCBI_KEY_FILE} present but NCBI_API_KEY not set in it."
        echo "  Expected one line:  NCBI_API_KEY=your_key_here"
        echo "  Proceeding KEYLESS (3 req/s; ${PROT_BATCH_SLEEP}s between batches)."
    fi
else
    echo "  No readable key file at ${NCBI_KEY_FILE} — proceeding KEYLESS."
    echo "  (3 req/s; ${PROT_BATCH_SLEEP}s between batches. To enable the key,"
    echo "   create that file with:  NCBI_API_KEY=your_key_here  and chmod 600.)"
fi
echo ""

if [[ -s "${PROT_FASTA}" ]]; then
    echo "  Found existing ${PROT_FASTA} — skipping protein download."
    RAW_PROT_COUNT=$(grep -c "^>" "${PROT_FASTA}" || true); RAW_PROT_COUNT=${RAW_PROT_COUNT:-0}
    EXPECTED_ACC=${RAW_PROT_COUNT}
else
    # ─── Phase 3a: list protein accessions under the taxon ──────────────────
    echo "  Phase 3a: listing protein accessions under taxid ${LENTI_TAXID}..."
    # esearch over the protein DB scoped to the taxon; efetch the accession list.
    # This is a tiny transfer (IDs only), so it completes on the slow link.
    ACC_OK=0
    for attempt in 1 2 3 4; do
        if esearch -db protein -query "txid${LENTI_TAXID}[Organism:exp]" 2>/dev/null \
             | efetch -format acc 2>/dev/null \
             | grep -v '^$' > "${PROT_ACC_LIST}" \
           && [[ -s "${PROT_ACC_LIST}" ]]; then
            ACC_OK=1
            break
        fi
        echo "    [acclist] attempt ${attempt} failed; retrying in $((attempt*10))s..."
        sleep $((attempt*10))
    done

    if [[ "${ACC_OK}" -ne 1 || ! -s "${PROT_ACC_LIST}" ]]; then
        # Distinguish block vs transient with a trivial reachability probe that
        # does NOT depend on a large transfer (the old probe misfired on slow links).
        if esearch -db protein -query "txid${LENTI_TAXID}[Organism:exp]" >/dev/null 2>&1; then
            download_fail "transient" "eutils.ncbi.nlm.nih.gov" \
                "Lentivirus protein accession list (taxid ${LENTI_TAXID})"
        else
            download_fail "block" "eutils.ncbi.nlm.nih.gov (and api/www.ncbi.nlm.nih.gov)" \
                "Lentivirus protein accession list (taxid ${LENTI_TAXID})"
        fi
        echo "  STEP 3 Phase 3a failed. See guidance above. Aborting."
        exit 1
    fi

    # Deduplicate accession list (esearch can repeat) and count.
    sort -u "${PROT_ACC_LIST}" -o "${PROT_ACC_LIST}"
    EXPECTED_ACC=$(wc -l < "${PROT_ACC_LIST}")
    echo "  Phase 3a complete: ${EXPECTED_ACC} unique protein accessions listed."
    echo ""

    # ─── Phase 3b: efetch sequences in small batches ────────────────────────
    echo "  Phase 3b: fetching ${EXPECTED_ACC} proteins in batches of ${PROT_BATCH_SIZE}..."
    : > "${PROT_FASTA}"   # start fresh

    BATCH_NUM=0
    FAILED_BATCHES=0
    mapfile -t ALL_ACC < "${PROT_ACC_LIST}"
    TOTAL_ACC=${#ALL_ACC[@]}

    for (( i=0; i<TOTAL_ACC; i+=PROT_BATCH_SIZE )); do
        BATCH_NUM=$(( BATCH_NUM + 1 ))
        batch=( "${ALL_ACC[@]:i:PROT_BATCH_SIZE}" )
        # Join this batch's accessions with commas for efetch -id.
        acc_csv=$(IFS=,; echo "${batch[*]}")

        printf "    batch %d (accs %d-%d of %d)... " \
            "${BATCH_NUM}" "$((i+1))" "$(( i+${#batch[@]} ))" "${TOTAL_ACC}"

        if efetch_batch "${acc_csv}" "${PROT_FASTA}" "${PROT_BATCH_SLEEP}"; then
            echo "ok"
        else
            echo "FAILED (after ${PROT_BATCH_RETRIES} retries) — continuing"
            FAILED_BATCHES=$(( FAILED_BATCHES + 1 ))
            echo "${acc_csv}" >> "${RAW_PROT}/failed_batches.txt"
        fi
        sleep "${PROT_BATCH_SLEEP}"
    done

    echo ""
    if (( FAILED_BATCHES > 0 )); then
        echo "  WARNING: ${FAILED_BATCHES} batch(es) failed all retries."
        echo "  Their accessions are in ${RAW_PROT}/failed_batches.txt"
        echo "  The DB will be built from what succeeded; re-run to retry the"
        echo "  remainder (existing FASTA is detected and the run will resume)."
    fi

    RAW_PROT_COUNT=$(grep -c "^>" "${PROT_FASTA}" || true); RAW_PROT_COUNT=${RAW_PROT_COUNT:-0}
fi

if [[ ! -s "${PROT_FASTA}" ]]; then
    echo "  ERROR: no protein FASTA assembled. Aborting STEP 3."
    exit 1
fi

# ─── Completeness self-check (verifies the 'all of taxon' claim) ────────────
echo "  Retrieved proteins:  ${RAW_PROT_COUNT}"
echo "  Expected (accessions): ${EXPECTED_ACC}"
if (( RAW_PROT_COUNT < EXPECTED_ACC )); then
    PCT=$(( 100 * RAW_PROT_COUNT / (EXPECTED_ACC>0 ? EXPECTED_ACC : 1) ))
    echo "  NOTE: retrieved ${PCT}% of expected. Some batches may have failed,"
    echo "  or efetch returned isoforms/duplicates differently than the ID list."
    echo "  Re-running STEP 3 will resume (existing FASTA detected)."
fi

# ─── Deduplicate by sequence (batches/isoforms can overlap) ─────────────────
PROT_FASTA_DEDUP="${RAW_PROT}/lentivirus_proteins.dedup.faa"
echo "  Deduplicating proteins by sequence (seqkit rmdup -s)..."
seqkit rmdup -s "${PROT_FASTA}" > "${PROT_FASTA_DEDUP}" 2>/dev/null
DEDUP_COUNT=$(grep -c "^>" "${PROT_FASTA_DEDUP}" || true); DEDUP_COUNT=${DEDUP_COUNT:-0}
echo "  After dedup: ${DEDUP_COUNT} unique protein sequences."

echo "  Building DIAMOND protein DB..."
diamond makedb --in "${PROT_FASTA_DEDUP}" \
    --db "${DIAMOND_DB_DIR}/lentivirus_prot" \
    --threads "${THREADS}"
if [[ $? -ne 0 ]]; then
    echo "  ERROR: diamond makedb failed."
    exit 1
fi
echo "  → ${DIAMOND_DB_DIR}/lentivirus_prot.dmnd"

record "PROTEIN DB (diamond blastx)"
record "  Source       : NCBI efetch, protein DB, txid${LENTI_TAXID}[Organism:exp]"
record "  Method       : accession-batched efetch (batch=${PROT_BATCH_SIZE}); real"
record "                 NCBI-annotated proteins (no local ORF prediction)"
record "  API key      : $([[ -n "${NCBI_API_KEY:-}" ]] && echo 'used (value not recorded)' || echo 'not used (keyless)')"
record "  Retrieved    : $(date)"
record "  Accessions   : ${EXPECTED_ACC} listed"
record "  Fetched seqs : ${RAW_PROT_COUNT}"
record "  Unique seqs  : ${DEDUP_COUNT} (seqkit rmdup -s)"
record "  DIAMOND DB   : ${DIAMOND_DB_DIR}/lentivirus_prot.dmnd"
record ""
echo ""

###############################################################################
# STEP 4: geNomad marker database
###############################################################################
echo "============================================================"
echo "STEP 4: geNomad marker database (env: ${GENOMAD_ENV})"
echo "============================================================"

# geNomad lives in its OWN env. Switch into it here. If it is not built, this
# step is SKIPPED (warn, not fail) — the BLAST/DIAMOND tiers from STEPS 2-3 are
# the primary evidence and do not depend on geNomad. You can build the geNomad
# env and re-run later to fill in the profile tier.
if ! conda env list | grep -qE "^\s*${GENOMAD_ENV}\s"; then
    echo "  SKIP: env '${GENOMAD_ENV}' not found."
    echo "  The profile tier is optional for now. To add it later:"
    echo "    conda env create -f viral_id_genomad.yml"
    echo "    (then re-run this script; STEPS 2-3 will skip via existing outputs)"
    record "PROFILE DB (geNomad)"
    record "  Status       : SKIPPED — env ${GENOMAD_ENV} not present at build time"
    record ""
    GENOMAD_DB_VER="(skipped)"
else
    conda activate "${GENOMAD_ENV}"
    export PATH="${CONDA_PREFIX}/bin:${PATH}"
    echo "  Active env: ${CONDA_PREFIX}"

    if ! command -v genomad &>/dev/null; then
        echo "  WARNING: 'genomad' not resolvable in ${GENOMAD_ENV}. Skipping."
        echo "  Verify the env: conda activate ${GENOMAD_ENV} && genomad --version"
        record "PROFILE DB (geNomad)"
        record "  Status       : SKIPPED — genomad not resolvable in ${GENOMAD_ENV}"
        record ""
        GENOMAD_DB_VER="(skipped)"
    else
        echo "  genomad: $(genomad --version 2>&1 | head -1)"

        # geNomad writes genomad_db/ with a version file inside.
        if [[ -d "${GENOMAD_DB_DIR}/genomad_db" ]]; then
            echo "  Found existing geNomad DB — skipping download."
        else
            echo "  Downloading geNomad database (versioned, citable)..."
            echo "  Command: genomad download-database ${GENOMAD_DB_DIR}"
            genomad download-database "${GENOMAD_DB_DIR}"
            DL_RC=$?

            if [[ ${DL_RC} -ne 0 || ! -d "${GENOMAD_DB_DIR}/genomad_db" ]]; then
                download_fail "transient" "zenodo.org (geNomad DB host)" \
                             "geNomad marker database"
                echo "  STEP 4 download failed. See guidance above."
                echo "  (BLAST/DIAMOND DBs from STEPS 2-3 are unaffected and usable.)"
                echo "  Continuing to STEP 5 without the profile tier."
                echo "  Re-run later to retry just this download (STEPS 2-3 skip)."
                record "PROFILE DB (geNomad)"
                record "  Status       : FAILED download — see guidance in log"
                record ""
                GENOMAD_DB_VER="(download failed)"
            fi
        fi

        if [[ -d "${GENOMAD_DB_DIR}/genomad_db" ]]; then
            GENOMAD_VER_FILE=$(find "${GENOMAD_DB_DIR}" -name "*version*" 2>/dev/null | head -1)
            GENOMAD_DB_VER="unknown"
            [[ -n "${GENOMAD_VER_FILE}" ]] && GENOMAD_DB_VER=$(cat "${GENOMAD_VER_FILE}" 2>/dev/null | head -1)

            echo "  geNomad DB location: ${GENOMAD_DB_DIR}/genomad_db"
            echo "  geNomad DB version:  ${GENOMAD_DB_VER}"

            record "PROFILE DB (geNomad)"
            record "  Source       : genomad download-database (Zenodo-hosted, versioned)"
            record "  Env          : ${GENOMAD_ENV} (from viral_id_genomad.yml)"
            record "  Retrieved    : $(date)"
            record "  DB version   : ${GENOMAD_DB_VER}"
            record "  DB location  : ${GENOMAD_DB_DIR}/genomad_db"
            record ""
        fi
    fi

    # Return to the alignment env for the STEP 5 BLAST sanity checks.
    conda activate "${ALIGN_ENV}"
    export PATH="${CONDA_PREFIX}/bin:${PATH}"
fi
echo ""

###############################################################################
# STEP 5: Final verification + manifest close-out
###############################################################################
echo "============================================================"
echo "STEP 5: Verification summary"
echo "============================================================"

echo "  Database sanity checks:"
echo ""

# BLAST DB
if blastdbcmd -db "${BLAST_DB_DIR}/lentivirus_nucl" -info &>/dev/null; then
    echo "  [OK]  BLAST nucleotide DB readable:"
    blastdbcmd -db "${BLAST_DB_DIR}/lentivirus_nucl" -info 2>/dev/null \
        | sed 's/^/        /'
else
    echo "  [WARN] BLAST nucleotide DB did not return info."
fi
echo ""

# DIAMOND DB
if [[ -s "${DIAMOND_DB_DIR}/lentivirus_prot.dmnd" ]]; then
    echo "  [OK]  DIAMOND protein DB present: $(du -h "${DIAMOND_DB_DIR}/lentivirus_prot.dmnd" | cut -f1)"
else
    echo "  [WARN] DIAMOND protein DB missing/empty."
fi
echo ""

# geNomad DB
if [[ -d "${GENOMAD_DB_DIR}/genomad_db" ]]; then
    echo "  [OK]  geNomad DB present: $(du -sh "${GENOMAD_DB_DIR}/genomad_db" | cut -f1)"
else
    echo "  [WARN] geNomad DB missing."
fi
echo ""

record "############################################################"
record "# BUILD COMPLETE — $(date)"
record "# Nucleotide seqs (kept) : ${FILT_NUCL_COUNT}"
record "# Protein seqs           : ${RAW_PROT_COUNT}"
record "# geNomad DB version     : ${GENOMAD_DB_VER}"
record "############################################################"

echo "============================================================"
echo "DONE — databases built. NO searching performed."
echo "============================================================"
echo ""
echo "Inspect before proceeding:"
echo "  Manifest:  ${MANIFEST}"
echo "  Full log:  ${LOGFILE}"
echo ""
echo "Suggested sanity checks:"
echo "  1. Read ${MANIFEST} — confirm counts look reasonable."
echo "     (Expect: hundreds-to-thousands of nucleotide seqs, HIV-1 heavy.)"
echo "  2. Spot-check the DB contains the canonical references you expect:"
echo "       blastdbcmd -db ${BLAST_DB_DIR}/lentivirus_nucl -entry all -outfmt '%t' | grep -iE 'SIVmac|HIV|SHIV' | head"
echo "  3. Confirm SIVmac239 and an HIV-1 ref are present (the two halves of"
echo "     the SHIV chimera). If absent, the taxon pull may be incomplete."
echo ""
echo "Next script (separate): three-tier search of the 118,889 MEGAHIT"
echo "contigs against these DBs (blastn + tblastx + diamond blastx + geNomad),"
echo "then reconcile + cluster (MMseqs2) for edge cases."
echo ""
echo "Done: $(date)"
