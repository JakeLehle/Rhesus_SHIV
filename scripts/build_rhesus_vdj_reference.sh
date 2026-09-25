#!/bin/bash
# ============================================================
# build_rhesus_vdj_reference.sh
#
# Builds a Cell Ranger VDJ reference for Rhesus Macaque
# (Macaca mulatta) from IMGT, then validates Zoey's custom
# TCR enrichment primers against the reference.
#
# Based on:
#   - Walsh et al. 2022 (J Immunol 208:762-771)
#   - aherman-umn/10x_Mmul guide
#   - 10x Genomics custom VDJ reference documentation
#
# Requirements:
#   - Cell Ranger v9.0.1
#   - BLAST+ (makeblastdb, blastn)
#   - Internet access (to download from IMGT)
#
# Usage:
#   bash build_rhesus_vdj_reference.sh
#
# Author: Jake Lehle
# Date: April 2026
# ============================================================

set -eo pipefail

# ============================================================
# CONFIGURATION — UPDATE THESE PATHS
# ============================================================

# Cell Ranger installation path
CELLRANGER_DIR="/master/jlehle/cellranger-9.0.1"

# Where to build the reference
WORK_DIR="/master/jlehle/WORKING/SC/fastq/FLEX/vdj_reference"

# Reference name
REF_NAME="vdj_rhesus_macaque_IMGT"

# Species for IMGT
SPECIES="Macaca mulatta"

# ============================================================
# ZOEY'S CUSTOM INNER ENRICHMENT PRIMERS
# These target the constant regions of all 4 rhesus TCR chains.
# Inner primers are what Cell Ranger checks against the reference.
# ============================================================

declare -A INNER_PRIMERS
INNER_PRIMERS[TCRA]="ATGCACGTCAGAATCCTTGCT"
INNER_PRIMERS[TCRB]="TCTGATGGCTCAAACACAGC"
INNER_PRIMERS[TCRD]="GGGAGAGACGACAATAGCAGGA"
INNER_PRIMERS[TCRG]="TGGAGGTTTGTTTCAGCAATGGA"

# ============================================================
# SETUP
# ============================================================

echo "============================================================"
echo " Rhesus Macaque VDJ Reference Builder"
echo " Cell Ranger v9.0.1"
echo "============================================================"
echo ""

mkdir -p "${WORK_DIR}"
cd "${WORK_DIR}"

# Source Cell Ranger environment
if [[ -f "${CELLRANGER_DIR}/sourceme.bash" ]]; then
    source "${CELLRANGER_DIR}/sourceme.bash"
    echo "[OK] Sourced Cell Ranger environment"
else
    echo "[WARN] sourceme.bash not found at ${CELLRANGER_DIR}"
    echo "       Making sure cellranger is on PATH..."
    export PATH="${CELLRANGER_DIR}:${PATH}"
fi

echo "[INFO] Cell Ranger version:"
cellranger --version
echo ""

# ============================================================
# STEP 1: PATCH fetch-imgt FOR RHESUS C-REGIONS
#
# The fetch-imgt script downloads V(D)J segments from IMGT.
# For rhesus, the C-region genes need a different IMGT version
# query (7.2 instead of 14.1) to pull TR and IG heavy chains.
#
# We copy the script and patch it rather than editing in-place.
# ============================================================

echo "============================================================"
echo " Step 1: Preparing fetch-imgt for Macaca mulatta"
echo "============================================================"

FETCH_IMGT_ORIG="${CELLRANGER_DIR}/lib/bin/fetch-imgt"

if [[ ! -f "${FETCH_IMGT_ORIG}" ]]; then
    echo "[ERROR] fetch-imgt not found at ${FETCH_IMGT_ORIG}"
    echo "        Check your Cell Ranger installation path."
    exit 1
fi

# Copy to working dir so we can patch safely
cp "${FETCH_IMGT_ORIG}" "${WORK_DIR}/fetch-imgt-rhesus"
chmod +x "${WORK_DIR}/fetch-imgt-rhesus"

# Patch the C-region query version for rhesus
# The original script has c_query = "14.1" which only works for
# IG kappa/lambda. Changing to "7.2" pulls TR and IG heavy C-regions.
if grep -q 'c_query = "14.1"' "${WORK_DIR}/fetch-imgt-rhesus"; then
    sed -i 's/c_query = "14.1"/c_query = "7.2"/' "${WORK_DIR}/fetch-imgt-rhesus"
    echo "[OK] Patched c_query from 14.1 to 7.2 for rhesus C-regions"
elif grep -q "c_query" "${WORK_DIR}/fetch-imgt-rhesus"; then
    echo "[WARN] c_query found but not matching expected pattern."
    echo "       You may need to manually edit fetch-imgt-rhesus"
    echo "       Change the c_query value to \"7.2\" for the C-region genes:"
    echo "         c_genes = [\"TRAC\", \"TRBC\", \"TRDC\", \"TRGC\", \"IGHC\"]"
    echo "         c_query = \"7.2\""
else
    echo "[WARN] c_query not found in fetch-imgt. The script format may"
    echo "       have changed in v9.0.1. Check manually and edit if needed."
    echo "       For rhesus, C-region genes need IMGT version 7.2"
fi
echo ""

# ============================================================
# STEP 2: FETCH IMGT SEQUENCES FOR MACACA MULATTA
# ============================================================

echo "============================================================"
echo " Step 2: Fetching IMGT sequences for ${SPECIES}"
echo "============================================================"

python3 "${WORK_DIR}/fetch-imgt-rhesus" \
    --genome "${REF_NAME}" \
    --species "${SPECIES}"

IMGT_FASTA="${WORK_DIR}/${REF_NAME}-mkvdjref-input.fasta"

if [[ ! -f "${IMGT_FASTA}" ]]; then
    echo "[ERROR] Expected output not found: ${IMGT_FASTA}"
    echo "        fetch-imgt may have failed. Check output above."
    exit 1
fi

# Count what we got
N_SEQS=$(grep -c "^>" "${IMGT_FASTA}")
echo ""
echo "[OK] Downloaded ${N_SEQS} V(D)J segment sequences"
echo ""

# Print summary of chain types
echo "Sequence breakdown by chain type:"
for CHAIN in TRAV TRBV TRDV TRGV TRAJ TRBJ TRDJ TRGJ TRBD TRAC TRBC TRDC TRGC \
             IGHV IGKV IGLV IGHJ IGKJ IGLJ IGHD IGHC IGKC IGLC; do
    COUNT=$(grep -c "${CHAIN}" "${IMGT_FASTA}" 2>/dev/null || echo 0)
    if [[ ${COUNT} -gt 0 ]]; then
        echo "  ${CHAIN}: ${COUNT}"
    fi
done
echo ""

# ============================================================
# STEP 3: BUILD THE VDJ REFERENCE WITH cellranger mkvdjref
# ============================================================

echo "============================================================"
echo " Step 3: Building VDJ reference with cellranger mkvdjref"
echo "============================================================"

# Remove existing reference dir if present (mkvdjref won't overwrite)
if [[ -d "${WORK_DIR}/${REF_NAME}" ]]; then
    echo "[WARN] Removing existing reference directory: ${REF_NAME}"
    rm -rf "${WORK_DIR}/${REF_NAME}"
fi

cellranger mkvdjref \
    --genome="${REF_NAME}" \
    --seqs="${IMGT_FASTA}"

if [[ -f "${WORK_DIR}/${REF_NAME}/fasta/regions.fa" ]]; then
    echo ""
    echo "[OK] VDJ reference built successfully at:"
    echo "     ${WORK_DIR}/${REF_NAME}/"
    echo ""
    echo "Reference contents:"
    tree "${WORK_DIR}/${REF_NAME}/" 2>/dev/null || ls -R "${WORK_DIR}/${REF_NAME}/"
else
    echo "[ERROR] Reference build failed — regions.fa not found"
    exit 1
fi
echo ""

# ============================================================
# STEP 4: VALIDATE INNER ENRICHMENT PRIMERS AGAINST REFERENCE
#
# Cell Ranger does exact string matching of inner primers against
# C-region sequences. If there's even one SNP, preflight fails.
# We BLAST the primers to check for mismatches and report any
# corrections needed.
# ============================================================

echo "============================================================"
echo " Step 4: Validating inner enrichment primers"
echo "============================================================"

REGIONS_FA="${WORK_DIR}/${REF_NAME}/fasta/regions.fa"

# Check if BLAST is available
if ! command -v makeblastdb &>/dev/null; then
    echo "[WARN] BLAST+ not found on PATH. Skipping primer validation."
    echo "       Install with: conda install -c bioconda blast"
    echo "       Then re-run this step manually."
    echo ""
    echo "       Alternatively, you can manually check primers by searching"
    echo "       for exact matches in: ${REGIONS_FA}"
else
    # Build BLAST database from reference
    BLAST_DB="${WORK_DIR}/vdj_blast_db"
    makeblastdb -in "${REGIONS_FA}" -dbtype nucl -out "${BLAST_DB}" -title "Rhesus VDJ" 2>/dev/null
    echo "[OK] BLAST database created"
    echo ""

    # Write primers to FASTA
    PRIMER_FASTA="${WORK_DIR}/inner_primers.fa"
    > "${PRIMER_FASTA}"
    for CHAIN in TCRA TCRB TCRD TCRG; do
        echo ">${CHAIN}_inner" >> "${PRIMER_FASTA}"
        echo "${INNER_PRIMERS[${CHAIN}]}" >> "${PRIMER_FASTA}"
    done

    # BLAST primers against reference
    echo "BLAST results for each inner enrichment primer:"
    echo "-----------------------------------------------"
    echo ""

    PRIMER_RESULTS="${WORK_DIR}/primer_validation_results.txt"
    > "${PRIMER_RESULTS}"

    ALL_PASS=true

    for CHAIN in TCRA TCRB TCRD TCRG; do
        PRIMER_SEQ="${INNER_PRIMERS[${CHAIN}]}"
        PRIMER_LEN=${#PRIMER_SEQ}

        echo "--- ${CHAIN} inner primer (${PRIMER_LEN}bp): ${PRIMER_SEQ} ---"

        # Write single primer to temp file
        TEMP_FA=$(mktemp)
        echo ">${CHAIN}_inner" > "${TEMP_FA}"
        echo "${PRIMER_SEQ}" >> "${TEMP_FA}"

        # BLAST with relaxed parameters for short sequences
        HITS=$(blastn -query "${TEMP_FA}" -db "${BLAST_DB}" \
            -word_size 7 -evalue 10 -dust no \
            -outfmt "6 sseqid pident length mismatch qseq sseq" \
            -max_target_seqs 5 2>/dev/null | head -5)

        if [[ -z "${HITS}" ]]; then
            echo "  [WARN] No BLAST hits found!"
            echo "         This primer may not match any C-region in the reference."
            echo "         Cell Ranger will likely fail preflight for this chain."
            ALL_PASS=false
            echo "${CHAIN}: NO HITS" >> "${PRIMER_RESULTS}"
        else
            echo "${HITS}" | while IFS=$'\t' read -r SUBJ PIDENT ALEN MISMATCH QSEQ SSEQ; do
                if [[ "${MISMATCH}" == "0" ]]; then
                    echo "  [OK] Perfect match to ${SUBJ} (${PIDENT}% identity, ${ALEN}bp aligned)"
                    echo "${CHAIN}: EXACT MATCH to ${SUBJ}" >> "${PRIMER_RESULTS}"
                else
                    echo "  [WARN] ${MISMATCH} mismatch(es) to ${SUBJ} (${PIDENT}% identity)"
                    echo "         Query:   ${QSEQ}"
                    echo "         Subject: ${SSEQ}"
                    echo "         Use the Subject sequence as the corrected primer."
                    echo "${CHAIN}: ${MISMATCH} MISMATCH(ES) to ${SUBJ}" >> "${PRIMER_RESULTS}"
                    echo "  Original: ${QSEQ}" >> "${PRIMER_RESULTS}"
                    echo "  Corrected: ${SSEQ}" >> "${PRIMER_RESULTS}"
                fi
            done
        fi
        echo ""
        rm -f "${TEMP_FA}"
    done

    echo "Validation summary saved to: ${PRIMER_RESULTS}"
    echo ""

    if ${ALL_PASS}; then
        echo "[OK] All primers have hits in the reference."
    else
        echo "[ACTION NEEDED] Some primers had no hits or mismatches."
        echo "  Check ${PRIMER_RESULTS} for details."
    fi
fi
echo ""

# ============================================================
# STEP 5: GENERATE PRIMER FILE FOR CELLRANGER MULTI
#
# Create the inner enrichment primer CSV that cellranger multi
# needs. If any primers needed correction from Step 4, update
# the sequences below.
# ============================================================

echo "============================================================"
echo " Step 5: Generating primer file for cellranger multi"
echo "============================================================"

PRIMER_CSV="${WORK_DIR}/inner_enrichment_primers.csv"
cat > "${PRIMER_CSV}" << 'PRIMERS_EOF'
gene_name,name,sequence
TCRA,TCRA_inner,ATGCACGTCAGAATCCTTGCT
TCRB,TCRB_inner,TCTGATGGCTCAAACACAGC
TCRD,TCRD_inner,GGGAGAGACGACAATAGCAGGA
TCRG,TCRG_inner,TGGAGGTTTGTTTCAGCAATGGA
PRIMERS_EOF

echo "[OK] Inner enrichment primer CSV written to:"
echo "     ${PRIMER_CSV}"
echo ""
echo "Contents:"
cat "${PRIMER_CSV}"
echo ""
echo "NOTE: If Step 4 reported mismatches, update the sequences in"
echo "      this file with the corrected (Subject) sequences before"
echo "      running cellranger multi."
echo ""

# ============================================================
# STEP 6: GENERATE TEMPLATE cellranger multi CONFIG
# ============================================================

echo "============================================================"
echo " Step 6: Template cellranger multi config"
echo "============================================================"

MULTI_CONFIG="${WORK_DIR}/multi_config_template.csv"
cat > "${MULTI_CONFIG}" << TEMPLATE_EOF
[gene-expression]
reference,/path/to/Mmul_10_GEX_reference
create-bam,true

[vdj]
reference,${WORK_DIR}/${REF_NAME}
inner-enrichment-primers,${PRIMER_CSV}

[libraries]
fastq_id,fastqs,feature_types
GEX_SAMPLE_ID,/path/to/gex_fastqs,Gene Expression
VDJ_SAMPLE_ID,/path/to/vdj_fastqs,VDJ-T
TEMPLATE_EOF

echo "[OK] Template multi config written to:"
echo "     ${MULTI_CONFIG}"
echo ""
echo "Zoey needs to update this with her actual paths before running:"
echo "  - GEX reference path"
echo "  - FASTQ paths for GEX and VDJ libraries"
echo "  - Sample IDs matching her FASTQ filenames"
echo ""

# ============================================================
# SUMMARY
# ============================================================

echo ""
echo "============================================================"
echo " DONE — Summary"
echo "============================================================"
echo ""
echo "VDJ reference:      ${WORK_DIR}/${REF_NAME}/"
echo "IMGT FASTA:         ${IMGT_FASTA}"
echo "Inner primers CSV:  ${PRIMER_CSV}"
echo "Multi config:       ${MULTI_CONFIG}"
echo "Primer validation:  ${WORK_DIR}/primer_validation_results.txt"
echo ""
echo "Next steps:"
echo "  1. Review primer_validation_results.txt"
echo "     - If mismatches found, update inner_enrichment_primers.csv"
echo "       with the corrected sequences"
echo "  2. Have Zoey update multi_config_template.csv with her paths"
echo "  3. Run:  cellranger multi --id=run_id --csv=multi_config.csv"
echo "  4. VDJ contigs can then be integrated with the annotated"
echo "     AnnData object from the GEX pipeline"
echo ""
echo "References:"
echo "  - Walsh et al. 2022 (doi:10.4049/jimmunol.2100824)"
echo "  - Brochu et al. 2020 (doi:10.4049/jimmunol.1901201)"
echo "  - https://github.com/ncsu-penglab/RhesusIgTCR"
echo "  - https://github.com/aherman-umn/10x_Mmul"
echo ""
