#!/bin/bash
#SBATCH -J CR_MULTI_SHIV
#SBATCH -o /master/jlehle/WORKING/LOGS/cellranger_multi_shiv_%A_%a.o.log
#SBATCH -e /master/jlehle/WORKING/LOGS/cellranger_multi_shiv_%A_%a.e.log
#SBATCH -t 5-00:00:00
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 24
#SBATCH --mem 800GB

###############################################################################
# Cell Ranger multi — Rhesus SHIV 5' scRNA-seq + CITE-seq + VDJ
#
# This script processes all three modalities together per library:
#   - GEX:     Gene expression (5' chemistry) → custom Mmul10 + SHIVAD8-EO ref
#   - PROTEIN: CITE-seq antibody capture → feature_reference.csv
#   - VDJ_RH:  TCR assembly (rhesus-specific primers) → IMGT VDJ reference
#
# CORRECTED Library → Sample Mapping (May 2026)
# ==============================================
# Verified against provider's demultiplexed per_sample_outs and
# Zoey's Group_H_CITEseq_Matrix_Info.xlsx.
#
# Libraries A-H correspond to PHYSICAL pooling for sequencing.
# Some libraries are MIXED-CONDITION pools (C, D, F) where animals
# from different timepoints or subspecies share a GEM well.
# Per-animal condition assignment happens downstream in Scanpy
# using the provider's genetic demultiplexing barcode assignments.
#
#   A = Pre PBMC                     (Chinese: 39272, 40702, 40707, 41861)
#   B = Pre LN                      (Chinese: 39272, 40702, 40707, 41861)
#   C = Pre PBMC [MIXED]            (Indian: 34315, 41903 + Chinese: 41862)
#   D = 21DPI+Pre LN [MIXED]        (21DPI: 39272, 40702, 40707 + Pre: 41862)
#   E = 21DPI PBMC                  (Chinese: 39272, 40702, 40707, 41861)
#   F = 21DPI+NonInf LN [MIXED]     (21DPI: 41861, 41862 + NonInf: 34315, 41903)
#   G = Necropsy PBMC               (Chinese: 40702, 40707, 41861, 41862)
#   H = Necropsy Mesenteric LN      (Chinese: 40702, 40707, 41861, 41862)
#
# Animal notes:
#   - 40707 = elite controller
#   - 34315, 41903 = Indian RM (non-infected controls only)
#   - 39272 = no necropsy samples (possibly sacrificed at 21 DPI)
#   - 41862 = pooled separately from other Chinese RM at Pre/21DPI
#
# NOTE on VDJ libraries:
#   Two VDJ library sets exist in the raw data:
#     - VDJ_LIBRARY_X     (index position 2, e.g. SI_TT_A2)
#     - VDJ_RH_LIBRARY_X  (index position 3, e.g. SI_TT_A3)
#   We use VDJ_RH (rhesus-specific TCR primers) with the custom
#   IMGT-based rhesus VDJ reference and inner enrichment primers.
#
# Usage:
#   # Test with one library first:
#   sbatch --array=0 CELLRANGER_MULTI_RH_SHIV.sh
#
#   # All 8 libraries:
#   sbatch --array=0-7 CELLRANGER_MULTI_RH_SHIV.sh
###############################################################################

set -euo pipefail

#--- Configuration ---#
CRDIR="/master/jlehle/cellranger-10.0.0"
FASTQ_BASE="/master/zwallis/WORKING/SC/10X_RAW"
OUTDIR="/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV/Analysis"
CONFIG_DIR="${OUTDIR}/multi_configs"

# References
GEX_REF="/master/jlehle/WORKING/SC/REF/Mmul10_SHIVAD8EO_v10"
GEX_REF_ALT="/master/jlehle/cellranger-10.0.0/Mmul10_SHIVAD8EO_v10"
VDJ_REF="/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV/vdj_reference/vdj_rhesus_macaque_IMGT"
INNER_PRIMERS="/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV/vdj_reference/inner_enrichment_primers.txt"
FEATURE_REF="/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV/feature_reference.csv"

#--- CORRECTED Library definitions ---#
# Order: A B C D E F G H (index 0-7)
LETTERS=(A B C D E F G H)
CONDITIONS=(
    "Pre_PBMC"
    "Pre_LN"
    "Pre_PBMC_IndianMixed"
    "21DPI_Pre_LN_Mixed"
    "21DPI_PBMC"
    "21DPI_NonInf_LN_Mixed"
    "Necropsy_PBMC"
    "Necropsy_LN"
)

# FASTQ sample prefixes (extracted from filenames before _S*_L*)
GEX_PREFIXES=(
    "GEX_LIBRARY_A-SI_TT_A1_23552KLT4"
    "GEX_LIBRARY_B-SI_TT_B1_23552KLT4"
    "GEX_LIBRARY_C-SI_TT_C1_23552KLT4"
    "GEX_LIBRARY_D-SI_TT_D1_23552KLT4"
    "GEX_LIBRARY_E-SI_TT_E1_23552KLT4"
    "GEX_LIBRARY_F-SI_TT_F1_23552KLT4"
    "GEX_LIBRARY_G-SI_TT_G1_23552KLT4"
    "GEX_LIBRARY_H-SI_TT_H1_23552KLT4"
)
PROTEIN_PREFIXES=(
    "PROTEIN_LIBRARY_A-SI_TN_A1_23552KLT4"
    "PROTEIN_LIBRARY_B-SI_TN_B1_23552KLT4"
    "PROTEIN_LIBRARY_C-SI_TN_C1_23552KLT4"
    "PROTEIN_LIBRARY_D-SI_TN_D1_23552KLT4"
    "PROTEIN_LIBRARY_E-SI_TN_E1_23552KLT4"
    "PROTEIN_LIBRARY_F-SI_TN_F1_23552KLT4"
    "PROTEIN_LIBRARY_G-SI_TN_G1_23552KLT4"
    "PROTEIN_LIBRARY_H-SI_TN_H1_23552KLT4"
)
VDJ_RH_PREFIXES=(
    "VDJ_RH_LIBRARY_A-SI_TT_A3_23552KLT4"
    "VDJ_RH_LIBRARY_B-SI_TT_B3_23552KLT4"
    "VDJ_RH_LIBRARY_C-SI_TT_C3_23552KLT4"
    "VDJ_RH_LIBRARY_D-SI_TT_D3_23552KLT4"
    "VDJ_RH_LIBRARY_E-SI_TT_E3_23552KLT4"
    "VDJ_RH_LIBRARY_F-SI_TT_F3_23552KLT4"
    "VDJ_RH_LIBRARY_G-SI_TT_G3_23552KLT4"
    "VDJ_RH_LIBRARY_H-SI_TT_H3_23552KLT4"
)

###############################################################################
# Preflight checks
###############################################################################

echo "=== Cell Ranger multi — Rhesus SHIV (CORRECTED MAPPING) ==="
echo ""

# Check GEX reference
if [ ! -d "${GEX_REF}" ]; then
    if [ -d "${GEX_REF_ALT}" ]; then
        echo "[INFO] GEX reference found in cellranger dir, not REF/."
        echo "       Using: ${GEX_REF_ALT}"
        GEX_REF="${GEX_REF_ALT}"
    else
        echo "ERROR: GEX reference not found at either location:"
        echo "  ${GEX_REF}"
        echo "  ${GEX_REF_ALT}"
        exit 1
    fi
fi
echo "[OK] GEX reference: ${GEX_REF}"

# Check VDJ reference
if [ ! -d "${VDJ_REF}" ]; then
    echo "ERROR: VDJ reference not found: ${VDJ_REF}"
    exit 1
fi
echo "[OK] VDJ reference: ${VDJ_REF}"

if [ ! -f "${INNER_PRIMERS}" ]; then
    echo "ERROR: Inner enrichment primers not found: ${INNER_PRIMERS}"
    exit 1
fi
echo "[OK] Inner primers: ${INNER_PRIMERS}"

# Check feature reference
if [ -z "${FEATURE_REF}" ]; then
    echo "[WARN] FEATURE_REF not set — skipping CITE-seq."
    INCLUDE_PROTEIN=false
elif [ ! -f "${FEATURE_REF}" ]; then
    echo "ERROR: Feature reference not found: ${FEATURE_REF}"
    exit 1
else
    echo "[OK] Feature reference: ${FEATURE_REF}"
    INCLUDE_PROTEIN=true
fi

# Check FASTQ directories
IDX=${SLURM_ARRAY_TASK_ID:-0}
LTR="${LETTERS[$IDX]}"

GEX_FASTQ_DIR="${FASTQ_BASE}/GEX_LIBRARY_${LTR}"
VDJ_FASTQ_DIR="${FASTQ_BASE}/VDJ_RH_LIBRARY_${LTR}"
PROTEIN_FASTQ_DIR="${FASTQ_BASE}/PROTEIN_LIBRARY_${LTR}"

for FDIR in "${GEX_FASTQ_DIR}" "${VDJ_FASTQ_DIR}"; do
    if [ ! -d "${FDIR}" ]; then
        echo "ERROR: FASTQ directory not found: ${FDIR}"
        exit 1
    fi
done
if ${INCLUDE_PROTEIN} && [ ! -d "${PROTEIN_FASTQ_DIR}" ]; then
    echo "ERROR: PROTEIN FASTQ directory not found: ${PROTEIN_FASTQ_DIR}"
    exit 1
fi

echo ""
echo "Processing library ${LTR}: ${CONDITIONS[$IDX]}"
echo ""

###############################################################################
# Generate multi config CSV
###############################################################################

mkdir -p "${CONFIG_DIR}" "${OUTDIR}"

SAMPLE_ID="SHIV_${CONDITIONS[$IDX]}"
CONFIG_CSV="${CONFIG_DIR}/multi_config_${LTR}_${CONDITIONS[$IDX]}.csv"

echo "Generating config: ${CONFIG_CSV}"

cat > "${CONFIG_CSV}" << CONFIGEOF
[gene-expression]
reference,${GEX_REF}
create-bam,true
expect-cells,10000

[vdj]
reference,${VDJ_REF}
inner-enrichment-primers,${INNER_PRIMERS}

CONFIGEOF

if ${INCLUDE_PROTEIN}; then
    cat >> "${CONFIG_CSV}" << CONFIGEOF2
[feature]
reference,${FEATURE_REF}

CONFIGEOF2
fi

echo "[libraries]" >> "${CONFIG_CSV}"
echo "fastq_id,fastqs,feature_types" >> "${CONFIG_CSV}"
echo "${GEX_PREFIXES[$IDX]},${GEX_FASTQ_DIR},Gene Expression" >> "${CONFIG_CSV}"
echo "${VDJ_RH_PREFIXES[$IDX]},${VDJ_FASTQ_DIR},VDJ-T" >> "${CONFIG_CSV}"

if ${INCLUDE_PROTEIN}; then
    echo "${PROTEIN_PREFIXES[$IDX]},${PROTEIN_FASTQ_DIR},Antibody Capture" >> "${CONFIG_CSV}"
fi

echo ""
echo "--- Config contents ---"
cat "${CONFIG_CSV}"
echo "--- End config ---"
echo ""

###############################################################################
# Run cellranger multi
###############################################################################

cd "${OUTDIR}"

if [ -d "${SAMPLE_ID}" ]; then
    echo "Removing previous output: ${SAMPLE_ID}/"
    rm -rf "${SAMPLE_ID}"
fi

echo "=== Running cellranger multi ==="
echo "Sample ID: ${SAMPLE_ID}"
echo ""

"${CRDIR}/cellranger" multi \
    --id="${SAMPLE_ID}" \
    --csv="${CONFIG_CSV}" \
    --localcores=${SLURM_CPUS_PER_TASK:-24} \
    --localmem=$((${SLURM_MEM_PER_NODE:-204800} / 1024))

echo ""
echo "=== Done: ${SAMPLE_ID} ==="
echo "Output: ${OUTDIR}/${SAMPLE_ID}/outs/"
echo ""
echo "Condition: ${CONDITIONS[$IDX]}"
echo "Library:   ${LTR}"
