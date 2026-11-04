#!/bin/bash
#SBATCH -J CR_MULTI_SHIV_v11
#SBATCH -o /master/jlehle/WORKING/LOGS/cellranger_multi_shiv_v11_%A_%a.o.log
#SBATCH -e /master/jlehle/WORKING/LOGS/cellranger_multi_shiv_v11_%A_%a.e.log
#SBATCH -t 5-00:00:00
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 24
#SBATCH --mem 800GB

###############################################################################
# 03 — Cell Ranger multi against Mmul10_SHIVAD8EO_v11
#
# Replaces CELLRANGER_MULTI_RH_SHIV.sh. Three modalities per library:
#   GEX      5' gene expression  -> Mmul10 + SHIVAD8EO_V4fullN (collapsed GTF)
#   PROTEIN  CITE-seq            -> feature_reference.csv
#   VDJ_RH   rhesus TCR primers  -> IMGT rhesus VDJ reference
#
# WHAT CHANGED, AND WHY
# ---------------------
# 1. chemistry,SC5P-PE-v3 IS NOW SET EXPLICITLY. The old config carried no
#    chemistry line, so Cell Ranger 10.0.0 auto-detected and settled on
#    SC5P-R2-v3, which reads R1 only for the cell barcode and UMI and throws
#    the rest away. SHIV_R1_Structure_Check.py measured R1 as
#    [CB 16][UMI 12][TSO 13] then cDNA, and shiv_read_length_check.py measured
#    R1 at 150bp in every GEX library, so roughly 109 bases of real transcript
#    per read were being discarded. Confirmed absent from both config.csv and
#    metrics_summary.csv on 2026-10-30.
#
# 2. THE REFERENCE IS v11, built by 02. It carries the 34bp cassette masked as
#    N at locus 5862 and ONE collapsed GTF feature spanning the whole contig.
#    v10's per-ORF GTF left 33.4% of the viral contig under two or more
#    features and 11.6% under none.
#
# 3. OUTPUT GOES TO A NEW TREE AND Analysis/ IS NEVER TOUCHED. The old script
#    did `rm -rf "${SAMPLE_ID}"` inside Analysis/ before every run. Doing that
#    now would destroy the v10 baseline the gate is measured against: the 3
#    SHIV-positive cells in H, its 29 UMI, and the metrics_summary.csv that
#    documents the missing chemistry row. Writing to pipeline_results/ instead
#    means the two runs coexist and can be compared directly, and this script
#    refuses to run if OUTDIR resolves anywhere inside Analysis/.
#
# 4. THE CHEMISTRY NAME IS CHECKED BEFORE A FIVE-DAY JOB STARTS, against the
#    installed Cell Ranger rather than against the documentation.
#
# 5. A POST-RUN REPORT pulls the chemistry row, the cell and UMI metrics, and
#    the strand split of viral reads. That last one matters: the collapsed
#    feature is plus-strand only and the STARsolo recovery needed
#    `soloStrand REVERSE`, so if Cell Ranger's 5' handling lands the other way
#    the viral count is zero for a reason that has nothing to do with the GTF.
#
# UNCHANGED, DELIBERATELY
# -----------------------
#   the corrected library-to-condition mapping, VDJ_RH (not VDJ), the IMGT
#   rhesus reference and inner enrichment primers, feature_reference.csv,
#   create-bam true (the barcode arm reads these BAMs), expect-cells 10000
#
# LIBRARY MAP (libraries are SEQUENCING POOLS, not conditions)
#   A = Pre PBMC                  Chinese: 39272, 40702, 40707, 41861
#   B = Pre LN                    Chinese: 39272, 40702, 40707, 41861
#   C = Pre PBMC         [MIXED]  Indian: 34315, 41903 + Chinese 41862
#   D = 21DPI + Pre LN   [MIXED]  21DPI: 39272, 40702, 40707 + Pre 41862
#   E = 21DPI PBMC                Chinese: 39272, 40702, 40707, 41861
#   F = 21DPI + NonInf LN [MIXED] 21DPI: 41861, 41862 + non-inf 34315, 41903
#   G = Necropsy PBMC             Chinese: 40702, 40707, 41861, 41862
#   H = Necropsy LN               Chinese: 40702, 40707, 41861, 41862
# Per-animal condition is assigned downstream from the provider's genetic
# demultiplexing, never from the library letter.
#
# USAGE
#   # THE GATE. H is index 7 and holds 422 of the 711 viral cell barcodes.
#   # G is index 6, Necropsy PBMC, documented zero viral signal across 17,019
#   # cells. Run them TOGETHER so the control moves with the signal: the GTF
#   # collapse raises counted reads in both, and only the ratio is informative.
#   sbatch --array=6,7 03_cellranger_multi.sh
#
#   # then, if the gate passes
#   sbatch --array=0-5 03_cellranger_multi.sh
#
# GATE CRITERIA, in the order each one tells you something
#   1. H SHIV-positive cells rise from 3. Below about 20 the GTF was not the
#      bottleneck and the full array is not justified.
#   2. H total SHIV UMI rises from 29. Separates more cells from more
#      molecules per cell.
#   3. G stays at or near zero. If G rises with H, the collapse is recovering
#      background as fast as signal and the clean-negatives claim needs
#      restating before any figure is drawn.
#   4. Reads spanning the cassette rise from 7. The only criterion that bears
#      on the barcode arm.
#   5. Host cell count and median UMI per cell stay sane. A large shift means
#      paired-end moved the QC thresholds and annotation needs rework.
#
# Author: Jake Lehle
# Date:   October 2026
###############################################################################

set -euo pipefail

#--- Configuration -----------------------------------------------------------
WORKDIR="/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV"
REFDIR="/master/jlehle/WORKING/SC/REF"
CRDIR="/master/jlehle/cellranger-10.0.0"
FASTQ_BASE="/master/zwallis/WORKING/SC/10X_RAW"

# New output tree. The v10 baseline lives in ${WORKDIR}/Analysis and is never
# written to by this script.
OUTDIR="${WORKDIR}/pipeline_results/03_cellranger_multi"
CONFIG_DIR="${OUTDIR}/multi_configs"
REPORT_DIR="${OUTDIR}/reports"
LEGACY_ANALYSIS="${WORKDIR}/Analysis"

# References
GEX_REF="${REFDIR}/Mmul10_SHIVAD8EO_v11"
VDJ_REF="${WORKDIR}/vdj_reference/vdj_rhesus_macaque_IMGT"
INNER_PRIMERS="${WORKDIR}/vdj_reference/inner_enrichment_primers.txt"
FEATURE_REF="${WORKDIR}/feature_reference.csv"

VIRAL_CONTIG="SHIVAD8EO_V4fullN"
VIRAL_GENE="SHIVAD8EO"

# THE change. Set explicitly so it can never silently auto-detect again.
CHEMISTRY="SC5P-PE-v3"
EXPECT_CELLS=10000
CREATE_BAM="true"          # the barcode arm reads these BAMs

#--- Library definitions (index 0-7) -----------------------------------------
LETTERS=(A B C D E F G H)
CONDITIONS=(
    "Pre_PBMC" "Pre_LN" "Pre_PBMC_IndianMixed" "21DPI_Pre_LN_Mixed"
    "21DPI_PBMC" "21DPI_NonInf_LN_Mixed" "Necropsy_PBMC" "Necropsy_LN"
)
GEX_PREFIXES=(
    "GEX_LIBRARY_A-SI_TT_A1_23552KLT4" "GEX_LIBRARY_B-SI_TT_B1_23552KLT4"
    "GEX_LIBRARY_C-SI_TT_C1_23552KLT4" "GEX_LIBRARY_D-SI_TT_D1_23552KLT4"
    "GEX_LIBRARY_E-SI_TT_E1_23552KLT4" "GEX_LIBRARY_F-SI_TT_F1_23552KLT4"
    "GEX_LIBRARY_G-SI_TT_G1_23552KLT4" "GEX_LIBRARY_H-SI_TT_H1_23552KLT4"
)
PROTEIN_PREFIXES=(
    "PROTEIN_LIBRARY_A-SI_TN_A1_23552KLT4" "PROTEIN_LIBRARY_B-SI_TN_B1_23552KLT4"
    "PROTEIN_LIBRARY_C-SI_TN_C1_23552KLT4" "PROTEIN_LIBRARY_D-SI_TN_D1_23552KLT4"
    "PROTEIN_LIBRARY_E-SI_TN_E1_23552KLT4" "PROTEIN_LIBRARY_F-SI_TN_F1_23552KLT4"
    "PROTEIN_LIBRARY_G-SI_TN_G1_23552KLT4" "PROTEIN_LIBRARY_H-SI_TN_H1_23552KLT4"
)
VDJ_RH_PREFIXES=(
    "VDJ_RH_LIBRARY_A-SI_TT_A3_23552KLT4" "VDJ_RH_LIBRARY_B-SI_TT_B3_23552KLT4"
    "VDJ_RH_LIBRARY_C-SI_TT_C3_23552KLT4" "VDJ_RH_LIBRARY_D-SI_TT_D3_23552KLT4"
    "VDJ_RH_LIBRARY_E-SI_TT_E3_23552KLT4" "VDJ_RH_LIBRARY_F-SI_TT_F3_23552KLT4"
    "VDJ_RH_LIBRARY_G-SI_TT_G3_23552KLT4" "VDJ_RH_LIBRARY_H-SI_TT_H3_23552KLT4"
)

hr () { echo; echo "=============================================================="; echo "$1"; echo "=============================================================="; }
die () { echo; echo "  STOP: $1"; exit 1; }

mkdir -p "${OUTDIR}" "${CONFIG_DIR}" "${REPORT_DIR}"

IDX=${SLURM_ARRAY_TASK_ID:-7}        # default to H, the gate library
LTR="${LETTERS[$IDX]}"
COND="${CONDITIONS[$IDX]}"
SAMPLE_ID="SHIV_${COND}"

hr "03 — Cell Ranger multi | library ${LTR} | ${COND}"
echo "  reference : ${GEX_REF}"
echo "  chemistry : ${CHEMISTRY}"
echo "  output    : ${OUTDIR}/${SAMPLE_ID}"


###############################################################################
hr "PREFLIGHT"
###############################################################################

# --- the v10 baseline must be unreachable from here -------------------------
case "$(readlink -m "${OUTDIR}")" in
    "$(readlink -m "${LEGACY_ANALYSIS}")"|"$(readlink -m "${LEGACY_ANALYSIS}")"/*)
        die "OUTDIR resolves inside ${LEGACY_ANALYSIS}. That tree is the v10
        baseline the gate is measured against and must not be overwritten." ;;
esac
echo "  [OK  ] OUTDIR is outside Analysis/; the v10 baseline is safe"
if [ -d "${LEGACY_ANALYSIS}/${SAMPLE_ID}" ]; then
    echo "         v10 baseline for ${LTR}: ${LEGACY_ANALYSIS}/${SAMPLE_ID}"
fi

# --- references -------------------------------------------------------------
[ -d "${GEX_REF}" ] || die "GEX reference not found: ${GEX_REF}
        Run 02_mkref_Mmul10_SHIVAD8EO.sh first."
[ -f "${GEX_REF}/reference.json" ] || die "GEX reference looks incomplete (no reference.json)"
if grep -qx "${VIRAL_CONTIG}" "${GEX_REF}/star/chrName.txt" 2>/dev/null; then
    echo "  [OK  ] viral contig ${VIRAL_CONTIG} present in the reference"
else
    die "viral contig ${VIRAL_CONTIG} is NOT in ${GEX_REF}/star/chrName.txt.
        This is the wrong reference, or 02 did not finish."
fi

[ -d "${VDJ_REF}" ]       || die "VDJ reference not found: ${VDJ_REF}"
[ -f "${INNER_PRIMERS}" ] || die "inner enrichment primers not found: ${INNER_PRIMERS}"
echo "  [OK  ] VDJ reference and inner primers"

if [ -f "${FEATURE_REF}" ]; then
    INCLUDE_PROTEIN=true
    echo "  [OK  ] feature reference: ${FEATURE_REF}"
else
    INCLUDE_PROTEIN=false
    echo "  WARN  feature reference missing; CITE-seq will be skipped"
fi

# --- chemistry name, checked against the install, not the docs --------------
# A wrong name fails Cell Ranger's own preflight in seconds rather than days,
# so this is a convenience check and deliberately warn-and-continue: the
# chemistry table may be compiled into the binary rather than sitting in a
# readable file, and absence of the string does not prove absence of support.
if grep -rqs --binary-files=text "${CHEMISTRY}" "${CRDIR}/lib" "${CRDIR}/bin" 2>/dev/null; then
    echo "  [OK  ] chemistry ${CHEMISTRY} found in the Cell Ranger install"
else
    echo "  WARN  could not find the string ${CHEMISTRY} under ${CRDIR}."
    echo "        It may be compiled in. Cell Ranger validates the name during"
    echo "        its own preflight, so a bad name costs seconds, not the job."
fi

# --- FASTQs -----------------------------------------------------------------
GEX_FASTQ_DIR="${FASTQ_BASE}/GEX_LIBRARY_${LTR}"
VDJ_FASTQ_DIR="${FASTQ_BASE}/VDJ_RH_LIBRARY_${LTR}"
PROTEIN_FASTQ_DIR="${FASTQ_BASE}/PROTEIN_LIBRARY_${LTR}"

for FDIR in "${GEX_FASTQ_DIR}" "${VDJ_FASTQ_DIR}"; do
    [ -d "${FDIR}" ] || die "FASTQ directory not found: ${FDIR}"
done
if ${INCLUDE_PROTEIN} && [ ! -d "${PROTEIN_FASTQ_DIR}" ]; then
    die "PROTEIN FASTQ directory not found: ${PROTEIN_FASTQ_DIR}"
fi
echo "  [OK  ] FASTQ directories present"

# Paired-end needs a real R1. Libraries B and H were repaired after transfer
# corruption (B to 82.3%, H to 18.1% for VDJ), so confirm R1 is still there.
N_R1=$(find "${GEX_FASTQ_DIR}" -name '*_R1_*.fastq.gz' | wc -l)
N_R2=$(find "${GEX_FASTQ_DIR}" -name '*_R2_*.fastq.gz' | wc -l)
echo "  GEX FASTQs: ${N_R1} R1 file(s), ${N_R2} R2 file(s)"
[ "${N_R1}" -gt 0 ] || die "no R1 FASTQs found. ${CHEMISTRY} aligns BOTH mates
        and cannot run without R1."
if [ "${N_R1}" -ne "${N_R2}" ]; then
    echo "  WARN  R1 and R2 file counts differ (${N_R1} vs ${N_R2})."
    echo "        10x FASTQs are generated in lockstep; check"
    echo "        diagnose_truncated_fastqs.sh output for library ${LTR}."
fi


###############################################################################
hr "CONFIG"
###############################################################################

CONFIG_CSV="${CONFIG_DIR}/multi_config_${LTR}_${COND}.csv"

cat > "${CONFIG_CSV}" <<CONFIGEOF
[gene-expression]
reference,${GEX_REF}
chemistry,${CHEMISTRY}
create-bam,${CREATE_BAM}
expect-cells,${EXPECT_CELLS}

[vdj]
reference,${VDJ_REF}
inner-enrichment-primers,${INNER_PRIMERS}

CONFIGEOF

if ${INCLUDE_PROTEIN}; then
    cat >> "${CONFIG_CSV}" <<CONFIGEOF2
[feature]
reference,${FEATURE_REF}

CONFIGEOF2
fi

{
    echo "[libraries]"
    echo "fastq_id,fastqs,feature_types"
    echo "${GEX_PREFIXES[$IDX]},${GEX_FASTQ_DIR},Gene Expression"
    echo "${VDJ_RH_PREFIXES[$IDX]},${VDJ_FASTQ_DIR},VDJ-T"
    ${INCLUDE_PROTEIN} && \
        echo "${PROTEIN_PREFIXES[$IDX]},${PROTEIN_FASTQ_DIR},Antibody Capture"
} >> "${CONFIG_CSV}"

echo "  wrote ${CONFIG_CSV}"
echo "  --- config ---"
sed 's/^/    /' "${CONFIG_CSV}"
echo "  --- end ---"

grep -q "^chemistry,${CHEMISTRY}$" "${CONFIG_CSV}" \
    || die "the chemistry line is missing from the generated config."
echo "  [OK  ] chemistry line present in the config"


###############################################################################
hr "RUN"
###############################################################################

cd "${OUTDIR}"

# Only ever removes inside OUTDIR, which the preflight proved is not Analysis/.
if [ -d "${SAMPLE_ID}" ]; then
    echo "  removing a previous v11 run of ${SAMPLE_ID} (inside OUTDIR only)..."
    rm -rf "${OUTDIR:?}/${SAMPLE_ID}"
fi

echo "  starting $(date -u +%Y-%m-%dT%H:%M:%SZ)"
"${CRDIR}/cellranger" multi \
    --id="${SAMPLE_ID}" \
    --csv="${CONFIG_CSV}" \
    --localcores="${SLURM_CPUS_PER_TASK:-24}" \
    --localmem=$(( ${SLURM_MEM_PER_NODE:-819200} / 1024 ))
echo "  finished $(date -u +%Y-%m-%dT%H:%M:%SZ)"


###############################################################################
hr "POST-RUN REPORT — library ${LTR}"
###############################################################################

OUTS="${OUTDIR}/${SAMPLE_ID}/outs"
REPORT="${REPORT_DIR}/${LTR}_${COND}_report.txt"
: > "${REPORT}"
rep () { echo "$1" | tee -a "${REPORT}"; }

rep "library            ${LTR}  (${COND})"
rep "reference          ${GEX_REF}"
rep "chemistry (config) ${CHEMISTRY}"
rep "finished           $(date -u +%Y-%m-%dT%H:%M:%SZ)"
rep ""

# --- 1. did the chemistry actually take -------------------------------------
MET=$(find "${OUTS}" -name 'metrics_summary.csv' 2>/dev/null | head -1 || true)
if [ -n "${MET}" ]; then
    rep "metrics_summary    ${MET}"
    # `cellranger multi` writes six columns:
    #   Category,Library Type,Grouped By,Group Name,Metric Name,Metric Value
    # and emits NO chemistry row at all, so its absence says nothing. The
    # decisive chemistry test is whether the BAM carries paired alignments,
    # done below. Filter on Library Type or the Antibody Capture rows get
    # reported as if they were gene expression, which the 2026-10-30 gate run
    # did.
    rep "chemistry row      not emitted by \`multi\`; see the paired-read check below"
    rep "gene expression metrics:"
    awk -F',' 'NR>1 && $2 ~ /Gene Expression/ {
                   printf "  %-44s %s\n", $5, $6 }' "${MET}" \
        | sed 's/^/  /' | tee -a "${REPORT}" || true
else
    rep "metrics_summary    NOT FOUND under ${OUTS}"
fi
rep ""

# --- 2. viral reads, and which strand they are on ---------------------------
# The collapsed feature is plus-strand only. If viral reads are plentiful but
# overwhelmingly antisense, the viral count will be zero for a strand reason
# and NOT because the GTF collapse failed. Read this before the gate criteria.
# The system samtools at /usr/local/bin is an ancient broken build. Prefer the
# conda one explicitly rather than warning after the fact, which is what the
# 2026-10-30 gate run did, returning NA for both strand counts.
SAMTOOLS=""
for CAND in "${CONDA_PREFIX:-}/bin/samtools" "$(command -v samtools || true)"; do
    case "${CAND}" in ""|/usr/local/bin/*) continue ;; esac
    [ -x "${CAND}" ] && { SAMTOOLS="${CAND}"; break; }
done
BAM=$(find "${OUTS}" -name 'sample_alignments.bam' 2>/dev/null | head -1 || true)
if [ -z "${SAMTOOLS}" ]; then
    rep "WARN no usable samtools. The system build is broken and is ignored."
    rep "     conda activate sc_pre && export PATH=\"\${CONDA_PREFIX}/bin:\${PATH}\""
fi
if [ -n "${SAMTOOLS}" ] && [ -n "${BAM}" ]; then
    rep "samtools           ${SAMTOOLS}"
    rep "bam                ${BAM}"
    FWD=$("${SAMTOOLS}" view -c -F 0x914 "${BAM}" "${VIRAL_CONTIG}" 2>/dev/null || echo NA)
    REV=$("${SAMTOOLS}" view -c -F 0x904 -f 0x10 "${BAM}" "${VIRAL_CONTIG}" 2>/dev/null || echo NA)
    # SC5P-PE-v3 aligns both mates; SC5P-R2-v3 aligns R2 only. Paired
    # alignments in the BAM are therefore direct evidence the config took.
    PTOT=$("${SAMTOOLS}" view -c -F 0x900 "${BAM}" 2>/dev/null || echo 0)
    PPAIR=$("${SAMTOOLS}" view -c -f 0x1 -F 0x900 "${BAM}" 2>/dev/null || echo 0)
    rep "primary reads      ${PTOT}"
    rep "paired alignments  ${PPAIR}"
    if [ "${PTOT}" -gt 0 ] && [ "${PPAIR}" -gt $(( PTOT / 2 )) ]; then
        rep "chemistry (actual) PAIRED-END, consistent with ${CHEMISTRY}"
    else
        rep "chemistry (actual) SINGLE-END. ${CHEMISTRY} did NOT take effect."
        rep "                   Do not compare viral counts against the"
        rep "                   baseline until this is resolved."
    fi
    rep ""
    rep "viral reads  fwd   ${FWD}"
    rep "viral reads  rev   ${REV}"
    if [ "${FWD}" != "NA" ] && [ "${REV}" != "NA" ]; then
        TOT=$(( FWD + REV ))
        rep "viral reads  total ${TOT}"
        if [ "${TOT}" -gt 0 ] && [ "${REV}" -gt $(( TOT * 4 / 5 )) ]; then
            rep ""
            rep "*** Over 80% of viral reads are on the reverse strand while the"
            rep "    GTF feature is plus-strand only. A near-zero viral count"
            rep "    here is a STRAND problem, not a GTF problem. Resolve this"
            rep "    before judging the gate."
        fi
    fi
else
    rep "viral strand check SKIPPED (samtools or BAM not found)"
fi
rep ""

# --- 3. the viral feature in the matrix -------------------------------------
H5=$(find "${OUTS}" -name 'sample_filtered_feature_bc_matrix.h5' 2>/dev/null | head -1 || true)
if [ -n "${H5}" ]; then
    rep "filtered matrix    ${H5}"
    rep "                   count the viral feature with:"
    rep "                     python -c \"import scanpy as sc; a=sc.read_10x_h5('${H5}');"
    rep "                     import numpy as np; m=a.var_names=='${VIRAL_GENE}';"
    rep "                     x=a[:,m].X; print('SHIV UMI', x.sum(), 'cells', (x>0).sum())\""
else
    rep "filtered matrix    NOT FOUND"
fi

rep ""
rep "v10 baseline for comparison: ${LEGACY_ANALYSIS}/${SAMPLE_ID}/outs"
rep "  H baseline was 3 SHIV-positive cells, 29 SHIV UMI, 7 cassette reads."
rep "  G baseline was zero viral signal across 17,019 cells."

echo
echo "  report written: ${REPORT}"

hr "DONE — ${SAMPLE_ID}"
cat <<EOF

  outs      : ${OUTS}
  report    : ${REPORT}
  baseline  : ${LEGACY_ANALYSIS}/${SAMPLE_ID}/outs   (untouched)

  If this was the gate (--array=6,7), read H and G together:
    H SHIV+ cells 3 -> ?      criterion 1, needs to clear ~20
    H SHIV UMI   29 -> ?      criterion 2
    G viral signal 0 -> ?     criterion 3, must stay near zero
  A rise in H that is matched proportionally in G is not a result.

  Then: sbatch --array=0-5 03_cellranger_multi.sh

  BEFORE stage 04, scanpy_analysis_script.py section 3 must become
      adata.var['shiv'] = adata.var_names == '${VIRAL_GENE}'
  with the case-insensitive fallback deleted, or host TAT is counted as virus
  with no error raised.

EOF
exit 0
