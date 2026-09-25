#!/bin/bash
#SBATCH --job-name=shiv_starsolo
#SBATCH --output=/master/jlehle/WORKING/LOGS/shiv_starsolo_%j.out
#SBATCH --error=/master/jlehle/WORKING/LOGS/shiv_starsolo_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --partition=normal

###############################################################################
# shiv_starsolo_recovery.sh
#
# Barcode-aware recovery of SHIV reads that CellRanger's STAR discarded.
#
# WHY: CellRanger drops genuine SHIV reads at the alignment LENGTH filter
#      (outFilter*OverLread = 0.66). minimap2 proved the reads are real
#      (Library H: 553 recovered, 91.4% genome breadth) but the FASTQ
#      extraction threw away CB/UMI, so those reads can't become cells.
#      STARsolo runs the ORIGINAL R1+R2, keeps barcodes, relaxes the length
#      filter, and emits a per-CELL SHIV matrix we merge into the strict host
#      counts in Scanpy.
#
# CONFIRMED barcode config (from check_sample_whitelist.py, empirical):
#   chemistry  : GEM-X 5' v3   (UMI length 12)
#   CB len 16, UMI len 12
#   whitelist  : 3M-5pgex-jan-2023.txt.gz   (100% match to observed barcodes)
#
# THREE DESIGN CHOICES (flagged for Jake — adjust at top if you disagree):
#   1. Collapsed GTF: ONE feature spanning the whole 10,262bp contig, so every
#      sense viral read counts regardless of which (overlapping) ORF it hits.
#   2. Relaxed length filters set to 0, with an absolute 30bp floor + 10% mismatch
#      cap. SHIV-only index => no host to cross-map to (the asymmetry).
#   3. soloStrand REVERSE. 10x 5' cDNA (R2) maps antisense to the +-strand gene.
#      Forward miscounted (458 aligned -> 7 gene UMIs); Reverse fixed it (-> ~337
#      gene reads, 31 SHIV+ cells in H, on top of CellRanger's original 27).
#      NOTE: aligned-read count / breadth are strand-INDEPENDENT, so they cannot
#      validate strand — only the Gene-counting rate (Solo Summary) can.
#
# RUN ORDER: validated on H, now running all 8. Negatives (A/B/C/G) are the
# false-positive control and should stay ~0 at the 2-UMI cell threshold.
#
# Author: Jake Lehle
# Date: June 2026
###############################################################################

set -euo pipefail
source ~/anaconda3/bin/activate
eval "$(conda shell.bash hook)"
conda activate sc_pre
export PATH="${CONDA_PREFIX}/bin:${PATH}"   # conda samtools (system one is broken)

# ─── Configuration (all adjustable knobs here) ───────────────────────────────
THREADS=${SLURM_CPUS_PER_TASK:-16}
WORKDIR="/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV"
ANALYSIS_DIR="${WORKDIR}/Analysis"
FASTQ_BASE="/master/zwallis/WORKING/SC/10X_RAW"
OUTDIR="${WORKDIR}/shiv_starsolo"
CR_COMBINED_FA="/master/jlehle/cellranger-10.0.0/Mmul10_SHIVAD8EO_v10/fasta/genome.fa"
SHIV_CONTIG="SHIVAD8EO"

# STAR binary: need a MODERN STAR (>=2.7.4) for the CellRanger-replication solo
# flags. CellRanger's bundled STAR is 2.7.2a (2019) and predates them — do NOT
# use it. Install once:  conda create -y -n star -c bioconda -c conda-forge star=2.7.11b
STAR_ENV_BIN="${HOME}/anaconda3/envs/star/bin/STAR"
STAR_BIN="$(command -v STAR || true)"
[[ -z "${STAR_BIN}" && -x "${STAR_ENV_BIN}" ]] && STAR_BIN="${STAR_ENV_BIN}"
[[ -z "${STAR_BIN}" ]] && STAR_BIN="${STAR_ENV_BIN}"   # default to the star env
MIN_STAR="2.7.4a"

# --- Confirmed barcode config (do not change without re-running the check) ---
CB_LEN=16
UMI_LEN=12
WHITELIST_GZ="/master/jlehle/cellranger-10.0.0/lib/python/cellranger/barcodes/3M-5pgex-jan-2023.txt.gz"

# --- Strand: REVERSE for 10x 5'. The cDNA read (R2) maps antisense to the
#     plus-strand gene, so Forward counted only the ~11% wrong-strand minority
#     (458 aligned -> 7 counted in H). Confirmed: R2 is ~89% minus-strand while
#     the collapsed SHIV gene is on +. Reverse counts the real reads. ---
STRAND="Reverse"

# --- Relaxed alignment filters (the whole point) ---
SCORE_MIN_OVER_L=0          # was 0.66 in CellRanger
MATCH_MIN_OVER_L=0          # was 0.66 in CellRanger
MATCH_MIN_ABS=30            # absolute floor: >=30 matched bases (anti-garbage)
MISMATCH_OVER_L=0.1         # allow up to 10% mismatch (covers the nM 9-10 tail)

# --- Libraries to process: letter -> CellRanger sample name ---
#     Strand fix validated on H (7 -> ~337 gene reads, 31 SHIV+ cells in merge).
#     Full run, all 8: negatives A B C G expected ~0 at the 2-UMI threshold,
#     infected D E F H carry reservoir-appropriate signal.
LIBS="A B C D E F G H"
declare -A LIB_NAME=(
    [A]="SHIV_Pre_PBMC"            [B]="SHIV_Pre_LN"
    [C]="SHIV_Pre_PBMC_IndianMixed" [D]="SHIV_21DPI_Pre_LN_Mixed"
    [E]="SHIV_21DPI_PBMC"          [F]="SHIV_21DPI_NonInf_LN_Mixed"
    [G]="SHIV_Necropsy_PBMC"       [H]="SHIV_Necropsy_LN"
)
# minimap2 recovery baselines (for the strand/recovery sanity check)
declare -A MM2_BASELINE=( [H]=553 [D]=35 [E]=48 )

mkdir -p "${OUTDIR}"/{ref,index,solo,reports}
cd "${OUTDIR}"
LOGFILE="${OUTDIR}/reports/starsolo_recovery.log"
exec > >(tee -a "${LOGFILE}") 2>&1

echo "============================================================"
echo "SHIV STARsolo Recovery — $(date)"
echo "STAR: ${STAR_BIN}"
if [[ ! -x "${STAR_BIN}" ]]; then
    echo "ERROR: no STAR at ${STAR_BIN}"
    echo "  install: conda create -y -n star -c bioconda -c conda-forge star=2.7.11b"
    exit 1
fi
STAR_VER="$(${STAR_BIN} --version)"
echo "STAR version: ${STAR_VER}"
# require >= 2.7.4a (CellRanger-replication solo flags); sort -V comparison
if [[ "$(printf '%s\n%s\n' "${MIN_STAR}" "${STAR_VER}" | sort -V | head -1)" != "${MIN_STAR}" ]]; then
    echo "ERROR: STAR ${STAR_VER} is too old; need >= ${MIN_STAR}."
    echo "  The bundled CellRanger STAR (2.7.2a) lacks soloCellFilter, 1MM_CR, etc."
    echo "  install: conda create -y -n star -c bioconda -c conda-forge star=2.7.11b"
    exit 1
fi
echo "============================================================"

# ─── Prep: SHIV-only fasta + collapsed single-feature GTF ────────────────────
SHIV_FA="${OUTDIR}/ref/SHIVAD8EO.fa"
SHIV_GTF="${OUTDIR}/ref/SHIVAD8EO_collapsed.gtf"

if [[ ! -f "${SHIV_FA}" ]]; then
    echo "[ref] extracting ${SHIV_CONTIG} from combined reference"
    samtools faidx "${CR_COMBINED_FA}" "${SHIV_CONTIG}" > "${SHIV_FA}"
    samtools faidx "${SHIV_FA}"
fi
CONTIG_LEN=$(cut -f2 "${SHIV_FA}.fai")
echo "[ref] ${SHIV_CONTIG} length = ${CONTIG_LEN} bp"

if [[ ! -f "${SHIV_GTF}" ]]; then
    echo "[ref] writing collapsed GTF (one feature over the whole contig)"
    # single gene/transcript/exon spanning 1..CONTIG_LEN on + strand
    printf '%s\tcollapsed\tgene\t1\t%s\t.\t+\t.\tgene_id "SHIVAD8EO"; gene_name "SHIV";\n' \
        "${SHIV_CONTIG}" "${CONTIG_LEN}" >  "${SHIV_GTF}"
    printf '%s\tcollapsed\ttranscript\t1\t%s\t.\t+\t.\tgene_id "SHIVAD8EO"; transcript_id "SHIVAD8EO_t"; gene_name "SHIV";\n' \
        "${SHIV_CONTIG}" "${CONTIG_LEN}" >> "${SHIV_GTF}"
    printf '%s\tcollapsed\texon\t1\t%s\t.\t+\t.\tgene_id "SHIVAD8EO"; transcript_id "SHIVAD8EO_t"; gene_name "SHIV"; exon_number "1";\n' \
        "${SHIV_CONTIG}" "${CONTIG_LEN}" >> "${SHIV_GTF}"
fi
cat "${SHIV_GTF}"

# ─── Whitelist: STARsolo needs a PLAIN-TEXT whitelist (gunzip if needed) ──────
WHITELIST="${OUTDIR}/ref/$(basename "${WHITELIST_GZ}" .gz)"
if [[ ! -f "${WHITELIST}" ]]; then
    echo "[ref] decompressing whitelist for STARsolo"
    zcat "${WHITELIST_GZ}" | cut -f1 > "${WHITELIST}"
fi
echo "[ref] whitelist: ${WHITELIST} ($(wc -l < "${WHITELIST}") barcodes)"

# ─── Build SHIV-only STAR index (small genome => reduced SAindexNbases) ───────
INDEX_DIR="${OUTDIR}/index"
VER_MARKER="${INDEX_DIR}/.star_version"
# STAR's recommended value for a ~10kb genome is 5 (it seg-faults / warns at 6).
SA_NBASES=5
BUILD_TAG="${STAR_VER}_SA${SA_NBASES}"
# Rebuild if missing OR built by a different STAR version / SA setting. The old
# 2.7.2a index cannot be loaded by a modern STAR. Build is sub-second here.
NEED_BUILD=1
if [[ -f "${INDEX_DIR}/SAindex" && -f "${VER_MARKER}" && "$(cat "${VER_MARKER}")" == "${BUILD_TAG}" ]]; then
    NEED_BUILD=0
fi
if [[ "${NEED_BUILD}" -eq 1 ]]; then
    rm -rf "${INDEX_DIR}"; mkdir -p "${INDEX_DIR}"
    echo "[index] building STAR index with ${STAR_VER} (genomeSAindexNbases=${SA_NBASES})"
    "${STAR_BIN}" --runMode genomeGenerate \
        --genomeDir "${INDEX_DIR}" \
        --genomeFastaFiles "${SHIV_FA}" \
        --sjdbGTFfile "${SHIV_GTF}" \
        --sjdbOverhang 149 \
        --genomeSAindexNbases ${SA_NBASES} \
        --runThreadN "${THREADS}"
    echo "${BUILD_TAG}" > "${VER_MARKER}"
else
    echo "[index] reusing index built by ${BUILD_TAG} at ${INDEX_DIR}"
fi
echo ""

# ─── Process each library ────────────────────────────────────────────────────
for L in ${LIBS}; do
    NAME="${LIB_NAME[$L]}"
    FQDIR="${FASTQ_BASE}/GEX_LIBRARY_${L}"
    PREFIX="${OUTDIR}/solo/${L}_${NAME}_"

    echo "============================================================"
    echo "Library ${L}: ${NAME}"
    echo "============================================================"

    if [[ ! -d "${FQDIR}" ]]; then
        echo "  [WARN] FASTQ dir not found, skipping: ${FQDIR}"
        continue
    fi

    # cDNA = R2, barcode = R1 ; matched lane order, comma-separated
    R2S=$(ls "${FQDIR}"/*_R2_*.fastq.gz 2>/dev/null | sort | paste -sd, -)
    R1S=$(ls "${FQDIR}"/*_R1_*.fastq.gz 2>/dev/null | sort | paste -sd, -)
    if [[ -z "${R2S}" || -z "${R1S}" ]]; then
        echo "  [WARN] missing R1/R2 in ${FQDIR}, skipping"
        continue
    fi
    echo "  R2 (cDNA):    ${R2S}"
    echo "  R1 (barcode): ${R1S}"

    # ── STARsolo: relaxed length filters, CellRanger-style barcode/UMI logic ──
    "${STAR_BIN}" \
        --runMode alignReads \
        --genomeDir "${INDEX_DIR}" \
        --readFilesIn "${R2S}" "${R1S}" \
        --readFilesCommand zcat \
        --runThreadN "${THREADS}" \
        --outFileNamePrefix "${PREFIX}" \
        --soloType CB_UMI_Simple \
        --soloCBwhitelist "${WHITELIST}" \
        --soloCBstart 1 --soloCBlen ${CB_LEN} \
        --soloUMIstart $((CB_LEN + 1)) --soloUMIlen ${UMI_LEN} \
        --soloBarcodeReadLength 0 \
        --soloStrand ${STRAND} \
        --soloFeatures Gene \
        --soloCellFilter None \
        --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
        --soloUMIdedup 1MM_CR \
        --soloUMIfiltering MultiGeneUMI_CR \
        --clipAdapterType CellRanger4 \
        --outFilterScoreMinOverLread ${SCORE_MIN_OVER_L} \
        --outFilterMatchNminOverLread ${MATCH_MIN_OVER_L} \
        --outFilterMatchNmin ${MATCH_MIN_ABS} \
        --outFilterMismatchNoverLmax ${MISMATCH_OVER_L} \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes NH HI AS nM CB UB \
        --limitBAMsortRAM 60000000000

    BAM="${PREFIX}Aligned.sortedByCoord.out.bam"
    samtools index "${BAM}"

    # ── Validation: aligned-read count + genome breadth vs minimap2 baseline ──
    ALN=$(samtools view -c -F 0x4 "${BAM}")
    BREADTH=$(samtools depth -a "${BAM}" | awk '{n++; if($3>0)c++} END{printf "%.1f", (n>0?100*c/n:0)}')
    HQ=$(samtools view -c -q 30 "${BAM}")
    echo ""
    echo "  ── RECOVERY (Library ${L}) ──"
    echo "    aligned reads (STARsolo): ${ALN}"
    echo "    MAPQ>=30:                 ${HQ}"
    echo "    genome breadth:           ${BREADTH}%"
    if [[ -n "${MM2_BASELINE[$L]:-}" ]]; then
        echo "    minimap2 baseline:        ${MM2_BASELINE[$L]} (strand sanity target)"
        echo "    >> if STARsolo aligned ~0 but minimap2 found ${MM2_BASELINE[$L]},"
        echo "       set STRAND=Reverse at top and rerun."
    fi

    # ── Solo matrix summary (the actual deliverable: per-cell SHIV counts) ──
    SUMMARY="${PREFIX}Solo.out/Gene/Summary.csv"
    RAWDIR="${PREFIX}Solo.out/Gene/raw"
    echo ""
    echo "  ── Solo summary ──"
    [[ -f "${SUMMARY}" ]] && sed 's/^/    /' "${SUMMARY}" || echo "    [no Summary.csv]"
    if [[ -d "${RAWDIR}" ]]; then
        # gzip-agnostic: count barcodes with >0 SHIV UMI from the raw matrix
        MTX="${RAWDIR}/matrix.mtx"; [[ -f "${MTX}.gz" ]] && MTX="${MTX}.gz"
        echo "    raw matrix: ${MTX}"
        echo "    (merge raw matrix barcodes with CellRanger filtered cells in Scanpy)"
    fi
    echo ""
done

echo "============================================================"
echo "DONE — $(date)"
echo ""
echo "NEXT:"
echo "  1. Confirm Library H recovery is in the minimap2 ballpark (strand OK)."
echo "  2. Add D E + Pre negatives to LIBS; rerun (negatives should stay ~0)."
echo "  3. In Scanpy: load each Solo.out/Gene/raw matrix, subset to the cell"
echo "     barcodes CellRanger already called, add SHIV as a count layer/obs."
echo "============================================================"
