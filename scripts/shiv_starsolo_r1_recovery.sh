#!/bin/bash
#SBATCH -J shiv_starsolo_r1
#SBATCH -o /master/jlehle/WORKING/LOGS/shiv_starsolo_r1_recovery.o.log
#SBATCH -e /master/jlehle/WORKING/LOGS/shiv_starsolo_r1_recovery.e.log
#SBATCH --partition=normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=48:00:00
# 48h, not 24. The R1 pass runs ~11,500 reads/min against the R2 pass's
# ~62,000, because R1 is listed TWICE in --readFilesIn and the gzip is
# therefore decompressed twice per library. At ~120M reads that is ~2.9h per
# library, so 8 libraries on one strand lands at ~23h: inside a 24h wall only
# if nothing slows down. Do not lose 23 hours of compute to a wall clock.
###############################################################################
# shiv_starsolo_r1_recovery.sh
#
# Aligns the cDNA in R1, which every run in this project has thrown away.
#
# WHY
# ---
# 10x 5' R1 is not a barcode-only read. It is
#     [CB 16][UMI 12][TSO 13] then cDNA, sense strand, from the molecule's 5' end
# and shiv_r1_structure_check.py measured exactly that on all 8 libraries: TSO
# at offset 28, cDNA from base 41, 109bp of it, poly-G 0.08-2.09%, entropy flat
# at ~1.90 bits across all five 20bp blocks. Real transcript, not dark-cycle
# artifact.
#
# --soloType CB_UMI_Simple reads nothing past base 28, so those 109bp have
# never been through any aligner. That matters because shiv_vpr_flank_scan.py
# found 11 of its 16 cassette hits in R1 and only 5 in R2. The half of the data
# that produced two thirds of the signal has never been aligned.
#
# WHY SINGLE-END R1 AND NOT PAIRED-END
# ------------------------------------
# The obvious approach is --soloBarcodeMate 1 with --readFilesIn R1 R2, which
# aligns both mates as a pair. It was rejected: STAR needs both mates to place,
# so a pair whose R2 fails but whose R1 would have succeeded is DROPPED. That
# can lose reads the existing R2-only run already finds, turning an additive
# experiment into a trade. Run R1 on its own and union on read id instead;
# nothing already recovered can be lost that way.
#
# No pre-trimming of 130GB of FASTQ is needed. R1 is passed TWICE:
#   first file   R1 with --clip5pNbases ${CLIP5}  -> the cDNA, aligned
#   second file  R1 unclipped                     -> the barcode, read at 1..28
#
# NOTE ON --clip5pNbases ARITY (2026-10-28). Under --soloType CB_UMI_Simple the
# barcode file is NOT a mate: --readFilesIn R1 R1 is ONE cDNA mate plus one
# barcode read. So --clip5pNbases takes ONE value. Passing two draws
#   "--clip5pNbases has to contain 1 values to match the number of mates"
# which is also confirmation that STAR accepts the same file listed twice.
#
# THAT RAISES ONE RISK WORTH CHECKING: if clipping were applied to the barcode
# read as well, CB would be read from base 42 and every barcode would fail the
# whitelist. CELL output below prints Solo.out/Barcodes.stats so that failure
# cannot pass silently; a collapsed match rate means clipping hit the barcode
# read and the FALLBACK at the bottom (pre-trim, no --clip5pNbases) is required.
#
# THE STRAND CALL IS THE RISK
# ---------------------------
# R1's cDNA runs SENSE (it reads out of the molecule's 5' end), while R2 is the
# antisense mate. So soloStrand should flip from Reverse to Forward. That is
# reasoning, not a measurement, and this exact parameter already cost this
# project once: on R2, Forward counted 7 gene reads where Reverse counted ~337.
# So MODE=smoke runs BOTH settings on a subsample and lets the data pick. Do not
# run the full pass until the smoke test has named a winner.
#
#   MODE=smoke sbatch shiv_starsolo_r1_recovery.sh     # ~20 min, library H
#   MODE=full  STRAND=Forward sbatch shiv_starsolo_r1_recovery.sh
#
# Author: Jake Lehle
# Date: October 2026
###############################################################################

set -euo pipefail
source ~/anaconda3/bin/activate
eval "$(conda shell.bash hook)"
conda activate sc_pre
export PATH="${CONDA_PREFIX}/bin:${PATH}"   # conda samtools; the system one is broken

# ─── Configuration ───────────────────────────────────────────────────────────
MODE="${MODE:-smoke}"                 # smoke | full
THREADS=${SLURM_CPUS_PER_TASK:-16}
WORKDIR="/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV"
FASTQ_BASE="/master/zwallis/WORKING/SC/10X_RAW"
REFDIR="${WORKDIR}/shiv_bc_reference"
OUTDIR="${WORKDIR}/shiv_starsolo_r1"

# V2allN: the pool-diff run showed 99-100% read overlap with the old reference,
# so this changes nothing about recovery, but it is the reference that matches
# the virus Brandon actually sent and the coordinates below assume it.
INDEX="${REFDIR}/idx_SHIVAD8EO_V2allN"
CONTIG="SHIVAD8EO_V2allN"
CASSETTE_START=5862
CASSETTE_END=5896

STAR_ENV_BIN="${HOME}/anaconda3/envs/star/bin/STAR"
STAR_BIN="$(command -v STAR || true)"
[[ -z "${STAR_BIN}" || ! -x "${STAR_BIN}" ]] && STAR_BIN="${STAR_ENV_BIN}"

# --- Read geometry, MEASURED by shiv_r1_structure_check.py (do not assume) ---
CB_LEN=16
UMI_LEN=12
TSO_LEN=13
CLIP5=$((CB_LEN + UMI_LEN + TSO_LEN))      # 41: CB + UMI + TSO
WHITELIST="${WORKDIR}/shiv_starsolo/ref/3M-5pgex-jan-2023.txt"

# --- Filters -----------------------------------------------------------------
# MATCH_MIN_ABS 25 -> 70, set by shiv_r1_alignment_audit.py on the 2026-10-28
# smoke run. At 25 the aligned population was 73% junk: median matched length
# 32 of 109 available bases, ~6 mismatches inside those 32, and a 46% reverse
# fraction, which is a coin flip and means no strand preference at all. The
# matched-length histogram was bimodal, a noise mode at 20-39 and a real mode
# at 70-109, so the cut goes between them.
#
# CUT ON MATCHED LENGTH, NOT ON MISMATCH RATE. This is the part that is easy to
# get backwards. A hypermutated read aligns across its WHOLE length carrying
# many mismatches; a spurious read aligns over a SHORT stretch. Tightening
# MISMATCH_OVER_L would discard exactly the APOBEC3-edited reads that alignment
# was adopted to find, while barely touching the noise. So the mismatch budget
# stays at 0.2 and the length floor does the work.
SCORE_MIN_OVER_L=0
MATCH_MIN_OVER_L=0
MATCH_MIN_ABS=70             # of 109bp cDNA; was 25, which admitted the noise mode
MISMATCH_OVER_L=0.2          # reference N is free; keep it for hypermutants

LIB_NAME_A="SHIV_Pre_PBMC";             LIB_NAME_B="SHIV_Pre_LN"
LIB_NAME_C="SHIV_Pre_PBMC_IndianMixed"; LIB_NAME_D="SHIV_21DPI_Pre_LN_Mixed"
LIB_NAME_E="SHIV_21DPI_PBMC";           LIB_NAME_F="SHIV_21DPI_NonInf_LN_Mixed"
LIB_NAME_G="SHIV_Necropsy_PBMC";        LIB_NAME_H="SHIV_Necropsy_LN"
libname() { local v="LIB_NAME_$1"; echo "${!v}"; }
group()   { case "$1" in D|E|F|H) echo INF ;; *) echo neg ;; esac; }

if [[ "${MODE}" == "smoke" ]]; then
    # A NEGATIVE IS NOT OPTIONAL. The first version of this smoke test ran
    # library H alone, which cannot detect the failure mode that matters most:
    # spurious alignment to a SHIV-only index, where no host sequence competes
    # for a read. H returned 3,623 aligned from 20M reads, 35x the per-million
    # rate of the R2 runs, with a 1.04x strand ratio where R2 gave 48x. A
    # strand tie is the signature of alignments with no strand preference,
    # i.e. noise. Only a library with no virus in it can tell the two apart.
    LIBS="${LIBS:-A H}"        # A = SHIV_Pre_PBMC, pre-infection: the floor
    STRANDS="Forward Reverse"
    READ_CAP=20000000          # first 20M reads; enough to separate a 50x effect
else
    LIBS="${LIBS:-A B C D E F G H}"
    # Forward, settled by the 2026-10-28 smoke test: H Solo-counted 119 on
    # Forward against 1 on Reverse. R1 reads sense out of the TSO while R2 is
    # the antisense mate, so this is the opposite of the R2 runs' Reverse.
    STRANDS="${STRAND:-Forward}"
    READ_CAP=-1                # -1 = no cap
    # If the 48h wall is still a worry on a busy partition, split it:
    #   LIBS="A B C D" MODE=full STRAND=Forward sbatch shiv_starsolo_r1_recovery.sh
    #   LIBS="E F G H" MODE=full STRAND=Forward sbatch shiv_starsolo_r1_recovery.sh
    # Outputs land in per-strand directories keyed by library, so two jobs
    # write to the same tree without colliding.
fi

mkdir -p "${OUTDIR}"/{solo,reports}
cd "${OUTDIR}"
LOGFILE="${OUTDIR}/reports/starsolo_r1_${MODE}.log"
exec > >(tee -a "${LOGFILE}") 2>&1

echo "============================================================"
echo "SHIV STARsolo R1-cDNA Recovery — $(date)"
echo "mode            : ${MODE}"
echo "STAR            : ${STAR_BIN}"
if [[ ! -x "${STAR_BIN}" ]]; then
    echo "ERROR: no STAR at ${STAR_BIN}"; exit 1
fi
echo "STAR version    : $(${STAR_BIN} --version)"
echo "index           : ${INDEX}"
echo "clip5pNbases    : ${CLIP5}  (CB ${CB_LEN} + UMI ${UMI_LEN} + TSO ${TSO_LEN})"
echo "strand(s)       : ${STRANDS}"
echo "libraries       : ${LIBS}"
echo "read cap        : ${READ_CAP}"
echo "============================================================"

if [[ ! -d "${INDEX}" ]]; then
    echo "ERROR: index not found: ${INDEX}"
    echo "  run shiv_build_bc_reference.py then build_and_test_indices.sh"
    exit 1
fi

SUMMARY="${OUTDIR}/reports/r1_recovery_${MODE}.tsv"
printf 'strand\tlib\tsample\tgroup\tinput\taligned\tmapq30\tbreadth_pct\tcassette_span\tsolo_counted\n' \
    > "${SUMMARY}"

for STR in ${STRANDS}; do
  echo ""
  echo "############################################################"
  echo "# soloStrand ${STR}"
  echo "############################################################"
  for L in ${LIBS}; do
    NAME="$(libname "${L}")"; GRP="$(group "${L}")"
    FQDIR="${FASTQ_BASE}/GEX_LIBRARY_${L}"
    PREFIX="${OUTDIR}/solo/${STR}/${L}_${NAME}_"
    mkdir -p "$(dirname "${PREFIX}")"

    echo ""
    echo "--- ${STR} / library ${L}: ${NAME} [${GRP}] ---"
    R1S=$(ls "${FQDIR}"/*_R1_*.fastq.gz 2>/dev/null | sort | paste -sd, -)
    if [[ -z "${R1S}" ]]; then
        echo "  [WARN] no R1 in ${FQDIR}, skipping"; continue
    fi

    CAP_ARG=()
    [[ "${READ_CAP}" != "-1" ]] && CAP_ARG=(--readMapNumber "${READ_CAP}")

    # R1 passed TWICE: mate 1 clipped past the TSO is the cDNA, mate 2
    # unclipped supplies CB+UMI. --clip5pNbases takes one value per mate.
    "${STAR_BIN}" \
        --runMode alignReads \
        --genomeDir "${INDEX}" \
        --readFilesIn "${R1S}" "${R1S}" \
        --readFilesCommand zcat \
        --runThreadN "${THREADS}" \
        --outFileNamePrefix "${PREFIX}" \
        "${CAP_ARG[@]}" \
        --clip5pNbases ${CLIP5} \
        --soloType CB_UMI_Simple \
        --soloCBwhitelist "${WHITELIST}" \
        --soloCBstart 1 --soloCBlen ${CB_LEN} \
        --soloUMIstart $((CB_LEN + 1)) --soloUMIlen ${UMI_LEN} \
        --soloBarcodeReadLength 0 \
        --soloStrand ${STR} \
        --soloFeatures Gene \
        --soloCellFilter None \
        --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
        --soloUMIdedup 1MM_CR \
        --soloUMIfiltering MultiGeneUMI_CR \
        --outFilterScoreMinOverLread ${SCORE_MIN_OVER_L} \
        --outFilterMatchNminOverLread ${MATCH_MIN_OVER_L} \
        --outFilterMatchNmin ${MATCH_MIN_ABS} \
        --outFilterMismatchNoverLmax ${MISMATCH_OVER_L} \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes NH HI AS nM CB UB \
        --limitBAMsortRAM 60000000000

    BAM="${PREFIX}Aligned.sortedByCoord.out.bam"
    samtools index "${BAM}"

    INPUT=$(awk -F'|' '/Number of input reads/{gsub(/[ \t]/,"",$2); print $2}' \
            "${PREFIX}Log.final.out" 2>/dev/null || echo NA)
    # -F 0x904 excludes unmapped, secondary AND supplementary. The earlier
    # -F 0x4 counted secondary alignments too and reported 3,623 where only
    # 2,686 were primary, inflating the apparent yield by 26%.
    ALN=$(samtools view -c -F 0x904 "${BAM}")
    ALN_ALL=$(samtools view -c -F 0x4 "${BAM}")
    SECONDARY=$((ALN_ALL - ALN))
    HQ=$(samtools view -c -q 30 -F 0x904 "${BAM}")
    BREADTH=$(samtools depth -a "${BAM}" \
              | awk '{n++; if($3>0)c++} END{printf "%.1f", (n>0?100*c/n:0)}')
    # One aligned block containing the whole cassette window. A CIGAR that
    # skips it with N reaches past the window while carrying no cassette base.
    SPAN=$(samtools view "${BAM}" "${CONTIG}:$((CASSETTE_START+1))-${CASSETTE_END}" 2>/dev/null \
           | awk -v s=${CASSETTE_START} -v e=${CASSETTE_END} '
               { r=$4-1; c=$6; ok=0;
                 while (match(c, /^[0-9]+[MIDNSHP=X]/)) {
                   n=substr(c,1,RLENGTH-1)+0; op=substr(c,RLENGTH,1);
                   if (op ~ /[M=X]/) { if (r<=s && r+n>=e) ok=1; r+=n }
                   else if (op ~ /[DN]/) r+=n;
                   c=substr(c,RLENGTH+1) }
                 if (ok) k++ }
               END{print k+0}')
    # Solo gene-counted reads: the strand call lives or dies here, not on the
    # aligned count. Aligned reads are strand-independent; counted ones are not.
    MTX="${PREFIX}Solo.out/Gene/raw/matrix.mtx"
    SOLO=$(awk 'NR>3{s+=$3} END{print s+0}' "${MTX}" 2>/dev/null || echo 0)

    # Did the barcode read survive? If --clip5pNbases reached the barcode file,
    # CB would be read from base 42 and nothing would match the whitelist. This
    # turns that into a loud failure instead of a quiet zero that looks like a
    # strand result.
    BCSTATS="${PREFIX}Solo.out/Barcodes.stats"
    if [[ -f "${BCSTATS}" ]]; then
        echo "  --- Barcodes.stats ---"
        sed 's/^/    /' "${BCSTATS}"
        # STAR 2.7.11b field names are yesWLmatch* / no*. An earlier version of
        # this guard looked for nExactMatch / nNoMatch, which exist in older
        # STAR, so it matched nothing and silently printed no percentage. The
        # raw stats dump above is what saved it. Match on the yes/no prefixes
        # instead, which are stable across versions.
        CBOK=$(awk '$1 ~ /^yes/ {s+=$2} END{print s+0}' "${BCSTATS}")
        CBNO=$(awk '$1 ~ /^no/  {s+=$2} END{print s+0}' "${BCSTATS}")
        TOTCB=$((CBOK + CBNO))
        if [[ ${TOTCB} -gt 0 ]]; then
            PCT=$(awk -v a=${CBOK} -v t=${TOTCB} 'BEGIN{printf "%.1f", 100*a/t}')
            echo "  CB whitelist match: ${PCT}%"
            BAD=$(awk -v p="${PCT}" 'BEGIN{print (p<20)?1:0}')
            if [[ "${BAD}" == "1" ]]; then
                echo "  *** CB match has COLLAPSED. --clip5pNbases is reaching the"
                echo "  *** barcode read. Stop; use the pre-trim FALLBACK below."
            fi
        fi
    else
        echo "  [WARN] no Barcodes.stats; cannot confirm the barcode read is intact"
    fi

    echo "  input ${INPUT}, aligned(primary) ${ALN}, secondary ${SECONDARY},"
    echo "  MAPQ30 ${HQ}, breadth ${BREADTH}%,"
    echo "  cassette-spanning ${SPAN}, Solo-counted ${SOLO}"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${STR}" "${L}" "${NAME}" "${GRP}" "${INPUT}" "${ALN}" "${HQ}" \
        "${BREADTH}" "${SPAN}" "${SOLO}" >> "${SUMMARY}"
  done
done

echo ""
echo "============================================================"
echo "DONE — $(date)"
echo "summary: ${SUMMARY}"
if [[ "${MODE}" == "smoke" ]]; then
cat <<'NOTE'

READ THE SOLO-COUNTED COLUMN, NOT THE ALIGNED COLUMN.
Aligned-read counts are strand-INDEPENDENT, so both settings will look
similar there. The strand call shows up in Solo counting: on R2 the wrong
setting gave 7 counted reads against ~337 for the right one. Whichever
strand wins by a wide margin is the answer; if they are close, something
other than strand is wrong and the full pass should wait.

Then:
  MODE=full STRAND=<winner> sbatch shiv_starsolo_r1_recovery.sh

FALLBACK if STAR refuses the same file twice in --readFilesIn:
  pre-trim R1 once per library, then pass the trimmed file as mate 1 and
  the original as mate 2, dropping --clip5pNbases:

    zcat R1.fastq.gz | awk 'NR%4==2||NR%4==0{$0=substr($0,42)}1' \
      | gzip > R1_cdna.fastq.gz

  (42 = CLIP5 + 1, because awk substr is 1-based. The quality line is
   trimmed with the sequence or STAR rejects the record.)
NOTE
else
cat <<'NOTE'

Next: union these R1 hits with the R2 results on read id, and with
flank_scan_hits.csv. The three instruments see different things:
  R2 alignment   mismatch-tolerant, but blind to R1
  R1 alignment   mismatch-tolerant, the half never aligned before
  flank scan     both reads, alignment-free, but exact-match only
A read found by alignment and not by the scan is most likely a
hypermutant, which is the population the exact-match probe cannot see.
NOTE
fi
echo "============================================================"
