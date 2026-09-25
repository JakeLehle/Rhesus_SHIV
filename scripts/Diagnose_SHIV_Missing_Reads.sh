#!/bin/bash
#SBATCH --job-name=shiv_divergence
#SBATCH --output=/master/jlehle/WORKING/LOGS/shiv_divergence_%j.out
#SBATCH --error=/master/jlehle/WORKING/LOGS/shiv_divergence_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --partition=normal

###############################################################################
# diagnose_shiv_divergence.sh
#
# Purpose: Assess SHIV reference divergence by extracting unmapped/SHIV-mapped
#          reads from cellranger multi BAMs, remapping permissively, assembling,
#          and comparing to MN816822.1.
#
# Strategy:
#   1. Extract reads already mapped to SHIV contig (any MAPQ)
#   2. Extract unmapped reads from same BAM
#   3. Remap unmapped reads against SHIV-only reference (minimap2, permissive)
#   4. Combine all SHIV-candidate reads
#   5. De novo assemble with SPAdes (--isolate mode for low-diversity viral)
#   6. BLAST assembled contigs against MN816822.1
#   7. Report divergence statistics
#
# Starting with Library G (SHIV_Necropsy_PBMC) — highest SHIV+ rate (40702)
#
# Author: Jake Lehle
# Date: May 26, 2026
###############################################################################

set -euo pipefail

# ─── Activate conda environment ─────────────────────────────────────────────
# System samtools (0.1.18) has broken libncurses; use conda sc_pre instead
source ~/anaconda3/bin/activate
conda activate sc_pre

# ─── Activate conda environment ─────────────────────────────────────────────
# System samtools at /usr/local/bin is broken (missing libncurses.so.5).
# SLURM jobs can have /usr/local/bin ahead of conda in PATH, so we force
# conda's bin to the front after activation.
eval "$(conda shell.bash hook)"
conda activate sc_pre
export PATH="${CONDA_PREFIX}/bin:${PATH}"
echo "Conda env: ${CONDA_DEFAULT_ENV:-none}"
echo "samtools:  $(which samtools)"

# ─── Configuration ──────────────────────────────────────────────────────────
THREADS=${SLURM_CPUS_PER_TASK:-16}
WORKDIR="/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV"
ANALYSIS_DIR="${WORKDIR}/Analysis"
OUTDIR="${WORKDIR}/shiv_divergence_analysis"

# Library G = SHIV_Necropsy_PBMC (highest SHIV+ rate, contains 40702)
LIB_NAME="SHIV_Necropsy_PBMC"
LIB_DIR="${ANALYSIS_DIR}/${LIB_NAME}/outs"

# BAM files from cellranger multi
SAMPLE_BAM="${LIB_DIR}/per_sample_outs/${LIB_NAME}/sample_alignments.bam"
UNASSIGNED_BAM="${LIB_DIR}/unassigned_alignments.bam"

# Reference SHIV FASTA (the one used for mkref, from GenBank)
SHIV_REF_FASTA="${WORKDIR}/shiv_divergence_analysis/MN816822.1.fasta"

# Full cellranger reference (to get contig name)
CR_REF="${WORKDIR}/../../cellranger-10.0.0/Mmul10_SHIVAD8EO_v10"
# Alternate path if the above doesn't resolve:
# CR_REF="/master/jlehle/cellranger-10.0.0/Mmul10_SHIVAD8EO_v10"

# ─── Setup ──────────────────────────────────────────────────────────────────
mkdir -p "${OUTDIR}"/{reads,assembly,blast,reports}
cd "${OUTDIR}"

LOGFILE="${OUTDIR}/reports/divergence_diagnostic.log"
exec > >(tee -a "${LOGFILE}") 2>&1

echo "============================================================"
echo "SHIV Divergence Diagnostic — $(date)"
echo "Library: ${LIB_NAME}"
echo "Threads: ${THREADS}"
echo "============================================================"

# ─── Dependency check ───────────────────────────────────────────────────────
echo ""
echo "[STEP 0] Checking dependencies..."

MISSING=()
for TOOL in samtools minimap2 blastn makeblastdb; do
    if ! command -v "${TOOL}" &>/dev/null; then
        MISSING+=("${TOOL}")
    fi
done

# Check for SPAdes (spades.py)
if ! command -v spades.py &>/dev/null; then
    # Try MEGAHIT as fallback
    if ! command -v megahit &>/dev/null; then
        MISSING+=("spades.py OR megahit")
        ASSEMBLER="none"
    else
        ASSEMBLER="megahit"
        echo "  SPAdes not found, will use MEGAHIT"
    fi
else
    ASSEMBLER="spades"
fi

if [[ ${#MISSING[@]} -gt 0 ]]; then
    echo ""
    echo "ERROR: Missing tools: ${MISSING[*]}"
    echo ""
    echo "Install with conda (in sc_pre or a new env):"
    echo "  conda install -c bioconda samtools minimap2 blast spades"
    exit 1
fi

echo "  samtools: $(samtools --version | head -1) [$(which samtools)]"
echo "  minimap2: $(minimap2 --version)"
echo "  blastn:   $(blastn -version | head -1)"
echo "  Assembler: ${ASSEMBLER}"

# ─── Step 1: Identify SHIV contig name ──────────────────────────────────────
echo ""
echo "[STEP 1] Identifying SHIV contig name from BAM header..."

# Check which BAM exists
if [[ -f "${SAMPLE_BAM}" ]]; then
    PRIMARY_BAM="${SAMPLE_BAM}"
    echo "  Using sample BAM: ${SAMPLE_BAM}"
elif [[ -f "${UNASSIGNED_BAM}" ]]; then
    # unassigned_alignments.bam shares the same reference header
    PRIMARY_BAM="${UNASSIGNED_BAM}"
    echo "  WARNING: sample_alignments.bam not found at expected path"
    echo "  Using unassigned BAM for header: ${PRIMARY_BAM}"
    echo "  (Will need to adjust read extraction steps)"
else
    echo "  ERROR: Cannot find BAM file. Checking available files..."
    find "${LIB_DIR}" -name "*.bam" -type f 2>/dev/null | head -20
    echo ""
    echo "  Please update SAMPLE_BAM path and re-run."
    exit 1
fi

# Extract SHIV contig name from BAM header
# Look for MN816822 or SHIV in the @SQ lines
SHIV_CONTIG=$(samtools view -H "${PRIMARY_BAM}" | \
    grep "^@SQ" | \
    awk -F'\t' '{for(i=1;i<=NF;i++) if($i ~ /^SN:/) print substr($i,4)}' | \
    grep -iE "MN816822|SHIV|shiv" | head -1)

if [[ -z "${SHIV_CONTIG}" ]]; then
    echo "  WARNING: Could not find SHIV contig by name. Listing all contigs..."
    samtools view -H "${PRIMARY_BAM}" | grep "^@SQ" | \
        awk -F'\t' '{for(i=1;i<=NF;i++) if($i ~ /^SN:/) print substr($i,4)}' | \
        tail -20  # Show the last 20 (SHIV would be appended after host chroms)
    echo ""
    echo "  Please identify the SHIV contig name from the list above,"
    echo "  set SHIV_CONTIG manually, and re-run."
    exit 1
fi

echo "  SHIV contig name: ${SHIV_CONTIG}"

# Get SHIV contig length from header
SHIV_LEN=$(samtools view -H "${PRIMARY_BAM}" | \
    grep "SN:${SHIV_CONTIG}" | \
    awk -F'\t' '{for(i=1;i<=NF;i++) if($i ~ /^LN:/) print substr($i,4)}')
echo "  SHIV contig length: ${SHIV_LEN} bp"

# ─── Step 2: Extract SHIV-mapped reads ──────────────────────────────────────
echo ""
echo "[STEP 2] Extracting reads mapped to SHIV contig..."

# All reads mapped to SHIV contig (any MAPQ)
samtools view -@ "${THREADS}" -b "${PRIMARY_BAM}" "${SHIV_CONTIG}" \
    > "${OUTDIR}/reads/shiv_mapped_all.bam"

SHIV_MAPPED=$(samtools view -c "${OUTDIR}/reads/shiv_mapped_all.bam")
echo "  Total reads mapped to SHIV contig: ${SHIV_MAPPED}"

# Breakdown by MAPQ
echo "  MAPQ breakdown:"
for Q in 0 1 10 20 30 60 255; do
    COUNT=$(samtools view -c -q ${Q} "${OUTDIR}/reads/shiv_mapped_all.bam")
    echo "    MAPQ >= ${Q}: ${COUNT}"
done

# Extract as FASTQ for assembly input
samtools view -@ "${THREADS}" -F 0x900 "${OUTDIR}/reads/shiv_mapped_all.bam" | \
    awk '{print "@"$1"\n"$10"\n+\n"$11}' \
    > "${OUTDIR}/reads/shiv_mapped.fastq"

SHIV_MAPPED_READS=$(grep -c "^@" "${OUTDIR}/reads/shiv_mapped.fastq" || echo "0")
echo "  Primary SHIV-mapped reads (for assembly): ${SHIV_MAPPED_READS}"

# ─── Step 3: Check SHIV mapping quality details ────────────────────────────
echo ""
echo "[STEP 3] Analyzing SHIV-mapped read quality..."

# Coverage across SHIV genome
samtools depth -a -r "${SHIV_CONTIG}" "${PRIMARY_BAM}" \
    > "${OUTDIR}/reports/shiv_coverage.txt"

# Summary stats
if [[ -s "${OUTDIR}/reports/shiv_coverage.txt" ]]; then
    awk '
    BEGIN { min=999999; max=0; sum=0; n=0; zero=0 }
    {
        sum += $3; n++
        if ($3 < min) min = $3
        if ($3 > max) max = $3
        if ($3 == 0) zero++
    }
    END {
        printf "  Coverage stats across %d positions:\n", n
        printf "    Mean depth: %.2f\n", (n>0 ? sum/n : 0)
        printf "    Min depth:  %d\n", min
        printf "    Max depth:  %d\n", max
        printf "    Zero-cov positions: %d (%.1f%%)\n", zero, (n>0 ? 100.0*zero/n : 0)
        printf "    Breadth (>=1x): %.1f%%\n", (n>0 ? 100.0*(n-zero)/n : 0)
    }' "${OUTDIR}/reports/shiv_coverage.txt"
else
    echo "  WARNING: No coverage data generated"
fi

# Mismatch analysis — NM tag from SHIV-mapped reads
echo ""
echo "  Edit distance (NM tag) distribution for SHIV-mapped reads:"
samtools view "${OUTDIR}/reads/shiv_mapped_all.bam" | \
    grep -oP 'NM:i:\K[0-9]+' | \
    sort -n | uniq -c | sort -rn | head -20 | \
    awk '{printf "    NM=%d: %d reads\n", $2, $1}'

# ─── Step 4: Extract unmapped reads ────────────────────────────────────────
echo ""
echo "[STEP 4] Extracting unmapped reads from BAM..."

# Flag 4 = unmapped
samtools view -@ "${THREADS}" -b -f 4 "${PRIMARY_BAM}" \
    > "${OUTDIR}/reads/unmapped.bam"

UNMAPPED_TOTAL=$(samtools view -c "${OUTDIR}/reads/unmapped.bam")
echo "  Total unmapped reads: ${UNMAPPED_TOTAL}"

# Convert to FASTQ (exclude secondary/supplementary)
samtools view -@ "${THREADS}" -F 0x900 "${OUTDIR}/reads/unmapped.bam" | \
    awk '{print "@"$1"\n"$10"\n+\n"$11}' \
    > "${OUTDIR}/reads/unmapped.fastq"

UNMAPPED_READS=$(grep -c "^@" "${OUTDIR}/reads/unmapped.fastq" || echo "0")
echo "  Unmapped reads (primary only, for remapping): ${UNMAPPED_READS}"

# Also check unassigned_alignments.bam if it exists
if [[ -f "${UNASSIGNED_BAM}" ]]; then
    echo ""
    echo "  Also checking unassigned_alignments.bam..."
    UNASSIGNED_UNMAPPED=$(samtools view -c -f 4 "${UNASSIGNED_BAM}")
    echo "  Unmapped reads in unassigned BAM: ${UNASSIGNED_UNMAPPED}"

    # Extract unmapped from unassigned too
    samtools view -@ "${THREADS}" -F 0x900 -f 4 "${UNASSIGNED_BAM}" | \
        awk '{print "@"$1"\n"$10"\n+\n"$11}' \
        >> "${OUTDIR}/reads/unmapped.fastq"

    UNMAPPED_READS=$(grep -c "^@" "${OUTDIR}/reads/unmapped.fastq" || echo "0")
    echo "  Combined unmapped reads: ${UNMAPPED_READS}"
fi

# ─── Step 5: Permissive remapping of unmapped reads to SHIV ────────────────
echo ""
echo "[STEP 5] Permissive remapping of unmapped reads against SHIV reference..."

# First, get the SHIV reference FASTA
# Extract from the cellranger reference, or download from GenBank
if [[ -f "${SHIV_REF_FASTA}" ]]; then
    echo "  Using existing SHIV reference: ${SHIV_REF_FASTA}"
elif [[ -d "${CR_REF}/fasta" ]]; then
    echo "  Extracting SHIV sequence from cellranger reference..."
    # The cellranger reference genome.fa has all contigs concatenated
    # Use samtools faidx to extract just the SHIV contig
    CR_FASTA="${CR_REF}/fasta/genome.fa"
    if [[ -f "${CR_FASTA}" ]]; then
        samtools faidx "${CR_FASTA}" "${SHIV_CONTIG}" > "${SHIV_REF_FASTA}"
        echo "  Extracted ${SHIV_CONTIG} to ${SHIV_REF_FASTA}"
    fi
else
    echo "  Downloading MN816822.1 from NCBI..."
    # Fallback: download directly
    mkdir -p "$(dirname "${SHIV_REF_FASTA}")"
    curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=MN816822.1&rettype=fasta&retmode=text" \
        > "${SHIV_REF_FASTA}" 2>/dev/null || {
        echo "  ERROR: Cannot download SHIV reference. Please provide ${SHIV_REF_FASTA}"
        echo "  You can extract it from the cellranger reference:"
        echo "    samtools faidx ${CR_REF}/fasta/genome.fa ${SHIV_CONTIG} > ${SHIV_REF_FASTA}"
        exit 1
    }
fi

# Index for minimap2
echo "  Building minimap2 index..."
minimap2 -d "${OUTDIR}/reads/shiv_ref.mmi" "${SHIV_REF_FASTA}" 2>/dev/null

# Permissive mapping: short reads mode, lower identity threshold
# -k10: shorter k-mer (default 15 for sr); catches more divergent reads
# -w5: smaller window
# -A1 -B1: reduced mismatch penalty (more permissive)
# --secondary=no: keep it clean
echo "  Remapping ${UNMAPPED_READS} unmapped reads against SHIV (permissive)..."
minimap2 -a -x sr \
    -k 10 -w 5 \
    -A 1 -B 1 -O 2,16 -E 1,0 \
    --secondary=no \
    -t "${THREADS}" \
    "${SHIV_REF_FASTA}" \
    "${OUTDIR}/reads/unmapped.fastq" 2>/dev/null | \
    samtools view -@ 4 -b -F 4 - | \
    samtools sort -@ 4 -o "${OUTDIR}/reads/unmapped_remapped_shiv.bam"

samtools index "${OUTDIR}/reads/unmapped_remapped_shiv.bam"

RECOVERED=$(samtools view -c "${OUTDIR}/reads/unmapped_remapped_shiv.bam")
echo "  Recovered SHIV-like reads from unmapped pool: ${RECOVERED}"

# Quality of recovered reads
if [[ "${RECOVERED}" -gt 0 ]]; then
    echo "  MAPQ distribution of recovered reads:"
    for Q in 0 1 10 20 30 60; do
        COUNT=$(samtools view -c -q ${Q} "${OUTDIR}/reads/unmapped_remapped_shiv.bam")
        echo "    MAPQ >= ${Q}: ${COUNT}"
    done

    # Coverage of recovered reads
    samtools depth -a "${OUTDIR}/reads/unmapped_remapped_shiv.bam" \
        > "${OUTDIR}/reports/recovered_coverage.txt"

    awk '
    BEGIN { sum=0; n=0; zero=0 }
    {
        sum += $3; n++
        if ($3 == 0) zero++
    }
    END {
        printf "  Recovered read coverage: mean=%.2f, breadth=%.1f%%\n",
            (n>0 ? sum/n : 0), (n>0 ? 100.0*(n-zero)/n : 0)
    }' "${OUTDIR}/reports/recovered_coverage.txt"

    # Extract recovered reads as FASTQ
    samtools view -F 0x900 "${OUTDIR}/reads/unmapped_remapped_shiv.bam" | \
        awk '{print "@"$1"_recovered\n"$10"\n+\n"$11}' \
        > "${OUTDIR}/reads/recovered_shiv.fastq"
fi

# ─── Step 6: Combine all SHIV reads and assemble ───────────────────────────
echo ""
echo "[STEP 6] Combining SHIV reads and running de novo assembly..."

cat "${OUTDIR}/reads/shiv_mapped.fastq" \
    "${OUTDIR}/reads/recovered_shiv.fastq" 2>/dev/null \
    > "${OUTDIR}/reads/all_shiv_reads.fastq"

TOTAL_SHIV_READS=$(grep -c "^@" "${OUTDIR}/reads/all_shiv_reads.fastq" || echo "0")
echo "  Total SHIV-candidate reads for assembly: ${TOTAL_SHIV_READS}"
echo "    - From original mapping: ${SHIV_MAPPED_READS}"
echo "    - Recovered from unmapped: ${RECOVERED:-0}"

if [[ "${TOTAL_SHIV_READS}" -lt 10 ]]; then
    echo ""
    echo "  WARNING: Very few SHIV reads (${TOTAL_SHIV_READS}). Assembly unlikely to succeed."
    echo "  Consider: (a) checking more libraries, (b) relaxing mapping further,"
    echo "  or (c) this may indicate the actual virus is very divergent from MN816822.1."
    echo ""
    echo "  Skipping assembly, but continuing with BLAST of individual reads."
fi

if [[ "${ASSEMBLER}" == "spades" && "${TOTAL_SHIV_READS}" -ge 10 ]]; then
    echo "  Running SPAdes (--isolate mode)..."
    spades.py \
        --isolate \
        -s "${OUTDIR}/reads/all_shiv_reads.fastq" \
        -o "${OUTDIR}/assembly/spades_out" \
        -t "${THREADS}" \
        -m 32 \
        --only-assembler \
        2>&1 | tail -20

    if [[ -f "${OUTDIR}/assembly/spades_out/contigs.fasta" ]]; then
        NCONTIGS=$(grep -c "^>" "${OUTDIR}/assembly/spades_out/contigs.fasta" || echo "0")
        echo "  SPAdes produced ${NCONTIGS} contigs"
        ASSEMBLY="${OUTDIR}/assembly/spades_out/contigs.fasta"

        # Basic assembly stats
        awk '/^>/ {if(len>0){n++;sum+=len; if(len>max)max=len; if(len<min||min==0)min=len} len=0; next}
             {len+=length($0)}
             END {if(len>0){n++;sum+=len; if(len>max)max=len; if(len<min||min==0)min=len}
                  printf "  Assembly stats: %d contigs, total %d bp, longest %d bp, shortest %d bp\n", n, sum, max, min
             }' "${ASSEMBLY}"
    else
        echo "  WARNING: SPAdes did not produce contigs"
        ASSEMBLY=""
    fi

elif [[ "${ASSEMBLER}" == "megahit" && "${TOTAL_SHIV_READS}" -ge 10 ]]; then
    echo "  Running MEGAHIT..."
    rm -rf "${OUTDIR}/assembly/megahit_out"  # MEGAHIT needs clean dir
    megahit \
        -r "${OUTDIR}/reads/all_shiv_reads.fastq" \
        -o "${OUTDIR}/assembly/megahit_out" \
        -t "${THREADS}" \
        --min-contig-len 200 \
        2>&1 | tail -10

    if [[ -f "${OUTDIR}/assembly/megahit_out/final.contigs.fa" ]]; then
        NCONTIGS=$(grep -c "^>" "${OUTDIR}/assembly/megahit_out/final.contigs.fa" || echo "0")
        echo "  MEGAHIT produced ${NCONTIGS} contigs"
        ASSEMBLY="${OUTDIR}/assembly/megahit_out/final.contigs.fa"
    else
        echo "  WARNING: MEGAHIT did not produce contigs"
        ASSEMBLY=""
    fi
else
    ASSEMBLY=""
fi

# ─── Step 7: BLAST against MN816822.1 ──────────────────────────────────────
echo ""
echo "[STEP 7] BLAST analysis against MN816822.1..."

# Build BLAST database from the reference SHIV
makeblastdb -in "${SHIV_REF_FASTA}" -dbtype nucl \
    -out "${OUTDIR}/blast/shiv_ref_db" 2>/dev/null

# BLAST assembled contigs if we have them
if [[ -n "${ASSEMBLY}" && -f "${ASSEMBLY}" ]]; then
    echo "  BLASTing assembled contigs against MN816822.1..."
    blastn \
        -query "${ASSEMBLY}" \
        -db "${OUTDIR}/blast/shiv_ref_db" \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" \
        -evalue 1e-5 \
        -max_target_seqs 1 \
        -num_threads "${THREADS}" \
        -out "${OUTDIR}/blast/contigs_vs_ref.tsv"

    echo ""
    echo "  BLAST results (contigs vs MN816822.1):"
    echo "  qseqid | %ident | aln_len | mismatches | q_start-end | s_start-end | qlen"
    echo "  -------------------------------------------------------------------"
    awk -F'\t' '{printf "  %-30s | %6.2f%% | %6d | %4d | %d-%d | %d-%d | %d\n",
        $1, $3, $4, $5, $7, $8, $9, $10, $13}' \
        "${OUTDIR}/blast/contigs_vs_ref.tsv" | head -30

    # Overall divergence summary
    echo ""
    echo "  === DIVERGENCE SUMMARY ==="
    awk -F'\t' '
    BEGIN { total_aln=0; total_mm=0; total_gap=0; n=0 }
    {
        total_aln += $4
        total_mm += $5
        total_gap += $6
        n++
    }
    END {
        if (total_aln > 0) {
            printf "  %d BLAST hits across %d aligned bp\n", n, total_aln
            printf "  Total mismatches: %d (%.2f%%)\n", total_mm, 100.0*total_mm/total_aln
            printf "  Total gaps: %d (%.2f%%)\n", total_gap, 100.0*total_gap/total_aln
            printf "  Overall identity: %.2f%%\n", 100.0*(total_aln-total_mm-total_gap)/total_aln
        } else {
            print "  No significant BLAST hits — virus may be highly divergent"
        }
    }' "${OUTDIR}/blast/contigs_vs_ref.tsv"

    # Also produce a full alignment for visual inspection
    echo ""
    echo "  Generating pairwise alignment for manual inspection..."
    blastn \
        -query "${ASSEMBLY}" \
        -db "${OUTDIR}/blast/shiv_ref_db" \
        -evalue 1e-5 \
        -max_target_seqs 1 \
        -num_threads "${THREADS}" \
        -out "${OUTDIR}/blast/contigs_vs_ref_alignment.txt"
    echo "  Full alignment: ${OUTDIR}/blast/contigs_vs_ref_alignment.txt"
fi

# ─── Step 8: Per-gene coverage analysis ────────────────────────────────────
echo ""
echo "[STEP 8] Per-gene coverage analysis..."

# SHIV genes from the GTF (approximate coordinates for SHIVAD8-EO MN816822.1)
# These are based on the GenBank annotation — verify against your GTF
cat > "${OUTDIR}/reports/shiv_genes.bed" << 'GENEBED'
MN816822.1	790	2310	gag	0	+
MN816822.1	2088	5174	pol	0	+
MN816822.1	5061	5531	vif	0	+
MN816822.1	5475	5756	vpx	0	+
MN816822.1	5822	6106	vpr	0	+
MN816822.1	6057	6270	tat_exon1	0	+
MN816822.1	6057	6270	rev_exon1	0	+
MN816822.1	6201	6270	vpu	0	+
MN816822.1	6313	8876	env	0	+
MN816822.1	8628	9232	nef	0	+
GENEBED

echo "  NOTE: Gene coordinates above are approximate for MN816822.1."
echo "  Verify against your actual GTF if numbers look off."
echo ""

# If the contig name in the BAM differs from MN816822.1, adjust the BED
if [[ "${SHIV_CONTIG}" != "MN816822.1" ]]; then
    sed -i "s/MN816822.1/${SHIV_CONTIG}/g" "${OUTDIR}/reports/shiv_genes.bed"
    echo "  (Adjusted BED contig name to: ${SHIV_CONTIG})"
fi

# Per-gene coverage from original mapping
if [[ -s "${OUTDIR}/reports/shiv_coverage.txt" ]]; then
    echo ""
    echo "  Per-gene mean coverage (original mapping):"
    while IFS=$'\t' read -r chrom start end gene score strand; do
        MEAN_COV=$(awk -v s="${start}" -v e="${end}" \
            '$2 >= s && $2 <= e { sum+=$3; n++ }
             END { if(n>0) printf "%.1f", sum/n; else print "0" }' \
            "${OUTDIR}/reports/shiv_coverage.txt")
        printf "    %-12s (%d-%d): mean depth = %s\n" "${gene}" "${start}" "${end}" "${MEAN_COV}"
    done < "${OUTDIR}/reports/shiv_genes.bed"
fi

# ─── Step 9: Summary report ────────────────────────────────────────────────
echo ""
echo "============================================================"
echo "DIAGNOSTIC SUMMARY"
echo "============================================================"
echo ""
echo "Library:                    ${LIB_NAME}"
echo "Primary BAM:                ${PRIMARY_BAM}"
echo "SHIV contig:                ${SHIV_CONTIG} (${SHIV_LEN} bp)"
echo ""
echo "Reads mapped to SHIV:       ${SHIV_MAPPED}"
echo "Unmapped reads (total):     ${UNMAPPED_TOTAL}"
echo "Recovered from unmapped:    ${RECOVERED:-0}"
echo "Total SHIV reads:           ${TOTAL_SHIV_READS}"
echo ""
if [[ -n "${ASSEMBLY}" ]]; then
    echo "Assembly:                   ${ASSEMBLY}"
    echo "Assembly contigs:           ${NCONTIGS}"
fi
echo ""
echo "Output directory:           ${OUTDIR}"
echo ""
echo "KEY FILES:"
echo "  ${OUTDIR}/reports/divergence_diagnostic.log   — This log"
echo "  ${OUTDIR}/reports/shiv_coverage.txt           — Per-base coverage"
echo "  ${OUTDIR}/blast/contigs_vs_ref.tsv            — BLAST tabular results"
echo "  ${OUTDIR}/blast/contigs_vs_ref_alignment.txt  — Full BLAST alignment"
echo "  ${OUTDIR}/reads/all_shiv_reads.fastq          — Combined SHIV reads"
if [[ -n "${ASSEMBLY}" ]]; then
    echo "  ${ASSEMBLY}                                    — Assembled contigs"
fi
echo ""
echo "NEXT STEPS:"
echo "  1. Review divergence summary — if >2-3% divergent, a custom ref is warranted"
echo "  2. If assembly is good, inspect contigs_vs_ref_alignment.txt for SNP hotspots"
echo "  3. Check if env has highest divergence (expected for SHIV in vivo evolution)"
echo "  4. If custom ref justified: build consensus, remake cellranger ref, re-run multi"
echo "  5. Consider running this on additional libraries (E=21DPI_PBMC) for comparison"
echo ""
echo "Done: $(date)"

