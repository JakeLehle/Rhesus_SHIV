#!/bin/bash
#SBATCH --job-name=shiv_denovo
#SBATCH --output=/master/jlehle/WORKING/LOGS/shiv_denovo_%j.out
#SBATCH --error=/master/jlehle/WORKING/LOGS/shiv_denovo_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --partition=normal

###############################################################################
# shiv_denovo_assembly.sh
#
# Purpose: Build de novo contigs from unmapped reads, then BLAST against
#          SHIVAD8EO to discover the actual in vivo SHIV sequence.
#
# Strategy:
#   Phase A — Characterize the 864 reads that DID map to SHIVAD8EO
#             (mismatch profiles, coverage gaps, per-gene stats)
#   Phase B — Extract ALL unmapped reads from Library H BAM
#   Phase C — De novo assemble unmapped reads (MEGAHIT, handles large pools)
#   Phase D — BLAST all assembled contigs against SHIVAD8EO reference
#   Phase E — Extract and analyze SHIV-like contigs
#
# Target: Library H (SHIV_Necropsy_LN) — 27 SHIV+ cells, 864 mapped reads
#
# Author: Jake Lehle
# Date: May 27, 2026
###############################################################################

set -euo pipefail
source ~/anaconda3/bin/activate

# ─── Activate conda ─────────────────────────────────────────────────────────
eval "$(conda shell.bash hook)"
conda activate sc_pre
export PATH="${CONDA_PREFIX}/bin:${PATH}"

# ─── Configuration ──────────────────────────────────────────────────────────
THREADS=${SLURM_CPUS_PER_TASK:-16}
WORKDIR="/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV"
ANALYSIS_DIR="${WORKDIR}/Analysis"
OUTDIR="${WORKDIR}/shiv_denovo_analysis"

# Library H = SHIV_Necropsy_LN (27 SHIV+ cells, 864 mapped reads)
LIB_NAME="SHIV_Necropsy_LN"
LIB_DIR="${ANALYSIS_DIR}/${LIB_NAME}/outs"
SAMPLE_BAM="${LIB_DIR}/per_sample_outs/${LIB_NAME}/sample_alignments.bam"
UNASSIGNED_BAM="${LIB_DIR}/unassigned_alignments.bam"

# SHIV contig name (confirmed from previous diagnostic)
SHIV_CONTIG="SHIVAD8EO"

# Reference SHIV FASTA — extract from cellranger reference
CR_REF="/master/jlehle/cellranger-10.0.0/Mmul10_SHIVAD8EO_v10"

# ─── Setup ──────────────────────────────────────────────────────────────────
mkdir -p "${OUTDIR}"/{reads,assembly,blast,reports}
cd "${OUTDIR}"

LOGFILE="${OUTDIR}/reports/denovo_assembly.log"
exec > >(tee -a "${LOGFILE}") 2>&1

echo "============================================================"
echo "SHIV De Novo Assembly — $(date)"
echo "Library: ${LIB_NAME}"
echo "Threads: ${THREADS}"
echo "============================================================"

# ─── Dependency check ───────────────────────────────────────────────────────
echo ""
echo "[CHECK] Dependencies..."
echo "  samtools: $(samtools --version | head -1) [$(which samtools)]"
echo "  minimap2: $(minimap2 --version)"
echo "  blastn:   $(blastn -version | head -1)"

# Prefer MEGAHIT for large unmapped read pools (faster, lower memory)
if command -v megahit &>/dev/null; then
    ASSEMBLER="megahit"
    echo "  megahit:  $(megahit --version 2>&1)"
elif command -v spades.py &>/dev/null; then
    ASSEMBLER="spades"
    echo "  spades:   $(spades.py --version 2>&1)"
else
    echo "  ERROR: Need megahit or spades. Install: conda install -c bioconda megahit"
    exit 1
fi

# Verify BAM exists
if [[ ! -f "${SAMPLE_BAM}" ]]; then
    echo "  ERROR: BAM not found: ${SAMPLE_BAM}"
    exit 1
fi
echo "  BAM:      ${SAMPLE_BAM}"

# ─── Get SHIV reference FASTA ──────────────────────────────────────────────
SHIV_REF="${OUTDIR}/SHIVAD8EO_reference.fasta"
if [[ ! -f "${SHIV_REF}" ]]; then
    echo ""
    echo "[PREP] Extracting SHIV reference sequence..."
    CR_FASTA="${CR_REF}/fasta/genome.fa"
    if [[ -f "${CR_FASTA}" ]]; then
        samtools faidx "${CR_FASTA}" "${SHIV_CONTIG}" > "${SHIV_REF}"
        echo "  Extracted ${SHIV_CONTIG} ($(grep -v '^>' "${SHIV_REF}" | wc -c | awk '{print $1}') bp)"
    else
        echo "  ERROR: Cannot find ${CR_FASTA}"
        echo "  Provide SHIV FASTA manually at: ${SHIV_REF}"
        exit 1
    fi
fi

###############################################################################
# PHASE A: Characterize the 864 SHIV-mapped reads
###############################################################################
echo ""
echo "============================================================"
echo "PHASE A: Characterize ${SHIV_CONTIG}-mapped reads"
echo "============================================================"

# Extract SHIV-mapped reads
samtools view -@ "${THREADS}" -b "${SAMPLE_BAM}" "${SHIV_CONTIG}" \
    > "${OUTDIR}/reads/shiv_mapped.bam"
samtools index "${OUTDIR}/reads/shiv_mapped.bam"

MAPPED_COUNT=$(samtools view -c "${OUTDIR}/reads/shiv_mapped.bam")
echo "  Reads on ${SHIV_CONTIG}: ${MAPPED_COUNT}"

# MAPQ distribution
echo ""
echo "  MAPQ distribution:"
samtools view "${OUTDIR}/reads/shiv_mapped.bam" | \
    awk '{print $5}' | sort -n | uniq -c | sort -rn | head -10 | \
    awk '{printf "    MAPQ %3d: %d reads\n", $2, $1}'

# Edit distance (NM tag) distribution
echo ""
echo "  Edit distance (NM tag) distribution:"
NM_DATA=$(samtools view "${OUTDIR}/reads/shiv_mapped.bam" | \
    grep -oP 'NM:i:\K[0-9]+' || true)
if [[ -n "${NM_DATA}" ]]; then
    echo "${NM_DATA}" | sort -n | uniq -c | sort -rn | head -15 | \
        awk '{printf "    NM=%2d: %d reads\n", $2, $1}'
else
    echo "    (No NM tags in BAM — cellranger may not include them)"
    echo "    Checking XS/AS alignment score tags instead..."
    samtools view "${OUTDIR}/reads/shiv_mapped.bam" | head -3 | \
        awk '{for(i=12;i<=NF;i++) print "    "$i}' || true
fi

# Per-base coverage
echo ""
echo "  Per-base coverage across SHIVAD8EO:"
samtools depth -a -r "${SHIV_CONTIG}" "${SAMPLE_BAM}" \
    > "${OUTDIR}/reports/shiv_mapped_coverage.txt"

awk '
BEGIN { min=999999; max=0; sum=0; n=0; zero=0 }
{
    sum += $3; n++
    if ($3 < min) min = $3
    if ($3 > max) max = $3
    if ($3 == 0) zero++
}
END {
    printf "    Positions: %d\n", n
    printf "    Mean depth: %.2f\n", (n>0 ? sum/n : 0)
    printf "    Max depth:  %d\n", max
    printf "    Zero-cov:   %d (%.1f%%)\n", zero, (n>0 ? 100.0*zero/n : 0)
    printf "    Breadth:    %.1f%%\n", (n>0 ? 100.0*(n-zero)/n : 0)
}' "${OUTDIR}/reports/shiv_mapped_coverage.txt"

# Per-gene coverage (using confirmed GTF coordinates)
echo ""
echo "  Per-gene coverage:"
declare -A GENE_COORDS=(
    ["gag"]="765-2297"
    ["pol"]="1940-5122"
    ["vif"]="5052-5696"
    ["vpx"]="5524-5862"
    ["vpr"]="5863-6168"
    ["tat"]="6014-8721"
    ["rev"]="6240-8905"
    ["vpu"]="6332-6577"
    ["env"]="6495-9047"
    ["nef"]="9060-9851"
)

for gene in gag pol vif vpx vpr tat rev vpu env nef; do
    coords="${GENE_COORDS[$gene]}"
    start="${coords%-*}"
    end="${coords#*-}"
    MEAN_COV=$(awk -v s="${start}" -v e="${end}" \
        '$2 >= s && $2 <= e { sum+=$3; n++ }
         END { if(n>0) printf "%.2f", sum/n; else print "0" }' \
        "${OUTDIR}/reports/shiv_mapped_coverage.txt")
    BREADTH=$(awk -v s="${start}" -v e="${end}" \
        '$2 >= s && $2 <= e { n++; if($3>0) cov++ }
         END { if(n>0) printf "%.1f", 100.0*cov/n; else print "0" }' \
        "${OUTDIR}/reports/shiv_mapped_coverage.txt")
    printf "    %-4s (%s): mean=%-6s breadth=%s%%\n" "${gene}" "${coords}" "${MEAN_COV}" "${BREADTH}"
done

# Save mapped reads as FASTQ (for optional later use)
samtools fastq -@ "${THREADS}" "${OUTDIR}/reads/shiv_mapped.bam" \
    > "${OUTDIR}/reads/shiv_mapped.fastq" 2>/dev/null
MAPPED_FQ=$(grep -c "^@" "${OUTDIR}/reads/shiv_mapped.fastq" || echo "0")
echo ""
echo "  SHIV-mapped reads saved: ${MAPPED_FQ}"

###############################################################################
# PHASE B: Extract ALL unmapped reads
###############################################################################
echo ""
echo "============================================================"
echo "PHASE B: Extract unmapped reads"
echo "============================================================"

# From sample BAM
echo "  Extracting unmapped reads from sample BAM..."
samtools view -@ "${THREADS}" -f 4 -F 0x900 "${SAMPLE_BAM}" | \
    awk '{print "@"$1"\n"$10"\n+\n"$11}' \
    > "${OUTDIR}/reads/unmapped_sample.fastq"

UNMAPPED_SAMPLE=$(( $(wc -l < "${OUTDIR}/reads/unmapped_sample.fastq") / 4 ))
echo "  Sample BAM unmapped reads: ${UNMAPPED_SAMPLE}"

# From unassigned BAM (if exists)
UNMAPPED_UNASSIGNED=0
if [[ -f "${UNASSIGNED_BAM}" ]]; then
    echo "  Extracting unmapped reads from unassigned BAM..."
    samtools view -@ "${THREADS}" -f 4 -F 0x900 "${UNASSIGNED_BAM}" | \
        awk '{print "@"$1"\n"$10"\n+\n"$11}' \
        > "${OUTDIR}/reads/unmapped_unassigned.fastq"
    UNMAPPED_UNASSIGNED=$(( $(wc -l < "${OUTDIR}/reads/unmapped_unassigned.fastq") / 4 ))
    echo "  Unassigned BAM unmapped reads: ${UNMAPPED_UNASSIGNED}"
fi

# Combine
cp "${OUTDIR}/reads/unmapped_sample.fastq" "${OUTDIR}/reads/all_unmapped.fastq"
if [[ -f "${OUTDIR}/reads/unmapped_unassigned.fastq" && -s "${OUTDIR}/reads/unmapped_unassigned.fastq" ]]; then
    cat "${OUTDIR}/reads/unmapped_unassigned.fastq" >> "${OUTDIR}/reads/all_unmapped.fastq"
fi

TOTAL_UNMAPPED=$(( $(wc -l < "${OUTDIR}/reads/all_unmapped.fastq") / 4 ))
echo "  Total unmapped reads for assembly: ${TOTAL_UNMAPPED}"

# Quick read length distribution
echo ""
echo "  Read length distribution (unmapped):"
awk 'NR%4==2 {print length($0)}' "${OUTDIR}/reads/all_unmapped.fastq" | \
    sort -n | uniq -c | sort -rn | head -5 | \
    awk '{printf "    %d bp: %d reads\n", $2, $1}'

###############################################################################
# PHASE C: De novo assembly of unmapped reads
###############################################################################
echo ""
echo "============================================================"
echo "PHASE C: De novo assembly (${ASSEMBLER})"
echo "============================================================"

if [[ "${ASSEMBLER}" == "megahit" ]]; then
    echo "  Running MEGAHIT on ${TOTAL_UNMAPPED} unmapped reads..."
    echo "  (This may take 30-90 minutes depending on read count)"

    # Clean previous run if exists
    rm -rf "${OUTDIR}/assembly/megahit_out"

    megahit \
        -r "${OUTDIR}/reads/all_unmapped.fastq" \
        -o "${OUTDIR}/assembly/megahit_out" \
        -t "${THREADS}" \
        --min-contig-len 200 \
        --k-min 21 \
        --k-max 99 \
        --k-step 10 \
        2>&1 | tail -20

    ASSEMBLY="${OUTDIR}/assembly/megahit_out/final.contigs.fa"

elif [[ "${ASSEMBLER}" == "spades" ]]; then
    echo "  Running SPAdes (--meta) on ${TOTAL_UNMAPPED} unmapped reads..."
    echo "  (This may take 1-3 hours)"

    spades.py \
        --meta \
        -s "${OUTDIR}/reads/all_unmapped.fastq" \
        -o "${OUTDIR}/assembly/spades_out" \
        -t "${THREADS}" \
        -m 120 \
        --only-assembler \
        2>&1 | tail -20

    ASSEMBLY="${OUTDIR}/assembly/spades_out/contigs.fasta"
fi

# Assembly stats
if [[ -f "${ASSEMBLY}" ]]; then
    NCONTIGS=$(grep -c "^>" "${ASSEMBLY}" || echo "0")
    echo ""
    echo "  Assembly complete: ${NCONTIGS} contigs"

    awk '/^>/ {if(len>0){n++; sum+=len; if(len>max)max=len; if(len<min||min==0)min=len; lens[n]=len}
              len=0; next}
         {len+=length($0)}
         END {if(len>0){n++; sum+=len; if(len>max)max=len; if(len<min||min==0)min=len; lens[n]=len}
              # N50
              asort(lens); cumul=0
              for(i=n;i>=1;i--){ cumul+=lens[i]; if(cumul>=sum/2){n50=lens[i]; break} }
              printf "  Total contigs:  %d\n", n
              printf "  Total length:   %d bp\n", sum
              printf "  Longest contig: %d bp\n", max
              printf "  Shortest:       %d bp\n", min
              printf "  N50:            %d bp\n", n50
         }' "${ASSEMBLY}"
else
    echo "  ERROR: Assembly failed — no contigs produced"
    echo "  Check ${OUTDIR}/assembly/ for logs"
    exit 1
fi

###############################################################################
# PHASE D: BLAST all contigs against SHIVAD8EO
###############################################################################
echo ""
echo "============================================================"
echo "PHASE D: BLAST all contigs against SHIVAD8EO reference"
echo "============================================================"

# Build BLAST database from reference SHIV
makeblastdb -in "${SHIV_REF}" -dbtype nucl \
    -out "${OUTDIR}/blast/shiv_ref_db" 2>/dev/null

# BLAST all contigs — keep all hits with reasonable e-value
echo "  BLASTing ${NCONTIGS} contigs against SHIVAD8EO..."
blastn \
    -query "${ASSEMBLY}" \
    -db "${OUTDIR}/blast/shiv_ref_db" \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" \
    -evalue 1e-5 \
    -max_target_seqs 1 \
    -num_threads "${THREADS}" \
    -out "${OUTDIR}/blast/all_contigs_vs_shiv.tsv"

NHITS=$(wc -l < "${OUTDIR}/blast/all_contigs_vs_shiv.tsv")
echo "  BLAST hits: ${NHITS}"

if [[ "${NHITS}" -eq 0 ]]; then
    echo ""
    echo "  WARNING: No contigs matched SHIVAD8EO at e-value < 1e-5"
    echo "  Trying more permissive search (e-value 1e-3, word_size 7)..."

    blastn \
        -query "${ASSEMBLY}" \
        -db "${OUTDIR}/blast/shiv_ref_db" \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" \
        -evalue 1e-3 \
        -word_size 7 \
        -max_target_seqs 1 \
        -num_threads "${THREADS}" \
        -out "${OUTDIR}/blast/all_contigs_vs_shiv.tsv"

    NHITS=$(wc -l < "${OUTDIR}/blast/all_contigs_vs_shiv.tsv")
    echo "  Permissive BLAST hits: ${NHITS}"
fi

###############################################################################
# PHASE E: Analyze SHIV-like contigs
###############################################################################
echo ""
echo "============================================================"
echo "PHASE E: Analyze SHIV-like contigs"
echo "============================================================"

if [[ "${NHITS}" -gt 0 ]]; then
    # Sort by alignment length (descending) to see best hits first
    sort -t$'\t' -k4 -rn "${OUTDIR}/blast/all_contigs_vs_shiv.tsv" \
        > "${OUTDIR}/blast/shiv_hits_sorted.tsv"

    echo ""
    echo "  Top SHIV-like contigs (sorted by alignment length):"
    echo "  ────────────────────────────────────────────────────────────────"
    printf "  %-30s %7s %7s %4s %12s %12s\n" "contig" "%ident" "aln_len" "mm" "ref_range" "contig_len"
    echo "  ────────────────────────────────────────────────────────────────"
    head -30 "${OUTDIR}/blast/shiv_hits_sorted.tsv" | \
        awk -F'\t' '{printf "  %-30s %6.1f%% %7d %4d %5d-%-5d %7d\n",
            $1, $3, $4, $5, $9, $10, $13}'

    # Extract SHIV-like contig sequences
    echo ""
    echo "  Extracting SHIV-like contig sequences..."
    SHIV_CONTIGS="${OUTDIR}/blast/shiv_like_contigs.fasta"
    awk -F'\t' '{print $1}' "${OUTDIR}/blast/shiv_hits_sorted.tsv" | \
        sort -u > "${OUTDIR}/blast/shiv_contig_ids.txt"

    N_SHIV_CONTIGS=$(wc -l < "${OUTDIR}/blast/shiv_contig_ids.txt")
    echo "  Unique SHIV-like contigs: ${N_SHIV_CONTIGS}"

    # Pull sequences from assembly
    python3 -c "
import sys
ids = set(open('${OUTDIR}/blast/shiv_contig_ids.txt').read().split())
printing = False
with open('${ASSEMBLY}') as f, open('${SHIV_CONTIGS}', 'w') as out:
    for line in f:
        if line.startswith('>'):
            name = line.split()[0][1:]
            printing = name in ids
        if printing:
            out.write(line)
print(f'  Wrote {len(ids)} contigs to ${SHIV_CONTIGS}')
"

    # Full pairwise alignment for manual inspection
    echo ""
    echo "  Generating full BLAST alignment..."
    blastn \
        -query "${SHIV_CONTIGS}" \
        -db "${OUTDIR}/blast/shiv_ref_db" \
        -evalue 1e-5 \
        -out "${OUTDIR}/blast/shiv_contigs_alignment.txt"

    # Overall divergence summary
    echo ""
    echo "  ════════════════════════════════════════════════════════"
    echo "  DIVERGENCE SUMMARY"
    echo "  ════════════════════════════════════════════════════════"
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
            printf "  Total BLAST HSPs:       %d\n", n
            printf "  Total aligned bp:       %d / 10262 (%.1f%% of ref)\n", total_aln, 100.0*total_aln/10262
            printf "  Total mismatches:       %d (%.2f%%)\n", total_mm, 100.0*total_mm/total_aln
            printf "  Total gaps:             %d (%.2f%%)\n", total_gap, 100.0*total_gap/total_aln
            printf "  Estimated identity:     %.2f%%\n", 100.0*(total_aln-total_mm-total_gap)/total_aln
        } else {
            print "  No significant alignments"
        }
    }' "${OUTDIR}/blast/shiv_hits_sorted.tsv"

    # Coverage of reference by SHIV-like contigs
    echo ""
    echo "  Reference coverage by SHIV-like contigs:"
    awk -F'\t' '
    BEGIN { for(i=1;i<=10262;i++) cov[i]=0 }
    {
        s = ($9 < $10) ? $9 : $10
        e = ($9 > $10) ? $9 : $10
        for(i=s;i<=e;i++) cov[i]=1
    }
    END {
        covered=0; for(i=1;i<=10262;i++) if(cov[i]) covered++
        printf "  Covered positions: %d / 10262 (%.1f%%)\n", covered, 100.0*covered/10262

        # Per-gene coverage
        split("gag,pol,vif,vpx,vpr,tat,rev,vpu,env,nef", genes, ",")
        split("765,1940,5052,5524,5863,6014,6240,6332,6495,9060", starts, ",")
        split("2297,5122,5696,5862,6168,8721,8905,6577,9047,9851", ends, ",")
        for(g=1;g<=10;g++){
            gc=0; gl=ends[g]-starts[g]+1
            for(i=starts[g];i<=ends[g];i++) if(cov[i]) gc++
            printf "    %-4s: %d/%d bp covered (%.0f%%)\n", genes[g], gc, gl, 100.0*gc/gl
        }
    }' "${OUTDIR}/blast/shiv_hits_sorted.tsv"

else
    echo ""
    echo "  No SHIV-like contigs found in assembly."
    echo "  This is unexpected given 864 mapped reads exist."
    echo ""
    echo "  Possible causes:"
    echo "    - Unmapped reads don't contain enough SHIV signal"
    echo "    - Assembly fragmented SHIV reads into very short contigs (<200bp)"
    echo "    - Try lowering --min-contig-len to 100"
fi

###############################################################################
# Summary
###############################################################################
echo ""
echo "============================================================"
echo "FINAL SUMMARY"
echo "============================================================"
echo ""
echo "Library:              ${LIB_NAME} (Necropsy LN)"
echo "SHIV-mapped reads:    ${MAPPED_COUNT}"
echo "Total unmapped reads: ${TOTAL_UNMAPPED}"
echo "Assembly contigs:     ${NCONTIGS}"
echo "SHIV-like contigs:    ${N_SHIV_CONTIGS:-0}"
echo "BLAST hits:           ${NHITS}"
echo ""
echo "KEY OUTPUT FILES:"
echo "  ${OUTDIR}/reports/denovo_assembly.log"
echo "  ${OUTDIR}/reports/shiv_mapped_coverage.txt"
echo "  ${OUTDIR}/blast/shiv_hits_sorted.tsv"
echo "  ${OUTDIR}/blast/shiv_like_contigs.fasta"
echo "  ${OUTDIR}/blast/shiv_contigs_alignment.txt"
echo "  ${ASSEMBLY}"
echo ""
echo "NEXT STEPS:"
echo "  1. Review shiv_contigs_alignment.txt for SNP positions"
echo "  2. If good coverage: build consensus from SHIV-like contigs"
echo "  3. Compare consensus to MN816822.1 — note divergence hotspots"
echo "  4. If warranted: rebuild cellranger ref with consensus, re-run multi"
echo "  5. Run same analysis on Library D and E for cross-animal comparison"
echo "     (Dr. Ling's convergent evolution question)"
echo ""
echo "Done: $(date)"
