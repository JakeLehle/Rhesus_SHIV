#!/bin/bash
#SBATCH -J MKREF_Mmul10_SHIVAD8EO_v11
#SBATCH -o /master/jlehle/WORKING/LOGS/mkref_Mmul10_SHIVAD8EO_v11.o.log
#SBATCH -e /master/jlehle/WORKING/LOGS/mkref_Mmul10_SHIVAD8EO_v11.e.log
#SBATCH -t 7-00:00:00
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 64
#SBATCH --mem 300GB

###############################################################################
# 02 — Build the Cell Ranger reference: Mmul10 + SHIVAD8EO_V4fullN
#
# Project: Rhesus SHIV scRNA-seq (5' GEX + CITE-seq + VDJ)
# Replaces MKREF_Mmul10_SHIVAD8EO.sh, which built v10.
#
# WHAT CHANGED FROM v10, AND WHY
# ------------------------------
# 1. THE SHIV FASTA AND GTF ARE NO LONGER GENERATED HERE. 01_build_shiv_fasta_gtf.py
#    builds and verifies them, and this script consumes them. That deletes the
#    GFF3 / local-GFF / GenBank-flat-file cascade and its two duplicated
#    converters, roughly 430 lines, because that converter IS the bug:
#    `parse_location` collapsed each GenBank join() to its outermost
#    coordinates, inflating *tat* to 2,708bp and *rev* to 2,666bp so they
#    swallowed *vpu* and most of *env*. Measured consequence, 2026-10-30:
#    33.4% of the viral contig sat under two or more features and 11.6% under
#    none, *tat*, *rev* and *vpu* counted zero reads in all three alignment
#    sources, and collapsing to one feature recovers 2.08x / 2.76x / 1.69x in
#    reads across Cell Ranger, the R1 pass and the R2 pass.
#
# 2. THE VIRAL CONTIG NOW CARRIES THE 34bp CASSETTE, fully masked as N at the
#    derived locus 5862. v10 used the unbarcoded parental, so a read spanning
#    the cassette could not align contiguously and the viral read pool was
#    depleted at exactly the locus the barcode arm needs.
#
# 3. THE vpr SANITY CHECK IS GONE. v10 warned if `vpr` was absent from the SHIV
#    GTF. The collapsed GTF has one feature named SHIVAD8EO and no per-ORF
#    features at all, so that check would now fire on a correct build. Removed
#    rather than loosened, so it cannot be misread as a real warning.
#
# 4. CONTIG MISMATCH IS NOW FATAL. v10 printed "This will cause mkref to fail"
#    and carried on. If the FASTA and GTF disagree, stop before mkref.
#
# 5. TRAILING-NEWLINE GUARD BEFORE CONCATENATION. `cat host.fna virus.fna`
#    silently corrupts the viral contig when the host FASTA has no trailing
#    newline, producing "...ACGT>SHIVAD8EO_V4fullN" on one line. v10 never
#    checked. This script normalises and then verifies the viral contig
#    survived the concatenation intact.
#
# 6. VERSION BUMPED TO v11 AND v10 IS NEVER TOUCHED. Every existing
#    Analysis/ output was built against v10 and it is the comparison baseline
#    for the rerun.
#
# 7. THE REFERENCE IS BUILT INTO ${REFDIR}, not into the Cell Ranger install
#    directory. v10 ran mkref from ${CRDIR}, which left the reference inside
#    the software install. CELLRANGER_MULTI_RH_SHIV.sh already prefers
#    ${REFDIR} and falls back to the install copy.
#
# KEPT FROM v10, DELIBERATELY
# ---------------------------
#   - the Mmul10 download, so this runs without a prestaged host genome
#   - the FULL unfiltered host GTF; mkgtf filtering dropped 40,307 genes to
#     29,349, Total Genes Detected 26,361 to 19,210, Median Genes per Cell
#     3,508 to 2,182
#   - the sed fixes for spaces in the host GTF source column, which make
#     mkref silently misparse lines
#   - the check that every transcript has at least one exon (v10 requirement)
#   - 64 threads, to avoid OOM on a large multi-contig reference
#
# AFTER THIS RUNS
# ---------------
#   03_cellranger_multi.sh points at ${REFDIR}/Mmul10_SHIVAD8EO_v11 and adds
#   chemistry,SC5P-PE-v3. Before stage 04, scanpy_analysis_script.py section 3
#   must become  adata.var['shiv'] = adata.var_names == 'SHIVAD8EO'  with the
#   case-insensitive fallback deleted, or host TAT (tyrosine aminotransferase)
#   is tagged as virus and the run does not error.
#
# Usage:  sbatch 02_mkref_Mmul10_SHIVAD8EO.sh
#
# Author: Jake Lehle
# Date:   October 2026
###############################################################################

set -euo pipefail

#--- Configuration -----------------------------------------------------------
WORKDIR="/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV"
REFDIR="/master/jlehle/WORKING/SC/REF"
CRDIR="/master/jlehle/cellranger-10.0.0"

COMBINED_GENOME="Mmul10_SHIVAD8EO_v11"
PREVIOUS_GENOME="Mmul10_SHIVAD8EO_v10"      # never deleted; the baseline

# Viral inputs, produced and verified by 01_build_shiv_fasta_gtf.py
SHIV_DIR="${WORKDIR}/pipeline_results/01_shiv_reference"
SHIV_CONTIG="SHIVAD8EO_V4fullN"
SHIV_FA="${SHIV_DIR}/${SHIV_CONTIG}.fa"
SHIV_GTF="${SHIV_DIR}/${SHIV_CONTIG}.gtf"

# Expected viral geometry. Re-asserted here so 02 cannot run on a stale or
# hand-edited file that 01 never blessed.
EXPECT_SHIV_LEN=10296
EXPECT_CASSETTE_N=34
EXPECT_LOCUS0=5862          # 0-based cassette start

# Host genome
HOST_BASE="GCF_003339765.1_Mmul_10_genomic"
HOST_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/339/765/GCF_003339765.1_Mmul_10"

THREADS=64
MEMGB=200                   # under the 300GB SLURM allocation

OUT_DIR="${WORKDIR}/pipeline_results/02_mkref"
mkdir -p "${OUT_DIR}" "${REFDIR}"

hr () { echo; echo "=============================================================="; echo "$1"; echo "=============================================================="; }
die () { echo; echo "  STOP: $1"; echo "        Nothing downstream was written."; exit 1; }

MANIFEST="${OUT_DIR}/reference_manifest.txt"
: > "${MANIFEST}"
note () { echo "$1" | tee -a "${MANIFEST}"; }

note "02_mkref  $(date -u +%Y-%m-%dT%H:%M:%SZ)"
note "genome    ${COMBINED_GENOME}"


###############################################################################
hr "STEP 1 — validate the viral inputs from 01 (hard stop on any failure)"
###############################################################################

[ -f "${SHIV_FA}" ]  || die "SHIV FASTA not found: ${SHIV_FA}
        Run 01_build_shiv_fasta_gtf.py first."
[ -f "${SHIV_GTF}" ] || die "SHIV GTF not found: ${SHIV_GTF}
        Run 01_build_shiv_fasta_gtf.py first."

echo "  FASTA : ${SHIV_FA}"
echo "  GTF   : ${SHIV_GTF}"

# These checks duplicate 01 on purpose. 01 proves the files were correct when
# written; this proves they are still correct now, which is not the same thing.
python3 - "${SHIV_FA}" "${SHIV_GTF}" "${SHIV_CONTIG}" \
          "${EXPECT_SHIV_LEN}" "${EXPECT_CASSETTE_N}" "${EXPECT_LOCUS0}" <<'PYEOF' \
    || die "viral input validation failed. Re-run 01 and read its verification.csv."
import re, sys, hashlib
fa, gtf, contig, exp_len, exp_n, exp_locus = sys.argv[1:7]
exp_len, exp_n, exp_locus = int(exp_len), int(exp_n), int(exp_locus)

FLANK5 = 'GGAGACCAGGACCTCCTCCTCCTCCCCCTCCAGGACTAGCATAA'
FLANK3 = 'ATGGAAGAAAGACCTCCAGAAAATGAAGGACCACAAAGGGAA'

names, seqs, cur = [], [], []
for line in open(fa):
    if line.startswith('>'):
        if cur:
            seqs.append(''.join(cur)); cur = []
        names.append(line[1:].split()[0])
    else:
        cur.append(line.strip())
if cur:
    seqs.append(''.join(cur))

fail = []
def ck(label, ok, detail=''):
    print(f"  [{'OK  ' if ok else 'FAIL'}] {label}" + (f"   {detail}" if detail else ''))
    if not ok:
        fail.append(label)

ck('FASTA holds exactly one contig', len(names) == 1, str(len(names)))
if len(names) != 1:
    sys.exit(1)
name, seq = names[0], seqs[0].upper()

ck('contig name matches the expected name', name == contig, name)
ck(f'length is {exp_len:,}', len(seq) == exp_len, f'{len(seq):,}')
ck('only A/C/G/T/N present', bool(re.fullmatch(r'[ACGTN]+', seq)))
ck(f'exactly {exp_n} N in the contig', seq.count('N') == exp_n, str(seq.count('N')))
runs = [(m.start(), m.end()) for m in re.finditer(r'N+', seq)]
ck('exactly one N run', len(runs) == 1, str(len(runs)))
if runs:
    s, e = runs[0]
    ck(f'N run starts at {exp_locus}', s == exp_locus, str(s))
    ck(f'N run is {exp_n} bp', e - s == exp_n, str(e - s))
    ck("5' Keele flank immediately precedes it",
       seq[s - len(FLANK5):s] == FLANK5)
    ck("3' Keele flank immediately follows it",
       seq[e:e + len(FLANK3)] == FLANK3)

rows = [l.rstrip('\n').split('\t') for l in open(gtf) if l.strip()]
ck('GTF has exactly 3 lines', len(rows) == 3, str(len(rows)))
ck('GTF lines have 9 fields', all(len(r) == 9 for r in rows))
ck('GTF feature types are gene/transcript/exon',
   [r[2] for r in rows] == ['gene', 'transcript', 'exon'])
ck('GTF seqname matches the FASTA contig', all(r[0] == name for r in rows))
ck('GTF span is 1..FASTA length',
   all(int(r[3]) == 1 and int(r[4]) == len(seq) for r in rows),
   f'{rows[0][3]}..{rows[0][4]} vs 1..{len(seq)}')
ck('no spaces in the GTF source column', all(' ' not in r[1] for r in rows))
ck('gene_name is SHIVAD8EO', all('gene_name "SHIVAD8EO"' in r[8] for r in rows))

print(f"\n  viral contig md5 : {hashlib.md5(seq.encode()).hexdigest()}")
if fail:
    print(f"\n  {len(fail)} check(s) failed: {fail}")
    sys.exit(1)
print("  viral inputs validated.")
PYEOF


note "shiv_fa   ${SHIV_FA}"
note "shiv_gtf  ${SHIV_GTF}"
note "shiv_md5  $(grep -v '^>' "${SHIV_FA}" | tr -d '\n' | md5sum | cut -d' ' -f1)"


###############################################################################
hr "STEP 2 — Mmul10 host genome (downloaded on demand)"
###############################################################################

cd "${REFDIR}"

if [ ! -f "${HOST_BASE}.fna" ]; then
    echo "  downloading host FASTA..."
    wget -q "${HOST_URL}/${HOST_BASE}.fna.gz" || die "host FASTA download failed"
    gunzip "${HOST_BASE}.fna.gz"
    echo "  done."
else
    echo "  host FASTA already present, skipping download."
fi

if [ ! -f "${HOST_BASE}.gtf" ]; then
    echo "  downloading host GTF..."
    wget -q "${HOST_URL}/${HOST_BASE}.gtf.gz" || die "host GTF download failed"
    gunzip "${HOST_BASE}.gtf.gz"
    echo "  done."
else
    echo "  host GTF already present, skipping download."
fi

echo "  host FASTA contigs : $(grep -c '^>' "${HOST_BASE}.fna")"
echo "  host GTF genes     : $(awk -F'\t' '$3=="gene"' "${HOST_BASE}.gtf" | wc -l)"
note "host_fna  ${REFDIR}/${HOST_BASE}.fna"
note "host_gtf  ${REFDIR}/${HOST_BASE}.gtf"

# No mkgtf filtering. Filtering dropped 40,307 genes to 29,349 on the Mtb/SIV
# project, Total Genes Detected 26,361 -> 19,210, Median Genes per Cell
# 3,508 -> 2,182.
echo "  using the FULL unfiltered host GTF (no mkgtf)."


###############################################################################
hr "STEP 3 — patch spaces out of the host GTF source column"
###############################################################################

# Spaces in column 2 make mkref misparse lines silently. These substitutions
# are idempotent, so re-running is harmless. Applied in place to avoid copying
# a multi-GB GTF. The viral GTF uses "custom" and 01 already asserts it has no
# spaces, so only the host file needs this.
for PAT in 'Protein Homology/Protein_Homology' \
           'RefSeq or Coverage/RefSeq_or_Coverage' \
           'Best Coverage/Best_Coverage' \
           'Best RefSeq/Best_RefSeq' \
           'Curated Genomic/Curated_Genomic'; do
    sed -i "s/${PAT}/g" "${HOST_BASE}.gtf"
done

REMAINING=$( { grep -vP '^#' "${HOST_BASE}.gtf" || true; } | awk -F'\t' 'NF>=9 && $2 ~ / /' | wc -l )
echo "  host GTF lines still carrying a space in the source column: ${REMAINING}"
if [ "${REMAINING}" -gt 0 ]; then
    echo "  offending source values:"
    grep -vP '^#' "${HOST_BASE}.gtf" | awk -F'\t' 'NF>=9 && $2 ~ / /{print $2}' \
        | sort -u | head -10 | sed 's/^/    /'
    die "unhandled spaces in the host GTF source column. Add a substitution above."
fi


###############################################################################
hr "STEP 4 — combine host and virus"
###############################################################################

COMBINED_FA="${REFDIR}/${COMBINED_GENOME}.fa"
COMBINED_GTF="${REFDIR}/${COMBINED_GENOME}.gtf"

# A host FASTA with no trailing newline makes `cat host virus` emit
# "...ACGT>SHIVAD8EO_V4fullN" on one line, which destroys the viral contig
# without any error. Normalise first, verify after.
if [ -n "$(tail -c 1 "${HOST_BASE}.fna")" ]; then
    echo "  host FASTA has NO trailing newline; adding one before concatenation."
    printf '\n' >> "${HOST_BASE}.fna"
else
    echo "  host FASTA ends with a newline, safe to concatenate."
fi
if [ -n "$(tail -c 1 "${HOST_BASE}.gtf")" ]; then
    echo "  host GTF has NO trailing newline; adding one."
    printf '\n' >> "${HOST_BASE}.gtf"
fi

cat "${HOST_BASE}.fna" "${SHIV_FA}"  > "${COMBINED_FA}"
cat "${HOST_BASE}.gtf" "${SHIV_GTF}" > "${COMBINED_GTF}"

echo "  combined FASTA : ${COMBINED_FA}"
echo "  combined GTF   : ${COMBINED_GTF}"
echo "  FASTA contigs  : $(grep -c '^>' "${COMBINED_FA}")"
echo "  GTF lines      : $(wc -l < "${COMBINED_GTF}")"
note "combined_fa  ${COMBINED_FA}"
note "combined_gtf ${COMBINED_GTF}"


###############################################################################
hr "STEP 5 — validate the combined reference before mkref"
###############################################################################

# 5a. The viral contig must have survived concatenation byte for byte.
echo "  checking the viral contig inside the combined FASTA..."
python3 - "${COMBINED_FA}" "${SHIV_FA}" "${SHIV_CONTIG}" <<'PYEOF' \
    || die "the viral contig did not survive concatenation."
import sys
comb, shiv, contig = sys.argv[1:4]

def grab(path, want):
    seq, on = [], False
    for line in open(path):
        if line.startswith('>'):
            if on:
                break
            on = (line[1:].split()[0] == want)
            continue
        if on:
            seq.append(line.strip())
    return ''.join(seq).upper()

a, b = grab(comb, contig), grab(shiv, contig)
if not a:
    print(f"  [FAIL] contig {contig} not found in the combined FASTA")
    sys.exit(1)
if a != b:
    print(f"  [FAIL] viral contig differs after concatenation "
          f"({len(a):,} vs {len(b):,} bp)")
    sys.exit(1)
print(f"  [OK  ] viral contig intact in the combined FASTA ({len(a):,} bp, "
      f"{a.count('N')} N)")
PYEOF


# 5b. Every GTF seqname must exist in the FASTA. Fatal, not a warning.
echo "  checking GTF seqnames against FASTA contigs..."
grep '^>' "${COMBINED_FA}" | sed 's/^>//' | cut -d' ' -f1 | sort -u > "${OUT_DIR}/.fa_contigs"
awk -F'\t' '!/^#/ && NF>=9 {print $1}' "${COMBINED_GTF}" | sort -u > "${OUT_DIR}/.gtf_contigs"
MISSING=$(comm -13 "${OUT_DIR}/.fa_contigs" "${OUT_DIR}/.gtf_contigs")
if [ -n "${MISSING}" ]; then
    echo "  GTF seqnames with no matching FASTA contig:"
    echo "${MISSING}" | head -10 | sed 's/^/    /'
    die "FASTA and GTF disagree on contig names. mkref would fail or misbuild."
fi
echo "  [OK  ] every GTF seqname is present in the FASTA"

# 5c. The viral contig must be annotated.
if ! grep -qP "^${SHIV_CONTIG}\t" "${COMBINED_GTF}"; then
    die "the viral contig ${SHIV_CONTIG} has no GTF entry in the combined file."
fi
echo "  [OK  ] viral contig is annotated in the combined GTF"

# 5d. Cell Ranger v10 requires every transcript to carry at least one exon.
echo "  checking every transcript has an exon..."
python3 - "${COMBINED_GTF}" <<'PYEOF' \
    || die "transcripts without exons. mkref will reject this GTF."
import re, sys
path = sys.argv[1]
with_exon, all_tx = set(), set()
for line in open(path):
    if line.startswith('#'):
        continue
    f = line.rstrip('\n').split('\t')
    if len(f) < 9:
        continue
    m = re.search(r'transcript_id "([^"]+)"', f[8])
    if not m:
        continue
    all_tx.add(m.group(1))
    if f[2] == 'exon':
        with_exon.add(m.group(1))
missing = all_tx - with_exon
if missing:
    print(f"  [FAIL] {len(missing):,} transcripts lack an exon entry")
    for t in sorted(missing)[:10]:
        print(f"         {t}")
    sys.exit(1)
print(f"  [OK  ] all {len(all_tx):,} transcripts have an exon")
PYEOF


rm -f "${OUT_DIR}/.fa_contigs" "${OUT_DIR}/.gtf_contigs"
echo
echo "  combined reference validated. Proceeding to mkref."


###############################################################################
hr "STEP 6 — cellranger mkref"
###############################################################################

# Built into REFDIR rather than into the Cell Ranger install directory.
cd "${REFDIR}"

# v10 is the comparison baseline for the rerun and is never removed here.
if [ -d "${REFDIR}/${PREVIOUS_GENOME}" ]; then
    echo "  previous reference ${PREVIOUS_GENOME} is present and will NOT be touched."
fi

if [ -d "${REFDIR}/${COMBINED_GENOME}" ]; then
    echo "  removing a previous build of ${COMBINED_GENOME}..."
    rm -rf "${REFDIR:?}/${COMBINED_GENOME}"
fi

# Cell Ranger also leaves a martian pipestance directory named after the --id
# it derives from --genome. If that survives, a re-run either refuses to start
# or tries to resume the old one. v10's script removed it and this must too.
if [ -d "${REFDIR}/mkref_${COMBINED_GENOME}" ]; then
    echo "  removing the previous mkref pipestance directory..."
    rm -rf "${REFDIR:?}/mkref_${COMBINED_GENOME}"
fi

# --output= needs coreutils 8.21+; fall back rather than abort under set -e
AVAIL_GB=$( { df -BG --output=avail "${REFDIR}" 2>/dev/null || df -BG "${REFDIR}"; } \
            | tail -1 | grep -oE '[0-9]+G' | head -1 | tr -dc '0-9' )
AVAIL_GB=${AVAIL_GB:-0}
echo "  free space in ${REFDIR}: ${AVAIL_GB}GB"
if [ "${AVAIL_GB}" -gt 0 ] && [ "${AVAIL_GB}" -lt 100 ]; then
    echo "  WARN  under 100GB free; mkref on a genome this size may run out."
fi

echo "  threads ${THREADS}, memgb ${MEMGB}"
echo "  starting $(date -u +%Y-%m-%dT%H:%M:%SZ)"

"${CRDIR}/cellranger" mkref \
    --nthreads="${THREADS}" \
    --memgb="${MEMGB}" \
    --genome="${COMBINED_GENOME}" \
    --fasta="${COMBINED_FA}" \
    --genes="${COMBINED_GTF}"

echo "  finished $(date -u +%Y-%m-%dT%H:%M:%SZ)"


###############################################################################
hr "STEP 7 — verify the built reference"
###############################################################################

REF_OUT="${REFDIR}/${COMBINED_GENOME}"
[ -d "${REF_OUT}" ] || die "mkref produced no ${REF_OUT} directory."

for f in reference.json star/SA star/Genome star/chrName.txt; do
    if [ -e "${REF_OUT}/${f}" ]; then
        echo "  [OK  ] ${f}"
    else
        echo "  [FAIL] missing ${f}"
        die "the built reference is incomplete."
    fi
done

if grep -qx "${SHIV_CONTIG}" "${REF_OUT}/star/chrName.txt"; then
    echo "  [OK  ] viral contig ${SHIV_CONTIG} is in the STAR index"
else
    echo "  [FAIL] viral contig not found in star/chrName.txt"
    die "the virus did not make it into the index."
fi

GENES_FILE="${REF_OUT}/genes/genes.gtf.gz"
if [ -f "${GENES_FILE}" ]; then
    N_SHIV=$( { zcat "${GENES_FILE}" || true; } | awk -F'\t' -v c="${SHIV_CONTIG}" '$1==c' | wc -l )
    echo "  [OK  ] viral GTF lines carried into the reference: ${N_SHIV}"
    if [ "${N_SHIV}" -lt 3 ]; then
        die "expected at least 3 viral GTF lines, found ${N_SHIV}."
    fi
    if zcat "${GENES_FILE}" | grep -q 'gene_name "SHIVAD8EO"'; then
        echo "  [OK  ] gene_name SHIVAD8EO present; it will be the matrix var_name"
    else
        echo "  WARN  gene_name SHIVAD8EO not found in the packaged GTF"
    fi
fi

note "reference ${REF_OUT}"
note "built     $(date -u +%Y-%m-%dT%H:%M:%SZ)"

hr "DONE — ${COMBINED_GENOME}"
cat <<EOF

  reference : ${REF_OUT}
  manifest  : ${MANIFEST}

  Baseline ${PREVIOUS_GENOME} was not modified.

  Cell Ranger leaves ${REFDIR}/mkref_${COMBINED_GENOME}/ behind, which
  holds its own copy of the reference under outs/. The reference above
  is the one to use; the pipestance directory is safe to delete once
  this run is confirmed good.

  NEXT, 03_cellranger_multi.sh:
    GEX_REF=${REF_OUT}
    add to [gene-expression]:   chemistry,SC5P-PE-v3
    archive Analysis/SHIV_Necropsy_LN/outs BEFORE the gate run; the current
    script rm -rf's prior output and that directory is the H baseline
    gate on  --array=7,6  (H and G) rather than H alone, so the control
    moves with the signal

  BEFORE stage 04, in scanpy_analysis_script.py section 3, replace the
  ten-ORF name list and its case-insensitive fallback with:

      adata.var['shiv'] = adata.var_names == 'SHIVAD8EO'

  The collapsed GTF means none of gag/pol/vif/vpx/vpr/tat/rev/vpu/env/nef
  exist as features. The fallback would then match host TAT, tyrosine
  aminotransferase, and report it as viral UMI without raising an error.

  ALSO VERIFY IN THE GATE: the collapsed feature is plus-strand only. The
  STARsolo recovery needed soloStrand REVERSE. Check the antisense fraction
  on ${SHIV_CONTIG} before reading the gate result, or a strand mismatch
  will look like a GTF failure.

EOF
exit 0
