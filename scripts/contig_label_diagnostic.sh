#!/bin/bash
# %% Cell 1 — parameters
set -o pipefail

DB_ROOT="/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV/viral_id_databases"
RESULTS="${DB_ROOT}/search_results"
BLASTN="${RESULTS}/blastn_contigs_vs_lentivirus.tsv"
DIAMOND="${RESULTS}/diamond_blastx_contigs_vs_lentivirus.tsv"
OUTDIR="${RESULTS}/label_discovery"

mkdir -p "${OUTDIR}"

# stitle is the LAST field in both tables. blastn has 16 cols, diamond has 16.
# We pull the trailing stitle by cutting from the 16th field onward, since the
# title itself contains spaces/commas but the leading 15 cols are clean tabs.

echo "=================================================="
echo "Label discovery across blastn + DIAMOND stitle"
echo "Started: $(date)"
echo "=================================================="

# %% Cell 2 — sanity: confirm column counts match expectation
echo ""
echo "### Column count check (should be 16 for both) ###"
echo "blastn  field count (first row): $(head -1 ${BLASTN} | awk -F'\t' '{print NF}')"
echo "diamond field count (first row): $(head -1 ${DIAMOND} | awk -F'\t' '{print NF}')"

# %% Cell 3 — extract stitle column (field 16) from each
echo ""
echo "### Extracting stitle fields ###"
cut -f16 "${BLASTN}"  > "${OUTDIR}/blastn_stitles.txt"
cut -f16 "${DIAMOND}" > "${OUTDIR}/diamond_stitles.txt"
echo "blastn stitles:  $(wc -l < ${OUTDIR}/blastn_stitles.txt) rows"
echo "diamond stitles: $(wc -l < ${OUTDIR}/diamond_stitles.txt) rows"

# %% Cell 4 — BLASTN: what organisms/labels appear?
# Grab the leading descriptor of each title. blastn titles tend to start with
# the organism phrase, so we look at the first several words to see the vocab.
echo ""
echo "### BLASTN: top 40 leading-3-word phrases in stitle ###"
awk '{print $1, $2, $3}' "${OUTDIR}/blastn_stitles.txt" \
    | sort | uniq -c | sort -rn | head -40

echo ""
echo "### BLASTN: keyword tallies (case-insensitive substring counts) ###"
for kw in "HIV-1" "human immunodeficiency virus 1" "HIV-2" "immunodeficiency virus 2" \
          "simian immunodeficiency" "SHIV" "simian-human" "SIVmac" "mac239" "mac251" \
          "smm" "sooty" "lentivirus" "endogenous"; do
    c=$(grep -icF "${kw}" "${OUTDIR}/blastn_stitles.txt")
    printf "  %-40s %d\n" "${kw}" "${c}"
done

# %% Cell 5 — DIAMOND: what subject proteins appear?
# diamond titles here are dominated by pdb chains and protein names. We want to
# see the protein-description vocabulary so we can separate true viral proteins
# from antibody/Ig/host contamination.
echo ""
echo "### DIAMOND: source-database prefix tally (pdb|, sp|, tr|, gb|, ref|, etc) ###"
awk '{print $1}' "${OUTDIR}/diamond_stitles.txt" \
    | sed -E 's/\|.*//' | sort | uniq -c | sort -rn | head -20

echo ""
echo "### DIAMOND: top 50 full stitle strings by frequency ###"
sort "${OUTDIR}/diamond_stitles.txt" | uniq -c | sort -rn | head -50

echo ""
echo "### DIAMOND: viral-protein keyword tallies ###"
for kw in "Gag" "Pol" "Env" "Nef" "Tat" "Rev" "Vif" "Vpr" "Vpu" "Vpx" \
          "integrase" "protease" "reverse transcriptase" "gp120" "gp41" "gp160" \
          "capsid" "matrix" "polyprotein"; do
    c=$(grep -icF "${kw}" "${OUTDIR}/diamond_stitles.txt")
    printf "  %-28s %d\n" "${kw}" "${c}"
done

echo ""
echo "### DIAMOND: host/antibody contamination keyword tallies ###"
for kw in "Fab" "light chain" "heavy chain" "immunoglobulin" "antibody" \
          "Chain H" "Chain L" "Chain A" "Chain B" "bnAb" "neutralizing" \
          "ADP-ribosylation" "kinase" "receptor" "CD4" "TCR" "MHC"; do
    c=$(grep -icF "${kw}" "${OUTDIR}/diamond_stitles.txt")
    printf "  %-28s %d\n" "${kw}" "${c}"
done

# %% Cell 6 — done
echo ""
echo "=================================================="
echo "Label discovery complete: $(date)"
echo "Inspect the tallies above to design the binning patterns."
echo "Raw stitle lists saved in: ${OUTDIR}"
echo "=================================================="
