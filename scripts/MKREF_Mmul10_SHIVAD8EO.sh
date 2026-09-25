#!/bin/bash
#SBATCH -J MKREF_Mmul10_SHIVAD8EO
#SBATCH -o /master/jlehle/WORKING/LOGS/mkref_Mmul10_SHIVAD8EO.o.log
#SBATCH -e /master/jlehle/WORKING/LOGS/mkref_Mmul10_SHIVAD8EO.e.log
#SBATCH -t 7-00:00:00
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 64
#SBATCH --mem 300GB

###############################################################################
# Build custom Cell Ranger reference: Mmul10 + SHIVAD8-EO
#
# Project:  Rhesus SHIV infection scRNA-seq (5' GEX + VDJ + CITE-seq)
# Virus:    SHIVAD8-EO (barcoded = SHIVAD8EOM; barcode between vpx/vpr)
# Adapted from: mkref_Mmul10_MtbCDC1551_SIVmac239_v10.sh
#
# Key lessons carried forward from the Mtb/SIV project:
#   - Use FULL unfiltered Rhesus GTF (mkgtf filtering dropped 11K genes)
#   - Every transcript needs >=1 exon entry (cellranger v10 requirement)
#   - Spaces in GTF source column break mkref — sed-replace them
#   - Verify FASTA contig names match GTF seqnames before mkref
#   - Use --expect-cells=10000 (not --force-cells) in downstream count
#   - Reduced thread count to 64 to avoid OOM with large multi-contig refs
#
# NOTE on barcodes: The 34bp SHIV barcodes in the 5'UTR of vpr do NOT
# need to be in the reference. Cell Ranger maps reads to the coding
# regions (including vpr). Barcode detection happens post-alignment by
# pulling virus-positive cell reads and matching the 5'UTR barcode
# sequences from Binhua's barcode catalog.
###############################################################################

set -euo pipefail

#--- Configuration ---#
REFDIR="/master/jlehle/WORKING/SC/REF"
CRDIR="/master/jlehle/cellranger-10.0.0"
COMBINED_GENOME="Mmul10_SHIVAD8EO_v10"
THREADS=64

###############################################################################
# SHIVAD8-EO GENOME SEQUENCE
#
# GenBank accession MN816822.1 — full proviral sequence of SHIV-AD8EO
# (SIVmac239 backbone + HIV-1 AD8 Env, CCR5-tropic)
# Reference: Gautam et al., PNAS 2012; Nishimura et al., J Virol 2010
#
# The FASTA is downloaded via NCBI efetch and the header is renamed
# to >SHIVAD8EO for clean contig naming in the Cell Ranger reference.
# The GTF is auto-generated from the NCBI GFF3 annotation.
#
# If you need to override with a local file instead, set SHIV_ACCESSION=""
# and place SHIVAD8EO.fna (header: >SHIVAD8EO) in ${REFDIR}/.
###############################################################################

SHIV_ACCESSION="MN816822.1"

mkdir -p "${REFDIR}"
cd "${REFDIR}"

###############################################################################
# 1. Download Rhesus Macaque (Mmul_10) genome and GTF from NCBI
###############################################################################
echo "=== Downloading Rhesus Macaque Mmul_10 ==="

if [ ! -f GCF_003339765.1_Mmul_10_genomic.fna ]; then
    wget -q https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/339/765/GCF_003339765.1_Mmul_10/GCF_003339765.1_Mmul_10_genomic.fna.gz
    gunzip GCF_003339765.1_Mmul_10_genomic.fna.gz
    echo "Rhesus FASTA downloaded."
else
    echo "Rhesus FASTA already exists, skipping."
fi

if [ ! -f GCF_003339765.1_Mmul_10_genomic.gtf ]; then
    wget -q https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/339/765/GCF_003339765.1_Mmul_10/GCF_003339765.1_Mmul_10_genomic.gtf.gz
    gunzip GCF_003339765.1_Mmul_10_genomic.gtf.gz
    echo "Rhesus GTF downloaded."
else
    echo "Rhesus GTF already exists, skipping."
fi

###############################################################################
# 2. Download or verify SHIVAD8-EO genome
###############################################################################
echo "=== Setting up SHIVAD8-EO genome ==="

if [ ! -f SHIVAD8EO.fna ]; then
    if [ -n "${SHIV_ACCESSION}" ]; then
        echo "Downloading SHIVAD8-EO from NCBI accession: ${SHIV_ACCESSION}"
        wget -q -O SHIVAD8EO_raw.fna \
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=${SHIV_ACCESSION}&rettype=fasta&retmode=text"

        # Validate we got a real FASTA
        if ! grep -q "^>" SHIVAD8EO_raw.fna; then
            echo "ERROR: Downloaded file does not appear to be a valid FASTA."
            echo "Check that SHIV_ACCESSION=${SHIV_ACCESSION} is correct."
            cat SHIVAD8EO_raw.fna
            exit 1
        fi

        # Rename the FASTA header to a clean contig name
        # (GenBank headers are long; Cell Ranger needs a simple name)
        awk '
        /^>/ {
            if (!done) {
                print ">SHIVAD8EO"
                done = 1
            }
            next
        }
        { print }
        ' SHIVAD8EO_raw.fna > SHIVAD8EO.fna
        rm -f SHIVAD8EO_raw.fna

        echo "SHIVAD8-EO FASTA downloaded and header cleaned."
        echo "Sequence length: $(grep -v '^>' SHIVAD8EO.fna | tr -d '\n' | wc -c) bp"
    else
        echo ""
        echo "================================================================="
        echo "  SHIVAD8-EO FASTA NOT FOUND AND NO ACCESSION PROVIDED"
        echo ""
        echo "  Either:"
        echo "    1. Set SHIV_ACCESSION at the top of this script, or"
        echo "    2. Place SHIVAD8EO.fna in ${REFDIR}/"
        echo "       (header must be: >SHIVAD8EO)"
        echo "================================================================="
        echo ""
        exit 1
    fi
else
    echo "SHIVAD8-EO FASTA already exists."
    echo "Sequence length: $(grep -v '^>' SHIVAD8EO.fna | tr -d '\n' | wc -c) bp"

    # Verify header is clean
    HEADER=$(head -1 SHIVAD8EO.fna)
    if [[ "${HEADER}" != ">SHIVAD8EO"* ]]; then
        echo "WARNING: FASTA header is '${HEADER}'"
        echo "Expected '>SHIVAD8EO'. Fixing..."
        sed -i '1s/^>.*/>SHIVAD8EO/' SHIVAD8EO.fna
        echo "Header fixed."
    fi
fi

###############################################################################
# 3. Generate SHIVAD8-EO GTF from GenBank annotation
#
#    SHIV is a chimeric virus — there will not be a standard NCBI GTF.
#    We build the GTF from the GenBank annotation.
#
#    Strategy (tries in order):
#      1. Download GFF3 from NCBI sviewer and convert
#      2. Use a local .gff/.gff3 file if present
#      3. Download GenBank flat file and parse gene coordinates
#
#    SHIV gene structure (SIVmac239 backbone + HIV-1 AD8 Env):
#      gag, pol, vif, vpx, vpr, tat, rev, vpu, env, nef
#    The barcode region (5'UTR of vpr, between vpx/vpr) is non-coding
#    and does NOT need GTF annotation — reads map to vpr CDS.
###############################################################################
echo "=== Setting up SHIVAD8-EO GTF ==="

if [ ! -f SHIVAD8EO.gtf ]; then
    # Try to get GFF from NCBI
    if [ -n "${SHIV_ACCESSION}" ]; then
        echo "Downloading GFF3 annotation for ${SHIV_ACCESSION}..."
        wget -q -O SHIVAD8EO_raw.gff3 \
            "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=${SHIV_ACCESSION}&report=gff3"

        if [ -s SHIVAD8EO_raw.gff3 ] && grep -q "gene" SHIVAD8EO_raw.gff3; then
            echo "GFF3 downloaded, converting to Cell Ranger-compatible GTF..."
        else
            echo "WARNING: GFF3 download failed or empty. Will generate manual GTF."
            rm -f SHIVAD8EO_raw.gff3
        fi
    fi

    if [ -f SHIVAD8EO_raw.gff3 ]; then
        # Convert GFF3 to Cell Ranger GTF
        # Handles the standard NCBI GFF3 format for viral genomes
        python3 << 'PYEOF'
import re
import sys

input_gff = "SHIVAD8EO_raw.gff3"
output_gtf = "SHIVAD8EO.gtf"
contig_name = "SHIVAD8EO"

print(f"Converting {input_gff} -> {output_gtf}")

genes = []  # list of (gene_id, gene_name, start, end, strand, feature_type)

with open(input_gff, 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        if len(fields) < 9:
            continue

        seqid, source, ftype, start, end, score, strand, phase, attrs = fields

        if ftype not in ('gene', 'CDS'):
            continue

        # Parse GFF3 attributes
        attr_dict = {}
        for a in attrs.split(';'):
            if '=' in a:
                k, v = a.split('=', 1)
                attr_dict[k.strip()] = v.strip()

        gene_id = attr_dict.get('ID', '')
        gene_name = attr_dict.get('gene', '') or attr_dict.get('Name', '') or attr_dict.get('locus_tag', '')

        if not gene_name and not gene_id:
            continue

        if not gene_name:
            gene_name = gene_id
        if not gene_id:
            gene_id = gene_name

        # Clean gene_id to remove GFF3 prefixes like "gene-" or "cds-"
        gene_id_clean = re.sub(r'^(gene|cds|rna)-', '', gene_id)
        if not gene_id_clean:
            gene_id_clean = gene_name

        genes.append({
            'ftype': ftype,
            'start': start,
            'end': end,
            'strand': strand,
            'gene_id': gene_id_clean,
            'gene_name': gene_name,
        })

# Deduplicate: for each gene_name, keep the widest span
gene_spans = {}
for g in genes:
    name = g['gene_name']
    if name not in gene_spans:
        gene_spans[name] = g.copy()
    else:
        # Expand coordinates to widest
        existing = gene_spans[name]
        if int(g['start']) < int(existing['start']):
            existing['start'] = g['start']
        if int(g['end']) > int(existing['end']):
            existing['end'] = g['end']

print(f"Found {len(gene_spans)} unique genes: {', '.join(sorted(gene_spans.keys()))}")

# Write Cell Ranger compatible GTF
# Each gene needs: gene, transcript, and exon entries (v10 requirement)
with open(output_gtf, 'w') as out:
    for gname in sorted(gene_spans.keys(), key=lambda x: int(gene_spans[x]['start'])):
        g = gene_spans[gname]
        gid = g['gene_id']
        tid = gid + "_t"

        gene_attrs = f'gene_id "{gid}"; gene_name "{gname}";'
        tx_attrs   = f'gene_id "{gid}"; transcript_id "{tid}"; gene_name "{gname}";'
        exon_attrs = f'gene_id "{gid}"; transcript_id "{tid}"; gene_name "{gname}"; exon_number "1";'

        out.write(f'{contig_name}\tSHIVAD8EO\tgene\t{g["start"]}\t{g["end"]}\t.\t{g["strand"]}\t.\t{gene_attrs}\n')
        out.write(f'{contig_name}\tSHIVAD8EO\ttranscript\t{g["start"]}\t{g["end"]}\t.\t{g["strand"]}\t.\t{tx_attrs}\n')
        out.write(f'{contig_name}\tSHIVAD8EO\texon\t{g["start"]}\t{g["end"]}\t.\t{g["strand"]}\t.\t{exon_attrs}\n')

print(f"Written: {output_gtf}")
PYEOF

        rm -f SHIVAD8EO_raw.gff3

    elif [ -f SHIVAD8EO.gff ] || [ -f SHIVAD8EO.gff3 ]; then
        # Local GFF file found — convert it
        LOCAL_GFF=$(ls SHIVAD8EO.gff SHIVAD8EO.gff3 2>/dev/null | head -1)
        echo "Found local GFF: ${LOCAL_GFF}, converting..."
        cp "${LOCAL_GFF}" SHIVAD8EO_raw.gff3
        # Re-run the same python converter above
        python3 << 'PYEOF2'
# (same converter as above, duplicated for the local file path)
import re

input_gff = "SHIVAD8EO_raw.gff3"
output_gtf = "SHIVAD8EO.gtf"
contig_name = "SHIVAD8EO"

genes = []
with open(input_gff, 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        if len(fields) < 9:
            continue
        seqid, source, ftype, start, end, score, strand, phase, attrs = fields
        if ftype not in ('gene', 'CDS'):
            continue
        attr_dict = {}
        for a in attrs.split(';'):
            if '=' in a:
                k, v = a.split('=', 1)
                attr_dict[k.strip()] = v.strip()
        gene_id = attr_dict.get('ID', '')
        gene_name = attr_dict.get('gene', '') or attr_dict.get('Name', '') or attr_dict.get('locus_tag', '')
        if not gene_name and not gene_id:
            continue
        if not gene_name: gene_name = gene_id
        if not gene_id: gene_id = gene_name
        gene_id_clean = re.sub(r'^(gene|cds|rna)-', '', gene_id)
        if not gene_id_clean: gene_id_clean = gene_name
        genes.append({'ftype': ftype, 'start': start, 'end': end,
                       'strand': strand, 'gene_id': gene_id_clean, 'gene_name': gene_name})

gene_spans = {}
for g in genes:
    name = g['gene_name']
    if name not in gene_spans:
        gene_spans[name] = g.copy()
    else:
        existing = gene_spans[name]
        if int(g['start']) < int(existing['start']): existing['start'] = g['start']
        if int(g['end']) > int(existing['end']): existing['end'] = g['end']

print(f"Found {len(gene_spans)} unique genes: {', '.join(sorted(gene_spans.keys()))}")

with open(output_gtf, 'w') as out:
    for gname in sorted(gene_spans.keys(), key=lambda x: int(gene_spans[x]['start'])):
        g = gene_spans[gname]
        gid = g['gene_id']; tid = gid + "_t"
        gene_attrs = f'gene_id "{gid}"; gene_name "{gname}";'
        tx_attrs   = f'gene_id "{gid}"; transcript_id "{tid}"; gene_name "{gname}";'
        exon_attrs = f'gene_id "{gid}"; transcript_id "{tid}"; gene_name "{gname}"; exon_number "1";'
        out.write(f'{contig_name}\tSHIVAD8EO\tgene\t{g["start"]}\t{g["end"]}\t.\t{g["strand"]}\t.\t{gene_attrs}\n')
        out.write(f'{contig_name}\tSHIVAD8EO\ttranscript\t{g["start"]}\t{g["end"]}\t.\t{g["strand"]}\t.\t{tx_attrs}\n')
        out.write(f'{contig_name}\tSHIVAD8EO\texon\t{g["start"]}\t{g["end"]}\t.\t{g["strand"]}\t.\t{exon_attrs}\n')
print(f"Written: {output_gtf}")
PYEOF2
        rm -f SHIVAD8EO_raw.gff3

    else
        # Last resort: download GenBank flat file and parse gene coordinates
        echo "GFF3 and local GFF not available. Trying GenBank flat file..."
        if [ -n "${SHIV_ACCESSION}" ]; then
            wget -q -O SHIVAD8EO_raw.gb \
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=${SHIV_ACCESSION}&rettype=gb&retmode=text"

            if [ -s SHIVAD8EO_raw.gb ] && grep -q "FEATURES" SHIVAD8EO_raw.gb; then
                echo "GenBank flat file downloaded. Parsing gene annotations..."

                python3 << 'PYEOF3'
import re

gb_file = "SHIVAD8EO_raw.gb"
output_gtf = "SHIVAD8EO.gtf"
contig_name = "SHIVAD8EO"

# Parse FEATURES section of GenBank flat file for gene and CDS entries
in_features = False
genes = {}  # gene_name -> {start, end, strand}
current_feature = None
current_quals = {}
current_location = ""

def parse_location(loc_str):
    """Parse GenBank location string to get start, end, strand."""
    strand = "+"
    if loc_str.startswith("complement("):
        strand = "-"
        loc_str = loc_str.replace("complement(", "").rstrip(")")
    # Handle join() for split genes — take outermost coordinates
    if "join(" in loc_str:
        loc_str = loc_str.replace("join(", "").rstrip(")")
    # Extract all numbers
    nums = [int(x) for x in re.findall(r'(\d+)', loc_str)]
    if len(nums) >= 2:
        return min(nums), max(nums), strand
    return None, None, None

def flush_feature():
    """Process accumulated feature data."""
    if current_feature in ("gene", "CDS") and current_location:
        gene_name = current_quals.get("gene", [""])[0]
        locus_tag = current_quals.get("locus_tag", [""])[0]
        name = gene_name or locus_tag
        if not name:
            return
        start, end, strand = parse_location(current_location)
        if start is None:
            return
        if name not in genes:
            genes[name] = {"start": start, "end": end, "strand": strand}
        else:
            # Expand to widest span (handles split genes like tat, rev)
            genes[name]["start"] = min(genes[name]["start"], start)
            genes[name]["end"] = max(genes[name]["end"], end)

with open(gb_file, "r") as f:
    for line in f:
        # Detect FEATURES section
        if line.startswith("FEATURES"):
            in_features = True
            continue
        if line.startswith("ORIGIN") or line.startswith("CONTIG"):
            flush_feature()
            in_features = False
            continue
        if not in_features:
            continue

        # New feature line (starts at column 6 with feature key)
        if len(line) > 5 and line[5] != " " and line[5] != "\n":
            flush_feature()
            parts = line.strip().split()
            if len(parts) >= 2:
                current_feature = parts[0]
                current_location = parts[1]
                current_quals = {}
            else:
                current_feature = None
        elif line.strip().startswith("/"):
            # Qualifier line
            qual = line.strip().lstrip("/")
            if "=" in qual:
                key, val = qual.split("=", 1)
                val = val.strip('"').strip()
                current_quals.setdefault(key, []).append(val)
        elif current_location and not line.strip().startswith("/"):
            # Continuation of location (multi-line join)
            current_location += line.strip()

# Flush last feature
flush_feature()

if not genes:
    print("ERROR: No genes found in GenBank flat file!")
    print("Check the downloaded file SHIVAD8EO_raw.gb")
    import sys; sys.exit(1)

print(f"Found {len(genes)} genes from GenBank: {', '.join(sorted(genes.keys()))}")

# Write Cell Ranger compatible GTF
with open(output_gtf, "w") as out:
    for gname in sorted(genes.keys(), key=lambda x: genes[x]["start"]):
        g = genes[gname]
        gid = gname
        tid = gid + "_t"

        gene_attrs = f'gene_id "{gid}"; gene_name "{gname}";'
        tx_attrs   = f'gene_id "{gid}"; transcript_id "{tid}"; gene_name "{gname}";'
        exon_attrs = f'gene_id "{gid}"; transcript_id "{tid}"; gene_name "{gname}"; exon_number "1";'

        out.write(f'{contig_name}\tSHIVAD8EO\tgene\t{g["start"]}\t{g["end"]}\t.\t{g["strand"]}\t.\t{gene_attrs}\n')
        out.write(f'{contig_name}\tSHIVAD8EO\ttranscript\t{g["start"]}\t{g["end"]}\t.\t{g["strand"]}\t.\t{tx_attrs}\n')
        out.write(f'{contig_name}\tSHIVAD8EO\texon\t{g["start"]}\t{g["end"]}\t.\t{g["strand"]}\t.\t{exon_attrs}\n')

print(f"Written: {output_gtf}")
PYEOF3
                rm -f SHIVAD8EO_raw.gb
            else
                echo ""
                echo "================================================================="
                echo "  FAILED TO DOWNLOAD GENBANK FLAT FILE"
                echo ""
                echo "  Could not retrieve annotation for ${SHIV_ACCESSION}"
                echo "  from either sviewer GFF3 or efetch GenBank format."
                echo ""
                echo "  Please provide a pre-made GTF at ${REFDIR}/SHIVAD8EO.gtf"
                echo "  with gene/transcript/exon entries for:"
                echo "  gag, pol, vif, vpx, vpr, tat, rev, vpu, env, nef"
                echo "================================================================="
                exit 1
            fi
        else
            echo ""
            echo "================================================================="
            echo "  NO SHIVAD8-EO ANNOTATION FOUND AND NO ACCESSION SET"
            echo ""
            echo "  Please provide one of:"
            echo "    - Set SHIV_ACCESSION in this script"
            echo "    - Place SHIVAD8EO.gff in ${REFDIR}/"
            echo "    - Place SHIVAD8EO.gtf in ${REFDIR}/"
            echo "================================================================="
            exit 1
        fi
    fi
else
    echo "SHIVAD8-EO GTF already exists."
fi

# Validate GTF contents
echo "--- SHIVAD8-EO GTF summary ---"
echo "Genes:       $(awk -F'\t' '$3 == "gene"' SHIVAD8EO.gtf | wc -l)"
echo "Transcripts: $(awk -F'\t' '$3 == "transcript"' SHIVAD8EO.gtf | wc -l)"
echo "Exons:       $(awk -F'\t' '$3 == "exon"' SHIVAD8EO.gtf | wc -l)"
echo "Gene names:  $(awk -F'\t' '$3 == "gene"' SHIVAD8EO.gtf | grep -oP 'gene_name "\K[^"]+' | tr '\n' ', ')"
echo ""

# Sanity check: vpr should be in there (reads mapping to vpr coding region
# is how we'll identify SHIV-positive cells before barcode extraction)
if grep -q 'vpr' SHIVAD8EO.gtf; then
    echo "OK: vpr gene found in GTF (critical for barcode detection strategy)."
else
    echo "WARNING: vpr gene NOT found in GTF!"
    echo "  vpr is required — SHIV barcode reads map to the vpr coding region."
    echo "  The 34bp barcodes sit in the 5'UTR upstream of vpr (between vpx/vpr)."
    echo "  Check the annotation source and verify gene names."
fi

###############################################################################
# 4. NO mkgtf filtering — use full Rhesus GTF
#
#    Lesson from Mtb/SIV project: mkgtf filtering removed pseudogenes,
#    snRNA, snoRNA, miRNA, tRNA, etc. reducing 40,307 -> 29,349 genes.
#    This caused Total Genes Detected to drop from 26,361 to 19,210 and
#    Median Genes per Cell from 3,508 to 2,182.
#    We use the full GTF to match standard processing.
###############################################################################
echo "=== Skipping mkgtf filtering (using full GTF) ==="
echo "Full Rhesus GTF genes: $(awk -F'\t' '$3 == "gene"' GCF_003339765.1_Mmul_10_genomic.gtf | wc -l)"

###############################################################################
# 5. Fix spaces in GTF source columns
#
#    Spaces in the GTF source column (column 2) cause cellranger mkref
#    to silently misparse lines. Known offenders from NCBI Rhesus GTF:
#    "Protein Homology", "Best Coverage", "Best RefSeq", etc.
###############################################################################
echo "=== Cleaning GTF source fields (replacing spaces with underscores) ==="

for GTF in GCF_003339765.1_Mmul_10_genomic.gtf SHIVAD8EO.gtf; do
    echo "Fixing spaces in source fields for ${GTF}..."
    sed -i 's/Protein Homology/Protein_Homology/g' "${GTF}"
    sed -i 's/RefSeq or Coverage/RefSeq_or_Coverage/g' "${GTF}"
    sed -i 's/Best Coverage/Best_Coverage/g' "${GTF}"
    sed -i 's/Best RefSeq/Best_RefSeq/g' "${GTF}"
    sed -i 's/Curated Genomic/Curated_Genomic/g' "${GTF}"
done

echo "--- Checking for remaining spaces in GTF source fields ---"
for GTF in GCF_003339765.1_Mmul_10_genomic.gtf SHIVAD8EO.gtf; do
    SPACE_COUNT=$(grep -vP '^#' "${GTF}" | awk '{
        match($0, /^[^ \t]+[ \t]+/)
        rest = substr($0, RLENGTH+1)
        match(rest, /^[^ \t]*[ \t]/)
        src = substr(rest, 1, RLENGTH-1)
        if (src ~ / /) print
    }' | wc -l)
    echo "${GTF}: ${SPACE_COUNT} lines still have spaces in source column"
done

###############################################################################
# 6. Combine FASTA and GTF files
###############################################################################
echo "=== Combining genomes ==="

COMBINED_FA="${REFDIR}/${COMBINED_GENOME}.fa"
COMBINED_GTF="${REFDIR}/${COMBINED_GENOME}.gtf"

cat GCF_003339765.1_Mmul_10_genomic.fna SHIVAD8EO.fna > "${COMBINED_FA}"
cat GCF_003339765.1_Mmul_10_genomic.gtf SHIVAD8EO.gtf > "${COMBINED_GTF}"

echo "Combined FASTA: ${COMBINED_FA}"
echo "Combined GTF:   ${COMBINED_GTF}"

echo "--- FASTA chromosome/contig counts ---"
grep -c "^>" "${COMBINED_FA}"
echo "--- GTF stats ---"
echo "Rhesus GTF genes: $(awk -F'\t' '$3 == "gene"' GCF_003339765.1_Mmul_10_genomic.gtf | wc -l)"
echo "SHIV GTF genes:   $(awk -F'\t' '$3 == "gene"' SHIVAD8EO.gtf | wc -l)"
echo "Combined lines:   $(wc -l < ${COMBINED_GTF})"

# Verify FASTA contig names match GTF seqnames
echo "--- Verifying FASTA contigs match GTF seqnames ---"
FA_CONTIGS=$(grep "^>" "${COMBINED_FA}" | sed 's/^>//' | cut -d' ' -f1 | sort -u)
GTF_CONTIGS=$(awk -F'\t' '!/^#/{print $1}' "${COMBINED_GTF}" | sort -u)
MISSING=$(comm -13 <(echo "${FA_CONTIGS}") <(echo "${GTF_CONTIGS}"))
if [ -n "${MISSING}" ]; then
    echo "WARNING: GTF has seqnames not found in FASTA:"
    echo "${MISSING}"
    echo ""
    echo "This will cause mkref to fail. Check contig naming."
    echo "FASTA contigs (first 5):"
    echo "${FA_CONTIGS}" | head -5
    echo "GTF contigs (first 5):"
    echo "${GTF_CONTIGS}" | head -5
else
    echo "All GTF seqnames found in FASTA. Good to go."
fi

# Verify no transcripts without exons (cellranger v10 requirement)
echo "--- Checking for transcripts without exons (v10 requirement) ---"
python3 -c "
import re

gtf_path = '${COMBINED_GTF}'
transcripts_with_exon = set()
all_transcript_ids = set()

with open(gtf_path, 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        if len(fields) < 9:
            continue
        m = re.search(r'transcript_id \"([^\"]+)\"', fields[8])
        if not m:
            continue
        tid = m.group(1)
        all_transcript_ids.add(tid)
        if fields[2] == 'exon':
            transcripts_with_exon.add(tid)

missing = all_transcript_ids - transcripts_with_exon
if missing:
    print(f'WARNING: {len(missing)} transcripts still lack exon entries!')
    for tid in sorted(list(missing))[:10]:
        print(f'  {tid}')
    if len(missing) > 10:
        print(f'  ... and {len(missing) - 10} more')
else:
    print(f'OK: All {len(all_transcript_ids)} transcripts have exon entries.')
"

###############################################################################
# 7. Run Cell Ranger v10.0.0 mkref
###############################################################################
echo "=== Running cellranger v10.0.0 mkref ==="

cd "${CRDIR}"
if [ -d "mkref_${COMBINED_GENOME}" ]; then
    echo "Removing previous mkref output directory..."
    rm -rf "mkref_${COMBINED_GENOME}"
fi
if [ -d "${COMBINED_GENOME}" ]; then
    echo "Removing previous reference directory..."
    rm -rf "${COMBINED_GENOME}"
fi

./cellranger mkref \
    --nthreads=${THREADS} \
    --genome="${COMBINED_GENOME}" \
    --fasta="${COMBINED_FA}" \
    --genes="${COMBINED_GTF}"

echo ""
echo "================================================================="
echo "  Reference built successfully: ${COMBINED_GENOME}"
echo "================================================================="
echo ""
echo "Copy to REF directory:"
echo "  cp -r ${CRDIR}/${COMBINED_GENOME} ${REFDIR}/"
echo ""
echo "For cellranger count (GEX), use:"
echo "  --transcriptome=${REFDIR}/${COMBINED_GENOME}"
echo "  --chemistry=auto"
echo "  --expect-cells=10000"
echo ""
echo "For cellranger vdj (TCR), use de novo assembly:"
echo "  cellranger vdj --id=<sample>_vdj \\"
echo "    --fastqs=<vdj_fastq_dir> \\"
echo "    --sample=<sample_prefix>"
echo "  (No --reference needed for de novo; annotate with IgBLAST after)"
echo ""
echo "Downstream barcode detection:"
echo "  1. Identify SHIV+ cells (vpr/env/gag expression > 0)"
echo "  2. Extract aligned reads for those cells from possorted_bam"
echo "  3. Match reads against Binhua barcode catalog (SHIVAD8EOM barcodes)"
echo "  4. Assign barcode lineage to each SHIV+ cell"
echo ""
exit 0
