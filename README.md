# Rhesus_SHIV pipeline

5' scRNA-seq + CITE-seq + TCR VDJ. 8 libraries (A–H), 7 animals, 3 timepoints,
2 tissues. Scanpy on the Texas Biomed HPC.

Everything in `pipeline/` is a numbered stage meant to be run in order, each
with a `run_<stage>.sh` sbatch submitter beside it.

```
pipeline/  01_build_shiv_fasta_gtf.py   02_mkref_Mmul10_SHIVAD8EO.sh
           03_cellranger_multi.sh       04_ingest_qc.py
           05_embed_annotate.py         06_viral_barcode.py
           07_object_report.py          PIPELINE.md
           run_04_ingest_qc.sh          run_05_embed_annotate.sh
           run_06_viral_barcode.sh      run_07_object_report.sh
scripts/   00_gtf_multiplicity_diagnostic.py   00b_refernce_comparison.py
           03b_gate_review.py                 04b_qc_loss_diagnostic.py
           04c_mito_threshold_sweep.py        05b_viral_and_label_audit.py
           run_03b_gate_evaluation.sh         run_05b_viral_and_label_audit.sh
           (plus the 80 legacy scripts)
```
 Everything in `scripts/` is
a diagnostic (number plus a letter suffix) or a legacy script kept for
reference. Diagnostics are read-only and never write into `pipeline_results/`
paths that a stage owns.

```
WORKDIR = /master/jlehle/WORKING/SC/fastq/Rhesus_SHIV
repo    = $WORKDIR/Rhesus_SHIV
outputs = $WORKDIR/pipeline_results/<stage>/
```

The old `Analysis/` tree is never written to. Stage 03 refuses to start if its
resolved output path falls inside it.

---

## Flow

```
    Brandon's SHIV AD8-EO parental FASTA
                 |
    ┌────────────▼──────────────────────────────────────────────┐
    │ 01  build_shiv_fasta_gtf.py                               │
    │     derive locus, insert 34-N cassette, write FASTA + GTF │
    └────────────┬──────────────────────────────────────────────┘
                 │  SHIVAD8EO_V4fullN.fa  (10,296 bp, 1 contig)
                 │  SHIVAD8EO_V4fullN.gtf (1 gene = whole contig)
                 │
    Mmul10 (downloaded)                                          
                 │                                               
    ┌────────────▼──────────────────────────────────────────────┐
    │ 02  mkref_Mmul10_SHIVAD8EO.sh                             │
    │     cellranger mkref, host + virus combined               │
    └────────────┬──────────────────────────────────────────────┘
                 │  REF/Mmul10_SHIVAD8EO_v11/
                 │
    FASTQ (GEX + ADT + VDJ)                                      
                 │                                               
    ┌────────────▼──────────────────────────────────────────────┐
    │ 03  cellranger_multi.sh                                   │
    │     per-library multi run, SC5P-PE-v3                     │
    └────────────┬──────────────────────────────────────────────┘
                 │  per_sample_outs/*/count/sample_filtered_*
                 │  vdj_t/filtered_contig_annotations.csv
                 │
    ┌────────────▼──────────────────────────────────────────────┐
    │ 04  ingest_qc.py                                          │
    │     merge 8 libraries, attach metadata, QC, doublets      │
    └────────────┬──────────────────────────────────────────────┘
                 │  adata_prefilter.h5ad  127,417 cells  (recovery source)
                 │  adata_qc.h5ad         115,402 cells  (the working object)
                 │
    ┌────────────▼──────────────────────────────────────────────┐
    │ 05  embed_annotate.py                                     │
    │     normalize, HVG, PCA, BBKNN, UMAP, Leiden              │
    │     Tier 1 + T/NK correction                              │
    │     VDJ attach, gamma-delta, CD3-dropout T recovery       │
    │     Tier 2 per-lineage subclustering by signature         │
    └────────────┬──────────────────────────────────────────────┘
                 │  adata_tier1_checkpoint.h5ad  (resume point)
                 │  adata_annotated.h5ad
                 │
                 │
    ┌────────────▼──────────────────────────────────────────────┐
    │ 06  viral_barcode.py                                      │
    │     pull reads aligned to the SHIV contig from the BAMs    │
    │     CIGAR-walk the cassette out of soft clips, R1 and R2   │
    │     collapse mates to molecules, dedup on (CB, UMI)        │
    │     negative floor, per-animal inventory, enrichment call  │
    └────────────┬──────────────────────────────────────────────┘
                 │  barcode_inventory.csv, detection_funnel.csv
                 │
    ┌────────────▼──────────────────────────────────────────────┐
    │ 07  object_report.py                                      │
    │     census every obs column and its provenance            │
    │     cell counts by library, animal, timepoint, tissue,    │
    │       Tier 1 and Tier 2, and the crossings of those       │
    │     state which contrasts the design supports             │
    │     NO FIGURES, by design                                 │
    └────────────┬──────────────────────────────────────────────┘
                 │  07_object_report.md  + one CSV per table
```

---

## Running the whole pipeline

Stages are numbered and run in order. Everything activates `sc_pre` and puts
`$CONDA_PREFIX/bin` first on `PATH`, because the system `samtools` at
`/usr/local/bin` is a broken build that otherwise wins.

```bash
conda activate sc_pre
cd $WORKDIR/Rhesus_SHIV/pipeline

python 01_build_shiv_fasta_gtf.py      # seconds; 30 checks must all pass
bash   02_mkref_Mmul10_SHIVAD8EO.sh    # ~2h, cellranger mkref
bash   03_cellranger_multi.sh          # the long one; per-library multi
sbatch run_04_ingest_qc.sh             # 8 cores, 300GB, ~8h
sbatch run_05_embed_annotate.sh        # 16 cores, 256GB, 4-8h (24h booked)
sbatch run_06_viral_barcode.sh         # 4 cores, 64GB, ~1 minute
sbatch run_07_object_report.sh         # 4 cores, 128GB, ~2 minutes
```

Logs land in `/master/jlehle/WORKING/LOGS/`. Every stage runs `python -u`, so
a log fills as the job runs rather than only at exit, and every submitter
echoes the settings that decide what the run actually did.

**Start here when reading the output:** `pipeline_results/07_object_report/07_object_report.md`.
It describes the finished object, which columns are populated, how the cells
break down every way that matters, and which contrasts the design can
actually support.

**Re-running 05 cheaply.** Set `RESUME_FROM_TIER1 = True` to skip the
embedding and reapply Tier 2 curation in minutes. The checkpoint carries a
recipe stamp of every setting that decides a label, and resuming under a
changed recipe is refused with a diff rather than applied silently.

**Renumbered 2026-11-02.** 05, 06 and 07 of the earlier plan (embed+Tier 1,
VDJ, Tier 2) are now one script, so viral barcode and summary moved from 08
and 09 down to 06 and 07.

---

## 01 — `01_build_shiv_fasta_gtf.py`

Builds the custom viral reference. Four ordered stages inside one script:
derive the locus, build the contig, verify it, emit the GTF.

| | |
|---|---|
| **in** | Brandon's parental SHIV AD8-EO FASTA (10,262 bp) |
| **out** | `SHIVAD8EO_V4fullN.fa`, `SHIVAD8EO_V4fullN.gtf` |

### Settings

```python
CONTIG_NAME           = 'SHIVAD8EO_V4fullN'
GENE_ID               = 'SHIVAD8EO_provirus'
GENE_NAME             = 'SHIVAD8EO'
CASSETTE              = 'N' * 34          # Brandon's spec
EXPECTED_LOCUS        = 5862
EXPECTED_PARENTAL_LEN = 10262
EXPECTED_BUILT_LEN    = 10296             # 10262 + 34
EXPECTED_CORE         = (5874, 5884)
```

### Checks: 30, all must pass

Built md5 `106f756c6f62ac89d39942d18db4d96c`. The independent cross-check
rebuild must come out IDENTICAL.

### The two decisions behind this stage

**One gene spanning the whole contig, not per-ORF features.** Cell Ranger
counts a read overlapping exactly one gene feature; a read overlapping two or
more is ambiguous and discarded; zero is intergenic and discarded. HIV ORFs
overlap heavily, so per-ORF features made most viral reads ambiguous and threw
them away. `00_gtf_multiplicity_diagnostic.py` measured the recovery from
collapsing to one feature: **cellranger 2.08x, starsolo R1 2.76x, starsolo R2
1.69x**. Collapsing cannot manufacture background, because the control
libraries have zero viral reads to begin with.

The GenBank source also had a parsing bug worth recording: `parse_location`
took min/max across all numbers in a `join()`, inflating *tat* to 2,708 bp and
*rev* to 2,666 bp so they swallowed *vpu* and most of *env*. Deriving the locus
independently, as 01 does, sidesteps it.

**34 N at the cassette, not ACGCGC + 22 N + GCGCGT.** The conserved flanks are
real sequence in the virus, so writing them into the reference is tempting.
It is wrong. Real bases in a reference window constrain where the aligner may
open a gap. 32 of the 1,290 catalog entries (2.5%) are shorter at this locus
and align with a deletion; with the flanks present, V2allN demonstrably forced
the ADE08M.11736 deletion into the register tag at position 5868, corrupting
the barcode call. All 34 positions stay N.

`00b_refernce_comparison.py` Cell 7 also settled which reference the existing
BAMs used: V4fullN appears in zero BAMs, V2allN in 10 R1 and 8 R2. That is why
every BAM had to be regenerated rather than reused.

---

## 02 — `02_mkref_Mmul10_SHIVAD8EO.sh`

| | |
|---|---|
| **in** | Mmul10 FASTA + GTF (downloaded by the script), 01's FASTA + GTF |
| **out** | `REF/Mmul10_SHIVAD8EO_v11/` |

The host reference downloads on demand, so this runs with no input staged.
The viral half is custom and must come from 01.

Cut from 675 to roughly 250 lines: the GFF3 / local-GFF / GenBank cascade and
both duplicated converters are gone. One path in, one path out.

### Verified on build

2,940 contigs; 40,307 host genes; 88,316 transcripts, all with exons.
The pipestance directory `mkref_${COMBINED_GENOME}` is removed on success.

### Two bash traps that bit here

`[ $? -eq 0 ] || die` is unreachable under `set -e`: the shell exits on the
failing command before the test runs. `|| die` must be attached to the command
itself. And `grep` returning 1 aborts under `pipefail`, so greps that may
legitimately match nothing are wrapped `{ cmd || true; }`.

---

## 03 — `03_cellranger_multi.sh`

| | |
|---|---|
| **in** | FASTQ for GEX, ADT and VDJ; `REF/Mmul10_SHIVAD8EO_v11`; `feature_reference.csv` |
| **out** | `pipeline_results/03_cellranger_multi/<library>/` |

### Settings

```
[gene-expression]
reference,/master/jlehle/WORKING/SC/REF/Mmul10_SHIVAD8EO_v11
chemistry,SC5P-PE-v3
create-bam,true
expect-cells,10000
```

All 8 libraries completed. Preflight resolves the output path with
`readlink -m` and refuses to run if it lands inside `Analysis/`.

### Correction to the rebuild plan

**The plan's chemistry claim was false and should be struck.** v10 was already
paired-end. The reasoning that produced the error: `multi` never emits a
chemistry row in `metrics_summary.csv`, so the absence of one proved nothing,
and the STARsolo R1 reads that looked skipped were in fact discarded by the
GTF multiplicity rule, not by the aligner. The chemistry line above is
explicit now for the collaborators' benefit, not because it changed anything.

`SESSION_CONTEXT.md` and the rebuild plan still state the chemistry claim as
fact in two places. Both need editing once the rerun is signed off.

Note for anyone reading logs: use the conda `samtools`. `/usr/local/bin/samtools`
on titan is broken and must be skipped.

---

## 04 — `04_ingest_qc.py`

| | |
|---|---|
| **in** | `pipeline_results/03_cellranger_multi/*/per_sample_outs/*/count/` |
| **out** | `adata_prefilter.h5ad`, `adata_qc.h5ad`, `qc_dropped_ledger.csv`, `shiv_whitelist.csv`, `qc_summary.csv` |
| **submit** | `sbatch run_04_ingest_qc.sh` (8 cores, 300GB, 8h) |

### Settings

```python
MT_THRESHOLD           = 20      # fixed. NO MAD on mito.
MAD_COUNTS_GENES       = 5       # two-sided, global, on counts and genes only
MIN_GENES              = 200
MIN_CELLS_PER_GENE     = 20
EXPECTED_DOUBLET_RATE  = 0.08
REMOVE_DOUBLETS        = False   # score only; see below
WHITELIST_VIRAL_CELLS  = True
MIN_WHITELIST_UMI      = 1
REUSE_PREFILTER        = True
SHIV_GENE              = 'SHIVAD8EO'
SHIV_GENE_ID           = 'SHIVAD8EO_provirus'
HOST_COLLISION         = ['TAT']
REQUIRE_SHIV_FEATURE   = True
```

### Cell order

```
CELL 2  metadata (31 per-animal entries, 7 animals)
CELL 3  ingestion          ── Cells 4/5 nested under `if not REUSED:`
CELL 5b verify viral tag   ── UNCONDITIONAL, runs on the reuse path too
CELL 6  prefilter object
CELL 7  QC
CELL 8  Scrublet            ── scoring only; REMOVE_DOUBLETS gates removal
CELL 9  save
```

CELL 5b is unconditional on purpose. An earlier version nested the viral
verification inside `if not REUSED:`, which raised `NameError: n_shiv` on the
reuse path and, worse, **silently skipped the host TAT guard**. 5b re-derives
`n_shiv` from the object, re-runs the TAT assertion, and hard-stops if no
viral feature is tagged.

### The host TAT collision

Rhesus tyrosine aminotransferase is `TAT`; viral *tat* differs only by case.
`HOST_COLLISION = ['TAT']` asserts that the host gene is present in the matrix
and is **not** tagged `shiv`. Verified in the current run: look-alike present,
none tagged.

### Why the mito filter is a fixed 20% and not MAD

`04b_qc_loss_diagnostic.py` Cell 6 is the decisive test: the same comparison
run *within* one library, so library and animal are separated.

| library | Indian removed | Chinese removed | gap |
|---|---|---|---|
| C | 53.5% | 16.9% | +36.6 pts |
| F | 91.4% | 15.6% | +75.8 pts |

The gap survives inside a single library, so it is an animal effect, not a
batch effect. The global 3-MAD bound landed at 5.51% while the Indian control
median is 7.58%, meaning most of the control arm was an outlier *by
construction*.

`04c_mito_threshold_sweep.py` scored 12 candidate rules on the same cells.
`fixed_20%` removes 7.8% Indian / 0.4% Chinese and keeps 117,461 cells, against
92,732 under the old rule. Per-sample MAD was rejected despite looking
symmetric: it is blind to a uniformly degraded sample, since if every cell in a
sample sits at 15% mito none of them is an outlier within it.

One thing to ignore in that sweep: the `asymmetry` column is a ratio of
fractions and is misleading when the absolute removals differ by an order of
magnitude.

### Viral whitelist

60 barcodes carry ≥1 viral UMI. 5 of them would have failed QC and are kept,
flagged `qc_would_have_dropped` so the exemption is auditable rather than
invisible. 59 survive the gene filter.

### Results of the current run

Measured with `REMOVE_DOUBLETS = True`, which is no longer the setting. The
re-run with removal off ends at **116,536 cells** rather than 115,402; every
other number here is unchanged, because doublet removal was the last step.

```
cells loaded (pre-filter)   127,417
after QC (+ doublets)       115,402  x 18,582 genes   -> 116,536 with removal off
mito filter removed           3,137  (2.5%)
  Indian controls             2,253 of 28,918  (7.8%)
  Chinese infected              368 of 91,164  (0.4%)
viral cells pre-filter             60
viral cells surviving              59
total viral UMI                   185
```

Per animal: 41861 22,068 / 41862 17,700 / 40702 17,004 / 40707 16,593 /
34315 15,311 / 41903 10,715 / 39272 10,440 / Unknown 5,571.

Per timepoint: Necropsy 35,474 / 21 DPI 33,542 / Non-Infected 26,026 /
Pre 14,789 / Unknown 5,571.

The control arm went from 7,753 cells to 26,026. That is the whole point of
the mito change.

### GTF fix, measured on library H

| | v10 | v11 | |
|---|---|---|---|
| viral UMI | 38 | 55 | 1.45x |
| SHIV+ cells | 6 | 10 | 1.67x |
| cells ≥1 UMI | 27 | 35 | |

Library G (control) has zero viral reads in both versions. Clean.

### Scrublet is unreliable here, so removal is OFF

```
A 3 (0.0%)   B 5 (0.1%)   C 1 (0.0%)   D 2 (0.0%)
E 0 (0.0%)   F 0 (0.0%)   G 1,121 (5.9%)   H 2 (0.0%)
total 1,134 = 0.97% against EXPECTED_DOUBLET_RATE = 0.08
```

Six of eight libraries return literally zero. Scrublet's automatic threshold is
picking a cut in the tail of a unimodal simulated-doublet distribution, which
is its known failure mode.

**Resolved 2026-11-02: `REMOVE_DOUBLETS = False`.** The first version of 04
hardcoded `adata = adata[~(preds & ~wl)].copy()` with no flag, which deleted
1,134 cells, 1,121 of them from library G alone. Removing 5.9% of one library
and ~0% of six others is not doublet removal, it is a batch artifact written
into the object before the embedding runs. Any threshold substituted in its
place would have been invented rather than measured.

`doublet_score` and `predicted_doublet` travel with the object. 05 handles
doublets at annotation as **evidence**: a subcluster is flagged
`QC: doublet-like` when it co-expresses markers of two lineages
(`QC_DOUBLET_COEXPR_FRAC = 0.50`) **and** carries a mean `doublet_score` at
least `QC_DOUBLET_SCORE_MULT = 1.5` times the lineage median. Flagged, never
deleted, so 06 can still recover a viral barcode from one.

If `04` was already run with removal on, re-run it. `REUSE_PREFILTER = True`
makes that a few minutes, since `adata_prefilter.h5ad` is untouched by this.

---

## 05 — change log

Every change since the first version, with the measurement that forced it.
The script carries the same list in its header. Newest first.

| date | change | what forced it |
|---|---|---|
| 11-03 | HVG exclusions applied to the **within-lineage** HVG, not just global | 48% of T and B cells were in subclusters defined by ribosomal or mitochondrial content |
| 11-03 | ribosomal, translation and MHC class I added to the exclusions | RIRA excludes the same sets; class I polymorphism corrupts counts, ribosomal genes are batch-prone |
| 11-03 | Monocyte `Intermediate` panel **deleted**, DC split into cDC1/cDC2/pDC/mregDC | the panel was an MHC-II gradient; it mislabelled 993 B cells and 142 cDC2 |
| 11-03 | `Effector/Cytotoxic` auto-splits on the measured granzyme skew | the skew separated 10,015 cells into TEM-like and TCM-like on its first run |
| 11-03 | `top_family_range` reported | "Treg, margin 0.005" read as meaningful when that panel's whole range is 0.100 |
| 11-03 | B panels replaced with RIRA modules; `InnateLike` added | *CD27*, *CD24*, *IGHD* absent from Mmul10 and unrecoverable by renaming |
| 11-03 | RIRA granzyme and S100A modules added as diagnostics | they cut across the assignment rather than competing with it |
| 11-03 | `CD3_MIN_COUNTS` 3 → 5 | measured background CD3+ was 12.11% at 3, 3.28% at 5 |
| 11-03 | NK-versus-v10 line stopped comparing across thresholds | it made a threshold change look like a biological shift |
| 11-02 | Tier 2 assignment moved from per-cell vote to subcluster mean | z-scoring gave uninformative panels equal voting power; 10 of 11 T subclusters Unresolved |
| 11-02 | CD3 gate generalized to every non-T label; marker support added | the gate touched only the two labels where the v10 collision was known |
| 11-02 | graded viral call; checkpoint recipe stamp | ≥1 viral UMI is not a detection; a stale checkpoint would silently reapply the old gate |

### The within-lineage HVG bug, 2026-11-03

`CELL 4` excluded `shiv` and `mt`. `CELL 12` excluded **nothing**. So the
global embedding was clean while Tier 2 subclustering was driven by whatever
technical gene varied most inside each lineage. Counting each subcluster's top
15 markers:

| subcluster | ribosomal | KEG06 (mito) | cells |
|---|---|---|---|
| T sub 0 | **13/15** | | 6,400 |
| T sub 4 | **12/15** | | 5,243 |
| T sub 2 | **11/15** | | 3,820 |
| T sub 8 | 2 | **9/15** | 6,749 |
| B sub 9 | **11/15** | | 4,803 |
| B sub 10 | **10/15** | | 6,458 |
| B sub 5 | **9/15** | | 2,732 |
| B sub 2 | 1 | **10/15** | 3,647 |

**39,852 of 83,004 T and B cells, 48%.** That is why the biological panels
scored near zero in those subclusters, and why `Naive/CM` collapsed from 0.290
to 0.061 between the two runs: those subclusters are defined by ribosomal
content and nothing in a naive panel competes with that.

The two mito-driven subclusters are the cleanest case. `KEG06_*` is the rhesus
mitochondrial genome (`KEG06_p01`–`p13` are the 13 protein-coding genes,
`KEG06_r01`/`r02` the two rRNAs). T sub 8 has a median `pct_counts_mt` of
7.46% and B sub 2 has 7.57%, against the Indian control animals' measured
median of **7.58%** (04c), and both subclusters are **100% Indian controls**.
They are the control animals grouped by sample quality.

**04 was never wrong about this.** Line 439 tags `KEG06_*` correctly as `mt`,
so `pct_counts_mt` and the fixed 20% filter were always computed on the right
gene set. Only the Tier 2 HVG step missed them.

One shared `hvg_exclusion()` helper now serves both cells so they cannot drift
apart again.

```python
EXCLUDE_FROM_HVG = ['shiv', 'mt']          # var-column flags from 04
EXCLUDE_HVG_PATTERNS = {
    'ribosomal':        r'^(RP[SL]|MRP[SL])',
    'ribosomal_fusion': r'^(FAU|UBA52)$',
    'translation':      r'^EEF1',
    'mhc_class_I':      r'^MAMU-(A|A3|AG|B|E|F)$',
}
EXCLUDE_MHC_CLASS_II = False
EXCLUDE_HVG_IN_LINEAGE = True
```

**MHC class II is deliberately kept.** RIRA's reason for dropping class I is
that high polymorphism makes the counts inaccurate in a single reference, and
that does not apply to class II. Class II is also genuinely informative here:
*MAMU-DRA* is in the Monocyte panel and MHC-II separates B cells from DCs.
`EXCLUDE_MHC_CLASS_II = True` drops it, at the cost of real signal. *B2M* is
kept for the same reason: it is the class I light chain but sits outside the
polymorphic region. Histones are kept, because that is real proliferation
signal and the `Proliferating` panel depends on it.

**These exclusions apply to HVG selection only.** Every gene stays in the
object and is still scored, plotted and available to `rank_genes_groups`, so
ribosomal and `KEG06_*` genes will still appear in the top-marker tables. They
just no longer define the axes. Verified that no panel gene is caught by any
pattern.

### Monocyte `Intermediate` deleted, 2026-11-03

The *MAMU-DRA* alias resolved and made the panel **worse**. *CD74*, *FCGR3*
and *MAMU-DRA* are all MHC-II pathway or Fc, so nothing in the panel is
monocyte-specific. Range 2.300, winning 7 of 12 subclusters: a lineage-wide
MHC-II gradient wearing a subset label. It mislabelled subcluster 4 (993
cells, top markers *MS4A1 CD79B BANK1 AFF3 FCRL1*, so B cells) at a margin of
0.845, and subcluster 11 (142 cells, *FCER1A* plus high MHC-II plus *CST3*, so
cDC2) at a margin of 2.356.

**RIRA found two monocyte populations in 41,199 rhesus myeloid cells, not
three:** "separation of CD14+ classical monocytes from CD16+ non-classical
monocytes". Forcing a third category with no monocyte-specific gene is what
produced those mislabels.

| panel | genes | note |
|---|---|---|
| Classical | *CD14 VCAN S100A8 S100A9* | RIRA found *VCAN* "more robust and specific" than *CD14* |
| NonClassical | *FCGR3 CX3CR1 ITGAX CDKN1C MS4A7* | RIRA: *CX3CR1* and *ITGAX* "more specific in scRNA-seq data" |
| Macrophage | *CD68 APOC1 C1QA FABP4 MRC1 CD163* | |
| cDC1 | *CLEC9A XCR1 CADM1* | |
| cDC2 | *CD1C FCER1A CLEC10A* | this is what subcluster 11 is |
| pDC | *IL3RA LILRA4* | RIRA also lists *GZMB* and *JCHAIN*; **omitted on purpose**, because Monocyte subcluster 8 is *GZMB*/*NKG7*/*GNLY*-high NK contamination and a *GZMB*-bearing pDC panel would capture it |
| mregDC | *CCR7 LAMP3 CCL22 FSCN1* | RIRA's DC3 / mature migratory DC |

### `Effector/Cytotoxic` sub-splits on the granzyme skew, 2026-11-03

The skew worked on its first run and separated one label into two real
populations:

| sub | GranzymeKM | GranzymeABH | skew | Perforin | S100A memory | tissue | reading |
|---|---|---|---|---|---|---|---|
| 5 | 0.352 | **1.229** | **+0.878** | 0.848 | 1.114 | PBMC +1.13 | TEM / terminal effector |
| 6 | **1.343** | -0.229 | **-1.572** | 0.478 | 0.473 | neutral | TCM / TTM-like |

Top markers confirm independently: subcluster 5 is *GZMB GZMH*, subcluster 6
is *GZMK GZMM*. Both were labelled `Effector/Cytotoxic`, 10,015 cells between
them, and `SUBTYPE_RENAME` could not tell them apart because it is keyed on
the provisional name. The alternative was to let it key on subcluster ID,
which is the thing removed from this script for being unstable across
rebuilds. **A split on a measured quantity is stable.**

`FAMILY_SUBSPLIT` thresholds the skew at ±0.30, well clear of both observed
values so it is not tuned to them. RIRA's mapping from sorted reference cells:
*GZMM* alone is naive CD8, *GZMM*+*GZMK* is TCM/TTM, *GZMB*+*GZMH* is TEM,
*GZMA*+*GZMB*+*GZMH* is CD16-positive NK, *GZMM*+*GZMK* is CD16-negative NK.

## 05 — `05_embed_annotate.py`

Consolidates six legacy scripts: `shiv_embeding.py`,
`Rhesus_T_NK_Cell_Correction.py`, `shiv_VDJ_integration.py`,
`shiv_tier1_tcr_recovery.py`, `shiv_tier2_diagnostic_v2.py` and
`shiv_tier2_annotation_write.py`.

| | |
|---|---|
| **in** | `pipeline_results/04_ingest_qc/adata_qc.h5ad`, plus `vdj_t/filtered_contig_annotations.csv` per library from 03 |
| **out** | `adata_annotated.h5ad`, `adata_tier1_checkpoint.h5ad`, `tier1_transition_matrix.csv`, `tier1_tcr_recovery_matrix.csv`, `tier1_marker_panels.csv`, `tier2_<lin>_subtype_table.csv`, `tier2_<lin>_subcluster_signature_means.csv`, `tier2_<lin>_top_markers.csv`, `tier2_<lin>_doublet_evidence.csv`, `tier2_<lin>_{animal,timepoint,tissue}_enrichment.csv`, `tier2_by_{timepoint,animal_id}.csv`, `adt_detection.csv`, `shiv_positive_cells.csv`, `figures/F05_*.pdf/.png` |
| **submit** | `sbatch run_05_embed_annotate.sh` (16 cores, 256GB, 24h) |

### Why these are one script and not four

Ordering, not convenience. TCR recovery adds CD3-dropout T cells back into the
T pool, and that pool has to be **complete before subclustering**, or the T
subclusters are built on a truncated population and every subtype fraction is
wrong. In the old chain each script wrote an intermediate object whose only
purpose was to be reloaded by the next, and each re-normalized defensively
because none could assume what the others had done. Normalization happens once.

The cost of merging is that iterating on Tier 2 would mean re-running the
embedding. `RESUME_FROM_TIER1 = True` loads `adata_tier1_checkpoint.h5ad` and
skips to Tier 2, finishing in minutes. Same pattern as `REUSE_PREFILTER` in 04.

`CELL 11b` runs unconditionally on both paths and re-derives everything Tier 2
needs from the object. That is the lesson from 04, where nesting a check inside
`if not REUSED:` raised a `NameError` and silently skipped the host *TAT* guard.

### Settings, carried forward unchanged

```python
# global embedding
COUNTS_LAYER     = 'counts';   TARGET_SUM      = 1e4
HVG_FLAVOR       = 'seurat';   N_TOP_HVG       = 3000
EXCLUDE_FROM_HVG = ['shiv', 'mt']        # var-column flags from 04
EXCLUDE_HVG_PATTERNS   = ribosomal, ribosomal_fusion, translation, mhc_class_I
EXCLUDE_HVG_IN_LINEAGE = True            # the 11-03 fix
N_PCS_PCA        = 50;         N_PCS_NEIGHBORS = 30
LEIDEN_RES       = 1.0;        RANDOM_STATE    = 42
BATCH_KEY        = 'library';  DO_PREBBKNN_UMAP = True

# Tier 1 and the T/NK correction
SCORE_THRESHOLD  = 0.1
CD3_MIN_COUNTS   = 3           # RAW counts. Was 1; see the audit below.
CD3_GATE_ALL_LABELS       = True
CD3_CALIBRATOR_LINEAGES   = ['Erythrocytes', 'Platelets']
REQUIRE_MARKER_SUPPORT    = True
MIN_MARKERS_DETECTED      = 1
NK_EXCLUSIVE     = ['KLRD1', 'KLRC3', 'NCAM1']
NK_SHARED        = ['GNLY', 'PRF1', 'GZMB', 'NKG7']     # informational only

# the graded viral call
VIRAL_BACKGROUND_LINEAGES = ['Erythrocytes', 'Platelets', 'Mast cells']
SHIV_CONFIDENT_MIN_UMI    = None      # None = derive the ceiling each run
SHIV_CONFIDENT_MIN_FRAC   = 0.01

# VDJ
GD_SUBTYPES      = ['TRG+TRD']                          # strict, paired
TCR_RECOVERY_QUARANTINE = ['B cells', 'Monocytes', 'Macrophages', 'DCs']

# Tier 2
TIER2_LINEAGES   = ['T cells', 'B cells', 'Monocytes']
MIN_LINEAGE_CELLS = 300
SUB_HVG = 2000;  SUB_PCS = 50;  SUB_NEIGH_PCS = 30;  SUB_RES = 1.0
MIN_SCORE_MARGIN   = 0.15                # replaced MIN_DOMINANT_FRAC
REQUIRE_POSITIVE_TOP = True
QC_MITO_SUBCLUSTER = 15.0
QC_DOUBLET_COEXPR_FRAC = 0.50;  QC_DOUBLET_SCORE_MULT = 1.5
SUBTYPE_RENAME = {}                                     # empty by default
```

Nothing is tuned in this stage. The object changed; the recipe did not.

### Step by step

| cell | step | settings / markers | in → out |
|---|---|---|---|
| 2 | load, or resume | hard-stop if obs/layers are from an older checkpoint | `adata_qc.h5ad` or `adata_tier1_checkpoint.h5ad` |
| 3 | stash raw counts, normalize | hard-stop if X is not raw and `layers['counts']` absent; `normalize_total(1e4)` → `log1p` | X → log-norm X |
| 4 | HVG | `flavor='seurat'`, `n_top_genes=3000`, virus and mito excluded | → `var['highly_variable']` |
| 5 | PCA | `n_comps=50`, `use_highly_variable=True` | → `obsm['X_pca']` |
| 5 | pre-correction UMAP | `n_pcs=30` | → `obsm['X_umap_prebbknn']` |
| 5 | BBKNN, UMAP, Leiden | `batch_key='library'`, `n_pcs=30`, `resolution=1.0`, seed 42 | → `obsm['X_umap']`, `obs['clusters']` |
| 6 | ADT CLR | per-cell `log1p(x) - mean(log1p(x))` | `obsm['protein_counts']` → `obsm['protein_clr']`, `obs['ADT_<tag>']` |
| 7 | Tier 1 baseline | 9 RIRA panels, `score_genes`, argmax, floor 0.1 | → `obs['tier1_celltype_uncorrected']` |
| 8 | T/NK correction | exclusive NK anchors + CD3 hard gate on raw counts | → `obs['tier1_celltype']`, `obs['nk_provisional']` |
| 9 | VDJ attach | scirpy `read_10x_vdj` + `chain_qc`; barcode suffix join; subset AIRR to GEX | → `obs['has_tcr']`, `receptor_subtype`, `chain_pairing`, `is_gamma_delta` |
| 9 | resolve provisional NK | TCR-positive → T; TCR-negative stays provisional | → `obs['nk_provisional']` narrowed |
| 10 | CD3-dropout T recovery | `has_tcr` AND not T AND not quarantined | → `obs['tier1_celltype']`, `obs['tcr_recovered']` |
| 11 | checkpoint | | → `adata_tier1_checkpoint.h5ad` |
| 11b | verify (unconditional) | re-derives `HAS_ADT`, `HAS_TCR`, `HAS_DBL` from the object | |
| 12 | Tier 2 per lineage | `SUB_HVG=2000` within lineage, PCA 50, BBKNN 30, sub-Leiden 1.0, signature z-score argmax, mode per subcluster | → `obs['sub_leiden']`, `tier2_subtype`, `tier2_dominant_frac`, `qc_exclude` |
| 13 | Tier 2 summary | | → `tier2_by_*.csv` |
| 14 | viral cells | | → `shiv_positive_cells.csv` |
| 15 | figures | 300 DPI, PDF + PNG, font 28–34, hex | → `figures/F05_*` |
| 16 | save | | → `adata_annotated.h5ad` |

### `BATCH_KEY = 'library'`, not `animal_id`

04b showed the dominant structure is per animal, which tempts you toward animal
as the batch key. That would delete the experiment. Animal is confounded with
the biology: Indian versus Chinese **is** infection status, and 40707 being the
elite controller is the centrepiece of most contrasts. Correcting on animal
removes exactly the signal the project exists to measure.

Library is the GEM well, the genuine technical batch. Libraries C and F each
carry both subspecies, so library correction does not flatten the subspecies
contrast.

`X_umap_prebbknn` exists so the batch effect can be seen rather than assumed.
If the uncorrected UMAP already mixes, BBKNN is doing little and the batch key
is worth revisiting.

### Virus and mito excluded from HVG

A 3,000-gene HVG set with the viral feature in it would let 59 infected cells
pull on the global embedding. Mitochondrial genes track dissociation stress
rather than cell identity. Both stay in the object and are still scored; they
just do not define the axes.

### Tier 1 marker panels

Nine lineages, `sc.tl.score_genes`, argmax, `Unassigned` below
`SCORE_THRESHOLD = 0.1`. Panels are RIRA (Bimber et al. 2025), rhesus-validated.

| lineage | markers |
|---|---|
| T cells | *CD3D CD3E CD3G* |
| B cells | *CD79A MS4A1 CD19 BANK1 PAX5* |
| NK cells | *KLRD1 KLRC3 GNLY PRF1 GZMB NKG7* (baseline; corrected below) |
| Monocytes | *S100A8 S100A9 VCAN CD14 LYZ* |
| Macrophages | *CD68 CD163 MRC1 MERTK APOC1* |
| DCs | *FLT3 CLEC9A CLEC10A IL3RA CLEC4C* |
| Mast cells | *CPA3 KIT HPGD ENPP3* |
| Erythrocytes | *HBB HBA1 HBA2* |
| Platelets | *PPBP PF4 GP9* |

`tier1_marker_panels.csv` records which markers were actually found. Rhesus
naming drops a few every time, so a quietly empty panel is a real risk and is
logged rather than left to produce a lineage with a flat zero score.

### The T/NK correction

The baseline NK panel is the defect. *GNLY PRF1 GZMB NKG7* is a shared
cytotoxic block expressed by cytotoxic T cells, and under argmax it overwhelms
the three-gene CD3 score. `Rhesus_T_NK_Cell_Diagnostic.py` quantified it on
the v10 object: **62.2% of Tier 1 NK calls were CD3-positive** with a median
NK-minus-T margin of 0.762, and **39.1% of Unassigned cells were CD3-positive
T cells**.

Four moves:

1. Rebuild the NK score on exclusive anchors only, *KLRD1 KLRC3 NCAM1* →
   `dx_nk_excl`.
2. Move the shared block to `dx_cyto_shared`, which **never** participates in
   assignment and exists only to be looked at.
3. Apply the CD3 hard gate at `CD3_MIN_COUNTS = 1` **raw** count. A
   CD3-positive cell cannot be NK, whatever the score says.
4. Recover CD3-positive `Unassigned` cells to T.

Both labels are computed and the full transition matrix is written to
`tier1_transition_matrix.csv`. The old label survives as
`tier1_celltype_uncorrected` so the correction is auditable.

**The CD3 hard gate before cytotoxic scoring is a standing convention on this
project.** Flat winner-take-all scoring leaks lymphocytes across lineages and
will do it again in any rebuilt object if the gate is not reapplied.

### What the 05b audit changed, 2026-11-02

`scripts/05b_viral_and_label_audit.py` is read-only and runs in minutes. It
found two defects in the first version of 05 and refuted a third that was only
a suspicion.

**`CD3_MIN_COUNTS` 1 to 3, calibrated rather than chosen.** At one count,
Erythrocytes were 20.2% CD3-positive and Platelets 13.0%. Neither can express
*CD3*, so ambient *CD3* mRNA alone makes 13 to 20% of **any** lineage
CD3-positive at one count. B cells sat at 19.1%, inside that background range,
which is why the apparent "6,948 CD3-positive B cells" dissolved on inspection.
The threshold sweep:

| CD3 ≥ | CD3+ cells in ungated labels |
|---|---|
| 1 | 9,432 |
| 2 | 3,660 |
| 3 | **1,831** |
| 5 | 561 |

05 now prints the background rate in the calibrator lineages at the chosen
threshold and warns above `CD3_BACKGROUND_TOLERANCE = 2.0%`, so the choice is
visibly a measurement.

**The gate now applies to every non-T label.** It previously touched only
`NK cells` and `Unassigned`, which were where the v10 collision happened to be
known. The gate exists because the scores are untrustworthy, and which label a
bad score produced is not a reason to trust it.

**A label must be backed by at least one of its own markers.** Assignment was
argmax over `score_genes` with a floor on the *score*, not on detection; in a
sparse cell both terms sit near zero and noise picks the winner. Measured cost,
per lineage:

| label | zero-markers | median detected | note |
|---|---|---|---|
| B cells | 0.0% | 3 of 5 | solid; *CD79A* 74.6%, *MS4A1* 77.5% |
| T cells | 10.8% | 2 of 3 | these are the TCR-recovered cells, by design |
| NK cells | 14.3% | 4 of 6 | thin |
| Macrophages | 0.0% | **1 of 5** | *CD68* 78% but *MERTK* 3.7% |
| DCs | 0.0% | **1 of 5** | *FLT3* 78% but *CLEC4C* 1.9% |
| Mast cells | 0.0% | **1 of 4** | *KIT* 77% but *CPA3* 11.9% |

A failing cell becomes `Unassigned`, never a different lineage, because failing
one panel is not evidence for another.

**Ordering matters and is the safeguard.** Marker support runs in CELL 8,
*before* TCR recovery in CELL 10. The recovered cells carry no *CD3* by
construction, which is the entire point of the recovery, and they account for
almost exactly the 10.8% of T cells showing zero markers. Re-applying the
support check after recovery would strip precisely the cells the recovery
exists to rescue, so it is never re-applied, and `marker_support_exempt`
records the exemption for any later stage that re-checks.

**31% B cells was NOT a bug, and the suspicion was mine.** Composition by
library is 7.2% B in A (PBMC) and 51.2% in H (LN), 23% against 39% by tissue.
That is textbook. The hypothesis that B cells were mass-mislabelled is refuted
by `pct_zero_markers = 0.0`.

**The checkpoint is recipe-stamped.** `uns['tier1_recipe']` records every
setting that decides a Tier 1 label. Resuming a checkpoint written under a
different recipe is refused with a diff, rather than silently reapplying the
old gate while Cell 1 claims otherwise.

### VDJ and TCR

scirpy 0.22.5, MuData-era API (`merge_with_tcr` is gone). VDJ barcodes are
`<bare>-1`; GEX barcodes are `<bare>-1-<lib>`, so the library suffix is
appended before joining. AIRR is then subset to cells present in GEX, which
drops the TCR cells whose GEX partner QC removed and keeps the modalities
aligned.

**Gamma-delta is called strictly from a paired `TRG+TRD` receptor_subtype.**
*TRG* persists in alpha-beta T cells and *TRD* recovery is sparse, so "any
*TRG* or *TRD*" overcounted gamma-delta at 2,982 cells. The strict call is a
floor, which is the honest direction to be wrong in. The loose count is printed
for transparency and never used.

**`nk_provisional` is resolved asymmetrically, on purpose.** A provisional NK
carrying a productive TCR is reassigned to T, and that is strong evidence. A
provisional NK with **no** TCR stays provisional, because TCR dropout means
absence is only supporting evidence for NK, never proof.

**CD3-dropout T recovery** flips a cell to T when it carries a TCR, is not
already T, and is not in the quarantine set. The quarantine is
`['B cells', 'Monocytes', 'Macrophages', 'DCs']`, because measured ambient
orphan-TCR rates there were B 9.3%, Monocytes 10.7%, Macrophages 14.4%, so a
lone chain is more likely soup than a receptor. `Unassigned`, `Erythrocytes`
and `Platelets` are **not** quarantined: their own *HBB*/*PPBP* signal is
itself ambient, so a detected receptor is the better evidence. On v10 the
Unassigned pool was 9,832 cells, 0% CD3-positive by RNA but 3,493 (35.5%)
carrying a called TCR.

### Tier 2 is assigned by signature, not by a hardcoded cluster map

`shiv_tier2_annotation_write.py` carried
`SUBTYPE_MAP = {'0': 'CD8 effector memory', '5': ..., '10': ...}`. Leiden
cluster IDs are not stable across a reference change, a cell-count change, or a
seed change, so that map silently mislabels everything after any rebuild, and
this rebuild changed all three.

Instead: z-score each signature across cells, argmax per cell, mode per
subcluster, and report `dominant_frac` so a split subcluster is visible rather
than hidden behind a majority vote. A subcluster below
`MIN_DOMINANT_FRAC = 0.40` is labelled `Unresolved`, which is a real answer and
not a failure: the signatures do not separate it at this resolution.

HVG are recomputed **inside** each lineage, because genes that separate T from
B are useless for separating naive T from effector T.

`rank_genes_groups` (Wilcoxon, top 15) runs alongside so the signature call can
be checked against genes nobody chose in advance.

#### Tier 2 panels

**T cells:** Naive *CCR7 SELL TCF7 LEF1 IL7R* · CentralMem *CCR7 SELL IL7R CD28
GPR183* · EffectorMem *GZMK CCL5 CST7 KLRG1* · Cytotoxic *GZMB GNLY PRF1 NKG7
GZMH FGFBP2* · Th1 *TBX21 IFNG CXCR3 STAT1* · Th17 *RORC IL17A CCR6 IL23R
KLRB1* · Tfh *CXCR5 PDCD1 BCL6 ICOS IL21* · Treg *FOXP3 IL2RA CTLA4 IKZF2* ·
Proliferating *MKI67 TOP2A*

**B cells:** Naive *TCL1A IGHD FCER2 IL4R* · Memory *CD27 TNFRSF13B AIM2
TXNIP* · Transitional *MME VPREB3 SOX4 CD24* · Plasma *MZB1 XBP1 PRDM1 JCHAIN
SDC1* · GerminalCenter *BCL6 AICDA RGS13 MEF2B* · Activated *CD69 CD83 NR4A1*

**Monocytes:** Classical *CD14 S100A8 S100A9 VCAN* · NonClassical *FCGR3
CDKN1C LILRB2 MS4A7* · Intermediate *CD74 FCGR3 HLA-DRA* · Macrophage *CD68
CD163 MRC1 MERTK APOC1 C1QA* · DC *FLT3 CLEC9A CLEC10A IL3RA*

#### Curation is explicit, not hardcoded

`SUBTYPE_RENAME` is **empty by default**, which means the names you see came
from the signatures rather than from a guess. Read the review tables, then
promote provisional names to curated ones by filling it in, keyed
`'<lineage>:<provisional>'`. The v10 adjudication, for reference:

| evidence | curated name |
|---|---|
| Cytotoxic + ~100% ADT-CD8 + LN-excluded + necropsy-high | CD8 terminal effector |
| EffectorMem + 90–95% ADT-CD8 + PBMC-enriched | CD8 effector memory |
| CentralMem + 77–82% ADT-CD4 | CD4 central memory |
| Naive, *LEF1*/*BACH2* | Naive/resting T |

#### Decisions carried as per-cell attributes, never as labels

Each of these was argued out in the v10 pass and lost. They stay per-cell:

- **CD4/CD8**: `cd48_adt` (CD8 CLR > CD4 CLR, the reliable call) and
  `cd48_rna`. The naive pool stays **one** label carrying the attribute.
- **Proliferation**: `is_proliferating` gate on *MKI67*/*TOP2A* plus the
  continuous `sig_Proliferating`. The clean proliferating cluster dissolved on
  the recovered pool, so this is per-cell.
- **CD4 activation**: `sig_Activation` on immediate-early genes. Partly a
  tissue-dissociation artifact, so it is a soft score and **never** a label.

#### QC rules are explicit thresholds, never cluster IDs

- `QC: high-mito` when the subcluster's **median** `pct_counts_mt` exceeds
  `QC_MITO_SUBCLUSTER = 15.0`.
- `QC: doublet-like` when at least `QC_DOUBLET_COEXPR_FRAC = 0.50` of cells
  co-express own and foreign lineage markers **and** the mean `doublet_score`
  is at least `QC_DOUBLET_SCORE_MULT = 1.5` times the lineage median. When
  `doublet_score` is unusable (which it currently is) co-expression alone
  carries the call, because Scrublet is the thing that failed, not the biology.

`qc_exclude` is a boolean and `qc_reason` records which rule fired. **Nothing
is deleted.** Downstream biology filters on `qc_exclude`; a QC-excluded cell
carrying a viral barcode is still the only evidence we have for that barcode.

### CITE-seq

CLR is per-cell `log1p(x)` minus the per-cell mean of `log1p(x)`, written to
`obsm['protein_clr']` and to `obs['ADT_<tag>']` for plotting. Tags are resolved
by exact match, then prefix-before-underscore, then prefix, because the feature
reference suffixes them inconsistently.

CD163 and CD169 are **structurally dead** and `adt_detection.csv` labels them
so. They are BD AbSeq 36-mers, incompatible with 10x 5' Antibody Capture: they
need a 3' poly-A and a BD PCR handle, so those molecules are never sequenced
into the ADT library. Settled, do not retest. The script warns if either shows
non-zero counts, since that would contradict the finding.

### Figures

300 DPI, PDF and PNG, font 28–34, all colors hex.

```
F05_umap_prebbknn_by_library      the batch effect, uncorrected
F05_umap_by_library               after BBKNN
F05_umap_by_animal / _by_timepoint
F05_umap_clusters                 Leiden 1.0
F05_umap_tier1                    corrected Tier 1
F05_umap_tier2                    Tier 2 subtypes
F05_umap_cd3_gate                 the hard gate
F05_umap_qc_exclude               flagged, not deleted
F05_umap_has_tcr
F05_umap_doublet_score
F05_umap_shiv_umi                 where the virus sits
F05_umap_ADT_<tag>                CD8 CD4 CX3CR1 CD28 CD69
F05_tier2_<lin>_subclusters       per lineage
```

These are diagnostic panels. Publication figures come later with the
`gen_mpl_labels` helper.

### The viral arm: what the audit cost, and what survives

**Per-cell viral attribution is not publishable.** The measurements:

- 59 cells carried ≥1 viral UMI, 179 UMI in total, median **1** molecule.
  38 of 59 carried exactly one.
- Three sat in lineages that cannot be productively infected: two
  Erythrocytes (no nucleus, no CD4) and one Mast cell. That is 5.1%
  known-false and a direct measurement that ambient viral RNA is present.
- **The strongest B cell (5 UMI, 0.53% of transcriptome) was weaker than the
  strongest erythrocyte (7 UMI, 0.58%).** That single comparison answers the
  B-cell question without any statistics.
- Only 2 of 59 cleared 1% of their transcriptome; none cleared 10%.
- No lineage beat its own abundance after correcting for the ten tested
  (B cells nominal p = 0.013, Bonferroni 0.132; Erythrocytes nominal
  p = 0.012, which is the same signal in an impossible category).

**The control arm was clean.** Zero viral cells in 26,026 uninfected Indian
control cells, and zero in library C. That rules out index hopping and a
mapping artifact as the explanation, which is genuinely good news about the
reference build.

**The background ceiling.** Above the 7 UMI maximum seen in the impossible
lineages sit four cells: three T cells and a monocyte, both legitimate
CCR5-tropic SHIV targets, and no B cells. The clearest is
`CATGCACAGAGCATTT-1-E` at 55 UMI and 5.2% of transcriptome, animal 41861,
21 DPI PBMC, with a productive TCR. That is a twentyfold jump over the next
cell and looks like a genuinely productively infected T cell.

The ceiling is derived from the data each run (`SHIV_CONFIDENT_MIN_UMI = None`)
and it is **fragile**: it rested on three cells, so one more background hit at
12 UMI would move it past two of the four survivors. 05 prints how many cells
it was derived from. Pin the value to freeze it.

**What survives.** The ambient viral RNA still came out of genuinely infected
cells somewhere in that animal, so **which barcodes are present in an animal,
and how many distinct ones**, remains a real measurement. What the ambient
destroys is the link from a molecule to a particular cell, not the existence
of the barcodes. So 06 is written on the barcode inventory, and the methods
should say that cell-level attribution was not attempted below the background
ceiling rather than disclaiming the arm wholesale.

Two things 05 now does so no figure can get this wrong by accident:
`shiv_call` is a three-level categorical (`none` / `background` / `confident`),
and `shiv_by_lineage.csv` prints detected and confident side by side with a
note that the detected column tracks abundance rather than tropism.

**Do not use `shiv_detected` for anything cell-level. Use `shiv_confident`.**

### Tier 2 assignment, rewritten 2026-11-03

Assignment is now from the **subcluster mean** of each signature, with
overlapping panels grouped into families and the gate on the margin to the
runner-up. `dominant_frac` is still computed and printed but **gates nothing**.

```python
MIN_SCORE_MARGIN     = 0.15     # top family mean minus runner-up
REQUIRE_POSITIVE_TOP = True     # a negative top score cannot name a subtype
REPORT_DOMINANT_FRAC = True     # homogeneity diagnostic, not a gate
PANEL_FAMILIES = {'T cells': {'Naive/CM': ['Naive', 'CentralMem'],
                              'Effector/Cytotoxic': ['EffectorMem', 'Cytotoxic']}}
```

#### Why the old statistic failed, which was not the reason first given

The first diagnosis blamed overlapping panels splitting the vote. That is part
of it but reconstructing the T population in simulation showed overlap alone
does not produce a 0.38 vote share. **The cause was the z-scoring.**

Dividing each panel by its SD across all cells is meant to make panels
comparable, but that SD is large for a panel with real signal and tiny for a
panel whose genes are barely expressed. Between-subcluster range of each T
panel:

```
Proliferating 1.701   Cytotoxic 1.393   EffectorMem 1.314   Naive 0.745
CentralMem 0.459  |  Treg 0.184   Th17 0.142   Tfh 0.107   Th1 0.090
```

*Th1*, *Th17*, *Tfh* and *Treg* never vary by more than 0.18 across the entire
T population, so they carry no information about which subcluster is which.
Z-scoring rescaled their pure noise up to the same unit variance as the real
signal and handed them equal voting power. They took 26% of the votes in a
subcluster whose Cytotoxic mean was 1.042.

`score_genes` already subtracts a control-gene mean, so a non-expressing panel
sits near zero by construction and the **raw means are directly comparable**
with no rescaling. A flat panel simply cannot out-rank a real one. It also
makes the log self-consistent: the "mean signature score per subcluster" table
is now the table the assignment is made from, which was not true before.

#### A panel-range floor was considered and rejected

Filtering out panels below a range threshold changes **not one label** once
assignment is on the mean, verified on all three lineages. It would have been
a second arbitrary threshold manufacturing `Unresolved` cells for no gain. The
range is printed as a diagnostic and gates nothing.

#### Families cannot be derived from gene overlap alone

| pair | shared genes | grouped | why |
|---|---|---|---|
| Naive + CentralMem | *CCR7 SELL IL7R* (3 of 5) | yes | nearly the same panel |
| EffectorMem + Cytotoxic | **none** | yes | cytotoxic program is co-expressed: subclusters 7 and 8 scored 0.80/0.80 and 1.04/1.01 while every other subcluster was negative on both |
| NonClassical + Intermediate | *FCGR3* | no | a panel defect to fix, not a family |
| any B pair | none | no | |

So the script prints the between-subcluster correlation of panel means as the
evidence, and families are declared in config with the reasoning in a comment.

#### `dominant_frac` is a homogeneity measure, not a confidence measure

A low value with a **large** margin means one population whose cells split
their vote between co-expressed panels. A low value with a **small** margin
means a genuinely mixed subcluster. Verified both ways in simulation: a 50/50
mix of naive and cytotoxic cells averages to nothing in the mean and shows up
only in the vote share, while a pure cytotoxic population with a near-tied
runner-up splits the vote while the mean stays unambiguous. Both statistics
are printed and the log names the subclusters on the bubble.

#### Predicted effect on the current labels

| lineage | Unresolved before | after | what changed |
|---|---|---|---|
| T cells | 9 of 11 | **3 of 11** | subclusters 7 and 8 become `Effector/Cytotoxic` (0.90 and 0.95 surface CD8, PBMC-enriched); 0, 1, 3, 5 become `Naive/CM` |
| Monocytes | 3 of 10 | **2 of 10** | subcluster 8, the 1,790-cell T population mislabelled `Monocytes DC`, correctly demoted to `Unresolved` at margin 0.023 |
| B cells | 5 of 10 | **6 of 10** | it refuses two calls it should not have made, including `B Naive` on a **negative** mean of -0.025 |

The ADT is the independent check and was not used in the assignment:
`Effector/Cytotoxic` lands on the two most CD8 subclusters at 0.90 and 0.95
surface CD8, `Naive/CM` on the CD4-leaning ones at 0.37 to 0.50.

Two honest caveats. B subclusters 0 and 1 keep `Memory` at margins of 0.156
and 0.151, which sit barely above `MIN_SCORE_MARGIN = 0.15`; they are on the
bubble by construction and the margin column plus `dominant_frac` is how you
see that. And **Monocyte subcluster 9 stays wrong**: a 582-cell B-cell
population labelled `Intermediate` at a margin of 1.202, because that is a
panel defect and no statistic fixes it.

### The panel defect: HLA-DRA does not exist in Mmul10

The 05 log said `Intermediate 2/3 MISSING ['HLA-DRA']` and nobody chased it.
The rhesus MHC class II alpha gene is **MAMU-DRA**, which appears throughout
the top-marker tables alongside *MAMU-DRB1*, *MAMU-DRB5*, *MAMU-DPA*,
*MAMU-DQA1* and *MAMU-DQB1*.

So the Monocyte Intermediate panel ran on *CD74* plus *FCGR3*, and *CD74* is
one of the most abundant transcripts in B cells. That is why "Intermediate"
scored 1.716 in a subcluster whose top markers are *CD79B*, *MS4A1*, *BANK1*
and *CD79A*. The panel had become an MHC-II score wearing a monocyte-subset
label. **A panel reduced to a non-specific gene produces confident wrong
answers, which is worse than `Unresolved`.**

05 now resolves panel genes through `GENE_ALIASES`, which maps the MHC genes
to their MAMU- equivalents and reports in the log which alias was used. Eight
genes remain missing and `scripts/05c_panel_audit.py` searches `var_names` for
candidates rather than guessing:

| gene | panel | priority |
|---|---|---|
| *HLA-DRA* | Monocyte Intermediate | fixed via alias |
| *CD27* | B Memory | high |
| *CD24* | B Transitional | high |
| *IGHD* | B Naive | high |
| *LILRB2* | Monocyte NonClassical | medium |
| *IL17A* | T Th17 | low, Th17 barely varies anyway |
| *PF4* | Platelets | low |
| *HBA2* | Erythrocytes | low |

*CD27*, *CD24* and *IGHD* are the three discriminating markers for naive,
memory and transitional B cells. Losing all three is why the B lineage cannot
separate, and no change to the statistic will fix that.

05c proposes and applies nothing. A wrong substitution is the failure it
exists to catch, and deciding *MAMU-DRA* stands in for *HLA-DRA* needs someone
who knows the biology rather than a regex.

### B panels replaced with RIRA modules, 2026-11-03

05c confirmed three alias substitutions and **rejected four**, because every
pattern candidate was a different gene:

| missing | candidates | verdict |
|---|---|---|
| *HLA-DRA* | *MAMU-DRA* | **accepted**, confirmed present |
| *PF4* | *PF4V1* | **accepted**, genuine platelet-expressed paralog |
| *HBA2* | *HBA* | **accepted**, second alpha-globin locus |
| *CD24* | *CD247*, *CD244*, *CD248*, *CD2* | rejected: CD3-zeta, 2B4, endosialin, CD2 |
| *CD27* | *CD276*, *CD274*, *CD2* | rejected: B7-H3, PD-L1, CD2 |
| *LILRB2* | *LILRB1*, *LILRB3* | rejected: LILR family is expanded and divergent in macaques, so 1:1 orthology cannot be assumed, and *LILRB1* is broader so it would make NonClassical **less** specific |
| *IGHD* | none | rejected; 25.2% of this matrix is LOC identifiers |
| *IL17A* | none | rejected, same reason |

So *CD27*, *CD24* and *IGHD* are **not recoverable by renaming**, and they
were the three discriminating markers for memory, transitional and naive B
cells. No change to the statistic fixes that.

**RIRA hit the same wall and solved it empirically.** From Mahyari et al. 2025
(Cell Genomics 5:100849, Figures 3E–3K): "while there were few
lineage-defining markers within mature B cells, we reasoned that scoring B
cells for enrichment of multi-gene modules would provide a classification
tool." The B panels are now their rhesus-validated modules:

| panel | genes | source |
|---|---|---|
| Naive | *TCL1A FCER2 CXCR4* | RIRA; *CXCR4* is their rhesus flow marker |
| Memory | *PLAC8 CR2 LY86 CD44* + *TNFRSF13B* | RIRA module; TACI kept from the old panel |
| **InnateLike** | *AHNAK CD1C SOX5* | RIRA; **new panel** |
| Transitional | *VPREB1 VPREB3 SOX4 MME* | RIRA pre-B markers |
| Plasma | *JCHAIN XBP1 IGHA1* + *MZB1 PRDM1 SDC1* | RIRA plus the old set |
| GerminalCenter | *RGS13 MEF2B AICDA BCL6* | RIRA, "consistent with published features" |

**InnateLike is the one to watch, and it was already visible in the data.**
RIRA describes CD40-low, PBMC-enriched innate-like B cells marked by *AHNAK*,
*CD1C* and *SOX5*. B subclusters 0 and 1 from the 2026-11-02 run carry exactly
that: *SOX5*, *AHNAK* and *CD1C* among their top markers, and both are
PBMC-enriched at +1.29 and +1.31 with LN at -2.00. Three independent lines
agree, so those two subclusters, which were getting a weak `Memory` call at a
margin of 0.15, are very likely innate-like B cells.

### Diagnostic panels: scored, printed, never assigned from

Same pattern and the same reasoning as `dx_cyto_shared` in Tier 1. A score
that is informative to look at is not therefore a label.

RIRA found the granzyme pattern to be the single most informative
discriminator among effector T and NK cells, above the canonical subset
markers: "the pattern of granzyme K/M versus granzyme A/B/H expression by
itself is informative at differentiating between many effector T and NK
subsets." From their sorted reference cells:

| granzyme profile | population |
|---|---|
| *GZMM* alone | naive CD8 |
| *GZMM* + *GZMK* | TCM / TTM |
| *GZMB* + *GZMH* | TEM |
| *GZMA* + *GZMB* + *GZMH* | CD16-positive NK |
| *GZMM* + *GZMK* | CD16-negative NK |

They also found S100A genes to be "a much more reliable signature of T cell
memory than genes encoding canonical protein markers".

```python
TIER2_DIAGNOSTIC_PANELS = {
    'T cells': {'dx_GranzymeKM':  ['GZMK', 'GZMM'],
                'dx_GranzymeABH': ['GZMA', 'GZMB', 'GZMH'],
                'dx_S100Amemory': ['S100A4','S100A6','S100A10','S100A11'],
                'dx_Perforin':    ['PRF1', 'GNLY']},
    'B cells': {'dx_Proliferating': ['MKI67', 'TOP2A']},
}
```

These stay out of the assignment because **they cut across it**. Cytotoxicity
is what dominates the clustering, so `Effector/Cytotoxic` answers "is this
cytotoxic" while the granzyme skew answers "TCM-like or TEM-like" *within*
that. A flat argmax cannot express two levels, and forcing it to would
recreate exactly the shared-program collision that broke the NK panel in
Tier 1.

05 prints a `granzyme skew` column, ABH minus KM, per subcluster. Positive
means TEM-like or CD16-positive NK; negative means TCM/TTM-like or
CD16-negative NK. That is how `Effector/Cytotoxic` gets promoted to a memory
stage through `SUBTYPE_RENAME`.

### Marker support: the cost of a stricter bar is reported, not taken

`MIN_MARKERS_DETECTED` stays at 1. At that bar it moved 365 cells, all NK,
because every lineage had 0.0% of cells at zero markers. The thin lineages are
thin differently: Macrophages, DCs, Mast, Erythrocytes and Platelets are each
called on a **median of one** detected marker. Raising the bar to 2 would
catch that and would also move a large number of real cells to `Unassigned`,
so `REPORT_ALT_MARKER_MIN = 2` prints the per-label cost and applies nothing.

### Superseded: the original Tier 2 diagnosis (kept for the record)

Measured on the 2026-11-02 run. T cells: **10 of 11 subclusters `Unresolved`**.
The cause is the statistic, not the biology.

Assignment takes the per-cell argmax over z-scored signatures, the mode per
subcluster, and requires `dominant_frac >= MIN_DOMINANT_FRAC = 0.40`. That
punishes panels that overlap. `Cytotoxic` and `EffectorMem` share the cytotoxic
program; `Naive` and `CentralMem` share *CCR7 SELL IL7R*. The per-cell argmax
splits between near-synonymous panels, so `dominant_frac` never clears the bar
even when the subcluster mean is unambiguous:

| T subcluster | Cytotoxic mean | EffectorMem mean | ADT CD8 | tissue | label given |
|---|---|---|---|---|---|
| 7 | **1.042** | 0.799 | 0.899 | PBMC | `Unresolved` (frac 0.381) |
| 8 | **1.014** | 0.803 | 0.946 | PBMC | `Unresolved` (frac 0.357) |

Those are CD8 cytotoxic effectors by every column and the rule refused to say
so. The B-cell calls are knife-edge in the other direction: Memory at 0.411,
0.411 and Naive at 0.410 all sit barely above 0.40, so they were labelled on a
hundredth of a fraction.

**Two subclusters are misassigned, not merely unresolved.** Inside the Monocyte
lineage:

| subcluster | n | label given | top markers | has_tcr |
|---|---|---|---|---|
| 8 | 1,790 | `Monocytes DC` | *SKAP1 CAMK4 BCL11B IL7R TNIK* | **52%** |
| 9 | 582 | `Monocytes Intermediate` | *CD79B MS4A1 BANK1 CD79A AFF3* | 3% |

Subcluster 8 is a T-cell population and subcluster 9 is a B-cell population.
Both passed marker support because abundant ambient myeloid transcripts (*LYZ*,
*S100A8/9*) clear a one-marker bar.

**Marker support as configured is nearly a no-op.** It moved 365 cells, all NK.
`pct_zero_markers` is 0.0 for every lineage except NK at 7.9%. The thin
lineages it was written to catch sit at a median of ONE detected marker
(Macrophages 1 of 5, DCs 1 of 5, Mast 1 of 4, Erythrocytes 1 of 2, Platelets
1 of 2) and `MIN_MARKERS_DETECTED = 1` does not touch them.

#### The proposed fix

1. Assign on the **subcluster mean** signature, not the per-cell argmax mode.
   The mean table already gives clean answers; the per-cell mode is the noisy
   step.
2. **Group overlapping panels into families** before the argmax, because the
   markers genuinely do not separate them: `{Naive, CentralMem}` →
   `Naive/CM`, `{EffectorMem, Cytotoxic}` → `Effector/Cytotoxic`, with Th1,
   Th17, Tfh, Treg and Proliferating standing alone. ADT corroborates the
   split: subclusters 7 and 8 are 0.90 and 0.95 ADT-CD8 while 3, 4 and 9 run
   0.21 to 0.37.
3. Gate on the **margin between the top two family means**, reporting
   `dominant_frac` as a diagnostic rather than using it as the gate.
4. Raise `MIN_MARKERS_DETECTED` to 2 so a lineage called on one ambient-prone
   gene fails, which is what sends the misassigned Monocyte subclusters to
   `Unassigned`.

Until this lands, `tier2_subtype` is a cluster identifier with a hint
attached, not a cell type.

### CD3 threshold: 3 was too low, now 5

The calibration that shipped with the fix immediately contradicted the value
chosen for it. Measured pre-gate on erythrocytes and platelets (1,495 cells):

| CD3 ≥ | background CD3+ | cells the gate moves |
|---|---|---|
| 1 | 35.79% | 9,297 |
| 2 | 20.60% | 3,609 |
| 3 | 12.11% | 1,805 |
| **5** | **3.28%** | **550** |
| 10 | 0.20% | 32 |

At 3 the gate was moving a population that is 12% ambient. 10 drives background
to nothing but leaves the gate inert. `CD3_MIN_COUNTS = 5` and
`CD3_BACKGROUND_TOLERANCE = 5.0`, the latter set to what is achievable on this
data rather than leaving a warning permanently lit. Raising the threshold
shifts work to the TCR evidence in CELL 9, which is the stronger claim anyway.

The "% of uncorrected NK calls that are CD3-positive" line also used to compare
a threshold-5 figure against the v10 benchmark of 62.2%, which was measured at
threshold 1. It now prints both and says which to compare.

### The animal-specificity result

`shiv_per_animal_within_library.csv`, which did not exist before this run.
Ambient is per GEM well, so under well-wide ambient the viral hits should
spread across the animals in a well in proportion to their cells. They do not:

| library | animal | observed | expected | fold | p |
|---|---|---|---|---|---|
| D | 40707 | 5 | 1.92 | 2.61 | 0.028 |
| E | 40707 | 8 | 1.99 | 4.01 | 0.0001 |
| F | 41862 | 4 | 0.84 | 4.74 | 0.002 |
| H | 40702 | 32 | 6.83 | 4.69 | 2e-5 |

Four libraries, one dominant animal each, and **a different animal in
different wells**. In library H, 32 of 34 viral cells are animal 40702 against
7 expected, at 1.5 UMI per cell, which is background-level burden concentrated
in one animal.

Well-wide ambient does not produce that. A plausible mechanism that does:
free viral RNA in each animal's own tube before pooling, adsorbing to that
animal's cells. Either way, the per-animal signal is a real correlate even
though per-cell attribution is not, and it is the strongest argument for
writing 06 around the per-animal barcode inventory.

Caveats that matter. The p-values are uncorrected and the tests are not
independent. A genuinely high-burden animal and a demultiplexing artifact are
indistinguishable in this table, and ~5% of cells are `Unknown` by demux.
40707 is the elite controller and shows the highest fold in two libraries,
which is counterintuitive, though both are 21 DPI when control may not yet be
established. And the four confident cells sit in four different animals,
including 41861 in library E where the detected cells are depleted.

### Two things still open

**The animal 40702 concentration.** 05b found library H holding 34 of 59
detected cells with roughly 32 of them in animal 40702. Nothing reported
per-animal cell counts *within* a library, so that was uninterpretable. 05 now
writes `shiv_per_animal_within_library.csv` with observed against expected
from each animal's share of its own well. Ambient is per GEM well, so under
background the hits should spread in proportion to cells. If 40702 beats its
share, that is where to look next, with the caveat that a genuinely
high-burden animal and a demultiplexing artifact are indistinguishable in that
table.

**Library C is 57.5% monocytes** against ≤15% elsewhere, with 19.8% T cells,
and `Monocytes Classical` was 4,232 of 4,557 in the Non-Infected arm. C is 83%
Indian control cells, whose mito median is 7.58%. That reads as degraded
samples with poor lymphocyte recovery. It matters more than the viral question:
the control arm is compositionally unlike the infected arm, so any
infected-versus-control contrast is confounded by composition before any
biology enters. Unresolved.

**One recoverable viral barcode.** `CGTCTTAGTCACAGCA-1-H`, 6 UMI in 603 total
(1.0% of transcriptome), removed by the gene filter rather than by QC, animal
`Unknown`, 27.4% mito. High fraction, dying cell, uncertain provenance.
Recoverable from `adata_prefilter.h5ad` if 06 wants it.

**Ignore *CD4* and *CCR5* as evidence in this dataset.** *CD4* mRNA is detected
in only 7.5% of T cells (3,685 of 49,228), so zero *CD4* in 27 B cells is what
dropout predicts (p ≈ 0.12). *CCR5* is detected in one cell in the entire
object. Neither column carries information here.

### Read these after 05

1. The **CD3 background calibration** printed in CELL 8. The rate in
   Erythrocytes and Platelets must be near zero at the chosen threshold, or
   the gate is moving ambient *CD3* into the T compartment.
2. `tier1_marker_support.csv`. B cells should pass cleanly; NK and the small
   myeloid lineages are where cells get reassigned to `Unassigned`.
3. `shiv_by_lineage.csv`, **the confident column only**.
4. `shiv_per_animal_within_library.csv`, for the 40702 question.
5. `F05_umap_prebbknn_by_library` against `F05_umap_by_library`.
6. `tier1_transition_matrix.csv`. The NK-to-T move is the correction working;
   v10 measured 62.2% of NK calls as CD3-positive and the rebuild measured
   58.0%, which is close enough to say nothing moved.
7. `tier2_*_subtype_table.csv`. Any `Unresolved` subcluster needs a decision,
   and `dominant_frac` near 0.40 means a split subcluster.
8. `tier2_*_doublet_evidence.csv`, since Scrublet is unreliable here and
   co-expression is carrying that call.

Then fill `SUBTYPE_RENAME` and re-run with `RESUME_FROM_TIER1 = True` to apply
the curated names in minutes.

---

## 06 — `06_viral_barcode.py`

Recovers Keele-style viral barcodes from reads **already aligned** to the SHIV
contig, and ends with an explicit verdict on whether targeted amplicon
enrichment is needed.

| | |
|---|---|
| **in** | `pipeline_results/03_cellranger_multi/*/outs/**/*.bam` (globbed) |
| | `pipeline_results/05_embed_annotate/adata_annotated.h5ad` (optional) |
| | `Binhua_FINAL_Barcode_*.xlsx` (validation only) |
| **out** | `pipeline_results/06_viral_barcode/` |
| **cost** | seconds to minutes, 4 cores, 64 GB. Not a long stage. |

### Why this replaced the FASTQ scan

The legacy approach (`shiv_bc_core_extract.py`, `shiv_vpr_flank_scan.py`,
`shiv_vpr_barcode_demux.py`) scanned raw FASTQ at full depth across all eight
libraries and took hours. That was never an alignment problem. Reads landed in
the BAM fine and were then discarded **at counting**, because they overlapped
two or more per-ORF GTF features. Stage 01 collapsing the contig to one gene
feature fixed the counting; the reads were always in the BAM.

So the FASTQ scan was solving a problem that no longer exists. Gone with it:
`MAX_READS_PER_LIB`, `PROGRESS_EVERY`, the smoke-test mode, the read caps, the
hit cache, and every host-genome cleavage-site scan. 06 reads a few thousand
SAM records.

Measured on library H, the richest of the eight (34 of the 59 viral cells):

| BAM | records | pairs | within 61 bp of the cassette |
|---|---|---|---|
| `sample_alignments.bam` | 862 | 431 | 13 |
| `unassigned_alignments.bam` | 90 | 45 | 0 |

Both BAMs are read, because `unassigned_alignments.bam` holds reads whose cell
barcode failed assignment, and a barcode recovered from one still proves the
barcode exists in that animal.

### R1 is in the BAM, and R1 and R2 are one molecule

431 R1 and 431 R2, exactly paired. 03 runs **SC5P-PE-v3** and Cell Ranger
aligns both mates in paired-end mode. The legacy scripts assert R1 is never
aligned; that was true of their STARsolo run (`CB_UMI_Simple` reads nothing
past base 28) and never of Cell Ranger multi in PE mode. Independent support:
the 00 GTF diagnostic measured `starsolo_r1` at **2.76x** fold recovery, so R1
reads did align to SHIV and were being lost at the GTF step.

R1 is a second look at the same molecule from the other end of the fragment,
which genuinely raises the per-molecule chance of seeing the cassette. It does
**not** double the molecule count. An earlier estimate that it would was wrong
by about 2x. Mates collapse on query name before anything is counted.

### Cassette geometry, verified not assumed

`pipeline_results/01_shiv_reference/cassette_context.txt` shows Brandon's
parental reference has **no cassette at all**: the flanks abut directly. The
34-N run is a landing pad stage 01 created, and Keele's published 44 bp and 42
bp flanks are verbatim the sequence either side of that junction.

```
cassette  [5862, 5896)   34bp, all N in the reference
  MluI L  [5862, 5868)    6bp  ACGCGC    cloning scar, NOT used
  reg L   [5868, 5874)    6bp  AGCGTA or TAGCTG
  CORE    [5874, 5884)   10bp  <- the variable barcode
  reg R   [5884, 5890)    6bp  CAGCTA or TACGCT
  MluI R  [5890, 5896)    6bp  GCGCGT    cloning scar, NOT used
```

**The MluI half-sites are out of the detection path on purpose.** They are
6-mers expected roughly every 4 kb by chance, which is why the old
`ACGCGC.{22}GCGCGT` regex returned three reads in total.

**The register markers stay in**, because they are what says which register the
core sits in: 707 of the 1,159 catalog entries are register 1 and 435 are
register 2. They are read off the alignment coordinate rather than searched for
genome-wide, so the specificity objection does not apply. They are masked in
the reference (everything is) but they are in the virus, present in 98.5% of
the catalog, and documented in Khanal et al. 2019.

### No canonicalisation between registers

`revcomp('AGCGTA') == 'TACGCT'`, so the two registers are the same cassette
read on opposite strands. But the catalog holds **1,158 distinct cores across
both**, so they are different barcodes rather than one barcode written two
ways, and the catalog stores each 34-mer as-written. 06 matches the extracted
window directly. Both orientations are tried and **which one hit is reported**,
so the first run tells us the convention instead of assuming it.

### The CIGAR walk, which is the crux

The cassette is N in the reference, so STAR prefers to end the alignment at the
N-run boundary and put the cassette bases into a **soft clip**. Soft clips
retain sequence in the BAM `SEQ` field; hard clips do not. `ref_to_query()`
projects clipped bases onto the reference positions they would have occupied
and tags every recovered base `aligned` or `projected`, so nothing downstream
mistakes a projection for an alignment. If most cores come out projected, the
log says so, because every one then rests on the no-indel-in-the-clip
assumption.

**SAM stores `SEQ` in reference orientation, including for reverse-strand
alignments**, so the CIGAR walk is identical for both strands and nothing is
reverse-complemented by FLAG. The design notes for this stage said to
reverse-complement minus-strand reads; that would have been a bug and it is
deliberately not done. The `reverse` column is logged and never used for
logic. What *does* read as a reverse complement is a cassette inserted in the
opposite register in the virus itself, and that is handled by trying both
markers.

`register_anchored_core()` finds the marker in the **raw read sequence**, with
the hit constrained to the query indices the projected window spans. Two
reasons it is not read out of the projected window: the window has a `.`
wherever a reference position has no query base, and stripping those dots
splices sequence across a deletion, so a core read from the stripped window can
be 10 bases that were never adjacent in the molecule. And a 6 bp marker found
anywhere in a 150 bp read would not be specific; the coordinate does the
localising and the marker only does the framing. This is the `canon_core` bug
from the legacy script, fixed: it trusted a fixed offset.

### Tested before delivery

The extraction functions were unit-tested against synthetic records covering a
trailing clip, a leading clip, a fully-spanned cassette, a minus-strand read, a
3 bp flank deletion, a 2 bp deletion inside the cassette, a decoy marker 90 bp
away, and a 6-base partial. All pass.

The whole stage was then run end to end against a synthetic BAM set with a
known answer: 7 planted cores across libraries D, E, F, H and one deliberate
decoy in negative library B, each core with a PCR duplicate sharing its
`(CB, UMI)`, plus one 6-base partial, plus background viral reads carrying no
cassette in all eight libraries. Every number came out right:

| check | expected | got |
|---|---|---|
| complete cores extracted | 21 (7 cores × 3 records) | 21 |
| partial held at tier 3 | 1 | 1 |
| pairs after mate collapse | 14 | 14 |
| molecules after `(CB, UMI)` dedup | 7 | 7 |
| decoy caught by the negative floor | 1 | 1 |
| in catalog | 5 of 7 | 5 of 7 |
| wrong-animal match flagged | 1 | 1 |

The wrong-animal case is the one worth noting: a core planted in library H
under animal 40702 that the catalog assigns to 40707 was flagged and printed,
which is the behaviour that matters, since agreement is weak evidence and
disagreement is strong.

### `MIN_CORE_BASES = 10`, not 6

The legacy scripts accepted partial cores down to 6 bases. Measured against the
catalog, the rate at which a random fragment matches a real core:

| bases visible | random match rate |
|---|---|
| 6 | **41.7%** |
| 7 | 12.9% |
| 8 | 3.4% |
| 9 | 0.88% |
| 10 | **0.11%** |

At six bases that is a coin flip, and it is how one barcode became three in the
2026-10-29 run (`ACCATT` / `ACCATTTCT` / `ACCATTTCTT`). Partials are reported
as tier 3 and never counted. The gate reads `MIN_CORE_BASES` directly and
warns loudly if it is ever set below the full core length.

### Confidence tiers

| tier | definition | counted |
|---|---|---|
| 1 | core ≥ `MIN_CORE_BASES`, register marker, **both** flanks ≥ 15 bp aligned | yes |
| 2 | same, **one** flank | yes |
| 3 | partial core or no register marker | **no**, reported only |

Tier 1 is trusted **even when the core is absent from the catalog**, for the
reason in the next section.

**A prediction that was wrong, recorded because it changes the methods.** The
first version of this section said tier 1 would be structurally empty, on the
reasoning that both flanks can only align if a read spans the cassette and STAR
would never align through 34 N. The 2026-11-04 run produced **3 tier-1
records**, so STAR does align through the N run. If tier 1 ever does come out
empty while tier 2 is populated, that means the aligner stopped at the boundary
on every read, and tier 2 is then the top tier available rather than a
degradation; the script prints that explanation when it happens.

### The catalog is validation, never a gate, and matching is exact

Brandon amplifies and mixes the whole stock before shipping, so the starting
pool is roughly 10,000 barcodes at unknown concentrations and nobody, including
Brandon, has the input list. The Binhua workbook is a **plasma readout at
specific draw dates**, 113 to 419 cores per animal across five animals
(39272, 40702, 40707, 41861, 41862). 34315 and 41903 are absent, which
independently confirms them as the uninfected controls.

Two consequences:

1. **A low catalog match rate is not evidence of false positives.** Most of
   what any animal carries was never in the plasma readout. 93% of catalog
   cores are animal-private, so the catalog is informative about *which* animal
   when it does hit, and uninformative when it misses.
2. **Error-tolerant matching is off entirely.** With 10,000 barcodes in a
   10-mer space each has about 0.29 neighbours one base away, so roughly a
   quarter of the real stock has a real `d=1` neighbour and a single sequencing
   error converts one real barcode into another. 6.2% of catalog cores already
   have a `d=1` neighbour inside the catalog alone. `EXACT_MATCH_ONLY = True`.

The animal cross-check (does the catalog animal agree with the demultiplexed
animal?) is a **sanity check, not a statistical test**. With 13
cassette-proximal records in the richest library there is no power. A core
matching the *wrong* animal is a genuine red flag; agreement is reassurance and
nothing more. An earlier claim that this was "the primary statistical test" was
wrong and is withdrawn.

### Molecules, not reads

Three collapses in order, because a repeat is only evidence of clonal expansion
if the observations are independent:

1. **mates to pair.** R1 and R2 of one query name are one observation.
2. **(CB, UMI) dedup.** Same cell, same UMI is one molecule amplified.
3. **truncation collapse.** A shorter core that is a prefix *or* suffix of a
   longer confirmed one merges into it. Suffix as well as prefix, because the
   two registers are reverse complements, so a truncated core can be compatible
   from either end depending on which register it came from.

### The negative floor

| class | why it is a floor |
|---|---|
| libraries A, B, C | pre-infection |
| library G | necropsy PBMC, measured zero viral reads in 05 |
| animals 34315, 41903 | uninfected Indian controls, absent from the catalog |

Any core called in a negative is a false positive **by construction**, so this
is the empirical error rate the tiers get judged against. Read
`negative_floor.csv` before believing any positive call.

### What this stage is actually for

The science question is how barcodes evolve over time, independent of the
starting mix. Viral reads exist only at **21 DPI** (D, E, F) and **Necropsy**
(H), so the only longitudinal statement available is whether a core appears at
both timepoints in the same animal. At this depth the answer is probably no,
and **that is a finding rather than a failure**: there is not enough coverage
to observe persistence even if it is occurring.

Per-cell viral attribution is already known not to be publishable (05b: the
strongest B cell carried less virus than the strongest erythrocyte). A barcode
recovered from a cell is evidence that the barcode exists in that **animal**,
not that the cell was infected. What survives is the per-animal inventory,
because the ambient RNA still came out of genuinely infected cells in that
animal.

**On enrichment.** The cassette is 34 bp of a 10,296 bp contig. Even with the
GTF fixed and both mates aligned, whole-transcriptome 5' capture puts a
vanishing fraction of reads across the barcode, and no analysis can create
coverage that was never sequenced. If 06 returns five or fewer distinct
barcodes, targeted VpxF1 amplicon enrichment on the same libraries is the only
route to a diversity or longitudinal result, and the log says so outright. The
limit is the library, not the method.

**Free observation, already in the 05 data.** At necropsy the virus is in LN
(library H, 34 viral cells) and absent from PBMC (library G, zero). Confounded
by depth, but worth recording.

### First real run, 2026-11-04, job 192594

**Two minutes of wall time reproduced the answer the 24-hour full-FASTQ scan
produced.** That is the result that matters for the handoff: the streamlined
method is not an approximation of the old one, it agrees with it.

| funnel stage | n |
|---|---|
| records aligned to the contig | 1,754 |
| pairs (molecules pre-UMI) | 877 |
| records within 60 bp of the cassette | 37 |
| records with a complete 10 bp core | 5 |
| records in tier 1 or 2 | 5 |
| molecules after `(CB, UMI)` dedup | 3 |
| molecules passing the negative floor | 3 |
| **distinct barcodes** | **3** |

Per library: A 0, B 0, C 2 (unassigned only, neither countable), D 344, E 352,
F 104, G 0, H 952. **Mates pair exactly in every library** (D 172/172, E
176/176, F 52/52, H 476/476), so R1 is in the BAM throughout and not just in H.

`negative_floor.csv`: **zero** molecules in any negative library or control
animal. There is no measurable false-positive floor at these tiers, which is
the good outcome and the thing to quote.

The catalog held 1,158 distinct cores across the five infected animals (113 to
420 entries each, 77–90% animal-private). Register composition came out 752
`AGCGTA` to 477 `TAGCTG`, with 10 entries in neither register, so the
register-marker assumption covers 99.2% of the catalog.

The 3 barcodes: `TCTGAAAGGT` (library E, tier 1, in catalog), `TATTCCCGGT`
(library H, tier 2, in catalog) and `TCGCCCTGTT` (library H, tier 1, not in
catalog).

### Two measured answers that settle the methods wording

Both were instrumented after job 192594 and answered by job 192595.

**1. The soft-clip projection recovered nothing.** All 5 countable cores came
out with all 10 bases **aligned** through the 34-N cassette; 0 were projected,
even in part. So STAR aligns through the N run and the CIGAR projection, which
the design treated as the crux of this stage, found no complete core. It is
insurance against a reference that could have forced clipping, not the
mechanism that produced the result. **Do not describe projection as
load-bearing in the methods.** It stays in the code because a different
reference build or aligner version could make it necessary, and it still
recovers the partial cores that populate tier 3.

**2. R1 earned its place, and the amount is measured.** `seen by R1 only: 1,
by R2 only: 0, by both: 2`. One of the three barcodes would have been **lost
outright** by an R2-only analysis, which is a third of the yield. That is the
justification for reading both mates, and it is a measurement rather than an
argument. The legacy scripts' assumption that R1 is never aligned would have
cost this barcode.

### Attribution, in three tiers of evidence

Job 192594 printed `unlinked` for 2 of the 3 barcodes, which conflated three
different situations that are not equally strong. `CELL 9` and `CELL 10` now
separate them and record which one applies in an `animal_source` column, so no
table or figure mixes them silently.

| `animal_source` | basis | strength |
|---|---|---|
| `cell` | the cell barcode is in `adata_annotated.h5ad` | strongest |
| `library` | the library contains exactly **one** animal in `obs`, so any molecule from it belongs to that animal regardless of whether its cell survived QC | as strong as the demultiplexing |
| `catalog` | the core is in the plasma catalog, which is 77–90% animal-private, so the catalog names one animal | **weakest** |
| `unattributed` | none of the above, or the core maps to several animals | none |

Single-animal libraries are **read from `obs` at run time**, not hardcoded, so
the script does not need to know the study design and cannot go stale against
it. A core matching several animals is left unattributed rather than guessed.

**`timepoint` and `tissue` are not in this hierarchy at all.** They are
properties of the *library*, so they are known whether or not a cell survived
QC, and printing them as `unlinked` threw away information the pipeline already
held. They now come from the library in every case. The only way either stays
`unknown` is a library that contributed molecules but has no cells in the
annotated object at all, which the script warns about by name, since that would
mean 04 or 05 failed to ingest it.

**A placeholder is not a value, and that cost a run.** Job 192748 still
printed `timepoint unknown` for the two library-H barcodes, even though H is
plainly Necropsy. The library fill requires exactly one distinct value, and
`dropna()` removes `NaN` but not the *string* `Unknown` that 05 writes for
cells which failed demultiplexing. So H's timepoints were
`{Necropsy, Unknown}`, length 2, and the fill declined. Tissue filled only
because H's tissue happens to be uniform. Placeholder strings are now stripped
before the uniqueness test, and the per-library animal count reports
placeholder cells separately rather than counting them as an extra animal.

**The circularity rule.** A catalog-attributed animal is *never* fed into the
animal cross-check. Using the catalog to assign the animal and then testing
whether the catalog agrees with the animal would manufacture a perfect score.
The cross-check runs only on `animal_source == 'cell'`.

### Two guards on the longitudinal claim

Persistence across timepoints is the headline longitudinal result, so the
weakest evidence tier must not be able to create it.

1. **A placeholder is not a timepoint.** A core seen at 21 DPI and again at an
   `unknown` timepoint has been seen at one known timepoint, not two. Counting
   the placeholder would invent persistence out of a missing value.
2. **Catalog attribution cannot supply the second timepoint.** The headline
   persistence line runs on `cell` and `library` evidence only. Candidates that
   appear only once catalog-attributed molecules are included are written to
   `barcode_persistence_catalog_only.csv` and labelled as not reportable: if
   the catalog supplied the animal, the catalog supplied the claim.

Both were verified against synthetic data in both directions, confirming that a
genuine persistence case on cell-plus-library evidence does still fire.

### The cell type call

`barcode_inventory.csv` now carries `tier1_celltype` and `tier2_subtype`,
available only for `cell`-attributed molecules since a molecule with no called
cell has no cell type by definition.

**This column is descriptive and is not evidence of infection.** 05b
established that per-cell viral attribution below `shiv_confident` is not
publishable, because the strongest B cell carried less virus than the strongest
erythrocyte. A barcode appearing in a T cell says the barcode is present in
that *animal*, and that this ambient molecule happened to land in a droplet
that also held a T cell. **Do not build a tropism claim on it.** The script
prints that caveat directly under the table so it travels with the output
rather than living only here.

### Bug found by job 192594 and fixed

The inventory printed **1 row while 3 barcodes existed**. `CELL 9` filled
`animal_id` with `'unlinked'` but left `timepoint` and `tissue` as `NaN`, and
`pandas.groupby` drops rows with a NaN key by default, so the two unattributed
barcodes were silently deleted from the headline table. They were in
`barcode_molecules.csv` the whole time. Fixed three ways: every linkage column
is filled, both groupbys pass `dropna=False`, and the inventory now asserts
that its molecule and core counts match what arrived and warns if they do not.

### Log hygiene fixed after the first run

Cell Ranger writes 40 BAMs per run and 24 of them are VDJ assembly BAMs with
**one contig per clonotype**. The near-miss warning printed those header lists
in full and buried the rest of the log under thousands of clonotype names.
VDJ BAMs are now skipped by path with a one-line count, and any remaining
near-miss list is truncated to four names plus a count. The submitter's
settings echo was also matching prose headings in the module docstring
("`MIN_CORE_BASES IS 10, NOT 6`") instead of the assignments, because its grep
pattern did not require an `=`.

### Outputs, in reading order

| file | what it answers |
|---|---|
| `detection_funnel.csv` | where the reads go and where they stop |
| `cassette_geometry.csv` | whether the core is readable at all |
| `negative_floor.csv` | the false-positive rate |
| `core_extractions.csv` | every record, every tier, every projection |
| `barcode_molecules.csv` | after mate collapse and UMI dedup |
| `barcode_inventory.csv` | per animal, timepoint, tissue, with evidence tier and cell type |
| `barcode_diversity.csv` | distinct cores per animal and timepoint, split by evidence tier |
| `barcode_persistence.csv` | the longitudinal question |
| `barcode_persistence_catalog_only.csv` | candidates that exist only on catalog evidence; **not reportable** |
| `animal_cross_check.csv` | catalog animal vs demultiplexed animal, cell-attributed only |
| `figures/F06_cassette_geometry` | where alignments stop, core completeness |

### Settings

```python
CONTIG            = 'SHIVAD8EO_V4fullN'   # checked against the BAM header
CASSETTE          = (5862, 5896)          # 0-based half-open
CORE              = (5874, 5884)
WINDOW_SLACK      = 8
REGISTER_MARKERS  = {'AGCGTA': 1, 'TAGCTG': 2}
USE_MLUI_AS_ANCHOR = False
MIN_CORE_BASES    = 10                    # was 6; see the table above
MIN_FLANK_ALIGNED = 15
EXACT_MATCH_ONLY  = True
PROJECT_SOFTCLIPS = True
INFECTED_LIBS     = {'D', 'E', 'F', 'H'}
NEGATIVE_LIBS     = {'A', 'B', 'C', 'G'}
CONTROL_ANIMALS   = ['34315', '41903']
EXCLUDE_FLAGS     = 0x900                 # secondary + supplementary
```

### Two implementation notes

**`samtools` via `subprocess`, not `pysam`.** `pysam` is not confirmed present
in `sc_pre` and there is no reason to add a build dependency to parse a few
thousand SAM lines. The submitter checks *which* `samtools` won the PATH race
and refuses to run on `/usr/local/bin/samtools`, which is a broken build on
titan.

**Paths are globbed and the contig name is read from the BAM header.** Nothing
is hardcoded. A guessed path already cost one round trip during development
(`per_sample_outs/<s>/count/sample_alignments.bam` does not exist; the real
path is `per_sample_outs/<s>/sample_alignments.bam`), and a silently renamed
contig is exactly how a believable zero gets produced. A near-miss in the
header is logged as a warning rather than passed over.

---

## 07 — `07_object_report.py`

A census of the finished object, written to be read by someone who did not
build it.

| | |
|---|---|
| **in** | `pipeline_results/05_embed_annotate/adata_annotated.h5ad` |
| | `pipeline_results/06_viral_barcode/barcode_inventory.csv` (optional) |
| **out** | `pipeline_results/07_object_report/07_object_report.md` + one CSV per table |
| **submit** | `sbatch run_07_object_report.sh` (4 cores, 128GB, ~2 minutes) |

It opens the object `backed='r'` and reads only `obs`, `var` and slot
metadata, so it never loads the matrix, cannot corrupt what it audits, and
does not need the memory to hold it.

### No figures, deliberately

**07 produces no plots and does not import matplotlib.** The figure set is a
separate decision that has not been made yet, and automating a figure nobody
has chosen produces work that gets thrown away. This report is the input to
that decision. The submitter states the same thing in its dependency check, so
a missing or broken matplotlib cannot stop the report.

### What it answers

1. **What is in here.** Every `obs` column with its dtype, cardinality, how
   much of it is populated, and **which stage wrote it**, inferred from name
   patterns. Columns matching no known pattern are listed rather than hidden,
   because an unrecognized column is either something the pipeline writes that
   this report does not know about or a leftover from an older run.
2. **How many cells**, by library, animal, timepoint, tissue, Tier 1 lineage
   and Tier 2 subtype, one factor at a time and then crossed.
3. **Which contrasts the design supports**, with the confound named where
   there is one.

### Placeholders are counted separately from missing values

05 writes `Unknown` for a cell that failed demultiplexing and `Unassigned` for
a lineage that passed no marker gate. Those are real cells with unresolved
labels, which is neither a `NaN` nor a usable value. `populated` therefore
counts them as **not** populated, because a column that is 100% non-null and
40% `Unknown` is not a usable column at face value. They are never dropped
from any table, and the report says how many there are for every factor.

### The three facts the design tables establish

**Library is not animal.** Every library is multiplexed across four or five
animals. Two consequences, and the second is the one that is easy to miss: a
per-library number is not a per-animal number, and ambient RNA pools across
the whole cell suspension before droplet formation, so the ambient transcripts
in a droplet come from every animal in the pool rather than from the animal
whose cell is in that droplet. That is why 06 cannot attribute an ambient
viral molecule to an animal by its cell barcode, and why the plasma catalog is
the *stronger* attribution route for an ambient molecule rather than the
weaker one.

**Composition differs sharply between libraries.** The report computes the
per-lineage share in each library and prints the spread. Library C is
monocyte-dominated, so a pooled per-cell differential test across libraries
will read that compositional difference as biology. Prefer a within-lineage
test or a per-animal pseudobulk, and report composition alongside every
cross-library contrast.

**Depth differs between libraries**, and depth drives detection. A lineage or
a transcript that looks absent from a shallow library may be under-sampled
there. This is the same confound behind the necropsy observation in 06: virus
present in LN and absent from PBMC is inseparable from the depth difference
between those two libraries.

### The answerable-questions section

For each candidate contrast the report states whether the cells exist on both
sides at or above `MIN_CELLS_PER_GROUP` (50) and `MIN_CELLS_PER_CELLTYPE`
(200), and names the confound. Nothing in it is a result; it is a map of what
a result could be built from. The two standing negatives are recorded there
permanently:

- **Cell-level viral tropism: not available.** 05b and 06 agree. Low-UMI viral
  signal is ambient, and the one 06 barcode attributable to a cell was in an
  erythrocyte, a lineage that cannot be infected. Not fixable by analysis.
- **Barcode diversity or lineage tracing: not available.** 3 barcodes across
  the whole experiment. Coverage, not method.

It also reports, per lineage, how many animals clear the threshold for
pseudobulk. **Per-animal pseudobulk is the analysis this design supports
best**: it treats the animal as the replicate, which is the unit the
experiment randomised, and it is immune to the compositional confound that
defeats a pooled per-cell test. Six lineages have all eight animals above
threshold; Mast cells have none.

### Two implementation notes on the report itself

**`subspecies` was missing from `FACTORS`** in the first version, which is how
the Chinese-versus-Indian confound went unreported until the second run. A
factor that is absent from `FACTORS` is invisible to every crossing in the
report, so adding one is not cosmetic.

**One pattern per concept.** Section 10 originally used a narrower TCR regex
than the column-provenance table, so it listed `has_tcr` and `chain_pairing`
but silently omitted `receptor_subtype`, `is_gamma_delta` and `tcr_recovered`.
Both now read the same pattern.

### Integrity checks

Three checks that look for the kind of defect nothing else in the pipeline
would catch, all three of which found something on the first real run:

**Byte-identical columns.** Candidates are grouped by an order-preserving hash
and then **confirmed by exact comparison**, and each group is labelled
`expected alias`, `REVIEW` or `BUG`, because two names for one quantity is
redundancy while two names for two different measurements is a defect.

The one real defect: `ADT_CD163_BDAbSeq` and `ADT_CD169_BDAbSeq` are identical
on every value, so one was written twice under two names and anybody reaching
for the second gets a silently wrong answer. Both are also entirely negative
after CLR (max −0.246), so they are duplicated *and* non-detected. **Check
`feature_reference.csv` and the ADT ingestion in 04 before using either as a
marker.** The other nine ADT columns are `_TotalSeqC` and look fine; only the
two `_BDAbSeq` ones are affected.

Several groups are benign aliases and are labelled as such rather than raised
as problems: `n_genes` with `n_genes_by_counts` (scanpy writes both),
`total_counts_shiv` with `shiv_umi_cellranger`, `shiv_whitelist` with
`shiv_detected` (identical by definition), `dx_NK_cells` with `dx_nk_old` (the
superseded NK panel, identical meaning the panel did not change), and
`tcr_recovered` with `marker_support_exempt` (05 assigns both from the same
boolean at lines 1524 and 1538).

### `outlier` does not mean what its name says in the annotated object

This came out of the duplicate check and is the most useful thing it found,
because it is a trap rather than a defect. `outlier == qc_would_have_dropped`
survived exact comparison, and the reason is **survivorship**: 04 computes
`outlier` from a 5-MAD test on three QC metrics, which flags a large fraction
of cells, and then drops `(outlier | mt_outlier) & ~shiv_whitelist`. Every
outlier that was not on the viral whitelist is therefore gone from the object,
so a surviving `outlier == True` cell must be whitelisted, which forces
`qc_would_have_dropped` true for exactly those cells.

So the handful of cells carrying `outlier = True` in `adata_annotated.h5ad`
are **not** "the QC outliers". They are the whitelisted cells QC would have
removed and that were deliberately rescued. The general outlier population
exists only in `adata_prefilter.h5ad`. Anybody treating this column as a QC
flag on the annotated object would be reasoning about a few rescued cells
while believing they were looking at all outliers. Use `qc_exclude` for
filtering; use `adata_prefilter.h5ad` if the question is about what QC
removed.

**A bug in this check, found and fixed.** The first version summed the
per-element hashes of each column. For a boolean column that sum collapses to
a function of the True count alone, so any two bool columns with the same
number of True values collided even when the True cells were entirely
different. The 2026-11-04 run therefore reported four spurious groups, one of
which wrongly grouped `shiv_confident` with `outlier` and
`qc_would_have_dropped`. All three simply had the same number of True cells.
The check now hashes the ordered sequence and then verifies every candidate
group by exact comparison, so the table cannot contain a false positive. If
`outlier` and `qc_would_have_dropped` survive that check they are worth a
look, since `qc_would_have_dropped` is `(outlier | mt_outlier) &
shiv_whitelist` and is therefore a strict subset of `outlier` by construction.

**Constant columns.** `doublet_whitelisted` is `False` for every cell, so the
doublet whitelist never matched anything. Worth knowing before relying on it.

**Per-animal composition outliers.** Section 6 reports the composition spread
per *library*, which cannot tell you whether a library is odd or one animal
inside it is, and that distinction decides whether the cause is batch or
biology. Animal **34315 is 52.2% monocytes while every other animal sits
between 4.6% and 9.7%**, a six-to-eleven-fold difference in one animal. It is
the whole reason library C reads 57.2% monocytes. A batch correction cannot
touch this, because it is not a batch effect.

### Three facts about the design that change what is worth testing

All three come out of the separability checks, and none were on the table
before 07 computed them.

**`subspecies` is perfectly confounded with infection status.** Every infected
animal is Chinese and both controls are Indian. An infected-versus-control
difference and a Chinese-versus-Indian difference are therefore *the same
comparison*, and no analysis of this dataset can tell them apart. Any such
result has to be reported with that stated rather than as an infection effect.
This is structural: adding cells cannot fix it.

**The controls share no timepoint with the infected animals.** Controls exist
only at `Non-Infected`; infected animals are at `Pre`, `21 DPI` and
`Necropsy`. So the contrast is not timepoint-matched and an infection effect
cannot be separated from a sampling-occasion effect.

**But libraries C and F each contain both arms.** Library F holds controls
34315 and 41903 alongside infected 41861 and 41862, all in LN. That supports a
**within-library infected-versus-control comparison that is free of the batch
confound** defeating the pooled version, and it is the cleanest contrast in the
dataset. Build on it before anything pooled. Sharing a library is not sharing a
timepoint, so check the timepoint composition inside it first.

### What the first run established about the object

| | |
|---|---|
| shape | 116,536 cells x 18,582 features, 1.63 GB |
| obs columns | 114, across 12 pipeline stages |
| features | 25.2% are `LOC` identifiers; 3,000 flagged highly variable |
| placeholders | 5,593 cells (4.8%) have no animal or timepoint; 11,990 (10.3%) no lineage |
| QC | 3,893 flagged `qc_exclude` (3.3%), 1,134 `predicted_doublet` (1.0%), none removed |
| depth | median `total_counts` ranges 865 (F) to 2,058 (A) across libraries |
| viral | 59 cells with any viral UMI, **4** `shiv_confident` |
| TCR | 44,303 cells (38.0%), of which 6,315 are off-lineage |
| ADT | 11 columns, 2 of them duplicates of each other |

**There is no clonotype ID column.** `chain_pairing` and `receptor_subtype`
describe each cell's chain configuration, which is not the same thing: neither
groups cells that share a receptor. So clonal expansion, clonal overlap
between timepoints and clonotype-to-phenotype links **cannot be computed from
this object as it stands**. 05 runs scirpy far enough to produce the chain
columns but never writes a clonotype definition. Recovering it needs a
`scirpy.tl.define_clonotypes` step added to 05 and a re-run, and it is worth
deciding on before any T-cell analysis is planned, because clonal structure is
usually the reason for adding VDJ in the first place.

### Settings

```python
FACTORS = ['library', 'animal_id', 'subspecies', 'timepoint', 'tissue',
           'tier1_celltype', 'tier2_subtype']
PLACEHOLDERS = {'unknown', 'unassigned', 'unlinked', 'none', 'nan', 'na',
                'n/a', '', 'nd', 'unresolved'}
CONTROL_ANIMALS = ['34315', '41903']
MIN_CELLS_PER_GROUP = 50
MIN_CELLS_PER_CELLTYPE = 200
```

### Two implementation notes

**Categorical columns.** `fillna` on a categorical refuses a value that is not
already a category, so placeholder detection goes through `astype(str)`
instead. This object is largely categorical and the naive version fails on the
first column.

**`.layers` carries a `None` key.** When an object is opened `backed='r'`,
`A.layers.keys()` includes `None` for `X` itself. It is not a layer, so it is
filtered out rather than reported as one.

---

## Diagnostics

Everything in `scripts/` with a number-plus-letter name is a read-only
diagnostic from this build. None of them write into a `pipeline_results/` path
a stage owns, and none are needed for a rerun; they exist because each one
settled a question that changed the pipeline.

| script | the question it settled |
|---|---|
| `00_gtf_multiplicity_diagnostic.py` | How many viral reads the per-ORF GTF was discarding. Measured 2.08x (cellranger), 2.76x (starsolo R1), 1.69x (starsolo R2) recovery from collapsing to one gene feature. This is why 01 emits a single-feature GTF. |
| `00b_refernce_comparison.py` | Whether the rebuilt reference matches an independent rebuild. md5 `106f756c6f62ac89d39942d18db4d96c`. |
| `03b_gate_review.py` | Whether the Cell Ranger cell-calling gate was sane per library. |
| `04b_qc_loss_diagnostic.py` | Where QC was losing cells. Found the control arm being decimated by a quantile mito filter. |
| `04c_mito_threshold_sweep.py` | What mito threshold to use. A fixed 20% took the control arm from 7,753 to 26,026 cells; the quantile version was removing real cells. |
| `05b_viral_and_label_audit.py` | Whether the viral signal was real. Found zero viral cells in 26,026 controls, 3 of 59 viral cells in lineages that cannot be infected, and the strongest B cell weaker than the strongest erythrocyte. This audit is why the viral arm is reported per animal and not per cell. |
| `05c_panel_audit.py` | Which marker genes actually exist in Mmul10. Found `HLA-DRA` absent, which had been silently running the Intermediate monocyte panel on `CD74 + FCGR3`; accepted three aliases and rejected four unrecoverable genes. |

`run_03b_gate_evaluation.sh` and `run_05b_viral_and_label_audit.sh` are the
submitters for the two that are slow enough to want one.

The remaining ~80 files in `scripts/` are legacy work kept for reference and
are not part of this pipeline. `scripts/README.md` covers those.

---

## Conventions

- `# %%` Spyder cells throughout. Cell 1 holds every tunable.
- Warn and continue over hard stops, except where continuing would produce a
  silently wrong answer (a missing viral feature, X that is not raw counts, no
  exclusive NK anchor).
- Figures: 300 DPI, PDF **and** PNG, all text 28–34 pt, all colors as hex.
- Gene names italicized when referring to the gene or mRNA, capitalized and
  upright for the protein.
- Diagnostics are read-only and live in `scripts/`.
- Every Python stage has a `run_<stage>.sh` submitter beside it. Each one
  activates `sc_pre`, puts `CONDA_PREFIX/bin` first on `PATH` (the system
  `samtools` at `/usr/local/bin` is a broken build that otherwise wins), runs
  `python -u` so the log fills as it goes, and echoes the settings that decide
  what the run actually did.
- QC-excluded cells are **flagged**, never deleted. A QC-excluded cell carrying
  a viral barcode is still the only evidence we have for that barcode.
- Cluster IDs never appear in a label map. Assignment is by signature with the
  margin reported, because Leiden IDs do not survive a rebuild.
- **Every threshold on a noisy measurement is calibrated against cells that
  cannot have the thing being measured.** Erythrocytes and Platelets cannot
  express *CD3*, so their CD3-positive rate sets the ambient floor for the gate.
  They also cannot be infected, so their viral burden sets the background
  ceiling for a viral call. Those lineages are kept and used as controls, never
  filtered out for being implausible.
- A derived threshold prints how many cells it was derived from, so its
  fragility is visible rather than implied.
- Checkpoints carry a recipe stamp of every setting that decides a label.
  Resuming under a changed recipe is refused with a diff, never applied
  silently.
