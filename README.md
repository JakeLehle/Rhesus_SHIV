# Rhesus SHIV

Single-cell multi-omics analysis of barcoded SHIV-AD8EO infection in rhesus macaques.
5' scRNA-seq + CITE-seq + TCR VDJ from Cell Ranger `multi`, analysed in Scanpy on the
Texas Biomed HPC.

Scientific lead and data owner: Zoey Wallis. Analysis: Jake Lehle.
Virus and barcode ground truth: Brandon Keele's group (barcoded SHIVAD8EOM);
per-animal plasma barcode catalog from Binhua.

---

## 1. Study design

Eight physical 10x libraries (A through H) across three timepoints and two tissues.
Libraries are **sequencing pools, not conditions**. C, D and F are mixed-condition
pools where animals from different timepoints or subspecies share a GEM well, so
per-animal condition assignment happens downstream in Scanpy using the provider's
genetic demultiplexing barcode assignments.

| Lib | Sample | Composition |
|-----|--------|-------------|
| A | `SHIV_Pre_PBMC` | Chinese: 39272, 40702, 40707, 41861 |
| B | `SHIV_Pre_LN` | Chinese: 39272, 40702, 40707, 41861 |
| C | `SHIV_Pre_PBMC_IndianMixed` | **mixed** Indian: 34315, 41903 + Chinese 41862 |
| D | `SHIV_21DPI_Pre_LN_Mixed` | **mixed** 21 DPI: 39272, 40702, 40707 + Pre 41862 |
| E | `SHIV_21DPI_PBMC` | Chinese: 39272, 40702, 40707, 41861 |
| F | `SHIV_21DPI_NonInf_LN_Mixed` | **mixed** 21 DPI: 41861, 41862 + non-inf 34315, 41903 |
| G | `SHIV_Necropsy_PBMC` | Chinese: 40702, 40707, 41861, 41862 |
| H | `SHIV_Necropsy_LN` | Chinese: 40702, 40707, 41861, 41862 (reservoir library) |

Five infected Chinese animals, two Indian non-infected controls. **40707 is the elite
controller** and the centrepiece of most contrasts. Two gaps that break requested
comparisons: 39272 has no necropsy at all, and 41862 has no 21 DPI PBMC because its
21 DPI sample sits in library F, which is lymph node. Six of Zoey's requested
contrasts are therefore impossible and are enumerated in the header of
`SHIV_DEG_comparison.py`.

Chemistry is **GEM-X 5' v3**, confirmed empirically rather than from product naming
(CB 16, UMI 12, whitelist `3M-5pgex-jan-2023.txt`). See `check_sample_whitelist.py`.

---

## 2. References: where they come from and where they go

This is the part worth reading before rerunning anything.

### Host + virus GEX reference

Built by **`MKREF_Mmul10_SHIVAD8EO.sh`**.

- Host: Mmul10, `GCF_003339765.1_Mmul_10_genomic.fna` + full unfiltered GTF,
  pulled from the NCBI FTP site by the script.
- Virus: **GenBank `MN816822.1`**, the SHIV-AD8EO proviral sequence, pulled by
  `efetch` and renamed to `>SHIVAD8EO`. **10,262 bp.**
- Staging directory: `/master/jlehle/WORKING/SC/REF/`
  (`SHIVAD8EO.fna`, `Mmul10_SHIVAD8EO_v10.fa`, `Mmul10_SHIVAD8EO_v10.gtf`)
- Final Cell Ranger reference: `/master/jlehle/cellranger-10.0.0/Mmul10_SHIVAD8EO_v10/`
  (some scripts reach it as `GEX_REF=/master/jlehle/WORKING/SC/REF/Mmul10_SHIVAD8EO_v10`
  with `GEX_REF_ALT` pointing at the cellranger copy)

Lessons carried forward from the Mtb/SIV project and encoded in the script: use the
full unfiltered GTF, every transcript needs at least one exon entry for Cell Ranger
v10, spaces in the GTF source column break `mkref`, verify FASTA contig names against
GTF seqnames before building, and use `--expect-cells` rather than `--force-cells`.

### A standing note in MKREF that turned out to matter

`MKREF_Mmul10_SHIVAD8EO.sh` records the decision that the 34 bp barcode does **not**
need to be in the reference, on the reasoning that Cell Ranger maps to coding regions
and barcode detection happens post-alignment. That is why the reference in use is the
unbarcoded parental. The consequence is traced in section 6: reads crossing the
cassette have nothing to align to there.

### Reference relationship to Keele's SHIVAD8EOM file

Keele's barcoded reference (`SHIVAD8EOM_reference.txt`, 9,563 bp) and
`MN816822.1` (10,262 bp) are the same virus in two different frames. Anchoring on the
two Khanal primers in both files pins the relationship exactly:

| Feature | Keele SHIVAD8EOM | MN816822.1 | offset |
|---|---|---|---|
| VpxF1 forward | 5192 | 5790 | +598 |
| 34 bp barcode cassette | 5264..5298 | **absent** | |
| vpr RT primer (revcomp) | 5533 | 6097 | +564 |

The upstream offset is +598 and the downstream offset is +564. The difference is
exactly 34 bp, which is the cassette present in Keele's file and missing from ours.
Total length also reconciles: 9,563 + 598 upstream + 135 downstream − 34 cassette
= 10,262. So the two references differ by LTR trimming at both ends plus the cassette,
with no detectable indel between the two primer anchors. Keele's file carries the
cassette as 34 literal `N`s, which cannot be aligned through, so it is a coordinate
and design document rather than a drop-in replacement reference.

### Other references and databases

| Reference | Built by | Location |
|---|---|---|
| Rhesus VDJ (IMGT, all four TCR loci) | `build_rhesus_vdj_reference.sh` | `${WORKDIR}/vdj_reference/vdj_rhesus_macaque_IMGT` |
| CITE-seq feature reference | hand-curated | `${WORKDIR}/feature_reference.csv` |
| SHIV-only STARsolo index + collapsed GTF | `shiv_starsolo_recovery.sh` | `${WORKDIR}/shiv_starsolo/ref/SHIVAD8EO.fa`, `SHIVAD8EO_collapsed.gtf` |
| Lentivirus BLAST nucleotide DB | `shiv_build_databases.sh` | `${WORKDIR}/viral_id_databases/blast_nucl_db/lentivirus_nucl` |
| Lentivirus DIAMOND protein DB | `shiv_build_databases.sh` | `${WORKDIR}/viral_id_databases/diamond_prot_db/lentivirus_prot.dmnd` |
| geNomad marker DB v1.9 | `shiv_build_databases.sh` | `${WORKDIR}/viral_id_databases/genomad_db` |
| Standalone MN816822.1 FASTA + GenBank | `Diagnose_SHIV_Missing_Reads.sh`, `shiv_candidate_contigs_vs_ref.sh` | `${WORKDIR}/shiv_divergence_analysis/MN816822.1.fasta` |

`WORKDIR` is `/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV` throughout.
Raw FASTQs are Zoey's, read-only, at `/master/zwallis/WORKING/SC/10X_RAW`.

When searching the lentivirus databases, grep on `"simian immunodeficiency"` and
`"HIV-1"`. A literal `SIVmac` returns almost nothing because records title themselves
"Simian immunodeficiency virus (Mac251)" and similar.

---

## 3. Object chain

```
Cell Ranger multi per-library outs
   └─ scanpy_analysis_script.py
         adata_qc.h5ad              animal_id, tissue, timepoint, subspecies, QC, raw ADT
         adata_pp.h5ad              + embedding, clusters, CLR ADT
         adata_annotated.h5ad       + Tier 1 / Tier 2 labels
   └─ shiv_scanpy_merge.py          attaches recovered SHIV signal onto an existing object
         shiv_host_merged.h5ad
   └─ Rhesus_T_NK_Cell_Correction.py
         shiv_host_tier1_corrected.h5ad
   └─ shiv_VDJ_integration.py
         shiv_host_tier1_tcr.h5ad
   └─ shiv_embeding.py
         shiv_host_final_embedded.h5ad      <- canonical object for figures and DEG
   └─ shiv_tier2_annotation_write.py
         shiv_host_tier2_annotated.h5ad
   └─ SHIV_tier2_myeloid_annotation_write.py
         shiv_host_tier2_myeloid_annotated.h5ad
```

`adata_annotated.h5ad` carries `shiv_positive` / `total_counts_shiv` from Cell Ranger,
which are the **undercount**. The recovered `shiv_umi` / `shiv_pos` attached by
`shiv_scanpy_merge.py` are the ones to use.

---

## 4. Scripts by effort

### 4.1 Reference construction

| Script | What it does |
|---|---|
| `MKREF_Mmul10_SHIVAD8EO.sh` | Builds the Mmul10 + SHIVAD8-EO Cell Ranger reference |
| `build_rhesus_vdj_reference.sh` | Builds the IMGT rhesus VDJ reference and validates Zoey's custom TCR enrichment primers against it |

### 4.2 FASTQ integrity

Libraries B and H arrived with transfer corruption. 10x FASTQs are generated in
lockstep, so read N in R1 is read N in R2/I1/I2 and damage is always at the tail.
Truncating every file to the shortest is therefore the correct paired-end repair.

| Script | What it does |
|---|---|
| `diagnose_truncated_fastqs.sh` | Counts reads across all files per library, reports mismatches |
| `repair_corrupt_fastqs.sh` | Salvages valid gzip stream, trims to last complete 4-line record, truncates all files to the minimum, recompresses, verifies. Backs up originals to `REPAIR_BACKUP/` |

Outcome: Library B repaired to 82.3%, Library H to 18.1% for VDJ.

### 4.3 Alignment and quantification

| Script | What it does |
|---|---|
| `CELLRANGER_MULTI_RH_SHIV.sh` | Cell Ranger `multi` across all three modalities per library (GEX, PROTEIN, VDJ_RH) with the corrected library-to-sample mapping |
| `check_sample_whitelist.py` | Confirms CB length, UMI length and whitelist empirically from the BAMs rather than from the chemistry name, resolving GEM-X vs v3 |
| `check_all_libraries_shiv.py` | Scans all 8 h5 matrices for non-zero SHIV counts; also catches host gene name collisions (`TAT` vs `tat`) |

### 4.4 The SHIV detection problem

The arc that consumed most of the project. SHIV detection came in at 46 cells out of
126K, and this group of scripts worked out why.

| Script | What it does |
|---|---|
| `diagnose_shiv_counts.py` | Works backwards from the count matrix to find where 46 SHIV+ cells came from when `idxstats` showed 0 reads on the SHIVAD8EO contig |
| `shiv_mapping_diagnostic.py` | Cross-library view of SHIV read availability and loss. Decomposes the SAM FLAG, classifies reads by LTR / single-gene / multi-gene-overlap position. Established that loss is **annotation** (overlapping ORFs plus unannotated LTRs), not STAR stringency |
| `Diagnose_SHIV_Missing_Reads.sh` | Extracts SHIV-mapped and unmapped reads, permissive minimap2 remap, SPAdes assembly, BLAST against MN816822.1 to test divergence |
| `shiv_targeted_remap.sh` | Permissive minimap2 remap of the unmapped pool for libraries D, E, H. Confirmed the reference is **not** divergent: 864 mapped reads cover 97.4% of the genome at `nM:i:0` |
| `shiv_starsolo_recovery.sh` | Barcode-aware recovery. Original R1+R2 through STARsolo against a SHIV-only index with relaxed length filters, collapsed single-feature GTF, `soloStrand REVERSE`. Lifted Library H from 27 to 31 SHIV+ cells |
| `shiv_merge_barcode_diagnostic.py` | Diagnosed the NovaSeq X barcode orientation mismatch: Cell Ranger reverse-complements the CB, plain-whitelist STARsolo does not, so the join key never matched |
| `shiv_scanpy_merge.py` | Metadata-preserving merge. Attaches recovered SHIV signal onto an existing pipeline object instead of re-reading raw cells, which is what the first version did wrong |

Conclusion: low SHIV detection is a **mapping and annotation recovery problem, not
viral divergence**, with a genuinely low viral transcription floor underneath it.

### 4.5 De novo assembly and viral identification

Reference-free route, in case a divergent quasispecies was escaping recruitment.

| Script | What it does |
|---|---|
| `shiv_denovo_assembly.sh` | Phases A to E: characterise the 864 mapped reads, extract unmapped, MEGAHIT assemble, BLAST against SHIVAD8EO, pull SHIV-like contigs |
| `shiv_full_denovo_assembly.sh` | Assembly only, all ~280M unmapped Library H reads, `--min-count 1`, `--presets meta-sensitive`, no pre-filter and no subsampling. 118,889 contigs, N50 356 bp, max 4,447 bp |
| `quick_inventory.py` | Read-only orientation pass over the STARsolo / unmapped / MEGAHIT read pools |
| `shiv_build_databases.sh` | Builds all three viral ID databases from scratch with a provenance manifest. No searching happens here |
| `LENTI_DB_DIAG.sh` | Probes lineage composition and redundancy depth of the built databases |
| `SHIV_BLASTN.sh` | blastn of all contigs, no `max_target_seqs` cap, e-value as the only selection |
| `SHIV_DIAMOND.sh` | DIAMOND blastx `--very-sensitive` for remote homology, no target cap |
| `contig_label_diagnostic.sh` | Parses the trailing `stitle` field out of the blastn and DIAMOND tables |
| `contig_full_diagnostic.py` | Full per-contig diagnostic over the search output |
| `troubleshoot_contig_diagnostic.py` | Profiling probe for where the diagnostic runtime was going |
| `shiv_candidate_contigs_vs_ref.sh` | Extracts the 706 candidate contigs and blasts them against a single-reference DB built from MN816822.1 |
| `shiv_summarize_candidate_contigs_vs_ref.py` | Builds the gene map from the GenBank annotation and summarises candidate hits by region |

Result: 14 contigs mapped to SHIVAD8-EO at 95 to 100% identity, which is what ruled
divergence out.

### 4.6 Viral barcode recovery (vpr cassette)

The attempt to link Keele barcodes to single cells. **This line is currently closed**;
see section 5.

Cassette design (Khanal et al. 2020, J Virol, confirmed against the Binhua workbook):
34 bp = 12 bp flank + 10 bp core + 12 bp flank, outer 6 bp being the MluI ligation
scar `ACGCGC ... GCGCGT`. That anchor is its own reverse complement at the outer
level, so a forward search catches either strand. Inner flank `AGCGTA` vs `TAGCTG`
says whether to reverse-complement the core to canonicalise it.

| Script | What it does |
|---|---|
| `shiv_vpr_barcode_feasibility.py`<br>`run_shiv_vpr_barcode_feasibility.sh` | Feasibility only. Locks the anchor geometry against the ground truth, locates the cassette site and the Khanal primers in the reference, counts anchor hits in BAMs and raw R2 across all 8 libraries |
| `shiv_vpr_read_diagnostic.py`<br>`SHIV_vpr_read_diagnostic.py` | Asks one question: do reads reach the insertion point and soft-clip the cassette the reference lacks? Decides whether reference surgery is worth building |
| `shiv_vpr_barcode_demux.py`<br>`run_shiv_vpr_barcode_demux.sh` | Treats the barcode as a feature-barcode problem rather than an alignment problem. Ceiling measurement over the complete raw read pool with negatives A/B/C/G as the false-positive floor |
| `shiv_vpr_barcode_linkage.py` | Extracts cores from vpr-window reads two ways and links barcode-positive cells to animal / tissue / timepoint / cell type |
| `shiv_vpr_barcode_verify.py` | Verification pass over the recovered set |
| `SHIV_vpr_barcode_gt_search.py`<br>`run_shiv_vpr_barcode_gt_search.sh` | Five matching tiers including partial-overlap and mismatch-tolerant search against the closed 1,245-member ground-truth set. Also fixes a silent parser bug that dropped 21 of 1,264 workbook entries |

Ground truth: `Binhua_FINAL_Barcode_4726.xlsx`, later renamed
`Binhua_FINAL_Barcode_4-7-26_Rhesus_SHIV_Barcodes.xlsx`. Five sheets, one per infected
Chinese animal, 1,264 entries in the `bc_seq` column.

Outcome: 3 barcoded reads across roughly 964M reads per library set. The permissive
tiers added single digits, which is the convergence that justified recommending
targeted VpxF1 amplicon enrichment instead.

### 4.7 CITE-seq / BD AbSeq recovery

The CD163 and CD169 tags returned zero counts. Three non-exclusive hypotheses were
tested: a pattern bug in `feature_reference.csv`, a PCR handle mismatch, or GEX
carryover of a poly(A)-tailed BD oligo.

| Script | What it does |
|---|---|
| `bd_abseq_probes.py` | Step 1. Builds the 16-probe set and the FASTQ manifest |
| `bd_abseq_scan.sh` | Step 2. Read-length report, full-file scan of every R2 for the BD probes, plus a subsampled TotalSeq-C control scan that makes a zero interpretable |
| `bd_abseq_summarize.py` | Step 3. Parses every hit line, locates each probe in each read, filters on the evidence side, writes a verdict |

Outcome: BD AbSeq 36-mer barcodes need a 3' poly-A plus a BD PCR handle and are
**completely incompatible** with 10x 5' Antibody Capture construction. Those molecules
are never sequenced into the ADT library. Myeloid Tier 2 annotation therefore proceeds
RNA-driven. Nine of eleven ADT tags are usable; the anti-CD8a tag cannot separate
CD8aa NK from CD3-dropout T cells, so it is treated as a surface map only.

### 4.8 Host processing pipeline

| Script | What it does |
|---|---|
| `scanpy_analysis_script.py` | The main pipeline, sections 0 to 12. Corrected library mapping, provider barcode-to-animal lookup, gene category tagging, MAD-based QC plus Scrublet, normalise / embed / BBKNN, Tier 1 and Tier 2 annotation, CITE-seq validation, SHIV+ identification, gene set scoring, export |
| `shiv_comprehensive_diagnostic.py`<br>`run_shiv_comprehensive_diagnostic.sh` | Read-only object-level diagnostic over four horizon items: the 109 paired-receptor cells not labelled T, clonotype feasibility for the 40707 CD8 story, plotting inputs, and low-purity clusters |
| `shiv_embeding.py` | PCA / neighbours / UMAP / Leiden plus ADT CLR on the VDJ-integrated object, reproducing the pipeline's section 5 recipe so the embedding stays consistent |
| `shiv_final_object_check.py` | Integrity check with consistency cross-checks: every gamma-delta cell is a T cell, no NK cell still carries a TCR, SHIV+ counts preserved, `animal_id` preserved |

### 4.9 Tier 1 annotation correction

| Script | What it does |
|---|---|
| `Rhesus_T_NK_Cell_Diagnostic.py` | Read-only. Quantifies the NK-vs-T misclassification before any logic changes. Found 62.2% of Tier-1 NK calls are CD3+ with a decisive median margin of 0.762, and 39.1% of Unassigned cells are CD3+ T cells |
| `Rhesus_T_NK_Cell_Correction.py` | Rebuilds the NK score on exclusive anchors only (*KLRD1*, *KLRC3*, *NCAM1*), applies a CD3 hard gate, recovers CD3+ Unassigned to T, flags CD3-negative retained NK as provisional. Computes both old and new labels and prints the full transition matrix |
| `shiv_tier1_tcr_recovery.py` | Recovers CD3-dropout T cells carrying a called TCR, with B and myeloid quarantined because those lineages carry documented ambient orphan-TCR rates of 9 to 14% |

The CD3 hard gate before cytotoxic scoring is now a standing convention on this
project. Flat winner-take-all scoring leaks lymphocytes across lineages.

### 4.10 VDJ / TCR

| Script | What it does |
|---|---|
| `shiv_VDJ_inventory.py` | Read-only inventory of what `vdj_t` actually produced per library, chain coverage across TRA/TRB/TRG/TRD, and the barcode join yield, before any merge was written |
| `shiv_VDJ_integration.py` | scirpy 0.22.5 MuData workflow. Attaches receptor status, resolves provisional NK, calls gamma-delta strictly from paired TRG+TRD |
| `shiv_tcr_abdg_breakout.py` | Alpha-beta vs gamma-delta per animal and timepoint, testing whether 40707 holds a steady gamma-delta fraction while the others spike then fall |

Gamma-delta is called from the paired TRG+TRD subtype only, which is a floor rather
than a count, because TRD dropout is real and TRG persists in alpha-beta T cells.

### 4.11 Tier 2 annotation

| Script | What it does |
|---|---|
| `shiv_tier2_diagnostic.py`<br>`shiv_tier2_diagnostic_v2.py`<br>`run_shiv_tier2_diagnostic.sh` | Per-lineage recluster with BBKNN, functional signature scoring, subcluster profiles and animal / timepoint / tissue enrichment. Prep for annotation, not annotation |
| `shiv_tier2_cluster9_adjudication.py` | Settles whether the CD79/MHC-II-high subcluster is residual B-doublet contamination or a genuine MHC-II-high activated T subtype |
| `shiv_tier2_annotation_write.py` | Writes final T subtype labels with every cluster's basis recorded. CD4/CD8 surface call, proliferation gate and CD4 activation score are carried as per-cell attributes rather than label names |
| `SHIV_tier2_myeloid_diagnostic.py` | v2 of the myeloid diagnostic. v1 was not annotatable because 62% of the subset was not myeloid biology and the panels were collinear (mean pairwise signature r = 0.41, PC1 explaining 49%). Adds a hard lineage-purity gate before subclustering plus QC signatures so artifact clusters self-identify |
| `SHIV_tier2_myeloid_annotation_write.py` | Writes myeloid labels onto the purity-gated subset and structurally refuses to label any cell whose Tier 1 call is not myeloid |
| `shiv_tier2_deg.py` | Cell-level DEG within fixed subtypes and within infected Chinese animals. Wilcoxon for ranking only; the evidence layer is the per-animal consistency panel |
| `Shiv_tier2_animal_enrichment.py` | Tier 2 animal enrichment heatmaps for both compartments on named subtypes, plus three UMAP embeddings |
| `Shiv_annotation_reference_collector.py` | Collects every marker panel and every cluster-to-label mapping into CSVs with provenance, parsing the source scripts with `ast` so the tables cannot drift from the code that produced the labels |

### 4.12 DEG and pathway enrichment

| Script | What it does |
|---|---|
| `SHIV_DEG_comparison.py` | Zoey's full DEG manifest across tissues, timepoints, animals, subspecies and Tier-2 populations, plus cross-tissue contrasts and a batch-controlled subspecies variant. Documents the six requested contrasts that cannot exist in the data |
| `SHIV_ORA_list_backfill.py` | Rebuilds the ORA gene lists from DEG tables already on disk, recomputing nothing, instead of rerunning hours of pseudobulk |
| `SHIV_ORA_verify_totals.py` | Verifies the backfill against the original run using the exact identity up + down == total at any `abs(log2FC) >= 1.5` gate |
| `SHIV_enrichment_prep.py` | Step 1 of 2. Connectivity probe, gene ID mapping and KEGG `mcc` gene set construction. Establishes that Enrichr has no macaque, so the route is KEGG `mcc` plus g:Profiler `mmulatta`, with no human ortholog conversion anywhere |
| `SHIV_enrichment_run.py` | Step 2 of 2. KEGG and GO over-representation per contrast, up and down tested separately, through two independent engines (g:Profiler server-side and a local hypergeometric validated against a worked example) |
| `SHIV_KEGG_local_diagnostic.py` | Explains why the local KEGG engine returned zero on all 35 lists while g:Profiler returned 51 terms, by counting through five candidate causes rather than assuming a bug |
| `SHIV_enrichment_figures.py` | KEGG bubble, GO term network and volcano. Python equivalents of the three R figures, recomputing no statistics. Uses fold enrichment rather than GeneRatio on x so up and down panels are comparable |

Pathway analysis here is over-representation on a gene list, not ranked-list GSEA, and
up and down genes are run as separate enrichments to keep direction.

### 4.13 Figures

| Script | What it does |
|---|---|
| `shiv_plot_figures.py`<br>`run_shiv_plot_figures.sh` | Cell types and Leiden clusters on the UMAP, the 9 usable surface markers, metadata overlays, composition by tissue and timepoint |
| `shiv_plot_figures_updated_T_cell_corrected.py` | Same figure set regenerated after the Tier 1 T/NK correction |

---

## 5. Current state

**Done and validated**

- Cell Ranger `multi` complete for all 8 libraries with the corrected mapping
- 126,252 cells loaded, 89,844 after QC
- Tier 1 annotation complete including recovery of CD3-dropout T cells
- Tier 2 T cell annotation complete; Tier 2 myeloid annotation complete, RNA-driven
- VDJ integrated, gamma-delta called strictly from paired TRG+TRD
- DEG comparison suite and KEGG/GO over-representation complete
- Figures delivered for Zoey's committee meeting

**Closed lines, recorded so they are not relitigated**

- *Viral divergence as the cause of low SHIV detection.* Ruled out: 97.4% genome
  breadth at zero mismatches, and 14 assembled contigs at 95 to 100% identity.
- *BD AbSeq CD163/CD169 recovery.* Ruled out as a reagent incompatibility with 10x 5'
  chemistry. Not recoverable in silico.
- *Single-cell linkage of the vpr viral barcode.* Closed 2026-09-25 as unsalvageable by
  in-silico recovery. The honest path is targeted VpxF1 amplicon enrichment, which
  needs new sequencing.

**Open**

- The vpr barcode line is being revisited after advice from Keele to anchor on the
  virus sequence flanking the cassette rather than on the cassette's own invariant
  half-sites. The flanks verify exactly against his reference and abut the 34 bp
  cassette with no off-by-one. Whether this changes the outcome depends on whether any
  reads reach the vpr window at all, which is the question
  `shiv_vpr_read_diagnostic.py` was built to answer.

---

## 6. Conventions

- **Diagnostic-first.** Read-only diagnostics before any object modification. Labels
  are only written in dedicated annotation passes after the tables have been reviewed.
- **Evidence funnel, not pre-filtering.** Run permissive, assign evidence tiers
  post-hoc, filter on the evidence side where the filtering is visible.
- **Modular and sequential.** Each script does one thing and writes to new files
  rather than overwriting. Every tunable parameter lives in Cell 1 at the top.
- **Warn and continue** rather than hard stops.
- **Dropped findings stay documented** in script headers with the reasoning.
- Spyder `# %%` cell format throughout; `matplotlib.use('Agg')` under SLURM.
- Figures saved as both PDF and PNG at 300 DPI, all text 28 to 34 pt, hex colours
  defined once and reused.
- Negative libraries are always carried at full depth as the false-positive floor.
  A zero is only interpretable next to a control that fires with the same machinery.

### Environment

- HPC: Titan / Zeus, SLURM `--partition=normal`, logs at `/master/jlehle/WORKING/LOGS/`
- Conda envs: `sc_pre` (primary, Python 3.11, scanpy / scirpy 0.22.5 / BBKNN / LIANA),
  `star` (STAR 2.7.11b for STARsolo), `viral_id`, `viral_id_genomad`
- Cell Ranger v10.0.0 for GEX (v9.0.1 for the VDJ reference build)

Two standing hazards, both learned the hard way:

- The system `samtools` at `/usr/local/bin/samtools` is ancient and broken. Always
  `source ~/anaconda3/bin/activate` then
  `export PATH="${CONDA_PREFIX}/bin:${PATH}"` in SLURM jobs.
- Cell Ranger's bundled STAR is too old for STARsolo flags. Use the separate `star`
  env.

### Rhesus naming gotchas

Mitochondrial genes use `KEG06_` locus tags, not `MT-`. MHC genes use `MAMU-`, not
`HLA-`. SIV/SHIV genes carry a `gene-` prefix on `var_names` with plain symbols
underneath. Marker panels are from RIRA (Bimber et al. 2025, Cell Genomics), which is
rhesus-validated.# Rhesus_SHIV
