#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
05 — Embedding, Tier 1, TCR recovery, Tier 2
=============================================

Consolidates five legacy scripts into one stage:

    shiv_embeding.py                  embedding
    Rhesus_T_NK_Cell_Correction.py    T/NK correction
    shiv_VDJ_integration.py           TCR attach, gamma-delta, NK resolution
    shiv_tier1_tcr_recovery.py        CD3-dropout T recovery
    shiv_tier2_diagnostic_v2.py       per-lineage subclustering and signatures
    shiv_tier2_annotation_write.py    subtype labels and per-cell attributes

    in   pipeline_results/04_ingest_qc/adata_qc.h5ad
    out  pipeline_results/05_embed_annotate/adata_annotated.h5ad

WHY THESE ARE ONE SCRIPT AND NOT FOUR
--------------------------------------
Ordering, not convenience. TCR recovery adds CD3-dropout T cells back into the
T pool, and the T pool has to be COMPLETE before subclustering, or the T
subclusters are built on a truncated population and every subtype fraction is
wrong. In the old chain each script wrote an intermediate object whose only
purpose was to be reloaded by the next one, and each re-normalized defensively
because none could assume what the others had done. Normalization happens once
here.

The cost of merging is that iterating on Tier 2 would mean re-running the
embedding. RESUME_FROM_TIER1 solves that: the object is checkpointed after
Tier 1 and TCR recovery, and a resumed run skips straight to Tier 2. Same
pattern as REUSE_PREFILTER in 04.

EMBEDDING PARAMETERS ARE CARRIED FORWARD UNCHANGED from the legacy scripts.
HVG seurat 3,000 global and 2,000 within-lineage, PCA 50, BBKNN 30, Leiden 1.0,
seed 42, SCORE_THRESHOLD 0.1, SUB_RES 1.0, GD_SUBTYPES TRG+TRD.

FOUR THINGS CHANGED after the 05b audit, all of them measurements rather than
preferences. Each is argued at its config entry in Cell 1:

  1. CD3_MIN_COUNTS 1 -> 3. At one count, Erythrocytes were 20.2% CD3-positive
     and Platelets 13.0%. Neither can express CD3, so a one-count gate was
     measuring ambient mRNA. B cells sat at 19.1%, inside that background
     range, which is why the "CD3-positive B cells" finding dissolved.

  2. The CD3 gate now applies to EVERY non-T label, not just NK and
     Unassigned. The gate exists because the scores are untrustworthy, and
     which label a bad score produced is not a reason to trust it.

  3. A label must be backed by at least one of its own markers. Assignment was
     argmax over score_genes with a floor on the SCORE, not on detection.
     B cells pass easily (0.0% zero-markers); NK at 14.3%, and Macrophages,
     DCs and Mast cells at a median of ONE marker each, are the thin calls.

  4. The viral call is graded. >=1 UMI is not a detection: of 59 such cells,
     38 carried exactly one molecule and three sat in lineages that cannot be
     infected. shiv_confident requires burden above a ceiling measured from
     those impossible lineages, or 1% of the transcriptome.

WHAT THE AUDIT ESTABLISHED, AND WHAT IT COSTS
----------------------------------------------
The strongest viral signal in any B cell (5 UMI, 0.53% of transcriptome) was
WEAKER than the strongest signal in an erythrocyte (7 UMI, 0.58%). Above the
background ceiling sat three T cells and a monocyte, both legitimate
CCR5-tropic targets, and no B cells. The control arm was clean: zero viral
cells in 26,026 uninfected cells, which rules out index hopping.

So per-cell attribution of these molecules is not publishable, and that is a
real loss. What survives is the barcode inventory per animal, because the
ambient RNA still came out of genuinely infected cells somewhere in that
animal. The cell-level link is what the ambient destroys, not the existence of
the barcodes. 06 should be written on that basis.

BATCH_KEY IS 'library', NOT 'animal_id'
----------------------------------------
04b showed the dominant structure is per animal, which tempts you toward
animal as the batch key. That would delete the experiment. Animal is
confounded with the biology: Indian versus Chinese IS infection status, and
40707 being the elite controller is the centrepiece of most contrasts.
Correcting on animal removes exactly the signal the project exists to measure.
Library is the GEM well, the genuine technical batch, and libraries C and F
each carry both subspecies so library correction does not flatten subspecies.

A pre-correction UMAP is saved as X_umap_prebbknn so the batch effect can be
seen rather than assumed.

TIER 2 IS ASSIGNED BY SIGNATURE, NOT BY A HARDCODED CLUSTER MAP
----------------------------------------------------------------
shiv_tier2_annotation_write.py carried SUBTYPE_MAP = {'0': 'CD8 effector
memory', '5': ..., '10': ...}. Leiden cluster IDs are not stable across a
reference change, a cell-count change, or a seed change, so that map silently
mislabels everything after any rebuild, and this rebuild changed all three.

Instead the subtype is assigned from the SUBCLUSTER MEAN of each signature,
with overlapping panels grouped into families and the gate on the margin to
the runner-up.

The first attempt at this used the per-cell argmax of z-scored signatures,
took the mode per subcluster, and gated on that vote share. It left 10 of 11 T
subclusters Unresolved, including two that are CD8 cytotoxic effectors by
every other measure. The cause was the z-scoring: dividing each panel by its
SD across cells rescales an uninformative panel's noise up to the same unit
variance as a real signal, so four T panels whose means never varied by more
than 0.18 across the whole population took 26% of the votes in a subcluster
whose Cytotoxic mean was 1.042.

score_genes already subtracts a control-gene mean, so the raw means are
comparable with no rescaling and a flat panel simply cannot out-rank a real
one. The vote share is still computed and printed as dominant_frac, but it
GATES NOTHING: it measures homogeneity, not confidence, and it is how you find
the subclusters on the bubble.

Your curated names from the v10 adjudication are NOT thrown away and NOT
hardcoded either. SUBTYPE_RENAME in Cell 1 is empty by default. Read the
review tables this script writes, then promote provisional names to curated
ones by filling it in. The curation stays explicit, auditable, and attached to
evidence instead of to cluster integers.

DOUBLETS ARE CARRIED, NOT REMOVED
----------------------------------
Scrublet's automatic threshold is failing on this dataset: six of eight
libraries return 0.0%, library G returns 5.9%, 1,134 total against an expected
~8%. 04 now has REMOVE_DOUBLETS = False. Any threshold substituted there would
be invented rather than measured. Doublets are handled here as evidence: a
subcluster is flagged QC: doublet-like when it co-expresses markers of two
lineages AND carries an elevated mean doublet_score. Flagged, never deleted,
so 06 can still recover a viral barcode from one.

QC-EXCLUDED CELLS ARE FLAGGED, NEVER DELETED
---------------------------------------------
qc_exclude is a boolean. Biological analysis downstream filters on it. Nothing
is dropped, because a QC-excluded cell carrying a viral barcode is still the
only evidence we have for that barcode.

CHANGE LOG
==========
Every change since the first version, with the measurement that forced it.
Each entry is argued in full at its own config entry in CELL 1.

2026-11-03  Exclusions now apply to the WITHIN-LINEAGE HVG, not just global
    CELL 4 excluded shiv and mt; CELL 12 excluded nothing. So Tier 2
    subclustering was driven by whatever technical gene varied most. Of each
    subcluster's top 15 markers: T sub 0 13/15 ribosomal, T sub 4 12/15,
    T sub 2 11/15, B sub 9 11/15, B sub 10 10/15, B sub 5 9/15, and T sub 8
    and B sub 2 were 9/15 and 10/15 KEG06_* (the rhesus mitochondrial genome).
    39,852 of 83,004 T and B cells, 48%, in subclusters defined by technical
    content. T sub 8 and B sub 2 had median pct_counts_mt of 7.46% and 7.57%
    against the Indian controls' 7.58% and were 100% Indian controls: the
    control animals grouped by sample quality. One shared hvg_exclusion()
    helper now serves both cells so they cannot drift apart again.

2026-11-03  Ribosomal, translation and MHC class I added to the exclusions
    Following RIRA, which excludes the same sets for the same stated reasons:
    class I polymorphism corrupts counts in a single reference, ribosomal
    genes are batch-prone. MHC class II is deliberately KEPT: the polymorphism
    argument does not apply to it, MAMU-DRA is in the Monocyte panel, and
    MHC-II separates B cells and DCs. Verified that no panel gene is caught.

2026-11-03  Monocyte 'Intermediate' panel DELETED, DC split into lineages
    Intermediate was CD74 + FCGR3 + HLA-DRA. Resolving HLA-DRA to MAMU-DRA
    made it worse: all three are then MHC-II or Fc, nothing monocyte-specific.
    Range 2.300, winning 7 of 12 subclusters. It mislabelled 993 B cells
    (MS4A1 CD79B BANK1) at margin 0.845 and 142 cDC2 (FCER1A) at margin 2.356.
    RIRA found TWO monocyte populations in 41,199 rhesus myeloid cells, not
    three. NonClassical gained CX3CR1 and ITGAX, which RIRA found more
    specific in scRNA-seq than FCGR3 alone.

2026-11-03  Effector/Cytotoxic auto-splits on the measured granzyme skew
    The skew separated two populations on its first run: sub 5 at +0.878
    (GZMB/GZMH, PBMC-enriched, perforin-high) and sub 6 at -1.572 (GZMK/GZMM).
    Both were one label of 10,015 cells and SUBTYPE_RENAME could not tell them
    apart, being keyed on the provisional name. Splitting on a MEASURED
    quantity is stable across rebuilds; keying on subcluster ID is not.

2026-11-03  top_family_range reported
    "Treg, margin 0.005" read as if it meant something when the Treg panel's
    whole between-subcluster range is 0.100. The row was correctly Unresolved;
    the column makes the reason visible.

2026-11-03  B panels replaced with RIRA rhesus-validated modules
    05c confirmed HLA-DRA->MAMU-DRA, PF4->PF4V1, HBA2->HBA and REJECTED four
    where every pattern candidate was a different gene: CD24->CD247 is
    CD3-zeta, CD27->CD276 is B7-H3, LILRB2->LILRB1/3 in a family that is
    expanded in macaques, IGHD and IL17A no candidate at all. CD27, CD24 and
    IGHD were the discriminating markers for memory, transitional and naive B
    cells, so renaming could not fix the lineage. RIRA hit the same wall and
    solved it with modules. A new InnateLike panel (AHNAK CD1C SOX5) was added
    and immediately confirmed: B subclusters 0 and 1, 6,601 cells, margins
    0.406 and 0.320, dominant_frac 0.638 and 0.609, PBMC-enriched at +1.36 and
    +1.32. Four independent lines agreeing.

2026-11-03  RIRA granzyme and S100A modules added as DIAGNOSTIC panels
    Scored and printed, never assigned from, same pattern as dx_cyto_shared in
    Tier 1. They cut across the assignment rather than competing with it:
    Effector/Cytotoxic answers "is this cytotoxic", the granzyme skew answers
    "TCM-like or TEM-like" within that. A flat argmax cannot express two
    levels, and forcing it to would recreate the shared-program collision that
    broke the NK panel.

2026-11-03  CD3_MIN_COUNTS 3 -> 5
    The calibration that shipped with the threshold immediately contradicted
    it. Measured background CD3+ in erythrocytes and platelets, which cannot
    express CD3: 35.79% at >=1, 20.60% at 2, 12.11% at 3, 3.28% at 5, 0.20%
    at 10. At 3 the gate was laundering a population that is 12% ambient. 10
    leaves it inert at 32 cells. Tolerance set to 5.0, which is what the data
    can actually meet rather than a warning left permanently lit.

2026-11-03  The NK-versus-v10 line stopped comparing across thresholds
    It printed a threshold-5 figure against v10's 62.2%, which was measured at
    threshold 1, making a threshold change look like a biological shift. Now
    prints both and says which is comparable. At >=1 this object gives 58.0%
    against v10's 62.2%, so nothing moved in the reference change.

2026-11-02  Tier 2 assignment moved from per-cell vote to subcluster mean
    The per-cell argmax of z-scored signatures, moded per subcluster and gated
    on dominant_frac >= 0.40, left 10 of 11 T subclusters Unresolved including
    two that are CD8 cytotoxic effectors by every other measure. The cause was
    the z-scoring: dividing each panel by its SD across cells rescales an
    uninformative panel's noise up to unit variance, so four T panels whose
    means never varied by more than 0.18 took 26% of the votes in a subcluster
    whose Cytotoxic mean was 1.042. A panel-range floor was considered and
    REJECTED: it changes not one label once assignment is on the mean.

2026-11-02  CD3 gate generalized to every non-T label; marker support added
    The gate previously touched only NK and Unassigned, which was where the
    v10 collision happened to be known. A label must now also detect at least
    one of its own markers. Marker support runs BEFORE TCR recovery, because
    recovered cells carry no CD3 by construction; the ordering is the
    safeguard and marker_support_exempt records it.

2026-11-02  Graded viral call; checkpoint recipe stamp
    >=1 viral UMI is not a detection: of 59 such cells, 38 carried exactly one
    molecule and three sat in lineages that cannot be infected. The strongest
    B cell (5 UMI, 0.53%) was weaker than the strongest erythrocyte (7 UMI,
    0.58%). shiv_confident requires burden above a ceiling measured from those
    lineages, or 1% of the transcriptome. The checkpoint carries a stamp of
    every Tier 1 setting; resuming under a changed recipe is refused with a
    diff rather than silently reapplying the old gate.

Environment:  conda activate sc_pre
Submit with:  sbatch run_05_embed_annotate.sh
Author: Jake Lehle
Date:   November 2026
"""

# %% ============================================================
# CELL 1 — CONFIGURATION (every tunable lives here)
# ============================================================
import os
os.environ.setdefault('MPLBACKEND', 'Agg')

import glob
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse, stats

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import warnings
warnings.filterwarnings('ignore')

WORKDIR = '/master/jlehle/WORKING/SC/fastq/Rhesus_SHIV'
IN_H5AD = os.path.join(WORKDIR, 'pipeline_results', '04_ingest_qc',
                       'adata_qc.h5ad')
MULTI_DIR = os.path.join(WORKDIR, 'pipeline_results', '03_cellranger_multi')
OUT_DIR = os.path.join(WORKDIR, 'pipeline_results', '05_embed_annotate')
CKPT_H5AD = os.path.join(OUT_DIR, 'adata_tier1_checkpoint.h5ad')
OUT_H5AD = os.path.join(OUT_DIR, 'adata_annotated.h5ad')
FIG_DIR = os.path.join(OUT_DIR, 'figures')

# --- resume ----------------------------------------------------------------
# True skips the embedding, Tier 1 and TCR recovery and loads the checkpoint.
# Use it when iterating on Tier 2 only. The checkpoint is written at CELL 11.
RESUME_FROM_TIER1 = False

# --- global embedding. Unchanged from shiv_embeding.py. --------------------
COUNTS_LAYER = 'counts'
TARGET_SUM = 1e4
HVG_FLAVOR = 'seurat'          # dispersion on log-norm X
N_TOP_HVG = 3000
EXCLUDE_FROM_HVG = ['shiv', 'mt']      # var-column flags set in 04

# Pattern-based HVG exclusions, added 2026-11-03. Applied to the global HVG in
# CELL 4 **and** the within-lineage HVG in CELL 12, which is the bug this
# fixes: CELL 4 excluded shiv and mt, CELL 12 excluded nothing, so Tier 2
# subclustering was driven by whatever technical gene happened to vary most.
#
# Measured on the 2026-11-03 run: of each subcluster's top 15 markers,
#
#     T sub 0   13/15 ribosomal        6,400 cells
#     T sub 4   12/15 ribosomal        5,243 cells
#     T sub 2   11/15 ribosomal        3,820 cells
#     T sub 8    9/15 KEG06 (mito)     6,749 cells
#     B sub 9   11/15 ribosomal        4,803 cells
#     B sub 10  10/15 ribosomal        6,458 cells
#     B sub 5    9/15 ribosomal        2,732 cells
#     B sub 2   10/15 KEG06 (mito)     3,647 cells
#
# which is 39,852 of 83,004 T and B cells, 48%, sitting in subclusters defined
# by technical content. That is why the biological panels scored near zero
# there, and why Naive/CM collapsed from 0.290 to 0.061 between runs.
#
# The two mito-driven subclusters are the clearest case. T sub 8 has a median
# pct_counts_mt of 7.46% and B sub 2 has 7.57%, against the Indian control
# animals' measured median of 7.58% (04c). Both subclusters are 100% Indian
# controls. They are the control animals grouped by sample quality, not by
# biology, and KEG06_* is the rhesus mitochondrial genome: KEG06_p01 to p13
# are the 13 protein-coding genes and KEG06_r01/r02 the two rRNAs. 04 tags
# them correctly as mt, so pct_counts_mt and the 20% filter were always right;
# only the Tier 2 HVG step missed them.
#
# RIRA does exactly this, and for the same stated reasons: "The top variable
# gene list was filtered to exclude genes that produce inaccurate counts, such
# as the HLA class I region, where high polymorphism across the species cannot
# be represented appropriately in a single reference genome, as well as gene
# sets empirically identified to be prone to sample-to-sample batch effects,
# such as ribosomal genes."
#
# MHC CLASS II IS DELIBERATELY NOT EXCLUDED. RIRA's reason for dropping class
# I is polymorphism making the counts inaccurate, which does not apply to
# class II, and class II is genuinely informative here: MAMU-DRA is in the
# Monocyte panel and MHC-II separates B cells and DCs. Set
# EXCLUDE_MHC_CLASS_II = True to drop it as well, but understand that it would
# remove real signal.
#
# These exclusions apply to HVG SELECTION ONLY. Every gene stays in the object
# and is still scored, plotted and available to rank_genes_groups, which is
# why KEG06_* and the ribosomal genes will still appear in the top-marker
# tables. They just no longer define the axes.
EXCLUDE_HVG_PATTERNS = {
    'ribosomal':       r'^(RP[SL]|MRP[SL])',   # RPS*, RPL*, RPLP*, mitoribo
    'ribosomal_fusion': r'^(FAU|UBA52)$',      # RPS30 and RPL40 fusion genes
    'translation':     r'^EEF1',               # EEF1A1 etc, same batch-prone class
    'mhc_class_I':     r'^MAMU-(A|A3|AG|B|E|F)$',
}
EXCLUDE_MHC_CLASS_II = False
EXCLUDE_HVG_IN_LINEAGE = True          # the fix: apply all of this in CELL 12 too
N_PCS_PCA = 50
N_PCS_NEIGHBORS = 30
LEIDEN_RES = 1.0
RANDOM_STATE = 42
BATCH_KEY = 'library'          # the GEM well. NOT animal_id. See the header.
DO_PREBBKNN_UMAP = True

# --- Tier 1. RIRA panels (Bimber et al. 2025), rhesus-validated. -----------
SCORE_THRESHOLD = 0.1
TIER1_MARKERS = {
    'T cells':      ['CD3D', 'CD3E', 'CD3G'],
    'B cells':      ['CD79A', 'MS4A1', 'CD19', 'BANK1', 'PAX5'],
    'NK cells':     ['KLRD1', 'KLRC3', 'GNLY', 'PRF1', 'GZMB', 'NKG7'],
    'Monocytes':    ['S100A8', 'S100A9', 'VCAN', 'CD14', 'LYZ'],
    'Macrophages':  ['CD68', 'CD163', 'MRC1', 'MERTK', 'APOC1'],
    'DCs':          ['FLT3', 'CLEC9A', 'CLEC10A', 'IL3RA', 'CLEC4C'],
    'Mast cells':   ['CPA3', 'KIT', 'HPGD', 'ENPP3'],
    'Erythrocytes': ['HBB', 'HBA1', 'HBA2'],
    'Platelets':    ['PPBP', 'PF4', 'GP9'],
}

# --- the T/NK correction ---------------------------------------------------
NK_OLD = ['KLRD1', 'KLRC3', 'GNLY', 'PRF1', 'GZMB', 'NKG7']   # the defect
NK_EXCLUSIVE = ['KLRD1', 'KLRC3', 'NCAM1']                    # the correction
NK_SHARED = ['GNLY', 'PRF1', 'GZMB', 'NKG7']                  # informational
CD3_GENES = ['CD3D', 'CD3E', 'CD3G']

# RAW counts. Raised 1 -> 3 after 05b, then 3 -> 5 after the real
# calibration ran. The synthetic test I sized this on said 3 would leave
# background near zero; the actual sweep on this object said otherwise, and
# the object wins.
#
# Measured on this object, pre-gate, in Erythrocytes and Platelets (1,495
# cells, none of which can express CD3):
#
#     threshold   background CD3+   cells the gate would move
#             1            35.79%                       9,297
#             2            20.60%                       3,609
#             3            12.11%                       1,805
#             5             3.28%                         550   <- chosen
#            10             0.20%                          32
#
# At 3 the gate was moving a population that is 12% ambient, which is not a
# gate, it is a laundering step. 10 drives background to nothing but leaves
# the gate nearly inert at 32 cells. 5 is the compromise, and the tolerance
# below is set to 5.0 to match what is actually achievable on this data
# rather than leaving a warning permanently lit, which only teaches people to
# ignore warnings.
#
# Raising the threshold shifts work from the CD3 gate to the TCR evidence in
# CELL 9, which is the better evidence anyway: a productive receptor is a
# stronger claim than three counts of CD3 in a well full of ambient CD3.
CD3_MIN_COUNTS = 5

# Apply the gate to EVERY non-T label, not just NK and Unassigned.
#
# The previous version gated only those two, because those were where the v10
# collision happened to be known. A CD3-positive cell that argmaxed to B cells,
# Monocytes, Macrophages, DCs, Mast, Erythrocytes or Platelets kept that label
# and was never examined. The gate exists because the scores are untrustworthy,
# and which label a bad score produced is not a reason to trust it.
CD3_GATE_ALL_LABELS = True

# Lineages that cannot express CD3 under any reading. Their CD3-positive rate
# at the chosen threshold IS the background rate, so the log prints it as a
# calibration. If it is not near zero, the threshold is too low and the gate is
# moving soup into the T compartment.
CD3_CALIBRATOR_LINEAGES = ['Erythrocytes', 'Platelets']
# 5.0, not 2.0. The measured floor on this object is 3.28% at threshold 5,
# and 2.0% is only reachable at threshold 10 where the gate stops doing
# anything. A tolerance that can never be met is a warning nobody reads.
CD3_BACKGROUND_TOLERANCE = 5.0      # percent; warn above this

# --- marker support: a label must be backed by its own markers -------------
# Assignment is argmax over score_genes with a floor on the SCORE, not on
# marker detection. score_genes subtracts a control-gene mean, so in a sparse
# cell both terms sit near zero and noise picks the winner.
#
# 05b measured how much this actually costs, and the answer is lineage
# specific. B cells are fine: 0.0% express zero B markers, median 3 of 5. The
# thin ones are NK at 14.3% zero-markers, and Macrophages, DCs and Mast cells
# which are each assigned on a median of ONE marker (CD68 78% but MERTK 3.7%;
# FLT3 78% but CLEC4C 1.9%; KIT 77% but CPA3 11.9%).
#
# A cell failing this becomes Unassigned, which is honest: nothing showed it
# was that lineage. It does NOT become a different lineage, because failing one
# panel is not evidence for another.
REQUIRE_MARKER_SUPPORT = True
MIN_MARKERS_DETECTED = 1
MARKER_SUPPORT_MIN_COUNTS = 1       # RAW counts for "detected"

# Report what this stricter bar would cost, and apply nothing. On the
# 2026-11-02 run the >=1 bar moved only 365 cells, all NK, because every
# lineage had 0.0% of cells at zero markers. The thin lineages are thin in a
# different way: Macrophages, DCs, Mast, Erythrocytes and Platelets are each
# called on a MEDIAN of one detected marker. Raising the bar to 2 would catch
# that, and it would also move a large number of real cells to Unassigned.
# That trade is a decision, so the log prints the cost per label and changes
# nothing. Set MIN_MARKERS_DETECTED = 2 to take it.
REPORT_ALT_MARKER_MIN = 2

# --- VDJ / TCR. Unchanged from shiv_VDJ_integration.py. --------------------
DO_VDJ = True
VDJ_SUBDIR = 'vdj_t'
VDJ_FILE = 'filtered_contig_annotations.csv'
LIB_NAME = {
    'A': 'SHIV_Pre_PBMC',             'B': 'SHIV_Pre_LN',
    'C': 'SHIV_Pre_PBMC_IndianMixed', 'D': 'SHIV_21DPI_Pre_LN_Mixed',
    'E': 'SHIV_21DPI_PBMC',           'F': 'SHIV_21DPI_NonInf_LN_Mixed',
    'G': 'SHIV_Necropsy_PBMC',        'H': 'SHIV_Necropsy_LN',
}
# Strict, paired. TRG persists in alpha-beta T cells and TRD recovery is
# sparse, so "any TRG or TRD" overcounted gamma-delta at 2,982 cells. The
# paired subtype is a floor, which is the honest direction to be wrong in.
GD_SUBTYPES = ['TRG+TRD']
GD_SUBTYPES_LOOSE = ['TRG+TRD', 'TRD', 'TRG']     # reported, never used

# TCR recovery quarantine. A lone TCR chain in these lineages is more likely
# ambient soup than a real receptor: measured orphan rates were B 9.3%,
# Monocytes 10.7%, Macrophages 14.4%. Unassigned, Erythrocytes and Platelets
# are NOT quarantined, because their own HBB/PPBP signal is itself ambient, so
# a detected receptor is the better evidence.
TCR_RECOVERY_QUARANTINE = ['B cells', 'Monocytes', 'Macrophages', 'DCs']

# --- Tier 2. Unchanged from shiv_tier2_diagnostic_v2.py. -------------------
TIER2_LINEAGES = ['T cells', 'B cells', 'Monocytes']
MIN_LINEAGE_CELLS = 300
SUB_HVG = 2000
SUB_PCS = 50
SUB_NEIGH_PCS = 30
SUB_RES = 1.0
N_TOP_MARKERS = 15
# --- Tier 2 assignment: subcluster mean, not per-cell vote ----------------
# The first version assigned by taking the per-cell argmax over z-scored
# signatures, the mode per subcluster, and requiring dominant_frac >= 0.40.
# On the 2026-11-02 run that left 10 of 11 T subclusters Unresolved, including
# two that are CD8 cytotoxic effectors by every other column.
#
# The cause was the z-scoring, not the panels. Dividing each panel by its SD
# across all cells is meant to make panels comparable, but that SD is large
# for a panel with real biological signal and tiny for a panel whose genes are
# barely expressed. Measured on the T lineage, the between-subcluster range of
# each panel was Proliferating 1.70, Cytotoxic 1.39, EffectorMem 1.31, Naive
# 0.75, CentralMem 0.46, then Treg 0.18, Th17 0.14, Tfh 0.11, Th1 0.09. The
# last four carry no information about which subcluster is which, and
# z-scoring rescaled their pure noise up to the same unit variance as the real
# signal, handing them equal voting power. They took 26% of the votes in a
# subcluster whose Cytotoxic mean was 1.042.
#
# score_genes already subtracts a control-gene mean, so a non-expressing panel
# sits near zero by construction and the RAW means are directly comparable.
# Assigning on the raw subcluster mean therefore needs no rescaling, and a
# flat panel cannot out-rank a real one because its mean is small. It also
# makes the log self-consistent: the "mean signature score per subcluster"
# table printed below is now the table the assignment is actually made from,
# which was not true before.
#
# A panel-range floor was considered and REJECTED. It changes not one label
# once assignment is on the mean, so it would have been a second arbitrary
# threshold producing Unresolved cells for no gain.
MIN_SCORE_MARGIN = 0.15        # top family mean minus runner-up
REQUIRE_POSITIVE_TOP = True    # a negative top score cannot name a subtype

# dominant_frac is still computed and reported, but it is NO LONGER A GATE. It
# measures HOMOGENEITY, not confidence: a subcluster can be unambiguously
# cytotoxic at the population level and still split its per-cell votes between
# two co-expressed panels, while a 50/50 mixture of two real populations would
# average to nothing in the mean and only show up in the vote share. The two
# statistics answer different questions and both are worth printing.
REPORT_DOMINANT_FRAC = True

# Panel families. Two panels go in a family when asking which of them a cell
# is cannot be answered from the markers available.
#
#   Naive + CentralMem share CCR7, SELL and IL7R, three of five genes each.
#   The panels are close to the same panel.
#
#   EffectorMem + Cytotoxic share ZERO genes. They are grouped because the
#   cytotoxic program is CO-EXPRESSED, which the data shows directly: T
#   subclusters 7 and 8 scored EffectorMem 0.80/0.80 and Cytotoxic 1.04/1.01
#   while every other subcluster was negative on both. The script prints the
#   between-subcluster correlation of the panel means as the evidence.
#
# B cells and Monocytes get no families: no B panels share a gene, and the one
# Monocyte overlap (NonClassical + Intermediate share FCGR3) is a panel defect
# to fix rather than a family to declare.
#
# These names are deliberately coarser than "CD4 central memory" and "CD8
# effector memory". Coarse and right beats specific and wrong; promote them
# through SUBTYPE_RENAME once the ADT and enrichment tables justify it.
PANEL_FAMILIES = {
    'T cells': {
        'Naive/CM': ['Naive', 'CentralMem'],
        'Effector/Cytotoxic': ['EffectorMem', 'Cytotoxic'],
    },
}

# Rhesus gene naming. A panel gene absent from Mmul10 under its human name
# silently shrinks the panel, and a panel reduced to a non-specific gene
# produces CONFIDENT WRONG answers, which is worse than Unresolved.
#
# The live example: the Monocyte Intermediate panel is CD74 + FCGR3 + HLA-DRA,
# and HLA-DRA does not exist in Mmul10. The rhesus MHC class II alpha gene is
# MAMU-DRA. So the panel ran on CD74 + FCGR3, and CD74 is one of the most
# abundant transcripts in B cells, which is why "Intermediate" scored 1.716 in
# a subcluster whose top markers are CD79B, MS4A1, BANK1 and CD79A.
#
# Each key resolves to the first alias present in the matrix, and the log
# reports which alias was used. Run scripts/05c_panel_audit.py to find
# candidates for the genes still missing (CD27, CD24, IGHD, IL17A, LILRB2,
# PF4, HBA2) rather than guessing at them here.
GENE_ALIASES = {
    'HLA-DRA':  ['MAMU-DRA'],      # confirmed present, 05c 2026-11-03
    'HLA-DRB1': ['MAMU-DRB1'],
    'HLA-DQA1': ['MAMU-DQA1'],
    'HLA-DPA1': ['MAMU-DPA'],
    'PF4':      ['PF4V1'],         # genuine paralog, platelet-expressed
    'HBA2':     ['HBA'],           # second alpha-globin locus; HBA1 already in panel
}

# REJECTED substitutions, recorded so nobody retries them. 05c proposed these
# by pattern matching and every one is a DIFFERENT GENE:
#
#   CD24   -> CD247 (CD3-zeta), CD244 (2B4/SLAMF4), CD248 (endosialin), CD2
#   CD27   -> CD276 (B7-H3), CD274 (PD-L1), CD2
#   LILRB2 -> LILRB1, LILRB3. The LILR family is expanded and divergent in
#             macaques, so 1:1 orthology cannot be assumed and LILRB1 has a
#             broader expression pattern that would make the NonClassical
#             panel LESS specific.
#   IGHD, IL17A -> no candidate found. Immunoglobulin and cytokine loci are
#             poorly annotated; both are likely LOC identifiers, and 25.2% of
#             this matrix is LOC.
#
# So CD27, CD24 and IGHD are NOT recoverable by renaming. They were the three
# discriminating markers for naive, memory and transitional B cells, which is
# why that lineage would not separate. The fix is not an alias, it is RIRA's
# rhesus-validated B modules below.
PANELS = {
    'T cells': {
        'Naive':         ['CCR7', 'SELL', 'TCF7', 'LEF1', 'IL7R'],
        'CentralMem':    ['CCR7', 'SELL', 'IL7R', 'CD28', 'GPR183'],
        'EffectorMem':   ['GZMK', 'CCL5', 'CST7', 'KLRG1'],
        'Cytotoxic':     ['GZMB', 'GNLY', 'PRF1', 'NKG7', 'GZMH', 'FGFBP2'],
        'Th1':           ['TBX21', 'IFNG', 'CXCR3', 'STAT1'],
        'Th17':          ['RORC', 'IL17A', 'CCR6', 'IL23R', 'KLRB1'],
        'Tfh':           ['CXCR5', 'PDCD1', 'BCL6', 'ICOS', 'IL21'],
        'Treg':          ['FOXP3', 'IL2RA', 'CTLA4', 'IKZF2'],
        'Proliferating': ['MKI67', 'TOP2A'],
    },
    # B-cell panels REPLACED 2026-11-03 with the RIRA rhesus-validated
    # modules (Mahyari et al. 2025, Cell Genomics 5:100849, Figures 3E-3K).
    #
    # The human-canonical versions could not work here. CD27, CD24 and IGHD
    # are absent from Mmul10 and are not recoverable by renaming (see the
    # rejected substitutions above), and those were the three discriminating
    # markers for memory, transitional and naive B cells. That is why five of
    # ten B subclusters came out Unresolved and the two that were labelled
    # sat on a margin of 0.15.
    #
    # RIRA hit the same wall and solved it empirically, scoring B cells for
    # modules rather than canonical markers: "while there were few
    # lineage-defining markers within mature B cells, we reasoned that scoring
    # B cells for enrichment of multi-gene modules would provide a
    # classification tool."
    #
    # InnateLike is a NEW panel and it is the one to watch. RIRA describes
    # CD40-low, PBMC-enriched innate-like B cells marked by AHNAK, CD1C and
    # SOX5. B subclusters 0 and 1 in the 2026-11-02 run already carry exactly
    # that: SOX5, AHNAK and CD1C among their top markers, and both are
    # PBMC-enriched at +1.29 and +1.31 with LN at -2.00. Three independent
    # lines agree, so those two subclusters are very likely innate-like B
    # rather than the weak "Memory" call they were getting.
    'B cells': {
        # TCL1A and FCER2 from RIRA; CXCR4 is RIRA's rhesus flow marker.
        # IGHD and IL4R dropped: IGHD absent, IL4R not RIRA-validated.
        'Naive':          ['TCL1A', 'FCER2', 'CXCR4'],
        # RIRA's memory module. TNFRSF13B (TACI) kept from the old panel as a
        # legitimate memory marker that IS present; AIM2 and TXNIP dropped as
        # the weakest and least specific of the old set.
        'Memory':         ['PLAC8', 'CR2', 'LY86', 'CD44', 'TNFRSF13B'],
        'InnateLike':     ['AHNAK', 'CD1C', 'SOX5'],
        # Pre-B / immature. RIRA uses VPREB1 and SOX4; VPREB3 and MME kept.
        'Transitional':   ['VPREB1', 'VPREB3', 'SOX4', 'MME'],
        # RIRA: JCHAIN, IGHA1, XBP1. MZB1, PRDM1, SDC1 kept from the old set.
        'Plasma':         ['JCHAIN', 'XBP1', 'MZB1', 'PRDM1', 'SDC1', 'IGHA1'],
        # RIRA: RGS13, MEF2B, AICDA, "all consistent with published features".
        'GerminalCenter': ['RGS13', 'MEF2B', 'AICDA', 'BCL6'],
        'Activated':      ['CD69', 'CD83', 'NR4A1'],
    },
    # Monocyte panels REBUILT 2026-11-03 on RIRA (Mahyari et al. 2025,
    # Figure 2B-2C), and 'Intermediate' is GONE.
    #
    # Intermediate was CD74 + FCGR3 + HLA-DRA. Resolving HLA-DRA to MAMU-DRA
    # made it worse, not better: all three genes are then MHC-II pathway or
    # Fc, so nothing in the panel is monocyte-specific. It won in 7 of 12
    # subclusters with a between-subcluster range of 2.300, which is a
    # lineage-wide MHC-II gradient wearing a subset label. It mislabelled
    # subcluster 4 (993 cells, top markers MS4A1 CD79B BANK1 AFF3 FCRL1, so B
    # cells) at a margin of 0.845, and subcluster 11 (142 cells, FCER1A plus
    # high MHC-II plus CST3, so cDC2) at a margin of 2.356.
    #
    # RIRA found TWO monocyte populations in 41,199 rhesus myeloid cells, not
    # three: "separation of CD14+ classical monocytes from CD16+ non-classical
    # monocytes". Forcing a third category with no monocyte-specific gene is
    # what produced those mislabels, so the category goes.
    #
    # NonClassical gains RIRA's better RNA markers: "non-classical monocytes
    # expressed both protein and RNA for the canonical marker CD16 (FCGR3);
    # however, additional RNA markers such as CX3CR1 and CD11c (ITGAX) were
    # more specific in scRNA-seq data." LILRB2 dropped: absent from Mmul10 and
    # not recoverable (see the rejected substitutions).
    #
    # The DC panel is split into RIRA's actual lineages, which is what
    # subcluster 11 is. pDC deliberately OMITS GZMB and JCHAIN even though
    # RIRA lists them, because Monocyte subcluster 8 is GZMB/NKG7/GNLY-high
    # NK contamination and a GZMB-bearing pDC panel would capture it. LILRA4
    # is the specific marker and carries the panel on its own merits.
    'Monocytes': {
        'Classical':    ['CD14', 'VCAN', 'S100A8', 'S100A9'],
        'NonClassical': ['FCGR3', 'CX3CR1', 'ITGAX', 'CDKN1C', 'MS4A7'],
        'Macrophage':   ['CD68', 'APOC1', 'C1QA', 'FABP4', 'MRC1', 'CD163'],
        'cDC1':         ['CLEC9A', 'XCR1', 'CADM1'],
        'cDC2':         ['CD1C', 'FCER1A', 'CLEC10A'],
        'pDC':          ['IL3RA', 'LILRA4'],
        'mregDC':       ['CCR7', 'LAMP3', 'CCL22', 'FSCN1'],
    },
}

# Panels that are SCORED AND PRINTED but never participate in assignment.
# Same pattern, and the same reasoning, as dx_cyto_shared in Tier 1: a score
# that is informative to look at is not therefore a label.
#
# RIRA (Mahyari et al. 2025) found the granzyme pattern to be the single most
# informative discriminator among effector T and NK cells, more so than the
# canonical subset markers: "the pattern of granzyme K/M versus granzyme A/B/H
# expression by itself is informative at differentiating between many effector
# T and NK subsets." Their mapping, from sorted reference cells:
#
#     GZMM alone          naive CD8
#     GZMM + GZMK         TCM / TTM
#     GZMB + GZMH         TEM
#     GZMA + GZMB + GZMH  CD16-positive NK
#     GZMM + GZMK         CD16-negative NK
#
# They also found S100A genes to be "a much more reliable signature of T cell
# memory than genes encoding canonical protein markers".
#
# These are NOT in the assignment because they cut across it. Cytotoxicity is
# what dominates the clustering, so Effector/Cytotoxic answers "is this
# cytotoxic" while the granzyme skew answers "TCM-like or TEM-like" WITHIN
# that. A flat argmax cannot express two levels, and forcing it to would
# recreate exactly the shared-program collision that broke the NK panel in
# Tier 1. Read them off the means table and promote through SUBTYPE_RENAME.
TIER2_DIAGNOSTIC_PANELS = {
    'T cells': {
        'dx_GranzymeKM':  ['GZMK', 'GZMM'],
        'dx_GranzymeABH': ['GZMA', 'GZMB', 'GZMH'],
        'dx_S100Amemory': ['S100A4', 'S100A6', 'S100A10', 'S100A11'],
        'dx_Perforin':    ['PRF1', 'GNLY'],
    },
    'B cells': {
        'dx_Proliferating': ['MKI67', 'TOP2A'],
    },
}

# Sub-split a family on a measured diagnostic, added 2026-11-03.
#
# The granzyme skew worked on the first run and it split one label into two
# real populations:
#
#     sub 5   GranzymeKM 0.352   GranzymeABH 1.229   skew +0.878   PBMC +1.13
#             top markers GZMB GZMH, Perforin 0.848, S100Amemory 1.114  -> TEM
#     sub 6   GranzymeKM 1.343   GranzymeABH -0.229  skew -1.572   neutral
#             top markers GZMK GZMM                                  -> TCM/TTM
#
# Both were labelled Effector/Cytotoxic, 10,015 cells between them, and
# SUBTYPE_RENAME could not tell them apart because it is keyed on the
# provisional name. The alternative was to let SUBTYPE_RENAME key on
# subcluster ID, which is the thing removed from this script for being
# unstable across rebuilds. A split on a MEASURED quantity is stable.
#
# RIRA's mapping from sorted reference cells: GZMM alone is naive CD8,
# GZMM+GZMK is TCM/TTM, GZMB+GZMH is TEM, GZMA+GZMB+GZMH is CD16-positive NK,
# GZMM+GZMK is CD16-negative NK. The threshold is well clear of both observed
# values, so it is not tuned to them.
FAMILY_SUBSPLIT = {
    'T cells': {
        'Effector/Cytotoxic': {
            'hi_panel': 'sig_dx_GranzymeABH',
            'lo_panel': 'sig_dx_GranzymeKM',
            'threshold': 0.30,
            'hi_label': 'Effector/Cytotoxic TEM-like',
            'lo_label': 'Effector/Cytotoxic TCM-like',
        },
    },
}

# Promote provisional names to curated ones AFTER reading the review tables.
# Keys are '<lineage>:<provisional>'. Empty by default on purpose: an empty
# dict means the labels you see came from the signatures, not from a guess.
# Your v10 adjudication, for reference when you fill this in:
#   Cytotoxic + ~100% ADT-CD8 + LN-excluded  -> 'CD8 terminal effector'
#   EffectorMem + 90-95% ADT-CD8 + PBMC      -> 'CD8 effector memory'
#   CentralMem + 77-82% ADT-CD4              -> 'CD4 central memory'
#   Naive                                    -> 'Naive/resting T'
SUBTYPE_RENAME = {}

# --- per-cell attributes, deliberately NOT labels --------------------------
# Each of these was argued out in the v10 pass and lost. They stay per-cell.
PROLIF_GENES = ['MKI67', 'TOP2A']
ACTIVATION_GENES = ['FOS', 'FOSB', 'JUN', 'JUNB', 'JUND',
                    'EGR1', 'NR4A1', 'DUSP1']   # immediate-early; soft score
CD8_GENES = ['CD8A', 'CD8B']
CD4_GENES = ['CD4']
POS_MIN = 1

# --- QC rules at Tier 2. Explicit thresholds, not cluster IDs. -------------
QC_MITO_SUBCLUSTER = 15.0      # median pct_counts_mt within the subcluster
QC_DOUBLET_PAIRS = {           # lineage: (own markers, foreign markers)
    'T cells':   (['CD3D', 'CD3E'], ['CD79A', 'CD79B', 'MS4A1']),
    'B cells':   (['CD79A', 'MS4A1'], ['CD3D', 'CD3E', 'CD3G']),
    'Monocytes': (['LYZ', 'CD14'], ['CD3D', 'CD3E', 'CD79A']),
}
QC_DOUBLET_COEXPR_FRAC = 0.50  # fraction of cells co-expressing both sides
QC_DOUBLET_SCORE_MULT = 1.5    # mean doublet_score vs the lineage median

# Everything that decides a Tier 1 label. Stamped onto the checkpoint at
# CELL 11 and checked on resume at CELL 2, so a checkpoint from an older run
# cannot be resumed into a config it does not match.
TIER1_RECIPE = {
    'SCORE_THRESHOLD': SCORE_THRESHOLD,
    'CD3_MIN_COUNTS': CD3_MIN_COUNTS,
    'CD3_GATE_ALL_LABELS': CD3_GATE_ALL_LABELS,
    'REQUIRE_MARKER_SUPPORT': REQUIRE_MARKER_SUPPORT,
    'MIN_MARKERS_DETECTED': MIN_MARKERS_DETECTED,
    'MARKER_SUPPORT_MIN_COUNTS': MARKER_SUPPORT_MIN_COUNTS,
    'TCR_RECOVERY_QUARANTINE': ','.join(TCR_RECOVERY_QUARANTINE),
    'NK_EXCLUSIVE': ','.join(NK_EXCLUSIVE),
}

# --- the graded viral call --------------------------------------------------
# 05b established that >=1 viral UMI is NOT a detection of infection. Of 59
# such cells, 38 carried exactly one molecule, the median was 1, and three sat
# in lineages that cannot be infected at all. The strongest B cell (5 UMI,
# 0.53% of transcriptome) was WEAKER than the strongest erythrocyte (7 UMI,
# 0.58%). So a bare >0 threshold manufactures a tropism story out of soup.
#
# The fix is a background ceiling measured from the cells that cannot be
# infected, which are the only internal negative control in the experiment.
# Above that ceiling in the audit: three T cells and a monocyte, both
# legitimate CCR5-tropic targets. No B cells. No impossible lineages.
#
# SHIV_CONFIDENT_MIN_UMI = None derives the ceiling from the data each run
# (max burden in VIRAL_BACKGROUND_LINEAGES, plus one). Set an integer to pin
# it. The derivation is honest but FRAGILE: it rested on three cells in the
# audit, so one more background hit at 12 UMI would move the ceiling past two
# of the four survivors. The log prints the ceiling and how many cells it was
# derived from, so the fragility is visible rather than implied.
VIRAL_BACKGROUND_LINEAGES = ['Erythrocytes', 'Platelets', 'Mast cells']
SHIV_CONFIDENT_MIN_UMI = None       # None = derive; int = pin
SHIV_CONFIDENT_MIN_FRAC = 0.01      # 1% of the transcriptome

# OR, not AND, and the two criteria answer different questions. Ambient
# deposits a roughly FIXED number of molecules per droplet, so absolute UMI is
# what compares against background regardless of depth. Viral FRACTION is what
# says "productively transcribing". A deep cell with 12 viral UMI is well above
# background even though its fraction is diluted to 0.28%; a shallow cell at 1%
# is at the literature level even with only 6 molecules. Requiring both would
# discard each of those for the wrong reason.
SHIV_CONFIDENT_REQUIRE_BOTH = False

# --- ADT -------------------------------------------------------------------
DO_ADT_CLR = True
ADT_OVERLAY = ['CD8', 'CD4', 'CX3CR1', 'CD28', 'CD69']
# CD163 and CD169 are BD AbSeq 36-mers, structurally incompatible with 10x 5'
# Antibody Capture: they need a 3' poly-A and a BD PCR handle, so those
# molecules are never sequenced into the ADT library. Settled, do not retest.
ADT_KNOWN_DEAD = ['CD163', 'CD169']

# --- metadata columns ------------------------------------------------------
ANIMAL_COL = 'animal_id'
TIME_COL = 'timepoint'
TISSUE_COL = 'tissue'
ELITE_CONTROLLER = '40707'
INDIAN = ['34315', '41903']
TIME_ORDER = ['Pre', '21 DPI', 'Necropsy', 'Non-Infected', 'Unknown']

# --- figures. All hex, 28-34 pt, PDF and PNG at 300 DPI. -------------------
FIG_DPI = 300
plt.rcParams.update({
    'font.size': 28, 'axes.titlesize': 34, 'axes.labelsize': 30,
    'xtick.labelsize': 28, 'ytick.labelsize': 28, 'legend.fontsize': 18,
    'figure.dpi': 100, 'pdf.fonttype': 42, 'ps.fonttype': 42,
})
TIER1_COLORS = {
    'T cells': '#D32F2F', 'B cells': '#1976D2', 'NK cells': '#7B1FA2',
    'Monocytes': '#F9A825', 'Macrophages': '#00897B', 'DCs': '#6D4C41',
    'Mast cells': '#455A64', 'Erythrocytes': '#8D6E63',
    'Platelets': '#C2185B', 'Unassigned': '#E0E0E0',
}
HEX_GRID = '#E0E0E0'
HEX_QC = '#BDBDBD'

os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)
sc.settings.figdir = FIG_DIR
sc.settings.n_jobs = int(os.environ.get('SLURM_CPUS_PER_TASK', 8))
pd.set_option('display.width', 240)
pd.set_option('display.max_columns', 60)
pd.set_option('display.max_rows', 200)


def hr(t):
    print('\n' + '=' * 94, flush=True)
    print(t, flush=True)
    print('=' * 94, flush=True)


def save_csv(df, name, index=True):
    df.to_csv(os.path.join(OUT_DIR, name), index=index, lineterminator='\n')
    print(f'  wrote {name}  ({df.shape[0]:,} rows)', flush=True)


def savefig(fig, stem):
    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(FIG_DIR, f'{stem}.{ext}'), dpi=FIG_DPI,
                    bbox_inches='tight')
    plt.close(fig)
    print(f'  wrote {stem}.pdf / .png', flush=True)


def present(genes, a):
    """Genes actually in this matrix. Rhesus naming drops a few every time."""
    return [g for g in genes if g in a.var_names]


def hvg_exclusion(a):
    """Features to keep OUT of HVG selection. Used by CELL 4 and CELL 12.

    One function so the global and within-lineage steps cannot drift apart,
    which is exactly how the Tier 2 subclustering ended up being driven by
    mitochondrial and ribosomal content while the global embedding was clean.
    """
    mask = np.zeros(a.n_vars, dtype=bool)
    detail = {}
    for flag in EXCLUDE_FROM_HVG:
        if flag in a.var.columns:
            m = a.var[flag].astype(bool).values
            mask |= m
            detail[flag] = int(m.sum())
    pats = dict(EXCLUDE_HVG_PATTERNS)
    if EXCLUDE_MHC_CLASS_II:
        pats['mhc_class_II'] = r'^MAMU-(D[MOPQR])'
    for name, pat in pats.items():
        # np.asarray, not .values. On pandas 3.x, Index.str.match returns a
        # bare ndarray which has no .values; on older pandas it returned
        # something that did. np.asarray is correct for both and for a Series.
        m = np.asarray(a.var_names.str.match(pat), dtype=bool)
        mask |= m
        detail[name] = int(m.sum())
    return mask, detail


def resolve_panel(genes, a):
    """Resolve a panel against the matrix, following GENE_ALIASES.

    Returns (resolved, used_alias, missing). A human gene name absent from
    Mmul10 is tried against its rhesus aliases before being given up on,
    because a silently shortened panel is how a non-specific gene ends up
    carrying a whole subtype call.
    """
    resolved, used, missing = [], {}, []
    for g in genes:
        if g in a.var_names:
            resolved.append(g)
            continue
        hit = next((alt for alt in GENE_ALIASES.get(g, [])
                    if alt in a.var_names), None)
        if hit:
            resolved.append(hit)
            used[g] = hit
        else:
            missing.append(g)
    return resolved, used, missing


def gset_sum(a, genes, layer=None):
    """Summed counts across a gene set. Zeros if none are present."""
    g = present(genes, a)
    if not g:
        return np.zeros(a.n_obs)
    m = a[:, g].layers[layer] if layer else a[:, g].X
    m = m.toarray() if sparse.issparse(m) else np.asarray(m)
    return m.sum(axis=1)


def resolve_adt(adt_names, tag):
    """ADT tags are inconsistently suffixed across the feature reference."""
    if tag in adt_names:
        return adt_names.index(tag)
    for i, n in enumerate(adt_names):
        if n.split('_')[0] == tag:
            return i
    for i, n in enumerate(adt_names):
        if n.startswith(tag):
            return i
    return None


def log2_enrichment(groups, values, clip=2.0):
    """Observed over expected per group x value, log2, clipped."""
    obs = pd.crosstab(groups, values)
    exp = np.outer(obs.sum(axis=1), obs.sum(axis=0)) / obs.values.sum()
    e = np.log2((obs.values + 1) / (exp + 1))
    return obs, pd.DataFrame(np.clip(e, -clip, clip),
                             index=obs.index, columns=obs.columns)


# %% ============================================================
# CELL 2 — LOAD, OR RESUME FROM THE TIER 1 CHECKPOINT
# ============================================================
hr('CELL 2 — load')

RESUMED = False
if RESUME_FROM_TIER1 and os.path.exists(CKPT_H5AD):
    adata = sc.read_h5ad(CKPT_H5AD)
    stamp = dict(adata.uns.get('tier1_recipe', {}))
    drift = {k: (stamp.get(k, '<absent>'), val)
             for k, val in TIER1_RECIPE.items()
             if str(stamp.get(k, '<absent>')) != str(val)}
    if drift:
        print(f'  STALE CHECKPOINT: {os.path.basename(CKPT_H5AD)} was written')
        print('  under a different Tier 1 recipe, so its labels do not match')
        print('  the config at the top of this file:')
        for k, (was, now) in drift.items():
            print(f'    {k}: checkpoint {was}  ->  config {now}')
        raise SystemExit(
            '  STOP: resuming would reapply the OLD gate and silently\n'
            '        produce labels that contradict Cell 1. Either set\n'
            '        RESUME_FROM_TIER1 = False to rebuild Tier 1, or delete\n'
            '        the checkpoint. Refusing rather than guessing.')
    RESUMED = True
    print(f'  RESUMED from {os.path.basename(CKPT_H5AD)}')
    print('  recipe stamp matches Cell 1; Tier 1 labels are current.')
    print('  skipping the embedding, Tier 1 and TCR recovery.')
    print('  Set RESUME_FROM_TIER1 = False to rebuild from adata_qc.h5ad.')
elif RESUME_FROM_TIER1:
    print('  RESUME_FROM_TIER1 is set but no checkpoint exists at')
    print(f'  {CKPT_H5AD}')
    print('  falling through to a full run.')

if not RESUMED:
    if not os.path.exists(IN_H5AD):
        raise SystemExit(f'not found: {IN_H5AD}\nRun 04_ingest_qc.py first.')
    adata = sc.read_h5ad(IN_H5AD)

print(f'  {adata.n_obs:,} cells x {adata.n_vars:,} genes')
print(f'  libraries : {sorted(adata.obs["library"].astype(str).unique())}')
print(f'  animals   : {sorted(adata.obs[ANIMAL_COL].astype(str).unique())}')

if 'predicted_doublet' in adata.obs.columns:
    nd = int(adata.obs['predicted_doublet'].sum())
    print(f'  predicted_doublet = True : {nd:,} cells, CARRIED not removed')
    if nd == 0:
        print('  WARN  zero predicted doublets. If 04 ran with')
        print('        REMOVE_DOUBLETS = True they were already deleted, and')
        print('        the per-library loss was 1,121 from G and ~0 elsewhere.')


# %% ============================================================
# CELL 3 — NORMALIZE (once, and only here)
# ============================================================
if not RESUMED:
    hr('CELL 3 — normalize')

    xmax = float(adata.X.max())
    head = adata.X[:20]
    head = head.toarray() if sparse.issparse(head) else np.asarray(head)
    looks_raw = xmax > 50 and np.allclose(head, np.round(head))

    if COUNTS_LAYER not in adata.layers:
        if not looks_raw:
            raise SystemExit(
                f'  STOP: X does not look like raw counts (max={xmax:.1f}) and\n'
                f'        layers["{COUNTS_LAYER}"] is absent. The CD3 hard gate\n'
                '        and the CD4/CD8 call both need raw counts and cannot\n'
                '        be reconstructed from normalized values.')
        adata.layers[COUNTS_LAYER] = adata.X.copy()
        print(f'  stored raw counts in layers["{COUNTS_LAYER}"] (max={xmax:.0f})')
    else:
        print(f'  layers["{COUNTS_LAYER}"] already present')

    sc.pp.normalize_total(adata, target_sum=TARGET_SUM)
    sc.pp.log1p(adata)
    print(f'  normalize_total(target_sum={TARGET_SUM:.0f}) then log1p')
    print(f'  X max now {float(adata.X.max()):.2f}')


# %% ============================================================
# CELL 4 — HIGHLY VARIABLE GENES (virus and mito excluded)
# ============================================================
if not RESUMED:
    hr('CELL 4 — HVG')

    exclude, detail = hvg_exclusion(adata)
    print(f'  excluded from HVG selection: {int(exclude.sum())} features')
    for k, v in detail.items():
        print(f'    {k:<18} {v:>5}')
    print('  HVG SELECTION ONLY. Every gene stays in the object and is still')
    print('  scored, plotted and ranked, so these will still appear in the')
    print('  top-marker tables. They just no longer define the axes.')
    print('  the virus is excluded so ~59 infected cells cannot pull on the')
    print('  global embedding; mito genes track dissociation stress, not')
    print('  cell identity. Both stay in the object and are still scored.')

    sub = adata[:, ~exclude].copy()
    if HVG_FLAVOR == 'seurat_v3':
        sc.pp.highly_variable_genes(sub, n_top_genes=N_TOP_HVG,
                                    flavor='seurat_v3', layer=COUNTS_LAYER)
    else:
        sc.pp.highly_variable_genes(sub, n_top_genes=N_TOP_HVG,
                                    flavor=HVG_FLAVOR)
    hvg = sub.var_names[sub.var['highly_variable']].tolist()
    del sub

    adata.var['highly_variable'] = adata.var_names.isin(hvg)
    n_hvg = int(adata.var['highly_variable'].sum())
    print(f'  HVGs selected: {n_hvg:,} '
          f'(flavor={HVG_FLAVOR}, n_top_genes={N_TOP_HVG})')
    if n_hvg < N_TOP_HVG * 0.8:
        print('  WARN  fewer HVGs than requested; check the gene filter in 04')


# %% ============================================================
# CELL 5 — PCA, BATCH CORRECTION, UMAP, CLUSTERS
# ============================================================
if not RESUMED:
    hr('CELL 5 — embedding')

    sc.tl.pca(adata, n_comps=N_PCS_PCA, use_highly_variable=True,
              random_state=RANDOM_STATE)
    print(f'  PCA: {N_PCS_PCA} components')

    if DO_PREBBKNN_UMAP:
        sc.pp.neighbors(adata, n_pcs=N_PCS_NEIGHBORS,
                        random_state=RANDOM_STATE)
        sc.tl.umap(adata, random_state=RANDOM_STATE)
        adata.obsm['X_umap_prebbknn'] = adata.obsm['X_umap'].copy()
        print('  pre-correction UMAP stored as obsm["X_umap_prebbknn"]')
        print('  compare it against the corrected one. If the uncorrected')
        print('  version already mixes, BBKNN is doing little and the batch')
        print('  key is worth revisiting rather than trusting.')

    print(f'  BBKNN on batch_key="{BATCH_KEY}" '
          f'({adata.obs[BATCH_KEY].nunique()} batches), '
          f'n_pcs={N_PCS_NEIGHBORS}')
    try:
        import bbknn
        bbknn.bbknn(adata, batch_key=BATCH_KEY, n_pcs=N_PCS_NEIGHBORS)
    except ImportError:
        print('  WARN  bbknn not installed; uncorrected neighbors instead')
        sc.pp.neighbors(adata, n_pcs=N_PCS_NEIGHBORS,
                        random_state=RANDOM_STATE)

    sc.tl.umap(adata, random_state=RANDOM_STATE)
    try:
        sc.tl.leiden(adata, key_added='clusters', resolution=LEIDEN_RES,
                     random_state=RANDOM_STATE, flavor='igraph',
                     n_iterations=2, directed=False)
    except TypeError:
        sc.tl.leiden(adata, key_added='clusters', resolution=LEIDEN_RES,
                     random_state=RANDOM_STATE)
    print(f'  Leiden: {adata.obs["clusters"].nunique()} clusters '
          f'at resolution {LEIDEN_RES}')
    print('\n  cells per cluster:')
    print(adata.obs['clusters'].value_counts().sort_index().to_string())


# %% ============================================================
# CELL 6 — ADT CLR NORMALIZATION
# ============================================================
if not RESUMED:
    hr('CELL 6 — CITE-seq')

    if DO_ADT_CLR and 'protein_counts' in adata.obsm:
        adt_names = list(adata.uns.get('adt_names', []))
        P = np.asarray(adata.obsm['protein_counts'], dtype=float)
        print(f'  {P.shape[1]} antibody tags, {P.shape[0]:,} cells')

        # CLR per cell: log1p(x) minus the per-cell mean of log1p(x).
        L = np.log1p(P)
        adata.obsm['protein_clr'] = L - L.mean(axis=1, keepdims=True)
        for i, nm in enumerate(adt_names):
            adata.obs[f'ADT_{nm}'] = adata.obsm['protein_clr'][:, i]
        print('  CLR written to obsm["protein_clr"] and obs["ADT_<tag>"]')

        det = pd.DataFrame({
            'tag': adt_names,
            'total_counts': P.sum(axis=0).astype(int),
            'cells_detected': (P > 0).sum(axis=0).astype(int),
            'pct_cells': (100 * (P > 0).mean(axis=0)).round(2),
        })
        det['status'] = np.where(det.tag.isin(ADT_KNOWN_DEAD),
                                 'structurally dead (BD AbSeq)', 'usable')
        print()
        print(det.to_string(index=False))
        save_csv(det, 'adt_detection.csv', index=False)

        alive = det[(det.tag.isin(ADT_KNOWN_DEAD)) & (det.total_counts > 0)]
        if len(alive):
            print('\n  NOTE  a BD AbSeq tag has non-zero counts, which')
            print('        contradicts the settled finding that they are')
            print('        never sequenced into the ADT library. Worth a look.')
    else:
        print('  no protein_counts in obsm; skipping')


# %% ============================================================
# CELL 7 — TIER 1 BASELINE SCORING
# ============================================================
if not RESUMED:
    hr('CELL 7 — Tier 1 baseline')

    rows = []
    for lin, genes in TIER1_MARKERS.items():
        found = present(genes, adata)
        rows.append({'lineage': lin, 'n_requested': len(genes),
                     'n_found': len(found), 'genes_found': ';'.join(found),
                     'genes_missing': ';'.join(
                         [g for g in genes if g not in found])})
        key = f"dx_{lin.replace(' ', '_')}"
        if found:
            sc.tl.score_genes(adata, gene_list=found, score_name=key,
                              random_state=RANDOM_STATE)
        else:
            adata.obs[key] = 0.0
            print(f'  WARN  {lin}: no markers present; score forced to 0')
    panels = pd.DataFrame(rows)
    print(panels[['lineage', 'n_requested', 'n_found',
                  'genes_missing']].to_string(index=False))
    save_csv(panels, 'tier1_marker_panels.csv', index=False)

    score_cols = [f"dx_{l.replace(' ', '_')}" for l in TIER1_MARKERS]
    S = adata.obs[score_cols].values
    lab = np.array(list(TIER1_MARKERS))[S.argmax(axis=1)]
    lab = np.where(S.max(axis=1) < SCORE_THRESHOLD, 'Unassigned', lab)
    adata.obs['tier1_celltype_uncorrected'] = pd.Categorical(lab)
    print(f'\n  argmax assignment, Unassigned below {SCORE_THRESHOLD}')
    print(adata.obs['tier1_celltype_uncorrected'].value_counts().to_string())


# %% ============================================================
# CELL 8 — THE T/NK CORRECTION
# ============================================================
if not RESUMED:
    hr('CELL 8 — T/NK correction')

    nk_excl = present(NK_EXCLUSIVE, adata)
    nk_shared = present(NK_SHARED, adata)
    cd3 = present(CD3_GENES, adata)
    print(f'  NK exclusive anchors : {nk_excl}')
    print(f'  shared cytotoxic     : {nk_shared}   (informational only)')
    print(f'  CD3 genes            : {cd3}')
    if not nk_excl:
        raise SystemExit('  STOP: no exclusive NK anchor present. The '
                         'correction cannot be made.')
    if not cd3:
        raise SystemExit('  STOP: no CD3 gene present. The hard gate cannot '
                         'be applied and T cells cannot be called.')

    sc.tl.score_genes(adata, gene_list=nk_excl, score_name='dx_nk_excl',
                      random_state=RANDOM_STATE)
    if nk_shared:
        sc.tl.score_genes(adata, gene_list=nk_shared,
                          score_name='dx_cyto_shared',
                          random_state=RANDOM_STATE)
    adata.obs['dx_nk_old'] = adata.obs['dx_NK_cells']

    # The CD3 hard gate, on RAW counts.
    cd3_counts = gset_sum(adata, cd3, layer=COUNTS_LAYER)
    adata.obs['cd3_counts'] = cd3_counts
    adata.obs['dx_cd3_pos'] = cd3_counts >= CD3_MIN_COUNTS
    print(f'\n  CD3 hard gate at >={CD3_MIN_COUNTS} raw count(s): '
          f'{int(adata.obs["dx_cd3_pos"].sum()):,} cells CD3-positive '
          f'({100 * adata.obs["dx_cd3_pos"].mean():.1f}%)')

    # Re-assign with the corrected NK score. The shared cytotoxic block is
    # deliberately NOT in this matrix. That is the whole fix.
    corr_cols = ['dx_nk_excl' if l == 'NK cells'
                 else f"dx_{l.replace(' ', '_')}" for l in TIER1_MARKERS]
    S2 = adata.obs[corr_cols].values
    lab2 = np.array(list(TIER1_MARKERS))[S2.argmax(axis=1)]
    lab2 = np.where(S2.max(axis=1) < SCORE_THRESHOLD, 'Unassigned', lab2)

    # --- CD3 threshold calibration, BEFORE the gate is applied -------------
    # Erythrocytes and Platelets cannot express CD3, so their CD3-positive
    # rate at any threshold is the background rate. Printed as a sweep so the
    # chosen threshold is visibly a measurement rather than a preference.
    cal = [c for c in CD3_CALIBRATOR_LINEAGES if (lab2 == c).any()]
    if cal:
        cal_mask = np.isin(lab2, cal)
        print(f'\n  CD3 background calibration on {cal} '
              f'({int(cal_mask.sum()):,} cells, none can express CD3):')
        print(f'    {"threshold":>9}  {"background CD3+":>15}  '
              f'{"CD3+ in other non-T labels":>27}')
        other = ~np.isin(lab2, cal + ['T cells', 'NK cells', 'Unassigned'])
        for thr in (1, 2, 3, 5, 10):
            bg = 100 * (cd3_counts[cal_mask] >= thr).mean()
            oth = int(((cd3_counts >= thr) & other).sum())
            mark = '  <- CD3_MIN_COUNTS' if thr == CD3_MIN_COUNTS else ''
            print(f'    {thr:>9}  {bg:>14.2f}%  {oth:>27,}{mark}')
        bg_now = 100 * (cd3_counts[cal_mask] >= CD3_MIN_COUNTS).mean()
        if bg_now > CD3_BACKGROUND_TOLERANCE:
            print(f'\n  WARN  background CD3+ is {bg_now:.2f}% at the chosen')
            print(f'        threshold, above the {CD3_BACKGROUND_TOLERANCE}% '
                  'tolerance. The gate is')
            print('        moving ambient CD3 into the T compartment. Raise')
            print('        CD3_MIN_COUNTS until this is near zero.')
        else:
            print(f'\n  background CD3+ at the chosen threshold: '
                  f'{bg_now:.2f}%  (OK)')
    else:
        print('\n  WARN  no calibrator lineage present; the CD3 threshold is')
        print('        uncalibrated and CD3_MIN_COUNTS is a guess.')

    # --- marker support ----------------------------------------------------
    # A cell must detect at least MIN_MARKERS_DETECTED of its assigned
    # lineage's own markers. Failing becomes Unassigned, never another
    # lineage: failing one panel is not evidence for a different one.
    if REQUIRE_MARKER_SUPPORT:
        n_det = np.zeros(adata.n_obs, dtype=int)
        print(f'\n  marker support: >={MIN_MARKERS_DETECTED} marker(s) of the '
              f'assigned lineage at >={MARKER_SUPPORT_MIN_COUNTS} raw count')
        rows = []
        for lin, genes in TIER1_MARKERS.items():
            m = lab2 == lin
            if not m.any():
                continue
            g = present(genes, adata)
            if not g:
                continue
            C = adata[m, g].layers[COUNTS_LAYER]
            C = C.toarray() if sparse.issparse(C) else np.asarray(C)
            d = (C >= MARKER_SUPPORT_MIN_COUNTS).sum(axis=1)
            n_det[m] = d
            rows.append({'label': lin, 'n_cells': int(m.sum()),
                         'n_markers': len(g),
                         'median_detected': float(np.median(d)),
                         'pct_zero': round(100 * (d == 0).mean(), 1),
                         'would_fail': int((d < MIN_MARKERS_DETECTED).sum())})
        msup = pd.DataFrame(rows)
        print(msup.to_string(index=False))
        save_csv(msup, 'tier1_marker_support.csv', index=False)

        # Report what a stricter bar WOULD cost, and apply nothing. Several
        # lineages are called on a median of ONE detected marker, which is
        # thin, but raising the bar moves real cells to Unassigned and that
        # trade is a decision rather than a fix. The column makes it visible
        # so it can be made deliberately.
        alt = REPORT_ALT_MARKER_MIN
        if alt and alt > MIN_MARKERS_DETECTED:
            extra = int(((n_det < alt) & (n_det >= MIN_MARKERS_DETECTED)
                         & (lab2 != 'Unassigned')).sum())
            print(f'\n  if MIN_MARKERS_DETECTED were {alt} instead of '
                  f'{MIN_MARKERS_DETECTED}, a further {extra:,} cells would')
            print('  move to Unassigned. NOT applied. Per label:')
            per = []
            for lin2 in TIER1_MARKERS:
                m2 = lab2 == lin2
                if not m2.any():
                    continue
                per.append({'label': lin2, 'n_cells': int(m2.sum()),
                            f'would_fail_at_{alt}': int((n_det[m2] < alt).sum()),
                            'pct': round(100 * (n_det[m2] < alt).mean(), 1)})
            print(pd.DataFrame(per).to_string(index=False))

        fails = (n_det < MIN_MARKERS_DETECTED) & (lab2 != 'Unassigned')
        adata.obs['tier1_markers_detected'] = n_det
        adata.obs['tier1_marker_support_failed'] = fails
        print(f'\n  moved to Unassigned for want of marker support: '
              f'{int(fails.sum()):,}')
        print('  05b measured where this bites: B cells 0.0% zero-markers')
        print('  (median 3 of 5) so they are untouched, while NK at 14.3%')
        print('  and Macrophages, DCs and Mast at a median of ONE marker are')
        print('  where the thin calls live.')
        lab2 = np.where(fails, 'Unassigned', lab2)
    else:
        adata.obs['tier1_markers_detected'] = -1
        adata.obs['tier1_marker_support_failed'] = False

    # --- the gate, now applied to every non-T label ------------------------
    cd3_pos = adata.obs['dx_cd3_pos'].values
    n_nk_to_t = int(((lab2 == 'NK cells') & cd3_pos).sum())
    lab2 = np.where((lab2 == 'NK cells') & cd3_pos, 'T cells', lab2)
    n_un_to_t = int(((lab2 == 'Unassigned') & cd3_pos).sum())
    lab2 = np.where((lab2 == 'Unassigned') & cd3_pos, 'T cells', lab2)

    n_other_to_t = 0
    if CD3_GATE_ALL_LABELS:
        other = (~np.isin(lab2, ['T cells', 'Unassigned'])) & cd3_pos
        if other.any():
            print('\n  CD3-positive cells in other non-T labels, by label:')
            print(pd.Series(lab2[other]).value_counts().to_string())
        n_other_to_t = int(other.sum())
        lab2 = np.where(other, 'T cells', lab2)

    adata.obs['tier1_celltype'] = pd.Categorical(lab2)
    adata.obs['nk_provisional'] = (lab2 == 'NK cells') & ~cd3_pos

    print(f'\n  CD3-positive cells called NK, moved to T : {n_nk_to_t:,}')
    print(f'  CD3-positive Unassigned recovered to T   : {n_un_to_t:,}')
    print(f'  CD3-positive OTHER labels moved to T     : {n_other_to_t:,}'
          f'   {"(gate generalized)" if CD3_GATE_ALL_LABELS else "(gate OFF)"}')
    print(f'  CD3-negative NK flagged provisional      : '
          f'{int(adata.obs["nk_provisional"].sum()):,}')
    print('\n  corrected Tier 1:')
    print(adata.obs['tier1_celltype'].value_counts().to_string())

    trans = pd.crosstab(adata.obs['tier1_celltype_uncorrected'],
                        adata.obs['tier1_celltype'])
    print('\n  transition matrix (rows = uncorrected, cols = corrected):')
    print(trans.to_string())
    save_csv(trans, 'tier1_transition_matrix.csv')

    old_nk = (adata.obs['tier1_celltype_uncorrected'] == 'NK cells').values
    if int(old_nk.sum()):
        # The v10 benchmark of 62.2% was measured at CD3 >= 1. Printing a
        # threshold-5 number against it compares two different measurements
        # and makes a threshold change look like a biological shift, which is
        # exactly the misreading this line was supposed to prevent. Report
        # both: the at-1 figure for comparability, and the at-threshold
        # figure for what the gate actually acts on.
        at1 = 100 * (cd3_counts[old_nk] >= 1).mean()
        atN = 100 * (cd3_counts[old_nk] >= CD3_MIN_COUNTS).mean()
        print('\n  UNCORRECTED NK calls that are CD3-positive:')
        print(f'    at CD3 >= 1                  : {at1:.1f}%   '
              '<- compare this one to v10')
        print(f'    at CD3 >= {CD3_MIN_COUNTS} (the gate)       : {atN:.1f}%')
        print('  v10 measured 62.2% AT A THRESHOLD OF 1. Compare like with')
        print('  like: a departure in the at-1 figure means something moved')
        print('  in the reference change. A difference in the at-threshold')
        print('  figure is just the threshold.')


# %% ============================================================
# CELL 9 — VDJ: ATTACH TCR, CALL GAMMA-DELTA, RESOLVE PROVISIONAL NK
# ============================================================
if not RESUMED:
    hr('CELL 9 — VDJ / TCR')

    adata.obs['has_tcr'] = False
    adata.obs['receptor_subtype'] = 'none'
    adata.obs['chain_pairing'] = 'none'

    if not DO_VDJ:
        print('  DO_VDJ = False; skipping')
    else:
        def find_vdj(sample):
            p = os.path.join(MULTI_DIR, sample, 'outs', 'per_sample_outs',
                             sample, VDJ_SUBDIR, VDJ_FILE)
            if os.path.exists(p):
                return p
            hits = glob.glob(os.path.join(MULTI_DIR, sample, 'outs', '**',
                                          VDJ_SUBDIR, VDJ_FILE),
                             recursive=True)
            return hits[0] if hits else None

        try:
            import scirpy as ir
            print(f'  scirpy {ir.__version__}')

            # GEX barcodes are '<bare>-1-<lib>'; VDJ are '<bare>-1'. Append
            # the library suffix so the two align, then subset AIRR to cells
            # present in GEX so the modalities stay aligned.
            parts = []
            for lib in sorted(LIB_NAME):
                path = find_vdj(LIB_NAME[lib])
                if path is None:
                    print(f'  library {lib}: no {VDJ_FILE} found, skipped')
                    continue
                a = ir.io.read_10x_vdj(path)
                a.obs_names = [f'{bc}-{lib}' for bc in a.obs_names]
                a.obs['library'] = lib
                parts.append(a)
                print(f'  library {lib}: {a.n_obs:,} TCR cells')

            if not parts:
                print('  WARN  no VDJ data found anywhere under')
                print(f'        {MULTI_DIR}')
                print('        has_tcr stays False for every cell, so the')
                print('        recovery in CELL 10 is a no-op and every')
                print('        provisional NK stays provisional.')
            else:
                import anndata as ad
                airr = ad.concat(parts, join='outer')
                airr.obs_names_make_unique()
                keep = airr.obs_names.intersection(adata.obs_names)
                print(f'\n  TCR cells total      : {airr.n_obs:,}')
                print(f'  with a GEX partner   : {len(keep):,} '
                      f'({100 * len(keep) / max(airr.n_obs, 1):.1f}%)')
                print('  the rest had their GEX partner removed by QC and are')
                print('  dropped here to keep the modalities aligned.')
                airr = airr[keep].copy()

                ir.tl.chain_qc(airr)
                for col in ('receptor_subtype', 'chain_pairing'):
                    if col in airr.obs.columns:
                        adata.obs[col] = (
                            airr.obs[col].astype(str)
                            .reindex(adata.obs_names, fill_value='none')
                            .values)
                adata.obs['has_tcr'] = adata.obs_names.isin(keep)

                print(f'\n  has_tcr: {int(adata.obs["has_tcr"].sum()):,} cells '
                      f'({100 * adata.obs["has_tcr"].mean():.1f}% of GEX)')
                print('\n  receptor_subtype distribution:')
                print(adata.obs['receptor_subtype'].value_counts().to_string())

                # Gamma-delta, strict and paired.
                adata.obs['is_gamma_delta'] = (
                    adata.obs['receptor_subtype'].isin(GD_SUBTYPES))
                loose = int(adata.obs['receptor_subtype']
                            .isin(GD_SUBTYPES_LOOSE).sum())
                print(f'\n  gamma-delta strict {GD_SUBTYPES}: '
                      f'{int(adata.obs["is_gamma_delta"].sum()):,}')
                print(f'  gamma-delta loose  {GD_SUBTYPES_LOOSE}: {loose:,}')
                print('  the strict call is used. TRG persists in alpha-beta T')
                print('  cells and TRD recovery is sparse, so loose overcounts.')
                print('  Strict is a floor, which is the honest direction to')
                print('  be wrong in.')

                # Resolve provisional NK. Asymmetric on purpose: a productive
                # TCR reassigns to T (strong evidence); no TCR keeps the cell
                # provisional, because TCR dropout means absence is only
                # SUPPORTING evidence for NK, never proof.
                prov = adata.obs['nk_provisional'].values
                resolve = prov & adata.obs['has_tcr'].values
                n_res = int(resolve.sum())
                lab = adata.obs['tier1_celltype'].astype(str).values.copy()
                lab[resolve] = 'T cells'
                adata.obs['tier1_celltype'] = pd.Categorical(lab)
                adata.obs['nk_provisional'] = prov & ~resolve
                print(f'\n  provisional NK with a TCR, reassigned to T: {n_res:,}')
                print(f'  provisional NK remaining (no TCR, stays NK): '
                      f'{int(adata.obs["nk_provisional"].sum()):,}')

        except ImportError:
            print('  WARN  scirpy not installed; TCR step skipped entirely.')
            print('        has_tcr stays False, CELL 10 is a no-op, and every')
            print('        provisional NK stays provisional.')


# %% ============================================================
# CELL 10 — TCR RECOVERY OF CD3-DROPOUT T CELLS
# ============================================================
if not RESUMED:
    hr('CELL 10 — CD3-dropout T recovery')

    adata.obs['tier1_celltype_preTCR'] = adata.obs['tier1_celltype'].copy()

    if not bool(adata.obs['has_tcr'].any()):
        print('  no TCR data attached; nothing to recover.')
        adata.obs['tcr_recovered'] = False
        adata.obs['marker_support_exempt'] = False
    else:
        cur = adata.obs['tier1_celltype'].astype(str).values
        flip = (adata.obs['has_tcr'].values
                & (cur != 'T cells')
                & ~pd.Series(cur).isin(TCR_RECOVERY_QUARANTINE).values)
        print(f'  quarantined lineages (ambient orphan-TCR risk): '
              f'{TCR_RECOVERY_QUARANTINE}')
        print('  measured orphan rates were B 9.3%, Monocytes 10.7%,')
        print('  Macrophages 14.4%, so a lone chain there is more likely soup')
        print('  than a receptor. Unassigned, Erythrocytes and Platelets are')
        print('  NOT quarantined: their own HBB/PPBP is itself ambient, so a')
        print('  detected receptor is the better evidence.')

        if flip.any():
            print('\n  flips by original label:')
            print(pd.Series(cur[flip]).value_counts().to_string())
        new = cur.copy()
        new[flip] = 'T cells'
        adata.obs['tier1_celltype'] = pd.Categorical(new)
        adata.obs['tcr_recovered'] = flip
        print(f'\n  recovered to T cells: {int(flip.sum()):,}')

        # These cells are T on TCR evidence and have NO CD3 by construction,
        # which is the entire point of the recovery: CD3 drops out hardest in
        # effector and cytotoxic T. 05b measured 10.8% of final T cells at
        # zero CD3 markers, which matched the recovered count almost exactly.
        #
        # So the marker-support check in CELL 8 must never be re-applied after
        # this cell, or it would strip exactly the cells this recovery exists
        # to rescue. The ordering is the safeguard (support runs before
        # recovery, once); this flag makes the exemption explicit so a later
        # stage that re-checks support honours it instead of rediscovering the
        # problem.
        adata.obs['marker_support_exempt'] = flip
        adata.obs.loc[flip, 'tier1_marker_support_failed'] = False
        print('  flagged marker_support_exempt: these are T on TCR evidence')
        print('  and carry no CD3 by design. Never re-apply the CELL 8')
        print('  marker-support check after this point.')

        trans2 = pd.crosstab(adata.obs['tier1_celltype_preTCR'],
                             adata.obs['tier1_celltype'])
        print('\n  transition matrix (rows = pre-TCR, cols = post-TCR):')
        print(trans2.to_string())
        save_csv(trans2, 'tier1_tcr_recovery_matrix.csv')

    print('\n  FINAL Tier 1:')
    print(adata.obs['tier1_celltype'].value_counts().to_string())


# %% ============================================================
# CELL 11 — CHECKPOINT
# ============================================================
if not RESUMED:
    hr('CELL 11 — checkpoint')
    # Stamp the Tier 1 recipe onto the object. A checkpoint written by an
    # earlier version of this script carries different labels, and resuming
    # from one would silently reapply the OLD gate while the config at the top
    # of this file says otherwise. CELL 2 refuses a stale checkpoint rather
    # than letting that happen.
    adata.uns['tier1_recipe'] = TIER1_RECIPE
    print('  stamped uns["tier1_recipe"]:')
    for k, val in TIER1_RECIPE.items():
        print(f'    {k} = {val}')
    adata.write(CKPT_H5AD)
    print(f'  wrote {os.path.basename(CKPT_H5AD)}  ({adata.n_obs:,} cells)')
    print('  Set RESUME_FROM_TIER1 = True to iterate on Tier 2 without')
    print('  re-running the embedding.')


# %% ============================================================
# CELL 11b — VERIFY WHAT TIER 2 NEEDS (runs on BOTH paths)
# ============================================================
# This cell is unconditional on purpose. In 04 an earlier version nested the
# viral verification inside `if not REUSED:`, which raised a NameError on the
# reuse path and, worse, SILENTLY SKIPPED the host TAT guard. Everything Tier 2
# depends on is re-derived from the object here rather than inherited from a
# variable that may never have been assigned.
hr('CELL 11b — verify')

required = ['tier1_celltype', 'clusters', 'library', ANIMAL_COL, TIME_COL]
missing = [c for c in required if c not in adata.obs.columns]
if missing:
    raise SystemExit(f'  STOP: obs is missing {missing}. The checkpoint is '
                     'from an older version of this script; re-run with '
                     'RESUME_FROM_TIER1 = False.')
if COUNTS_LAYER not in adata.layers:
    raise SystemExit(f'  STOP: layers["{COUNTS_LAYER}"] is absent. The '
                     'CD4/CD8 call needs raw counts.')

adt_names = list(adata.uns.get('adt_names', []))
HAS_ADT = 'protein_clr' in adata.obsm and len(adt_names) > 0
HAS_TCR = bool(adata.obs.get('has_tcr', pd.Series(False)).any())
HAS_DBL = ('doublet_score' in adata.obs.columns
           and float(np.nanmax(adata.obs['doublet_score'])) > 0)

print(f'  tier1_celltype   : {adata.obs["tier1_celltype"].nunique()} labels')
print(f'  clusters         : {adata.obs["clusters"].nunique()}')
print(f'  ADT available    : {HAS_ADT} ({len(adt_names)} tags)')
print(f'  TCR available    : {HAS_TCR}')
print(f'  doublet scores   : {HAS_DBL}')
if not HAS_ADT:
    print('  WARN  no ADT. The CD4/CD8 call falls back to RNA only, which is')
    print('        much noisier, and ADT_* overlays are skipped.')
if not HAS_DBL:
    print('  WARN  no usable doublet_score. The QC: doublet-like rule will')
    print('        rely on marker co-expression alone, which is weaker.')

if 'pct_counts_mt' not in adata.obs.columns:
    print('  WARN  pct_counts_mt absent; the QC: high-mito rule is skipped.')


# %% ============================================================
# CELL 12 — TIER 2: PER-LINEAGE SUBCLUSTERING AND SIGNATURES
# ============================================================
hr('CELL 12 — Tier 2 subclustering')

tier2_global = adata.obs['tier1_celltype'].astype(str).copy()
adata.obs['sub_leiden'] = 'NA'
adata.obs['tier2_subtype'] = tier2_global.values
adata.obs['tier2_dominant_frac'] = np.nan
adata.obs['qc_exclude'] = False
adata.obs['qc_reason'] = ''

for lin in TIER2_LINEAGES:
    panel = PANELS.get(lin, {})
    if not panel:
        print(f'\n  {lin}: no panel defined, skipped')
        continue

    mask = (adata.obs['tier1_celltype'].astype(str) == lin).values
    n = int(mask.sum())
    print(f'\n{"-" * 94}\n  {lin}: {n:,} cells')
    if n < MIN_LINEAGE_CELLS:
        print(f'  fewer than MIN_LINEAGE_CELLS={MIN_LINEAGE_CELLS}; '
              'left at the Tier 1 label.')
        continue

    sub = adata[mask].copy()
    tag = lin.split()[0]

    # --- within-lineage recluster. HVG are recomputed INSIDE the lineage,
    #     because genes that separate T from B are useless for separating
    #     naive T from effector T.
    #
    # THE EXCLUSIONS APPLY HERE TOO. They did not before 2026-11-03, and that
    # single omission is what put 48% of T and B cells into subclusters
    # defined by ribosomal or mitochondrial content. CELL 4 was clean the
    # whole time, which is why the global embedding looked fine while Tier 2
    # did not.
    if EXCLUDE_HVG_IN_LINEAGE:
        sx, sdetail = hvg_exclusion(sub)
        print(f'  excluded from within-lineage HVG: {int(sx.sum())} features '
              + ', '.join(f'{k} {v}' for k, v in sdetail.items() if v))
        shv = sub[:, ~sx].copy()
        sc.pp.highly_variable_genes(shv, n_top_genes=SUB_HVG, flavor='seurat')
        keep_hvg = shv.var_names[shv.var['highly_variable']].tolist()
        del shv
        sub.var['highly_variable'] = sub.var_names.isin(keep_hvg)
    else:
        sc.pp.highly_variable_genes(sub, n_top_genes=SUB_HVG, flavor='seurat')
    print(f'  HVG within lineage: {int(sub.var["highly_variable"].sum()):,}')
    sc.tl.pca(sub, n_comps=min(SUB_PCS, sub.n_obs - 1),
              use_highly_variable=True, random_state=RANDOM_STATE)
    try:
        import bbknn
        bbknn.bbknn(sub, batch_key=BATCH_KEY, n_pcs=SUB_NEIGH_PCS)
    except Exception as e:                                          # noqa: BLE001
        print(f'  WARN  bbknn unavailable ({type(e).__name__}); '
              'uncorrected neighbors')
        sc.pp.neighbors(sub, n_pcs=SUB_NEIGH_PCS, random_state=RANDOM_STATE)
    sc.tl.umap(sub, random_state=RANDOM_STATE)
    try:
        sc.tl.leiden(sub, key_added='sub_leiden', resolution=SUB_RES,
                     random_state=RANDOM_STATE, flavor='igraph',
                     n_iterations=2, directed=False)
    except TypeError:
        sc.tl.leiden(sub, key_added='sub_leiden', resolution=SUB_RES,
                     random_state=RANDOM_STATE)
    print(f'  sub-Leiden at {SUB_RES}: '
          f'{sub.obs["sub_leiden"].nunique()} subclusters')

    # --- signature scores, resolving rhesus gene names via GENE_ALIASES
    score_cols = []
    for name, genes in panel.items():
        g, used, miss = resolve_panel(genes, sub)
        if g:
            sc.tl.score_genes(sub, gene_list=g, score_name=f'sig_{name}',
                              random_state=RANDOM_STATE)
            score_cols.append(f'sig_{name}')
            line = f'    {name:<14} {len(g)}/{len(genes)}'
            if used:
                line += '   ALIAS ' + ', '.join(f'{k}->{v}'
                                                for k, v in used.items())
            if miss:
                line += f'   MISSING {miss}'
            print(line)
        else:
            print(f'    {name:<14} WARN no markers present, skipped')
    if not score_cols:
        print('  no signature could be scored; left at the Tier 1 label.')
        continue
    # Diagnostic panels: scored and printed, never assigned from.
    diag_cols = []
    for name, genes in TIER2_DIAGNOSTIC_PANELS.get(lin, {}).items():
        g, used, miss = resolve_panel(genes, sub)
        if g:
            sc.tl.score_genes(sub, gene_list=g, score_name=f'sig_{name}',
                              random_state=RANDOM_STATE)
            diag_cols.append(f'sig_{name}')
            print(f'    {name:<16} {len(g)}/{len(genes)}  DIAGNOSTIC ONLY'
                  + (f'   MISSING {miss}' if miss else ''))
        else:
            print(f'    {name:<16} WARN no markers present, skipped')

    if any(len(resolve_panel(g, sub)[2]) for g in panel.values()):
        print('    NOTE a shortened panel can be carried by whichever gene')
        print('    survives. If that gene is not lineage-specific the subtype')
        print('    call will be confidently wrong, which is worse than')
        print('    Unresolved. scripts/05c_panel_audit.py searches var_names')
        print('    for rhesus names of the missing genes.')

    if diag_cols:
        dmeans = sub.obs.groupby('sub_leiden',
                                 observed=True)[diag_cols].mean()
        dmeans.columns = [c.replace('sig_dx_', '') for c in dmeans.columns]
        print('\n  DIAGNOSTIC scores per subcluster (never assigned from):')
        print(dmeans.round(3).to_string())
        save_csv(dmeans.round(4), f'tier2_{tag}_diagnostic_means.csv')
        if {'GranzymeKM', 'GranzymeABH'} <= set(dmeans.columns):
            skew = (dmeans['GranzymeABH'] - dmeans['GranzymeKM']).round(3)
            print('\n  granzyme skew (ABH minus KM), the RIRA discriminator:')
            print(skew.to_string())
            print('  positive => TEM-like / CD16+ NK; negative => TCM/TTM-like')
            print('  or CD16-negative NK. Use this to promote the')
            print('  Effector/Cytotoxic family into a memory stage.')

    means = sub.obs.groupby('sub_leiden', observed=True)[score_cols].mean()
    means.columns = [c.replace('sig_', '') for c in means.columns]
    print('\n  mean signature score per subcluster:')
    print(means.round(3).to_string())
    save_csv(means.round(4), f'tier2_{tag}_subcluster_signature_means.csv')

    # --- which panels actually distinguish anything, reported not filtered
    # A panel whose mean barely moves between subclusters cannot tell them
    # apart. This is printed so a subtype carried by a flat panel is visible.
    # It gates NOTHING: once assignment is on the mean, a flat panel cannot
    # out-rank a real one anyway, and a floor here would only be a second
    # arbitrary threshold manufacturing Unresolved cells.
    rng_ = (means.max() - means.min()).sort_values(ascending=False)
    print('\n  between-subcluster range of each panel '
          '(how much signal it carries):')
    print(rng_.round(3).to_string())
    print('  Diagnostic only. A small range means that panel cannot')
    print('  distinguish these subclusters, whatever it scores.')

    # --- which panels move together, the evidence behind PANEL_FAMILIES
    if means.shape[0] >= 3 and means.shape[1] >= 2:
        corr = means.corr().round(2)
        print('\n  correlation of panel means across subclusters:')
        print(corr.to_string())
        print('  Two panels correlating near 1 cannot be assigned')
        print('  independently. That is what PANEL_FAMILIES encodes.')

    # --- build family means -----------------------------------------------
    fams = PANEL_FAMILIES.get(lin, {})
    grouped, members = {}, set()
    for fname, parts in fams.items():
        have = [p for p in parts if p in means.columns]
        if len(have) >= 1:
            grouped[fname] = means[have].mean(axis=1)
            members.update(have)
    for c in means.columns:
        if c not in members:
            grouped[c] = means[c]
    FAM = pd.DataFrame(grouped)
    if fams:
        print('\n  panel families applied: '
              + '; '.join(f'{k} = {"+".join(v)}' for k, v in fams.items()))
        print('  family mean per subcluster:')
        print(FAM.round(3).to_string())

    # --- assign on the subcluster mean, gate on the margin ----------------
    # score_genes subtracts a control-gene mean, so a non-expressing panel
    # sits near zero by construction and these raw means are comparable with
    # no rescaling. This is the table printed above, so the log and the
    # assignment now agree.
    rows = []
    for cl in FAM.index:
        srt = FAM.loc[cl].sort_values(ascending=False)
        top, topv = srt.index[0], float(srt.iloc[0])
        second = srt.index[1] if len(srt) > 1 else None
        margin = topv - float(srt.iloc[1]) if len(srt) > 1 else np.inf
        ok = margin >= MIN_SCORE_MARGIN and (topv > 0 or not REQUIRE_POSITIVE_TOP)
        rows.append({'top_family': top, 'top_mean': round(topv, 3),
                     'runner_up': second, 'margin': round(margin, 3),
                     'provisional_subtype': top if ok else 'Unresolved'})
    tbl = pd.DataFrame(rows, index=FAM.index)
    tbl['n_cells'] = sub.obs['sub_leiden'].value_counts()

    # How much signal the WINNING panel actually carries across subclusters.
    # Without it, "Treg, margin 0.005" reads as if it meant something when
    # the Treg panel's whole range is 0.100, i.e. noise. The row is still
    # correctly Unresolved; this column makes the reason visible.
    fam_range = (FAM.max() - FAM.min())
    tbl['top_family_range'] = tbl.top_family.map(fam_range).round(3)

    # --- sub-split a family on a measured diagnostic -----------------------
    split = FAMILY_SUBSPLIT.get(lin, {})
    for fam_name, cfg in split.items():
        hi, lo = cfg['hi_panel'], cfg['lo_panel']
        if hi not in sub.obs.columns or lo not in sub.obs.columns:
            print(f'\n  WARN  cannot sub-split {fam_name}: '
                  f'{hi} or {lo} was not scored')
            continue
        gh = sub.obs.groupby('sub_leiden', observed=True)[hi].mean()
        gl = sub.obs.groupby('sub_leiden', observed=True)[lo].mean()
        skew = (gh - gl)
        hit = tbl.provisional_subtype == fam_name
        if not hit.any():
            continue
        tbl['subsplit_skew'] = skew.round(3)
        for cl in tbl.index[hit]:
            v = float(skew.get(cl, 0.0))
            if v >= cfg['threshold']:
                tbl.loc[cl, 'provisional_subtype'] = cfg['hi_label']
            elif v <= -cfg['threshold']:
                tbl.loc[cl, 'provisional_subtype'] = cfg['lo_label']
        print(f'\n  sub-split {fam_name} on ({hi} minus {lo}) '
              f'at +/-{cfg["threshold"]}:')
        print(tbl.loc[hit, ['n_cells', 'subsplit_skew',
                            'provisional_subtype']].to_string())
        print('  A split on a measured diagnostic is stable across rebuilds;')
        print('  keying SUBTYPE_RENAME on subcluster ID would not be.')

    # --- dominant_frac, reported as HOMOGENEITY and gating nothing --------
    if REPORT_DOMINANT_FRAC:
        Z = sub.obs[score_cols].copy()
        Z = (Z - Z.mean()) / (Z.std() + 1e-9)
        sub.obs['provisional_subtype'] = np.asarray(
            Z.idxmax(axis=1).str.replace('sig_', '', regex=False))
        gb = sub.obs.groupby('sub_leiden', observed=True)['provisional_subtype']
        tbl['vote_winner'] = gb.agg(lambda x: x.value_counts().index[0])
        tbl['dominant_frac'] = gb.agg(
            lambda x: round(x.value_counts(normalize=True).iloc[0], 3))

    n_unres = int((tbl.provisional_subtype == 'Unresolved').sum())
    print(f'\n  subtype per subcluster (mean-based, Unresolved below margin '
          f'{MIN_SCORE_MARGIN}'
          + (' or on a negative top score' if REQUIRE_POSITIVE_TOP else '')
          + '):')
    cols_show = ['n_cells', 'top_family', 'top_mean', 'top_family_range',
                 'runner_up', 'margin', 'provisional_subtype']
    if REPORT_DOMINANT_FRAC:
        cols_show += ['vote_winner', 'dominant_frac']
    print(tbl[cols_show].to_string())

    if REPORT_DOMINANT_FRAC:
        print('\n  dominant_frac is the per-cell vote share and it GATES')
        print('  NOTHING. It measures homogeneity: a low value with a large')
        print('  margin means one population whose cells split their vote')
        print('  between co-expressed panels, while a low value with a small')
        print('  margin means a genuinely mixed subcluster. Those are the')
        print('  ones on the bubble and this column is how you find them.')
        mixed = tbl[(tbl.dominant_frac < 0.40) & (tbl.margin < MIN_SCORE_MARGIN)]
        if len(mixed):
            print(f'  on the bubble ({len(mixed)}): '
                  f'{", ".join(str(i) for i in mixed.index)}')

    if n_unres:
        print(f'\n  {n_unres} subcluster(s) Unresolved. That is a real answer,')
        print('  not a failure: the signatures do not separate them at this')
        print('  resolution. Lower SUB_RES or fix the panels before naming.')

    # --- data-driven top markers, so the signature call can be checked
    #     against genes nobody chose in advance
    try:
        sc.tl.rank_genes_groups(sub, groupby='sub_leiden', method='wilcoxon',
                                n_genes=N_TOP_MARKERS, use_raw=False)
        nm = sub.uns['rank_genes_groups']['names']
        top = pd.DataFrame({cl: [nm[cl][i] for i in range(N_TOP_MARKERS)]
                            for cl in sub.obs['sub_leiden'].cat.categories})
        print(f'\n  top {N_TOP_MARKERS} markers per subcluster:')
        print(top.to_string(index=False))
        save_csv(top, f'tier2_{tag}_top_markers.csv', index=False)
    except Exception as e:                                          # noqa: BLE001
        print(f'  WARN  rank_genes_groups skipped: {type(e).__name__}: {e}')

    # --- per-cell attributes that deliberately stay OFF the labels
    sub.obs['prolif_counts'] = gset_sum(sub, PROLIF_GENES, layer=COUNTS_LAYER)
    sub.obs['is_proliferating'] = sub.obs['prolif_counts'] >= POS_MIN
    act = present(ACTIVATION_GENES, sub)
    if act:
        sc.tl.score_genes(sub, gene_list=act, score_name='sig_Activation',
                          random_state=RANDOM_STATE)
        print('\n  sig_Activation scored on immediate-early genes. This is')
        print('  partly a tissue-dissociation artifact, so it stays a soft')
        print('  score and NEVER becomes a label.')

    if lin == 'T cells':
        c8 = gset_sum(sub, CD8_GENES, layer=COUNTS_LAYER)
        c4 = gset_sum(sub, CD4_GENES, layer=COUNTS_LAYER)
        sub.obs['cd48_rna'] = np.where(
            (c8 >= POS_MIN) & (c4 < POS_MIN), 'CD8',
            np.where((c4 >= POS_MIN) & (c8 < POS_MIN), 'CD4',
                     np.where((c4 >= POS_MIN) & (c8 >= POS_MIN),
                              'double-pos', 'double-neg')))
        print('\n  CD4/CD8 by RNA:')
        print(sub.obs['cd48_rna'].value_counts().to_string())

        if HAS_ADT:
            i8 = resolve_adt(adt_names, 'CD8')
            i4 = resolve_adt(adt_names, 'CD4')
            if i8 is not None and i4 is not None:
                clr = sub.obsm['protein_clr']
                sub.obs['cd48_adt'] = np.where(
                    clr[:, i8] > clr[:, i4], 'CD8', 'CD4')
                print('\n  CD4/CD8 by ADT (relative CLR, the reliable call):')
                print(sub.obs['cd48_adt'].value_counts().to_string())
                print('\n  ADT CD8 fraction per subcluster:')
                print((sub.obs.groupby('sub_leiden', observed=True)['cd48_adt']
                       .apply(lambda s: round((s == 'CD8').mean(), 3))
                       ).to_string())
            else:
                print('\n  WARN  CD4 or CD8 ADT tag not resolvable; '
                      'RNA call only.')

        # The naive pool stays ONE label carrying CD4/CD8 as an attribute,
        # rather than being split into CD4-naive and CD8-naive here.
        print('\n  CD4/CD8 is carried as a per-cell attribute, not folded')
        print('  into the subtype name.')

    # --- QC rules. Explicit thresholds on measurable quantities, never
    #     cluster IDs.
    qc_lab = {}
    if 'pct_counts_mt' in sub.obs.columns:
        mt_med = sub.obs.groupby('sub_leiden',
                                 observed=True)['pct_counts_mt'].median()
        for cl, v in mt_med.items():
            if v > QC_MITO_SUBCLUSTER:
                qc_lab[cl] = 'QC: high-mito'
        print(f'\n  median pct_counts_mt per subcluster '
              f'(flagged above {QC_MITO_SUBCLUSTER}%):')
        print(mt_med.round(2).to_string())

    own, foreign = QC_DOUBLET_PAIRS.get(lin, ([], []))
    if own and foreign:
        o = gset_sum(sub, own, layer=COUNTS_LAYER) >= POS_MIN
        f = gset_sum(sub, foreign, layer=COUNTS_LAYER) >= POS_MIN
        sub.obs['_coexpr'] = o & f
        co = sub.obs.groupby('sub_leiden', observed=True)['_coexpr'].mean()
        if HAS_DBL:
            ds = sub.obs.groupby('sub_leiden',
                                 observed=True)['doublet_score'].mean()
            ds_med = float(sub.obs['doublet_score'].median())
            thr = max(ds_med * QC_DOUBLET_SCORE_MULT, 1e-9)
        else:
            ds = pd.Series(np.nan, index=co.index)
            thr = np.inf
        dbl = pd.DataFrame({'coexpr_frac': co.round(3),
                            'mean_doublet_score': ds.round(4)})
        print(f'\n  doublet evidence per subcluster '
              f'(own {own} with foreign {foreign}):')
        print(dbl.to_string())
        for cl in co.index:
            hit_co = co[cl] >= QC_DOUBLET_COEXPR_FRAC
            hit_ds = HAS_DBL and ds[cl] >= thr
            # Co-expression alone is enough when doublet scores are unusable,
            # since Scrublet is the thing that failed, not the biology.
            if hit_co and (hit_ds or not HAS_DBL):
                qc_lab[cl] = 'QC: doublet-like'
        save_csv(dbl, f'tier2_{tag}_doublet_evidence.csv')

    # --- assemble the final subtype per subcluster
    final = {}
    for cl in tbl.index:
        if cl in qc_lab:
            final[cl] = qc_lab[cl]
            continue
        p = tbl.loc[cl, 'provisional_subtype']
        final[cl] = SUBTYPE_RENAME.get(f'{lin}:{p}', f'{tag} {p}')
    tbl['final_subtype'] = pd.Series(final)
    tbl['qc_flag'] = tbl.index.map(lambda c: c in qc_lab)
    print('\n  FINAL subtype per subcluster:')
    print(tbl.to_string())
    save_csv(tbl, f'tier2_{tag}_subtype_table.csv')
    if not SUBTYPE_RENAME:
        print('\n  SUBTYPE_RENAME is empty, so these names came from the')
        print('  signatures rather than from curation. Read the tables above,')
        print('  then fill SUBTYPE_RENAME in Cell 1 to promote them.')

    # --- enrichment, to check a subtype is not just one animal
    keep = ~sub.obs[ANIMAL_COL].astype(str).isin(['Unknown'])
    if keep.sum():
        _, e = log2_enrichment(sub.obs.loc[keep, 'sub_leiden'],
                               sub.obs.loc[keep, ANIMAL_COL].astype(str))
        print('\n  log2 enrichment by animal (clipped +/-2):')
        print(e.round(2).to_string())
        save_csv(e.round(3), f'tier2_{tag}_animal_enrichment.csv')
    for col in (TIME_COL, TISSUE_COL):
        if col in sub.obs.columns:
            _, e = log2_enrichment(sub.obs['sub_leiden'],
                                   sub.obs[col].astype(str))
            print(f'\n  log2 enrichment by {col}:')
            print(e.round(2).to_string())
            save_csv(e.round(3), f'tier2_{tag}_{col}_enrichment.csv')

    # --- write back onto the global object
    idx = adata.obs_names.get_indexer(sub.obs_names)
    adata.obs.iloc[idx, adata.obs.columns.get_loc('sub_leiden')] = \
        [f'{tag}_{c}' for c in sub.obs['sub_leiden'].astype(str)]
    mapped = sub.obs['sub_leiden'].astype(str).map(final).values
    adata.obs.iloc[idx, adata.obs.columns.get_loc('tier2_subtype')] = mapped
    adata.obs.iloc[idx, adata.obs.columns.get_loc('tier2_dominant_frac')] = \
        sub.obs['sub_leiden'].astype(str).map(tbl['dominant_frac']).values

    qc_cells = np.asarray(pd.Series(mapped).str.startswith('QC:'), dtype=bool)
    adata.obs.iloc[idx, adata.obs.columns.get_loc('qc_exclude')] = qc_cells
    adata.obs.iloc[idx, adata.obs.columns.get_loc('qc_reason')] = \
        np.where(qc_cells, mapped, '')

    def propagate(src, src_col, dst_col):
        """Write a per-lineage column into the global obs at the right rows.

        src is passed in rather than closed over, because the loop ends with
        `del sub` and a closure would break the moment this function were
        called anywhere else. The fill value matches the source dtype:
        initialising a bool column with '' and then assigning bools into it
        produces an object column holding both, which anndata will not always
        write to h5ad.
        """
        col = src.obs[src_col]
        if dst_col not in adata.obs.columns:
            k = col.dtype.kind
            adata.obs[dst_col] = (np.nan if k == 'f'
                                  else False if k == 'b'
                                  else 0 if k in 'iu'
                                  else '')
        adata.obs.iloc[idx, adata.obs.columns.get_loc(dst_col)] = col.values

    for c in ('provisional_subtype', 'is_proliferating', 'sig_Activation',
              'cd48_rna', 'cd48_adt', 'prolif_counts'):
        if c in sub.obs.columns:
            propagate(sub, c, c)

    # Signature columns are prefixed with the lineage tag on the global
    # object. Without it T-cell 'Naive' and B-cell 'Naive' both resolve to
    # obs['sig_Naive'] and overwrite each other's rows, leaving one column
    # that silently means two different things.
    for c in score_cols + diag_cols:
        propagate(sub, c, f'{tag}_{c}')
    print(f'\n  signature columns written as {tag}_sig_* on the global object'
          + (f' ({len(diag_cols)} of them diagnostic)' if diag_cols else ''))

    # per-lineage UMAP, since the global one cannot show subclusters
    try:
        fig, ax = plt.subplots(figsize=(16, 14))
        sc.pl.umap(sub, color='sub_leiden', ax=ax, show=False,
                   legend_loc='on data', frameon=False, size=12)
        ax.set_title(f'{lin}: sub-Leiden {SUB_RES}', fontsize=34)
        savefig(fig, f'F05_tier2_{tag}_subclusters')
    except Exception as e:                                          # noqa: BLE001
        print(f'  WARN  subcluster UMAP: {type(e).__name__}: {e}')

    del sub


# %% ============================================================
# CELL 13 — TIER 2 SUMMARY
# ============================================================
hr('CELL 13 — Tier 2 summary')

adata.obs['tier2_subtype'] = adata.obs['tier2_subtype'].astype(str)
print('  cells per Tier 2 subtype:')
print(adata.obs['tier2_subtype'].value_counts().to_string())

n_qc = int(adata.obs['qc_exclude'].sum())
print(f'\n  qc_exclude = True: {n_qc:,} cells '
      f'({100 * adata.obs["qc_exclude"].mean():.1f}%)')
if n_qc:
    print(adata.obs.loc[adata.obs['qc_exclude'],
                        'qc_reason'].value_counts().to_string())
print('  FLAGGED, not deleted. Downstream biology filters on qc_exclude; a')
print('  QC-excluded cell carrying a viral barcode is still the only')
print('  evidence we have for that barcode.')

for col in (TIME_COL, ANIMAL_COL):
    t = pd.crosstab(adata.obs['tier2_subtype'], adata.obs[col])
    save_csv(t, f'tier2_by_{col}.csv')
print()
print(pd.crosstab(adata.obs['tier2_subtype'],
                  adata.obs[TIME_COL]).to_string())

if HAS_TCR:
    tt = pd.crosstab(adata.obs['tier2_subtype'], adata.obs['has_tcr'])
    print('\n  TCR detection per subtype:')
    print(tt.to_string())
    save_csv(tt, 'tier2_tcr_by_subtype.csv')


# %% ============================================================
# CELL 14 — VIRAL CELLS ON THE ANNOTATED OBJECT
# ============================================================
hr('CELL 14 — viral cells')

if 'shiv_umi_cellranger' not in adata.obs.columns:
    print('  no shiv_umi_cellranger column; skipping')
else:
    v = adata.obs['shiv_umi_cellranger'].values.astype(float)
    if 'total_counts' in adata.obs.columns:
        tot = adata.obs['total_counts'].values.astype(float)
    else:
        Xc = adata.layers[COUNTS_LAYER]
        tot = np.asarray(Xc.sum(axis=1)).ravel()
    adata.obs['shiv_detected'] = v > 0
    adata.obs['viral_frac'] = np.where(tot > 0, v / tot, 0.0)
    n = int(adata.obs['shiv_detected'].sum())
    print(f'  cells with >=1 viral UMI : {n}')
    print(f'  total viral UMI          : {v.sum():,.0f}')

    # --- the background ceiling, from cells that cannot be infected --------
    bg_lin = [c for c in VIRAL_BACKGROUND_LINEAGES
              if (adata.obs['tier1_celltype'].astype(str) == c).any()]
    bg_mask = adata.obs['tier1_celltype'].astype(str).isin(bg_lin).values
    bg_hits = v[bg_mask & (v > 0)]
    print(f'\n  background lineages present: {bg_lin} '
          f'({int(bg_mask.sum()):,} cells)')
    print('  These cannot be productively infected. An erythrocyte has no')
    print('  nucleus and no CD4, so their viral burden IS the background.')

    if SHIV_CONFIDENT_MIN_UMI is not None:
        ceiling = int(SHIV_CONFIDENT_MIN_UMI)
        print(f'  ceiling PINNED at {ceiling} UMI by SHIV_CONFIDENT_MIN_UMI')
    elif len(bg_hits):
        ceiling = int(bg_hits.max()) + 1
        print(f'  background hits: {len(bg_hits)} cell(s) at '
              f'{sorted(bg_hits.astype(int))} UMI')
        print(f'  ceiling DERIVED at {ceiling} UMI (background max + 1)')
        print(f'  FRAGILE: derived from {len(bg_hits)} cell(s). One more')
        print('  background hit above this moves the ceiling and changes')
        print('  which cells survive. Pin SHIV_CONFIDENT_MIN_UMI to freeze it.')
    else:
        ceiling = 2
        print('  no viral UMI in the background lineages, so no ceiling could')
        print(f'  be measured. Falling back to {ceiling} UMI, which is a')
        print('  guess, not a measurement. Treat confident calls with care.')

    hi_umi = v >= ceiling
    hi_frac = adata.obs['viral_frac'].values >= SHIV_CONFIDENT_MIN_FRAC
    conf = (hi_umi & hi_frac) if SHIV_CONFIDENT_REQUIRE_BOTH else (hi_umi | hi_frac)
    adata.obs['shiv_confident'] = conf & adata.obs['shiv_detected'].values
    adata.obs['shiv_call'] = np.where(
        adata.obs['shiv_confident'], 'confident',
        np.where(adata.obs['shiv_detected'], 'background', 'none'))
    n_conf = int(adata.obs['shiv_confident'].sum())
    print(f'\n  confident call: >={ceiling} UMI '
          f'{"AND" if SHIV_CONFIDENT_REQUIRE_BOTH else "OR"} '
          f'>={SHIV_CONFIDENT_MIN_FRAC:.0%} of transcriptome')
    print(f'    by burden only   : {int((hi_umi & ~hi_frac & (v > 0)).sum())}')
    print(f'    by fraction only : {int((hi_frac & ~hi_umi & (v > 0)).sum())}')
    print(f'    by both          : {int((hi_umi & hi_frac).sum())}')
    print(f'    CONFIDENT total  : {n_conf} of {n} detected')
    print(f'    background       : {n - n_conf}')

    if not n:
        print('\n  nothing detected; skipping the rest of this cell.')
    else:
        # --- burden distribution, the measurement that settles it ----------
        vp = v[v > 0]
        fp = adata.obs['viral_frac'].values[v > 0]
        print(f'\n  burden per detected cell: min {vp.min():.0f}  median '
              f'{np.median(vp):.0f}  mean {vp.mean():.2f}  max {vp.max():.0f}')
        print(f'  viral fraction: median {np.median(fp):.2e}  max {fp.max():.2e}')
        bands = pd.DataFrame({
            'min_viral_UMI': [1, 2, 3, 5, 10, 25],
            'n_cells': [int((vp >= b).sum()) for b in (1, 2, 3, 5, 10, 25)]})
        print('\n  cells at or above each burden:')
        print(bands.to_string(index=False))
        n_one = int((vp == 1).sum())
        print(f'\n  {n_one} of {n} detected cells carry EXACTLY ONE molecule.')
        print('  One molecule is a Poisson draw from background, not a')
        print('  detection of infection. That is why shiv_confident exists.')

        cols = ['library', ANIMAL_COL, TIME_COL, 'shiv_umi_cellranger',
                'viral_frac', 'shiv_call', 'tier1_celltype', 'tier2_subtype',
                'clusters', 'sub_leiden', 'qc_exclude']
        if TISSUE_COL in adata.obs.columns:
            cols.insert(3, TISSUE_COL)
        for extra in ('has_tcr', 'cd3_counts', 'tier1_markers_detected'):
            if extra in adata.obs.columns:
                cols.append(extra)
        viral = adata.obs.loc[adata.obs['shiv_detected'], cols].copy()
        viral.index.name = 'barcode'
        viral = viral.sort_values('shiv_umi_cellranger', ascending=False)
        save_csv(viral.reset_index(), 'shiv_positive_cells.csv', index=False)

        print('\n  by Tier 1 lineage, DETECTED vs CONFIDENT:')
        byl = pd.DataFrame({
            'n_cells': adata.obs['tier1_celltype'].value_counts(),
            'detected': viral['tier1_celltype'].value_counts(),
            'confident': adata.obs.loc[adata.obs['shiv_confident'],
                                       'tier1_celltype'].value_counts(),
        }).fillna(0).astype(int)
        print(byl.to_string())
        save_csv(byl, 'shiv_by_lineage.csv')
        print('\n  The confident column is the one to read. The detected')
        print('  column is dominated by single-molecule background and will')
        print('  track lineage abundance rather than tropism.')

        if n_conf:
            print(f'\n  the {n_conf} confident cell(s):')
            print(viral[viral.shiv_call == 'confident'].to_string())

        # --- the weak ambient null, kept but labelled honestly -------------
        frac_l = adata.obs['tier1_celltype'].value_counts(normalize=True)
        cmp = pd.DataFrame({
            'observed_detected': viral['tier1_celltype'].value_counts(),
            'pct_of_object': (100 * frac_l).round(2)}).fillna(0)
        cmp['expected'] = (frac_l * n).round(2)
        cmp['fold'] = (cmp.observed_detected / cmp.expected).round(2)
        print('\n  detected cells against equal-per-cell contamination:')
        print(cmp.sort_values('pct_of_object', ascending=False).to_string())
        save_csv(cmp, 'shiv_vs_ambient_null.csv')
        print('\n  CAVEAT on that null, and it is a real one. It assumes every')
        print('  cell has the same chance of picking up a viral molecule.')
        print('  05b measured that assumption failing: Erythrocytes carried')
        print('  viral UMI at 0.61% of cells against 0.049% for everything')
        print('  else, a 12-fold excess in a lineage that cannot be infected.')
        print('  Ambient tracks how much SOUP a droplet holds, and a damaged')
        print('  or empty-ish droplet holds mostly soup. So treat this table')
        print('  as a rough orientation and the burden ceiling as the test.')

        # --- per animal WITHIN each library --------------------------------
        # 05b found library H holding 34 of 59 detected cells with nearly all
        # of them in animal 40702. Without per-animal cell counts inside each
        # library that is uninterpretable: it could be tropism, a genuinely
        # high-burden animal, or 40702 simply being most of H. This is that
        # test. Ambient is per GEM well, so under background the hits should
        # spread across the animals in a well in proportion to their cells.
        print('\n  per animal WITHIN each library (ambient is per GEM well,')
        print('  so background spreads across a well in proportion to cells):')
        rows = []
        for lib in sorted(adata.obs['library'].astype(str).unique()):
            lm = (adata.obs['library'].astype(str) == lib).values
            nvl = int(adata.obs['shiv_detected'].values[lm].sum())
            if nvl == 0:
                continue
            ser = adata.obs.loc[lm, ANIMAL_COL].astype(str)
            for an, cnt in ser.value_counts().items():
                am = lm & (adata.obs[ANIMAL_COL].astype(str) == an).values
                o = int(adata.obs['shiv_detected'].values[am].sum())
                oc = int(adata.obs['shiv_confident'].values[am].sum())
                share = cnt / lm.sum()
                rows.append({
                    'library': lib, 'animal': an, 'n_cells': int(cnt),
                    'pct_of_library': round(100 * share, 1),
                    'viral_cells_in_library': nvl,
                    'expected': round(share * nvl, 2), 'observed': o,
                    'fold': round(o / (share * nvl), 2) if share > 0 else np.nan,
                    'confident': oc,
                    'viral_UMI': float(v[am].sum()),
                })
        pa = pd.DataFrame(rows)
        if len(pa):
            pa['p_greater'] = [
                round(stats.binomtest(int(r.observed), int(r.viral_cells_in_library),
                                      r.pct_of_library / 100,
                                      alternative='greater').pvalue, 4)
                for r in pa.itertuples()]
            print(pa.sort_values(['library', 'observed'],
                                 ascending=[True, False]).to_string(index=False))
            save_csv(pa, 'shiv_per_animal_within_library.csv', index=False)
            sig = pa[(pa.p_greater < 0.05) & (pa.observed >= 3)]
            if len(sig):
                print('\n  animals carrying MORE than their share of a well:')
                print(sig[['library', 'animal', 'observed', 'expected',
                           'fold', 'p_greater']].to_string(index=False))
                print('  Uncorrected p, and these tests are not independent')
                print('  of each other. Read it as where to look, not as a')
                print('  result. An animal genuinely carrying more virus and')
                print('  a demultiplexing artifact look identical here.')
            else:
                print('\n  No animal carries more than its share of its own')
                print('  well. That is consistent with background, and it')
                print('  answers the 40702 concentration in library H: it')
                print('  tracks how much of H that animal is.')
        else:
            print('  no per-animal breakdown available')

        print('\n  FOR THE RECORD, because 06 and any figure depend on it:')
        print('  per-cell attribution of these molecules is NOT publishable.')
        print('  The ambient rate is high enough to put virus in erythrocytes.')
        print('  What survives is the barcode inventory per animal, since the')
        print('  ambient RNA still came out of real infected cells somewhere')
        print('  in that animal. Use shiv_confident for anything cell-level,')
        print('  and say in the methods that cell-level attribution was not')
        print('  attempted below the background ceiling.')


# %% ============================================================
# CELL 15 — FIGURES
# ============================================================
hr('CELL 15 — figures')


def umap_panel(color, stem, title, basis='X_umap',
               legend_loc='right margin', palette=None):
    try:
        fig, ax = plt.subplots(figsize=(16, 14))
        sc.pl.embedding(adata, basis=basis, color=color, ax=ax, show=False,
                        legend_loc=legend_loc, frameon=False, size=3,
                        palette=palette)
        ax.set_title(title, fontsize=34)
        savefig(fig, stem)
    except Exception as e:                                          # noqa: BLE001
        print(f'  WARN  {stem}: {type(e).__name__}: {e}')


t1_pal = [TIER1_COLORS.get(c, '#000000')
          for c in adata.obs['tier1_celltype'].astype('category').cat.categories]

if 'X_umap_prebbknn' in adata.obsm:
    umap_panel(BATCH_KEY, 'F05_umap_prebbknn_by_library',
               f'Before BBKNN, by {BATCH_KEY}', basis='X_umap_prebbknn')
umap_panel(BATCH_KEY, 'F05_umap_by_library', f'After BBKNN, by {BATCH_KEY}')
umap_panel(ANIMAL_COL, 'F05_umap_by_animal', 'By animal')
umap_panel(TIME_COL, 'F05_umap_by_timepoint', 'By timepoint')
umap_panel('clusters', 'F05_umap_clusters', f'Leiden {LEIDEN_RES}',
           legend_loc='on data')
umap_panel('tier1_celltype', 'F05_umap_tier1', 'Tier 1, corrected',
           legend_loc='on data', palette=t1_pal)
umap_panel('tier2_subtype', 'F05_umap_tier2', 'Tier 2 subtypes',
           legend_loc='on data')
umap_panel('dx_cd3_pos', 'F05_umap_cd3_gate', 'CD3 hard gate')
umap_panel('qc_exclude', 'F05_umap_qc_exclude', 'QC-excluded (flagged only)')
if HAS_TCR:
    umap_panel('has_tcr', 'F05_umap_has_tcr', 'TCR detected')
if HAS_DBL:
    umap_panel('doublet_score', 'F05_umap_doublet_score', 'Scrublet score')
if 'shiv_detected' in adata.obs.columns:
    umap_panel('shiv_umi_cellranger', 'F05_umap_shiv_umi', 'Viral UMI per cell')
    # The graded call, so no figure can show 38 single-molecule cells as
    # infected just because a bare >0 threshold called them positive.
    umap_panel('shiv_call', 'F05_umap_shiv_call',
               'Viral call: none / background / confident')
for tag in ADT_OVERLAY:
    i = resolve_adt(adt_names, tag) if HAS_ADT else None
    if i is not None:
        umap_panel(f'ADT_{adt_names[i]}', f'F05_umap_ADT_{tag}',
                   f'ADT {tag} (CLR)')


# %% ============================================================
# CELL 16 — SAVE
# ============================================================
hr('CELL 16 — save')

adata.write(OUT_H5AD)
print(f'  wrote adata_annotated.h5ad  ({adata.n_obs:,} cells x '
      f'{adata.n_vars:,} genes)')

hr('DONE — 05_embed_annotate.py')
print(f'  outputs: {OUT_DIR}')
print('    adata_annotated.h5ad            the object 06 continues from')
print('    adata_tier1_checkpoint.h5ad     resume point, recipe-stamped')
print('    tier1_transition_matrix.csv     uncorrected vs corrected')
print('    tier1_marker_support.csv        which labels are actually backed')
print('    tier1_tcr_recovery_matrix.csv   pre-TCR vs post-TCR')
print('    tier2_<lineage>_subtype_table.csv  the review table to read')
print('    shiv_by_lineage.csv             detected vs confident per lineage')
print('    shiv_per_animal_within_library.csv  the 40702 question')
print()
print('  Read these five things before 06:')
print('    1. The CD3 background calibration in CELL 8. Background CD3+ in')
print('       Erythrocytes and Platelets should be near zero at the chosen')
print('       threshold. If it is not, the gate is moving soup into T.')
print('    2. tier1_marker_support.csv. B cells should pass cleanly; NK and')
print('       the small myeloid lineages are where cells get reassigned.')
print('    3. shiv_by_lineage.csv, the CONFIDENT column only. The detected')
print('       column is single-molecule background and tracks abundance.')
print('    4. shiv_per_animal_within_library.csv. This answers whether the')
print('       library-H concentration in animal 40702 is real or is just')
print('       how much of H that animal is.')
print('    5. tier2_*_subtype_table.csv. Unresolved subclusters need a')
print('       decision; dominant_frac near 0.40 means a split subcluster.')
print()
print('  Then fill SUBTYPE_RENAME in Cell 1 and re-run with')
print('  RESUME_FROM_TIER1 = True to apply the curated names in minutes.')
print('  The checkpoint is recipe-stamped, so a resume after changing any')
print('  Tier 1 setting is refused rather than silently reapplying the old')
print('  gate.')
print()
print('  ON THE VIRAL ARM, so it is written down and not just in a chat log:')
print('  do not use per-cell viral attribution below shiv_confident in any')
print('  figure or claim. The ambient rate puts virus in erythrocytes. The')
print('  barcode inventory per animal is the part that survives, because the')
print('  ambient RNA still came from real infected cells in that animal.')
print()
print('  Next: 06_viral_barcode.py')
