# CMLRE Deep-Sea eDNA Pipeline — Full Architecture

**Project:** SIH PS DJS_26_SW_15 — AI-driven taxonomy identification and biodiversity
assessment from deep-sea eDNA, without heavy reliance on reference databases.

**Status key used throughout:** ✅ Implemented & verified | 🔶 Implemented, not yet
merged/verified | ⏳ Designed, not built | 🅿️ Parked / side-branch

---

## 0. High-level pipeline flow

```
Raw FASTQ (paired-end)
        │
        ▼
[Stage 0] DADA2 Preprocessing (R)                              ✅
        │  filterAndTrim → learnErrors → dada → mergePairs →
        │  makeSequenceTable → removeBimeraDenovo
        ▼
ASV table (sequence + abundance) + FASTA
        │
        ▼
[Stage 1] Embedding Backbone (DNABERT-S, pretrained)            ✅
        │  patched to bypass FlashAttention/Triton incompatibility
        │  (pure PyTorch attention fallback)
        ▼
embeddings.npy (one 768-dim vector per ASV)
        │
        ├─────────────────────────────────────────────┐
        ▼                                              ▼
[Stage 2] Reference Index Build (one-time, offline)    │        ✅
        │  PR2 FASTA + taxonomy → same DNABERT-S       │
        │  model → FAISS IVF index                     │
        ▼                                              │
reference.faiss + reference_metadata.json               │
        │                                              │
        │                                              ▼
        │                              [Stage 3] Clustering (HDBSCAN)   ✅ (v2)
        │                                      │  L2-normalize → UMAP →  🔶 (v3, UMAP added)
        │                                      │  HDBSCAN                
        │                                      │  + Stability scoring    🔶 (v4, in progress)
        │                                      │  + Bootstrap reproducibility 🔶 (v4)
        │                                      ▼
        │                              Cluster-labeled ASVs +
        │                              per-cluster representative
        │                              + Stability_Score +
        │                              Bootstrap_Reproducibility_Pct
        │                                      │
        └──────────────────┬───────────────────┘
                            ▼
        [Stage 4] Species Annotation (per cluster representative only)   ✅
                │  PRIMARY: FAISS top-5 nearest-neighbor search           
                │    + k-NN genus-level agreement cross-check             🔶 (crossover #3)
                │  SECONDARY (optional): BLAST cross-validation           ✅ (optional flag)
                │  Empirically-calibrated confidence thresholds           ✅ (crossover #1)
                ▼
        Cluster_ID → Best_Match, Similarity_Pct, Status,
                     kNN_Agreement_Pct
                            │
                            ▼
        [Stage 5] Abundance + Biodiversity Metrics                       ✅
                │  Per-cluster read totals, % of sample
                │  Shannon, Simpson, Chao1 richness estimator
                ▼
        [Stage 6] Final Report                                          ✅
                │  final_report.csv (main deliverable)
                │  candidate_novel_taxa.csv
                │  biodiversity_summary.json
                ▼
        [Stage 7] Dashboard / Visualization                              ⏳
                (UMAP Explorer, QC/Disagreement Review pages —
                 referenced by teammate's checklist, not yet built/seen)
```

---

## 1. Stage-by-stage detail

### Stage 0 — DADA2 Preprocessing ✅
- **Script:** `run_dada2.R`
- **Input:** raw paired-end FASTQ (real deep-sea eDNA sample, e.g. ERR3466765,
  CMLRE/eDNAbyss project, ENA accession PRJEB33873/PRJEB37673)
- **Steps:** `filterAndTrim()` (maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE) →
  `learnErrors()` (forward + reverse, run once, ~55 min each unbatched on CPU) →
  `dada()` denoising → `mergePairs()` → `makeSequenceTable()` →
  `removeBimeraDenovo()` (chimera removal)
- **Output:** ASV table (unique sequence + read count), FASTA of unique
  chimera-free sequences
- **Verified real numbers (single sample):** 200,056 raw reads → 166,184 post
  N-filter → 103,268 successfully merged → 67,865 total reads in final ASV
  table (2,192 ASVs) after chimera removal
- **Known scaling note:** for multi-sample runs, `learnErrors()` should be run
  ONCE on pooled samples (standard DADA2 workflow), not per-sample — not yet
  implemented, single-sample only so far

### Stage 1 — Embedding Backbone (DNABERT-S, pretrained) ✅
- **Script:** `s2.py`
- **Model:** DNABERT-S, pretrained, local weights (no fine-tuning, in the
  adopted primary pipeline)
- **Fix required to run at all:** `transformers==4.35.2` + `tokenizers==0.14.1`
  pinned (modern `transformers` defaults to meta-device model construction,
  incompatible with DNABERT-S's custom MosaicBERT-style `__init__` code,
  causing `RuntimeError: Tensor on device meta is not on the expected device
  cpu!`). Requires a full kernel restart after installing these versions —
  a pip install alone does not take effect on an already-imported library.
- **Also required:** manual monkey-patch bypassing FlashAttention/Triton
  entirely, replaced with a pure PyTorch attention implementation
  (`pytorch_flash_attn()`), since the original Triton kernel path crashes on
  this GPU/environment combination independent of the meta-device issue.
- **Output:** `embeddings.npy` — one 768-dimensional vector per ASV
- **Known performance issue:** if run unbatched (one sequence at a time),
  processing 240,000 reference sequences was estimated at ~172 hours.
  Batching is required for practical runtime; documented as "10-30x speedup"
  expected but not confirmed fixed in the reviewed `s2.py` version.

### Stage 2 — Reference Index Build (one-time, offline) ✅
- **Script:** `build_reference_index.py`
- **Input:** PR2 v5.1.1 SSU reference FASTA + taxonomy TSV (240,201 sequences,
  capped/sampled to 240,000)
- **Process:** embeds every reference sequence with the SAME DNABERT-S model
  used for query embeddings (critical — query and reference embeddings must
  share the same vector space), builds a FAISS IVF (or HNSW) approximate
  nearest-neighbor index
- **Output:** `reference.faiss`, `reference_metadata.json` (seq_id → lineage
  string), `reference_config.json`
- **Run once per reference database version** — not re-embedded per sample
- **Verified:** real run completed successfully using the transformers
  version fix; produced a working 240,000-sequence, 768-dim IVF index

### Stage 3 — Clustering (HDBSCAN) ✅ base / 🔶 crossover-enhanced
- **Script:** `recluster_improved.py`
- **Input:** `embeddings.npy` (full per-ASV embeddings) — does NOT re-run
  DNABERT-S, reuses existing embeddings
- **v1 (original, via `s2.py`):** raw HDBSCAN on unnormalized embeddings —
  598 clusters, high "Unclustered" noise bucket, 3.4% confident-match rate
  downstream
- **v2 (improved):** L2-normalize embeddings before clustering (Euclidean
  distance → effectively cosine distance), `min_samples` set lower than
  `min_cluster_size` for more lenient clustering — improved to 1,024
  clusters, confident-match rate rose to 7.6%
- **v3 (crossover feature — UMAP before HDBSCAN):** 🔶 added a UMAP
  dimensionality-reduction step (768-dim → ~30-dim) before HDBSCAN, to
  counter the "curse of dimensionality" degrading density-based clustering
  in high-dim space. Original full-dim embeddings preserved separately for
  FAISS annotation (only the clustering step operates in UMAP-reduced
  space). Real run result: 55 clusters, 143 unclustered on the full 19,857-ASV
  set with these particular UMAP settings — a substantial change from v2,
  not yet fully validated against v2 as strictly "better" (fewer, tighter
  clusters vs. many small ones — tradeoff not yet analyzed).
- **v4 (crossover feature — stability + bootstrap reproducibility):** 🔶
  in progress. Two additions:
  - **Stability scoring:** surfaces HDBSCAN's native
    `cluster_persistence_` score per cluster (previously computed
    internally but unused/undisplayed)
  - **Bootstrap reproducibility:** re-runs UMAP+HDBSCAN on N=10 random 80%
    subsamples of the data; for each original cluster, checks whether a
    majority-overlapping cluster reappears in each bootstrap run; reports
    `Bootstrap_Reproducibility_Pct` per cluster. Intended to give the
    "candidate novel taxa" claim actual reproducibility evidence rather
    than resting on a single clustering run.
  - **Known cost:** bootstrap step re-runs the full UMAP+HDBSCAN pipeline
    10 additional times: real added runtime, roughly ~10x a single run.

### Stage 4 — Species Annotation ✅ base / 🔶 crossover-enhanced
- **Script:** `step3_annotate_species.py`
- **Design principle (deliberate, do not change):** annotation runs ONLY on
  one representative sequence per cluster, not on every individual read —
  directly serves the problem statement's "reduce computational time"
  requirement.
- **PRIMARY — FAISS embedding-similarity search:**
  - Original (v1): top-1 nearest-neighbor search only, thresholds applied
    directly to cosine similarity score
  - **Crossover feature #1 — Confidence calibration:** ✅ thresholds
    (Confident ≥0.80, Uncertain ≥0.65, else Novel) empirically calibrated
    against the ACTUAL observed similarity score distribution on the full
    240K PR2 database (max ~85%, 95th percentile ~80%, median ~72%) —
    not naive BLAST-style 95%/80% cutoffs, which were tried first and
    incorrectly flagged 93% of clusters as "novel" before recalibration.
  - **Crossover feature #3 — k-NN top-5 cross-check:** 🔶 implemented.
    Extended search from top-1 to top-5 nearest neighbors. Computes
    genus-level agreement (not full-lineage-string match, which
    over-penalized legitimate `_sp.` unresolved-species placeholder
    entries) between the top-1 hit and its 4 next-nearest neighbors. A
    "Confident match" whose neighbors mostly disagree at genus level gets
    downgraded to "Low confidence / possible divergent (neighbor
    disagreement)". Real verified result on actual data: of 80 raw
    top-1 "Confident match" clusters, 46 were downgraded under this check
    (both full-string and genus-level comparison gave nearly identical
    downgrade counts — 47 vs 46 — indicating the disagreement is genuine,
    not a placeholder-naming artifact).
- **SECONDARY (optional) — BLAST cross-validation:** ✅ implemented as an
  opt-in flag (`--enable-blast-crossvalidation`), not run by default.
  Independent secondary check against the same or a different reference
  database, using traditional alignment rather than embedding similarity.
- **Output columns:** `Cluster_ID`, `Best_Match` (full lineage string),
  `Similarity_Pct`, `Status`, `kNN_Agreement_Pct`, optionally
  `BLAST_Best_Match` / `BLAST_Similarity_Pct` / `BLAST_Status` if enabled

### Stage 5 — Abundance + Biodiversity Metrics ✅
- **Function:** `compute_abundance()`, `compute_diversity_metrics()` in
  `step3_annotate_species.py`
- Per-cluster total reads and % of sample
- Shannon diversity index, Simpson diversity index, Chao1 richness
  estimator (using singleton/doubleton counts), computed on
  taxonomically-assigned (non-"Unclustered") reads only
- Real verified output (one full pipeline run): species richness 1,023,
  Shannon 5.86, Simpson 0.992, Chao1 1,023.0, 144,776 total reads used

### Stage 6 — Final Report ✅
- **Outputs:** `final_report.csv` (main deliverable — merges annotation +
  abundance + ASV counts per cluster), `candidate_novel_taxa.csv` (subset
  where Status == "Candidate novel taxon"), `biodiversity_summary.json`
- Real verified run: 1,024 clusters total → 80 raw "Confident match" (7.8%,
  closely matching a previously reported 7.6% headline figure), 214 "Candidate
  novel taxon" (20.9%, closely matching a previously reported 21.9% figure),
  730 "Low confidence / possible divergent" (71.3% — the largest bucket,
  and the primary target of the k-NN crossover feature)

### Stage 7 — Dashboard / Visualization ⏳
- Referenced in a teammate-provided implementation checklist as already
  existing ("UMAP Explorer," "QC/Disagreement Review" pages) — not yet
  reviewed or confirmed built as of this document. Open item.
- **Planned/implied scope (from crossover doc + checklist):**
  - UMAP 2D projection view (separate from the ~30-dim UMAP used for
    clustering — a dedicated 2D projection purely for visualization)
  - Per-cluster stability/reproducibility display
  - Within-cluster disagreement flagging view (sample member reads, check
    agreement with centroid label) — mentioned in teammate's checklist as
    "Critical," not yet implemented anywhere in reviewed code
  - Benchmark comparison view (this pipeline vs. BLAST/QIIME2 baseline)

---

## 2. Side-branch — Fine-tuned DNABERT-S 🅿️ (parked, not merged)

**Status: fine-tuning run completed successfully; NOT yet validated or
integrated into the main pipeline. Treated as an experimental alternative,
not a replacement, until proven better on real numbers.**

- **Objective:** adapt DNABERT-S's embedding space specifically to this
  project's taxonomic domain (PR2 clade structure, 18S/COI marker genes),
  rather than relying on off-the-shelf pretrained embeddings with no
  domain adaptation.
- **Method:** triplet loss / contrastive fine-tuning
  - Anchor = a training read
  - Positive = another read from the same known genus (PR2)
  - Negative = a **hard negative** — a read from a different genus within
    the SAME family (taxonomically close but distinct), not a random
    unrelated sequence. Hard negatives sharpen fine-grained separation
    more effectively than easy/random negatives.
  - Loss: `max(0, d(anchor,pos) - d(anchor,neg) + margin)` (standard
    triplet margin loss, margin=1.0)
- **Training data:** PR2 stratified train split (Epic A4) — genera
  thresholded to ≥6 sequences for a mechanically valid split (238,009 rows,
  1,501 genera after dropping 925 rare genera); unresolved/organelle
  placeholder labels (`_X`, `:mito`, etc.) deliberately included as real
  classes.
- **Real training run (completed):** 30,000 triplets, 5 epochs, batch
  size 4 (documented GPU-OOM-safe ceiling — not increased), single T4 GPU
  (dual-GPU DataParallel explicitly tried previously and found SLOWER, not
  adopted), ~9.5-10.4 hour estimated/actual runtime, fp32 (no mixed
  precision — deliberately avoided due to a previously-documented AMP-related
  crash from leftover GPU memory on an earlier attempt)
- **Checkpointing:** saved every epoch (5 full checkpoints) plus every 1,000
  steps mid-epoch (additional safety net for a multi-hour run) — many
  checkpoint folders exist; the canonical usable output is the final,
  no-suffix `dnabert-s-finetuned` folder (post-epoch-5 weights)
- **Known incompatibility, fixed identically to the pretrained pipeline:**
  same meta-device `RuntimeError` on model load in any fresh Kaggle session
  without the `transformers==4.35.2`/`tokenizers==0.14.1` pin + kernel
  restart; same FlashAttention/Triton monkey-patch required.
- **What re-validating this requires (not yet done):**
  1. Re-embed the ASV sample set using the fine-tuned model (Stage 1, rerun)
  2. Re-embed the PR2 reference database using the fine-tuned model
     (Stage 2, rerun — same ~172-hour unbatched risk applies if batching
     was never fixed in `s2.py`)
  3. Re-cluster (Stage 3) and re-annotate (Stage 4) on the new embeddings
  4. Compare confident-match %, novel %, Shannon index, and cluster count
     side-by-side against the pretrained-model baseline (v2/v3/v4 results)
  5. Adopt as primary ONLY if genuinely better/equal on real numbers — per
     explicit team decision, a working pretrained pipeline with real
     numbers is not to be replaced on a theoretical/unvalidated
     improvement alone

---

## 3. Original two-model architecture — deprecated, superseded ⏳/🅿️

**Historical note, included for completeness — this was the FIRST design
attempted, before the cluster-first pipeline was adopted. It never produced
working output. Documented here so the "ideal pipeline" picture is
complete, not because it is in use.**

- **Design:** classify-first, not cluster-first
  - **Model 1:** supervised classifier — cross-entropy loss, temperature
    scaling for calibrated confidence, class-weighted/focal loss for rare
    taxa, k-NN cross-check, per-rank (genus vs. species) confidence
    thresholds
  - **Model 2:** unsupervised clustering, applied only to low-confidence
    leftovers from Model 1 (the reverse order of the adopted pipeline,
    which clusters FIRST then annotates representatives)
- **Why abandoned as the primary path:** repeatedly blocked on the same
  meta-device loading crash during initial development, then GPU OOM
  crashes during training attempts, DataParallel found slower than single
  GPU, batch size capped at 4, a crash at `model.to(DEVICE)` from leftover
  GPU memory after a prior AMP failure — never reached a working, complete
  training run under this original design.
- **Explicit decision:** do NOT reintroduce classify-first ordering or
  per-read (as opposed to per-cluster-representative) annotation into the
  adopted pipeline — both current design choices (cluster-first,
  representative-only annotation) are considered deliberate improvements
  over this original design, not compromises.

---

## 4. What's fully verified vs. designed-but-unbuilt (honest summary)

| Component | Status | Evidence |
|---|---|---|
| DADA2 preprocessing | ✅ Verified | Real run, real numbers, single sample |
| DNABERT-S embedding (pretrained) | ✅ Verified | Real embeddings.npy produced |
| FAISS reference index build | ✅ Verified | Real 240K-sequence index built |
| HDBSCAN clustering v2 (L2-norm) | ✅ Verified | Real numbers: 1024 clusters |
| UMAP clustering v3 | 🔶 Run once | Real numbers produced, not yet compared/validated as better than v2 |
| Confidence calibration | ✅ Verified | Built into working annotation, real thresholds |
| k-NN top-5 cross-check | ✅ Verified | Real downgrade counts on real data (46-47 clusters) |
| BLAST secondary cross-check | ✅ Built | Present as opt-in flag, not run in the verified end-to-end pass |
| Abundance + diversity metrics | ✅ Verified | Real Shannon/Simpson/Chao1 numbers |
| Stability scoring | 🔶 Code written | Not yet run/verified on real data |
| Bootstrap reproducibility | 🔶 Code written | Not yet run/verified on real data (real runtime cost, ~10x a single clustering run) |
| Fine-tuned DNABERT-S | 🅿️ Trained | Model exists, NOT yet re-validated through the full pipeline |
| Dashboard (UMAP Explorer, QC/Disagreement Review) | ⏳ Unbuilt/unseen | Referenced by teammate's checklist only |
| Multi-sample DADA2 (pooled error learning) | ⏳ Not implemented | Currently single-sample only |
| Original two-model architecture | 🅿️ Abandoned | Never produced working output, superseded |

---

## 5. Reference: threshold/parameter values actually used

| Parameter | Value | Where |
|---|---|---|
| DADA2 `maxEE` | c(2,2) | Stage 0 |
| DADA2 `truncQ` | 2 | Stage 0 |
| HDBSCAN `min_cluster_size` | 5 | Stage 3 |
| HDBSCAN `min_samples` | 2 | Stage 3 |
| UMAP `n_components` | 30 | Stage 3 (v3+) |
| UMAP `n_neighbors` | 15 | Stage 3 (v3+) |
| UMAP `min_dist` | 0.0 | Stage 3 (v3+) |
| Bootstrap runs (N) | 10 | Stage 3 (v4) |
| Bootstrap subsample fraction | 80% | Stage 3 (v4) |
| FAISS confident threshold | 0.80 (80% similarity) | Stage 4 |
| FAISS uncertain threshold | 0.65 (65% similarity) | Stage 4 |
| k-NN cross-check K | 5 | Stage 4 |
| k-NN agreement comparison level | Genus (not full lineage string) | Stage 4 |
| Fine-tuning triplets | 30,000 | Side-branch |
| Fine-tuning epochs | 5 | Side-branch |
| Fine-tuning batch size | 4 | Side-branch |
| Fine-tuning learning rate | 1e-5 | Side-branch |
| Fine-tuning margin (triplet loss) | 1.0 | Side-branch |
